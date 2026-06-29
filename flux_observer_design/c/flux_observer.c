#include "flux_observer.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct FoComplex {
    float re;
    float im;
} FoComplex;

typedef struct FoDerived {
    float ls_h;
    float lr_h;
    float det_l;
    float cs0;
    float cs1;
    float cr0;
    float cr1;
} FoDerived;

static FoComplex fo_c(float re, float im)
{
    FoComplex z;
    z.re = re;
    z.im = im;
    return z;
}

static FoComplex fo_c_add(FoComplex a, FoComplex b)
{
    return fo_c(a.re + b.re, a.im + b.im);
}

static FoComplex fo_c_sub(FoComplex a, FoComplex b)
{
    return fo_c(a.re - b.re, a.im - b.im);
}

static FoComplex fo_c_mul(FoComplex a, FoComplex b)
{
    return fo_c(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
}

static FoComplex fo_c_scale(FoComplex a, float k)
{
    return fo_c(k * a.re, k * a.im);
}

static float fo_c_abs2(FoComplex a)
{
    return a.re * a.re + a.im * a.im;
}

static FoComplex fo_c_div(FoComplex a, FoComplex b)
{
    float den = fo_c_abs2(b);
    return fo_c((a.re * b.re + a.im * b.im) / den, (a.im * b.re - a.re * b.im) / den);
}

static int fo_is_finite_positive(float x)
{
    return isfinite(x) && (x > 0.0f);
}

static FluxObserverStatus fo_fetch_config(FluxObserver *observer, FluxObserverMotorConfig *cfg)
{
    if ((observer == NULL) || (cfg == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (observer->api.get_motor_config == NULL) {
        return FLUX_OBSERVER_ERR_API;
    }
    if (observer->api.get_motor_config(observer->api.user, cfg) != 0) {
        return FLUX_OBSERVER_ERR_API;
    }
    if (!fo_is_finite_positive(cfg->rs_ohm) ||
        !fo_is_finite_positive(cfg->rr_ohm) ||
        !fo_is_finite_positive(cfg->lls_h) ||
        !fo_is_finite_positive(cfg->llr_h) ||
        !fo_is_finite_positive(cfg->lm_h) ||
        (cfg->pole_pairs == 0u) ||
        !fo_is_finite_positive(cfg->control_period_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static FluxObserverStatus fo_derive(const FluxObserverMotorConfig *cfg, FoDerived *d)
{
    float det_l;
    if ((cfg == NULL) || (d == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    d->ls_h = cfg->lls_h + cfg->lm_h;
    d->lr_h = cfg->llr_h + cfg->lm_h;
    det_l = d->ls_h * d->lr_h - cfg->lm_h * cfg->lm_h;
    if (!isfinite(det_l) || (fabsf(det_l) < 1.0e-12f)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    d->det_l = det_l;
    d->cs0 = d->lr_h / det_l;
    d->cs1 = -cfg->lm_h / det_l;
    d->cr0 = -cfg->lm_h / det_l;
    d->cr1 = d->ls_h / det_l;
    return FLUX_OBSERVER_OK;
}

static FoComplex fo_current_from_flux(float c0, float c1, FoComplex psi_s, FoComplex psi_r)
{
    return fo_c(c0 * psi_s.re + c1 * psi_r.re, c0 * psi_s.im + c1 * psi_r.im);
}

static FluxObserverStatus fo_observer_h(
    const FluxObserver *observer,
    const FluxObserverMotorConfig *cfg,
    const FoDerived *d,
    float omega_r_rad_s,
    float omega_k_rad_s,
    FoComplex h[2])
{
    FoComplex a00;
    FoComplex a01;
    FoComplex a10;
    FoComplex a11;
    FoComplex adj00;
    FoComplex adj01;
    FoComplex adj10;
    FoComplex adj11;
    FoComplex row20;
    FoComplex row21;
    FoComplex lhs00;
    FoComplex lhs01;
    FoComplex lhs10;
    FoComplex lhs11;
    FoComplex rhs0;
    FoComplex rhs1;
    FoComplex det_lhs;
    FoComplex det_a;
    FoComplex trace_a;
    float bandwidth;
    float ratio;
    float desired_trace;
    float desired_det;

    if ((observer == NULL) || (cfg == NULL) || (d == NULL) || (h == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    bandwidth = observer->observer_bandwidth_rad_s;
    ratio = observer->observer_pole_ratio;
    if (!fo_is_finite_positive(bandwidth) || !fo_is_finite_positive(ratio)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    a00 = fo_c(-cfg->rs_ohm * d->cs0, -omega_k_rad_s);
    a01 = fo_c(-cfg->rs_ohm * d->cs1, 0.0f);
    a10 = fo_c(-cfg->rr_ohm * d->cr0, 0.0f);
    a11 = fo_c(-cfg->rr_ohm * d->cr1, -(omega_k_rad_s - omega_r_rad_s));

    adj00 = a11;
    adj01 = fo_c_scale(a01, -1.0f);
    adj10 = fo_c_scale(a10, -1.0f);
    adj11 = a00;

    row20 = fo_c_add(fo_c_scale(adj00, d->cs0), fo_c_scale(adj10, d->cs1));
    row21 = fo_c_add(fo_c_scale(adj01, d->cs0), fo_c_scale(adj11, d->cs1));

    lhs00 = fo_c(d->cs0, 0.0f);
    lhs01 = fo_c(d->cs1, 0.0f);
    lhs10 = row20;
    lhs11 = row21;

    desired_trace = -(1.0f + ratio) * bandwidth;
    desired_det = ratio * bandwidth * bandwidth;
    trace_a = fo_c_add(a00, a11);
    det_a = fo_c_sub(fo_c_mul(a00, a11), fo_c_mul(a01, a10));
    rhs0 = fo_c_sub(trace_a, fo_c(desired_trace, 0.0f));
    rhs1 = fo_c_sub(det_a, fo_c(desired_det, 0.0f));

    det_lhs = fo_c_sub(fo_c_mul(lhs00, lhs11), fo_c_mul(lhs01, lhs10));
    if (fo_c_abs2(det_lhs) < 1.0e-18f) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    h[0] = fo_c_div(fo_c_sub(fo_c_mul(rhs0, lhs11), fo_c_mul(lhs01, rhs1)), det_lhs);
    h[1] = fo_c_div(fo_c_sub(fo_c_mul(lhs00, rhs1), fo_c_mul(rhs0, lhs10)), det_lhs);
    return FLUX_OBSERVER_OK;
}

static void fo_store_h_matrix(FluxObserver *observer, const FoComplex h[2])
{
    observer->h[0][0] = h[0].re;
    observer->h[0][1] = -h[0].im;
    observer->h[1][0] = h[0].im;
    observer->h[1][1] = h[0].re;
    observer->h[2][0] = h[1].re;
    observer->h[2][1] = -h[1].im;
    observer->h[3][0] = h[1].im;
    observer->h[3][1] = h[1].re;
}

void FluxObserver_Init(FluxObserver *observer, FluxObserverApi api)
{
    if (observer == NULL) {
        return;
    }
    memset(observer, 0, sizeof(*observer));
    observer->api = api;
    observer->observer_bandwidth_rad_s = 2200.0f;
    observer->observer_pole_ratio = 1.55f;
}

void FluxObserver_SetPolePlacement(FluxObserver *observer, float bandwidth_rad_s, float pole_ratio)
{
    if (observer == NULL) {
        return;
    }
    if (fo_is_finite_positive(bandwidth_rad_s)) {
        observer->observer_bandwidth_rad_s = bandwidth_rad_s;
    }
    if (fo_is_finite_positive(pole_ratio)) {
        observer->observer_pole_ratio = pole_ratio;
    }
}

FluxObserverStatus FluxObserver_ResetFlux(
    FluxObserver *observer,
    float psi_sd_wb,
    float psi_sq_wb,
    float psi_rd_wb,
    float psi_rq_wb)
{
    if (observer == NULL) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    observer->psi_sd_wb = psi_sd_wb;
    observer->psi_sq_wb = psi_sq_wb;
    observer->psi_rd_wb = psi_rd_wb;
    observer->psi_rq_wb = psi_rq_wb;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus FluxObserver_ResetFromCurrents(
    FluxObserver *observer,
    float isd_a,
    float isq_a,
    float ird_a,
    float irq_a)
{
    FluxObserverMotorConfig cfg;
    FoDerived d;
    FluxObserverStatus status;

    status = fo_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    status = fo_derive(&cfg, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    observer->last_config = cfg;
    observer->psi_sd_wb = d.ls_h * isd_a + cfg.lm_h * ird_a;
    observer->psi_sq_wb = d.ls_h * isq_a + cfg.lm_h * irq_a;
    observer->psi_rd_wb = cfg.lm_h * isd_a + d.lr_h * ird_a;
    observer->psi_rq_wb = cfg.lm_h * isq_a + d.lr_h * irq_a;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus FluxObserver_Step(
    FluxObserver *observer,
    const FluxObserverInput *input,
    FluxObserverOutput *output)
{
    FluxObserverMotorConfig cfg;
    FoDerived d;
    FoComplex h[2];
    FoComplex psi_s;
    FoComplex psi_r;
    FoComplex i_s_hat;
    FoComplex i_r_hat;
    FoComplex innovation;
    FoComplex a00;
    FoComplex a01;
    FoComplex a10;
    FoComplex a11;
    FoComplex v_s;
    FoComplex dx_s;
    FoComplex dx_r;
    float omega_r;
    float omega_k;
    FluxObserverStatus status;

    if ((observer == NULL) || (input == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    status = fo_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    status = fo_derive(&cfg, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    omega_r = (float)cfg.pole_pairs * input->omega_m_rad_s;
    omega_k = omega_r + input->omega_slip_rad_s;

    status = fo_observer_h(observer, &cfg, &d, omega_r, omega_k, h);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    fo_store_h_matrix(observer, h);
    observer->last_config = cfg;

    psi_s = fo_c(observer->psi_sd_wb, observer->psi_sq_wb);
    psi_r = fo_c(observer->psi_rd_wb, observer->psi_rq_wb);
    i_s_hat = fo_current_from_flux(d.cs0, d.cs1, psi_s, psi_r);
    innovation = fo_c(input->isd_a - i_s_hat.re, input->isq_a - i_s_hat.im);

    a00 = fo_c(-cfg.rs_ohm * d.cs0, -omega_k);
    a01 = fo_c(-cfg.rs_ohm * d.cs1, 0.0f);
    a10 = fo_c(-cfg.rr_ohm * d.cr0, 0.0f);
    a11 = fo_c(-cfg.rr_ohm * d.cr1, -(omega_k - omega_r));
    v_s = fo_c(input->vsd_v, input->vsq_v);

    dx_s = fo_c_add(fo_c_add(fo_c_mul(a00, psi_s), fo_c_mul(a01, psi_r)), v_s);
    dx_s = fo_c_add(dx_s, fo_c_mul(h[0], innovation));
    dx_r = fo_c_add(fo_c_mul(a10, psi_s), fo_c_mul(a11, psi_r));
    dx_r = fo_c_add(dx_r, fo_c_mul(h[1], innovation));

    psi_s = fo_c_add(psi_s, fo_c_scale(dx_s, cfg.control_period_s));
    psi_r = fo_c_add(psi_r, fo_c_scale(dx_r, cfg.control_period_s));

    observer->psi_sd_wb = psi_s.re;
    observer->psi_sq_wb = psi_s.im;
    observer->psi_rd_wb = psi_r.re;
    observer->psi_rq_wb = psi_r.im;

    if (output != NULL) {
        i_s_hat = fo_current_from_flux(d.cs0, d.cs1, psi_s, psi_r);
        i_r_hat = fo_current_from_flux(d.cr0, d.cr1, psi_s, psi_r);
        output->psi_sd_wb = psi_s.re;
        output->psi_sq_wb = psi_s.im;
        output->psi_rd_wb = psi_r.re;
        output->psi_rq_wb = psi_r.im;
        output->isd_hat_a = i_s_hat.re;
        output->isq_hat_a = i_s_hat.im;
        output->ird_hat_a = i_r_hat.re;
        output->irq_hat_a = i_r_hat.im;
        output->omega_r_rad_s = omega_r;
        output->omega_k_rad_s = omega_k;
        memcpy(output->h, observer->h, sizeof(observer->h));
    }

    return FLUX_OBSERVER_OK;
}

FluxObserverStatus FluxObserver_GetLastH(
    const FluxObserver *observer,
    float h[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS])
{
    if ((observer == NULL) || (h == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    memcpy(h, observer->h, sizeof(observer->h));
    return FLUX_OBSERVER_OK;
}

const char *FluxObserver_StatusString(FluxObserverStatus status)
{
    switch (status) {
    case FLUX_OBSERVER_OK:
        return "OK";
    case FLUX_OBSERVER_ERR_NULL:
        return "NULL pointer";
    case FLUX_OBSERVER_ERR_API:
        return "API error";
    case FLUX_OBSERVER_ERR_PARAM:
        return "parameter error";
    case FLUX_OBSERVER_ERR_SINGULAR:
        return "singular observer equation";
    default:
        return "unknown error";
    }
}
