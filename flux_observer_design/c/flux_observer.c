#include "flux_observer.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct FoDerived {
    float ls_h;
    float lr_h;
    float det_l;
    float cs0;
    float cs1;
    float cr0;
    float cr1;
} FoDerived;

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

static void fo_primary_current_from_flux(
    const FoDerived *d,
    float psi_sd,
    float psi_sq,
    float psi_rd,
    float psi_rq,
    float *isd,
    float *isq)
{
    *isd = d->cs0 * psi_sd + d->cs1 * psi_rd;
    *isq = d->cs0 * psi_sq + d->cs1 * psi_rq;
}

static void fo_secondary_current_from_flux(
    const FoDerived *d,
    float psi_sd,
    float psi_sq,
    float psi_rd,
    float psi_rq,
    float *ird,
    float *irq)
{
    *ird = d->cr0 * psi_sd + d->cr1 * psi_rd;
    *irq = d->cr0 * psi_sq + d->cr1 * psi_rq;
}

static FluxObserverStatus fo_solve4(float aug[4][5], float x[4])
{
    unsigned int col;
    unsigned int row;
    unsigned int pivot_row;

    for (col = 0u; col < 4u; ++col) {
        float pivot_abs = fabsf(aug[col][col]);
        pivot_row = col;
        for (row = col + 1u; row < 4u; ++row) {
            float candidate = fabsf(aug[row][col]);
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot_row = row;
            }
        }
        if (!isfinite(pivot_abs) || (pivot_abs < 1.0e-18f)) {
            return FLUX_OBSERVER_ERR_SINGULAR;
        }
        if (pivot_row != col) {
            unsigned int k;
            for (k = col; k < 5u; ++k) {
                float tmp = aug[col][k];
                aug[col][k] = aug[pivot_row][k];
                aug[pivot_row][k] = tmp;
            }
        }
        {
            float pivot = aug[col][col];
            unsigned int k;
            for (k = col; k < 5u; ++k) {
                aug[col][k] /= pivot;
            }
        }
        for (row = 0u; row < 4u; ++row) {
            if (row != col) {
                float factor = aug[row][col];
                unsigned int k;
                for (k = col; k < 5u; ++k) {
                    aug[row][k] -= factor * aug[col][k];
                }
            }
        }
    }

    for (row = 0u; row < 4u; ++row) {
        x[row] = aug[row][4];
        if (!isfinite(x[row])) {
            return FLUX_OBSERVER_ERR_SINGULAR;
        }
    }
    return FLUX_OBSERVER_OK;
}

static FluxObserverStatus fo_observer_H(
    const FluxObserver *observer,
    const FluxObserverMotorConfig *cfg,
    const FoDerived *d,
    float omega_r_rad_s,
    float omega_k_rad_s,
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS])
{
    float a00_re;
    float a00_im;
    float a01_re;
    float a01_im;
    float a10_re;
    float a10_im;
    float a11_re;
    float a11_im;
    float row20_re;
    float row20_im;
    float row21_re;
    float row21_im;
    float desired_trace;
    float desired_det;
    float rhs0_re;
    float rhs0_im;
    float det_re;
    float det_im;
    float rhs1_re;
    float rhs1_im;
    float g[4];
    float aug[4][5];
    FluxObserverStatus status;

    if ((observer == NULL) || (cfg == NULL) || (d == NULL) || (H == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!fo_is_finite_positive(observer->observer_bandwidth_rad_s) ||
        !fo_is_finite_positive(observer->observer_pole_ratio)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    a00_re = -cfg->rs_ohm * d->cs0;
    a00_im = -omega_k_rad_s;
    a01_re = -cfg->rs_ohm * d->cs1;
    a01_im = 0.0f;
    a10_re = -cfg->rr_ohm * d->cr0;
    a10_im = 0.0f;
    a11_re = -cfg->rr_ohm * d->cr1;
    a11_im = -(omega_k_rad_s - omega_r_rad_s);

    row20_re = d->cs0 * a11_re - d->cs1 * a10_re;
    row20_im = d->cs0 * a11_im - d->cs1 * a10_im;
    row21_re = -d->cs0 * a01_re + d->cs1 * a00_re;
    row21_im = -d->cs0 * a01_im + d->cs1 * a00_im;

    desired_trace = -(1.0f + observer->observer_pole_ratio) * observer->observer_bandwidth_rad_s;
    desired_det = observer->observer_pole_ratio * observer->observer_bandwidth_rad_s * observer->observer_bandwidth_rad_s;
    rhs0_re = (a00_re + a11_re) - desired_trace;
    rhs0_im = a00_im + a11_im;
    det_re = (a00_re * a11_re - a00_im * a11_im) - (a01_re * a10_re - a01_im * a10_im);
    det_im = (a00_re * a11_im + a00_im * a11_re) - (a01_re * a10_im + a01_im * a10_re);
    rhs1_re = det_re - desired_det;
    rhs1_im = det_im;

    aug[0][0] = d->cs0;
    aug[0][1] = 0.0f;
    aug[0][2] = d->cs1;
    aug[0][3] = 0.0f;
    aug[0][4] = rhs0_re;

    aug[1][0] = 0.0f;
    aug[1][1] = d->cs0;
    aug[1][2] = 0.0f;
    aug[1][3] = d->cs1;
    aug[1][4] = rhs0_im;

    aug[2][0] = row20_re;
    aug[2][1] = -row20_im;
    aug[2][2] = row21_re;
    aug[2][3] = -row21_im;
    aug[2][4] = rhs1_re;

    aug[3][0] = row20_im;
    aug[3][1] = row20_re;
    aug[3][2] = row21_im;
    aug[3][3] = row21_re;
    aug[3][4] = rhs1_im;

    status = fo_solve4(aug, g);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    H[0][0] = g[0];
    H[0][1] = -g[1];
    H[1][0] = g[1];
    H[1][1] = g[0];
    H[2][0] = g[2];
    H[2][1] = -g[3];
    H[3][0] = g[3];
    H[3][1] = g[2];
    return FLUX_OBSERVER_OK;
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
    float isd_hat;
    float isq_hat;
    float ird_hat;
    float irq_hat;
    float innovation_d;
    float innovation_q;
    float dx_sd;
    float dx_sq;
    float dx_rd;
    float dx_rq;
    float omega_r;
    float omega_k;
    float omega_slip;
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
    omega_slip = omega_k - omega_r;

    status = fo_observer_H(observer, &cfg, &d, omega_r, omega_k, observer->H);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    observer->last_config = cfg;

    fo_primary_current_from_flux(
        &d,
        observer->psi_sd_wb,
        observer->psi_sq_wb,
        observer->psi_rd_wb,
        observer->psi_rq_wb,
        &isd_hat,
        &isq_hat);
    innovation_d = input->isd_a - isd_hat;
    innovation_q = input->isq_a - isq_hat;

    fo_secondary_current_from_flux(
        &d,
        observer->psi_sd_wb,
        observer->psi_sq_wb,
        observer->psi_rd_wb,
        observer->psi_rq_wb,
        &ird_hat,
        &irq_hat);

    dx_sd = input->vsd_v - cfg.rs_ohm * isd_hat + omega_k * observer->psi_sq_wb;
    dx_sq = input->vsq_v - cfg.rs_ohm * isq_hat - omega_k * observer->psi_sd_wb;
    dx_rd = -cfg.rr_ohm * ird_hat + omega_slip * observer->psi_rq_wb;
    dx_rq = -cfg.rr_ohm * irq_hat - omega_slip * observer->psi_rd_wb;

    dx_sd += observer->H[0][0] * innovation_d + observer->H[0][1] * innovation_q;
    dx_sq += observer->H[1][0] * innovation_d + observer->H[1][1] * innovation_q;
    dx_rd += observer->H[2][0] * innovation_d + observer->H[2][1] * innovation_q;
    dx_rq += observer->H[3][0] * innovation_d + observer->H[3][1] * innovation_q;

    observer->psi_sd_wb += cfg.control_period_s * dx_sd;
    observer->psi_sq_wb += cfg.control_period_s * dx_sq;
    observer->psi_rd_wb += cfg.control_period_s * dx_rd;
    observer->psi_rq_wb += cfg.control_period_s * dx_rq;

    if (output != NULL) {
        fo_primary_current_from_flux(
            &d,
            observer->psi_sd_wb,
            observer->psi_sq_wb,
            observer->psi_rd_wb,
            observer->psi_rq_wb,
            &isd_hat,
            &isq_hat);
        fo_secondary_current_from_flux(
            &d,
            observer->psi_sd_wb,
            observer->psi_sq_wb,
            observer->psi_rd_wb,
            observer->psi_rq_wb,
            &ird_hat,
            &irq_hat);
        output->psi_sd_wb = observer->psi_sd_wb;
        output->psi_sq_wb = observer->psi_sq_wb;
        output->psi_rd_wb = observer->psi_rd_wb;
        output->psi_rq_wb = observer->psi_rq_wb;
        output->isd_hat_a = isd_hat;
        output->isq_hat_a = isq_hat;
        output->ird_hat_a = ird_hat;
        output->irq_hat_a = irq_hat;
        output->omega_r_rad_s = omega_r;
        output->omega_k_rad_s = omega_k;
        memcpy(output->H, observer->H, sizeof(observer->H));
    }

    return FLUX_OBSERVER_OK;
}

FluxObserverStatus FluxObserver_GetLastH(
    const FluxObserver *observer,
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS])
{
    if ((observer == NULL) || (H == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    memcpy(H, observer->H, sizeof(observer->H));
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
