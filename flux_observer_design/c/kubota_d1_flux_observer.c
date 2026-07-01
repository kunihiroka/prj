#include "kubota_d1_flux_observer.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct KubotaD1Derived {
    float sigma;
    float l_sigma_h;
    float rho;
    float a_r11;
    float a_r12;
    float a_i12;
    float a_r21;
    float a_r22;
    float a_i22;
    float b1;
    float g1;
    float g2;
    float g3;
    float g4;
} KubotaD1Derived;

static int kd1_is_finite_positive(float x)
{
    return isfinite(x) && (x > 0.0f);
}

static FluxObserverStatus kd1_fetch_config(KubotaD1FluxObserver *observer, FluxObserverMotorConfig *cfg)
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
    if (!kd1_is_finite_positive(cfg->rs_ohm) ||
        !kd1_is_finite_positive(cfg->rr_ohm) ||
        !kd1_is_finite_positive(cfg->lls_h) ||
        !kd1_is_finite_positive(cfg->llr_h) ||
        !kd1_is_finite_positive(cfg->lm_h) ||
        (cfg->pole_pairs == 0u) ||
        !kd1_is_finite_positive(cfg->control_period_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static FluxObserverStatus kd1_derive(
    const FluxObserverMotorConfig *cfg,
    float omega_m_e,
    float k_gain,
    KubotaD1Derived *d)
{
    float ls_h;
    float lr_h;
    float det_l;
    float tau_r;
    float alpha_rel;

    if ((cfg == NULL) || (d == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!isfinite(k_gain) || (k_gain <= 0.0f)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    ls_h = cfg->lls_h + cfg->lm_h;
    lr_h = cfg->llr_h + cfg->lm_h;
    det_l = ls_h * lr_h - cfg->lm_h * cfg->lm_h;
    if (!isfinite(det_l) || (fabsf(det_l) < 1.0e-12f)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    d->sigma = 1.0f - cfg->lm_h * cfg->lm_h / (ls_h * lr_h);
    d->l_sigma_h = det_l / lr_h;
    d->rho = cfg->lm_h / lr_h;
    tau_r = lr_h / cfg->rr_ohm;
    if (!kd1_is_finite_positive(d->sigma) ||
        !kd1_is_finite_positive(tau_r) ||
        !kd1_is_finite_positive(d->l_sigma_h)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    d->a_r11 = -(cfg->rs_ohm / (d->sigma * ls_h) + cfg->rr_ohm * (1.0f - d->sigma) / (d->sigma * lr_h));
    d->a_r12 = cfg->lm_h / (d->sigma * ls_h * lr_h * tau_r);
    d->a_i12 = -cfg->lm_h * omega_m_e / (d->sigma * ls_h * lr_h);
    d->a_r21 = cfg->lm_h / tau_r;
    d->a_r22 = -1.0f / tau_r;
    d->a_i22 = omega_m_e;
    d->b1 = 1.0f / (d->sigma * ls_h);

    alpha_rel = -d->sigma * ls_h * lr_h / cfg->lm_h;
    d->g1 = (k_gain - 1.0f) * (d->a_r11 + d->a_r22);
    d->g2 = (k_gain - 1.0f) * d->a_i22;
    d->g3 = (k_gain * k_gain - 1.0f) * (-alpha_rel * d->a_r11 + d->a_r21);
    d->g3 += alpha_rel * (k_gain - 1.0f) * (d->a_r11 + d->a_r22);
    d->g4 = alpha_rel * (k_gain - 1.0f) * d->a_i22;

    if (!isfinite(d->a_r11) || !isfinite(d->a_r12) || !isfinite(d->a_i12) ||
        !isfinite(d->a_r21) || !isfinite(d->a_r22) || !isfinite(d->a_i22) ||
        !isfinite(d->b1) || !isfinite(d->g1) || !isfinite(d->g2) ||
        !isfinite(d->g3) || !isfinite(d->g4)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static float kd1_clamp_denominator(float denominator, float min_abs)
{
    if (fabsf(denominator) >= min_abs) {
        return denominator;
    }
    return (denominator >= 0.0f) ? min_abs : -min_abs;
}

static void kd1_fill_output(
    const KubotaD1FluxObserver *observer,
    const KubotaD1Derived *d,
    float omega_m_e,
    float omega_k,
    KubotaD1FluxObserverOutput *output)
{
    if (output == NULL) {
        return;
    }
    output->isd_hat_a = observer->isd_hat_a;
    output->isq_hat_a = observer->isq_hat_a;
    output->phi_r_hat_wb = observer->phi_r_hat_wb;
    output->psi_sd_hat_wb = d->l_sigma_h * observer->isd_hat_a + d->rho * observer->phi_r_hat_wb;
    output->psi_sq_hat_wb = d->l_sigma_h * observer->isq_hat_a;
    output->omega_m_e_rad_s = omega_m_e;
    output->omega_k_rad_s = omega_k;
    output->omega_slip_rad_s = omega_k - omega_m_e;
    output->k_gain = observer->k_gain;
    output->g1 = d->g1;
    output->g2 = d->g2;
    output->g3 = d->g3;
    output->g4 = d->g4;
}

void KubotaD1FluxObserver_Init(KubotaD1FluxObserver *observer, FluxObserverApi api)
{
    if (observer == NULL) {
        return;
    }
    memset(observer, 0, sizeof(*observer));
    observer->api = api;
    observer->k_gain = 1.2f;
    observer->min_denominator_wb = 1.0e-6f;
}

void KubotaD1FluxObserver_SetK(KubotaD1FluxObserver *observer, float k_gain)
{
    if (observer == NULL) {
        return;
    }
    if (!isfinite(k_gain) || (k_gain <= 0.0f)) {
        return;
    }
    observer->k_gain = k_gain;
}

void KubotaD1FluxObserver_SetMinDenominator(KubotaD1FluxObserver *observer, float min_denominator_wb)
{
    if (observer == NULL) {
        return;
    }
    if (!kd1_is_finite_positive(min_denominator_wb)) {
        return;
    }
    observer->min_denominator_wb = min_denominator_wb;
}

FluxObserverStatus KubotaD1FluxObserver_Reset(
    KubotaD1FluxObserver *observer,
    float isd_hat_a,
    float isq_hat_a,
    float phi_r_hat_wb)
{
    if (observer == NULL) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!isfinite(isd_hat_a) || !isfinite(isq_hat_a) || !isfinite(phi_r_hat_wb)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    observer->isd_hat_a = isd_hat_a;
    observer->isq_hat_a = isq_hat_a;
    observer->phi_r_hat_wb = phi_r_hat_wb;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus KubotaD1FluxObserver_ResetFromCurrents(
    KubotaD1FluxObserver *observer,
    float isd_a,
    float isq_a)
{
    FluxObserverMotorConfig cfg;
    FluxObserverStatus status;

    status = kd1_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    observer->last_config = cfg;
    observer->isd_hat_a = isd_a;
    observer->isq_hat_a = isq_a;
    observer->phi_r_hat_wb = cfg.lm_h * isd_a;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus KubotaD1FluxObserver_Step(
    KubotaD1FluxObserver *observer,
    const KubotaD1FluxObserverInput *input,
    KubotaD1FluxObserverOutput *output)
{
    FluxObserverMotorConfig cfg;
    KubotaD1Derived d;
    float omega_m_e;
    float ed;
    float eq;
    float denominator;
    float omega_k;
    float disd;
    float disq;
    float dphi_r;
    FluxObserverStatus status;

    if ((observer == NULL) || (input == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    status = kd1_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    omega_m_e = (float)cfg.pole_pairs * input->omega_m_rad_s;
    status = kd1_derive(&cfg, omega_m_e, observer->k_gain, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    ed = observer->isd_hat_a - input->isd_a;
    eq = observer->isq_hat_a - input->isq_a;
    denominator = kd1_clamp_denominator(observer->phi_r_hat_wb, observer->min_denominator_wb);
    omega_k = d.a_i22 + (d.a_r21 * observer->isq_hat_a + d.g4 * ed + d.g3 * eq) / denominator;

    disd =
        d.a_r11 * observer->isd_hat_a +
        omega_k * observer->isq_hat_a +
        d.a_r12 * observer->phi_r_hat_wb +
        d.b1 * input->usd_v +
        d.g1 * ed -
        d.g2 * eq;
    disq =
        d.a_r11 * observer->isq_hat_a -
        omega_k * observer->isd_hat_a +
        d.a_i12 * observer->phi_r_hat_wb +
        d.b1 * input->usq_v +
        d.g2 * ed +
        d.g1 * eq;
    dphi_r =
        d.a_r21 * observer->isd_hat_a +
        d.a_r22 * observer->phi_r_hat_wb +
        d.g3 * ed -
        d.g4 * eq;

    observer->isd_hat_a += cfg.control_period_s * disd;
    observer->isq_hat_a += cfg.control_period_s * disq;
    observer->phi_r_hat_wb += cfg.control_period_s * dphi_r;
    observer->last_config = cfg;

    if (!isfinite(observer->isd_hat_a) || !isfinite(observer->isq_hat_a) || !isfinite(observer->phi_r_hat_wb)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    kd1_fill_output(observer, &d, omega_m_e, omega_k, output);
    return FLUX_OBSERVER_OK;
}
