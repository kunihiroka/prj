#include "sled23_flux_observer.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct Sled23Derived {
    float lr_h;
    float det_l;
    float l_sigma_h;
    float l_m_inv_gamma_h;
    float rr_inv_gamma_ohm;
    float r_sigma_ohm;
    float alpha_rad_s;
} Sled23Derived;

static int sled_is_finite_positive(float x)
{
    return isfinite(x) && (x > 0.0f);
}

static FluxObserverStatus sled_fetch_config(Sled23FluxObserver *observer, FluxObserverMotorConfig *cfg)
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
    if (!sled_is_finite_positive(cfg->rs_ohm) ||
        !sled_is_finite_positive(cfg->rr_ohm) ||
        !sled_is_finite_positive(cfg->lls_h) ||
        !sled_is_finite_positive(cfg->llr_h) ||
        !sled_is_finite_positive(cfg->lm_h) ||
        (cfg->pole_pairs == 0u) ||
        !sled_is_finite_positive(cfg->control_period_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static FluxObserverStatus sled_derive(const FluxObserverMotorConfig *cfg, Sled23Derived *d)
{
    float ls_h;
    float lm_over_lr;

    if ((cfg == NULL) || (d == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    ls_h = cfg->lls_h + cfg->lm_h;
    d->lr_h = cfg->llr_h + cfg->lm_h;
    d->det_l = ls_h * d->lr_h - cfg->lm_h * cfg->lm_h;
    if (!isfinite(d->det_l) || (fabsf(d->det_l) < 1.0e-12f)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    lm_over_lr = cfg->lm_h / d->lr_h;
    d->l_sigma_h = d->det_l / d->lr_h;
    d->l_m_inv_gamma_h = cfg->lm_h * cfg->lm_h / d->lr_h;
    d->rr_inv_gamma_ohm = cfg->rr_ohm * lm_over_lr * lm_over_lr;
    d->r_sigma_ohm = cfg->rs_ohm + d->rr_inv_gamma_ohm;
    d->alpha_rad_s = d->rr_inv_gamma_ohm / d->l_m_inv_gamma_h;
    if (!sled_is_finite_positive(d->l_sigma_h) ||
        !sled_is_finite_positive(d->l_m_inv_gamma_h) ||
        !sled_is_finite_positive(d->rr_inv_gamma_ohm) ||
        !sled_is_finite_positive(d->r_sigma_ohm) ||
        !sled_is_finite_positive(d->alpha_rad_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static float sled_clamp_denominator(float denominator, float min_abs)
{
    if (fabsf(denominator) >= min_abs) {
        return denominator;
    }
    return (denominator >= 0.0f) ? min_abs : -min_abs;
}

static void sled_fill_output(
    const Sled23FluxObserver *observer,
    const Sled23Derived *d,
    float omega_m_e,
    float omega_s,
    float b,
    float gamma,
    float k1,
    float k2,
    Sled23FluxObserverOutput *output)
{
    if (output == NULL) {
        return;
    }
    output->isd_hat_a = observer->isd_hat_a;
    output->isq_hat_a = observer->isq_hat_a;
    output->psi_r_hat_wb = observer->psi_r_hat_wb;
    output->psi_sd_hat_wb = observer->psi_r_hat_wb + d->l_sigma_h * observer->isd_hat_a;
    output->psi_sq_hat_wb = d->l_sigma_h * observer->isq_hat_a;
    output->omega_m_e_rad_s = omega_m_e;
    output->omega_s_rad_s = omega_s;
    output->omega_slip_rad_s = omega_s - omega_m_e;
    output->alpha_rad_s = d->alpha_rad_s;
    output->alpha_i_rad_s = observer->alpha_i_rad_s;
    output->b_rad_s = b;
    output->gamma_rad_s = gamma;
    output->k1 = k1;
    output->k2 = k2;
}

void Sled23FluxObserver_Init(Sled23FluxObserver *observer, FluxObserverApi api)
{
    if (observer == NULL) {
        return;
    }
    memset(observer, 0, sizeof(*observer));
    observer->api = api;
    observer->alpha_i_rad_s = 2200.0f;
    observer->zeta_inf = 0.2f;
    observer->fixed_b_rad_s = 0.0f;
    observer->use_fixed_b = 0u;
    observer->min_denominator_wb = 1.0e-6f;
}

void Sled23FluxObserver_SetDesign(
    Sled23FluxObserver *observer,
    float alpha_i_rad_s,
    float zeta_inf)
{
    if (observer == NULL) {
        return;
    }
    if (!sled_is_finite_positive(alpha_i_rad_s) || !isfinite(zeta_inf) || (zeta_inf < 0.0f)) {
        return;
    }
    observer->alpha_i_rad_s = alpha_i_rad_s;
    observer->zeta_inf = zeta_inf;
    observer->use_fixed_b = 0u;
}

void Sled23FluxObserver_SetFixedB(
    Sled23FluxObserver *observer,
    float alpha_i_rad_s,
    float b_rad_s)
{
    if (observer == NULL) {
        return;
    }
    if (!sled_is_finite_positive(alpha_i_rad_s) || !sled_is_finite_positive(b_rad_s)) {
        return;
    }
    observer->alpha_i_rad_s = alpha_i_rad_s;
    observer->fixed_b_rad_s = b_rad_s;
    observer->use_fixed_b = 1u;
}

void Sled23FluxObserver_SetMinDenominator(Sled23FluxObserver *observer, float min_denominator_wb)
{
    if (observer == NULL) {
        return;
    }
    if (!sled_is_finite_positive(min_denominator_wb)) {
        return;
    }
    observer->min_denominator_wb = min_denominator_wb;
}

FluxObserverStatus Sled23FluxObserver_Reset(
    Sled23FluxObserver *observer,
    float isd_hat_a,
    float isq_hat_a,
    float psi_r_hat_wb)
{
    if (observer == NULL) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!isfinite(isd_hat_a) || !isfinite(isq_hat_a) || !isfinite(psi_r_hat_wb)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    observer->isd_hat_a = isd_hat_a;
    observer->isq_hat_a = isq_hat_a;
    observer->psi_r_hat_wb = psi_r_hat_wb;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus Sled23FluxObserver_ResetFromCurrents(
    Sled23FluxObserver *observer,
    float isd_a,
    float isq_a)
{
    FluxObserverMotorConfig cfg;
    Sled23Derived d;
    FluxObserverStatus status;

    status = sled_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    status = sled_derive(&cfg, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    observer->last_config = cfg;
    observer->isd_hat_a = isd_a;
    observer->isq_hat_a = isq_a;
    observer->psi_r_hat_wb = d.l_m_inv_gamma_h * isd_a;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus Sled23FluxObserver_Step(
    Sled23FluxObserver *observer,
    const Sled23FluxObserverInput *input,
    Sled23FluxObserverOutput *output)
{
    FluxObserverMotorConfig cfg;
    Sled23Derived d;
    float omega_m_e;
    float b;
    float gamma;
    float k_den;
    float k1;
    float k2;
    float eisd;
    float eisq;
    float denominator;
    float omega_s;
    float disd;
    float disq;
    float dpsi_r;
    FluxObserverStatus status;

    if ((observer == NULL) || (input == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    status = sled_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    status = sled_derive(&cfg, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    omega_m_e = (float)cfg.pole_pairs * input->omega_m_rad_s;
    b = observer->use_fixed_b ? observer->fixed_b_rad_s : (2.0f * observer->zeta_inf * fabsf(omega_m_e) + d.alpha_rad_s);
    gamma = observer->alpha_i_rad_s - d.alpha_rad_s;
    k_den = d.alpha_rad_s * d.alpha_rad_s + omega_m_e * omega_m_e;
    if (!sled_is_finite_positive(k_den) || !sled_is_finite_positive(observer->alpha_i_rad_s) || !sled_is_finite_positive(b)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    k1 = b * d.alpha_rad_s / k_den;
    k2 = b * omega_m_e / k_den;

    eisd = input->isd_a - observer->isd_hat_a;
    eisq = input->isq_a - observer->isq_hat_a;
    denominator = sled_clamp_denominator(
        observer->psi_r_hat_wb - d.l_sigma_h * eisd,
        observer->min_denominator_wb);
    omega_s = omega_m_e +
        (d.rr_inv_gamma_ohm * input->isq_a +
         k2 * observer->alpha_i_rad_s * d.l_sigma_h * eisd -
         gamma * d.l_sigma_h * eisq) /
            denominator;

    disd =
        (d.alpha_rad_s * observer->psi_r_hat_wb -
         d.r_sigma_ohm * observer->isd_hat_a +
         omega_s * d.l_sigma_h * observer->isq_hat_a +
         input->usd_v +
         d.l_sigma_h * (gamma * eisd - omega_m_e * eisq)) /
        d.l_sigma_h;
    disq =
        (-omega_m_e * observer->psi_r_hat_wb -
         d.r_sigma_ohm * observer->isq_hat_a -
         omega_s * d.l_sigma_h * observer->isd_hat_a +
         input->usq_v +
         d.l_sigma_h * (gamma * eisq + omega_m_e * eisd)) /
        d.l_sigma_h;
    dpsi_r =
        -d.alpha_rad_s * observer->psi_r_hat_wb +
        d.rr_inv_gamma_ohm * observer->isd_hat_a +
        (k1 * observer->alpha_i_rad_s - gamma) * d.l_sigma_h * eisd -
        (omega_s - omega_m_e) * d.l_sigma_h * eisq;

    observer->isd_hat_a += cfg.control_period_s * disd;
    observer->isq_hat_a += cfg.control_period_s * disq;
    observer->psi_r_hat_wb += cfg.control_period_s * dpsi_r;
    observer->last_config = cfg;

    if (!isfinite(observer->isd_hat_a) || !isfinite(observer->isq_hat_a) || !isfinite(observer->psi_r_hat_wb)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }

    sled_fill_output(observer, &d, omega_m_e, omega_s, b, gamma, k1, k2, output);
    return FLUX_OBSERVER_OK;
}
