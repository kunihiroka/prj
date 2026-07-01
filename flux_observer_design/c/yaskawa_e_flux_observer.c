#include "yaskawa_e_flux_observer.h"

#include <math.h>
#include <stddef.h>
#include <string.h>

typedef struct YaskawaEDerived {
    float sigma;
    float epsilon;
    float a11;
    float a12_r;
    float a12_j;
    float a21;
    float a22_r;
    float a22_j;
    float b1;
    float g1;
    float g2;
    float g3;
    float g4;
    float k1_sched;
    float k2_sched;
} YaskawaEDerived;

static int ye_is_finite_positive(float x)
{
    return isfinite(x) && (x > 0.0f);
}

static float ye_absf(float x)
{
    return (x >= 0.0f) ? x : -x;
}

static FluxObserverStatus ye_fetch_config(YaskawaEFluxObserver *observer, FluxObserverMotorConfig *cfg)
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
    if (!ye_is_finite_positive(cfg->rs_ohm) ||
        !ye_is_finite_positive(cfg->rr_ohm) ||
        !ye_is_finite_positive(cfg->lls_h) ||
        !ye_is_finite_positive(cfg->llr_h) ||
        !ye_is_finite_positive(cfg->lm_h) ||
        (cfg->pole_pairs == 0u) ||
        !ye_is_finite_positive(cfg->control_period_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static float ye_ramp_down_gain(float abs_speed, float w1, float w2)
{
    if (abs_speed <= w1) {
        return 1.0f;
    }
    if (abs_speed >= w2) {
        return 0.0f;
    }
    if (w2 <= w1) {
        return 0.0f;
    }
    return (w2 - abs_speed) / (w2 - w1);
}

static FluxObserverStatus ye_derive(
    const FluxObserverMotorConfig *cfg,
    const YaskawaEFluxObserver *observer,
    float omega_m_e,
    YaskawaEDerived *d)
{
    const float two_pi = 6.2831853071795864769f;
    float ls_h;
    float lr_h;
    float det_l;
    float coef;
    float t2;
    float g1_base;
    float rated_omega_m_e;
    float wx1;
    float wx2;
    float wx3;
    float wx4;
    float abs_omega;

    if ((cfg == NULL) || (observer == NULL) || (d == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    ls_h = cfg->lls_h + cfg->lm_h;
    lr_h = cfg->llr_h + cfg->lm_h;
    det_l = ls_h * lr_h - cfg->lm_h * cfg->lm_h;
    if (!isfinite(det_l) || (fabsf(det_l) < 1.0e-12f)) {
        return FLUX_OBSERVER_ERR_SINGULAR;
    }
    d->sigma = 1.0f - cfg->lm_h * cfg->lm_h / (ls_h * lr_h);
    d->epsilon = d->sigma * ls_h * lr_h / cfg->lm_h;
    if (!ye_is_finite_positive(d->sigma) || !ye_is_finite_positive(d->epsilon)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    coef = cfg->lm_h / (d->sigma * ls_h * lr_h);
    d->a11 = -(cfg->rs_ohm + cfg->rr_ohm * cfg->lm_h * cfg->lm_h / (lr_h * lr_h)) / (d->sigma * ls_h);
    d->a12_r = coef * cfg->rr_ohm / lr_h;
    d->a12_j = -coef * omega_m_e;
    d->a21 = cfg->lm_h * cfg->rr_ohm / lr_h;
    d->a22_r = -d->epsilon * d->a12_r;
    d->a22_j = -d->epsilon * d->a12_j;
    d->b1 = 1.0f / (d->sigma * ls_h);

    rated_omega_m_e = (float)cfg->pole_pairs * observer->rated_speed_rpm * two_pi / 60.0f;
    wx1 = ye_absf(rated_omega_m_e) * observer->wx1_ratio;
    wx2 = ye_absf(rated_omega_m_e) * observer->wx2_ratio;
    wx3 = ye_absf(rated_omega_m_e) * observer->wx3_ratio;
    wx4 = ye_absf(rated_omega_m_e) * observer->wx4_ratio;
    abs_omega = ye_absf(omega_m_e);
    d->k1_sched = ye_ramp_down_gain(abs_omega, wx1, wx2);
    d->k2_sched = ye_ramp_down_gain(abs_omega, wx3, wx4);

    t2 = lr_h / cfg->rr_ohm;
    if (!ye_is_finite_positive(t2)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    g1_base = -2.0f / t2 + cfg->rs_ohm / (d->sigma * ls_h) + cfg->rr_ohm / (d->sigma * lr_h);
    d->g1 = g1_base * d->k1_sched;
    d->g2 = 0.0f;
    d->g3 = (-d->epsilon * d->g1 - d->epsilon * cfg->rr_ohm / lr_h + cfg->rs_ohm * lr_h / cfg->lm_h) *
            d->k2_sched;
    d->g4 = (-d->epsilon * omega_m_e) * d->k2_sched;

    if (!isfinite(d->a11) || !isfinite(d->a12_r) || !isfinite(d->a12_j) ||
        !isfinite(d->a21) || !isfinite(d->a22_r) || !isfinite(d->a22_j) ||
        !isfinite(d->b1) || !isfinite(d->g1) || !isfinite(d->g2) ||
        !isfinite(d->g3) || !isfinite(d->g4) ||
        !isfinite(d->k1_sched) || !isfinite(d->k2_sched)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    return FLUX_OBSERVER_OK;
}

static void ye_derivative(
    const YaskawaEDerived *d,
    const YaskawaEFluxObserverInput *input,
    const YaskawaEFluxObserver *observer,
    float *dis_alpha,
    float *dis_beta,
    float *dphi_alpha,
    float *dphi_beta)
{
    float e_alpha;
    float e_beta;
    float j_e_alpha;
    float j_e_beta;
    float j_phi_alpha;
    float j_phi_beta;

    e_alpha = observer->is_alpha_hat_a - input->is_alpha_a;
    e_beta = observer->is_beta_hat_a - input->is_beta_a;
    j_e_alpha = -e_beta;
    j_e_beta = e_alpha;
    j_phi_alpha = -observer->phi_r_beta_hat_wb;
    j_phi_beta = observer->phi_r_alpha_hat_wb;

    *dis_alpha =
        d->a11 * observer->is_alpha_hat_a +
        d->a12_r * observer->phi_r_alpha_hat_wb +
        d->a12_j * j_phi_alpha +
        d->b1 * input->us_alpha_v +
        d->g1 * e_alpha +
        d->g2 * j_e_alpha;
    *dis_beta =
        d->a11 * observer->is_beta_hat_a +
        d->a12_r * observer->phi_r_beta_hat_wb +
        d->a12_j * j_phi_beta +
        d->b1 * input->us_beta_v +
        d->g1 * e_beta +
        d->g2 * j_e_beta;
    *dphi_alpha =
        d->a21 * observer->is_alpha_hat_a +
        d->a22_r * observer->phi_r_alpha_hat_wb +
        d->a22_j * j_phi_alpha +
        d->g3 * e_alpha +
        d->g4 * j_e_alpha;
    *dphi_beta =
        d->a21 * observer->is_beta_hat_a +
        d->a22_r * observer->phi_r_beta_hat_wb +
        d->a22_j * j_phi_beta +
        d->g3 * e_beta +
        d->g4 * j_e_beta;
}

static void ye_fill_output(
    const YaskawaEFluxObserver *observer,
    const YaskawaEDerived *d,
    float omega_m_e,
    YaskawaEFluxObserverOutput *output)
{
    if (output == NULL) {
        return;
    }
    output->is_alpha_hat_a = observer->is_alpha_hat_a;
    output->is_beta_hat_a = observer->is_beta_hat_a;
    output->phi_r_alpha_hat_wb = observer->phi_r_alpha_hat_wb;
    output->phi_r_beta_hat_wb = observer->phi_r_beta_hat_wb;
    output->omega_m_e_rad_s = omega_m_e;
    output->g1 = d->g1;
    output->g2 = d->g2;
    output->g3 = d->g3;
    output->g4 = d->g4;
    output->k1_sched = d->k1_sched;
    output->k2_sched = d->k2_sched;
    output->sigma = d->sigma;
    output->epsilon = d->epsilon;
}

void YaskawaEFluxObserver_Init(YaskawaEFluxObserver *observer, FluxObserverApi api)
{
    if (observer == NULL) {
        return;
    }
    memset(observer, 0, sizeof(*observer));
    observer->api = api;
    observer->rated_speed_rpm = 5000.0f;
    observer->wx1_ratio = 0.10f;
    observer->wx2_ratio = 0.15f;
    observer->wx3_ratio = 0.30f;
    observer->wx4_ratio = 0.50f;
}

void YaskawaEFluxObserver_SetRatedSpeed(YaskawaEFluxObserver *observer, float rated_speed_rpm)
{
    if (observer == NULL) {
        return;
    }
    if (!ye_is_finite_positive(rated_speed_rpm)) {
        return;
    }
    observer->rated_speed_rpm = rated_speed_rpm;
}

void YaskawaEFluxObserver_SetSchedule(
    YaskawaEFluxObserver *observer,
    float wx1_ratio,
    float wx2_ratio,
    float wx3_ratio,
    float wx4_ratio)
{
    if (observer == NULL) {
        return;
    }
    if (!isfinite(wx1_ratio) || !isfinite(wx2_ratio) || !isfinite(wx3_ratio) || !isfinite(wx4_ratio)) {
        return;
    }
    if ((wx1_ratio < 0.0f) || (wx1_ratio > wx2_ratio) || (wx2_ratio > wx3_ratio) ||
        (wx3_ratio > wx4_ratio) || (wx4_ratio <= 0.0f)) {
        return;
    }
    observer->wx1_ratio = wx1_ratio;
    observer->wx2_ratio = wx2_ratio;
    observer->wx3_ratio = wx3_ratio;
    observer->wx4_ratio = wx4_ratio;
}

FluxObserverStatus YaskawaEFluxObserver_Reset(
    YaskawaEFluxObserver *observer,
    float is_alpha_hat_a,
    float is_beta_hat_a,
    float phi_r_alpha_hat_wb,
    float phi_r_beta_hat_wb)
{
    if (observer == NULL) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!isfinite(is_alpha_hat_a) || !isfinite(is_beta_hat_a) ||
        !isfinite(phi_r_alpha_hat_wb) || !isfinite(phi_r_beta_hat_wb)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }
    observer->is_alpha_hat_a = is_alpha_hat_a;
    observer->is_beta_hat_a = is_beta_hat_a;
    observer->phi_r_alpha_hat_wb = phi_r_alpha_hat_wb;
    observer->phi_r_beta_hat_wb = phi_r_beta_hat_wb;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus YaskawaEFluxObserver_ResetFromCurrents(
    YaskawaEFluxObserver *observer,
    float is_alpha_a,
    float is_beta_a)
{
    FluxObserverMotorConfig cfg;
    FluxObserverStatus status;

    status = ye_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    observer->last_config = cfg;
    observer->is_alpha_hat_a = is_alpha_a;
    observer->is_beta_hat_a = is_beta_a;
    observer->phi_r_alpha_hat_wb = cfg.lm_h * is_alpha_a;
    observer->phi_r_beta_hat_wb = cfg.lm_h * is_beta_a;
    return FLUX_OBSERVER_OK;
}

FluxObserverStatus YaskawaEFluxObserver_Step(
    YaskawaEFluxObserver *observer,
    const YaskawaEFluxObserverInput *input,
    YaskawaEFluxObserverOutput *output)
{
    FluxObserverMotorConfig cfg;
    YaskawaEDerived d;
    float omega_m_e;
    float dis_alpha;
    float dis_beta;
    float dphi_alpha;
    float dphi_beta;
    FluxObserverStatus status;

    if ((observer == NULL) || (input == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    if (!isfinite(input->us_alpha_v) || !isfinite(input->us_beta_v) ||
        !isfinite(input->is_alpha_a) || !isfinite(input->is_beta_a) ||
        !isfinite(input->omega_m_rad_s)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    status = ye_fetch_config(observer, &cfg);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    omega_m_e = (float)cfg.pole_pairs * input->omega_m_rad_s;
    status = ye_derive(&cfg, observer, omega_m_e, &d);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    ye_derivative(&d, input, observer, &dis_alpha, &dis_beta, &dphi_alpha, &dphi_beta);
    observer->is_alpha_hat_a += cfg.control_period_s * dis_alpha;
    observer->is_beta_hat_a += cfg.control_period_s * dis_beta;
    observer->phi_r_alpha_hat_wb += cfg.control_period_s * dphi_alpha;
    observer->phi_r_beta_hat_wb += cfg.control_period_s * dphi_beta;
    observer->last_config = cfg;

    if (!isfinite(observer->is_alpha_hat_a) || !isfinite(observer->is_beta_hat_a) ||
        !isfinite(observer->phi_r_alpha_hat_wb) || !isfinite(observer->phi_r_beta_hat_wb)) {
        return FLUX_OBSERVER_ERR_PARAM;
    }

    ye_fill_output(observer, &d, omega_m_e, output);
    return FLUX_OBSERVER_OK;
}
