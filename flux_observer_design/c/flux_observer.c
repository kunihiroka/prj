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

static int fo_is_valid_observer_pole(float x)
{
    return isfinite(x) && (x < 0.0f);
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

static FluxObserverStatus fo_solve16(float aug[16][17], float x[16])
{
    unsigned int col;
    unsigned int row;
    unsigned int pivot_row;

    for (col = 0u; col < 16u; ++col) {
        float pivot_abs = fabsf(aug[col][col]);
        pivot_row = col;
        for (row = col + 1u; row < 16u; ++row) {
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
            for (k = col; k < 17u; ++k) {
                float tmp = aug[col][k];
                aug[col][k] = aug[pivot_row][k];
                aug[pivot_row][k] = tmp;
            }
        }
        {
            float pivot = aug[col][col];
            unsigned int k;
            for (k = col; k < 17u; ++k) {
                aug[col][k] /= pivot;
            }
        }
        for (row = 0u; row < 16u; ++row) {
            if (row != col) {
                float factor = aug[row][col];
                unsigned int k;
                for (k = col; k < 17u; ++k) {
                    aug[row][k] -= factor * aug[col][k];
                }
            }
        }
    }

    for (row = 0u; row < 16u; ++row) {
        x[row] = aug[row][16];
        if (!isfinite(x[row])) {
            return FLUX_OBSERVER_ERR_SINGULAR;
        }
    }
    return FLUX_OBSERVER_OK;
}

static void fo_build_state_matrices(
    const FluxObserverMotorConfig *cfg,
    const FoDerived *d,
    float omega_r_rad_s,
    float omega_k_rad_s,
    float A[4][4],
    float C[2][4])
{
    float omega_slip = omega_k_rad_s - omega_r_rad_s;

    memset(A, 0, 16u * sizeof(float));
    A[0][0] = -cfg->rs_ohm * d->cs0;
    A[0][1] = omega_k_rad_s;
    A[0][2] = -cfg->rs_ohm * d->cs1;
    A[1][0] = -omega_k_rad_s;
    A[1][1] = -cfg->rs_ohm * d->cs0;
    A[1][3] = -cfg->rs_ohm * d->cs1;
    A[2][0] = -cfg->rr_ohm * d->cr0;
    A[2][2] = -cfg->rr_ohm * d->cr1;
    A[2][3] = omega_slip;
    A[3][1] = -cfg->rr_ohm * d->cr0;
    A[3][2] = -omega_slip;
    A[3][3] = -cfg->rr_ohm * d->cr1;

    memset(C, 0, 8u * sizeof(float));
    C[0][0] = d->cs0;
    C[0][2] = d->cs1;
    C[1][1] = d->cs0;
    C[1][3] = d->cs1;
}

static FluxObserverStatus fo_build_target_matrix(const FluxObserver *observer, float F[4][4])
{
    unsigned int i;

    if ((observer == NULL) || (F == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    memset(F, 0, 16u * sizeof(float));
    if (observer->gain_design == FLUX_OBSERVER_GAIN_FOUR_POLE) {
        for (i = 0u; i < FLUX_OBSERVER_POLE_COUNT; ++i) {
            if (!fo_is_valid_observer_pole(observer->observer_poles_rad_s[i])) {
                return FLUX_OBSERVER_ERR_PARAM;
            }
            F[i][i] = observer->observer_poles_rad_s[i];
        }
        return FLUX_OBSERVER_OK;
    }

    if (observer->gain_design == FLUX_OBSERVER_GAIN_HORI_5_3) {
        float alpha = observer->hori_alpha_rad_s;
        float beta = observer->hori_beta_rad_s;
        if (!fo_is_finite_positive(alpha) || !fo_is_finite_positive(beta)) {
            return FLUX_OBSERVER_ERR_PARAM;
        }

        F[0][0] = -alpha;
        F[0][1] = -beta;
        F[1][0] = beta;
        F[1][1] = -alpha;
        F[2][2] = -alpha;
        F[2][3] = -beta;
        F[3][2] = beta;
        F[3][3] = -alpha;
        return FLUX_OBSERVER_OK;
    }

    return FLUX_OBSERVER_ERR_PARAM;
}

static FluxObserverStatus fo_build_distribution_matrix(const FluxObserver *observer, float G[4][2])
{
    if ((observer == NULL) || (G == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }

    memset(G, 0, 8u * sizeof(float));
    if (observer->gain_design == FLUX_OBSERVER_GAIN_FOUR_POLE) {
        G[0][0] = 1.0f;
        G[1][1] = 1.0f;
        G[2][0] = 1.0f;
        G[3][1] = 1.0f;
        return FLUX_OBSERVER_OK;
    }

    if (observer->gain_design == FLUX_OBSERVER_GAIN_HORI_5_3) {
        G[0][0] = 1.0f;
        G[1][1] = 1.0f;
        G[2][1] = 1.0f;
        G[3][0] = 1.0f;
        return FLUX_OBSERVER_OK;
    }

    return FLUX_OBSERVER_ERR_PARAM;
}

static FluxObserverStatus fo_observer_H(
    const FluxObserver *observer,
    const FluxObserverMotorConfig *cfg,
    const FoDerived *d,
    float omega_r_rad_s,
    float omega_k_rad_s,
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS])
{
    float A[4][4];
    float C[2][4];
    float F[4][4];
    float G[4][2];
    float rhs[4][4];
    float aug16[16][17];
    float t_vec[16];
    float T[4][4];
    unsigned int i;
    unsigned int j;
    unsigned int k;
    FluxObserverStatus status;

    if ((observer == NULL) || (cfg == NULL) || (d == NULL) || (H == NULL)) {
        return FLUX_OBSERVER_ERR_NULL;
    }
    status = fo_build_target_matrix(observer, F);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }
    status = fo_build_distribution_matrix(observer, G);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    fo_build_state_matrices(cfg, d, omega_r_rad_s, omega_k_rad_s, A, C);

    for (i = 0u; i < 4u; ++i) {
        for (j = 0u; j < 4u; ++j) {
            rhs[i][j] = G[i][0] * C[0][j] + G[i][1] * C[1][j];
        }
    }

    memset(aug16, 0, sizeof(aug16));
    for (i = 0u; i < 4u; ++i) {
        for (j = 0u; j < 4u; ++j) {
            unsigned int eq = i + 4u * j;
            for (k = 0u; k < 4u; ++k) {
                unsigned int unknown = i + 4u * k;
                aug16[eq][unknown] += A[k][j];
            }
            for (k = 0u; k < 4u; ++k) {
                unsigned int unknown = k + 4u * j;
                aug16[eq][unknown] -= F[i][k];
            }
            aug16[eq][16] = rhs[i][j];
        }
    }

    status = fo_solve16(aug16, t_vec);
    if (status != FLUX_OBSERVER_OK) {
        return status;
    }

    for (i = 0u; i < 4u; ++i) {
        for (j = 0u; j < 4u; ++j) {
            T[i][j] = t_vec[i + 4u * j];
        }
    }

    for (j = 0u; j < 2u; ++j) {
        float aug4[4][5];
        float h_col[4];
        for (i = 0u; i < 4u; ++i) {
            for (k = 0u; k < 4u; ++k) {
                aug4[i][k] = T[i][k];
            }
            aug4[i][4] = G[i][j];
        }
        status = fo_solve4(aug4, h_col);
        if (status != FLUX_OBSERVER_OK) {
            return status;
        }
        for (i = 0u; i < 4u; ++i) {
            H[i][j] = h_col[i];
        }
    }

    return FLUX_OBSERVER_OK;
}

void FluxObserver_Init(FluxObserver *observer, FluxObserverApi api)
{
    if (observer == NULL) {
        return;
    }
    memset(observer, 0, sizeof(*observer));
    observer->api = api;
    FluxObserver_SetPolePlacement(observer, 2200.0f, 2.0f);
}

void FluxObserver_SetPolePlacement(FluxObserver *observer, float bandwidth_rad_s, float fastest_ratio)
{
    if (observer == NULL) {
        return;
    }
    if (!fo_is_finite_positive(bandwidth_rad_s)) {
        return;
    }
    if (!fo_is_finite_positive(fastest_ratio) || (fastest_ratio <= 1.55f)) {
        fastest_ratio = 2.0f;
    }
    observer->gain_design = FLUX_OBSERVER_GAIN_FOUR_POLE;
    observer->observer_poles_rad_s[0] = -bandwidth_rad_s;
    observer->observer_poles_rad_s[1] = -1.25f * bandwidth_rad_s;
    observer->observer_poles_rad_s[2] = -1.55f * bandwidth_rad_s;
    observer->observer_poles_rad_s[3] = -fastest_ratio * bandwidth_rad_s;
}

void FluxObserver_SetObserverPoles(
    FluxObserver *observer,
    float pole0_rad_s,
    float pole1_rad_s,
    float pole2_rad_s,
    float pole3_rad_s)
{
    if (observer == NULL) {
        return;
    }
    if (!fo_is_valid_observer_pole(pole0_rad_s) ||
        !fo_is_valid_observer_pole(pole1_rad_s) ||
        !fo_is_valid_observer_pole(pole2_rad_s) ||
        !fo_is_valid_observer_pole(pole3_rad_s)) {
        return;
    }
    observer->observer_poles_rad_s[0] = pole0_rad_s;
    observer->observer_poles_rad_s[1] = pole1_rad_s;
    observer->observer_poles_rad_s[2] = pole2_rad_s;
    observer->observer_poles_rad_s[3] = pole3_rad_s;
    observer->gain_design = FLUX_OBSERVER_GAIN_FOUR_POLE;
}

void FluxObserver_SetHori53PolePlacement(
    FluxObserver *observer,
    float alpha_rad_s,
    float beta_rad_s)
{
    if (observer == NULL) {
        return;
    }
    if (!fo_is_finite_positive(alpha_rad_s) || !fo_is_finite_positive(beta_rad_s)) {
        return;
    }
    observer->hori_alpha_rad_s = alpha_rad_s;
    observer->hori_beta_rad_s = beta_rad_s;
    observer->gain_design = FLUX_OBSERVER_GAIN_HORI_5_3;
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
