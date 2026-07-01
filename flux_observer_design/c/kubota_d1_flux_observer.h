#ifndef KUBOTA_D1_FLUX_OBSERVER_H
#define KUBOTA_D1_FLUX_OBSERVER_H

#include "flux_observer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct KubotaD1FluxObserverInput {
    float usd_v;
    float usq_v;
    float isd_a;
    float isq_a;
    float omega_m_rad_s;
} KubotaD1FluxObserverInput;

typedef struct KubotaD1FluxObserverOutput {
    float isd_hat_a;
    float isq_hat_a;
    float phi_r_hat_wb;
    float psi_sd_hat_wb;
    float psi_sq_hat_wb;
    float omega_m_e_rad_s;
    float omega_k_rad_s;
    float omega_slip_rad_s;
    float k_gain;
    float g1;
    float g2;
    float g3;
    float g4;
} KubotaD1FluxObserverOutput;

typedef struct KubotaD1FluxObserver {
    FluxObserverApi api;
    float isd_hat_a;
    float isq_hat_a;
    float phi_r_hat_wb;
    float k_gain;
    float min_denominator_wb;
    FluxObserverMotorConfig last_config;
} KubotaD1FluxObserver;

void KubotaD1FluxObserver_Init(KubotaD1FluxObserver *observer, FluxObserverApi api);
void KubotaD1FluxObserver_SetK(KubotaD1FluxObserver *observer, float k_gain);
void KubotaD1FluxObserver_SetMinDenominator(KubotaD1FluxObserver *observer, float min_denominator_wb);
FluxObserverStatus KubotaD1FluxObserver_Reset(
    KubotaD1FluxObserver *observer,
    float isd_hat_a,
    float isq_hat_a,
    float phi_r_hat_wb);
FluxObserverStatus KubotaD1FluxObserver_ResetFromCurrents(
    KubotaD1FluxObserver *observer,
    float isd_a,
    float isq_a);
FluxObserverStatus KubotaD1FluxObserver_Step(
    KubotaD1FluxObserver *observer,
    const KubotaD1FluxObserverInput *input,
    KubotaD1FluxObserverOutput *output);

#ifdef __cplusplus
}
#endif

#endif
