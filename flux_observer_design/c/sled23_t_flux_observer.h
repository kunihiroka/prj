#ifndef SLED23_T_FLUX_OBSERVER_H
#define SLED23_T_FLUX_OBSERVER_H

#include "flux_observer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Sled23TFluxObserverInput {
    float usd_v;
    float usq_v;
    float isd_a;
    float isq_a;
    float omega_m_rad_s;
} Sled23TFluxObserverInput;

typedef struct Sled23TFluxObserverOutput {
    float isd_hat_a;
    float isq_hat_a;
    float phi_r_hat_wb;
    float psi_sd_hat_wb;
    float psi_sq_hat_wb;
    float omega_m_e_rad_s;
    float omega_s_rad_s;
    float omega_slip_rad_s;
    float alpha_rad_s;
    float alpha_i_rad_s;
    float b_rad_s;
    float gamma_rad_s;
    float k1;
    float k2;
} Sled23TFluxObserverOutput;

typedef struct Sled23TFluxObserver {
    FluxObserverApi api;
    float isd_hat_a;
    float isq_hat_a;
    float phi_r_hat_wb;
    float alpha_i_rad_s;
    float zeta_inf;
    float fixed_b_rad_s;
    uint8_t use_fixed_b;
    float min_denominator_wb;
    FluxObserverMotorConfig last_config;
} Sled23TFluxObserver;

void Sled23TFluxObserver_Init(Sled23TFluxObserver *observer, FluxObserverApi api);
void Sled23TFluxObserver_SetDesign(
    Sled23TFluxObserver *observer,
    float alpha_i_rad_s,
    float zeta_inf);
void Sled23TFluxObserver_SetFixedB(
    Sled23TFluxObserver *observer,
    float alpha_i_rad_s,
    float b_rad_s);
void Sled23TFluxObserver_SetMinDenominator(Sled23TFluxObserver *observer, float min_denominator_wb);
FluxObserverStatus Sled23TFluxObserver_Reset(
    Sled23TFluxObserver *observer,
    float isd_hat_a,
    float isq_hat_a,
    float phi_r_hat_wb);
FluxObserverStatus Sled23TFluxObserver_ResetFromCurrents(
    Sled23TFluxObserver *observer,
    float isd_a,
    float isq_a);
FluxObserverStatus Sled23TFluxObserver_Step(
    Sled23TFluxObserver *observer,
    const Sled23TFluxObserverInput *input,
    Sled23TFluxObserverOutput *output);

#ifdef __cplusplus
}
#endif

#endif
