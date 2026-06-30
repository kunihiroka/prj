#ifndef SLED23_FLUX_OBSERVER_H
#define SLED23_FLUX_OBSERVER_H

#include "flux_observer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Sled23FluxObserverInput {
    float usd_v;
    float usq_v;
    float isd_a;
    float isq_a;
    float omega_m_rad_s;
} Sled23FluxObserverInput;

typedef struct Sled23FluxObserverOutput {
    float isd_hat_a;
    float isq_hat_a;
    float psi_r_hat_wb;
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
} Sled23FluxObserverOutput;

typedef struct Sled23FluxObserver {
    FluxObserverApi api;
    float isd_hat_a;
    float isq_hat_a;
    float psi_r_hat_wb;
    float alpha_i_rad_s;
    float zeta_inf;
    float fixed_b_rad_s;
    uint8_t use_fixed_b;
    float min_denominator_wb;
    FluxObserverMotorConfig last_config;
} Sled23FluxObserver;

void Sled23FluxObserver_Init(Sled23FluxObserver *observer, FluxObserverApi api);
void Sled23FluxObserver_SetDesign(
    Sled23FluxObserver *observer,
    float alpha_i_rad_s,
    float zeta_inf);
void Sled23FluxObserver_SetFixedB(
    Sled23FluxObserver *observer,
    float alpha_i_rad_s,
    float b_rad_s);
void Sled23FluxObserver_SetMinDenominator(Sled23FluxObserver *observer, float min_denominator_wb);
FluxObserverStatus Sled23FluxObserver_Reset(
    Sled23FluxObserver *observer,
    float isd_hat_a,
    float isq_hat_a,
    float psi_r_hat_wb);
FluxObserverStatus Sled23FluxObserver_ResetFromCurrents(
    Sled23FluxObserver *observer,
    float isd_a,
    float isq_a);
FluxObserverStatus Sled23FluxObserver_Step(
    Sled23FluxObserver *observer,
    const Sled23FluxObserverInput *input,
    Sled23FluxObserverOutput *output);

#ifdef __cplusplus
}
#endif

#endif
