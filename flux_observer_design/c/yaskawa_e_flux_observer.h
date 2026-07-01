#ifndef YASKAWA_E_FLUX_OBSERVER_H
#define YASKAWA_E_FLUX_OBSERVER_H

#include "flux_observer.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct YaskawaEFluxObserverInput {
    float us_alpha_v;
    float us_beta_v;
    float is_alpha_a;
    float is_beta_a;
    float omega_m_rad_s;
} YaskawaEFluxObserverInput;

typedef struct YaskawaEFluxObserverOutput {
    float is_alpha_hat_a;
    float is_beta_hat_a;
    float phi_r_alpha_hat_wb;
    float phi_r_beta_hat_wb;
    float omega_m_e_rad_s;
    float g1;
    float g2;
    float g3;
    float g4;
    float k1_sched;
    float k2_sched;
    float sigma;
    float epsilon;
} YaskawaEFluxObserverOutput;

typedef struct YaskawaEFluxObserver {
    FluxObserverApi api;
    float is_alpha_hat_a;
    float is_beta_hat_a;
    float phi_r_alpha_hat_wb;
    float phi_r_beta_hat_wb;
    float rated_speed_rpm;
    float wx1_ratio;
    float wx2_ratio;
    float wx3_ratio;
    float wx4_ratio;
    FluxObserverMotorConfig last_config;
} YaskawaEFluxObserver;

void YaskawaEFluxObserver_Init(YaskawaEFluxObserver *observer, FluxObserverApi api);
void YaskawaEFluxObserver_SetRatedSpeed(YaskawaEFluxObserver *observer, float rated_speed_rpm);
void YaskawaEFluxObserver_SetSchedule(
    YaskawaEFluxObserver *observer,
    float wx1_ratio,
    float wx2_ratio,
    float wx3_ratio,
    float wx4_ratio);
FluxObserverStatus YaskawaEFluxObserver_Reset(
    YaskawaEFluxObserver *observer,
    float is_alpha_hat_a,
    float is_beta_hat_a,
    float phi_r_alpha_hat_wb,
    float phi_r_beta_hat_wb);
FluxObserverStatus YaskawaEFluxObserver_ResetFromCurrents(
    YaskawaEFluxObserver *observer,
    float is_alpha_a,
    float is_beta_a);
FluxObserverStatus YaskawaEFluxObserver_Step(
    YaskawaEFluxObserver *observer,
    const YaskawaEFluxObserverInput *input,
    YaskawaEFluxObserverOutput *output);

#ifdef __cplusplus
}
#endif

#endif
