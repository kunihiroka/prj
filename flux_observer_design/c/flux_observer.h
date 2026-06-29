#ifndef FLUX_OBSERVER_H
#define FLUX_OBSERVER_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FLUX_OBSERVER_H_ROWS (4u)
#define FLUX_OBSERVER_H_COLS (2u)

typedef enum FluxObserverStatus {
    FLUX_OBSERVER_OK = 0,
    FLUX_OBSERVER_ERR_NULL = -1,
    FLUX_OBSERVER_ERR_API = -2,
    FLUX_OBSERVER_ERR_PARAM = -3,
    FLUX_OBSERVER_ERR_SINGULAR = -4
} FluxObserverStatus;

typedef struct FluxObserverMotorConfig {
    float rs_ohm;
    float rr_ohm;
    float lls_h;
    float llr_h;
    float lm_h;
    uint16_t pole_pairs;
    float control_period_s;
} FluxObserverMotorConfig;

typedef int (*FluxObserverGetMotorConfigFn)(void *user, FluxObserverMotorConfig *config);

typedef struct FluxObserverApi {
    void *user;
    FluxObserverGetMotorConfigFn get_motor_config;
} FluxObserverApi;

typedef struct FluxObserverInput {
    float vsd_v;
    float vsq_v;
    float isd_a;
    float isq_a;
    float omega_m_rad_s;
    float omega_slip_rad_s;
} FluxObserverInput;

typedef struct FluxObserverOutput {
    float psi_sd_wb;
    float psi_sq_wb;
    float psi_rd_wb;
    float psi_rq_wb;
    float isd_hat_a;
    float isq_hat_a;
    float ird_hat_a;
    float irq_hat_a;
    float omega_r_rad_s;
    float omega_k_rad_s;
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS];
} FluxObserverOutput;

typedef struct FluxObserver {
    FluxObserverApi api;
    float psi_sd_wb;
    float psi_sq_wb;
    float psi_rd_wb;
    float psi_rq_wb;
    float observer_bandwidth_rad_s;
    float observer_pole_ratio;
    FluxObserverMotorConfig last_config;
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS];
} FluxObserver;

void FluxObserver_Init(FluxObserver *observer, FluxObserverApi api);
void FluxObserver_SetPolePlacement(FluxObserver *observer, float bandwidth_rad_s, float pole_ratio);
FluxObserverStatus FluxObserver_ResetFlux(
    FluxObserver *observer,
    float psi_sd_wb,
    float psi_sq_wb,
    float psi_rd_wb,
    float psi_rq_wb);
FluxObserverStatus FluxObserver_ResetFromCurrents(
    FluxObserver *observer,
    float isd_a,
    float isq_a,
    float ird_a,
    float irq_a);
FluxObserverStatus FluxObserver_Step(
    FluxObserver *observer,
    const FluxObserverInput *input,
    FluxObserverOutput *output);
FluxObserverStatus FluxObserver_GetLastH(
    const FluxObserver *observer,
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS]);
const char *FluxObserver_StatusString(FluxObserverStatus status);

#ifdef __cplusplus
}
#endif

#endif
