#ifndef FLUX_OBSERVER_H
#define FLUX_OBSERVER_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define FLUX_OBSERVER_H_ROWS (4u)
#define FLUX_OBSERVER_H_COLS (2u)
#define FLUX_OBSERVER_POLE_COUNT (4u)

typedef enum FluxObserverStatus {
    FLUX_OBSERVER_OK = 0,
    FLUX_OBSERVER_ERR_NULL = -1,
    FLUX_OBSERVER_ERR_API = -2,
    FLUX_OBSERVER_ERR_PARAM = -3,
    FLUX_OBSERVER_ERR_SINGULAR = -4
} FluxObserverStatus;

typedef enum FluxObserverGainDesign {
    FLUX_OBSERVER_GAIN_FOUR_POLE = 0,
    FLUX_OBSERVER_GAIN_HORI_5_3 = 1
} FluxObserverGainDesign;

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
    float last_isd_a;
    float last_isq_a;
    uint8_t has_last_current;
    float observer_poles_rad_s[FLUX_OBSERVER_POLE_COUNT];
    FluxObserverGainDesign gain_design;
    float hori_alpha_rad_s;
    float hori_beta_rad_s;
    FluxObserverMotorConfig last_config;
    float H[FLUX_OBSERVER_H_ROWS][FLUX_OBSERVER_H_COLS];
} FluxObserver;

void FluxObserver_Init(FluxObserver *observer, FluxObserverApi api);
void FluxObserver_SetPolePlacement(FluxObserver *observer, float bandwidth_rad_s, float fastest_ratio);
void FluxObserver_SetObserverPoles(
    FluxObserver *observer,
    float pole0_rad_s,
    float pole1_rad_s,
    float pole2_rad_s,
    float pole3_rad_s);
/* Uses Hori 1986 section 5.3 k1/k2 alpha-beta pole placement. This mode
 * observes rotor flux with the paper's model-correction observer; stator
 * flux is reconstructed from measured stator current and estimated rotor flux. */
void FluxObserver_SetHori53PolePlacement(
    FluxObserver *observer,
    float alpha_rad_s,
    float beta_rad_s);
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
