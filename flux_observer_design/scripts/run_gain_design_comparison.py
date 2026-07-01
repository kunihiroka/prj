"""Compare the implemented flux-observer gain-design methods.

Methods:

    A_four_pole_sylvester:
        Real four-state observer with four real poles and Sylvester equation.

    B_hori53_paper:
        Hori 1986 section 5.3 model-correction flux observer with paper
        k1/k2 alpha/beta pole placement.

    C_sled23_t_type:
        SLED 2023 Appendix A observer in estimated rotor-flux coordinates,
        written with the T-type physical rotor flux.

    D1_kubota_k1p2_t_type:
        Kubota/Matsuse k-times natural-pole observer, transformed to the
        estimated rotor-flux d-axis coordinates.

    E_yaskawa_gain_scheduled:
        Yaskawa/Takase 2023 adaptive rotor-flux observer gain scheduling
        based on Popov hyperstability and the Kalman-Yakubovich lemma.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

import run_flux_observer_evaluation as fo


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"

METHOD_A = "A_four_pole_sylvester"
METHOD_B = "B_hori53_paper"
METHOD_C = "C_sled23_t_type"
METHOD_D1 = "D1_kubota_k1p2_t_type"
METHOD_E = "E_yaskawa_gain_scheduled"
METHODS = [METHOD_D1, METHOD_C, METHOD_E, METHOD_A, METHOD_B]

SLED_ALPHA_I = 1000.0
SLED_ZETA_INF = 0.4
KUBOTA_D1_K = 1.2
YASKAWA_E_RATED_SPEED_RPM = 5000.0
YASKAWA_E_WX_RATIOS = (0.10, 0.15, 0.30, 0.50)


@dataclass
class MethodResult:
    method: str
    op: fo.OperatingPoint
    case: fo.ErrorCase
    t: np.ndarray
    current_error_norm: np.ndarray
    flux_error_norm: np.ndarray
    stable: bool
    metrics: dict[str, float]


@dataclass(frozen=True)
class TTypeDerived:
    l_sigma: float
    rho: float
    alpha: float
    c_r: float
    r_sigma: float


@dataclass(frozen=True)
class KubotaD1Derived:
    a_r11: float
    a_r12: float
    a_i12: float
    a_r21: float
    a_r22: float
    a_i22: float
    b1: float
    g1: float
    g2: float
    g3: float
    g4: float


@dataclass(frozen=True)
class YaskawaEDerived:
    sigma: float
    epsilon: float
    a11: np.ndarray
    a12: np.ndarray
    a21: np.ndarray
    a22: np.ndarray
    b1: np.ndarray
    g1: float
    g2: float
    g3: float
    g4: float
    k1_sched: float
    k2_sched: float


def convergence_time(t: np.ndarray, err: np.ndarray, threshold: float) -> float:
    ok = np.abs(err) <= threshold
    if not np.any(ok):
        return float("nan")
    suffix_ok = np.flip(np.cumprod(np.flip(ok).astype(int))).astype(bool)
    idx = np.argmax(suffix_ok)
    if suffix_ok[idx]:
        return float(t[idx])
    return float("nan")


def result_from_flux_evaluation(
    method: str,
    true_params: fo.MotorParams,
    op: fo.OperatingPoint,
    case: fo.ErrorCase,
) -> MethodResult:
    if method == METHOD_A:
        result = fo.simulate_case(true_params, op, case, gain_design=fo.GAIN_DESIGN_FOUR_POLE)
    elif method == METHOD_B:
        result = fo.simulate_case(
            true_params,
            op,
            case,
            gain_design=fo.GAIN_DESIGN_HORI_5_3,
            hori_alpha=2200.0,
            hori_beta=500.0,
        )
    else:
        raise ValueError(method)

    flux_err = result.est_flux[:, 2:4] - result.true_flux[:, 2:4]
    current_err = result.est_current[:, 0:2] - result.true_current[:, 0:2]
    return MethodResult(
        method=method,
        op=op,
        case=case,
        t=result.t,
        current_error_norm=np.linalg.norm(current_err, axis=1),
        flux_error_norm=np.linalg.norm(flux_err, axis=1),
        stable=result.stable,
        metrics={
            "psi_r_rms_wb": result.metrics["psi_r_rms_wb"],
            "psi_r_rel_pct": result.metrics["psi_r_rel_pct"],
            "is_rms_a": result.metrics["is_rms_a"],
            "is_rel_pct": result.metrics["is_rel_pct"],
            "psi_r_conv_time_ms": result.metrics["psi_r_conv_time_ms"],
            "is_conv_time_ms": result.metrics["is_conv_time_ms"],
        },
    )


def t_type_derived(params: fo.MotorParams) -> TTypeDerived:
    rho = params.lm / params.lr
    return TTypeDerived(
        l_sigma=params.det_l / params.lr,
        rho=rho,
        alpha=params.rr / params.lr,
        c_r=params.rr * params.lm / params.lr,
        r_sigma=params.rs + params.rr * rho * rho,
    )


def sled23_t_step_with_observer_params(
    obs_params: fo.MotorParams,
    x_hat: np.ndarray,
    i_meas: np.ndarray,
    u_meas: np.ndarray,
    omega_m_e: float,
    dt: float,
) -> tuple[np.ndarray, float]:
    """SLED Appendix A written with the T-type physical rotor flux phi_r.

    The inverse-gamma SLED paper state psi_R is related to the physical
    T-equivalent rotor flux by psi_R = (M / Lr) phi_r.  This function keeps
    phi_r as the state so that all reported flux errors are physical rotor
    flux errors and comparable with methods A/B.
    """
    d = t_type_derived(obs_params)
    isd_hat, isq_hat, phi_hat = x_hat
    eisd = i_meas[0] - isd_hat
    eisq = i_meas[1] - isq_hat

    b = 2.0 * SLED_ZETA_INF * abs(omega_m_e) + d.alpha
    gamma = SLED_ALPHA_I - d.alpha
    k_den = d.alpha * d.alpha + omega_m_e * omega_m_e
    k1 = b * d.alpha / k_den
    k2 = b * omega_m_e / k_den

    det_over_m = obs_params.det_l / obs_params.lm
    denominator = phi_hat - det_over_m * eisd
    if abs(denominator) < 1.0e-9:
        denominator = math.copysign(1.0e-9, denominator if denominator != 0.0 else 1.0)
    omega_s = omega_m_e + (
        d.c_r * i_meas[1] + k2 * SLED_ALPHA_I * det_over_m * eisd - gamma * det_over_m * eisq
    ) / denominator

    disd = (
        d.alpha * d.rho * phi_hat
        - d.r_sigma * isd_hat
        + omega_s * d.l_sigma * isq_hat
        + u_meas[0]
        + d.l_sigma * (gamma * eisd - omega_m_e * eisq)
    ) / d.l_sigma
    disq = (
        -omega_m_e * d.rho * phi_hat
        - d.r_sigma * isq_hat
        - omega_s * d.l_sigma * isd_hat
        + u_meas[1]
        + d.l_sigma * (gamma * eisq + omega_m_e * eisd)
    ) / d.l_sigma
    dphi = (
        -d.alpha * phi_hat
        + d.c_r * isd_hat
        + (k1 * SLED_ALPHA_I - gamma) * det_over_m * eisd
        - (omega_s - omega_m_e) * det_over_m * eisq
    )

    return x_hat + dt * np.array([disd, disq, dphi], dtype=float), omega_s


def kubota_d1_derived(params: fo.MotorParams, omega_r_e: float, k_gain: float) -> KubotaD1Derived:
    sigma = 1.0 - params.lm * params.lm / (params.ls * params.lr)
    tau_r = params.lr / params.rr
    a_r11 = -(params.rs / (sigma * params.ls) + params.rr * (1.0 - sigma) / (sigma * params.lr))
    a_r12 = params.lm / (sigma * params.ls * params.lr * tau_r)
    a_i12 = -params.lm * omega_r_e / (sigma * params.ls * params.lr)
    a_r21 = params.lm / tau_r
    a_r22 = -1.0 / tau_r
    a_i22 = omega_r_e
    alpha_rel = -sigma * params.ls * params.lr / params.lm
    g1 = (k_gain - 1.0) * (a_r11 + a_r22)
    g2 = (k_gain - 1.0) * a_i22
    g3 = (k_gain * k_gain - 1.0) * (-alpha_rel * a_r11 + a_r21)
    g3 += alpha_rel * (k_gain - 1.0) * (a_r11 + a_r22)
    g4 = alpha_rel * (k_gain - 1.0) * a_i22
    return KubotaD1Derived(
        a_r11=a_r11,
        a_r12=a_r12,
        a_i12=a_i12,
        a_r21=a_r21,
        a_r22=a_r22,
        a_i22=a_i22,
        b1=1.0 / (sigma * params.ls),
        g1=g1,
        g2=g2,
        g3=g3,
        g4=g4,
    )


def kubota_d1_step_with_observer_params(
    obs_params: fo.MotorParams,
    x_hat: np.ndarray,
    i_meas: np.ndarray,
    u_meas: np.ndarray,
    omega_r_e: float,
    dt: float,
    k_gain: float = KUBOTA_D1_K,
) -> tuple[np.ndarray, float]:
    d = kubota_d1_derived(obs_params, omega_r_e, k_gain)
    isd_hat, isq_hat, phi_hat = x_hat
    ed = isd_hat - i_meas[0]
    eq = isq_hat - i_meas[1]
    denominator = phi_hat
    if abs(denominator) < 1.0e-9:
        denominator = math.copysign(1.0e-9, denominator if denominator != 0.0 else 1.0)

    omega_k = d.a_i22 + (d.a_r21 * isq_hat + d.g4 * ed + d.g3 * eq) / denominator
    disd = (
        d.a_r11 * isd_hat
        + omega_k * isq_hat
        + d.a_r12 * phi_hat
        + d.b1 * u_meas[0]
        + d.g1 * ed
        - d.g2 * eq
    )
    disq = (
        d.a_r11 * isq_hat
        - omega_k * isd_hat
        + d.a_i12 * phi_hat
        + d.b1 * u_meas[1]
        + d.g2 * ed
        + d.g1 * eq
    )
    dphi = d.a_r21 * isd_hat + d.a_r22 * phi_hat + d.g3 * ed - d.g4 * eq
    return x_hat + dt * np.array([disd, disq, dphi], dtype=float), omega_k


def j_matrix_2() -> np.ndarray:
    return np.array([[0.0, -1.0], [1.0, 0.0]], dtype=float)


def rotate_dq_to_stationary(theta: float, v_dq: np.ndarray) -> np.ndarray:
    c = math.cos(theta)
    s = math.sin(theta)
    return np.array([c * v_dq[0] - s * v_dq[1], s * v_dq[0] + c * v_dq[1]], dtype=float)


def ramp_down_gain(abs_speed: float, w1: float, w2: float) -> float:
    if abs_speed <= w1:
        return 1.0
    if abs_speed >= w2:
        return 0.0
    if w2 <= w1:
        return 0.0
    return float((w2 - abs_speed) / (w2 - w1))


def yaskawa_e_derived(
    params: fo.MotorParams,
    omega_m_e: float,
    rated_speed_rpm: float = YASKAWA_E_RATED_SPEED_RPM,
    wx_ratios: tuple[float, float, float, float] = YASKAWA_E_WX_RATIOS,
) -> YaskawaEDerived:
    """Yaskawa/Takase 2023 Eq. (11) in stationary alpha-beta coordinates."""
    sigma = 1.0 - params.lm * params.lm / (params.ls * params.lr)
    epsilon = sigma * params.ls * params.lr / params.lm
    ident = np.eye(2, dtype=float)
    j = j_matrix_2()
    a11_scalar = -(params.rs + params.rr * params.lm * params.lm / (params.lr * params.lr)) / (sigma * params.ls)
    a11 = a11_scalar * ident
    a12 = params.lm / (sigma * params.ls * params.lr) * (params.rr / params.lr * ident - omega_m_e * j)
    a21 = params.lm * params.rr / params.lr * ident
    a22 = -epsilon * a12
    b1 = 1.0 / (sigma * params.ls) * ident

    rated_omega_m_e = params.pole_pairs * rated_speed_rpm * 2.0 * math.pi / 60.0
    wx1, wx2, wx3, wx4 = [abs(rated_omega_m_e) * r for r in wx_ratios]
    abs_omega = abs(omega_m_e)
    k1_sched = ramp_down_gain(abs_omega, wx1, wx2)
    k2_sched = ramp_down_gain(abs_omega, wx3, wx4)
    t2 = params.lr / params.rr
    g1_base = -2.0 / t2 + params.rs / (sigma * params.ls) + params.rr / (sigma * params.lr)
    g1 = g1_base * k1_sched
    g2 = 0.0
    g3 = (-epsilon * g1 - epsilon * params.rr / params.lr + params.rs * params.lr / params.lm) * k2_sched
    g4 = (-epsilon * omega_m_e) * k2_sched
    return YaskawaEDerived(
        sigma=sigma,
        epsilon=epsilon,
        a11=a11,
        a12=a12,
        a21=a21,
        a22=a22,
        b1=b1,
        g1=g1,
        g2=g2,
        g3=g3,
        g4=g4,
        k1_sched=k1_sched,
        k2_sched=k2_sched,
    )


def yaskawa_e_step_with_observer_params(
    obs_params: fo.MotorParams,
    x_hat: np.ndarray,
    i_meas: np.ndarray,
    u_meas: np.ndarray,
    omega_m_e: float,
    dt: float,
) -> np.ndarray:
    return x_hat + dt * yaskawa_e_derivative_with_observer_params(
        obs_params=obs_params,
        x_hat=x_hat,
        i_meas=i_meas,
        u_meas=u_meas,
        omega_m_e=omega_m_e,
    )


def yaskawa_e_derivative_with_observer_params(
    obs_params: fo.MotorParams,
    x_hat: np.ndarray,
    i_meas: np.ndarray,
    u_meas: np.ndarray,
    omega_m_e: float,
) -> np.ndarray:
    d = yaskawa_e_derived(obs_params, omega_m_e)
    i_hat = x_hat[0:2]
    phi_hat = x_hat[2:4]
    e = i_hat - i_meas
    g_upper = d.g1 * np.eye(2, dtype=float) + d.g2 * j_matrix_2()
    g_lower = d.g3 * np.eye(2, dtype=float) + d.g4 * j_matrix_2()
    di = d.a11 @ i_hat + d.a12 @ phi_hat + d.b1 @ u_meas + g_upper @ e
    dphi = d.a21 @ i_hat + d.a22 @ phi_hat + g_lower @ e
    return np.array([di[0], di[1], dphi[0], dphi[1]], dtype=float)


def result_from_sled23_t(
    true_params: fo.MotorParams,
    op: fo.OperatingPoint,
    case: fo.ErrorCase,
    dt: float = 1.0e-5,
    t_end: float = 0.12,
    seed: int = 1,
) -> MethodResult:
    obs_params = fo.apply_param_scale(true_params, case.param_scale)
    ss = fo.steady_operating_point(true_params, op)
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    i_true = np.array(ss["i_s"], dtype=float)
    u_true = np.array(ss["v_s"], dtype=float)
    phi_true = float(np.array(ss["psi_r"], dtype=float)[0])
    x_hat = np.array([0.6 * i_true[0], 1.4 * i_true[1], 0.5 * phi_true], dtype=float)

    current_error_norm = np.zeros(n, dtype=float)
    flux_error_norm = np.zeros(n, dtype=float)
    stable = True
    rng = np.random.default_rng(seed)
    v_offset = np.array([case.voltage_offset_d, case.voltage_offset_q], dtype=float)
    i_offset = np.array([case.current_offset_d, case.current_offset_q], dtype=float)

    for k in range(n):
        current_error_norm[k] = float(np.linalg.norm(x_hat[0:2] - i_true))
        flux_error_norm[k] = abs(x_hat[2] - phi_true)
        if k == n - 1:
            break
        noise = np.zeros(2, dtype=float)
        if case.current_noise_rms > 0.0:
            noise = rng.normal(0.0, case.current_noise_rms, size=2)
        u_meas = (1.0 + case.voltage_gain) * u_true + v_offset
        i_meas = (1.0 + case.current_gain) * i_true + i_offset + noise
        x_hat, _ = sled23_t_step_with_observer_params(
            obs_params=obs_params,
            x_hat=x_hat,
            i_meas=i_meas,
            u_meas=u_meas,
            omega_m_e=float(ss["omega_r"]),
            dt=dt,
        )
        if not np.all(np.isfinite(x_hat)) or np.max(np.abs(x_hat)) > 1.0e6:
            stable = False
            current_error_norm = current_error_norm[: k + 1]
            flux_error_norm = flux_error_norm[: k + 1]
            t = t[: k + 1]
            break

    start = max(0, int(0.8 * len(t)))
    is_mag = max(float(np.linalg.norm(i_true)), 1.0e-12)
    psi_mag = max(abs(phi_true), 1.0e-12)
    psi_r_rms = float(np.sqrt(np.mean(flux_error_norm[start:] ** 2)))
    is_rms = float(np.sqrt(np.mean(current_error_norm[start:] ** 2)))
    return MethodResult(
        method=METHOD_C,
        op=op,
        case=case,
        t=t,
        current_error_norm=current_error_norm,
        flux_error_norm=flux_error_norm,
        stable=stable,
        metrics={
            "psi_r_rms_wb": psi_r_rms,
            "psi_r_rel_pct": 100.0 * psi_r_rms / psi_mag,
            "is_rms_a": is_rms,
            "is_rel_pct": 100.0 * is_rms / is_mag,
            "psi_r_conv_time_ms": 1000.0 * convergence_time(t, flux_error_norm, 0.01 * psi_mag),
            "is_conv_time_ms": 1000.0 * convergence_time(t, current_error_norm, 1.0),
        },
    )


def result_from_kubota_d1(
    true_params: fo.MotorParams,
    op: fo.OperatingPoint,
    case: fo.ErrorCase,
    dt: float = 1.0e-5,
    t_end: float = 0.12,
    seed: int = 1,
    k_gain: float = KUBOTA_D1_K,
) -> MethodResult:
    obs_params = fo.apply_param_scale(true_params, case.param_scale)
    ss = fo.steady_operating_point(true_params, op)
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    i_true = np.array(ss["i_s"], dtype=float)
    u_true = np.array(ss["v_s"], dtype=float)
    phi_true = float(np.array(ss["psi_r"], dtype=float)[0])
    x_hat = np.array([0.6 * i_true[0], 1.4 * i_true[1], 0.5 * phi_true], dtype=float)

    current_error_norm = np.zeros(n, dtype=float)
    flux_error_norm = np.zeros(n, dtype=float)
    stable = True
    rng = np.random.default_rng(seed)
    v_offset = np.array([case.voltage_offset_d, case.voltage_offset_q], dtype=float)
    i_offset = np.array([case.current_offset_d, case.current_offset_q], dtype=float)

    for k in range(n):
        current_error_norm[k] = float(np.linalg.norm(x_hat[0:2] - i_true))
        flux_error_norm[k] = abs(x_hat[2] - phi_true)
        if k == n - 1:
            break
        noise = np.zeros(2, dtype=float)
        if case.current_noise_rms > 0.0:
            noise = rng.normal(0.0, case.current_noise_rms, size=2)
        u_meas = (1.0 + case.voltage_gain) * u_true + v_offset
        i_meas = (1.0 + case.current_gain) * i_true + i_offset + noise
        x_hat, _ = kubota_d1_step_with_observer_params(
            obs_params=obs_params,
            x_hat=x_hat,
            i_meas=i_meas,
            u_meas=u_meas,
            omega_r_e=float(ss["omega_r"]),
            dt=dt,
            k_gain=k_gain,
        )
        if not np.all(np.isfinite(x_hat)) or np.max(np.abs(x_hat)) > 1.0e6:
            stable = False
            current_error_norm = current_error_norm[: k + 1]
            flux_error_norm = flux_error_norm[: k + 1]
            t = t[: k + 1]
            break

    start = max(0, int(0.8 * len(t)))
    is_mag = max(float(np.linalg.norm(i_true)), 1.0e-12)
    phi_mag = max(abs(phi_true), 1.0e-12)
    psi_r_rms = float(np.sqrt(np.mean(flux_error_norm[start:] ** 2)))
    is_rms = float(np.sqrt(np.mean(current_error_norm[start:] ** 2)))
    return MethodResult(
        method=METHOD_D1,
        op=op,
        case=case,
        t=t,
        current_error_norm=current_error_norm,
        flux_error_norm=flux_error_norm,
        stable=stable,
        metrics={
            "psi_r_rms_wb": psi_r_rms,
            "psi_r_rel_pct": 100.0 * psi_r_rms / phi_mag,
            "is_rms_a": is_rms,
            "is_rel_pct": 100.0 * is_rms / is_mag,
            "psi_r_conv_time_ms": 1000.0 * convergence_time(t, flux_error_norm, 0.01 * phi_mag),
            "is_conv_time_ms": 1000.0 * convergence_time(t, current_error_norm, 1.0),
        },
    )


def result_from_yaskawa_e(
    true_params: fo.MotorParams,
    op: fo.OperatingPoint,
    case: fo.ErrorCase,
    dt: float = 1.0e-5,
    t_end: float = 0.12,
    seed: int = 1,
) -> MethodResult:
    obs_params = fo.apply_param_scale(true_params, case.param_scale)
    ss = fo.steady_operating_point(true_params, op)
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    i_dq = np.array(ss["i_s"], dtype=float)
    u_dq = np.array(ss["v_s"], dtype=float)
    phi_dq = np.array(ss["psi_r"], dtype=float)
    phi_mag = max(float(np.linalg.norm(phi_dq)), 1.0e-12)
    i_mag = max(float(np.linalg.norm(i_dq)), 1.0e-12)
    omega_k = float(ss["omega_k"])
    omega_m_e = float(ss["omega_r"])
    x_hat = np.array([0.6 * i_dq[0], 1.4 * i_dq[1], 0.5 * phi_dq[0], 0.25 * phi_mag], dtype=float)

    current_error_norm = np.zeros(n, dtype=float)
    flux_error_norm = np.zeros(n, dtype=float)
    stable = True
    rng = np.random.default_rng(seed)
    v_offset = np.array([case.voltage_offset_d, case.voltage_offset_q], dtype=float)
    i_offset = np.array([case.current_offset_d, case.current_offset_q], dtype=float)

    def sampled_vectors(t_sample: float, noise_vec: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        theta_sample = omega_k * t_sample
        i_sample = rotate_dq_to_stationary(theta_sample, i_dq)
        u_sample = rotate_dq_to_stationary(theta_sample, u_dq)
        phi_sample = rotate_dq_to_stationary(theta_sample, phi_dq)
        # ErrorCase offsets are defined in the same dq coordinates as the other methods.
        # Rotate them into the paper's stationary alpha-beta implementation.
        i_off_sample = rotate_dq_to_stationary(theta_sample, i_offset)
        v_off_sample = rotate_dq_to_stationary(theta_sample, v_offset)
        i_meas_sample = (1.0 + case.current_gain) * i_sample + i_off_sample + noise_vec
        u_meas_sample = (1.0 + case.voltage_gain) * u_sample + v_off_sample
        return i_sample, phi_sample, i_meas_sample, u_meas_sample

    for k in range(n):
        theta = omega_k * t[k]
        i_true = rotate_dq_to_stationary(theta, i_dq)
        phi_true = rotate_dq_to_stationary(theta, phi_dq)
        current_error_norm[k] = float(np.linalg.norm(x_hat[0:2] - i_true))
        flux_error_norm[k] = float(np.linalg.norm(x_hat[2:4] - phi_true))
        if k == n - 1:
            break

        noise = np.zeros(2, dtype=float)
        if case.current_noise_rms > 0.0:
            noise = rng.normal(0.0, case.current_noise_rms, size=2)
        _, _, i_meas_1, u_meas_1 = sampled_vectors(t[k], noise)
        _, _, i_meas_2, u_meas_2 = sampled_vectors(t[k] + 0.5 * dt, noise)
        _, _, i_meas_4, u_meas_4 = sampled_vectors(t[k] + dt, noise)
        k1 = yaskawa_e_derivative_with_observer_params(obs_params, x_hat, i_meas_1, u_meas_1, omega_m_e)
        k2 = yaskawa_e_derivative_with_observer_params(
            obs_params, x_hat + 0.5 * dt * k1, i_meas_2, u_meas_2, omega_m_e
        )
        k3 = yaskawa_e_derivative_with_observer_params(
            obs_params, x_hat + 0.5 * dt * k2, i_meas_2, u_meas_2, omega_m_e
        )
        k4 = yaskawa_e_derivative_with_observer_params(obs_params, x_hat + dt * k3, i_meas_4, u_meas_4, omega_m_e)
        x_hat = x_hat + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        if not np.all(np.isfinite(x_hat)) or np.max(np.abs(x_hat)) > 1.0e6:
            stable = False
            current_error_norm = current_error_norm[: k + 1]
            flux_error_norm = flux_error_norm[: k + 1]
            t = t[: k + 1]
            break

    start = max(0, int(0.8 * len(t)))
    psi_r_rms = float(np.sqrt(np.mean(flux_error_norm[start:] ** 2)))
    is_rms = float(np.sqrt(np.mean(current_error_norm[start:] ** 2)))
    return MethodResult(
        method=METHOD_E,
        op=op,
        case=case,
        t=t,
        current_error_norm=current_error_norm,
        flux_error_norm=flux_error_norm,
        stable=stable,
        metrics={
            "psi_r_rms_wb": psi_r_rms,
            "psi_r_rel_pct": 100.0 * psi_r_rms / phi_mag,
            "is_rms_a": is_rms,
            "is_rel_pct": 100.0 * is_rms / i_mag,
            "psi_r_conv_time_ms": 1000.0 * convergence_time(t, flux_error_norm, 0.01 * phi_mag),
            "is_conv_time_ms": 1000.0 * convergence_time(t, current_error_norm, 1.0),
        },
    )


def run_all() -> list[MethodResult]:
    params = fo.MotorParams()
    ops = fo.base_operating_points(params)
    cases = fo.evaluation_cases()
    results: list[MethodResult] = []
    for op in ops:
        print(f"Evaluating {op.name}")
        for idx, case in enumerate(cases):
            results.append(result_from_kubota_d1(params, op, case, seed=4000 + idx))
            results.append(result_from_sled23_t(params, op, case, seed=3000 + idx))
            results.append(result_from_yaskawa_e(params, op, case, seed=5000 + idx))
            results.append(result_from_flux_evaluation(METHOD_A, params, op, case))
            results.append(result_from_flux_evaluation(METHOD_B, params, op, case))
    return results


def save_combined_summary(results: list[MethodResult]) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / "gain_design_comparison_summary.csv"
    keys = list(results[0].metrics.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "operating_point", "case", "stable", *keys])
        for result in results:
            writer.writerow(
                [
                    result.method,
                    result.op.name,
                    result.case.name,
                    result.stable,
                    *[result.metrics[key] for key in keys],
                ]
            )
    return path


def save_nominal_table(results: list[MethodResult]) -> Path:
    path = DATA_DIR / "gain_design_nominal_comparison.csv"
    nominal = [r for r in results if r.case.name == "nominal"]
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "method",
                "operating_point",
                "stable",
                "psi_r_rel_pct",
                "is_rms_a",
                "psi_r_conv_time_ms",
                "is_conv_time_ms",
            ]
        )
        for r in nominal:
            writer.writerow(
                [
                    r.method,
                    r.op.name,
                    r.stable,
                    r.metrics["psi_r_rel_pct"],
                    r.metrics["is_rms_a"],
                    r.metrics["psi_r_conv_time_ms"],
                    r.metrics["is_conv_time_ms"],
                ]
            )
    return path


def save_worst_error_table(results: list[MethodResult]) -> Path:
    path = DATA_DIR / "gain_design_worst_error_comparison.csv"
    op_name = "5000rpm_-220Nm"
    rows: list[list[object]] = []
    for method in METHODS:
        method_rows = [r for r in results if r.method == method and r.op.name == op_name and r.case.name != "nominal"]
        worst_param = max(
            [r for r in method_rows if r.case.param_scale],
            key=lambda r: r.metrics["psi_r_rel_pct"],
        )
        worst_sensor = max(
            [r for r in method_rows if not r.case.param_scale],
            key=lambda r: r.metrics["psi_r_rel_pct"],
        )
        rows.append(
            [
                method,
                "parameter",
                worst_param.case.name,
                worst_param.stable,
                worst_param.metrics["psi_r_rel_pct"],
                worst_param.metrics["is_rms_a"],
            ]
        )
        rows.append(
            [
                method,
                "voltage_current",
                worst_sensor.case.name,
                worst_sensor.stable,
                worst_sensor.metrics["psi_r_rel_pct"],
                worst_sensor.metrics["is_rms_a"],
            ]
        )
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["method", "error_group", "worst_case", "stable", "psi_r_rel_pct", "is_rms_a"])
        writer.writerows(rows)
    return path


def plot_nominal_convergence(results: list[MethodResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nominal = [r for r in results if r.case.name == "nominal"]
    ops = [r.op.name for r in nominal if r.method == METHOD_D1]
    fig, axes = plt.subplots(len(ops), 2, figsize=(13.5, 9.8), sharex=False)
    colors = {
        METHOD_D1: "#6f4e7c",
        METHOD_C: "#0072b2",
        METHOD_E: "#009e73",
        METHOD_A: "#000000",
        METHOD_B: "#e69f00",
    }
    labels = {
        METHOD_D1: "D1 Kubota k=1.2",
        METHOD_C: "C T-type SLED23",
        METHOD_E: "E Yaskawa GS",
        METHOD_A: "A four-pole",
        METHOD_B: "B Hori 5.3",
    }
    for row, op_name in enumerate(ops):
        for method in METHODS:
            r = [x for x in nominal if x.op.name == op_name and x.method == method][0]
            t_ms = r.t * 1000.0
            axes[row, 0].semilogy(t_ms, np.maximum(r.current_error_norm, 1.0e-15), label=labels[method], color=colors[method])
            axes[row, 1].semilogy(t_ms, np.maximum(r.flux_error_norm, 1.0e-18), label=labels[method], color=colors[method])
        axes[row, 0].set_ylabel(f"{op_name}\ncurrent error [A]")
        axes[row, 1].set_ylabel("rotor flux error [Wb]")
        for ax in axes[row, :]:
            ax.grid(True, which="both", alpha=0.35)
            ax.legend(fontsize=8)
    axes[-1, 0].set_xlabel("time [ms]")
    axes[-1, 1].set_xlabel("time [ms]")
    fig.suptitle("Nominal convergence comparison of observer gain-design methods")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    path = OUT_DIR / "gain_design_nominal_convergence.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_worst_error_summary(results: list[MethodResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    op_name = "5000rpm_-220Nm"
    methods = METHODS
    labels = ["D1 Kubota", "C T-type", "E Yaskawa", "A four-pole", "B Hori 5.3"]
    param_psi = []
    sensor_psi = []
    param_i = []
    sensor_i = []
    for method in methods:
        rows = [r for r in results if r.method == method and r.op.name == op_name and r.case.name != "nominal"]
        param = max([r for r in rows if r.case.param_scale], key=lambda r: r.metrics["psi_r_rel_pct"])
        sensor = max([r for r in rows if not r.case.param_scale], key=lambda r: r.metrics["psi_r_rel_pct"])
        param_psi.append(param.metrics["psi_r_rel_pct"])
        sensor_psi.append(sensor.metrics["psi_r_rel_pct"])
        param_i.append(param.metrics["is_rms_a"])
        sensor_i.append(sensor.metrics["is_rms_a"])

    x = np.arange(len(methods))
    width = 0.35
    fig, axes = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    axes[0].bar(x - width / 2, param_psi, width, label="worst parameter error")
    axes[0].bar(x + width / 2, sensor_psi, width, label="worst voltage/current error")
    axes[1].bar(x - width / 2, param_i, width, label="worst parameter error")
    axes[1].bar(x + width / 2, sensor_i, width, label="worst voltage/current error")
    axes[0].set_ylabel("rotor flux RMS error [%]")
    axes[1].set_ylabel("stator current RMS error [A]")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels)
    for ax in axes:
        ax.grid(True, axis="y", alpha=0.35)
        ax.legend(fontsize=8)
    fig.suptitle("Worst error comparison at 5000 r/min, -220 Nm")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    path = OUT_DIR / "gain_design_worst_error_summary.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def print_report_tables(results: list[MethodResult]) -> None:
    print("\nNominal comparison")
    for r in [x for x in results if x.case.name == "nominal"]:
        print(
            f"{r.method:24s} {r.op.name:17s} stable={r.stable} "
            f"psi={r.metrics['psi_r_rel_pct']:.4g}% "
            f"is={r.metrics['is_rms_a']:.4g}A "
            f"tpsi={r.metrics['psi_r_conv_time_ms']:.4g}ms "
            f"ti={r.metrics['is_conv_time_ms']:.4g}ms"
        )
    print("\nWorst errors at 5000 rpm, -220 Nm")
    op_name = "5000rpm_-220Nm"
    for method in METHODS:
        rows = [r for r in results if r.method == method and r.op.name == op_name and r.case.name != "nominal"]
        param = max([r for r in rows if r.case.param_scale], key=lambda r: r.metrics["psi_r_rel_pct"])
        sensor = max([r for r in rows if not r.case.param_scale], key=lambda r: r.metrics["psi_r_rel_pct"])
        print(method, "param", param.case.name, param.metrics)
        print(method, "sensor", sensor.case.name, sensor.metrics)


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = run_all()
    paths = [
        save_combined_summary(results),
        save_nominal_table(results),
        save_worst_error_table(results),
        plot_nominal_convergence(results),
        plot_worst_error_summary(results),
    ]
    print_report_tables(results)
    for path in paths:
        print(f"Saved: {path}")


if __name__ == "__main__":
    main()
