"""Rotating-dq full-order flux observer design and evaluation.

The observer state is

    x = [psi_s, psi_r]^T

where psi_s and psi_r are complex dq vectors.  This is equivalent to the
four-real-state vector

    [psi_sd, psi_sq, psi_rd, psi_rq]^T.

Only stator current is used for correction.  Stator and rotor currents are
calculated from the estimated stator and rotor fluxes.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Iterable

import csv
import math

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"


@dataclass(frozen=True)
class MotorParams:
    rs: float = 0.00762
    rr: float = 0.008041
    lls: float = 0.0000419
    llr: float = 0.0000419
    lm: float = 0.0001608
    pole_pairs: int = 4

    @property
    def ls(self) -> float:
        return self.lls + self.lm

    @property
    def lr(self) -> float:
        return self.llr + self.lm

    @property
    def det_l(self) -> float:
        return self.ls * self.lr - self.lm * self.lm

    @property
    def sigma_ls(self) -> float:
        return self.det_l / self.lr


@dataclass(frozen=True)
class OperatingPoint:
    name: str
    speed_rpm: float
    torque_nm: float
    id_ref: float


@dataclass(frozen=True)
class ErrorCase:
    name: str
    param_scale: dict[str, float] | None = None
    voltage_gain: float = 0.0
    voltage_offset_d: float = 0.0
    voltage_offset_q: float = 0.0
    current_gain: float = 0.0
    current_offset_d: float = 0.0
    current_offset_q: float = 0.0
    current_noise_rms: float = 0.0


@dataclass
class SimResult:
    op: OperatingPoint
    case: ErrorCase
    t: np.ndarray
    true_flux: np.ndarray
    est_flux: np.ndarray
    true_current: np.ndarray
    est_current: np.ndarray
    poles: np.ndarray
    stable: bool
    metrics: dict[str, float]


def apply_param_scale(params: MotorParams, scale: dict[str, float] | None) -> MotorParams:
    if not scale:
        return params
    values = {
        "rs": params.rs,
        "rr": params.rr,
        "lls": params.lls,
        "llr": params.llr,
        "lm": params.lm,
        "pole_pairs": params.pole_pairs,
    }
    for key, factor in scale.items():
        if key not in values:
            raise KeyError(f"unknown parameter {key}")
        values[key] = values[key] * factor
    return MotorParams(**values)


def current_rows(params: MotorParams) -> tuple[np.ndarray, np.ndarray]:
    """Return complex rows that calculate stator and rotor current from flux."""
    d = params.det_l
    c_s = np.array([params.lr / d, -params.lm / d], dtype=complex)
    c_r = np.array([-params.lm / d, params.ls / d], dtype=complex)
    return c_s, c_r


def state_matrix(params: MotorParams, omega_r: float, omega_k: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return A, B, C_s for the complex rotating-dq model.

    Flux equations:

        d psi_s / dt = v_s - Rs i_s - j omega_k psi_s
        d psi_r / dt =     - Rr i_r - j(omega_k - omega_r) psi_r
    """
    c_s, c_r = current_rows(params)
    a = np.zeros((2, 2), dtype=complex)
    a[0, :] = -params.rs * c_s
    a[0, 0] += -1j * omega_k
    a[1, :] = -params.rr * c_r
    a[1, 1] += -1j * (omega_k - omega_r)
    b = np.array([1.0 + 0j, 0.0 + 0j], dtype=complex)
    return a, b, c_s


def observer_gain_by_pole_placement(
    a: np.ndarray, c: np.ndarray, desired_poles: Iterable[complex]
) -> np.ndarray:
    """Place the two complex observer poles of A - l*C.

    The equivalent real dq observer has the placed poles and their conjugates.
    For a 2-state complex SISO observer, trace and determinant matching give a
    linear equation in l:

        C*l = trace(A) - trace_desired
        C*adj(A)*l = det(A) - det_desired
    """
    p = np.array(list(desired_poles), dtype=complex)
    if p.shape != (2,):
        raise ValueError("exactly two complex poles are required")

    desired_trace = p[0] + p[1]
    desired_det = p[0] * p[1]
    adj = np.array([[a[1, 1], -a[0, 1]], [-a[1, 0], a[0, 0]]], dtype=complex)
    lhs = np.vstack([c, c @ adj])
    rhs = np.array([np.trace(a) - desired_trace, np.linalg.det(a) - desired_det], dtype=complex)
    gain = np.linalg.solve(lhs, rhs)
    return gain


def real_observer_gain(gain_complex: np.ndarray) -> np.ndarray:
    """Return the 4x2 real observer gain equivalent to the complex gain."""
    l_s, l_r = gain_complex
    return np.array(
        [
            [l_s.real, -l_s.imag],
            [l_s.imag, l_s.real],
            [l_r.real, -l_r.imag],
            [l_r.imag, l_r.real],
        ],
        dtype=float,
    )


def desired_observer_poles(observer_bandwidth: float = 2200.0) -> np.ndarray:
    """Observer pole pair in the complex two-state model."""
    return np.array([-observer_bandwidth, -1.55 * observer_bandwidth], dtype=complex)


def steady_operating_point(params: MotorParams, op: OperatingPoint) -> dict[str, complex | float]:
    omega_m = op.speed_rpm * 2.0 * math.pi / 60.0
    omega_r = params.pole_pairs * omega_m
    id_s = op.id_ref
    torque_coeff = 1.5 * params.pole_pairs * params.lm * params.lm / params.lr
    iq_s = op.torque_nm / (torque_coeff * id_s)
    slip = params.rr / params.lr * iq_s / id_s
    omega_k = omega_r + slip

    psi_r = params.lm * id_s + 0j
    psi_s = params.ls * id_s + 1j * params.sigma_ls * iq_s
    i_s = id_s + 1j * iq_s
    _, c_r = current_rows(params)
    i_r = c_r @ np.array([psi_s, psi_r], dtype=complex)
    v_s = params.rs * i_s + 1j * omega_k * psi_s
    torque = 1.5 * params.pole_pairs * (psi_s.real * i_s.imag - psi_s.imag * i_s.real)

    return {
        "omega_m": omega_m,
        "omega_r": omega_r,
        "omega_k": omega_k,
        "slip": slip,
        "psi_s": psi_s,
        "psi_r": psi_r,
        "i_s": i_s,
        "i_r": i_r,
        "v_s": v_s,
        "torque": torque,
    }


def all_currents(params: MotorParams, x: np.ndarray) -> tuple[complex, complex]:
    c_s, c_r = current_rows(params)
    return c_s @ x, c_r @ x


def rms_complex(x: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.abs(x) ** 2)))


def convergence_time(t: np.ndarray, err: np.ndarray, threshold: float) -> float:
    ok = np.abs(err) <= threshold
    if not np.any(ok):
        return float("nan")
    # First time after which all remaining samples stay below threshold.
    suffix_ok = np.flip(np.cumprod(np.flip(ok).astype(int))).astype(bool)
    idx = np.argmax(suffix_ok)
    if suffix_ok[idx]:
        return float(t[idx])
    return float("nan")


def simulate_case(
    true_params: MotorParams,
    op: OperatingPoint,
    case: ErrorCase,
    dt: float = 1.0e-5,
    t_end: float = 0.12,
    seed: int = 1,
) -> SimResult:
    obs_params = apply_param_scale(true_params, case.param_scale)
    ss = steady_operating_point(true_params, op)
    omega_r = float(ss["omega_r"])
    omega_k = float(ss["omega_k"])
    v_true = complex(ss["v_s"])

    a_true, b_true, _ = state_matrix(true_params, omega_r, omega_k)
    a_obs, b_obs, c_obs = state_matrix(obs_params, omega_r, omega_k)
    poles = desired_observer_poles()
    gain = observer_gain_by_pole_placement(a_obs, c_obs, poles)
    placed_poles = np.linalg.eigvals(a_obs - np.outer(gain, c_obs))

    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    x_true = np.array([complex(ss["psi_s"]), complex(ss["psi_r"])], dtype=complex)
    x_hat = np.array(
        [
            complex(ss["psi_s"]) * (0.25 + 0.20j),
            complex(ss["psi_r"]) * (1.45 - 0.25j),
        ],
        dtype=complex,
    )

    true_flux = np.zeros((n, 2), dtype=complex)
    est_flux = np.zeros((n, 2), dtype=complex)
    true_current = np.zeros((n, 2), dtype=complex)
    est_current = np.zeros((n, 2), dtype=complex)

    rng = np.random.default_rng(seed)
    stable = True

    v_offset = case.voltage_offset_d + 1j * case.voltage_offset_q
    i_offset = case.current_offset_d + 1j * case.current_offset_q

    for k in range(n):
        is_true, ir_true = all_currents(true_params, x_true)
        is_hat, ir_hat = all_currents(obs_params, x_hat)

        true_flux[k, :] = x_true
        est_flux[k, :] = x_hat
        true_current[k, :] = [is_true, ir_true]
        est_current[k, :] = [is_hat, ir_hat]

        if k == n - 1:
            break

        noise = 0j
        if case.current_noise_rms > 0.0:
            noise = rng.normal(0.0, case.current_noise_rms) + 1j * rng.normal(
                0.0, case.current_noise_rms
            )
        v_meas = (1.0 + case.voltage_gain) * v_true + v_offset
        i_meas = (1.0 + case.current_gain) * is_true + i_offset + noise

        # True plant is integrated as well, although the chosen input is a steady-state input.
        dx_true = a_true @ x_true + b_true * v_true
        x_true = x_true + dt * dx_true

        innovation = i_meas - (c_obs @ x_hat)
        dx_hat = a_obs @ x_hat + b_obs * v_meas + gain * innovation
        x_hat = x_hat + dt * dx_hat

        if not np.all(np.isfinite(x_hat)) or np.max(np.abs(x_hat)) > 10.0:
            stable = False
            true_flux = true_flux[: k + 1, :]
            est_flux = est_flux[: k + 1, :]
            true_current = true_current[: k + 1, :]
            est_current = est_current[: k + 1, :]
            t = t[: k + 1]
            break

    start = max(0, int(0.8 * len(t)))
    flux_err = est_flux - true_flux
    curr_err = est_current - true_current
    psi_s_mag = max(rms_complex(true_flux[start:, 0]), 1e-12)
    psi_r_mag = max(rms_complex(true_flux[start:, 1]), 1e-12)
    is_mag = max(rms_complex(true_current[start:, 0]), 1e-12)
    ir_mag = max(rms_complex(true_current[start:, 1]), 1e-12)

    metrics = {
        "psi_s_rms_wb": rms_complex(flux_err[start:, 0]),
        "psi_r_rms_wb": rms_complex(flux_err[start:, 1]),
        "is_rms_a": rms_complex(curr_err[start:, 0]),
        "ir_rms_a": rms_complex(curr_err[start:, 1]),
        "psi_s_rel_pct": 100.0 * rms_complex(flux_err[start:, 0]) / psi_s_mag,
        "psi_r_rel_pct": 100.0 * rms_complex(flux_err[start:, 1]) / psi_r_mag,
        "is_rel_pct": 100.0 * rms_complex(curr_err[start:, 0]) / is_mag,
        "ir_rel_pct": 100.0 * rms_complex(curr_err[start:, 1]) / ir_mag,
        "psi_r_conv_time_ms": 1000.0
        * convergence_time(t, flux_err[:, 1], 0.01 * max(np.abs(true_flux[-1, 1]), 1e-12)),
        "is_conv_time_ms": 1000.0 * convergence_time(t, curr_err[:, 0], 1.0),
        "max_abs_flux_wb": float(np.max(np.abs(est_flux))),
    }

    return SimResult(
        op=op,
        case=case,
        t=t,
        true_flux=true_flux,
        est_flux=est_flux,
        true_current=true_current,
        est_current=est_current,
        poles=placed_poles,
        stable=stable,
        metrics=metrics,
    )


def base_operating_points(params: MotorParams) -> list[OperatingPoint]:
    id_ref = 550.0 * math.sqrt(3.0) / math.sqrt(2.0)
    return [
        OperatingPoint("5000rpm_+220Nm", 5000.0, +220.0, id_ref),
        OperatingPoint("5000rpm_-220Nm", 5000.0, -220.0, id_ref),
        OperatingPoint("1000rpm_-220Nm", 1000.0, -220.0, id_ref),
    ]


def evaluation_cases() -> list[ErrorCase]:
    cases = [ErrorCase("nominal")]
    for key in ["rs", "rr", "lm", "lls", "llr"]:
        for pct in [-20, -10, -5, 5, 10, 20]:
            cases.append(ErrorCase(f"{key}_{pct:+d}pct", param_scale={key: 1.0 + pct / 100.0}))
    cases.extend(
        [
            ErrorCase("voltage_gain_+5pct", voltage_gain=+0.05),
            ErrorCase("voltage_gain_-5pct", voltage_gain=-0.05),
            ErrorCase("voltage_offset_+1Vd_-1Vq", voltage_offset_d=+1.0, voltage_offset_q=-1.0),
            ErrorCase("current_gain_+1pct", current_gain=+0.01),
            ErrorCase("current_gain_-1pct", current_gain=-0.01),
            ErrorCase("current_offset_+1Ad_-1Aq", current_offset_d=+1.0, current_offset_q=-1.0),
            ErrorCase("current_noise_1Arms", current_noise_rms=1.0),
        ]
    )
    return cases


def save_summary_csv(results: list[SimResult]) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / "evaluation_summary.csv"
    metric_keys = list(results[0].metrics.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["operating_point", "case", "stable", *metric_keys])
        for r in results:
            writer.writerow([r.op.name, r.case.name, r.stable, *[r.metrics[k] for k in metric_keys]])
    return path


def plot_nominal_convergence(results: list[SimResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    nominal = [r for r in results if r.case.name == "nominal"]
    fig, axes = plt.subplots(3, 2, figsize=(13, 9), sharex=False)
    for row, r in enumerate(nominal):
        t_ms = r.t * 1000.0
        flux_err_s = np.abs(r.est_flux[:, 0] - r.true_flux[:, 0])
        flux_err_r = np.abs(r.est_flux[:, 1] - r.true_flux[:, 1])
        curr_err_s = np.abs(r.est_current[:, 0] - r.true_current[:, 0])
        curr_err_r = np.abs(r.est_current[:, 1] - r.true_current[:, 1])
        ax = axes[row, 0]
        ax.semilogy(t_ms, flux_err_s, label="primary flux")
        ax.semilogy(t_ms, flux_err_r, label="secondary flux")
        ax.grid(True, which="both", alpha=0.35)
        ax.set_ylabel(f"{r.op.name}\nflux error [Wb]")
        ax.legend(loc="upper right")
        ax = axes[row, 1]
        ax.semilogy(t_ms, curr_err_s, label="primary current")
        ax.semilogy(t_ms, curr_err_r, label="secondary current")
        ax.grid(True, which="both", alpha=0.35)
        ax.set_ylabel("current error [A]")
        ax.legend(loc="upper right")
    axes[-1, 0].set_xlabel("time [ms]")
    axes[-1, 1].set_xlabel("time [ms]")
    fig.suptitle("Nominal observer convergence from initial flux error")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    path = OUT_DIR / "nominal_convergence.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_pole_map(params: MotorParams, ops: list[OperatingPoint]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 5))
    for op in ops:
        ss = steady_operating_point(params, op)
        a, _, c = state_matrix(params, float(ss["omega_r"]), float(ss["omega_k"]))
        g = observer_gain_by_pole_placement(a, c, desired_observer_poles())
        natural = np.linalg.eigvals(a)
        placed = np.linalg.eigvals(a - np.outer(g, c))
        ax.plot(natural.real, natural.imag, "x", ms=8, label=f"{op.name} natural")
        ax.plot(placed.real, placed.imag, "o", ms=5, label=f"{op.name} observer")
    ax.axvline(0.0, color="k", lw=0.8)
    ax.grid(True, alpha=0.35)
    ax.set_xlabel("real part [1/s]")
    ax.set_ylabel("imaginary part [rad/s]")
    ax.set_title("Frozen-time complex observer pole placement")
    ax.legend(fontsize=7, ncols=2)
    fig.tight_layout()
    path = OUT_DIR / "observer_pole_map.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_parameter_error_sweep(results: list[SimResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    op_name = "5000rpm_-220Nm"
    params = ["rs", "rr", "lm", "lls", "llr"]
    pct_values = [-20, -10, -5, 5, 10, 20]
    fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    for param in params:
        xs = []
        psi_r = []
        is_err = []
        for pct in pct_values:
            name = f"{param}_{pct:+d}pct"
            match = [r for r in results if r.op.name == op_name and r.case.name == name][0]
            xs.append(pct)
            psi_r.append(match.metrics["psi_r_rel_pct"])
            is_err.append(match.metrics["is_rms_a"])
        axes[0].plot(xs, psi_r, marker="o", label=param)
        axes[1].plot(xs, is_err, marker="o", label=param)
    axes[0].set_ylabel("secondary flux RMS error [%]")
    axes[1].set_ylabel("primary current RMS error [A]")
    axes[1].set_xlabel("observer parameter error [%]")
    for ax in axes:
        ax.grid(True, alpha=0.35)
        ax.legend(ncols=5, fontsize=8)
    fig.suptitle("Parameter error sensitivity at 5000 r/min, -220 Nm")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    path = OUT_DIR / "parameter_error_sweep.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_sensor_error_summary(results: list[SimResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    op_name = "5000rpm_-220Nm"
    sensor_names = [
        "nominal",
        "voltage_gain_+5pct",
        "voltage_gain_-5pct",
        "voltage_offset_+1Vd_-1Vq",
        "current_gain_+1pct",
        "current_gain_-1pct",
        "current_offset_+1Ad_-1Aq",
        "current_noise_1Arms",
    ]
    rows = [r for r in results if r.op.name == op_name and r.case.name in sensor_names]
    rows.sort(key=lambda r: sensor_names.index(r.case.name))
    labels = [r.case.name.replace("_", "\n") for r in rows]
    psi_r = [r.metrics["psi_r_rel_pct"] for r in rows]
    is_err = [r.metrics["is_rms_a"] for r in rows]
    x = np.arange(len(rows))
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    axes[0].bar(x, psi_r, color="#4c78a8")
    axes[1].bar(x, is_err, color="#f58518")
    axes[0].set_ylabel("secondary flux RMS error [%]")
    axes[1].set_ylabel("primary current RMS error [A]")
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(labels, rotation=0, fontsize=8)
    for ax in axes:
        ax.grid(True, axis="y", alpha=0.35)
    fig.suptitle("Voltage/current measurement error sensitivity at 5000 r/min, -220 Nm")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    path = OUT_DIR / "sensor_error_summary.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_nominal_waveform(results: list[SimResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    r = [item for item in results if item.op.name == "5000rpm_-220Nm" and item.case.name == "nominal"][0]
    t_ms = r.t * 1000.0
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharex=True)
    axes[0, 0].plot(t_ms, r.true_flux[:, 0].real, label="psi_sd true")
    axes[0, 0].plot(t_ms, r.est_flux[:, 0].real, "--", label="psi_sd obs")
    axes[0, 0].plot(t_ms, r.true_flux[:, 0].imag, label="psi_sq true")
    axes[0, 0].plot(t_ms, r.est_flux[:, 0].imag, "--", label="psi_sq obs")
    axes[0, 1].plot(t_ms, r.true_flux[:, 1].real, label="psi_rd true")
    axes[0, 1].plot(t_ms, r.est_flux[:, 1].real, "--", label="psi_rd obs")
    axes[0, 1].plot(t_ms, r.true_flux[:, 1].imag, label="psi_rq true")
    axes[0, 1].plot(t_ms, r.est_flux[:, 1].imag, "--", label="psi_rq obs")
    axes[1, 0].plot(t_ms, r.true_current[:, 0].real, label="isd true")
    axes[1, 0].plot(t_ms, r.est_current[:, 0].real, "--", label="isd obs")
    axes[1, 0].plot(t_ms, r.true_current[:, 0].imag, label="isq true")
    axes[1, 0].plot(t_ms, r.est_current[:, 0].imag, "--", label="isq obs")
    axes[1, 1].plot(t_ms, r.true_current[:, 1].real, label="ird true")
    axes[1, 1].plot(t_ms, r.est_current[:, 1].real, "--", label="ird obs")
    axes[1, 1].plot(t_ms, r.true_current[:, 1].imag, label="irq true")
    axes[1, 1].plot(t_ms, r.est_current[:, 1].imag, "--", label="irq obs")
    titles = ["primary flux", "secondary flux", "primary current", "secondary current"]
    for ax, title in zip(axes.ravel(), titles):
        ax.set_title(title)
        ax.grid(True, alpha=0.35)
        ax.legend(fontsize=8, ncols=2)
    axes[1, 0].set_xlabel("time [ms]")
    axes[1, 1].set_xlabel("time [ms]")
    fig.suptitle("Nominal waveform convergence at 5000 r/min, -220 Nm")
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    path = OUT_DIR / "nominal_waveform_5000rpm_regen.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def print_key_tables(results: list[SimResult]) -> None:
    print("\nNominal cases")
    print("op, psi_r_rel_pct, is_rms_a, psi_r_conv_ms, is_conv_ms, stable")
    for r in results:
        if r.case.name == "nominal":
            print(
                f"{r.op.name}, {r.metrics['psi_r_rel_pct']:.6g}, "
                f"{r.metrics['is_rms_a']:.6g}, {r.metrics['psi_r_conv_time_ms']:.3g}, "
                f"{r.metrics['is_conv_time_ms']:.3g}, {r.stable}"
            )

    print("\nWorst cases by secondary flux relative error")
    sorted_rows = sorted(results, key=lambda r: r.metrics["psi_r_rel_pct"], reverse=True)
    for r in sorted_rows[:10]:
        print(
            f"{r.op.name:17s} {r.case.name:28s} "
            f"psi_r={r.metrics['psi_r_rel_pct']:.3f}% "
            f"is={r.metrics['is_rms_a']:.3f}A stable={r.stable}"
        )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    params = MotorParams()
    ops = base_operating_points(params)
    cases = evaluation_cases()

    results: list[SimResult] = []
    for op in ops:
        ss = steady_operating_point(params, op)
        print(
            f"{op.name}: torque={ss['torque']:.3f} Nm, "
            f"omega_k={ss['omega_k']:.3f} rad/s, slip={ss['slip']:.3f} rad/s"
        )
        for idx, case in enumerate(cases):
            results.append(simulate_case(params, op, case, seed=1000 + idx))

    summary_path = save_summary_csv(results)
    figs = [
        plot_nominal_waveform(results),
        plot_nominal_convergence(results),
        plot_pole_map(params, ops),
        plot_parameter_error_sweep(results),
        plot_sensor_error_summary(results),
    ]
    print_key_tables(results)
    print(f"\nSaved summary: {summary_path}")
    for fig in figs:
        print(f"Saved figure:  {fig}")


if __name__ == "__main__":
    main()
