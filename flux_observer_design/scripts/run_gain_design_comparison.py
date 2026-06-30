"""Compare the three implemented flux-observer gain-design methods.

Methods:

    A_four_pole_sylvester:
        Real four-state observer with four real poles and Sylvester equation.

    B_alpha_beta_sylvester:
        Real four-state observer with alpha/beta pole blocks and Sylvester
        equation.

    C_sled23_closed_form:
        SLED 2023 Appendix A observer in estimated rotor-flux coordinates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

import run_flux_observer_evaluation as fo
import run_sled23_appendix_observer_demo as sled


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"

METHOD_A = "A_four_pole_sylvester"
METHOD_B = "B_alpha_beta_sylvester"
METHOD_C = "C_sled23_closed_form"

SLED_ALPHA_I = 1000.0
SLED_ZETA_INF = 0.4


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


def convergence_time(t: np.ndarray, err: np.ndarray, threshold: float) -> float:
    ok = np.abs(err) <= threshold
    if not np.any(ok):
        return float("nan")
    suffix_ok = np.flip(np.cumprod(np.flip(ok).astype(int))).astype(bool)
    idx = np.argmax(suffix_ok)
    if suffix_ok[idx]:
        return float(t[idx])
    return float("nan")


def result_from_sylvester(
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


def sled23_step_with_observer_params(
    true_params: fo.MotorParams,
    obs_params: fo.MotorParams,
    x_hat: np.ndarray,
    i_true: np.ndarray,
    i_meas: np.ndarray,
    u_meas: np.ndarray,
    omega_m: float,
    dt: float,
) -> tuple[np.ndarray, float]:
    inv = sled.inverse_gamma_params(obs_params)
    isd_hat, isq_hat, psi_hat = x_hat
    eisd = i_meas[0] - isd_hat
    eisq = i_meas[1] - isq_hat

    b = 2.0 * SLED_ZETA_INF * abs(omega_m) + inv.alpha
    gamma = SLED_ALPHA_I - inv.alpha
    k_den = inv.alpha * inv.alpha + omega_m * omega_m
    k1 = b * inv.alpha / k_den
    k2 = b * omega_m / k_den

    denominator = psi_hat - inv.l_sigma * eisd
    if abs(denominator) < 1.0e-9:
        denominator = math.copysign(1.0e-9, denominator if denominator != 0.0 else 1.0)
    omega_s = omega_m + (
        inv.rr * i_meas[1] + k2 * SLED_ALPHA_I * inv.l_sigma * eisd - gamma * inv.l_sigma * eisq
    ) / denominator

    disd = (
        inv.alpha * psi_hat
        - inv.r_sigma * isd_hat
        + omega_s * inv.l_sigma * isq_hat
        + u_meas[0]
        + inv.l_sigma * (gamma * eisd - omega_m * eisq)
    ) / inv.l_sigma
    disq = (
        -omega_m * psi_hat
        - inv.r_sigma * isq_hat
        - omega_s * inv.l_sigma * isd_hat
        + u_meas[1]
        + inv.l_sigma * (gamma * eisq + omega_m * eisd)
    ) / inv.l_sigma
    dpsi = (
        -inv.alpha * psi_hat
        + inv.rr * isd_hat
        + (k1 * SLED_ALPHA_I - gamma) * inv.l_sigma * eisd
        - (omega_s - omega_m) * inv.l_sigma * eisq
    )

    _ = true_params  # Keeps the signature explicit: true and observer params are intentionally separated.
    return x_hat + dt * np.array([disd, disq, dpsi], dtype=float), omega_s


def result_from_sled23(
    true_params: fo.MotorParams,
    op: fo.OperatingPoint,
    case: fo.ErrorCase,
    dt: float = 1.0e-5,
    t_end: float = 0.12,
    seed: int = 1,
) -> MethodResult:
    obs_params = fo.apply_param_scale(true_params, case.param_scale)
    ss = sled.steady_rotor_flux_frame(true_params, op)
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    i_true = np.array([ss["id"], ss["iq"]], dtype=float)
    u_true = np.array([ss["usd"], ss["usq"]], dtype=float)
    psi_true = float(ss["psi_r"])
    x_hat = np.array([0.6 * ss["id"], 1.4 * ss["iq"], 0.5 * psi_true], dtype=float)

    current_error_norm = np.zeros(n, dtype=float)
    flux_error_norm = np.zeros(n, dtype=float)
    stable = True
    rng = np.random.default_rng(seed)
    v_offset = np.array([case.voltage_offset_d, case.voltage_offset_q], dtype=float)
    i_offset = np.array([case.current_offset_d, case.current_offset_q], dtype=float)

    for k in range(n):
        current_error_norm[k] = float(np.linalg.norm(x_hat[0:2] - i_true))
        flux_error_norm[k] = abs(x_hat[2] - psi_true)
        if k == n - 1:
            break
        noise = np.zeros(2, dtype=float)
        if case.current_noise_rms > 0.0:
            noise = rng.normal(0.0, case.current_noise_rms, size=2)
        u_meas = (1.0 + case.voltage_gain) * u_true + v_offset
        i_meas = (1.0 + case.current_gain) * i_true + i_offset + noise
        x_hat, _ = sled23_step_with_observer_params(
            true_params=true_params,
            obs_params=obs_params,
            x_hat=x_hat,
            i_true=i_true,
            i_meas=i_meas,
            u_meas=u_meas,
            omega_m=float(ss["omega_m"]),
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
    psi_mag = max(abs(psi_true), 1.0e-12)
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


def run_all() -> list[MethodResult]:
    params = fo.MotorParams()
    ops = fo.base_operating_points(params)
    cases = fo.evaluation_cases()
    results: list[MethodResult] = []
    for op in ops:
        print(f"Evaluating {op.name}")
        for idx, case in enumerate(cases):
            results.append(result_from_sylvester(METHOD_A, params, op, case))
            results.append(result_from_sylvester(METHOD_B, params, op, case))
            results.append(result_from_sled23(params, op, case, seed=3000 + idx))
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
    for method in [METHOD_A, METHOD_B, METHOD_C]:
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
    ops = [r.op.name for r in nominal if r.method == METHOD_A]
    fig, axes = plt.subplots(len(ops), 2, figsize=(13, 9), sharex=False)
    colors = {METHOD_A: "#1f77b4", METHOD_B: "#d62728", METHOD_C: "#2ca02c"}
    labels = {METHOD_A: "A four-pole", METHOD_B: "B alpha/beta", METHOD_C: "C SLED23"}
    for row, op_name in enumerate(ops):
        for method in [METHOD_A, METHOD_B, METHOD_C]:
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
    fig.suptitle("Nominal convergence comparison of three observer gain-design methods")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    path = OUT_DIR / "gain_design_nominal_convergence.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def plot_worst_error_summary(results: list[MethodResult]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    op_name = "5000rpm_-220Nm"
    methods = [METHOD_A, METHOD_B, METHOD_C]
    labels = ["A four-pole", "B alpha/beta", "C SLED23"]
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
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
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
    for method in [METHOD_A, METHOD_B, METHOD_C]:
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
