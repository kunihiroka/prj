"""Demo of the SLED 2023 Appendix A full-order observer.

The observer is implemented in the estimated rotor-flux coordinates. Its state
is the estimated stator current and the estimated rotor-flux magnitude:

    x_hat = [isd_hat, isq_hat, psi_R_hat].

The equations follow Appendix A, (19a)-(19d), of:

    Tiitinen, Hinkkanen, Harnefors,
    "Speed-Adaptive Full-Order Observer Revisited: Closed-Form Design for
    Induction Motor Drives", IEEE SLED 2023.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import math

import matplotlib.pyplot as plt
import numpy as np

from run_flux_observer_evaluation import MotorParams, OperatingPoint, base_operating_points


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "figures"
DATA_DIR = ROOT / "data"


@dataclass(frozen=True)
class InverseGammaParams:
    l_sigma: float
    lm: float
    rr: float
    r_sigma: float
    alpha: float


@dataclass
class Sled23Result:
    op: OperatingPoint
    t: np.ndarray
    true_state: np.ndarray
    est_state: np.ndarray
    omega_s: np.ndarray
    stable: bool
    metrics: dict[str, float]


def inverse_gamma_params(params: MotorParams) -> InverseGammaParams:
    lr = params.lr
    l_sigma = params.det_l / lr
    lm_inv_gamma = params.lm * params.lm / lr
    rr_inv_gamma = params.rr * (params.lm / lr) ** 2
    return InverseGammaParams(
        l_sigma=l_sigma,
        lm=lm_inv_gamma,
        rr=rr_inv_gamma,
        r_sigma=params.rs + rr_inv_gamma,
        alpha=rr_inv_gamma / lm_inv_gamma,
    )


def steady_rotor_flux_frame(params: MotorParams, op: OperatingPoint) -> dict[str, float]:
    inv = inverse_gamma_params(params)
    omega_m = params.pole_pairs * op.speed_rpm * 2.0 * math.pi / 60.0
    id_s = op.id_ref
    torque_coeff = 1.5 * params.pole_pairs * params.lm * params.lm / params.lr
    iq_s = op.torque_nm / (torque_coeff * id_s)
    psi_r = inv.lm * id_s
    omega_s = omega_m + inv.rr * iq_s / psi_r

    # Steady-state input voltage in the estimated rotor-flux coordinates.
    usd = inv.r_sigma * id_s - omega_s * inv.l_sigma * iq_s - inv.alpha * psi_r
    usq = omega_m * psi_r + inv.r_sigma * iq_s + omega_s * inv.l_sigma * id_s
    return {
        "omega_m": omega_m,
        "id": id_s,
        "iq": iq_s,
        "psi_r": psi_r,
        "omega_s": omega_s,
        "usd": usd,
        "usq": usq,
    }


def sled23_step(
    params: MotorParams,
    x_hat: np.ndarray,
    i_meas: np.ndarray,
    u_s: np.ndarray,
    omega_m: float,
    alpha_i: float,
    zeta_inf: float,
    dt: float,
) -> tuple[np.ndarray, float]:
    inv = inverse_gamma_params(params)
    isd_hat, isq_hat, psi_hat = x_hat
    eisd = i_meas[0] - isd_hat
    eisq = i_meas[1] - isq_hat

    b = 2.0 * zeta_inf * abs(omega_m) + inv.alpha
    gamma = alpha_i - inv.alpha
    k_den = inv.alpha * inv.alpha + omega_m * omega_m
    k1 = b * inv.alpha / k_den
    k2 = b * omega_m / k_den

    denominator = psi_hat - inv.l_sigma * eisd
    if abs(denominator) < 1e-9:
        denominator = math.copysign(1e-9, denominator if denominator != 0.0 else 1.0)
    omega_s = omega_m + (
        inv.rr * i_meas[1] + k2 * alpha_i * inv.l_sigma * eisd - gamma * inv.l_sigma * eisq
    ) / denominator

    disd = (
        inv.alpha * psi_hat
        - inv.r_sigma * isd_hat
        + omega_s * inv.l_sigma * isq_hat
        + u_s[0]
        + inv.l_sigma * (gamma * eisd - omega_m * eisq)
    ) / inv.l_sigma
    disq = (
        -omega_m * psi_hat
        - inv.r_sigma * isq_hat
        - omega_s * inv.l_sigma * isd_hat
        + u_s[1]
        + inv.l_sigma * (gamma * eisq + omega_m * eisd)
    ) / inv.l_sigma
    dpsi = (
        -inv.alpha * psi_hat
        + inv.rr * isd_hat
        + (k1 * alpha_i - gamma) * inv.l_sigma * eisd
        - (omega_s - omega_m) * inv.l_sigma * eisq
    )
    return x_hat + dt * np.array([disd, disq, dpsi], dtype=float), omega_s


def simulate_case(
    params: MotorParams,
    op: OperatingPoint,
    dt: float = 1.0e-5,
    t_end: float = 0.08,
    alpha_i: float = 1000.0,
    zeta_inf: float = 0.4,
) -> Sled23Result:
    ss = steady_rotor_flux_frame(params, op)
    n = int(round(t_end / dt)) + 1
    t = np.linspace(0.0, t_end, n)
    true_state = np.tile(np.array([ss["id"], ss["iq"], ss["psi_r"]], dtype=float), (n, 1))
    est_state = np.zeros((n, 3), dtype=float)
    omega_s = np.zeros(n, dtype=float)
    x_hat = np.array([0.6 * ss["id"], 1.4 * ss["iq"], 0.5 * ss["psi_r"]], dtype=float)
    stable = True

    for k in range(n):
        est_state[k, :] = x_hat
        omega_s[k] = ss["omega_s"] if k == 0 else omega_s[k - 1]
        if k == n - 1:
            break
        x_hat, omega_s[k] = sled23_step(
            params=params,
            x_hat=x_hat,
            i_meas=np.array([ss["id"], ss["iq"]], dtype=float),
            u_s=np.array([ss["usd"], ss["usq"]], dtype=float),
            omega_m=ss["omega_m"],
            alpha_i=alpha_i,
            zeta_inf=zeta_inf,
            dt=dt,
        )
        if not np.all(np.isfinite(x_hat)) or np.max(np.abs(x_hat)) > 1.0e6:
            stable = False
            est_state = est_state[: k + 1, :]
            true_state = true_state[: k + 1, :]
            omega_s = omega_s[: k + 1]
            t = t[: k + 1]
            break

    err = est_state - true_state
    start = max(0, int(0.8 * len(t)))
    metrics = {
        "id_rms_a": float(np.sqrt(np.mean(err[start:, 0] ** 2))),
        "iq_rms_a": float(np.sqrt(np.mean(err[start:, 1] ** 2))),
        "psi_r_rms_wb": float(np.sqrt(np.mean(err[start:, 2] ** 2))),
        "omega_s_final_rad_s": float(omega_s[-1]),
    }
    return Sled23Result(op=op, t=t, true_state=true_state, est_state=est_state, omega_s=omega_s, stable=stable, metrics=metrics)


def save_summary(results: list[Sled23Result]) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / "sled23_appendix_observer_summary.csv"
    keys = list(results[0].metrics.keys())
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["operating_point", "stable", *keys])
        for result in results:
            writer.writerow([result.op.name, result.stable, *[result.metrics[k] for k in keys]])
    return path


def plot_results(results: list[Sled23Result]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(len(results), 2, figsize=(12, 8), sharex=False)
    for row, result in enumerate(results):
        t_ms = result.t * 1000.0
        err = result.est_state - result.true_state
        axes[row, 0].plot(t_ms, err[:, 0], label="id error")
        axes[row, 0].plot(t_ms, err[:, 1], label="iq error")
        axes[row, 0].grid(True, alpha=0.35)
        axes[row, 0].set_ylabel(f"{result.op.name}\ncurrent [A]")
        axes[row, 0].legend(loc="upper right", fontsize=8)
        axes[row, 1].plot(t_ms, err[:, 2], label="psi_R error")
        axes[row, 1].plot(t_ms, result.omega_s - result.omega_s[-1], label="omega_s error")
        axes[row, 1].grid(True, alpha=0.35)
        axes[row, 1].set_ylabel("flux [Wb], speed [rad/s]")
        axes[row, 1].legend(loc="upper right", fontsize=8)
    axes[-1, 0].set_xlabel("time [ms]")
    axes[-1, 1].set_xlabel("time [ms]")
    fig.suptitle("SLED 2023 Appendix A observer convergence")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    path = OUT_DIR / "sled23_appendix_observer_convergence.png"
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return path


def main() -> None:
    params = MotorParams()
    results = [simulate_case(params, op) for op in base_operating_points(params)]
    summary = save_summary(results)
    figure = plot_results(results)
    for result in results:
        print(result.op.name, result.stable, result.metrics)
    print(f"Saved summary: {summary}")
    print(f"Saved figure:  {figure}")


if __name__ == "__main__":
    main()
