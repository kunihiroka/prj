from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT_ROOTS = Path("outputs/full_closed_loop_eigenvalues.png")
OUT_SWEEP = Path("outputs/full_closed_loop_speed_torque_sweep.png")
OUT_TEXT = Path("outputs/full_closed_loop_eigenvalues.txt")

RS = 0.00762
RR = 0.008041
LLS = 0.0000419
LLR = 0.0000419
LM = 0.0001583
LS = LLS + LM
LR = LLR + LM
POLE_PAIRS = 4
SIGMA = 1.0 - LM * LM / (LS * LR)
SIGMA_LS = SIGMA * LS
KP = SIGMA_LS * 1000.0
KI = RS / SIGMA_LS
KSLIP = RR * LM / LR
TR = LR / RR
ID_BASE = 550.0 * math.sqrt(3.0) / math.sqrt(2.0)


def currents_from_flux(x: np.ndarray) -> np.ndarray:
    psid, psiq, phid, phiq = x[:4]
    det = LS * LR - LM * LM
    ids = (LR * psid - LM * phid) / det
    iqs = (LR * psiq - LM * phiq) / det
    idr = (-LM * psid + LS * phid) / det
    iqr = (-LM * psiq + LS * phiq) / det
    return np.array([ids, iqs, idr, iqr], dtype=float)


def operating_flux(id_ref: float, iq_ref: float) -> np.ndarray:
    phi = LM * id_ref
    i_rd = 0.0
    i_rq = -LM / LR * iq_ref
    psi_sd = LS * id_ref + LM * i_rd
    psi_sq = LS * iq_ref + LM * i_rq
    return np.array([psi_sd, psi_sq, phi, 0.0], dtype=float)


def slip_command(flux: np.ndarray, phihat: float, iq_ref: float, flux_source: str) -> float:
    if flux_source == "id_ref":
        phi_for_slip = phihat
    elif flux_source == "id_feedback":
        phi_for_slip = phihat
    elif flux_source == "actual_flux":
        phi_for_slip = flux[2]
    else:
        raise ValueError(flux_source)
    return KSLIP * iq_ref / phi_for_slip


def required_voltage(flux: np.ndarray, speed_rpm: float, iq_ref: float, flux_source: str) -> tuple[float, float, float, float]:
    ids, iqs, _, _ = currents_from_flux(flux)
    phi_hat = LM * ids
    slip = slip_command(flux, phi_hat, iq_ref, flux_source)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip
    psid, psiq, _, _ = flux
    vd = RS * ids - omega_e * psiq
    vq = RS * iqs + omega_e * psid
    return vd, vq, slip, omega_e


def decoupling_voltage(
    flux: np.ndarray,
    phihat: float,
    speed_rpm: float,
    id_ref: float,
    iq_ref: float,
    mode: str,
    flux_source: str,
) -> np.ndarray:
    if mode == "none":
        return np.zeros(2)
    slip_cmd = slip_command(flux, phihat, iq_ref, flux_source)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip_cmd
    if mode == "command":
        phi_dec = LM * id_ref
    elif mode == "flux_feedback":
        phi_dec = flux[2]
    else:
        raise ValueError(mode)
    vd_ff = -omega_e * SIGMA_LS * iq_ref
    vq_ff = omega_e * (SIGMA_LS * id_ref + (LM / LR) * phi_dec)
    return np.array([vd_ff, vq_ff], dtype=float)


def initial_state(speed_rpm: float, iq_ref: float, ff_mode: str, flux_source: str) -> np.ndarray:
    id_ref = ID_BASE
    flux = operating_flux(id_ref, iq_ref)
    phihat = LM * id_ref
    vd_req, vq_req, _, _ = required_voltage(flux, speed_rpm, iq_ref, flux_source)
    ff = decoupling_voltage(flux, phihat, speed_rpm, id_ref, iq_ref, ff_mode, flux_source)
    z = np.array([vd_req, vq_req]) - ff
    return np.r_[flux, z, phihat]


def rhs(x: np.ndarray, speed_rpm: float, iq_ref: float, ff_mode: str, flux_source: str) -> np.ndarray:
    id_ref = ID_BASE
    flux = x[:4]
    z = x[4:6]
    phihat = x[6]
    ids, iqs, idr, iqr = currents_from_flux(flux)
    e = np.array([id_ref - ids, iq_ref - iqs], dtype=float)
    ff = decoupling_voltage(flux, phihat, speed_rpm, id_ref, iq_ref, ff_mode, flux_source)
    v = KP * e + z + ff
    slip_cmd = slip_command(flux, phihat, iq_ref, flux_source)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip_cmd
    psid, psiq, phid, phiq = flux
    dflux = np.array(
        [
            v[0] - RS * ids + omega_e * psiq,
            v[1] - RS * iqs - omega_e * psid,
            -RR * idr + slip_cmd * phiq,
            -RR * iqr - slip_cmd * phid,
        ],
        dtype=float,
    )
    dz = KP * KI * e
    if flux_source == "id_feedback":
        phihat_input = LM * ids
    else:
        phihat_input = LM * id_ref
    dphihat = (phihat_input - phihat) / TR
    return np.r_[dflux, dz, dphihat]


def jacobian(x0: np.ndarray, speed_rpm: float, iq_ref: float, ff_mode: str, flux_source: str) -> np.ndarray:
    n = len(x0)
    a = np.zeros((n, n))
    for k in range(n):
        h = 1e-6 * max(1.0, abs(x0[k]))
        xp = x0.copy()
        xm = x0.copy()
        xp[k] += h
        xm[k] -= h
        a[:, k] = (rhs(xp, speed_rpm, iq_ref, ff_mode, flux_source) - rhs(xm, speed_rpm, iq_ref, ff_mode, flux_source)) / (2.0 * h)
    return a


def eigen_summary(speed_rpm: float, iq_ref: float, ff_mode: str, flux_source: str = "id_ref") -> tuple[np.ndarray, np.ndarray]:
    x0 = initial_state(speed_rpm, iq_ref, ff_mode, flux_source)
    residual = rhs(x0, speed_rpm, iq_ref, ff_mode, flux_source)
    a = jacobian(x0, speed_rpm, iq_ref, ff_mode, flux_source)
    return np.linalg.eigvals(a), residual


def format_eigs(eigs: np.ndarray) -> list[str]:
    lines = []
    for lam in sorted(eigs, key=lambda z: z.real, reverse=True):
        freq = abs(lam.imag) / (2.0 * math.pi)
        damping = -lam.real / max(abs(lam), 1e-30)
        lines.append(f"{lam.real:+10.4f} {lam.imag:+10.4f}j   f={freq:8.3f} Hz   zeta={damping:+7.4f}")
    return lines


def plot_roots(results: dict[str, np.ndarray]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    colors = {
        "motoring / no decoupling": "#0072B2",
        "regen / no decoupling": "#E69F00",
        "motoring / flux-feedback decoupling": "#56B4E9",
        "regen / flux-feedback decoupling": "#D55E00",
    }
    for ax in axes:
        for label, eigs in results.items():
            ax.scatter(eigs.real, eigs.imag / (2.0 * math.pi), s=44, label=label, color=colors.get(label), alpha=0.9)
        ax.axvline(0, color="black", lw=0.9)
        ax.axhline(0, color="black", lw=0.7)
        ax.grid(True, color="#cbd5e1", alpha=0.85)
        ax.set_xlabel("real part [1/s]")
        ax.set_ylabel("imaginary frequency [Hz]")
    axes[0].set_title("All eigenvalues")
    axes[1].set_title("Low-frequency modes")
    axes[1].set_xlim(-130, 15)
    axes[1].set_ylim(-12, 12)
    axes[1].legend(fontsize=8, loc="lower right")
    fig.suptitle("Full closed-loop linearized eigenvalues at 5000 r/min, |torque|=220 Nm equivalent")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    OUT_ROOTS.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_ROOTS, dpi=170)
    print(OUT_ROOTS)


def sweep_max_real(ff_mode: str, flux_source: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    speeds = np.array([1000, 2000, 3000, 4000, 5000, 6000], dtype=float)
    torque_fracs = np.array([0.25, 0.5, 0.75, 1.0], dtype=float)
    mot = np.zeros((len(torque_fracs), len(speeds)))
    reg = np.zeros_like(mot)
    for i, frac in enumerate(torque_fracs):
        for j, speed in enumerate(speeds):
            iq = ID_BASE * frac
            mot[i, j] = max(eigen_summary(speed, iq, ff_mode, flux_source)[0].real)
            reg[i, j] = max(eigen_summary(speed, -iq, ff_mode, flux_source)[0].real)
    return speeds, torque_fracs, mot, reg


def plot_sweep() -> None:
    speeds, fracs, mot, reg = sweep_max_real("flux_feedback", "id_ref")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True, constrained_layout=True)
    for ax, data, title in [(axes[0], mot, "motoring"), (axes[1], reg, "regen")]:
        im = ax.imshow(data, origin="lower", aspect="auto", cmap="coolwarm", vmin=-80, vmax=80)
        ax.set_title(f"max Re(lambda): {title}")
        ax.set_xticks(np.arange(len(speeds)), [f"{int(s)}" for s in speeds])
        ax.set_yticks(np.arange(len(fracs)), [f"{f:.2f}" for f in fracs])
        ax.set_xlabel("speed [r/min]")
        ax.set_ylabel("|torque| / rated")
        for i in range(len(fracs)):
            for j in range(len(speeds)):
                ax.text(j, i, f"{data[i,j]:.1f}", ha="center", va="center", fontsize=8)
    fig.colorbar(im, ax=axes.ravel().tolist(), label="max real eigenvalue [1/s]", shrink=0.86)
    fig.suptitle("Closed-loop eigenvalue sweep with flux-feedback decoupling")
    OUT_SWEEP.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SWEEP, dpi=170)
    print(OUT_SWEEP)


def main() -> None:
    cases = {
        "motoring / no decoupling / id_ref slip": (5000.0, ID_BASE, "none", "id_ref"),
        "regen / no decoupling / id_ref slip": (5000.0, -ID_BASE, "none", "id_ref"),
        "regen / no decoupling / id_feedback slip": (5000.0, -ID_BASE, "none", "id_feedback"),
        "regen / no decoupling / actual-flux slip": (5000.0, -ID_BASE, "none", "actual_flux"),
        "motoring / flux-feedback decoupling / id_ref slip": (5000.0, ID_BASE, "flux_feedback", "id_ref"),
        "regen / flux-feedback decoupling / id_ref slip": (5000.0, -ID_BASE, "flux_feedback", "id_ref"),
        "regen / flux-feedback decoupling / id_feedback slip": (5000.0, -ID_BASE, "flux_feedback", "id_feedback"),
        "regen / flux-feedback decoupling / actual-flux slip": (5000.0, -ID_BASE, "flux_feedback", "actual_flux"),
    }
    results: dict[str, np.ndarray] = {}
    text_lines = []
    for label, (speed, iq, ff_mode, flux_source) in cases.items():
        eigs, residual = eigen_summary(speed, iq, ff_mode, flux_source)
        results[label] = eigs
        text_lines.append(f"\n[{label}]")
        text_lines.append(f"max residual = {np.max(np.abs(residual)):.6e}")
        text_lines.append(f"max real eigenvalue = {np.max(eigs.real):+.6f} 1/s")
        text_lines.extend(format_eigs(eigs))
    OUT_TEXT.parent.mkdir(parents=True, exist_ok=True)
    OUT_TEXT.write_text("\n".join(text_lines), encoding="utf-8")
    print(OUT_TEXT)
    print("\n".join(text_lines))
    plot_roots(results)
    plot_sweep()


if __name__ == "__main__":
    main()
