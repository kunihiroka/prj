from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT = Path("outputs/full_current_loop_speed_sweep_idref.png")
OUT_TXT = Path("outputs/full_current_loop_speed_sweep_idref.txt")

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
ID0 = 550.0 * math.sqrt(3.0) / math.sqrt(2.0)


def currents_from_flux(flux: np.ndarray) -> np.ndarray:
    psid, psiq, phid, phiq = flux
    det = LS * LR - LM * LM
    ids = (LR * psid - LM * phid) / det
    iqs = (LR * psiq - LM * phiq) / det
    idr = (-LM * psid + LS * phid) / det
    iqr = (-LM * psiq + LS * phiq) / det
    return np.array([ids, iqs, idr, iqr], dtype=float)


def operating_flux(iq_ref: float) -> np.ndarray:
    phi = LM * ID0
    psi_sd = LS * ID0
    psi_sq = LS * iq_ref + LM * (-LM / LR * iq_ref)
    return np.array([psi_sd, psi_sq, phi, 0.0], dtype=float)


def slip_command(phihat: float, iq_ref: float) -> float:
    return KSLIP * iq_ref / phihat


def decoupling_voltage(flux: np.ndarray, phihat: float, iq_ref: float, speed_rpm: float) -> np.ndarray:
    id_ref = ID0
    slip = slip_command(phihat, iq_ref)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip
    phi_dec = flux[2]
    vd_ff = -omega_e * SIGMA_LS * iq_ref
    vq_ff = omega_e * (SIGMA_LS * id_ref + (LM / LR) * phi_dec)
    return np.array([vd_ff, vq_ff], dtype=float)


def required_voltage(flux: np.ndarray, phihat: float, iq_ref: float, speed_rpm: float) -> np.ndarray:
    ids, iqs, _, _ = currents_from_flux(flux)
    slip = slip_command(phihat, iq_ref)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip
    psid, psiq, _, _ = flux
    return np.array([RS * ids - omega_e * psiq, RS * iqs + omega_e * psid], dtype=float)


def plant_rhs(x: np.ndarray, u_pi: np.ndarray, iq_ref: float, speed_rpm: float) -> np.ndarray:
    flux = x[:4]
    phihat = x[4]
    ids, iqs, idr, iqr = currents_from_flux(flux)
    slip = slip_command(phihat, iq_ref)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip
    v = u_pi + decoupling_voltage(flux, phihat, iq_ref, speed_rpm)
    psid, psiq, phid, phiq = flux
    dflux = np.array(
        [
            v[0] - RS * ids + omega_e * psiq,
            v[1] - RS * iqs - omega_e * psid,
            -RR * idr + slip * phiq,
            -RR * iqr - slip * phid,
        ],
        dtype=float,
    )
    dphihat = (LM * ID0 - phihat) / TR
    return np.r_[dflux, dphihat]


def output_current(x: np.ndarray) -> np.ndarray:
    ids, iqs, _, _ = currents_from_flux(x[:4])
    return np.array([ids, iqs], dtype=float)


def operating_point(iq_ref: float, speed_rpm: float) -> tuple[np.ndarray, np.ndarray]:
    flux = operating_flux(iq_ref)
    phihat = LM * ID0
    x0 = np.r_[flux, phihat]
    u0 = required_voltage(flux, phihat, iq_ref, speed_rpm) - decoupling_voltage(flux, phihat, iq_ref, speed_rpm)
    return x0, u0


def linearize_plant(iq_ref: float, speed_rpm: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x0, u0 = operating_point(iq_ref, speed_rpm)
    n = len(x0)
    a = np.zeros((n, n))
    b = np.zeros((n, 2))
    c = np.zeros((2, n))
    for k in range(n):
        h = 1e-6 * max(1.0, abs(x0[k]))
        xp = x0.copy()
        xm = x0.copy()
        xp[k] += h
        xm[k] -= h
        a[:, k] = (plant_rhs(xp, u0, iq_ref, speed_rpm) - plant_rhs(xm, u0, iq_ref, speed_rpm)) / (2.0 * h)
        c[:, k] = (output_current(xp) - output_current(xm)) / (2.0 * h)
    for k in range(2):
        h = 1e-4
        up = u0.copy()
        um = u0.copy()
        up[k] += h
        um[k] -= h
        b[:, k] = (plant_rhs(x0, up, iq_ref, speed_rpm) - plant_rhs(x0, um, iq_ref, speed_rpm)) / (2.0 * h)
    return a, b, c


def ci_response(w: float) -> complex:
    return KP * (1.0 + KI / (1j * w))


def loop_response(iq_ref: float, speed_rpm: float, omega: np.ndarray) -> np.ndarray:
    a, b, c = linearize_plant(iq_ref, speed_rpm)
    eye = np.eye(a.shape[0])
    loops = np.zeros((len(omega), 2, 2), dtype=complex)
    for i, w in enumerate(omega):
        plant = c @ np.linalg.solve(1j * w * eye - a, b)
        loops[i] = ci_response(w) * plant
    return loops


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def main() -> None:
    speeds = [1000.0, 3000.0, 5000.0, 6000.0]
    colors = ["#0072B2", "#009E73", "#E69F00", "#D55E00"]
    f = np.logspace(-1, 4, 1600)
    w = 2.0 * math.pi * f

    fig, axes = plt.subplots(3, 2, figsize=(15, 11), sharex="col")
    fig.suptitle("Speed sweep of MIMO current-loop return ratio, id_ref slip, regen torque", fontsize=14, fontweight="bold")
    lines = []
    for speed, color in zip(speeds, colors):
        loops = loop_response(-ID0, speed, w)
        det_i_l = np.linalg.det(np.eye(2)[None, :, :] + loops)
        sigmas = np.array([np.linalg.svd(l, compute_uv=False) for l in loops])
        eigvals = np.array([np.linalg.eigvals(l) for l in loops])
        label = f"{speed:.0f} r/min"

        axes[0, 0].semilogx(f, 20 * np.log10(np.maximum(sigmas[:, 0], 1e-30)), color=color, lw=2.0, label=f"{label} sigma_max")
        axes[0, 0].semilogx(f, 20 * np.log10(np.maximum(sigmas[:, -1], 1e-30)), color=color, lw=1.2, alpha=0.55)
        axes[1, 0].semilogx(f, 20 * np.log10(np.maximum(np.abs(det_i_l), 1e-30)), color=color, lw=2.0, label=label)
        axes[2, 0].semilogx(f, phase_deg(det_i_l), color=color, lw=2.0, label=label)
        axes[1, 1].plot(det_i_l.real, det_i_l.imag, color=color, lw=2.0, label=label)
        for mode in range(2):
            axes[0, 1].plot(eigvals[:, mode].real, eigvals[:, mode].imag, color=color, lw=1.5, alpha=0.8)

        dist = np.min(np.abs(det_i_l))
        idx = int(np.argmin(np.abs(det_i_l)))
        lines.append(f"{label}: min |det(I+L)|={dist:.4f} at {f[idx]:.3f} Hz, det={det_i_l[idx].real:+.4f}{det_i_l[idx].imag:+.4f}j")

    axes[0, 0].axhline(0, color="black", lw=0.8)
    axes[0, 0].set_title("Singular values of return ratio L_i")
    axes[0, 0].set_ylabel("gain [dB]")
    axes[0, 0].legend(fontsize=8)

    axes[1, 0].axhline(0, color="black", lw=0.8)
    axes[1, 0].set_title("Bode gain of det(I + L_i)")
    axes[1, 0].set_ylabel("|det(I+L_i)| [dB]")
    axes[1, 0].legend(fontsize=8)

    axes[2, 0].axhline(0, color="black", lw=0.8)
    axes[2, 0].set_title("Bode phase of det(I + L_i)")
    axes[2, 0].set_ylabel("phase [deg]")
    axes[2, 0].set_xlabel("frequency [Hz]")
    axes[2, 0].legend(fontsize=8)

    axes[0, 1].plot([-1], [0], marker="x", color="#d62728", ms=10, mew=2.2, label="-1 + j0")
    axes[0, 1].axhline(0, color="black", lw=0.8)
    axes[0, 1].axvline(0, color="black", lw=0.8)
    axes[0, 1].set_title("Eigenloci of L_i")
    axes[0, 1].set_xlabel("Re{eig(L_i)}")
    axes[0, 1].set_ylabel("Im{eig(L_i)}")
    axes[0, 1].set_xlim(-3.0, 3.0)
    axes[0, 1].set_ylim(-3.0, 3.0)
    axes[0, 1].legend(fontsize=8)

    axes[1, 1].plot([0], [0], marker="x", color="#d62728", ms=10, mew=2.2, label="origin")
    axes[1, 1].axhline(0, color="black", lw=0.8)
    axes[1, 1].axvline(0, color="black", lw=0.8)
    axes[1, 1].set_title("MIMO Nyquist: det(I + L_i)")
    axes[1, 1].set_xlabel("Re{det(I+L_i)}")
    axes[1, 1].set_ylabel("Im{det(I+L_i)}")
    axes[1, 1].legend(fontsize=8)

    axes[2, 1].axis("off")
    axes[2, 1].text(0.02, 0.95, "\n".join(lines), va="top", family="monospace", fontsize=10)

    for ax in [axes[0, 0], axes[1, 0], axes[2, 0], axes[0, 1], axes[1, 1]]:
        ax.grid(True, which="both", color="#cbd5e1", alpha=0.85)
    for ax in axes[:, 0]:
        ax.set_xlim(f[0], f[-1])
        ax.tick_params(axis="x", which="both", labelbottom=True)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170)
    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)
    print(OUT_TXT)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
