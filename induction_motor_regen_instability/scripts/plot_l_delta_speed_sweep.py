from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT = Path("outputs/l_delta_speed_sweep_bode.png")

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
ID0 = 550.0 * math.sqrt(3.0) / math.sqrt(2.0)
KSLIP = RR * LM / LR


def currents_from_flux(x: np.ndarray) -> np.ndarray:
    psid, psiq, phid, phiq = x
    det = LS * LR - LM * LM
    id_s = (LR * psid - LM * phid) / det
    iq_s = (LR * psiq - LM * phiq) / det
    id_r = (-LM * psid + LS * phid) / det
    iq_r = (-LM * psiq + LS * phiq) / det
    return np.array([id_s, iq_s, id_r, iq_r], dtype=float)


def operating_state(speed_rpm: float, iq0: float) -> tuple[np.ndarray, float, float]:
    phi = LM * ID0
    i_rd = 0.0
    i_rq = -LM / LR * iq0
    psi_sd = LS * ID0 + LM * i_rd
    psi_sq = LS * iq0 + LM * i_rq
    psi_rd = phi
    psi_rq = 0.0
    slip = RR * LM * iq0 / (LR * phi)
    omega_r = POLE_PAIRS * speed_rpm * 2.0 * math.pi / 60.0
    omega_e = omega_r + slip
    return np.array([psi_sd, psi_sq, psi_rd, psi_rq], dtype=float), slip, omega_e


def derivatives(x: np.ndarray, u: np.ndarray, slip: float, omega_e: float) -> np.ndarray:
    ids, iqs, idr, iqr = currents_from_flux(x)
    vd, vq = u
    psid, psiq, phid, phiq = x
    return np.array(
        [
            vd - RS * ids + omega_e * psiq,
            vq - RS * iqs - omega_e * psid,
            -RR * idr + slip * phiq,
            -RR * iqr - slip * phid,
        ],
        dtype=float,
    )


def linearize(speed_rpm: float, iq0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x0, slip, omega_e = operating_state(speed_rpm, iq0)
    n = 4
    a = np.zeros((n, n))
    b = np.zeros((n, 2))
    for i in range(n):
        h = 1e-7 * max(1.0, abs(x0[i]))
        xp = x0.copy()
        xm = x0.copy()
        xp[i] += h
        xm[i] -= h
        a[:, i] = (derivatives(xp, np.zeros(2), slip, omega_e) - derivatives(xm, np.zeros(2), slip, omega_e)) / (2.0 * h)
    for j in range(2):
        h = 1e-4
        up = np.zeros(2)
        um = np.zeros(2)
        up[j] = h
        um[j] = -h
        b[:, j] = (derivatives(x0, up, slip, omega_e) - derivatives(x0, um, slip, omega_e)) / (2.0 * h)

    det = LS * LR - LM * LM
    c = np.array(
        [
            [LR / det, 0.0, -LM / det, 0.0],
            [0.0, LR / det, 0.0, -LM / det],
            [0.0, 0.0, 1.0, 0.0],
        ],
        dtype=float,
    )
    return a, b, c


def gm_response(speed_rpm: float, iq0: float, omega: np.ndarray) -> np.ndarray:
    a, b, c = linearize(speed_rpm, iq0)
    eye = np.eye(a.shape[0])
    out = np.zeros((len(omega), 3, 2), dtype=complex)
    for k, w in enumerate(omega):
        out[k] = c @ np.linalg.solve(1j * w * eye - a, b)
    return out


def ci_response(omega: np.ndarray) -> np.ndarray:
    return KP * (1.0 + KI / (1j * omega))


def s_omega(iq0: float) -> np.ndarray:
    phi = LM * ID0
    return np.array([0.0, KSLIP / phi, -KSLIP * iq0 / (phi * phi)], dtype=float)


def l_delta(speed_rpm: float, iq0: float, omega: np.ndarray) -> np.ndarray:
    gm = gm_response(speed_rpm, iq0, omega)
    ci = ci_response(omega)
    m_delta_e = np.array([-iq0, ID0], dtype=float)
    s = s_omega(iq0)
    hdot = np.zeros(len(omega), dtype=complex)
    for k in range(len(omega)):
        hdot[k] = -s @ gm[k] @ (ci[k] * m_delta_e)
    return hdot / (1j * omega)


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def main() -> None:
    speeds = [1000.0, 3000.0, 5000.0, 6000.0]
    f = np.logspace(-1, 4, 1400)
    w = 2.0 * math.pi * f
    colors = ["#0072B2", "#009E73", "#E69F00", "#D55E00"]

    fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex=True)
    fig.suptitle("Speed sweep of reduced axis-error loop L_delta", fontsize=14, fontweight="bold")

    for speed, color in zip(speeds, colors):
        mot = l_delta(speed, ID0, w)
        reg = l_delta(speed, -ID0, w)
        axes[0, 0].semilogx(f, 20 * np.log10(np.maximum(np.abs(mot), 1e-30)), color=color, lw=2.0, label=f"{speed:.0f} r/min")
        axes[1, 0].semilogx(f, phase_deg(mot), color=color, lw=2.0)
        axes[0, 1].semilogx(f, 20 * np.log10(np.maximum(np.abs(reg), 1e-30)), color=color, lw=2.0, label=f"{speed:.0f} r/min")
        axes[1, 1].semilogx(f, phase_deg(reg), color=color, lw=2.0)

    axes[0, 0].set_title("Motoring: gain")
    axes[1, 0].set_title("Motoring: phase")
    axes[0, 1].set_title("Regeneration: gain")
    axes[1, 1].set_title("Regeneration: phase")
    for ax in axes[0, :]:
        ax.axhline(0, color="black", lw=0.8)
        ax.set_ylabel("|L_delta| [dB]")
    for ax in axes[1, :]:
        ax.axhline(-180, color="#64748b", ls=":", lw=1.1)
        ax.axhline(180, color="#64748b", ls=":", lw=1.1)
        ax.set_ylabel("phase [deg]")
        ax.set_xlabel("frequency [Hz]")
    for ax in axes.ravel():
        ax.grid(True, which="both", color="#cbd5e1", alpha=0.85)
        ax.set_xlim(f[0], f[-1])
        ax.tick_params(axis="x", which="both", labelbottom=True)
    axes[0, 0].legend(fontsize=8)
    axes[0, 1].legend(fontsize=8)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170)
    print(OUT)


if __name__ == "__main__":
    main()
