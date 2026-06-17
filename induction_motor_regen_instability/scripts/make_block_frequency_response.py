from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUT = Path("outputs/block_frequency_response_mechanism.png")
OUT_NYQUIST = Path("outputs/l_delta_nyquist.png")
OUT_CURRENT_CLOSED = Path("outputs/l_delta_current_loop_closure.png")


RS = 0.00762
RR = 0.008041
LLS = 0.0000419
LLR = 0.0000419
LM = 0.0001583
LS = LLS + LM
LR = LLR + LM
POLE_PAIRS = 4
SPEED_RPM = 5000.0
OMEGA_M = SPEED_RPM * 2.0 * math.pi / 60.0
OMEGA_R = POLE_PAIRS * OMEGA_M
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


def operating_state(iq0: float) -> tuple[np.ndarray, float, float]:
    phi = LM * ID0
    i_rd = 0.0
    i_rq = -LM / LR * iq0
    psi_sd = LS * ID0 + LM * i_rd
    psi_sq = LS * iq0 + LM * i_rq
    psi_rd = phi
    psi_rq = 0.0
    slip = RR * LM * iq0 / (LR * phi)
    omega_e = OMEGA_R + slip
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


def linearize(iq0: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x0, slip, omega_e = operating_state(iq0)
    n = 4
    A = np.zeros((n, n))
    B = np.zeros((n, 2))
    for i in range(n):
        h = 1e-7 * max(1.0, abs(x0[i]))
        xp = x0.copy()
        xm = x0.copy()
        xp[i] += h
        xm[i] -= h
        A[:, i] = (derivatives(xp, np.zeros(2), slip, omega_e) - derivatives(xm, np.zeros(2), slip, omega_e)) / (2.0 * h)
    for j in range(2):
        h = 1e-4
        up = np.zeros(2)
        um = np.zeros(2)
        up[j] = h
        um[j] = -h
        B[:, j] = (derivatives(x0, up, slip, omega_e) - derivatives(x0, um, slip, omega_e)) / (2.0 * h)

    det = LS * LR - LM * LM
    C = np.array(
        [
            [LR / det, 0.0, -LM / det, 0.0],  # id
            [0.0, LR / det, 0.0, -LM / det],  # iq
            [0.0, 0.0, 1.0, 0.0],  # phi_d as flux magnitude perturbation near phi_q=0
        ],
        dtype=float,
    )
    return A, B, C


def gm_response(iq0: float, omega: np.ndarray) -> np.ndarray:
    A, B, C = linearize(iq0)
    eye = np.eye(A.shape[0])
    out = np.zeros((len(omega), 3, 2), dtype=complex)
    for k, w in enumerate(omega):
        out[k] = C @ np.linalg.solve(1j * w * eye - A, B)
    return out


def ci_response(omega: np.ndarray) -> np.ndarray:
    return KP * (1.0 + KI / (1j * omega))


def s_omega(iq0: float) -> np.ndarray:
    phi = LM * ID0
    return np.array([0.0, KSLIP / phi, -KSLIP * iq0 / (phi * phi)], dtype=float)


def delta_loop(iq0: float, omega: np.ndarray) -> np.ndarray:
    gm = gm_response(iq0, omega)
    ci = ci_response(omega)
    m_delta_e = np.array([-iq0, ID0], dtype=float)
    s = s_omega(iq0)
    hdot = np.zeros(len(omega), dtype=complex)
    for k in range(len(omega)):
        hdot[k] = -s @ gm[k] @ (ci[k] * m_delta_e)
    return hdot / (1j * omega)


def delta_loop_with_current_feedback(iq0: float, omega: np.ndarray) -> np.ndarray:
    gm = gm_response(iq0, omega)
    ci = ci_response(omega)
    m_delta_e = np.array([-iq0, ID0], dtype=complex)
    s = s_omega(iq0).astype(complex)
    p_i = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=complex)
    eye2 = np.eye(2, dtype=complex)
    hdot = np.zeros(len(omega), dtype=complex)
    for k in range(len(omega)):
        current_loop = eye2 + p_i @ gm[k] * ci[k]
        e = np.linalg.solve(current_loop, m_delta_e)
        x_m = gm[k] @ (ci[k] * e)
        hdot[k] = -s @ x_m
    return hdot / (1j * omega)


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def add_direction_arrows(ax: plt.Axes, z: np.ndarray, color: str) -> None:
    for frac in (0.08, 0.22, 0.45, 0.70):
        idx = int(frac * (len(z) - 2))
        ax.annotate(
            "",
            xy=(z[idx + 1].real, z[idx + 1].imag),
            xytext=(z[idx].real, z[idx].imag),
            arrowprops=dict(arrowstyle="->", color=color, lw=1.4, shrinkA=0, shrinkB=0),
        )


def add_frequency_markers(ax: plt.Axes, f: np.ndarray, z: np.ndarray, color: str, labels: tuple[float, ...]) -> None:
    for target in labels:
        idx = int(np.argmin(np.abs(f - target)))
        ax.plot(z[idx].real, z[idx].imag, marker="o", ms=4.5, color=color)
        ax.annotate(
            f"{target:g} Hz",
            xy=(z[idx].real, z[idx].imag),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7,
            color=color,
        )


def plot_nyquist(f: np.ndarray, l_mot: np.ndarray, l_reg: np.ndarray) -> None:
    mot_color = "#0072B2"
    reg_color = "#E69F00"
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8))
    fig.suptitle("Nyquist plot of the composed axis-error loop L_delta", fontsize=14, fontweight="bold")

    panels = [
        ("Zoom around -1", (-2.5, 1.0), (-1.75, 1.75), True),
        ("Medium range", (-12.0, 3.0), (-7.5, 7.5), True),
        (f"Full range, f={f[0]:.1f} to {f[-1]:.0f} Hz", None, None, False),
    ]
    for ax, (title, xlim, ylim, equal_aspect) in zip(axes, panels):
        ax.plot(l_mot.real, l_mot.imag, color=mot_color, lw=2.0, label="motoring")
        ax.plot(l_mot.real, -l_mot.imag, color=mot_color, lw=1.0, alpha=0.35)
        ax.plot(l_reg.real, l_reg.imag, color=reg_color, lw=2.0, ls="--", label="regen")
        ax.plot(l_reg.real, -l_reg.imag, color=reg_color, lw=1.0, ls="--", alpha=0.35)
        add_direction_arrows(ax, l_mot, mot_color)
        add_direction_arrows(ax, l_reg, reg_color)
        add_frequency_markers(ax, f, l_mot, mot_color, (0.1, 1.0, 3.0, 10.0))
        add_frequency_markers(ax, f, l_reg, reg_color, (0.1, 1.0, 3.0, 10.0))
        ax.plot([-1], [0], marker="x", ms=10, mew=2.2, color="#d62728", label="-1 + j0")
        ax.axhline(0, color="black", lw=0.8)
        ax.axvline(0, color="black", lw=0.8)
        ax.grid(True, color="#cbd5e1", alpha=0.85)
        ax.set_xlabel("Re{L_delta(jw)}")
        ax.set_ylabel("Im{L_delta(jw)}")
        if equal_aspect:
            ax.set_aspect("equal", adjustable="box")
        ax.legend(fontsize=8)
        ax.set_title(title)
        if xlim is not None and ylim is not None:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        else:
            real_min = min(np.min(l_mot.real), np.min(l_reg.real), -1.0)
            real_max = max(np.max(l_mot.real), np.max(l_reg.real), 1.0)
            imag_abs = max(np.max(np.abs(l_mot.imag)), np.max(np.abs(l_reg.imag)), 1.0)
            pad = 0.08 * max(real_max - real_min, 2.0)
            ax.set_xlim(real_min - pad, real_max + pad)
            ax.set_ylim(-imag_abs * 1.08, imag_abs * 1.08)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    OUT_NYQUIST.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_NYQUIST, dpi=170)
    print(OUT_NYQUIST)


def plot_current_loop_closure(f: np.ndarray, l_mot_open: np.ndarray, l_reg_open: np.ndarray, l_mot_cl: np.ndarray, l_reg_cl: np.ndarray) -> None:
    mot_color = "#0072B2"
    reg_color = "#E69F00"
    fig, axes = plt.subplots(2, 2, figsize=(15, 8), sharex="col")
    fig.suptitle("Effect of closing the d/q current loop inside L_delta", fontsize=14, fontweight="bold")

    series = [
        ("motoring open", l_mot_open, mot_color, "-"),
        ("motoring current-closed", l_mot_cl, mot_color, ":"),
        ("regen open", l_reg_open, reg_color, "--"),
        ("regen current-closed", l_reg_cl, reg_color, "-."),
    ]
    for label, loop, color, ls in series:
        axes[0, 0].semilogx(f, 20 * np.log10(np.maximum(np.abs(loop), 1e-30)), color=color, ls=ls, lw=2.0, label=label)
        axes[1, 0].semilogx(f, phase_deg(loop), color=color, ls=ls, lw=2.0, label=label)
    axes[0, 0].axhline(0, color="black", lw=0.8)
    axes[0, 0].set_title("Bode gain")
    axes[0, 0].set_ylabel("|L_delta| [dB]")
    axes[1, 0].set_title("Bode phase")
    axes[1, 0].set_ylabel("phase [deg]")
    axes[1, 0].set_xlabel("frequency [Hz]")

    for label, loop, color, ls in series:
        axes[0, 1].plot(loop.real, loop.imag, color=color, ls=ls, lw=2.0, label=label)
        axes[1, 1].plot(loop.real, loop.imag, color=color, ls=ls, lw=2.0, label=label)
    for ax, title, xlim, ylim in [
        (axes[0, 1], "Nyquist around -1", (-2.5, 1.0), (-1.75, 1.75)),
        (axes[1, 1], "Nyquist medium range", (-12.0, 3.0), (-7.5, 7.5)),
    ]:
        ax.plot([-1], [0], marker="x", ms=10, mew=2.2, color="#d62728", label="-1 + j0")
        ax.axhline(0, color="black", lw=0.8)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("Re{L_delta(jw)}")
        ax.set_ylabel("Im{L_delta(jw)}")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")

    for ax in axes.ravel():
        ax.grid(True, which="both", color="#cbd5e1", alpha=0.85)
        ax.legend(fontsize=7)
        if ax.get_xscale() == "log":
            ax.set_xlim(f[0], f[-1])
            ax.tick_params(axis="x", which="both", labelbottom=True)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    OUT_CURRENT_CLOSED.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_CURRENT_CLOSED, dpi=170)
    print(OUT_CURRENT_CLOSED)


def main() -> None:
    f = np.logspace(-1, 4, 1000)
    w = 2.0 * math.pi * f
    ci = ci_response(w) / KP
    gm_mot = gm_response(ID0, w)
    gm_reg = gm_response(-ID0, w)
    l_mot = delta_loop(ID0, w)
    l_reg = delta_loop(-ID0, w)
    l_mot_cl = delta_loop_with_current_feedback(ID0, w)
    l_reg_cl = delta_loop_with_current_feedback(-ID0, w)

    fig, axes = plt.subplots(4, 2, figsize=(15, 13), sharex="col")
    fig.suptitle("Small-signal frequency responses of the axis-error loop blocks", fontsize=15, fontweight="bold")

    axes[0, 0].semilogx(f, 20 * np.log10(np.abs(ci)), color="#1f6feb", lw=2.2)
    axes[0, 0].set_ylabel("|Ci/Kp| [dB]")
    axes[0, 0].set_title("Current PI block vs frequency")
    axes[1, 0].semilogx(f, phase_deg(ci), color="#1f6feb", lw=2.2)
    axes[1, 0].set_ylabel("phase [deg]")
    axes[1, 0].set_xlabel("frequency [Hz]")
    axes[1, 0].axhline(-90, color="#94a3b8", ls=":", lw=1.0)
    axes[1, 0].axhline(0, color="black", lw=0.8)

    mot_color = "#0072B2"  # Okabe-Ito blue
    reg_color = "#E69F00"  # Okabe-Ito orange
    d_color = "#0072B2"
    q_color = "#CC79A7"

    for label, gm, color, ls in [
        ("id/vd", gm_mot[:, 0, 0], d_color, "-"),
        ("iq/vq", gm_mot[:, 1, 1], q_color, "--"),
    ]:
        axes[0, 1].semilogx(f, 20 * np.log10(np.maximum(np.abs(gm), 1e-30)), color=color, ls=ls, label=label)
        axes[1, 1].semilogx(f, phase_deg(gm), color=color, ls=ls, label=label)
    axes[0, 1].set_title("Gm motoring: voltage -> current")
    axes[0, 1].set_ylabel("gain [A/V dB]")
    axes[1, 1].set_ylabel("phase [deg]")
    axes[1, 1].set_xlabel("frequency [Hz]")
    axes[0, 1].legend(fontsize=8)
    axes[1, 1].legend(fontsize=8)

    for label, gm, color, ls in [
        ("id/vd", gm_reg[:, 0, 0], reg_color, "-"),
        ("iq/vq", gm_reg[:, 1, 1], "#D55E00", "--"),
    ]:
        axes[2, 0].semilogx(f, 20 * np.log10(np.maximum(np.abs(gm), 1e-30)), color=color, ls=ls, label=label)
        axes[3, 0].semilogx(f, phase_deg(gm), color=color, ls=ls, label=label)
    axes[2, 0].set_title("Gm regeneration: voltage -> current")
    axes[2, 0].set_ylabel("gain [A/V dB]")
    axes[3, 0].set_ylabel("phase [deg]")
    axes[3, 0].set_xlabel("frequency [Hz]")
    axes[2, 0].legend(fontsize=8)
    axes[3, 0].legend(fontsize=8)

    for label, loop, color, ls in [("motoring", l_mot, mot_color, "-"), ("regen", l_reg, reg_color, "--")]:
        axes[2, 1].semilogx(f, 20 * np.log10(np.maximum(np.abs(loop), 1e-30)), color=color, ls=ls, lw=2.2, label=label)
        axes[3, 1].semilogx(f, phase_deg(loop), color=color, ls=ls, lw=2.2, label=label)
    axes[2, 1].axhline(0, color="black", lw=0.8)
    axes[2, 1].set_title("Composed axis-error loop L_delta")
    axes[2, 1].set_ylabel("|L_delta| [dB]")
    axes[3, 1].set_ylabel("phase [deg]")
    axes[3, 1].set_xlabel("frequency [Hz]")
    axes[2, 1].legend()
    axes[3, 1].legend()

    for ax in axes.ravel():
        if ax.has_data():
            ax.grid(True, which="both", color="#cbd5e1", alpha=0.85)
            ax.set_xlim(f[0], f[-1])
            ax.tick_params(axis="x", which="both", labelbottom=True)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170)
    print(OUT)
    plot_nyquist(f, l_mot, l_reg)
    plot_current_loop_closure(f, l_mot, l_reg, l_mot_cl, l_reg_cl)


if __name__ == "__main__":
    main()
