from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import make_block_frequency_response as base


OUT = Path("outputs/phi_contribution_bode.png")
OUT_TXT = Path("outputs/phi_contribution_bode.txt")


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def contribution(iq0: float, omega: np.ndarray) -> dict[str, np.ndarray]:
    gm = base.gm_response(iq0, omega)
    ci = base.ci_response(omega)
    # Unit axis error: Δe = M_delta_e * Δδ, Δδ=1 rad.
    m_delta_e = np.array([-iq0, base.ID0], dtype=complex)
    vd = ci * m_delta_e[0]
    vq = ci * m_delta_e[1]
    phi_from_vd = gm[:, 2, 0] * vd
    phi_from_vq = gm[:, 2, 1] * vq
    phi_total = phi_from_vd + phi_from_vq
    return {
        "vd": vd,
        "vq": vq,
        "phi_from_vd": phi_from_vd,
        "phi_from_vq": phi_from_vq,
        "phi_total": phi_total,
    }


def bode(ax_mag, ax_phase, f, y, label, color, ls="-", lw=1.9):
    ax_mag.semilogx(f, 20.0 * np.log10(np.maximum(np.abs(y), 1e-30)), color=color, ls=ls, lw=lw, label=label)
    ax_phase.semilogx(f, phase_deg(y), color=color, ls=ls, lw=lw, label=label)


def at_freq_summary(f: np.ndarray, data: dict[str, np.ndarray], target: float, prefix: str) -> list[str]:
    idx = int(np.argmin(np.abs(f - target)))

    def ph(z: complex) -> float:
        return math.degrees(math.atan2(z.imag, z.real))

    return [
        f"{prefix} at {f[idx]:.3f} Hz:",
        f"  Δphi via vd = {abs(data['phi_from_vd'][idx]):.3g} Wb/rad, phase={ph(data['phi_from_vd'][idx]):+.1f} deg",
        f"  Δphi via vq = {abs(data['phi_from_vq'][idx]):.3g} Wb/rad, phase={ph(data['phi_from_vq'][idx]):+.1f} deg",
        f"  Δphi total  = {abs(data['phi_total'][idx]):.3g} Wb/rad, phase={ph(data['phi_total'][idx]):+.1f} deg",
    ]


def vector_panel(ax, data: dict[str, np.ndarray], idx: int, title: str):
    a = data["phi_from_vd"][idx]
    b = data["phi_from_vq"][idx]
    total = data["phi_total"][idx]
    ax.arrow(0, 0, a.real, a.imag, length_includes_head=True, head_width=0.003, color="#0072B2", lw=2.0)
    ax.arrow(a.real, a.imag, b.real, b.imag, length_includes_head=True, head_width=0.003, color="#D55E00", lw=2.0, linestyle="--")
    ax.arrow(0, 0, total.real, total.imag, length_includes_head=True, head_width=0.003, color="#000000", lw=2.4)
    ax.text(a.real, a.imag, " via vd", color="#0072B2", fontsize=9)
    ax.text((a + b).real, (a + b).imag, " total", color="#000000", fontsize=9)
    ax.axhline(0, color="black", lw=0.8)
    ax.axvline(0, color="black", lw=0.8)
    ax.set_title(title)
    ax.set_xlabel("Re{Δphi_d} [Wb/rad]")
    ax.set_ylabel("Im{Δphi_d} [Wb/rad]")
    pts = np.array([0 + 0j, a, a + b, total])
    xmin, xmax = pts.real.min(), pts.real.max()
    ymin, ymax = pts.imag.min(), pts.imag.max()
    pad_x = max(0.01, 0.2 * (xmax - xmin + 1e-12))
    pad_y = max(0.01, 0.2 * (ymax - ymin + 1e-12))
    ax.set_xlim(xmin - pad_x, xmax + pad_x)
    ax.set_ylim(ymin - pad_y, ymax + pad_y)
    ax.set_aspect("equal", adjustable="box")


def main() -> None:
    f = np.logspace(-1, 4, 1800)
    w = 2.0 * math.pi * f
    mot = contribution(base.ID0, w)
    reg = contribution(-base.ID0, w)

    fig, axes = plt.subplots(3, 2, figsize=(15, 11))
    fig.suptitle("Δphi_d contributions in the axis-error loop: Gφd·Δvd and Gφq·Δvq", fontsize=15, fontweight="bold")

    bode(axes[0, 0], axes[1, 0], f, mot["phi_from_vd"], "via vd: Gφd·Δvd", "#0072B2", "-")
    bode(axes[0, 0], axes[1, 0], f, mot["phi_from_vq"], "via vq: Gφq·Δvq", "#D55E00", "--")
    bode(axes[0, 0], axes[1, 0], f, mot["phi_total"], "total Δphi_d", "#000000", "-", lw=2.4)
    axes[0, 0].set_title("Motoring +220 Nm: Δδ -> Δphi_d decomposition")
    axes[0, 0].set_ylabel("gain [Wb/rad dB]")
    axes[1, 0].set_ylabel("phase [deg]")
    axes[1, 0].set_xlabel("frequency [Hz]")

    bode(axes[0, 1], axes[1, 1], f, reg["phi_from_vd"], "via vd: Gφd·Δvd", "#0072B2", "-")
    bode(axes[0, 1], axes[1, 1], f, reg["phi_from_vq"], "via vq: Gφq·Δvq", "#D55E00", "--")
    bode(axes[0, 1], axes[1, 1], f, reg["phi_total"], "total Δphi_d", "#000000", "-", lw=2.4)
    axes[0, 1].set_title("Regeneration -220 Nm: Δδ -> Δphi_d decomposition")
    axes[0, 1].set_ylabel("gain [Wb/rad dB]")
    axes[1, 1].set_ylabel("phase [deg]")
    axes[1, 1].set_xlabel("frequency [Hz]")

    target = 3.42
    idx = int(np.argmin(np.abs(f - target)))
    vector_panel(axes[2, 0], mot, idx, f"Motoring vector sum at {f[idx]:.2f} Hz")
    vector_panel(axes[2, 1], reg, idx, f"Regen vector sum at {f[idx]:.2f} Hz")

    for ax in axes.ravel():
        ax.grid(True, which="both", color="#cbd5e1", alpha=0.85)
        handles, labels = ax.get_legend_handles_labels()
        if handles:
            ax.legend(fontsize=8)
    for ax in axes[:2, :].ravel():
        ax.set_xlim(f[0], f[-1])
        ax.tick_params(axis="x", which="both", labelbottom=True)
    for ax in axes[0, :]:
        ax.axhline(0, color="black", lw=0.8)
    for ax in axes[1, :]:
        ax.axhline(0, color="black", lw=0.8)
        ax.axhline(-180, color="#94a3b8", lw=0.9, ls=":")
        ax.axhline(-90, color="#94a3b8", lw=0.9, ls=":")

    lines = at_freq_summary(f, mot, target, "motoring") + [""] + at_freq_summary(f, reg, target, "regen")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(OUT, dpi=170)
    OUT_TXT.write_text("\n".join(lines), encoding="utf-8")
    print(OUT)
    print(OUT_TXT)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
