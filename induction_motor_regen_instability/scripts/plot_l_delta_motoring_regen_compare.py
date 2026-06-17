from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import make_block_frequency_response as base


OUT = Path("outputs/l_delta_motoring_regen_compare.png")
OUT_TXT = Path("outputs/l_delta_motoring_regen_compare.txt")


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def zero_db_crossings(f: np.ndarray, loop: np.ndarray) -> list[tuple[float, float]]:
    mag_db = 20.0 * np.log10(np.maximum(np.abs(loop), 1e-30))
    ph = phase_deg(loop)
    out: list[tuple[float, float]] = []
    for i in range(len(f) - 1):
        if (mag_db[i] == 0.0) or (mag_db[i] * mag_db[i + 1] < 0.0):
            a = abs(mag_db[i])
            b = abs(mag_db[i + 1])
            r = a / (a + b) if a + b > 0 else 0.0
            fc = f[i] + (f[i + 1] - f[i]) * r
            pc = ph[i] + (ph[i + 1] - ph[i]) * r
            out.append((fc, pc))
    return out


def main() -> None:
    f = np.logspace(-1, 4, 1800)
    w = 2.0 * math.pi * f
    l_mot = base.delta_loop(base.ID0, w)
    l_reg = base.delta_loop(-base.ID0, w)

    mot_color = "#0072B2"
    reg_color = "#D55E00"

    fig, axes = plt.subplots(2, 2, figsize=(15, 9))
    fig.suptitle(
        "Axis-error to axis-error loop L_delta: motoring vs regeneration",
        fontsize=15,
        fontweight="bold",
    )

    for label, loop, color, ls in [
        ("motoring +220 Nm", l_mot, mot_color, "-"),
        ("regen -220 Nm", l_reg, reg_color, "--"),
    ]:
        axes[0, 0].semilogx(
            f,
            20.0 * np.log10(np.maximum(np.abs(loop), 1e-30)),
            color=color,
            ls=ls,
            lw=2.2,
            label=label,
        )
        axes[1, 0].semilogx(f, phase_deg(loop), color=color, ls=ls, lw=2.2, label=label)
        axes[0, 1].plot(loop.real, loop.imag, color=color, ls=ls, lw=2.0, label=label)
        axes[1, 1].plot(loop.real, loop.imag, color=color, ls=ls, lw=2.0, label=label)

    axes[0, 0].axhline(0, color="black", lw=0.8)
    axes[0, 0].set_title("Bode gain")
    axes[0, 0].set_ylabel("|L_delta| [dB]")
    axes[0, 0].legend(fontsize=9)

    axes[1, 0].axhline(-180, color="#64748b", lw=0.9, ls=":")
    axes[1, 0].axhline(0, color="black", lw=0.8)
    axes[1, 0].set_title("Bode phase")
    axes[1, 0].set_xlabel("frequency [Hz]")
    axes[1, 0].set_ylabel("phase [deg]")
    axes[1, 0].legend(fontsize=9)

    for ax, title, xlim, ylim in [
        (axes[0, 1], "Nyquist around -1", (-2.5, 1.0), (-1.75, 1.75)),
        (axes[1, 1], "Nyquist medium range", (-12.0, 3.0), (-7.5, 7.5)),
    ]:
        ax.plot([-1], [0], marker="x", color="#cc3311", ms=10, mew=2.2, label="-1 + j0")
        ax.axhline(0, color="black", lw=0.8)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(title)
        ax.set_xlabel("Re{L_delta(jw)}")
        ax.set_ylabel("Im{L_delta(jw)}")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal", adjustable="box")
        ax.legend(fontsize=8)

    lines = []
    for label, loop in [("motoring +220 Nm", l_mot), ("regen -220 Nm", l_reg)]:
        crossings = zero_db_crossings(f, loop)
        if crossings:
            for fc, ph in crossings:
                lines.append(f"{label}: 0 dB crossing at {fc:.3f} Hz, phase={ph:.1f} deg")
        else:
            lines.append(f"{label}: no 0 dB crossing")
        idx = int(np.argmin(np.abs(loop + 1.0)))
        lines.append(
            f"{label}: nearest to -1 at {f[idx]:.3f} Hz, "
            f"L={loop[idx].real:+.3f}{loop[idx].imag:+.3f}j, "
            f"distance={abs(loop[idx] + 1.0):.3f}"
        )

    axes[0, 0].text(
        0.02,
        0.05,
        "\n".join(lines),
        transform=axes[0, 0].transAxes,
        fontsize=8.5,
        family="monospace",
        va="bottom",
        bbox=dict(facecolor="white", edgecolor="#cbd5e1", alpha=0.85),
    )

    for ax in axes.ravel():
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
