from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import plot_full_current_loop_return_ratio as rr


OUT = Path("outputs/full_current_loop_motoring_regen_compare.png")
OUT_TXT = Path("outputs/full_current_loop_motoring_regen_compare.txt")


def phase_deg(z: np.ndarray) -> np.ndarray:
    return np.unwrap(np.angle(z)) * 180.0 / math.pi


def main() -> None:
    f = np.logspace(-1, 4, 1600)
    w = 2.0 * math.pi * f
    cases = [
        ("motoring +220 Nm", rr.ID0, "#0072B2", "-"),
        ("regen -220 Nm", -rr.ID0, "#D55E00", "--"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(15, 11), sharex="col")
    fig.suptitle(
        "Motoring vs regeneration: MIMO current-loop return ratio, id_ref slip, 5000 r/min",
        fontsize=14,
        fontweight="bold",
    )
    lines: list[str] = []
    for label, iq_ref, color, ls in cases:
        loops = rr.loop_response(iq_ref, "flux_feedback", "id_ref", w)
        sigmas = np.array([np.linalg.svd(l, compute_uv=False) for l in loops])
        det_i_l = np.linalg.det(np.eye(2)[None, :, :] + loops)
        eigvals = np.array([np.linalg.eigvals(l) for l in loops])
        axes[0, 0].semilogx(
            f,
            20.0 * np.log10(np.maximum(sigmas[:, 0], 1e-30)),
            color=color,
            ls=ls,
            lw=2.0,
            label=f"{label} sigma_max",
        )
        axes[0, 0].semilogx(
            f,
            20.0 * np.log10(np.maximum(sigmas[:, -1], 1e-30)),
            color=color,
            ls=ls,
            lw=1.2,
            alpha=0.55,
            label=f"{label} sigma_min",
        )
        axes[1, 0].semilogx(
            f,
            20.0 * np.log10(np.maximum(np.abs(det_i_l), 1e-30)),
            color=color,
            ls=ls,
            lw=2.0,
            label=label,
        )
        axes[2, 0].semilogx(f, phase_deg(det_i_l), color=color, ls=ls, lw=2.0, label=label)
        axes[1, 1].plot(det_i_l.real, det_i_l.imag, color=color, ls=ls, lw=2.0, label=label)
        for mode in range(2):
            axes[0, 1].plot(eigvals[:, mode].real, eigvals[:, mode].imag, color=color, ls=ls, lw=1.7)

        idx = int(np.argmin(np.abs(det_i_l)))
        lines.append(
            f"{label}: min |det(I+L)|={abs(det_i_l[idx]):.4f} "
            f"at {f[idx]:.3f} Hz, det={det_i_l[idx].real:+.4f}{det_i_l[idx].imag:+.4f}j"
        )

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

    axes[0, 1].plot([-1], [0], marker="x", color="#cc3311", ms=10, mew=2.2, label="-1 + j0")
    axes[0, 1].axhline(0, color="black", lw=0.8)
    axes[0, 1].axvline(0, color="black", lw=0.8)
    axes[0, 1].set_title("Eigenloci of L_i")
    axes[0, 1].set_xlabel("Re{eig(L_i)}")
    axes[0, 1].set_ylabel("Im{eig(L_i)}")
    axes[0, 1].set_xlim(-3.0, 3.0)
    axes[0, 1].set_ylim(-3.0, 3.0)
    axes[0, 1].legend(fontsize=8)

    axes[1, 1].plot([0], [0], marker="x", color="#cc3311", ms=10, mew=2.2, label="origin")
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
