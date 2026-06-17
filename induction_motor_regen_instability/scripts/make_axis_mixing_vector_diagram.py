from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc, FancyArrowPatch


OUT = Path("outputs/axis_mixing_motoring_regen_vector_diagram.png")


def arrow(ax, start, end, color, label="", lw=2.0, mutation_scale=13, ls="-"):
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation_scale,
        linewidth=lw,
        color=color,
        linestyle=ls,
    )
    ax.add_patch(patch)
    if label:
        x = 0.5 * (start[0] + end[0])
        y = 0.5 * (start[1] + end[1])
        ax.text(x, y, label, color=color, fontsize=9, ha="center", va="center")


def unit(theta):
    return math.cos(theta), math.sin(theta)


def scale(v, k):
    return v[0] * k, v[1] * k


def add(v, w):
    return v[0] + w[0], v[1] + w[1]


def draw_case(ax, title, iq_sign, delta_sign):
    delta = delta_sign * math.radians(16)
    id_ctrl = 1.15
    iq_ctrl = iq_sign * 0.95

    d_real = (1.0, 0.0)
    q_real = (0.0, 1.0)
    d_ctrl = unit(delta)
    q_ctrl = unit(delta + math.pi / 2.0)

    i_id = scale(d_ctrl, id_ctrl)
    i_iq = scale(q_ctrl, iq_ctrl)
    i_vec = add(i_id, i_iq)
    iq_to_d_real = iq_ctrl * q_ctrl[0]
    approx = -delta * iq_ctrl

    ax.axhline(0, color="#94a3b8", lw=1.0)
    ax.axvline(0, color="#94a3b8", lw=1.0)
    arrow(ax, (0, 0), scale(d_real, 1.55), "#0f172a", "d real", lw=1.7)
    arrow(ax, (0, 0), scale(q_real, 1.35), "#0f172a", "q real", lw=1.7)
    arrow(ax, (0, 0), scale(d_ctrl, 1.45), "#64748b", "d ctrl", lw=1.4, ls="--")
    arrow(ax, (0, 0), scale(q_ctrl, 1.18), "#64748b", "q ctrl", lw=1.4, ls="--")

    arrow(ax, (0, 0), i_id, "#1f6feb", r"$i_d$ on ctrl d", lw=2.1)
    arrow(ax, (0, 0), i_iq, "#d73a49", r"$i_q$ on ctrl q", lw=2.4)
    arrow(ax, (0, 0), i_vec, "#7c3aed", r"$i_s$", lw=2.4)

    # Projection of iq component onto real d-axis.
    proj_end = (iq_to_d_real, 0.0)
    color = "#188038" if iq_to_d_real > 0 else "#b91c1c"
    arrow(ax, (0, -0.17), (iq_to_d_real, -0.17), color, r"$-\delta i_q$", lw=3.0)
    ax.plot([i_iq[0], iq_to_d_real], [i_iq[1], 0.0], color=color, lw=1.0, ls=":")

    arc_radius = 0.38
    theta1 = 0 if delta_sign > 0 else math.degrees(delta)
    theta2 = math.degrees(delta) if delta_sign > 0 else 0
    ax.add_patch(Arc((0, 0), arc_radius, arc_radius, theta1=theta1, theta2=theta2, color="#f97316", lw=2))
    ax.text(0.27, 0.08 * delta_sign, r"$\delta$", color="#f97316", fontsize=12)

    direction = "increase" if iq_to_d_real > 0 else "decrease"
    ax.text(
        -1.45,
        -1.35,
        "\n".join(
            [
                title,
                rf"$i_q {'>' if iq_sign > 0 else '<'} 0,\ \delta {'>' if delta_sign > 0 else '<'} 0$",
                rf"$i_{{d,real}}\simeq i_{{d,ctrl}}-\delta i_q$",
                rf"$-\delta i_q={approx:+.2f}$ -> d-real {direction}",
            ]
        ),
        fontsize=10,
        ha="left",
        va="bottom",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#f8fafc", edgecolor="#cbd5e1"),
    )

    ax.set_aspect("equal")
    ax.set_xlim(-1.65, 1.65)
    ax.set_ylim(-1.45, 1.45)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=12, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_color("#cbd5e1")


def main():
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Axis-misalignment current mixing: sign of q-axis current reverses d-axis projection",
        fontsize=15,
        fontweight="bold",
    )
    draw_case(axes[0, 0], "Motoring: iq > 0, delta > 0", iq_sign=1, delta_sign=1)
    draw_case(axes[1, 0], "Motoring: iq > 0, delta < 0", iq_sign=1, delta_sign=-1)
    draw_case(axes[0, 1], "Regen: iq < 0, delta > 0", iq_sign=-1, delta_sign=1)
    draw_case(axes[1, 1], "Regen: iq < 0, delta < 0", iq_sign=-1, delta_sign=-1)

    fig.text(
        0.5,
        0.025,
        r"Projection rule: $i_{d,real}\simeq i_{d,ctrl}-\delta i_q$. "
        r"Therefore the same $\delta$ gives opposite $i_d$ error in motoring and regeneration.",
        ha="center",
        fontsize=12,
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=170)
    print(OUT)


if __name__ == "__main__":
    main()
