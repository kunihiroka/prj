import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch


def add_box(ax, xy, w, h, text, fc="#ffffff", ec="#2f4a6d", lw=1.6, fs=10, weight="normal"):
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
    )
    ax.add_patch(box)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, weight=weight)
    return box


def add_arrow(ax, p1, p2, color="#2f4a6d", lw=1.6, text=None, tpos=0.5, dy=0.0):
    arr = FancyArrowPatch(
        p1,
        p2,
        arrowstyle="-|>",
        mutation_scale=12,
        linewidth=lw,
        color=color,
        shrinkA=4,
        shrinkB=4,
    )
    ax.add_patch(arr)
    if text:
        x = p1[0] + (p2[0] - p1[0]) * tpos
        y = p1[1] + (p2[1] - p1[1]) * tpos + dy
        ax.text(x, y, text, ha="center", va="center", fontsize=8.5, color=color)
    return arr


def main():
    plt.rcParams["font.family"] = ["Yu Gothic", "Meiryo", "DejaVu Sans"]
    fig, ax = plt.subplots(figsize=(15.5, 8.8), dpi=180)
    ax.set_xlim(0, 15.5)
    ax.set_ylim(0, 8.8)
    ax.axis("off")
    fig.patch.set_facecolor("#f6f8fb")
    ax.set_facecolor("#f6f8fb")

    ax.text(
        0.35,
        8.35,
        "MIMO current-loop return ratio model",
        fontsize=20,
        weight="bold",
        color="#18263a",
    )
    ax.text(
        0.36,
        8.02,
        "PI入力でループを切り、滑り・磁束・角度・座標変換を含む実効プラント G_i(s) を線形化する",
        fontsize=11,
        color="#526274",
    )

    # Outer current loop
    add_box(ax, (0.45, 6.45), 1.75, 0.7, "電流指令\n$i_d^*, i_q^*$", fc="#eaf2ff")
    sum_c = Circle((2.95, 6.80), 0.28, edgecolor="#2f4a6d", facecolor="#ffffff", linewidth=1.6)
    ax.add_patch(sum_c)
    ax.text(2.88, 6.91, "+", fontsize=12, ha="center", va="center")
    ax.text(3.03, 6.66, "-", fontsize=12, ha="center", va="center")
    add_box(ax, (3.55, 6.45), 1.55, 0.7, "電流誤差\n$e_d,e_q$", fc="#ffffff")
    add_box(ax, (5.55, 6.32), 1.95, 0.96, "PI制御器\n$C_i(s)=K_p(1+K_i/s)$", fc="#fff5df", ec="#9c6b12")
    add_box(ax, (8.0, 6.45), 1.65, 0.7, "PI電圧\n$v_{d,PI},v_{q,PI}$", fc="#ffffff")
    add_box(ax, (10.15, 6.22), 2.55, 1.15, "実効プラント\n$G_i(s)$\nPI電圧 -> 検出dq電流", fc="#eaf7ef", ec="#1d6f42", weight="bold")
    add_box(ax, (13.3, 6.45), 1.75, 0.7, "検出電流\n$i_{d,meas}, i_{q,meas}$", fc="#ffffff")

    add_arrow(ax, (2.2, 6.8), (2.67, 6.8))
    add_arrow(ax, (3.23, 6.8), (3.55, 6.8))
    add_arrow(ax, (5.1, 6.8), (5.55, 6.8))
    add_arrow(ax, (7.5, 6.8), (8.0, 6.8))
    add_arrow(ax, (9.65, 6.8), (10.15, 6.8))
    add_arrow(ax, (12.7, 6.8), (13.3, 6.8))
    add_arrow(ax, (14.18, 6.45), (14.18, 5.85), color="#7a8797")
    add_arrow(ax, (14.18, 5.85), (2.95, 5.85), color="#7a8797", text="負帰還", tpos=0.88, dy=0.17)
    add_arrow(ax, (2.95, 5.85), (2.95, 6.52), color="#7a8797")

    # Detail container
    detail = FancyBboxPatch(
        (0.45, 0.55),
        14.6,
        5.05,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.2,
        edgecolor="#aab8c9",
        facecolor="#ffffff",
    )
    ax.add_patch(detail)
    ax.text(0.75, 5.25, "$G_i(s)$ に含めた内部経路", fontsize=14, weight="bold", color="#18263a")

    # Internal path boxes
    add_box(ax, (0.85, 4.05), 1.35, 0.62, "$v_{PI}$", fc="#f7fbff")
    add_box(ax, (2.75, 3.88), 1.65, 0.96, "非干渉項\n$v_{ff}(x)$", fc="#fff5df", ec="#9c6b12")
    add_box(ax, (5.0, 4.05), 1.55, 0.62, "印加電圧\n$v_d,v_q$", fc="#f7fbff")
    add_box(ax, (7.0, 3.75), 2.4, 1.2, "誘導機状態方程式\n$\\psi_{sd},\\psi_{sq}$\n$\\psi_{rd},\\psi_{rq}$", fc="#eaf7ef", ec="#1d6f42")
    add_box(ax, (10.0, 4.22), 1.65, 0.72, "実dq電流\n$i_{d,a},i_{q,a}$", fc="#f7fbff")
    add_box(ax, (10.0, 3.15), 1.65, 0.72, "実ロータ磁束\n$\\Phi_a,\\theta_{flux}$", fc="#f7fbff")
    add_box(ax, (12.75, 4.08), 1.75, 0.82, "dq検出\n座標変換", fc="#eef1ff", ec="#4c5eb8")
    add_box(ax, (10.0, 2.05), 1.65, 0.78, "滑り計算\n$\\omega_{slip,cmd}$", fc="#eef1ff", ec="#4c5eb8")
    add_box(ax, (7.25, 2.05), 1.8, 0.78, "磁束推定器\n$\\hat\\phi$", fc="#eef1ff", ec="#4c5eb8")
    add_box(ax, (12.75, 2.05), 1.75, 0.78, "制御角\n$\\theta_{ctrl}$", fc="#eef1ff", ec="#4c5eb8")
    add_box(ax, (12.75, 0.95), 1.75, 0.74, "検出dq電流\n$i_{d,meas},i_{q,meas}$", fc="#f7fbff")

    add_arrow(ax, (2.2, 4.36), (2.75, 4.36))
    add_arrow(ax, (4.4, 4.36), (5.0, 4.36))
    add_arrow(ax, (6.55, 4.36), (7.0, 4.36))
    add_arrow(ax, (9.4, 4.58), (10.0, 4.58))
    add_arrow(ax, (9.4, 3.55), (10.0, 3.51))
    add_arrow(ax, (11.65, 4.58), (12.75, 4.52), text="電流", dy=0.16)
    add_arrow(ax, (11.65, 3.51), (12.75, 4.28), text="実磁束角", dy=0.22)
    add_arrow(ax, (14.5, 4.49), (14.85, 4.49))
    ax.text(14.95, 4.49, "出力 y", fontsize=9, va="center", color="#2f4a6d")

    add_arrow(ax, (10.82, 4.22), (8.15, 2.83), color="#7a8797", text="$i_d$源", tpos=0.55, dy=0.12)
    add_arrow(ax, (9.05, 2.44), (10.0, 2.44), color="#7a8797", text="$\\hat\\phi$", tpos=0.45, dy=0.15)
    add_arrow(ax, (10.82, 3.15), (10.82, 2.83), color="#7a8797", text="$\\Phi_a$", tpos=0.45, dy=0.14)
    add_arrow(ax, (11.65, 2.44), (12.75, 2.44), color="#7a8797")
    add_arrow(ax, (13.62, 2.83), (13.62, 4.08), color="#7a8797", text="$\\theta_{ctrl}$", tpos=0.55, dy=0.16)
    add_arrow(ax, (13.62, 4.08), (13.62, 1.69), color="#7a8797")

    # Linearization equations
    eq_box = FancyBboxPatch(
        (0.75, 0.78),
        8.9,
        1.0,
        boxstyle="round,pad=0.04,rounding_size=0.05",
        linewidth=1.2,
        edgecolor="#c6d0dc",
        facecolor="#f8fafc",
    )
    ax.add_patch(eq_box)
    ax.text(
        1.0,
        1.38,
        "$\\Delta i_{meas}(s)=G_i(s)\\Delta v_{PI}(s)$     "
        "$L_i(s)=C_i(s)G_i(s)$     "
        "$\\det(I+L_i(j\\omega))\\to 0$ で低減衰/不安定境界に接近",
        fontsize=12,
        color="#18263a",
        va="center",
    )
    ax.text(
        1.0,
        1.03,
        "$G_i(s)$ は2入力2出力。非干渉、磁束推定、滑り計算、角度生成、dq変換を含めて線形化する。",
        fontsize=9.5,
        color="#526274",
        va="center",
    )

    out = "outputs/mimo_return_ratio_block_diagram.png"
    fig.savefig(out, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(out)


if __name__ == "__main__":
    main()
