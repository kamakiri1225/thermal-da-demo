"""逆距離加重（IDW）補間の概念図：どの点からどの点へ、どう重みを付けるか.

左：CFDはセル中心に、FEMは節点に値を持つ（位置が一致しない）
右：1つの節点について、近いk個のセル中心から距離の逆数で加重平均する
出力: docs/img/idw_concept.png
再現: OMP_NUM_THREADS=4 python3 run/make_idw_concept_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle
IMG = os.path.join(ROOT, "docs", "img")


def panel_left(ax):
    ax.set_title("① CFDとFEMは「値を持つ場所」が違う", fontsize=13.5, weight="bold")
    for i in range(4):
        for j in range(3):
            ax.add_patch(Rectangle((i, j), 1, 1, fc="#eaf1f8", ec="#9fb6cc", lw=1.4))
            ax.plot(i + .5, j + .5, "o", color="#2e6fb8", ms=11, zorder=3)
    for i in range(5):
        for j in range(4):
            ax.plot(i, j, "o", mfc="none", mec="#c0392b", mew=2.4, ms=13, zorder=4)
    ax.plot([], [], "o", color="#2e6fb8", ms=11, label="CFDのセル中心（温度を持つ）")
    ax.plot([], [], "o", mfc="none", mec="#c0392b", mew=2.4, ms=13, label="FEMの節点（温度が必要）")
    ax.legend(fontsize=11, loc="upper center", bbox_to_anchor=(.5, -.04), frameon=False)
    ax.text(2, 3.42, "有限体積法（FVM）はセルの中心、\n有限要素法（FEM）は節点に値を置く\n→ 同じ形状でも位置がずれる",
            ha="center", fontsize=11.5, color="#33475b")
    ax.set_xlim(-.7, 4.7); ax.set_ylim(-1.0, 4.35); ax.set_aspect("equal"); ax.axis("off")


def panel_right(ax):
    ax.set_title("② 近い $k$ 個のセルから「距離の逆数」で加重平均", fontsize=13.5, weight="bold")
    node = np.array([0.0, 0.0])
    cells = np.array([[-0.9, 0.5], [1.1, 0.35], [0.35, -1.15], [-1.35, -1.0]])
    temps = [25.0, 23.0, 24.0, 22.0]
    d = np.linalg.norm(cells - node, axis=1)
    w = 1 / d; wn = w / w.sum()
    for c, dd, ww, T in zip(cells, d, wn, temps):
        ax.plot([node[0], c[0]], [node[1], c[1]], "-", color="#2e6fb8",
                lw=1.0 + 9 * ww, alpha=.55, zorder=1)
        ax.plot(*c, "o", color="#2e6fb8", ms=14, zorder=3)
        off = np.sign(c) * 0.3
        ax.annotate(f"$T$={T:.0f}℃\n$d$={dd:.2f}\n$w$={ww:.2f}",
                    xy=c, xytext=c + off * [1.25, 1.1], fontsize=10.5, ha="center",
                    color="#1b2430",
                    bbox=dict(boxstyle="round,pad=0.3", fc="#eaf1f8", ec="#9fb6cc"))
        mid = (node + c) / 2
        ax.text(mid[0], mid[1], f"{dd:.2f}", fontsize=9.5, color="#2e6fb8",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=.85))
    ax.plot(*node, "*", color="#c0392b", ms=30, zorder=5)
    Tp = float(np.sum(wn * np.array(temps)))
    ax.annotate(f"節点の温度\n$T_p$ = {Tp:.2f} ℃", xy=node, xytext=(0.05, 1.28),
                fontsize=12.5, weight="bold", color="#c0392b", ha="center",
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.8))
    ax.text(0, -2.25,
            r"$T_p=\dfrac{\sum_c w_c\,T_c}{\sum_c w_c},\qquad w_c=\dfrac{1}{|x_p-x_c|}$",
            fontsize=14, ha="center")
    ax.text(0, -2.85, "近いセルほど重みが大きい（距離が半分なら重みは2倍）\n"
                      "節点がセル中心と一致したら、その値をそのまま使う",
            fontsize=11, ha="center", color="#33475b")
    ax.set_xlim(-2.6, 2.6); ax.set_ylim(-3.2, 1.9); ax.set_aspect("equal"); ax.axis("off")


def main():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(14.4, 6.0),
                                 gridspec_kw={"width_ratios": [1.0, 1.25]})
    panel_left(a0); panel_right(a1)
    fig.suptitle("OpenFOAM（セル中心）→ FrontISTR（節点）への温度の渡し方：逆距離加重（IDW）",
                 fontsize=14.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(IMG, "idw_concept.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
