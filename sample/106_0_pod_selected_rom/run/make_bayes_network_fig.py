"""本研究の構造をベイジアンネットワーク（DAG）として描き、観測の価値を相互情報量で示す.

左: 発熱量Q・放熱h → 温度T → 変位u → 観測y という因果の向きを持つDAG
右: 推定対象ごとの相互情報量 I(X;y)。Δの順位と一致することを示す
出力: docs/img/bayes_network.png
再現: OMP_NUM_THREADS=4 python3 run/make_bayes_network_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
IMG = os.path.join(ROOT, "docs", "img")

B = np.array([[1.36, 0.12, 1.2], [0.12, 1.04, 0.4], [1.2, 0.4, 4.0]])
W = np.array([0.5, 0.2, 0.0]); RT = 0.09; RU = 0.01


def node(ax, x, y, w, h, txt, fc, fs=11.5, ec="#33475b"):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.05", fc=fc, ec=ec, lw=1.8))
    ax.text(x, y, txt, ha="center", va="center", fontsize=fs, weight="bold")


def arrow(ax, p, q, col="#33475b", lw=2.2, ls="-", lbl=None, dx=0, dy=0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=18,
                 color=col, lw=lw, linestyle=ls, shrinkA=2, shrinkB=2))
    if lbl:
        ax.text((p[0] + q[0]) / 2 + dx, (p[1] + q[1]) / 2 + dy, lbl,
                fontsize=9.5, color=col, ha="center", style="italic")


def dag(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    # パラメータ（推定したい未知数）
    node(ax, 1.5, 8.2, 2.2, 1.0, "発熱量 $Q$", "#fde3d6")
    node(ax, 1.5, 6.4, 2.2, 1.0, "放熱 $h$", "#fde3d6")
    # 状態
    node(ax, 5.0, 7.3, 2.6, 1.3, "温度場 $T$\n(PODで $a_1,a_2$)", "#dceafb")
    # 変位
    node(ax, 8.4, 7.3, 2.2, 1.0, "変位 $u$", "#e8f4e0")
    # 観測
    node(ax, 5.0, 3.4, 2.6, 1.0, "温度観測 $y_T$", "#f5f5b0", 11)
    node(ax, 8.4, 3.4, 2.2, 1.0, "変位観測 $y_u$", "#f5f5b0", 11)
    # 物理の矢印
    arrow(ax, (2.6, 8.2), (3.8, 7.6), lbl="発熱", dy=.28)
    arrow(ax, (2.6, 6.4), (3.8, 7.0), lbl="放熱", dy=-.30)
    arrow(ax, (6.3, 7.3), (7.3, 7.3), lbl="$W=K_s^{-1}H_T$", dy=.35)
    # 観測の矢印
    arrow(ax, (5.0, 6.65), (5.0, 3.9), col="#8a6d00", lw=1.8, ls="--",
          lbl="＋ノイズ $r_T$", dx=1.05)
    arrow(ax, (8.4, 6.8), (8.4, 3.9), col="#8a6d00", lw=1.8, ls="--",
          lbl="＋ノイズ $r_u$", dx=1.05)
    # 情報の逆流
    ax.add_patch(FancyArrowPatch((4.0, 3.2), (1.6, 5.8), arrowstyle="-|>",
                 mutation_scale=20, color="#c0392b", lw=2.6, linestyle=(0, (5, 3)),
                 connectionstyle="arc3,rad=0.28"))
    ax.text(1.9, 3.9, "観測から\nパラメータへ\n情報が逆流\n（データ同化）",
            fontsize=10.5, color="#c0392b", weight="bold", ha="center")
    ax.text(5.0, 1.3, "矢印＝因果（物理）。同化はこの矢印を逆にたどって未知数を絞り込む",
            ha="center", fontsize=10.5, color="dimgray")
    ax.set_title("① 本研究の構造をベイジアンネットワークで描くと", fontsize=13, weight="bold")


def mi_panel(ax):
    cands = [("温度点1", np.array([1., 0, 0]), RT), ("温度点2", np.array([0, 1., 0]), RT),
             ("変位 $u$", W, RU)]
    tgts = [("$T_1$", np.array([1., 0, 0])), ("$T_2$", np.array([0, 1., 0])),
            ("$Q$", np.array([0, 0, 1.])), ("全温度場", None)]
    cols = ["#3b74b8", "#2e9e5b", "#c0392b"]
    width = 0.26; xs = np.arange(len(tgts))
    for j, (cn, hh, r) in enumerate(cands):
        vals = []
        for tn, c in tgts:
            vy = hh @ B @ hh + r
            if c is None:
                D = sum((np.eye(3)[i] @ B @ hh) ** 2 for i in (0, 1)) / vy
                V = B[0, 0] + B[1, 1]
            else:
                D = (c @ B @ hh) ** 2 / vy; V = c @ B @ c
            vals.append(-0.5 * np.log(max(1e-12, 1 - D / V)) / np.log(2))
        ax.bar(xs + (j - 1) * width, vals, width, color=cols[j], label=cn)
        for x, v in zip(xs + (j - 1) * width, vals):
            ax.text(x, v, f"{v:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(xs); ax.set_xticklabels([t[0] for t in tgts], fontsize=12)
    ax.set_ylabel("相互情報量 $I(X;y)$ [bit]")
    ax.set_title("② 観測の価値＝相互情報量\n対象ごとに勝つ観測が変わる（$\\Delta$ の順位と一致）",
                 fontsize=13, weight="bold")
    ax.legend(fontsize=10.5); ax.grid(alpha=.3, axis="y")


def main():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(15.0, 6.2),
                                 gridspec_kw={"width_ratios": [1.25, 1.0]})
    dag(a0); mi_panel(a1)
    fig.tight_layout()
    out = os.path.join(IMG, "bayes_network.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
