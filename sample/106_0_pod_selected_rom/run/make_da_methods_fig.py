"""「データ同化とは」スライド用：手法の位置づけマップ＋EnKF vs PFの実測比較.

左: 手法マップ（線形/非線形 × ガウス/非ガウス）と今回の選択
右: 103_1で実施したEnKF vs PFの比較（同一真値・同一観測）
出力: docs/img/da_methods_map.png
再現: OMP_NUM_THREADS=4 python3 run/make_da_methods_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
IMG = os.path.join(ROOT, "docs", "img")

# 103_1/results/compare_summary.yaml の実測値
ENKF = dict(n=60, rmse_all=0.0761, rmse_unobs=0.0794)
PF   = dict(n=400, rmse_all=0.4580, rmse_unobs=0.5134, ess_min=1.02)


def box(ax, x, y, w, h, title, body, fc, ec, bold=False):
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.05",
                 fc=fc, ec=ec, lw=3.0 if bold else 1.6))
    ax.text(x, y + h/2 - 0.42, title, ha="center", va="top",
            fontsize=13 if bold else 12, weight="bold",
            color="#8a1c1c" if bold else "#1b2430")
    ax.text(x, y - 0.18, body, ha="center", va="center", fontsize=10.2, color="#333")


def method_map(ax):
    ax.set_xlim(0, 10); ax.set_ylim(0, 10); ax.axis("off")
    ax.text(5, 9.6, "データ同化の代表的な手法", ha="center", fontsize=13.5, weight="bold")
    box(ax, 2.7, 7.6, 4.6, 2.3, "OI（最適内挿）",
        "共分散 $B$ を固定して使う\n最も軽い・仕組みが分かりやすい",
        "#ffe9d6", "#c0392b", bold=True)
    box(ax, 7.6, 7.6, 4.2, 2.3, "カルマンフィルタ",
        "線形・ガウスなら最適\n共分散を毎回更新",
        "#eef2f6", "#8899aa")
    box(ax, 2.7, 4.6, 4.6, 2.3, "アンサンブルEnKF",
        "分身の散らばりから共分散を作る\n高次元の場に強い",
        "#ffe9d6", "#c0392b", bold=True)
    box(ax, 7.6, 4.6, 4.2, 2.3, "粒子フィルタ PF",
        "非ガウス・多峰でもOK\n高次元だと粒子が縮退",
        "#eef2f6", "#8899aa")
    box(ax, 5.0, 1.7, 5.4, 1.9, "変分法（3D/4D-Var）",
        "時間窓全体を最適化・随伴が要る",
        "#eef2f6", "#8899aa")
    ax.text(0.15, 6.05, "本研究で使用", fontsize=11.5, weight="bold", color="#c0392b",
            rotation=90, va="center")
    ax.plot([0.42, 0.42], [3.4, 8.8], color="#c0392b", lw=4, solid_capstyle="round")
    ax.text(5, 0.35,
            "熱伝導は線形・誤差はガウスとみなせる → カルマン系（OI・EnKF）で十分",
            ha="center", fontsize=11.5, weight="bold", color="#14459c")


def evidence(ax):
    labels = ["全セル温度\nRMSE", "未観測点\nRMSE"]
    x = np.arange(2); w = 0.36
    b1 = ax.bar(x - w/2, [ENKF["rmse_all"], ENKF["rmse_unobs"]], w,
                color="#c0392b", label=f"EnKF（{ENKF['n']}メンバー）")
    b2 = ax.bar(x + w/2, [PF["rmse_all"], PF["rmse_unobs"]], w,
                color="#8899aa", label=f"PF（{PF['n']}粒子）")
    for bb in (b1, b2):
        for r in bb:
            ax.text(r.get_x() + r.get_width()/2, r.get_height(),
                    f"{r.get_height():.3f}", ha="center", va="bottom",
                    fontsize=11, weight="bold")
    ax.axhline(0.3, color="#2e9e5b", ls="--", lw=2)
    ax.text(1.46, 0.315, "観測ノイズ 0.3 K", fontsize=10.5, color="#2e9e5b",
            ha="right", weight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=11.5)
    ax.set_ylabel("最終RMSE [K]")
    ax.set_ylim(0, 0.62)
    ax.set_title("実際に両方やってみた（同一の真値・観測）\n"
                 "PFは粒子6.7倍でもEnKFに6倍及ばず",
                 fontsize=13, weight="bold")
    ax.legend(fontsize=11, loc="upper left"); ax.grid(alpha=.3, axis="y")
    ax.text(0.5, -0.255,
            f"PFの有効粒子数は最小 {PF['ess_min']:.2f} 個まで低下（400粒子中）\n"
            "＝ 約2万次元の温度場では粒子が1個に縮退（次元の呪い）",
            transform=ax.transAxes, ha="center", fontsize=11,
            weight="bold", color="#8a1c1c",
            bbox=dict(boxstyle="round,pad=0.45", fc="#fdecea", ec="#c0392b"))


def main():
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(15.2, 6.4),
                                 gridspec_kw={"width_ratios": [1.28, 1.0]})
    method_map(a0); evidence(a1)
    fig.subplots_adjust(left=0.04, right=0.98, top=0.93, bottom=0.17, wspace=0.22)
    out = os.path.join(IMG, "da_methods_map.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
