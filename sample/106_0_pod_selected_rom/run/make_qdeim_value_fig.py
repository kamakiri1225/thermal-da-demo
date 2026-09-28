"""Q-DEIMの価値を示す2枚組の図（発表スライドS4/S5用）.

左: 1点ずつ選ぶたびに残差ノルムが減る様子（貪欲選択の進行）
右: ランダム5点200通りの復元RMSE分布と、Q-DEIMの位置
出力: docs/img/qdeim_value.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_qdeim_value_fig.py
"""
from __future__ import annotations
import os, sys, importlib.util
import numpy as np
from scipy.linalg import qr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
IMG = os.path.join(ROOT, 'docs', 'img'); R = 5


def main():
    ts, C, X = m.load_snapshots()
    Xc = X - X.mean(axis=1)[:, None]
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    Ur = U[:, :R].copy()

    # --- 左: 貪欲選択の残差減衰 ---
    W = Ur.copy(); picked = []; res = []
    for k in range(R):
        n = np.linalg.norm(W, axis=1); p = int(np.argmax(n))
        picked.append(p); res.append(float(n[p]))
        v = W[p] / np.linalg.norm(W[p]); W -= np.outer(W @ v, v)
    res_after = float(np.linalg.norm(W, axis=1).max())

    def err(pts):
        a = np.linalg.pinv(Ur[pts, :]) @ Xc[pts, :]
        return float(np.sqrt(((Ur @ a - Xc) ** 2).mean()))
    e_q = err(picked)
    rng = np.random.default_rng(0)
    rs = np.array([err(list(rng.choice(len(C), R, replace=False))) for _ in range(200)])

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(13.2, 5.0))
    xs = np.arange(1, R + 1)
    b = ax0.bar(xs, res, color="#3b74b8")
    for xi, v in zip(xs, res):
        ax0.text(xi, v, f"{v:.4f}", ha="center", va="bottom", fontsize=10, weight="bold")
    ax0.bar([R + 1], [res_after], color="#bbbbbb")
    ax0.text(R + 1, res_after, "0.0000", ha="center", va="bottom", fontsize=10, weight="bold")
    ax0.set_xticks(list(xs) + [R + 1])
    ax0.set_xticklabels([f"{i}点目" for i in xs] + ["5点\n選んだ後"], fontsize=10)
    ax0.set_ylabel("残差ノルム（未説明の情報量）", fontsize=11)
    ax0.set_title("1点選ぶごとに「分からない情報」が減る\n5点で残差ゼロ＝5モードを張り切った",
                  fontsize=12.5, weight="bold")
    ax0.grid(alpha=.3, axis="y")

    ax1.hist(rs, bins=40, color="#e58f2a", alpha=.75, label="ランダム5点（200通り）")
    ax1.axvline(e_q, color="#c0392b", lw=3.0, label=f"Q-DEIM  {e_q:.5f} K")
    ax1.axvline(np.median(rs), color="#555", lw=2.0, ls="--",
                label=f"ランダム中央値  {np.median(rs):.5f} K")
    ax1.set_xscale("log")
    ax1.set_xlabel("5点から全20,696セルを復元したRMSE [K]（対数軸）")
    ax1.set_ylabel("回数")
    ax1.set_title(f"勘で置くと最悪 {rs.max():.2f} K（Q-DEIMの約{rs.max()/e_q:.0f}倍）\n"
                  "Q-DEIMの価値＝大外れを確実に避けること", fontsize=12.5, weight="bold")
    ax1.legend(fontsize=10); ax1.grid(alpha=.3)
    fig.tight_layout()
    out = os.path.join(IMG, "qdeim_value.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == '__main__':
    main()
