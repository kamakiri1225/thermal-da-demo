"""σ の意味を散布図で見せる（blog_008 §5-2）.

4000 本の ROM 計算（t=150 秒）で、横軸＝候補センサの値（ノイズ込み）、縦軸＝反り A−O を描く。
センサの値が分かったときに残る反りの縦のばらつき（回帰直線からのずれの標準偏差）が σ。
  σ² = Var(d) − Cov(d, y)² / (Var(y) + r)
出力: docs/img/sigma_scatter.png
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from dacore import rom_general as rg

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
plt.rcParams.update({"font.size": 15})


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], 5); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]
    UPp = np.linalg.pinv(U[pod, :])
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    w = Wall[nn((0.028, 0, 0.1005)), 2] - Wall[nn((-0.028, 0, 0.1005)), 2]
    cands = [("① 変位計 S1（上面、横方向 Ux）", Wall[nn((0.0097, -0.0362, 0.1005)), 0], "µm", "#C0392B",
              "反りと強く一緒に動く\n→ 値が分かれば反りがほぼ決まる"),
             ("② 変位計 B（中段、上下方向 Uz）", Wall[nn((0.028, 0, 0.075)), 2], "µm", "#E07A1F",
              "反りとの関係が弱い\n→ 値が分かっても反りは絞れない"),
             ("③ 変位計（底面近く、Uz）", Wall[nn((-0.0375, 0, 0.005)), 2], "µm", "#8E44AD",
              "ほとんど動かず、ノイズ 0.3 µm に埋もれる\n→ 値が分かっても反りは絞れない"),
             ("④ 温度計 P4（反対側の下部）", np.eye(5)[4], "K", "#2E6FD8",
              "反りと無関係に上がり下がりする\n→ 1本だけでは反りは分からない")]
    rng = np.random.default_rng(20260930); N = 4000
    Tp = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N, 5)); qp = rng.uniform(0.3, 1.8, N)
    hp = np.clip(rng.normal(0.02, 0.01, N), 1e-3, 0.1)
    Tp = rg.integrate_ensemble(Tp, C, Km, hp, qp, heat, 0.0, 150.0, 2.0)
    dT = Tp - Tp.mean(0); dd = dT @ w; s0 = dd.std(ddof=1)
    rn = np.random.default_rng(1)
    fig, axs = plt.subplots(1, 4, figsize=(22, 6.4), sharey=True)
    for ax, (title, h, unit, col, note) in zip(axs, cands):
        y = dT @ h + rn.normal(0, 0.3, N)                      # センサの値（ノイズ 0.3 込み）
        cov = np.cov(dd, y)[0, 1]; vy = y.var(ddof=1); b = cov / vy
        sig = np.sqrt(dd.var(ddof=1) - cov ** 2 / vy); rho = cov / np.sqrt(vy * dd.var(ddof=1))
        ax.scatter(y, dd, s=4, alpha=.25, color=col, edgecolor="none")
        xs = np.linspace(np.percentile(y, 1), np.percentile(y, 99), 50)
        ax.plot(xs, b * xs, color="#14213D", lw=2.5)
        ax.fill_between(xs, b * xs - sig, b * xs + sig, color="#14213D", alpha=.12)
        ax.set_title(title, fontsize=16)
        ax.set_xlabel(f"センサの値（平均からのずれ）[{unit}]")
        ax.text(0.03, 0.97, f"相関 {rho:+.2f}\n測った後の幅 σ = {sig:.2f} µm", transform=ax.transAxes, va="top",
                fontsize=16, weight="bold", color=col, bbox=dict(fc="white", ec="none", alpha=.85))
        ax.text(0.03, 0.03, note, transform=ax.transAxes, va="bottom", fontsize=14, color="#333",
                bbox=dict(fc="white", ec="none", alpha=.85))
        ax.grid(alpha=.3)
    axs[0].set_ylabel("反り A−O（平均からのずれ）[µm]")
    fig.suptitle(f"σ ＝ センサの値が分かった後に残る、反りの縦の幅（測る前の幅は ±{s0:.2f} µm。4000 本の ROM 計算、t = 150 秒）",
                 fontsize=19, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(IMG, "sigma_scatter.png"), dpi=100, facecolor="white"); plt.close(fig)
    print("wrote sigma_scatter.png")


if __name__ == "__main__":
    main()
