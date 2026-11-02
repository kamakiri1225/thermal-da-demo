"""アンサンブルカルマンフィルタ（EnKF）の計算のイメージ図。数字は本番の計算（ROM＋EnKF、60 メンバー）から取る。

  条件：15 W・0〜300 秒加熱の OpenFOAM の結果を正解にし、温度計 P2・P0 の値（ノイズ 0.3 K）を 30 秒ごとに同化。
        状態は代表 5 点の温度＋発熱の倍率 q＋放熱 h（holdout_AO_realsolver.py と同じ設定、seed は 1 つ目）。
  図に描くのは P2（ヒータに一番近い代表点）の温度。
  左：3 つの時刻での「予報の散らばり（点線）→ 観測 → 補正後の散らばり（実線）」と、次の時刻への 1 期先予測（矢印）
  右：メンバーの時刻歴（灰）、平均（赤）、観測（赤い点）。30 秒ごとに補正されて縦に跳ぶ
  python3 presentation/make_enkf_image_figure.py → docs/img/enkf_image.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "run"))
from holdout_AO_realsolver import rg, enkf_update, integ_on, NPT, IQ, IH, NAUG, N_ENS, SIG_T, INFL, SEEDS, TN, SCHED
from dacore import plots as _p  # noqa: F401  （日本語フォント）
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from scipy.stats import gaussian_kde

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
K = 273.15; P = 2            # 描く代表点（P2）
ROWS = [1, 2, 3]             # 左の図に描くサイクル（30・60・90 秒）
NSHOW = 30                   # 右の図に描くメンバーの本数（計算は 60）


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"])
    pod = np.load(os.path.join(RES, "qdeim_points.npz"))["cell_idx"]
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; T5 = z["Tfield"][:, pod]
    on = SCHED["learned"]; seed = SEEDS[0]
    rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
    Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
    Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
    Rd = np.diag([SIG_T ** 2] * len(TN))
    tl = [0.0]; traj = [Z[:, P].copy()]          # 右の図の線（6 秒ごと）
    fc, an, obs = {}, {}, {}
    tp = 0.0
    for ci, tb in enumerate(t[1:], 1):
        for a, b in zip(np.arange(tp, tb, 6.0), np.arange(tp, tb, 6.0) + 6.0):
            Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, a, min(b, tb))
            tl.append(min(b, tb)); traj.append(Z[:, P].copy())
        tp = tb
        fc[ci] = Z[:, P].copy()
        y = T5[ci][TN] + ro.normal(0, SIG_T, len(TN)); obs[ci] = y[TN.index(P)]
        Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Z[:, TN])
        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
        an[ci] = Z[:, P].copy(); tl.append(tb); traj.append(Z[:, P].copy())
    tl = np.array(tl); traj = np.array(traj) - K
    for ci in ROWS:
        print(f"t={t[ci]:.0f}s 予報 平均 {fc[ci].mean()-K:.2f} 幅(標準偏差) {fc[ci].std():.2f}  観測 {obs[ci]-K:.2f}  "
              f"補正後 平均 {an[ci].mean()-K:.2f} 幅 {an[ci].std():.3f}  正解 {T5[ci][P]-K:.2f}")

    plt.rcParams.update({"font.size": 16, "axes.labelsize": 18, "xtick.labelsize": 16, "ytick.labelsize": 16})
    fig = plt.figure(figsize=(16, 8.4))
    # ---- 左：3 つの時刻の散らばり（下から上へ時間が進む） ----
    axs = [fig.add_axes([0.07, 0.07 + 0.245 * r, 0.40, 0.19]) for r in range(3)]
    lo = min(min(fc[c].min(), an[c].min()) for c in ROWS) - K - 0.5
    hi = max(max(fc[c].max(), an[c].max()) for c in ROWS) - K + 0.5
    xs = np.linspace(lo, hi, 400)
    for ax, ci in zip(axs, ROWS):
        f = fc[ci] - K; a = an[ci] - K; yo = obs[ci] - K
        kf = gaussian_kde(f)(xs); ka = gaussian_kde(a)(xs); sc = 1.0 / ka.max()
        ax.plot(xs, kf * sc, "--", color="#555", lw=2)
        ax.plot(xs, ka * sc, "-", color="#C0392B", lw=2.6)
        ax.plot(f, np.full_like(f, 0.04), "o", mfc="white", mec="#555", ms=7, mew=1.2)
        ax.plot(a, np.full_like(a, 0.04), "o", color="#C0392B", ms=6, alpha=.8)
        ax.plot([yo], [0.04], "D", color="#2E8B57", ms=13, zorder=5)
        ax.errorbar([yo], [0.04], xerr=[[0.3], [0.3]], color="#2E8B57", lw=2, capsize=5, zorder=4)
        ax.set_xlim(lo, hi); ax.set_ylim(0, 1.15); ax.set_yticks([])
        for sp in ("top", "right", "left"): ax.spines[sp].set_visible(False)
        ax.text(-0.02, 0.5, f"t = {t[ci]:.0f} 秒", transform=ax.transAxes, ha="right", va="center", fontsize=18)
        ax.text(0.99, 0.95, f"予報の幅 ±{f.std():.1f} K → 補正後 ±{a.std():.2f} K", transform=ax.transAxes,
                ha="right", va="top", fontsize=15)
        if ax is not axs[0]: ax.tick_params(labelbottom=False)
    axs[0].set_xlabel("P2 の温度 [°C]")
    # 補正後のメンバー → 次の時刻の予報（1 期先予測）を、何本かの矢印で結ぶ
    for r in range(2):
        a = an[ROWS[r]] - K; f = fc[ROWS[r + 1]] - K
        for m in range(0, N_ENS, 6):
            fig.add_artist(ConnectionPatch((a[m], 0.04), (f[m], 0.04), "data", "data", axesA=axs[r], axesB=axs[r + 1],
                                           arrowstyle="-|>", color="#4A90D9", lw=0.9, ls="--", alpha=.7, mutation_scale=10))
    fig.text(0.27, 0.80, "散らばり（確率密度）の変化", ha="center", va="bottom", fontsize=19)
    # 凡例（左上）
    from matplotlib.lines import Line2D
    hs = [Line2D([], [], ls="--", color="#555", marker="o", mfc="white", mec="#555", label="予報（ROM で 30 秒先へ）"),
          Line2D([], [], ls="", marker="D", color="#2E8B57", ms=11, label="観測（温度計 ±0.3 K）"),
          Line2D([], [], ls="-", color="#C0392B", marker="o", label="補正後（カルマンゲインで寄せる）"),
          Line2D([], [], ls="--", color="#4A90D9", label="1 期先予測（各メンバーを ROM で進める）")]
    hs.append(Line2D([], [], color="#999", lw=1.2, label="メンバー（右）"))
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.99), ncol=3, frameon=False, fontsize=15, columnspacing=1.5)
    # ---- 右：時刻歴 ----
    ax = fig.add_axes([0.56, 0.07, 0.42, 0.72])
    show = np.linspace(0, N_ENS - 1, NSHOW).astype(int)
    for m in show: ax.plot(tl, traj[:, m], "-", color="#999", lw=0.9, alpha=.75)
    ax.plot(tl, traj.mean(1), "-", color="#C0392B", lw=3, label="アンサンブルの平均（推定値）")
    ax.plot(t, T5[:, P] - K, "--", color="black", lw=2, label="正解（OpenFOAM）")
    ax.plot(t[1:], np.array([obs[c] for c in range(1, len(t))]) - K, "D", color="#2E8B57", ms=10, label="観測（30 秒ごと）", zorder=5)
    ax.set_xlim(0, 300); ax.set_ylim(15, 34); ax.set_xlabel("時刻 [秒]"); ax.set_ylabel("P2 の温度 [°C]")
    ax.grid(alpha=.3); ax.legend(loc="upper left", fontsize=15, framealpha=.95)
    ax.set_title(f"メンバーの時刻歴（60 本のうち {NSHOW} 本）", fontsize=19, pad=10)
    fig.savefig(os.path.join(IMG, "enkf_image.png"), dpi=120, bbox_inches="tight", facecolor="white")
    print("wrote docs/img/enkf_image.png")


if __name__ == "__main__":
    main()
