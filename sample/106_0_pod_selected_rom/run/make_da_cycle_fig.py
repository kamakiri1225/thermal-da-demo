"""データ同化のサイクル（予報→解析）を実際のEnKF計算から描く（スライド「データ同化とは」用）.

横軸=時刻、縦軸=温度。予報値・解析値・予報分散（スプレッド）・観測を1枚に示す。
出力: docs/img/da_cycle.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_da_cycle_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
IMG = os.path.join(ROOT, "docs", "img"); RES = os.path.join(ROOT, "results")
KC = 273.15; NPT = 5; IQ = 5; IH = 6; NAUG = 7
DT = 2.0; OBS_DT = 60.0; T_END = 420.0; N_ENS = 40; SIG_T = 0.35
P = 2   # 表示する点（ヒータ側）


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; K = rg.tri_to_matrix(d["K_upper"], NPT)
    h_t = float(d["h"]); heat = int(d["heat_node"])
    rng = np.random.default_rng(7); rng_o = np.random.default_rng(11)

    cyc = np.arange(OBS_DT, T_END + 1e-9, OBS_DT)
    # 真値
    T = np.full(NPT, rg.T_AIR_K); tt = [0.0]; tr_all = [T.copy()]
    for a, b in zip(np.r_[0, cyc[:-1]], cyc):
        ts, tr = rg.integrate_single(T, C, K, h_t, 1.0, heat, a, b, DT)
        tt += list(ts[1:]); tr_all += list(tr[1:]); T = tr[-1]
    tt = np.array(tt); truth = np.array(tr_all)[:, P] - KC

    # アンサンブル（わざと低めの初期値から始める）
    Z = np.zeros((N_ENS, NAUG))
    Z[:, :NPT] = rng.normal(rg.T_AIR_K - 1.2, 1.1, (N_ENS, NPT))
    Z[:, IQ] = rng.uniform(0.45, 1.5, N_ENS)
    Z[:, IH] = np.clip(rng.normal(0.018, 0.006, N_ENS), 1e-3, 0.1)

    segs = []          # (時刻列, 平均, 標準偏差) 予報区間ごと
    obs_t, obs_y = [], []
    ana_t, ana_m = [0.0], [Z[:, P].mean() - KC]
    tp = 0.0
    for ci, tb in enumerate(cyc, 1):
        n = int(round((tb - tp) / DT))
        ts_seg = np.linspace(tp, tb, n + 1)
        mu = [Z[:, P].mean() - KC]; sd = [Z[:, P].std()]
        Zs = Z.copy()
        for k in range(n):
            Zs[:, :NPT] = rg.integrate_ensemble(Zs[:, :NPT], C, K, Zs[:, IH], Zs[:, IQ],
                                                heat, ts_seg[k], ts_seg[k + 1], DT)
            mu.append(Zs[:, P].mean() - KC); sd.append(Zs[:, P].std())
        segs.append((ts_seg, np.array(mu), np.array(sd)))
        Z = Zs
        # 観測と更新
        ytrue = truth[np.argmin(abs(tt - tb))] + KC
        y = np.array([ytrue + rng_o.normal(0, SIG_T)])
        obs_t.append(tb); obs_y.append(y[0] - KC)
        Z = enkf_update(Z, y, None, np.diag([SIG_T**2]), rng, inflation=1.02, Yf=Z[:, [P]])
        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
        ana_t.append(tb); ana_m.append(Z[:, P].mean() - KC)
        tp = tb

    fig, ax = plt.subplots(figsize=(12.6, 6.0))
    ax.plot(tt, truth, color="k", lw=3.4, alpha=.78, label="真値（本当はこう）", zorder=2)
    for i, (ts_seg, mu, sd) in enumerate(segs):
        ax.fill_between(ts_seg, mu - sd, mu + sd, color="#3b74b8", alpha=.16, lw=0,
                        label="予報の分散（ばらつき）" if i == 0 else None, zorder=1)
        ax.plot(ts_seg, mu, color="#3b74b8", lw=2.4, ls="--",
                label="予報値（シミュレーション）" if i == 0 else None, zorder=3)
        # 解析への跳ね
        if i < len(ana_t) - 1:
            ax.plot([ts_seg[-1], ts_seg[-1]], [mu[-1], ana_m[i + 1]],
                    color="#c0392b", lw=2.6, zorder=4)
    ax.plot(ana_t[1:], ana_m[1:], "o", color="#c0392b", ms=11, zorder=6,
            label="解析値（同化後）")
    ax.errorbar(obs_t, obs_y, yerr=SIG_T, fmt="s", color="#2e9e5b", ms=9, capsize=5,
                lw=2, zorder=5, label=f"観測（ノイズ ±{SIG_T} K）")

    i0 = 2
    ts0, mu0, sd0 = segs[i0]
    j = len(ts0) // 2
    ax.annotate("予報ステップ\n（モデルで進める・分散が増える）",
                xy=(ts0[j], mu0[j] + sd0[j]),
                xytext=(212, 22.4), fontsize=12, color="#3b74b8", weight="bold",
                ha="center", arrowprops=dict(arrowstyle="-|>", color="#3b74b8", lw=1.8))
    ax.annotate("解析ステップ\n（観測で引き戻す・分散が減る）",
                xy=(60, 21.6), xytext=(118, 19.6),
                fontsize=12, color="#c0392b", weight="bold", ha="center",
                arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=1.8))
    ax.set_xlabel("時刻 [s]", fontsize=13)
    ax.set_ylabel("温度（ヒータ側の代表点）[℃]", fontsize=13)
    ax.set_title("データ同化のサイクル：予報で進め、観測が来たら解析で引き戻す（実際のEnKF計算）",
                 fontsize=14, weight="bold")
    ax.legend(fontsize=11, loc="upper left", ncol=2, framealpha=.95)
    ax.grid(alpha=.3); ax.set_xlim(0, T_END)
    fig.tight_layout()
    out = os.path.join(IMG, "da_cycle.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
