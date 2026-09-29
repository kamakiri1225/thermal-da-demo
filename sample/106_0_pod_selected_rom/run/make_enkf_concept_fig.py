"""EnKFの図解：アンサンブルとは何か、なぜ温度観測から発熱量が当たるのか.

A: メンバー（少しずつ違う初期値・パラメータ）を前進させた温度履歴と観測
B: サイクル0 の散布図 (温度, Q) … 無相関の雲
C: サイクル10 の散布図 … 斜めに伸びる＝相関が生まれた＝共分散ができた
出力: docs/img/enkf_concept.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_enkf_concept_fig.py
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
RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
KC = 273.15; NPT = 5; IQ = 5; IH = 6; NAUG = 7
DT = 2.0; OBS_DT = 30.0; T_END = 330.0; N = 60; SIG = 0.30; P = 2


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; K = rg.tri_to_matrix(d["K_upper"], NPT)
    h_t = float(d["h"]); heat = int(d["heat_node"])
    rng = np.random.default_rng(20260913); rng_o = np.random.default_rng(7)

    cyc = np.arange(OBS_DT, T_END + 1e-9, OBS_DT)
    T = np.full(NPT, rg.T_AIR_K); tt = [0.0]; tr = [T.copy()]
    for a, b in zip(np.r_[0, cyc[:-1]], cyc):
        ts, x = rg.integrate_single(T, C, K, h_t, 1.0, heat, a, b, DT)
        tt += list(ts[1:]); tr += list(x[1:]); T = x[-1]
    tt = np.array(tt); truth = np.array(tr)[:, P] - KC

    Z = np.zeros((N, NAUG))
    Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N, NPT))
    Z[:, IQ] = rng.uniform(0.3, 1.8, N)
    Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N), 1e-3, 0.1)
    snap0 = (Z[:, P].copy() - KC, Z[:, IQ].copy() * 15)

    hist = [Z[:, P].copy() - KC]; ht = [0.0]; tp = 0.0
    obs_t, obs_y = [], []
    for ci, tb in enumerate(cyc, 1):
        n = int(round((tb - tp) / DT)); seg = np.linspace(tp, tb, n + 1)
        for k in range(n):
            Z[:, :NPT] = rg.integrate_ensemble(Z[:, :NPT], C, K, Z[:, IH], Z[:, IQ], heat,
                                               seg[k], seg[k + 1], DT)
            hist.append(Z[:, P].copy() - KC); ht.append(seg[k + 1])
        yv = truth[np.argmin(abs(tt - tb))] + KC + rng_o.normal(0, SIG)
        obs_t.append(tb); obs_y.append(yv - KC)
        Z = enkf_update(Z, np.array([yv]), None, np.diag([SIG**2]), rng, inflation=1.02, Yf=Z[:, [P]])
        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
        hist.append(Z[:, P].copy() - KC); ht.append(tb)
        tp = tb
    hist = np.array(hist); ht = np.array(ht)
    snapN = (Z[:, P].copy() - KC, Z[:, IQ].copy() * 15)

    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(16.2, 5.2),
                                     gridspec_kw={"width_ratios": [1.5, 1, 1]})
    # A: アンサンブル
    for m in range(N):
        a0.plot(ht, hist[:, m], color="#3b74b8", lw=.7, alpha=.28)
    a0.plot(ht, hist.mean(1), color="#14459c", lw=3.0, label="アンサンブル平均")
    a0.plot(tt, truth, color="k", lw=3.0, alpha=.8, label="真値")
    a0.errorbar(obs_t, obs_y, yerr=SIG, fmt="s", color="#2e9e5b", ms=8, capsize=4,
                lw=1.8, zorder=5, label="観測")
    a0.plot([], [], color="#3b74b8", lw=1, alpha=.6, label=f"各メンバー（{N}本）")
    a0.set_xlabel("時刻 [s]"); a0.set_ylabel("温度 [℃]")
    a0.set_title("A：少しずつ違う条件のメンバーを前進させる\n"
                 "（初期温度・発熱量・放熱をばらつかせる）", fontsize=12.5, weight="bold")
    a0.legend(fontsize=9.5, loc="lower right"); a0.grid(alpha=.3)

    # サイクル番号だけだと何秒か分からないので、実時刻を併記する（観測間隔 OBS_DT=30 s）
    for ax, (Tv, Qv), ttl in [(a1, snap0, "B：サイクル0  t = 0 s（前進前）"),
                              (a2, snapN,
                               f"C：サイクル{len(cyc)}  t = {cyc[-1]:.0f} s（{len(cyc)}回前進後）")]:
        ax.scatter(Tv, Qv, s=34, color="#3b74b8", alpha=.72, edgecolors="w", lw=.5)
        r = np.corrcoef(Tv, Qv)[0, 1]
        if abs(r) > 0.25:
            c = np.polyfit(Tv, Qv, 1)
            xs = np.linspace(Tv.min(), Tv.max(), 10)
            ax.plot(xs, np.polyval(c, xs), color="#c0392b", lw=2.6)
        ax.axhline(15, color="k", ls="--", lw=1.6)
        ax.text(ax.get_xlim()[0], 15.25, " 真値 Q=15 W", fontsize=9.5, va="bottom")
        ax.set_xlabel("メンバーの温度 [℃]"); ax.set_ylabel("メンバーの発熱量 Q [W]")
        ax.set_title(f"{ttl}\n相関 = {r:+.2f}", fontsize=12.5, weight="bold",
                     color="#c0392b" if abs(r) > .25 else "#1b2430")
        ax.grid(alpha=.3)
    a2.text(.04, .04, "斜めに伸びた＝\n温度とQが連動\n→ これが共分散",
            transform=a2.transAxes, fontsize=11, weight="bold", color="#c0392b",
            bbox=dict(boxstyle="round,pad=0.35", fc="#fdecea", ec="#c0392b"))
    fig.suptitle("EnKFの中身：メンバーの「ばらつき」そのものが共分散になる ― "
                 "温度しか測らなくても、連動しているので発熱量が当たる",
                 fontsize=13.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.9])
    out = os.path.join(IMG, "enkf_concept.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print("wrote", out)
    print(f"  サイクル0 相関 {np.corrcoef(*snap0)[0,1]:+.3f} → "
          f"サイクル{len(cyc)} 相関 {np.corrcoef(*snapN)[0,1]:+.3f}")


if __name__ == "__main__":
    main()
