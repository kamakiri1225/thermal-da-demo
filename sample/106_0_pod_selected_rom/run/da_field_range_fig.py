"""データ同化で、全体の温度の範囲（最高・最低）と分布がどれだけ正解に近づくかを見せる（今の ROM）。
  上：全セルの最高・最低温度の時刻歴　下：150 秒の全セルの温度（正解 vs 推定）
  da_timeseries_fig.py と同じ同化の設定。

  比べるもの（同化の設定は holdout_unobserved_cells.py と同じ：60 メンバー、30 秒ごと、乱数 5 通りの平均）
    同化なし            ：初期温度・発熱・放熱が分からないまま（初期のばらつきの平均から）ROM で進めただけ
    温度 2 点で同化      ：温度計 P2・P0
    温度 2 点＋変位 2 点 ：＋目的指向で選んだ変位計 S1・S2
  描く量：測らない場所の温度（上面・反対側 T3）、測らない反り A−O（上面 2 点の上下の変位の差）
出力：docs/img/da_field_range.png
再現：OMP_NUM_THREADS=4 python3 run/da_field_range_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import holdout_AO_realsolver as H
from holdout_AO_realsolver import (rg, enkf_update, integ_on, build_idw_mapping, cylinder_mesh, NPT, IQ, IH, NAUG,
                                   N_ENS, SIG_T, SIG_U, INFL, SEEDS, TN, SCHED, RES, IMG)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt

PLACE = (-0.028, 0.0, 0.0975)          # 上面・反対側（T3）。温度計のない場所
A_XYZ = (0.028, 0.0, 0.1005); O_XYZ = (-0.028, 0.0, 0.1005)
CASES = [("learned", "15 W・0〜300 s"), ("q25", "25 W・0〜300 s"), ("intermittent", "15 W 間欠加熱")]


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); mp = mean[pod]
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); um = op["u_mean"]; coords = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    mesh = cylinder_mesh.build_cylinder_mesh(H.NR, H.NTH, H.NZ, H.R_IN, H.R_OUT, H.H_)
    node_ids = [n for n, _ in mesh["nodes"]]; tmap = build_idw_mapping(Cc, coords)
    near = lambda x: int(np.linalg.norm(coords - np.array(x), axis=1).argmin())
    sel = [(near(x), c) for x, c in H.SEL]; iA, iO = near(A_XYZ), near(O_XYZ)
    cell = int(np.linalg.norm(Cc - np.array(PLACE), axis=1).argmin())
    warp = lambda T5: (um[iA, 2] + Wall[iA, 2] @ (T5 - mp)) - (um[iO, 2] + Wall[iO, 2] @ (T5 - mp))
    tcell = lambda T5: mean[cell] + U[cell] @ (UPp @ (T5 - mp))

    plt.rcParams.update({"font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
    fig, axs = plt.subplots(2, 3, figsize=(16, 10.5), gridspec_kw=dict(hspace=0.32, wspace=0.27))
    for j, (case, lab) in enumerate(CASES):
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]; Tf = z["Tfield"]; T5 = Tf[:, pod]
        utrue = H.truth_disp_all(case, coords, node_ids, tmap, Tf, t); on = SCHED[case]
        res = {}
        for key in ("none", "temp", "sel"):
            Ts, Ws = [], []
            for seed in SEEDS:
                rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
                Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
                Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
                obs = sel if key == "sel" else []
                Rd = np.diag([SIG_T ** 2] * len(TN) + [SIG_U ** 2] * len(obs)); tp = 0.0
                m0 = Z[:, :NPT].mean(0); Tr = [m0]; Wr = [warp(m0)]
                for ci, tb in enumerate(t[1:], 1):
                    Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, tp, tb); tp = tb
                    if key != "none":
                        yv = list(T5[ci][TN]); Yf = Z[:, TN]
                        for (i, c) in obs:
                            yv.append(utrue[ci, i, c]); Yf = np.column_stack([Yf, um[i, c] + (Z[:, :NPT] - mp) @ Wall[i, c]])
                        y = np.array(yv) + ro.normal(0, np.sqrt(np.diag(Rd)))
                        Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Yf)
                        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
                    zm = Z[:, :NPT].mean(0); Tr.append(zm); Wr.append(warp(zm))
                Ts.append(Tr); Ws.append(Wr)
            T5m = np.mean(Ts, 0); res[key] = (mean[None] + ((T5m - mp) @ UPp.T) @ U.T - rg.T_AIR_K + 20, np.mean(Ws, 0))
        Ftrue = Tf - rg.T_AIR_K + 20; k150 = int(np.argmin(abs(t - 150)))
        sty = {"none": ("#999999", "--", "同化なし（初期温度・発熱が分からないまま）"), "temp": ("#1f77b4", "-", "温度 2 点で同化"),
               "sel": ("#C0392B", "-", "温度 2 点＋変位 2 点で同化（本研究）")}
        ax = axs[0, j]
        ax.fill_between(t, Ftrue.min(1), Ftrue.max(1), color="k", alpha=.15, lw=0, label="正解の範囲（全セルの最低〜最高）")
        for key, (c, ls, l) in sty.items():
            ax.plot(t, res[key][0].max(1), ls, color=c, lw=2.4, label=l); ax.plot(t, res[key][0].min(1), ls, color=c, lw=2.4)
        ax.axvline(150, color="k", lw=1, ls=":"); ax.grid(alpha=.3); ax.set_title(lab); ax.set_xlabel("時刻 [s]"); ax.set_xticks([0, 200, 400, 600])
        ax = axs[1, j]; rs = np.random.default_rng(0).choice(Ftrue.shape[1], 3000, replace=False)
        lo, hi = Ftrue[k150].min() - 0.5, Ftrue[k150].max() + 0.5
        for key, (c, ls, l) in sty.items():
            e = res[key][0][k150]; ax.plot(Ftrue[k150, rs], e[rs], ".", color=c, ms=3, alpha=.5)
            rm = float(np.sqrt(((e - Ftrue[k150]) ** 2).mean()))
            ax.text(0.03, {"none": 0.95, "temp": 0.86, "sel": 0.77}[key], f"RMSE {rm:.2f} K", color=c, transform=ax.transAxes, va="top", fontsize=15)
        ax.plot([lo, hi], [lo, hi], "k-", lw=1); ax.set_xlim(lo, hi); ax.set_ylim(lo, max(hi, res["none"][0][k150].max() + 0.3))
        ax.grid(alpha=.3); ax.set_xlabel("正解の温度 [°C]（150 秒、全セル）")
        print(case, "150 s RMSE", {k: round(float(np.sqrt(((res[k][0][k150] - Ftrue[k150]) ** 2).mean())), 3) for k in sty}, flush=True)
    axs[0, 0].set_ylabel("全体の最高・最低温度 [°C]"); axs[1, 0].set_ylabel("推定した温度 [°C]")
    h_, l_ = axs[0, 0].get_legend_handles_labels()
    from matplotlib.patches import Patch
    fig.legend(h_, l_, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.02), frameon=False)
    fig.savefig(os.path.join(IMG, "da_field_range.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote docs/img/da_field_range.png")


if __name__ == "__main__":
    main()
