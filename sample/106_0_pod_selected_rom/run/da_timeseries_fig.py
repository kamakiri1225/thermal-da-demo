"""データ同化でどれだけ良くなるかを、温度と変位の時刻歴で見せる（今の ROM：代表 5 点の熱回路）。

  比べるもの（同化の設定は holdout_unobserved_cells.py と同じ：60 メンバー、30 秒ごと、乱数 5 通りの平均）
    同化なし            ：初期温度・発熱・放熱が分からないまま（初期のばらつきの平均から）ROM で進めただけ
    温度 2 点で同化      ：温度計 P2・P0
    温度 2 点＋変位 2 点 ：＋目的指向で選んだ変位計 S1・S2
  描く量：測らない場所の温度（上面・反対側 T3）、測らない反り A−O（上面 2 点の上下の変位の差）
出力：docs/img/da_timeseries.png
再現：OMP_NUM_THREADS=4 python3 run/da_timeseries_fig.py
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
    uz = lambda i, T5: um[i, 2] + Wall[i, 2] @ (T5 - mp)
    figAO, axAO = plt.subplots(2, 3, figsize=(16, 9.2), sharex=True, gridspec_kw=dict(hspace=0.14, wspace=0.25))

    plt.rcParams.update({"font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.2), sharex=True, gridspec_kw=dict(hspace=0.14, wspace=0.25))
    for j, (case, lab) in enumerate(CASES):
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]; Tf = z["Tfield"]; T5 = Tf[:, pod]
        utrue = H.truth_disp_all(case, coords, node_ids, tmap, Tf, t); on = SCHED[case]
        res = {}
        for key in ("none", "temp", "sel"):
            Ts, Ws, As, Os = [], [], [], []
            for seed in SEEDS:
                rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
                Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
                Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
                obs = sel if key == "sel" else []
                Rd = np.diag([SIG_T ** 2] * len(TN) + [SIG_U ** 2] * len(obs)); tp = 0.0
                m0 = Z[:, :NPT].mean(0); Tr = [tcell(m0)]; Wr = [warp(m0)]; Ar = [uz(iA, m0)]; Or = [uz(iO, m0)]
                for ci, tb in enumerate(t[1:], 1):
                    Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, tp, tb); tp = tb
                    if key != "none":
                        yv = list(T5[ci][TN]); Yf = Z[:, TN]
                        for (i, c) in obs:
                            yv.append(utrue[ci, i, c]); Yf = np.column_stack([Yf, um[i, c] + (Z[:, :NPT] - mp) @ Wall[i, c]])
                        y = np.array(yv) + ro.normal(0, np.sqrt(np.diag(Rd)))
                        Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Yf)
                        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
                    zm = Z[:, :NPT].mean(0); Tr.append(tcell(zm)); Wr.append(warp(zm)); Ar.append(uz(iA, zm)); Or.append(uz(iO, zm))
                Ts.append(Tr); Ws.append(Wr); As.append(Ar); Os.append(Or)
            res[key] = (np.mean(Ts, 0) - rg.T_AIR_K + 20, np.mean(Ws, 0), np.mean(As, 0), np.mean(Os, 0))
        Ttrue = Tf[:, cell] - rg.T_AIR_K + 20; Wtrue = utrue[:, iA, 2] - utrue[:, iO, 2]
        sty = {"none": ("#999999", "--", "同化なし（初期温度・発熱が分からないまま）"), "temp": ("#1f77b4", "-", "温度 2 点で同化"),
               "sel": ("#C0392B", "-", "温度 2 点＋変位 2 点で同化（本研究）")}
        for row, (tru, idx) in enumerate(((Ttrue, 0), (Wtrue, 1))):
            ax = axs[row, j]
            ax.axvspan(0, 300 if case != "intermittent" else 150, color="#FFF4B8", lw=0)
            if case == "intermittent": ax.axvspan(300, 450, color="#FFF4B8", lw=0)
            ax.plot(t, tru, "-", color="k", lw=7, alpha=.22, label="正解（OpenFOAM＋FrontISTR）")
            if row == 0 and j == 0: ax.text(0.97, 0.04, "上面・反対側（温度計なし）", transform=ax.transAxes, ha="right", fontsize=14)
            for key, (c, ls, l) in sty.items():
                ax.plot(t, res[key][idx], ls, color=c, lw=2.6, marker="o" if key != "none" else None, ms=4, label=l)
            ax.grid(alpha=.3)
        axs[0, j].set_title(lab); axs[1, j].set_xlabel("時刻 [s]"); axs[1, j].set_xticks([0, 200, 400, 600])
        for row, (i_, idx) in enumerate(((iA, 2), (iO, 3))):
            ax = axAO[row, j]
            ax.axvspan(0, 300 if case != "intermittent" else 150, color="#FFF4B8", lw=0)
            if case == "intermittent": ax.axvspan(300, 450, color="#FFF4B8", lw=0)
            ax.plot(t, utrue[:, i_, 2], "-", color="k", lw=7, alpha=.22, label="正解（OpenFOAM＋FrontISTR）")
            for key, (c, ls, l) in sty.items():
                ax.plot(t, res[key][idx], ls, color=c, lw=2.6, marker="o" if key != "none" else None, ms=4, label=l)
            ax.grid(alpha=.3)
        axAO[0, j].set_title(lab); axAO[1, j].set_xlabel("時刻 [s]"); axAO[1, j].set_xticks([0, 200, 400, 600])
        print(case, "温度の誤差（同化なし/温度/温度＋変位、30〜300s 平均）",
              [round(float(np.abs(res[k][0][1:11] - Ttrue[1:11]).mean()), 3) for k in ("none", "temp", "sel")],
              "反りの誤差", [round(float(np.abs(res[k][1][1:11] - Wtrue[1:11]).mean()), 3) for k in ("none", "temp", "sel")], flush=True)
    axs[0, 0].set_ylabel("測らない点の温度 [°C]"); axs[1, 0].set_ylabel("測らない反り A−O [µm]")
    h_, l_ = axs[0, 0].get_legend_handles_labels()
    from matplotlib.patches import Patch
    fig.legend(h_ + [Patch(fc="#FFF4B8")], l_ + ["ヒータ ON"], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.02), frameon=False)
    fig.savefig(os.path.join(IMG, "da_timeseries.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote docs/img/da_timeseries.png")
    axAO[0, 0].set_ylabel("A の上下の変位 [µm]"); axAO[1, 0].set_ylabel("O の上下の変位 [µm]")
    axAO[0, 0].text(0.97, 0.04, "A：上面・ヒータ側", transform=axAO[0, 0].transAxes, ha="right", fontsize=14)
    axAO[1, 0].text(0.97, 0.04, "O：上面・反対側", transform=axAO[1, 0].transAxes, ha="right", fontsize=14)
    h_, l_ = axAO[0, 0].get_legend_handles_labels()
    figAO.legend(h_ + [Patch(fc="#FFF4B8")], l_ + ["ヒータ ON"], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.02), frameon=False)
    figAO.savefig(os.path.join(IMG, "da_AO.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(figAO)
    print("wrote docs/img/da_AO.png")


if __name__ == "__main__":
    main()
