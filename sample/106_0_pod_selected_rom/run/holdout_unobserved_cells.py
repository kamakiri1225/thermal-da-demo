"""【未観測点の温度】データ同化で直した代表5点の温度から、gappy-POD で全 20,696 セルの温度を戻し、
観測していない場所の温度を OpenFOAM（正解）と比べる。

  データ同化は run/holdout_AO_realsolver.py と同じ（ROM＋EnKF、60 メンバー、30 秒ごと、5 seed）。
  観測：温度2点（P2・P0）のみ ／ 温度2点＋選んだ変位2点（上面 S1・S2）
  全体の温度：T = 平均場 + Φ (Φ_P)^+ (T_P − 平均場_P)   （Φ：POD モード、Φ_P：代表5点の行）
  比べる場所：温度計のない代表点 P4、代表点でない 3 セル（上面のヒータ側・反対側、反対側の中段）、全セルの RMSE

出力: results/holdout_unobserved_cells.json, docs/img/holdout_unobserved_cells.png, docs/img/holdout_field_rmse.png
再現: OMP_NUM_THREADS=4 python3 run/holdout_unobserved_cells.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import holdout_AO_realsolver as H
from holdout_AO_realsolver import (rg, enkf_update, integ_on, build_idw_mapping, cylinder_mesh, NPT, IQ, IH, NAUG,
                                   N_ENS, SIG_T, SIG_U, INFL, SEEDS, TN, SCHED, CASES, RES, IMG)
from dacore import plots as _p   # noqa: F401  （日本語フォント）
import matplotlib.pyplot as plt

PLACES = [("P4（代表点、温度計なし）", None),
          ("上面・ヒータ側", (0.028, 0.0, 0.0975)),
          ("上面・反対側", (-0.028, 0.0, 0.0975)),
          ("反対側・中段", (-0.0287, 0.0, 0.050))]
CFG = [("温度 2 点のみ", "temp", "#9AA5B1"), ("温度 2 点＋選んだ変位 2 点（本研究）", "sel", "#C0392B")]


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :])
    op = np.load(os.path.join(RES, "dispop_allnodes.npz"))
    um = op["u_mean"]; coords = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    mesh = cylinder_mesh.build_cylinder_mesh(H.NR, H.NTH, H.NZ, H.R_IN, H.R_OUT, H.H_)
    node_ids = [n for n, _ in mesh["nodes"]]
    temperature_map = build_idw_mapping(Cc, coords)
    near_node = lambda x: int(np.linalg.norm(coords - np.array(x), axis=1).argmin())
    sel = [(near_node(x), c) for x, c in H.SEL]
    mp = mean[pod]
    cells = [int(pod[4])] + [int(np.linalg.norm(Cc - np.array(x), axis=1).argmin()) for _, x in PLACES[1:]]
    names = [n for n, _ in PLACES]
    observed = set(int(pod[i]) for i in TN)
    assert not observed & set(cells)

    def recon(T5):                                   # gappy-POD：代表5点 → 全セル
        return mean[None, :] + ((T5 - mp) @ UPp.T) @ U.T

    out = {"note": "データ同化（5 seed 平均）で直した代表5点の温度から gappy-POD で全セルを戻し、OpenFOAM と比べた。"
                   "誤差は加熱期30〜300秒の平均 [K]", "cells": {n: Cc[c].tolist() for n, c in zip(names, cells)}, "cases": {}}
    hist = {}
    for case, lab in CASES:
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz"))
        t = z["times"]; Tf = z["Tfield"]; T5 = Tf[:, pod]
        utrue = H.truth_disp_all(case, coords, node_ids, temperature_map, Tf, t)
        cyc = t[1:]; on = SCHED[case]
        res = {}; hist[case] = dict(t=t, truth=Tf[:, cells])
        for nm, key, _ in CFG:
            obs = [] if key == "temp" else sel
            ws = [Wall[i, c] for i, c in obs]; u0 = [um[i, c] for i, c in obs]
            pts = []; rm = []
            for seed in SEEDS:
                rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
                Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
                Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
                Rd = np.diag([SIG_T ** 2] * len(TN) + [SIG_U ** 2] * len(ws))
                T5m = [np.full(NPT, np.nan)]; tp = 0.0
                for ci, tb in enumerate(cyc, 1):
                    Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, tp, tb); tp = tb
                    yv = list(T5[ci][TN]); Yf = Z[:, TN]
                    for (i, c), w, b0 in zip(obs, ws, u0):
                        yv.append(utrue[ci, i, c]); Yf = np.column_stack([Yf, b0 + (Z[:, :NPT] - mp) @ w])
                    y = np.array(yv) + ro.normal(0, np.sqrt(np.diag(Rd)))
                    Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Yf)
                    Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
                    T5m.append(Z[:, :NPT].mean(0))
                T5m = np.array(T5m)
                full = recon(T5m[1:])                                   # (20, 20696)
                pts.append(np.vstack([np.full(len(cells), np.nan), full[:, cells]]))
                rm.append(np.concatenate([[np.nan], np.sqrt(((full - Tf[1:]) ** 2).mean(1))]))
            pts = np.mean(pts, 0); rm = np.mean(rm, 0)
            heat_on = (t > 0) & (t <= 300)
            res[nm] = dict(field_rmse_K=float(np.nanmean(rm[heat_on])),
                           cells_abs_err_K={n: float(np.nanmean(np.abs(pts[heat_on, j] - Tf[heat_on, cells[j]])))
                                            for j, n in enumerate(names)})
            hist[case][key] = pts; hist[case][key + "_rmse"] = rm
            print(f"{case:13s} {nm:28s} 全セル {res[nm]['field_rmse_K']:.3f} K  "
                  + "  ".join(f"{n.split('（')[0]} {v:.3f}" for n, v in res[nm]['cells_abs_err_K'].items()), flush=True)
        out["cases"][case] = res
    json.dump(out, open(os.path.join(RES, "holdout_unobserved_cells.json"), "w"), ensure_ascii=False, indent=1)
    np.savez(os.path.join(RES, "holdout_unobserved_cells.npz"),
             **{f"{c}_{k}": v for c, h in hist.items() for k, v in h.items()})
    plot(hist, names, out)


def plot(hist, names, out):
    plt.rcParams.update({'font.size': 15, 'axes.titlesize': 16, 'axes.labelsize': 16, 'xtick.labelsize': 15,
                         'ytick.labelsize': 15, 'legend.fontsize': 15, 'axes.unicode_minus': False})
    fig, axs = plt.subplots(3, 4, figsize=(16, 10), sharex=True)
    for r, (case, lab) in enumerate(CASES):
        h = hist[case]; t = h["t"]
        for j, n in enumerate(names):
            ax = axs[r, j]; ax.axvspan(0, 300, color="#FDEBD0", alpha=.4, lw=0)
            ax.plot(t, h["truth"][:, j] - 273.15, "o-", color="k", lw=2.2, ms=3.5)
            for nm, key, col in CFG:
                ax.plot(t[1:], h[key][1:, j] - 273.15, "-", color=col, lw=2)
            ax.grid(alpha=.3)
            if r == 0: ax.set_title(f"T{j + 1}：{n}")
            if j == 0: ax.set_ylabel(f"{lab.split('（')[0]}\n温度 [°C]")
            if r == 2: ax.set_xlabel("時刻 [s]")
    hd = [plt.Line2D([], [], color="k", marker="o", lw=2.2, ms=3.5)] + [plt.Line2D([], [], color=c, lw=2) for _, _, c in CFG]
    lb = ["正解（OpenFOAM）"] + [f"データ同化 → gappy-POD：{nm}" for nm, _, _ in CFG]
    fig.legend(hd, lb, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.04))
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(os.path.join(IMG, "holdout_unobserved_cells.png"), dpi=100, bbox_inches="tight"); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True)
    ymax = max(np.nanmax(hist[c][k + "_rmse"]) for c, _ in CASES for _, k, _ in CFG)
    for ax, (case, lab) in zip(axs, CASES):
        h = hist[case]; t = h["t"]; ax.axvspan(0, 300, color="#FDEBD0", alpha=.4, lw=0)
        for nm, key, col in CFG:
            v = out["cases"][case][nm]["field_rmse_K"]
            ax.plot(t[1:], h[key + "_rmse"][1:], "-o", color=col, lw=2, ms=4, label=f"{nm.split('（')[0]}（{v:.2f} K）")
        ax.set_title(lab.replace("（", "\n（")); ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3); ax.set_ylim(0, ymax * 1.05)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), fontsize=14)
    axs[0].set_ylabel("全 20,696 セルの温度の RMSE [K]")
    fig.tight_layout(); fig.savefig(os.path.join(IMG, "holdout_field_rmse.png"), dpi=100, bbox_inches="tight"); plt.close(fig)
    print("wrote holdout_unobserved_cells.png, holdout_field_rmse.png")


def replot():
    z = np.load(os.path.join(RES, "holdout_unobserved_cells.npz"))
    out = json.load(open(os.path.join(RES, "holdout_unobserved_cells.json")))
    hist = {c: {k.split("_", 1)[1]: z[k] for k in z.files if k.startswith(c + "_")} for c, _ in CASES}
    plot(hist, list(out["cells"]), out)


if __name__ == "__main__":
    replot() if "--replot" in sys.argv else main()
