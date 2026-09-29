"""温度センサの選択基準を2つ比べる：発熱感度 dT/dQ vs 熱感度 W（FrontISTR由来）.

問い：温度センサ1点をどう選ぶべきか。
  基準A  dT/dQ（発熱→温度の感度）   … 「よく温まる点」
  基準B  |d(Uz(A)-Uz(O))/dT_i|      … 「その温度が変形に効く点」＝FrontISTRの熱感度W由来
それぞれの最大点・最小点で同化し、温度と変位の推定精度を比較する。

出力: results/sensitivity_criterion.json, docs/img/sensitivity_criterion.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensitivity_criterion_compare.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NPT, IQ, IH, NAUG = 5, 5, 6, 7
DT = 2.0; OBS_DT = 30.0; T_END = 600.0; N_ENS = 60; SIG_T = 0.30; INFL = 1.02
SEEDS = [20260913, 20260914, 20260915, 20260916, 20260917]


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Kmat = rg.tri_to_matrix(d["K_upper"], NPT)
    heat = int(d["heat_node"]); h_true = float(d["h"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]
    UPp = np.linalg.pinv(U[pod, :])
    do = np.load(os.path.join(RES, "disp_operator.npz"))
    uzm = do["uz_mean"]; D = do["Dmode"]

    # --- 2つの感度 ---
    wAO = (D[0] - D[1]) @ UPp                     # 熱感度W由来 [µm/K]
    b0 = rg.integrate_single(np.full(NPT, rg.T_AIR_K), C, Kmat, h_true, 1.0, heat, 0, 300, DT)[1][-1]
    b1 = rg.integrate_single(np.full(NPT, rg.T_AIR_K), C, Kmat, h_true, 1.1, heat, 0, 300, DT)[1][-1]
    sQ = (b1 - b0) / 0.1                          # 発熱感度 [K/W]
    print("熱感度 |d(A-O)/dT| [µm/K] :", np.round(np.abs(wAO), 3))
    print("発熱感度 dT/dQ    [K/W]   :", np.round(sQ, 2))

    cfg = [("基準A: dT/dQ 最大",  int(np.argmax(sQ)),          "#3b74b8"),
           ("基準A: dT/dQ 最小",  int(np.argmin(sQ)),          "#9db8d6"),
           ("基準B: 熱感度W 最大", int(np.argmax(np.abs(wAO))), "#c0392b"),
           ("基準B: 熱感度W 最小", int(np.argmin(np.abs(wAO))), "#e8a99f")]

    cyc = np.arange(OBS_DT, T_END + 1e-9, OBS_DT); tg = np.r_[0, cyc]
    T = np.full(NPT, rg.T_AIR_K); Ttr = [T.copy()]
    for a, b in zip(np.r_[0, cyc[:-1]], cyc):
        _, tr = rg.integrate_single(T, C, Kmat, h_true, 1.0, heat, a, b, DT)
        T = tr[-1]; Ttr.append(T.copy())
    Ttr = np.array(Ttr)

    def qoi(T5):
        a = (T5 - mean[pod]) @ UPp.T
        u = uzm + a @ D.T
        return u[..., 0] - u[..., 1]
    qt = qoi(Ttr); ht = (tg > 0) & (tg <= 300)

    def run(node, seed):
        rng = np.random.default_rng(seed); rng_o = np.random.default_rng(seed + 7)
        Z = np.zeros((N_ENS, NAUG))
        Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
        Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS)
        Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
        rec = [Z[:, :NPT].mean(0).copy()]; tp = 0.0
        Rd = np.diag([SIG_T**2])
        for ci, tb in enumerate(cyc, 1):
            Z = Z.copy()
            Z[:, :NPT] = rg.integrate_ensemble(Z[:, :NPT], C, Kmat, Z[:, IH], Z[:, IQ], heat, tp, tb, DT)
            tp = tb
            y = np.array([Ttr[ci][node]]) + rng_o.normal(0, SIG_T, 1)
            Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Z[:, [node]])
            Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
            rec.append(Z[:, :NPT].mean(0).copy())
        return np.array(rec)

    out = {"wAO": wAO.tolist(), "dTdQ": sQ.tolist(), "results": []}
    print(f"\n{'選択基準':<22}{'点':>4}{'温度RMSE[K]':>13}{'変位差誤差[µm]':>15}")
    print("-" * 56)
    names, eTs, eUs, cols = [], [], [], []
    for nm, node, col in cfg:
        recs = np.mean([run(node, s) for s in SEEDS], axis=0)
        eT = float(np.sqrt(((recs - Ttr) ** 2).mean(axis=1)[ht].mean()))
        eU = float(np.abs(qoi(recs) - qt)[ht].mean())
        print(f"{nm:<22}{'P'+str(node):>4}{eT:>13.3f}{eU:>15.3f}")
        out["results"].append(dict(criterion=nm, node=node, rmse_T=eT, err_U=eU))
        names.append(f"{nm}\n(P{node})"); eTs.append(eT); eUs.append(eU); cols.append(col)

    fig, (a0, a1) = plt.subplots(1, 2, figsize=(14.0, 5.4))
    for ax, vals, ylab, ttl in [(a0, eTs, "加熱期 全5点温度RMSE [K]", "温度の推定精度"),
                                (a1, eUs, "加熱期 変位差 Uz(A)−Uz(O) 誤差 [µm]", "未観測の変位差の推定精度")]:
        b = ax.bar(range(4), vals, color=cols, edgecolor="k", lw=.6)
        for r, v in zip(b, vals):
            ax.text(r.get_x() + r.get_width() / 2, v, f"{v:.3f}", ha="center", va="bottom",
                    fontsize=11, weight="bold")
        ax.set_xticks(range(4)); ax.set_xticklabels(names, fontsize=9.5)
        ax.set_ylabel(ylab); ax.set_title(ttl, fontsize=13, weight="bold"); ax.grid(alpha=.3, axis="y")
    fig.suptitle("温度センサ1点の選択基準を比べる：発熱感度 dT/dQ と 熱感度 W（FrontISTR由来）\n"
                 "どちらの基準でも「感度が高い点＝良い」とは限らない（5 seed平均・EnKF）",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.88])
    o = os.path.join(IMG, "sensitivity_criterion.png")
    fig.savefig(o, dpi=140); plt.close(fig)
    with open(os.path.join(RES, "sensitivity_criterion.json"), "w") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\nwrote", o)


if __name__ == "__main__":
    main()
