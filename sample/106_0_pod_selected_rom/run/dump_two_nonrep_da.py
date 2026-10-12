"""blog_009：σ で選んだ2本で、データ同化の1サイクルの中身を数字で書き出す.

selected_pair_unobserved_check.py の run() と同じ計算（15 W、seed 20260913）を1サイクルずつ追い、
予報・観測・共分散・カルマンゲイン・更新量を保存する。
出力: results/two_nonrep_da_dump.json
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import rom_general as rg
from improve_heating_schedule import integrate as integ_on

RES = os.path.join(ROOT, "results")
NPT, IQ, IH, NAUG, N_ENS, INFL = 5, 5, 6, 7, 60, 1.02
SEED = 20260913
COMPI = {"Ux": 0, "Uy": 1, "Uz": 2}


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"]); h_true = float(d["h"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); mp = mean[pod]; G = U @ UPp
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); um = op["u_mean"]; co = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    nc = lambda x: int(np.linalg.norm(Cc - np.array(x), axis=1).argmin())
    iA, iO = nn((0.028, 0, 0.1005)), nn((-0.028, 0, 0.1005))
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; Tf = z["Tfield"]
    ut = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]
    js = json.load(open(os.path.join(RES, "two_nonrep_selection.json"), encoding="utf-8"))
    on = lambda tt: (tt < 300.0)
    out = {"note": "15 W、seed 20260913、真値＝OpenFOAM＋FrontISTR。温度は K、変位は µm。q は発熱の倍率（Q＝15q W）、h は ROM の放熱 [W/K]。",
           "h_true_WperK": h_true, "configs": {}}
    for key in ["T2", "U2"]:
        picks = js["modes"][key]
        rows, offs, truth_of, labels = [], [], [], []
        for p in picks:
            x = tuple(np.array(p["xyz_mm"]) / 1000)
            if p["kind"] == "温度":
                j = nc(x); rows.append(G[j]); offs.append(mean[j]); truth_of.append(("T", j))
                labels.append(f"温度 {tuple(p['xyz_mm'])}")
            else:
                i = nn(x); c = COMPI[p["comp"]]; rows.append(Wall[i, c]); offs.append(um[i, c]); truth_of.append(("U", i, c))
                labels.append(f"変位 {p['comp']} {tuple(p['xyz_mm'])}")
        Hm = np.array(rows); off = np.array(offs); Rd = np.eye(2) * 0.09
        rng = np.random.default_rng(SEED); ro = np.random.default_rng(SEED + 7)
        Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
        Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
        cyc = []; tp = 0.0
        init = {"T5_mean_C": (Z[:, :NPT].mean(0) - 273.15).round(2).tolist(), "q_mean": round(float(Z[:, IQ].mean()), 3),
                "h_mean": round(float(Z[:, IH].mean()), 4)}
        for ci, tb in enumerate(t[1:], 1):
            Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, tp, tb); tp = tb
            Yf = off + (Z[:, :NPT] - mp) @ Hm.T                        # 予報観測（60×2）
            ytrue = np.array([Tf[ci, a[1]] if a[0] == "T" else ut[ci, a[1], a[2]] for a in truth_of])
            y = ytrue + ro.normal(0, 0.3, 2)                           # 観測（真値＋ノイズ）
            # ---- enkf_update と同じ手順を展開 ----
            zbar0 = Z.mean(0); Zi = zbar0 + INFL * (Z - zbar0)
            yb0 = Yf.mean(0); Yi = yb0 + INFL * (Yf - yb0)
            dZ = Zi - Zi.mean(0); dY = Yi - Yi.mean(0)
            Czy = dZ.T @ dY / (N_ENS - 1); Cyy = dY.T @ dY / (N_ENS - 1); S = Cyy + Rd
            K = Czy @ np.linalg.inv(S)
            eps = rng.multivariate_normal(np.zeros(2), Rd, size=N_ENS)
            innov = (y[None, :] + eps) - Yi
            Za = Zi + innov @ np.linalg.solve(S, Czy.T)
            Za[:, IQ] = np.clip(Za[:, IQ], 0, 3); Za[:, IH] = np.clip(Za[:, IH], 1e-4, 0.2)
            T5t = Tf[ci, pod]
            rec = {"t": float(tb),
                   "forecast_T5_mean_C": (Zi[:, :NPT].mean(0) - 273.15).round(3).tolist(),
                   "forecast_q": round(float(Zi[:, IQ].mean()), 3), "forecast_h": round(float(Zi[:, IH].mean()), 4),
                   "Yf_mean": Yi.mean(0).round(3).tolist(), "y_obs": y.round(3).tolist(), "y_true": ytrue.round(3).tolist(),
                   "innov_mean": innov.mean(0).round(3).tolist(),
                   "Cyy": Cyy.round(4).tolist(), "S": S.round(4).tolist(), "Czy": Czy.round(4).tolist(), "K": K.round(4).tolist(),
                   "dz_mean": (Za.mean(0) - Zi.mean(0)).round(4).tolist(),
                   "analysis_T5_mean_C": (Za[:, :NPT].mean(0) - 273.15).round(3).tolist(),
                   "analysis_q": round(float(Za[:, IQ].mean()), 3), "analysis_h": round(float(Za[:, IH].mean()), 4),
                   "truth_T5_C": (T5t - 273.15).round(3).tolist(),
                   "truth_AO": round(float(ut[ci, iA, 2] - ut[ci, iO, 2]), 3)}
            dm = Za[:, :NPT].mean(0) - mp
            rec["est_AO"] = round(float((um[iA, 2] + Wall[iA, 2] @ dm) - (um[iO, 2] + Wall[iO, 2] @ dm)), 3)
            rec["est_UzA"] = round(float(um[iA, 2] + Wall[iA, 2] @ dm), 3); rec["truth_UzA"] = round(float(ut[ci, iA, 2]), 3)
            cyc.append(rec); Z = Za
        out["configs"][key] = {"labels": labels, "H_rows": Hm.round(4).tolist(), "offset": off.round(3).tolist(), "init": init, "cycles": cyc}
        c1 = cyc[0]; print(key, labels); print(" t=30 Yf", c1["Yf_mean"], "y", c1["y_obs"], "K", np.array(c1["K"]).round(3).tolist())
        print(" q traj", [c["analysis_q"] for c in cyc[:10]], "AO est/truth", [(c["est_AO"], c["truth_AO"]) for c in cyc[:10:3]])
    json.dump(out, open(os.path.join(RES, "two_nonrep_da_dump.json"), "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
