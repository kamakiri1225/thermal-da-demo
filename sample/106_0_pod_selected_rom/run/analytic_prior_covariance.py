"""4,000 本の乱数計算の代わりに、温度のばらつき B(t) を式で直接計算する（共分散の伝播）.

ROM（代表5点の熱回路）
    C dT/dt = (K − diag(ΣK) − h I) T + h T_air 1 + on(t) q Q_rated e_heat
分からない値：初期温度 T0（5点それぞれ一様乱数 −3〜+12 K）、発熱の倍率 q（一様 0.3〜1.8）、放熱 h（正規 0.02±0.01）
7個の値 x = (T, q, h) をまとめ、平均の軌道のまわりで一次近似すると
    dP/dt = F(t) P + P F(t)^T,   F = [[A(h̄), g, −C^{-1}(T̄ − T_air)], [0, 0, 0], [0, 0, 0]]
P(0) = diag(15²/12 ×5, 1.5²/12, 0.01²)。B(t) = P(t) の左上 5×5。
乱数を使わないので、結果はぶれず、本数を決める必要もない。

比較：乱数 4,000 本・40,000 本の B、選定の点数（反り A−O の σ、30〜300 秒平均）と選ばれる場所。
出力: results/analytic_prior_covariance.json
再現: OMP_NUM_THREADS=4 python3 run/analytic_prior_covariance.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import rom_general as rg

RES = os.path.join(ROOT, "results")
TIMES = np.arange(30.0, 300.0 + 1e-9, 30.0); DT = 2.0


def setup():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], 5); heat = int(d["heat_node"])
    return C, Km, heat


def analytic_B(C, Km, heat):
    n = 5; Ci = 1.0 / C
    Tair = rg.T_AIR_K; Qr = rg.HEATER_RATED_W
    mT0 = Tair + 4.5; mq = 1.05; mh = 0.02                   # 一様 −3〜+12 K の平均、一様 0.3〜1.8 の平均
    P = np.diag([15.0**2 / 12] * n + [1.5**2 / 12, 0.01**2])
    T = np.full(n, mT0)
    Ksum = Km.sum(1)
    def A(h): return (Km - np.diag(Ksum) - h * np.eye(n)) * Ci[:, None]
    g = np.zeros(n); g[heat] = Qr * Ci[heat]
    def f_mean(T):   # 平均の軌道（加熱中 on=1）
        return A(mh) @ T + mh * Tair * Ci + g * mq
    def F_of(T):
        F = np.zeros((7, 7)); F[:5, :5] = A(mh); F[:5, 5] = g; F[:5, 6] = -Ci * (T - Tair); return F
    def dP(P, T):
        F = F_of(T); return F @ P + P @ F.T
    Bs = []; t = 0.0
    for tb in TIMES:
        steps = int(round((tb - t) / DT))
        for _ in range(steps):   # RK4 で平均と共分散を一緒に進める
            k1T = f_mean(T); k1P = dP(P, T)
            T2 = T + .5 * DT * k1T; k2T = f_mean(T2); k2P = dP(P + .5 * DT * k1P, T2)
            T3 = T + .5 * DT * k2T; k3T = f_mean(T3); k3P = dP(P + .5 * DT * k2P, T3)
            T4 = T + DT * k3T; k4T = f_mean(T4); k4P = dP(P + DT * k3P, T4)
            T = T + DT / 6 * (k1T + 2 * k2T + 2 * k3T + k4T); P = P + DT / 6 * (k1P + 2 * k2P + 2 * k3P + k4P)
        t = tb; Bs.append(P[:5, :5].copy())
    return Bs


def mc_B(C, Km, heat, N, seed=20260930):
    rng = np.random.default_rng(seed)
    Tp = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N, 5)); qp = rng.uniform(0.3, 1.8, N)
    hp = np.clip(rng.normal(0.02, 0.01, N), 1e-3, 0.1)
    Bs = []; tp = 0.0
    for tb in TIMES:
        Tp = rg.integrate_ensemble(Tp, C, Km, hp, qp, heat, tp, tb, DT); tp = tb; Bs.append(np.cov(Tp.T))
    return Bs


def scores(Bs, w, Hfix, Hc, r=0.09):
    """温度計 P2・P0 を固定し、候補を1本足したときの反りの σ（10時刻平均）"""
    s = np.zeros(len(Hc))
    for B in Bs:
        S = Hfix @ B @ Hfix.T + r * np.eye(len(Hfix))
        P = B - B @ Hfix.T @ np.linalg.solve(S, Hfix @ B)
        HP = Hc @ P
        s += np.sqrt(np.maximum(w @ P @ w - (HP @ w) ** 2 / (np.einsum("ij,ij->i", HP, Hc) + r), 0))
    return s / len(Bs)


def main():
    C, Km, heat = setup()
    t0 = time.time(); Ba = analytic_B(C, Km, heat); t_an = time.time() - t0
    t0 = time.time(); B4 = mc_B(C, Km, heat, 4000); t4 = time.time() - t0
    t0 = time.time(); B40 = mc_B(C, Km, heat, 40000); t40 = time.time() - t0
    t0 = time.time(); B400 = mc_B(C, Km, heat, 400); t400 = time.time() - t0
    rel = lambda X, Y: float(np.mean([np.linalg.norm(x - y) / np.linalg.norm(y) for x, y in zip(X, Y)]))
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]
    W = np.einsum("nck,kj->ncj", op["D"], np.linalg.pinv(U[pod, :]))
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    iA, iO = nn((0.028, 0, 0.1005)), nn((-0.028, 0, 0.1005)); w = W[iA, 2] - W[iO, 2]
    far = (np.linalg.norm(co - co[iA], axis=1) > 0.012) & (np.linalg.norm(co - co[iO], axis=1) > 0.012)
    Wf = W.reshape(-1, 5); idx = np.where(np.repeat(far, 3))[0]; Hc = Wf[idx]
    Hfix = np.eye(5)[[2, 0]]
    res = {}
    for name, Bs in [("analytic", Ba), ("mc400", B400), ("mc4000", B4), ("mc40000", B40)]:
        sc = scores(Bs, w, Hfix, Hc); k = int(sc.argmin()); i, c = divmod(int(idx[k]), 3)
        res[name] = {"best_xyz_mm": np.round(co[i] * 1000, 1).tolist(), "best_comp": "xyz"[c], "best_sigma": float(sc[k]),
                     "sigma_none": float(np.mean([np.sqrt(w @ B @ w) for B in Bs])), "scores": sc}
    def rankcorr(a, b):
        ra = np.argsort(np.argsort(a)); rb = np.argsort(np.argsort(b)); return float(np.corrcoef(ra, rb)[0, 1])
    top = lambda a, n=100: set(np.argsort(a)[:n])
    out = {"note": "B(t)：温度のばらつき（30〜300秒の10時刻）。比較の基準は乱数 40,000 本。採点は温度計 P2・P0 固定＋変位計1本の反り σ。",
           "time_s": {"analytic": t_an, "mc400": t400, "mc4000": t4, "mc40000": t40},
           "B150_diag": {"analytic": np.diag(Ba[4]).round(3).tolist(), "mc400": np.diag(B400[4]).round(3).tolist(),
                         "mc4000": np.diag(B4[4]).round(3).tolist(), "mc40000": np.diag(B40[4]).round(3).tolist()},
           "B_rel_diff_vs_mc40000": {"analytic": rel(Ba, B40), "mc400": rel(B400, B40), "mc4000": rel(B4, B40)},
           "selection": {k: {kk: vv for kk, vv in v.items() if kk != "scores"} for k, v in res.items()},
           "rank_corr_vs_mc40000": {k: rankcorr(res[k]["scores"], res["mc40000"]["scores"]) for k in ["analytic", "mc400", "mc4000"]},
           "top100_overlap_vs_mc40000": {k: len(top(res[k]["scores"]) & top(res["mc40000"]["scores"])) for k in ["analytic", "mc400", "mc4000"]}}
    json.dump(out, open(os.path.join(RES, "analytic_prior_covariance.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "note"}, ensure_ascii=False, indent=1))


if __name__ == "__main__":
    main()
