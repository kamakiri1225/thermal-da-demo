"""EnKFが「温度と発熱量Qの相関」を自動で作る様子を追う（blog_003 §4の検証）.

固定BのOIでは、温度とQを独立にサンプリングすると Cov(T,Q)=0 になり、
観測が来てもQが動かない。EnKFは前進積分するだけで相関が立ち上がる、
という差をアンサンブルの標本相関係数で可視化する。

出力: 標準出力（サイクルごとの corr(T_i, Q)）
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/enkf_correlation_demo.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES = os.path.join(ROOT, 'results')
NPT, IQ, IH, NAUG = 5, 5, 6, 7
N_ENS = 60; DT = 2.0; OBS_DT = 30.0; T_END = 600.0; SIG_T = 0.30
TSENS = [2]          # 温度センサはP2の1点だけ


def main():
    np.set_printoptions(suppress=True, precision=3, floatmode='fixed')
    d = np.load(os.path.join(RES, 'rom_calibrated_pod.npz'))
    C = d['C']; K = rg.tri_to_matrix(d['K_upper'], NPT)
    h_true = float(d['h']); heat = int(d['heat_node'])
    rng = np.random.default_rng(20260913); rng_o = np.random.default_rng(7)

    cyc = np.arange(OBS_DT, T_END + 1e-9, OBS_DT)
    T = np.full(NPT, rg.T_AIR_K); Ttr = [T.copy()]
    for a, b in zip(np.r_[0, cyc[:-1]], cyc):
        _, tr = rg.integrate_single(T, C, K, h_true, 1.0, heat, a, b, DT)
        T = tr[-1]; Ttr.append(T.copy())
    Ttr = np.array(Ttr)

    # 温度とQを「独立に」サンプリングする（ここが出発点）
    Z = np.zeros((N_ENS, NAUG))
    Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
    Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS)
    Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)

    def corr_TQ(Z):
        return np.array([np.corrcoef(Z[:, i], Z[:, IQ])[0, 1] for i in range(NPT)])

    print("サイクル0（独立サンプリング直後）: corr(T,Q) =", corr_TQ(Z))
    print("  → ほぼ0。この状態では観測が来てもQは動かせない")

    tp = 0.0
    for ci, tb in enumerate(cyc, 1):
        Z = Z.copy()
        Z[:, :NPT] = rg.integrate_ensemble(Z[:, :NPT], C, K, Z[:, IH], Z[:, IQ], heat, tp, tb, DT)
        tp = tb
        if ci in (1, 2, 5, 10):
            print(f"サイクル{ci}（t={tb:.0f}s, 更新前）: corr(T,Q) =", corr_TQ(Z))
        Rd = np.diag([SIG_T ** 2] * len(TSENS))
        y = Ttr[ci][TSENS] + rng_o.normal(0, SIG_T, len(TSENS))
        Z = enkf_update(Z, y, None, Rd, rng, inflation=1.02, Yf=Z[:, TSENS])
        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)

    print(f"\n最終: Q = {Z[:, IQ].mean()*15:.2f} W (真値15) / "
          f"h = {Z[:, IH].mean():.5f} W/K (真値{h_true:.5f}) / "
          f"温度RMSE = {np.sqrt(((Z[:, :NPT].mean(0)-Ttr[-1])**2).mean()):.4f} K")


if __name__ == '__main__':
    main()
