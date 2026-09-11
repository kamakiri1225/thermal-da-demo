"""DUMPWパッチ版FrontISTRがダンプした K, H から感度行列 W=K⁻¹H の全体を構築する.

前提: run/dumpw_cylinder.py 実行済みで openfoam/dumpw_case/ に
  dump_matrix_1_0.mm (剛性行列K, 15120x15120, MatrixMarket)
  H_matrix.mtx       (熱ひずみ荷重行列H, 15120x5040)
がある。K を splu 分解し 504列ずつブロック求解、uz行(2::3)だけ抜き出して保存する。

出力:
  results/Wz_full.npy          … Wz[5040節点, 5040温度] [m/K] float32 (~100MB, git管理外)
  results/W_full_uz_rowsens.npy … 行感度 Σ_j|Wz[i,j]|×1e6 [µm/K]
検証: Wz[qa]-Wz[qo] は DUMPW出力 Wdiff_z と最大相対差 1.9e-7 で一致(同一行列由来)。
実行: ROOT(105_0)で `OMP_NUM_THREADS=4 python3 run/build_W_from_dumps.py` (約2分)
"""
import os, sys, time
import numpy as np
import scipy.io as sio
from scipy.sparse.linalg import splu

HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
CASE=os.path.join(ROOT,"openfoam","dumpw_case"); RES=os.path.join(ROOT,"results")

t0=time.time()
K=sio.mmread(os.path.join(CASE,"dump_matrix_1_0.mm")).tocsc()
H=sio.mmread(os.path.join(CASE,"H_matrix.mtx")).tocsc()
print(f"read {time.time()-t0:.1f}s  K{K.shape} H{H.shape}",flush=True)
t0=time.time(); lu=splu(K); print(f"LU {time.time()-t0:.1f}s",flush=True)
n_t=H.shape[1]; Wz=np.zeros((n_t,n_t),dtype=np.float32)
t0=time.time()
for c0 in range(0,n_t,504):
    c1=min(c0+504,n_t)
    X=lu.solve(np.asarray(H[:,c0:c1].todense()))
    Wz[:,c0:c1]=X[2::3,:]          # uz成分のみ
    print(f"  cols {c1}/{n_t} ({time.time()-t0:.0f}s)",flush=True)
sens=np.abs(Wz.astype(np.float64)).sum(axis=1)*1e6
np.save(os.path.join(RES,"W_full_uz_rowsens.npy"),sens)
np.save(os.path.join(RES,"Wz_full.npy"),Wz)
print(f"row-sens(uz): min={sens.min():.3g} max={sens.max():.3g} µm/K",flush=True)
