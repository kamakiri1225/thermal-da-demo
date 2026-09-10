"""ROM 5点温度 → 上面変位2点 の線形演算子 M を構築する（IDW→FrontISTR）。

FrontISTR は線形静解析なので、上面変位は温度場の線形汎関数。ROM の5ノード温度から
IDW でメッシュ節点温度を作り FrontISTR に通す一連は、5つの単位温度モードの応答を
並べた行列 M（2×5, uz=M(T_node - T_ref)）で厳密に表せる。M[:,i] はノード i を +1K
にしたときの上面2点の Uz。

出力: config/M_frontistr_operator.npy（単位: m/K）
使い方: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/compute_M_operator.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, fistr_case
from dacore.node_locations import NODE_XYZ
from fem.fem_obs import MATERIAL

NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
OBS_XYZ=np.array([[0.028,0,0.1005],[-0.028,0,0.1005]])  # 上面 heater側 / opposite側

def main():
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    romn=np.array(list(NODE_XYZ.values()))
    W=np.zeros((len(coords),5))
    for i,p in enumerate(coords):
        d=np.linalg.norm(romn-p,axis=1)
        if d.min()<1e-9: W[i,d.argmin()]=1
        else: w=1/d**2; W[i]=w/w.sum()
    oi=[int(np.linalg.norm(coords-p,axis=1).argmin()) for p in OBS_XYZ]
    work=os.path.join(ROOT,"openfoam","M_probe"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    M=np.zeros((2,5))
    for i in range(5):
        Tmesh=Tref+(W@np.eye(5)[i])   # ノードiを+1K
        fistr_case.write_cnt(Path(work),node_ids,Tmesh,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        M[:,i]=[disp[node_ids[j]][2] for j in oi]
        print(f"[M] node {i}: uz応答 = {np.round(M[:,i]*1e6,3)} µm/K",flush=True)
    out=os.path.join(ROOT,"config","M_frontistr_operator.npy"); np.save(out,M)
    print(f"[M] wrote {os.path.relpath(out,ROOT)} (m/K)")

if __name__=="__main__": main()
