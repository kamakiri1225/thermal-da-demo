"""全FEM節点について 温度→変位 の写像を作る（FrontISTR 6回のみ）.

disp_operator.npz は上面A/Oの2点ぶんしか持っていないが、FrontISTR は1回の解析で
全節点の変位を返す。したがって同じ6回（平均場＋PODモード5本）で、全節点×3成分の
    u = u_mean + D @ a,   a = U_rP^+ (T_5 - mean_P)
を作れる。これで任意の候補点の変位を追加FEMなしで評価できる。

出力: results/dispop_allnodes.npz  (u_mean (n,3), D (n,3,5), coords (n,3), node_ids)
再現: OMP_NUM_THREADS=4 python3 run/build_full_disp_operator.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results")
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
OUT=os.path.join(RES,"dispop_allnodes.npz")


def main():
    if os.path.exists(OUT): print("already exists:",OUT); return
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([x for _n,x in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    _,near=cKDTree(Cc).query(coords)
    work=Path(ROOT,"openfoam","dispop_all"); work.mkdir(parents=True,exist_ok=True)
    fistr_case.write_mesh(work,NR,NTH,NZ,R_IN,R_OUT,H,young_modulus=MATERIAL["young_modulus_Pa"],
        poisson_ratio=MATERIAL["poisson_ratio"],density=MATERIAL["density_kg_m3"],
        thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(work)
    def run_field(Tfield):
        fistr_case.write_cnt(work,node_ids,Tfield,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(work); disp=fistr_case.read_displacement(work)
        return np.array([disp[n] for n in node_ids])*1e6      # (n,3) µm
    u_mean=run_field(mean[near]); print("[all] mean field done",flush=True)
    r=U.shape[1]; D=np.zeros((len(node_ids),3,r))
    for k in range(r):
        D[:,:,k]=run_field(Tref+U[near,k]); print(f"[all] mode {k+1}/{r}",flush=True)
    np.savez_compressed(OUT,u_mean=u_mean,D=D,coords=coords,node_ids=np.array(node_ids))
    print("wrote",OUT,u_mean.shape,D.shape)


if __name__=="__main__": main()
