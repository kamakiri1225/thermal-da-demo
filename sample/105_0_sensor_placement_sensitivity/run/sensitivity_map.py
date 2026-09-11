"""変位観測2点の「熱感度分布」∂uz/∂T(x) を FrontISTR で実計算する.

どの場所の温度が上がると、観測している上面変位(uz_heater / uz_opp)が
どれだけ動くか ― を円筒全体の分布として求める。方法はパッチ摂動:
メッシュを周方向 NTH_P × 軸方向 NZ_P のパッチに分け、各パッチの節点に
+1 K を与えて FrontISTR を実行し、2観測点の変位応答を記録する。
線形なのでこれが厳密な感度(パッチ内一様摂動に対する)になる。

出力:
  results/sensitivity_uz.npz      … パッチ感度 (2, NTH_P, NZ_P) [µm/K] とメタ
  docs/img/sensitivity_map.png    … 円筒面上の感度分布(2パネル)
使い方:
  OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensitivity_map.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np

HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL

NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
NTH_P, NZ_P = 12, 10                    # パッチ分割(周方向×軸方向)
OBS_XYZ=np.array([[0.028,0,H],[-0.028,0,H]])   # uz_heater / uz_opp (上面)

def main():
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    theta=np.arctan2(coords[:,1],coords[:,0])%(2*np.pi)
    z=coords[:,2]
    ith=np.minimum((theta/(2*np.pi)*NTH_P).astype(int),NTH_P-1)
    iz=np.minimum((z/H*NZ_P).astype(int),NZ_P-1)
    oi=[int(np.linalg.norm(coords-p,axis=1).argmin()) for p in OBS_XYZ]

    work=os.path.join(ROOT,"openfoam","sens_probe"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))

    S=np.zeros((2,NTH_P,NZ_P))          # 感度 [m/K→後でµm]
    npatch=np.zeros((NTH_P,NZ_P),dtype=int)
    k=0
    for a in range(NTH_P):
        for b in range(NZ_P):
            sel=(ith==a)&(iz==b)
            npatch[a,b]=sel.sum()
            T=np.full(len(coords),Tref); T[sel]+=1.0     # このパッチだけ+1K
            fistr_case.write_cnt(Path(work),node_ids,T,reference_temperature=Tref,
                young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
                thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
            fistr_case.run_fistr(Path(work))
            disp=fistr_case.read_displacement(Path(work))
            for j in range(2): S[j,a,b]=disp[node_ids[oi[j]]][2]
            k+=1
            if k%12==0: print(f"[sens] {k}/{NTH_P*NZ_P} patches",flush=True)
    S_um=S*1e6
    os.makedirs(os.path.join(ROOT,"results"),exist_ok=True)
    np.savez(os.path.join(ROOT,"results","sensitivity_uz.npz"),
             S_um=S_um,npatch=npatch,NTH_P=NTH_P,NZ_P=NZ_P,
             note="S_um[j,a,b]=パッチ(a=theta,b=z)を+1KしたときのUz応答[µm]。j=0:heater側,1:opp側")
    print("[sens] wrote results/sensitivity_uz.npz")
    print(f"[sens] uz_heater感度: min={S_um[0].min():.4f} max={S_um[0].max():.4f} µm/K")
    print(f"[sens] uz_opp   感度: min={S_um[1].min():.4f} max={S_um[1].max():.4f} µm/K")

    # ---- 可視化: 円筒展開図(θ-z) 2パネル ----
    from dacore import plots as _p  # フォント
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(13,5),sharey=True)
    ext=[0,360,0,H*1000]
    vmax=max(abs(S_um).max(),1e-9)
    for j,(ax,name) in enumerate(zip(axs,["uz_heater (ヒータ側上面)","uz_opp (反対側上面)"])):
        im=ax.imshow(S_um[j].T,origin="lower",aspect="auto",extent=ext,
                     cmap="coolwarm",vmin=-vmax,vmax=vmax)
        ax.set_title(f"∂{name.split()[0]}/∂T の分布")
        ax.set_xlabel("周方向角度 θ [deg]  (0°=ヒータ中心)")
        if j==0: ax.set_ylabel("高さ z [mm]")
        ax.axvline(0,color='orange',lw=2,alpha=0.6); ax.axvline(360,color='orange',lw=2,alpha=0.6)
        ax.text(5,H*1000*0.94,"↑ヒータ側",color='darkorange',fontsize=12)
        ax.text(175,H*1000*0.94,"反対側",color='navy',fontsize=12)
        fig.colorbar(im,ax=ax,label="感度 [µm/K]")
    fig.suptitle("変位観測2点の熱感度分布（パッチ+1K → FrontISTR、円筒展開図）",fontsize=17)
    fig.tight_layout()
    out=os.path.join(ROOT,"docs","img","sensitivity_map.png")
    fig.savefig(out,dpi=140)
    import subprocess; subprocess.run(["cp",out,os.path.join(ROOT,"presentation","assets","sensitivity_map.png")],check=False)
    print(f"[sens] wrote {os.path.relpath(out,ROOT)}")

if __name__=="__main__": main()
