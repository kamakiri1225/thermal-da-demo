"""観測に使わない「別の2点間変位差（未観測QoI）」をDA後の状態から評価する.

観測に使う量：温度2点（高 dT/dQ）＋ 変位2点＝上面A/O(+X/-X, z=100.5mm, W選定)。
評価する量（観測に使わない）：中高さ z=75mm の +X/-X 上向き変位差 Uz(B)-Uz(C)。
  → 同じ非対称曲げの物理だが、センサを置いていない別位置。

手順：
 1. ROM EnKF で「温度2点＋変位2点(A/O)」を同化し、解析後の5点温度履歴 recT を得る。
 2. 未観測対 B/C の変位演算子を FrontISTR で作る（平均場＋5モードの6回）。
 3. DA後の recT→gappy-POD全場復元→B/C変位差 を計算し、真値・同化なしと比較。
→ 観測したQoIだけでなく未観測QoIも真値へ寄れば、同化後の温度場が大域的に正しい
  ＝デジタルツインとして任意箇所の熱変形を予測できる、の検証になる。

出力: docs/img/oi_unobserved_qoi.png（106へ）
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/eval_unobserved_qoi.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
K=273.15; NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
# 観測する上面A/O（+X/-X, z=H）と、観測しない中高さ対 B/C（+X/-X, z=0.75H）
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])
ZUN=0.075
B_XYZ=np.array([0.028,0,ZUN]); C_XYZ=np.array([-0.028,0,ZUN])


def build_op_for(targets, U, mean, Cc, tag):
    """targets(list of xyz) の上向き変位演算子を作る。戻り uz_mean(n,), Dmode(n,r) [µm]。"""
    cache=os.path.join(RES,f"dispop_{tag}.npz")
    if os.path.exists(cache):
        d=np.load(cache); return d["uz_mean"], d["Dmode"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    qidx=[int(np.linalg.norm(coords-t,axis=1).argmin()) for t in targets]
    for t,q in zip(targets,qidx):
        print(f"[unobs] target {np.round(t*1000,1)}mm -> node {node_ids[q]} at {np.round(coords[q]*1000,1)}mm")
    _,near=cKDTree(Cc).query(coords)
    work=os.path.join(ROOT,"openfoam",f"dispop_{tag}"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    def run_field(Tfield):
        fistr_case.write_cnt(Path(work),node_ids,Tfield,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        return np.array([disp[node_ids[q]][2] for q in qidx])*1e6   # µm, 上向きUz
    r=U.shape[1]
    uz_mean=run_field(mean[near])
    Dmode=np.zeros((len(targets),r))
    for k in range(r):
        Dmode[:,k]=run_field(Tref+U[near,k]); print(f"[unobs] {tag} mode {k+1}/{r} FrontISTR",flush=True)
    np.savez(cache,uz_mean=uz_mean,Dmode=Dmode)
    return uz_mean, Dmode


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    xyz=d["xyz"]; heat_node=int(d["heat_node"]); C=d["C"]
    Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    pod_cells=kv["cell_idx"]; UP_pinv=np.linalg.pinv(U[pod_cells,:])

    # 観測用A/O演算子（既存キャッシュ）と、未観測B/C演算子（新規）
    do=np.load(os.path.join(RES,"disp_operator.npz")); uzA,DmodeA=do["uz_mean"],do["Dmode"]
    uzB,DmodeB=build_op_for([B_XYZ,C_XYZ],U,mean,Cc,"unobs_z75")

    def a_of(T5): return (T5-mean[pod_cells])@UP_pinv.T
    def dispAO(T5): a=a_of(T5); return uzA+a@DmodeA.T            # 観測QoI用(上面A/O)
    def qoiAO(T5):  u=dispAO(T5); return u[...,0]-u[...,1]
    def qoiBC(T5):  a=a_of(T5); u=uzB+a@DmodeB.T; return u[...,0]-u[...,1]  # 未観測QoI(中高さB/C)

    # dT/dQ で温度2点を選定
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens=(pert-base)/0.1; hi2=list(np.argsort(sens)[::-1][:2])
    print("[unobs] dT/dQ=",np.round(sens,2)," 温度観測点(高dT/dQ)=P",hi2)

    # 真値トラジェクトリ
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tgrid=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(nodes,use_disp,seed):
        """解析後の5点温度履歴 recT[nt,5] を返す（nodes=Noneで同化なし）。"""
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        if nodes is not None:
            Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
        recT=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Tn=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
            Z=Z.copy(); Z[:,:NPT]=Tn; tp=tb
            if nodes is not None:
                yv=list(Ttr[ci][nodes]); Yf=Z[:,nodes]
                if use_disp:
                    yv=yv+list(dispAO(Ttr[ci])); Yf=np.column_stack([Yf,dispAO(Z[:,:NPT])])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            recT.append(Z[:,:NPT].mean(0).copy())
        return np.array(recT)

    # 未観測QoI(B/C) と 観測QoI(A/O) を、真値/同化なし/同化(温度2+変位A/O) で評価（5seed平均）
    recT_free=np.mean([run(None,False,s) for s in SEEDS],axis=0)
    recT_da  =np.mean([run(hi2,True ,s) for s in SEEDS],axis=0)
    bc_true=qoiBC(Ttr); bc_free=qoiBC(recT_free); bc_da=qoiBC(recT_da)
    ao_true=qoiAO(Ttr); ao_free=qoiAO(recT_free); ao_da=qoiAO(recT_da)
    ht=(tgrid>0)&(tgrid<=300)
    eB_free=np.abs(bc_free-bc_true)[ht].mean(); eB_da=np.abs(bc_da-bc_true)[ht].mean()
    eA_free=np.abs(ao_free-ao_true)[ht].mean(); eA_da=np.abs(ao_da-ao_true)[ht].mean()
    print(f"[unobs] 未観測QoI(B/C z75) 加熱期平均誤差[µm]: 同化なし={eB_free:.2f} 同化後={eB_da:.2f}")
    print(f"[unobs] 観測QoI(A/O top)   加熱期平均誤差[µm]: 同化なし={eA_free:.2f} 同化後={eA_da:.2f}")

    fig,axes=plt.subplots(1,2,figsize=(14.5,5.8))
    for ax,(ttl,qt,qf,qd,ef,ed,obs) in zip(axes,[
        ("未観測QoI：中高さ z=75mm の +X/−X 上向きUz差（センサ無し）",bc_true,bc_free,bc_da,eB_free,eB_da,False),
        ("観測QoI：上面 A/O(+X/−X, z=100.5mm) の Uz差（変位センサあり）",ao_true,ao_free,ao_da,eA_free,eA_da,True)]):
        ax.axvspan(0,300,color="orange",alpha=0.06)
        ax.plot(tgrid,qt,"-",color="k",lw=4,alpha=0.35,label="真値")
        ax.plot(tgrid,qf,":",color="tab:gray",lw=2.2,label=f"同化なし  誤差{ef:.1f}µm")
        ax.plot(tgrid,qd,"--",color="tab:red",lw=2.2,marker="o",ms=4,label=f"温度2点+変位A/Oで同化  誤差{ed:.1f}µm")
        ax.set_title(ttl,fontsize=11,color=("black" if obs else "dimgray"))
        ax.set_xlabel("time [s]"); ax.set_ylabel("2点間 上向きUz差 [µm]"); ax.grid(alpha=0.3); ax.legend(fontsize=9)
    fig.suptitle("観測に使わない別の変位差(未観測QoI)もDA後に真値へ寄るか ― ROMで検証\n"
                 "温度2点(高dT/dQ)+変位A/Oのみ同化。未観測の中高さ変位差(左)を、同化後の全場復元から評価。",
                 fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.92])
    out=os.path.join(IMG,"oi_unobserved_qoi.png"); fig.savefig(out,dpi=140); plt.close(fig)
    import json
    json.dump({"unobserved_pair_z_mm":ZUN*1000,"temp_obs_nodes":list(map(int,hi2)),
               "unobs_QoI_heat_MAE_um":{"no_DA":float(eB_free),"DA":float(eB_da)},
               "obs_QoI_heat_MAE_um":{"no_DA":float(eA_free),"DA":float(eA_da)}},
              open(os.path.join(RES,"unobserved_qoi.json"),"w"),ensure_ascii=False,indent=2)
    print("[unobs] wrote docs/img/oi_unobserved_qoi.png, results/unobserved_qoi.json")


if __name__=="__main__": main()
