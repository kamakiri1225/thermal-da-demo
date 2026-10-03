"""【切り口B の実ソルバ検証】A・O を一度も測らず、別の場所の変位2点で A・O を当てられるか.

blog_008 §5 は真値も ROM の双子実験だった。ここでは真値を実ソルバにする。
  温度の真値 ＝ OpenFOAM(CHT) の固体温度場（results/limit_truth_<case>.npz）
  変位の真値 ＝ その温度場を各時刻 FrontISTR に入れて計算した全節点の変位
同化モデルは ROM のまま（観測演算子も ROM 由来の D）。
観測：温度2点(P2+P0) ＋ 変位2点（置き場所を変える）。A・O は評価専用で一度も観測しない。

出力: results/holdout_AO_realsolver.json, results/truth_disp_all_<case>.npz（FrontISTR真値のキャッシュ）
      docs/img/holdout_AO_realsolver.png
再現: OMP_NUM_THREADS=4 python3 run/holdout_AO_realsolver.py
"""
from __future__ import annotations
import os, sys, json
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from improve_heating_schedule import integrate as integ_on
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TN=[2,0]; COMP=["Ux","Uy","Uz"]
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H_=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
A_XYZ=(0.028,0.0,H_); O_XYZ=(-0.028,0.0,H_)
SEL=[((0.0097,-0.0362,0.1005),0),((-0.0346,0.0144,0.1005),2)]   # blog_008 §5 で選定した2点
BC =[((0.028,0.0,0.075),2),((-0.028,0.0,0.075),2)]
LOW=[((0.0375,0.0,0.005),2),((-0.0375,0.0,0.005),2)]
SCHED={"learned":lambda t:(t<300.0),"q25":lambda t:(t<300.0),
       "intermittent":lambda t:(t<150.0)|((t>=300.0)&(t<450.0))}
QS={"learned":1.0,"q25":25/15,"intermittent":1.0}
CASES=[("learned","15 W（ROM の係数を決めた計算）"),("q25","25 W（新しく計算した条件）"),
       ("intermittent","15 W 間欠加熱（新しく計算した条件）")]


def truth_disp_all(case,coords,node_ids,near,Tf,times):
    """OpenFOAM の温度場を各時刻 FrontISTR に入れ、全節点の変位を得る（キャッシュ）"""
    cache=os.path.join(RES,f"truth_disp_all_{case}.npz")
    if os.path.exists(cache): return np.load(cache)["u"]
    work=Path(ROOT,"openfoam","limit",f"fistr_all_{case}"); work.mkdir(parents=True,exist_ok=True)
    fistr_case.write_mesh(work,NR,NTH,NZ,R_IN,R_OUT,H_,young_modulus=MATERIAL["young_modulus_Pa"],
        poisson_ratio=MATERIAL["poisson_ratio"],density=MATERIAL["density_kg_m3"],
        thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(work)
    u=np.zeros((len(times),len(node_ids),3))
    for k,t in enumerate(times):
        fistr_case.write_cnt(work,node_ids,Tf[k][near],reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(work); d=fistr_case.read_displacement(work)
        u[k]=np.array([d[n] for n in node_ids])*1e6
        print(f"  [{case}] FrontISTR t={t:g}s",flush=True)
    np.savez_compressed(cache,u=u); return u


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; Dall=op["D"]; coords=op["coords"]
    Wall=np.einsum("nck,kj->ncj",Dall,UPp)
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H_)
    node_ids=[n for n,_ in mesh["nodes"]]
    _,near=cKDTree(Cc).query(coords)
    near_node=lambda x:int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    iA=near_node(A_XYZ); iO=near_node(O_XYZ)
    sel=[(near_node(x),c) for x,c in SEL]; bc=[(near_node(x),c) for x,c in BC]; low=[(near_node(x),c) for x,c in LOW]
    mp=mean[pod]
    out={"note":"真値＝OpenFOAMの温度場＋それをFrontISTRに入れた全節点変位。観測＝温度2点(P2+P0)＋変位2点。"
                "A・Oは一度も観測しない。誤差は加熱期30〜300秒・5 seed平均 [µm]","cases":{}}
    fig,axs=plt.subplots(1,3,figsize=(16.5,5.0))
    for j,(case,lab) in enumerate(CASES):
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz"))
        t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]
        utrue=truth_disp_all(case,coords,node_ids,near,Tf,t)      # (時刻, 節点, 3) µm
        aot=utrue[:,iA,2]-utrue[:,iO,2]
        cyc=t[1:]; on=SCHED[case]
        def run(obs,seed):
            ws=[Wall[i,c] for i,c in obs]; u0=[um[i,c] for i,c in obs]
            rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
            Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
            Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
            Rd=np.diag([SIG_T**2]*len(TN)+[SIG_U**2]*len(ws))
            wA=Wall[iA,2]; wO=Wall[iO,2]; aA=um[iA,2]; aO=um[iO,2]
            eAO=[];eA=[];rec=[np.nan]; tp=0.0
            for ci,tb in enumerate(cyc,1):
                Z=Z.copy(); Z[:,:NPT]=integ_on(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,on,tp,tb); tp=tb
                yv=list(T5[ci][TN]); Yf=Z[:,TN]
                for (i,c),w,b0 in zip(obs,ws,u0):
                    yv.append(utrue[ci,i,c]); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
                y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                m=Z[:,:NPT].mean(0); dm=m-mp
                est=(aA+wA@dm)-(aO+wO@dm); rec.append(est)
                if tb<=300:
                    eAO.append(abs(est-aot[ci])); eA.append(abs((aA+wA@dm)-utrue[ci,iA,2]))
            return np.mean(eAO),np.mean(eA),np.array(rec)
        cfgs=[("変位なし（温度2点のみ）",[]),("A・O 自身を観測（循環・参考）",[(iA,2),(iO,2)]),
              ("B・C（中段 z=75 mm）",bc),("選定2点（上面）",sel),("底面近く z=5 mm（悪い例）",low)]
        res={}; traj={}
        for nm,obs in cfgs:
            r=[run(obs,s) for s in SEEDS]
            res[nm]=dict(AO_um=float(np.mean([x[0] for x in r])),AO_sd=float(np.std([x[0] for x in r],ddof=1)),
                         A_um=float(np.mean([x[1] for x in r])))
            traj[nm]=np.mean([x[2] for x in r],axis=0)
            print(f"{case:13s} {nm:26s} 反りA−O {res[nm]['AO_um']:.3f}±{res[nm]['AO_sd']:.3f}  Uz(A) {res[nm]['A_um']:.3f} µm",flush=True)
        out["cases"][case]=res
        ax=axs[j]; ax.axvspan(0,300,color="#FDEBD0",alpha=.40)
        ax.plot(t,aot,"o-",color="k",lw=2.4,ms=4,label="真値（OpenFOAM＋FrontISTR）")
        for nm,c in [("変位なし（温度2点のみ）","#9AA5B1"),("B・C（中段 z=75 mm）","#E67E22"),("選定2点（上面）","#2E8B57")]:
            ax.plot(t,traj[nm],"-",color=c,lw=2,label=f"{nm}（{res[nm]['AO_um']:.3f} µm）")
        ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3); ax.set_title(lab,fontsize=11.5)
        if j==0: ax.set_ylabel("反り A−O [µm]"); ax.legend(fontsize=8.5,loc="upper left")
    json.dump(out,open(os.path.join(RES,"holdout_AO_realsolver.json"),"w"),ensure_ascii=False,indent=1)
    fig.suptitle("A・O を一度も観測せずに反りを当てる ― 真値＝OpenFOAM＋FrontISTR（5 seed 平均）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"holdout_AO_realsolver.png"),dpi=150); plt.close(fig)
    print("wrote holdout_AO_realsolver.png")


if __name__=="__main__": main()
