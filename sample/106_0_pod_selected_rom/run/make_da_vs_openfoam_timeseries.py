"""blog_007 用：OpenFOAM の温度（真値）と、同化した温度を時刻歴で比較する.

run/limit_realsolver_truth.py と同じ設定で同化を回し、
  黒＝OpenFOAM の温度（その点のセル値、真値）
  灰＝同化なし
  緑＝温度2点のみ（P2+P4）
  赤＝温度2点＋変位2点（A/O）
を、代表5点のうち3点について3条件ぶん描く。

出力: docs/img/da_vs_openfoam_timeseries.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_da_vs_openfoam_timeseries.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
CASES=[("learned","15 W・0〜300 s（モデル作成に使った条件）"),
       ("q25","25 W・0〜300 s（使っていない条件）"),
       ("intermittent","15 W 間欠加熱（使っていない条件）")]
SHOW=[(2,"P2（ヒータ横・観測点）"),(1,"P1（反対側・未観測）"),(3,"P3（底面・未観測）")]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"disp_operator.npz")); ub=op["uz_mean"]; Dm=op["Dmode"]
    disp=lambda T5: ub+((T5-mean[pod])@UPp.T)@Dm.T
    fig,axs=plt.subplots(len(SHOW),len(CASES),figsize=(16,10.5),sharex=True)
    for j,(case,lab) in enumerate(CASES):
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]; uz=z["uz"]
        cyc=t[1:]
        def run(nodes,use_disp,seed):
            rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
            Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
            Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
            if nodes is not None:
                Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
            rec=[Z[:,:NPT].mean(0).copy()]; tp=0.0
            for ci,tb in enumerate(cyc,1):
                Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
                if nodes is not None:
                    yv=list(T5[ci][nodes]); Yf=Z[:,nodes]
                    if use_disp: yv+=list(uz[ci]); Yf=np.column_stack([Yf,disp(Z[:,:NPT])])
                    y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
                    Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                    Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                rec.append(Z[:,:NPT].mean(0).copy())
            return np.array(rec)
        runs={}
        for nm,nd,ud in [("同化なし",None,False),("温度2点（P2+P4）",[2,4],False),("温度2点＋変位2点",[2,4],True)]:
            runs[nm]=np.mean([run(nd,ud,s) for s in SEEDS],axis=0)
        print("[da-ts]",case,"done",flush=True)
        for i,(pi,pname) in enumerate(SHOW):
            ax=axs[i,j]; ax.axvspan(0,300,color="#FDEBD0",alpha=.40)
            ax.plot(t,T5[:,pi]-273.15,"o-",color="k",lw=2.6,ms=4,label="OpenFOAM（真値）")
            for nm,col,ls in [("同化なし","#9AA5B1","--"),("温度2点（P2+P4）","#2E8B57","-"),("温度2点＋変位2点","#C0392B","-")]:
                ax.plot(t,runs[nm][:,pi]-273.15,ls,color=col,lw=1.9,label=nm)
            ax.grid(alpha=.3)
            if j==0: ax.set_ylabel(f"{pname}\n温度 [℃]",fontsize=10.5)
            if i==0: ax.set_title(lab,fontsize=12)
            if i==len(SHOW)-1: ax.set_xlabel("時刻 [s]")
            ax.set_ylim(16,32 if case=="q25" else 28)
    axs[0,0].legend(fontsize=9,loc="lower right")
    fig.suptitle("OpenFOAM の温度（黒）と、同化した温度の時刻歴（5 seed 平均）\n観測は温度2点 P2+P4（＋変位2点 A/O）。P1・P3 は観測していない点",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig(os.path.join(IMG,"da_vs_openfoam_timeseries.png"),dpi=150); plt.close(fig)
    print("wrote da_vs_openfoam_timeseries.png")


if __name__=="__main__": main()
