"""OpenCAE用：温度1点→温度2点→温度2点＋変位2点 でEnKF推定精度が上がるのを示す.

さらに変位2点を「高熱感度W」と「低熱感度W」で置いた差も示す。
ROM(POD選定5点)＋拡大状態EnKF。変位オペレータ D はキャッシュ(dispop_*.npz)を利用。

構成:
  1. 温度1点（ヒータ側の高感度点）
  2. 温度2点（ヒータ側＋反対側）
  3. 温度2点＋変位2点（低熱感度W）
  4. 温度2点＋変位2点（高熱感度W）
評価: 全5点温度の加熱期RMSE、未観測の変位差 Uz(A)-Uz(O) の誤差。

出力: docs/img/opencae_sensor_progression.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/opencae_sensor_progression.py
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
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    hi=np.load(os.path.join(RES,"dispop_hiW.npz")); lo=np.load(os.path.join(RES,"dispop_loW.npz"))
    do=np.load(os.path.join(RES,"disp_operator.npz")); uzAO=do["uz_mean"]; Dao=do["Dmode"]

    # 温度センサ: dT/dQ 最大(ヒータ側)と反対側(最小)の2点
    b0=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat,0,300,DT)[1][-1]
    b1=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat,0,300,DT)[1][-1]
    sQ=(b1-b0)/0.1; T1=int(np.argmax(sQ)); T2=int(np.argmin(sQ))   # 高感度・反対側
    print(f"[prog] 温度センサ P{T1}(高感度) と P{T2}(反対側)")

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    def dispA(T5,uzm,D): a=(T5-mean[pod])@UPp.T; return uzm+a@D.T
    def qoiAO(T5): a=(T5-mean[pod])@UPp.T; u=uzAO+a@Dao.T; return u[...,0]-u[...,1]

    def run(tsens, dispkind, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        m=0 if dispkind is None else 2
        Rd=np.diag([SIG_T**2]*len(tsens)+[SIG_U**2]*m)
        rmse=[np.sqrt(((Z[:,:NPT].mean(0)-Ttr[0])**2).mean())]; recT=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][tsens]); Yf=Z[:,tsens]
            if dispkind is not None:
                op=hi if dispkind=="hi" else lo
                yv+=list(dispA(Ttr[ci],op["uz_mean"],op["D"]))
                Yf=np.column_stack([Yf,dispA(Z[:,:NPT],op["uz_mean"],op["D"])])
            y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            rmse.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean())); recT.append(Z[:,:NPT].mean(0).copy())
        return np.array(rmse),np.array(recT)

    cfgs=[("温度1点",       [T1],      None, "tab:gray"),
          ("温度2点",       [T1,T2],   None, "tab:green"),
          ("温度2点＋変位2点(低W)",[T1,T2],"lo","tab:orange"),
          ("温度2点＋変位2点(高W)",[T1,T2],"hi","tab:blue")]
    ht=(tg>0)&(tg<=300); qoi_true=qoiAO(Ttr)
    res={}
    for name,ts,dk,col in cfgs:
        rs=[run(ts,dk,s) for s in SEEDS]
        rmse=np.mean([r[0] for r in rs],axis=0); recT=np.mean([r[1] for r in rs],axis=0)
        eT=rmse[ht].mean(); eU=np.abs(qoiAO(recT)-qoi_true)[ht].mean()
        res[name]=dict(rmse=rmse,eT=eT,eU=eU,col=col,recT=recT)
        print(f"[prog] {name:20s} 温度RMSE={eT:.3f}K  変位差誤差={eU:.3f}µm")

    names=[c[0] for c in cfgs]; cols=[c[3] for c in cfgs]
    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(13.5,5.6))
    b=ax0.bar(range(4),[res[n]["eT"] for n in names],color=cols)
    for bar,n in zip(b,names): ax0.text(bar.get_x()+bar.get_width()/2,res[n]["eT"],f"{res[n]['eT']:.3f}",ha="center",va="bottom",fontsize=10,weight="bold")
    ax0.set_xticks(range(4)); ax0.set_xticklabels(names,fontsize=9.5,rotation=8)
    ax0.set_ylabel("加熱期 全5点温度RMSE [K]"); ax0.set_title("温度の推定精度",fontsize=13,weight="bold"); ax0.grid(alpha=.3,axis="y")
    b=ax1.bar(range(4),[res[n]["eU"] for n in names],color=cols)
    for bar,n in zip(b,names): ax1.text(bar.get_x()+bar.get_width()/2,res[n]["eU"],f"{res[n]['eU']:.3f}",ha="center",va="bottom",fontsize=10,weight="bold")
    ax1.set_xticks(range(4)); ax1.set_xticklabels(names,fontsize=9.5,rotation=8)
    ax1.set_ylabel("加熱期 変位差 Uz(A)-Uz(O) 誤差 [µm]"); ax1.set_title("未観測の変位差の推定精度",fontsize=13,weight="bold"); ax1.grid(alpha=.3,axis="y")
    fig.suptitle("観測を増やすと精度が上がる（ROM＋EnKF, 5seed平均）：温度1点→2点→＋変位。\n"
                 "変位は『高熱感度W』の点に置くほど効く（低Wはあまり効かない）",fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.9])
    out=os.path.join(IMG,"opencae_sensor_progression.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[prog] wrote",out)


if __name__=="__main__": main()
