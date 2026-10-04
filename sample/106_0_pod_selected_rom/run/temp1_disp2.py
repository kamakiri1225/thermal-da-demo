"""温度1点＋変位2点でもデータ同化でき、精度が上がるかを確かめる（blog_008 §0-8）.

観測の組み合わせを変えて、同じ EnKF を回す。真値は OpenFOAM（温度）＋FrontISTR（変位）。
  温度1点（P2）／温度1点（P2）＋変位2点（選定 S1・S2）／温度2点（P2・P4）／温度2点＋変位2点
評価：全 20,696 セルの温度 RMSE、観測していない反り A−O、推定した発熱量（加熱期 30〜300 秒、5 seed 平均）
条件：15 W・25 W（0〜300 秒）、15 W 間欠加熱（ヒータの ON/OFF を ROM に与える）

出力: results/temp1_disp2.json, docs/img/temp1_disp2.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/temp1_disp2.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from improve_heating_schedule import integrate as integ_on
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
SCHED={"learned":lambda t:(t<300.0),"q25":lambda t:(t<300.0),
       "intermittent":lambda t:(t<150.0)|((t>=300.0)&(t<450.0))}
QTRUE={"learned":15.0,"q25":25.0,"intermittent":15.0}
CASES=[("learned","15 W"),("q25","25 W"),("intermittent","15 W 間欠加熱")]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod]); mp=mean[pod]; Af=U@UPp
    op=np.load(os.path.join(RES,"dispop_allnodes.npz")); um=op["u_mean"]; co=op["coords"]
    Wall=np.einsum("nck,kj->ncj",op["D"],UPp)
    near=lambda x:int(np.linalg.norm(co-np.array(x),axis=1).argmin())
    iA,iO=near((0.028,0,0.1005)),near((-0.028,0,0.1005))
    sel=[(near((0.0097,-0.0362,0.1005)),0),(near((-0.0346,0.0144,0.1005)),2)]
    cfgs=[("温度1点（P2）",[2],[]),("温度1点（P2）＋変位2点",[2],sel),
          ("温度2点（P2・P4）",[2,4],[]),("温度2点（P2・P4）＋変位2点",[2,4],sel)]
    out={"note":"真値＝OpenFOAM＋FrontISTR。変位2点は blog_008 §5 で選んだ上面の S1（Ux）・S2（Uz）。"
                "評価は加熱期30〜300秒・5 seed平均。間欠加熱はヒータの ON/OFF を ROM に与えた","cases":{}}
    for case,lab in CASES:
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]
        ut=np.load(os.path.join(RES,f"truth_disp_all_{case}.npz"))["u"]; aot=ut[:,iA,2]-ut[:,iO,2]
        on=SCHED[case]; res={}
        for nm,tn,obs in cfgs:
            ws=[Wall[i,c] for i,c in obs]; u0=[um[i,c] for i,c in obs]; r=[]
            for seed in SEEDS:
                rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
                Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
                Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
                Rd=np.diag([SIG_T**2]*len(tn)+[SIG_U**2]*len(ws))
                eF=[];eAO=[];eQ=[]; tp=0.0
                for ci,tb in enumerate(t[1:],1):
                    Z=Z.copy(); Z[:,:NPT]=integ_on(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,on,tp,tb); tp=tb
                    yv=list(T5[ci][tn]); Yf=Z[:,tn]
                    for (i,c),w,b0 in zip(obs,ws,u0):
                        yv.append(ut[ci,i,c]); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
                    y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
                    Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                    Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                    if tb<=300:
                        m=Z[:,:NPT].mean(0); dm=m-mp
                        eF.append(np.sqrt(((mean+Af@dm-Tf[ci])**2).mean()))
                        eAO.append(abs((um[iA,2]+Wall[iA,2]@dm)-(um[iO,2]+Wall[iO,2]@dm)-aot[ci]))
                        eQ.append(15*Z[:,IQ].mean())
                r.append([np.mean(eF),np.mean(eAO),eQ[-1]])
            r=np.array(r); res[nm]=dict(field_K=float(r[:,0].mean()),AO_um=float(r[:,1].mean()),Q_W_at300=float(r[:,2].mean()))
            print(f"{lab:10s} {nm:22s} 温度場 {res[nm]['field_K']:.3f} K  反り {res[nm]['AO_um']:.3f} µm  Q(300s) {res[nm]['Q_W_at300']:.2f} W",flush=True)
        out["cases"][case]=dict(label=lab,Q_true=QTRUE[case],configs=res)
    json.dump(out,open(os.path.join(RES,"temp1_disp2.json"),"w"),ensure_ascii=False,indent=1)
    names=[c[0] for c in cfgs]; cols=["#9AA5B1","#C0392B","#5D6D7E","#7B241C"]
    fig,axs=plt.subplots(1,2,figsize=(15,5.8))
    x=np.arange(3); w=.2
    for ax,key,ttl,un in [(axs[0],"field_K","全 20,696 セルの温度の誤差","K"),(axs[1],"AO_um","観測していない反り A−O の誤差","µm")]:
        for k,nm in enumerate(names):
            v=[out["cases"][c]["configs"][nm][key] for c,_ in CASES]
            ax.bar(x+(k-1.5)*w,v,w,color=cols[k],label=nm)
            for i,val in enumerate(v): ax.text(x[i]+(k-1.5)*w,val,f"{val:.2f}",ha="center",va="bottom",fontsize=10)
        ax.set_xticks(x); ax.set_xticklabels([l for _,l in CASES]); ax.set_ylabel(f"誤差 [{un}]")
        ax.set_title(ttl,fontsize=14); ax.grid(axis="y",alpha=.3)
    axs[0].legend(fontsize=11)
    fig.suptitle("温度1点に変位2点を足すと、どれだけ良くなるか（真値＝OpenFOAM＋FrontISTR、5 seed 平均、加熱期 30〜300 秒）",fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"temp1_disp2.png"),dpi=140); plt.close(fig)
    print("wrote temp1_disp2.png")


if __name__=="__main__": main()
