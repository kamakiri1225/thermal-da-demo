"""間欠加熱で同化が合わなかった問題の改善を試す.

原因：ROM のヒータが「0〜300 s に一定の発熱」と決め打ち（rom_general.heater_power_W）。
  基準     … いまの方法（0〜300 s だけ q×15 W）
  方法1    … ヒータの ON/OFF 時刻を与える（NC が運転スケジュールを知っている想定）。大きさ q は推定
  方法2    … ON/OFF も分からない想定。発熱は常に q×15 W とし、q を毎回少しずつ変わってよい値（ランダムウォーク）として推定
           予報のたびに各メンバーの q に N(0, σ_rw²) を足す（σ_rw=0.1, 0.2, 0.4。q=1 が 15 W）
真値：OpenFOAM＋FrontISTR（results/limit_truth_<case>.npz）。同化の他の設定は run/limit_realsolver_truth.py と同じ。
温度センサ P2+P4、変位計 A/O。5 seed 平均。

出力: results/improve_heating_schedule.json, docs/img/improve_*.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/improve_heating_schedule.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]; NODES=[2,4]
SCHED={"learned":lambda t:(t<300.0),"q25":lambda t:(t<300.0),"intermittent":lambda t:(t<150.0)|((t>=300.0)&(t<450.0))}
QTRUE={"learned":lambda t:15.0*(t<300),"q25":lambda t:25.0*(t<300),"intermittent":lambda t:15.0*((t<150)|((t>=300)&(t<450)))}
LAB={"learned":"15 W・0〜300 s（ROM の係数を決めた計算）","q25":"25 W・0〜300 s（新しく計算した条件）","intermittent":"15 W 間欠加熱（新しく計算した条件）"}
METHODS=[("基準（0〜300 sに一定と仮定）","fixed",0.0),("方法1：ON/OFF時刻を与える","known",0.0),
         ("方法2：発熱を時間で変わる値に（σ=0.1）","rw",0.1),("方法2：発熱を時間で変わる値に（σ=0.2）","rw",0.2),("方法2：発熱を時間で変わる値に（σ=0.4）","rw",0.4)]


def integrate(T,C,Km,h,q,heat,on,t0,t1):
    n_ens,n=T.shape; Ks=Km.sum(1); M=np.broadcast_to(Km,(n_ens,n,n))/C[None,:,None]; M=M.copy()
    idx=np.arange(n); M[:,idx,idx]=-(Ks[None,:]+h[:,None])/C[None,:]; b=h[:,None]*rg.T_AIR_K/C[None,:]
    steps=max(1,int(round((t1-t0)/DT))); dt=(t1-t0)/steps
    def f(X,t):
        H=np.zeros_like(X); H[:,heat]=float(on(t))*q*rg.HEATER_RATED_W/C[heat]
        return np.einsum("eij,ej->ei",M,X)+b+H
    t=t0
    for _ in range(steps):
        k1=f(T,t);k2=f(T+.5*dt*k1,t+.5*dt);k3=f(T+.5*dt*k2,t+.5*dt);k4=f(T+dt*k3,t+dt); T=T+dt/6*(k1+2*k2+2*k3+k4); t+=dt
    return T


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz")); C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz")); U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); op=np.load(os.path.join(RES,"disp_operator.npz")); ub=op["uz_mean"]; D=op["Dmode"]
    disp=lambda X: ub+((X-mean[pod])@UPp.T)@D.T; field=lambda X: mean+U@(UPp@(X-mean[pod]))
    reg=json.load(open(os.path.join(RES,"baseline_regression_vs_da.json")))["cases"]
    out={"note":"温度P2+P4＋変位A/O。誤差は30〜600 s平均、5 seed平均。Q_mae は推定Qと真の発熱（ON/OFF込み）の差の平均","cases":{}}; traj={}
    for case in ["learned","q25","intermittent"]:
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]; uz=z["uz"]; aot=uz[:,0]-uz[:,1]
        res={}; traj[case]={"t":t,"truth":aot}
        for name,kind,srw in METHODS:
            on={"fixed":lambda tt:(tt<300.0),"known":SCHED[case],"rw":lambda tt:True}[kind]
            eAO=[];eF=[];eQ=[];AOs=[];Qs=[]
            for seed in SEEDS:
                rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7); rq=np.random.default_rng(seed+99)
                Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
                Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
                Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2); tp=0.0; a=[];f_=[];q_=[];AO=[np.nan];Q=[15*Z[:,IQ].mean()]
                for ci,tb in enumerate(t[1:],1):
                    Z=Z.copy(); Z[:,:NPT]=integrate(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,on,tp,tb); tp=tb
                    if kind=="rw": Z[:,IQ]=np.clip(Z[:,IQ]+rq.normal(0,srw,N_ENS),0,3)
                    Yf=np.column_stack([Z[:,NODES],disp(Z[:,:NPT])])
                    y=np.r_[T5[ci][NODES],uz[ci]]+ro.normal(0,np.sqrt(np.diag(Rd)))
                    Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                    Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                    m=Z[:,:NPT].mean(0); u=disp(m); qW=15*Z[:,IQ].mean()*float(on(tb-1e-6))
                    a.append(abs((u[0]-u[1])-aot[ci])); f_.append(np.sqrt(((field(m)-Tf[ci])**2).mean())); q_.append(abs(qW-QTRUE[case](tb-1e-6)))
                    AO.append(u[0]-u[1]); Q.append(qW)
                eAO.append(np.mean(a)); eF.append(np.mean(f_)); eQ.append(np.mean(q_)); AOs.append(AO); Qs.append(Q)
            res[name]=dict(AO_um=float(np.mean(eAO)),AO_um_sd=float(np.std(eAO,ddof=1)),field_K=float(np.mean(eF)),Q_mae_W=float(np.mean(eQ)))
            traj[case][name]=(np.nanmean(AOs,0),np.mean(Qs,0))
            print(f"{case:12s} {name:28s} 反り {res[name]['AO_um']:.3f}±{res[name]['AO_um_sd']:.3f} µm  温度場 {res[name]['field_K']:.3f} K  発熱の誤差 {res[name]['Q_mae_W']:.2f} W",flush=True)
        res["回帰式（温度P2+P4）"]=dict(AO_um=reg[case]["P2+P4"]["regression_all"])
        out["cases"][case]=res
    json.dump(out,open(os.path.join(RES,"improve_heating_schedule.json"),"w"),ensure_ascii=False,indent=1)
    best_rw=min([n for n,k,_ in METHODS if k=="rw"],key=lambda n:out["cases"]["intermittent"][n]["AO_um"])
    show=[("基準（0〜300 sに一定と仮定）","#2E8B57","--"),("方法1：ON/OFF時刻を与える","#C0392B","-"),(best_rw,"#2E6FD8","-")]
    fig,axs=plt.subplots(2,3,figsize=(16,7.6))
    for j,case in enumerate(["learned","q25","intermittent"]):
        tr=traj[case]; t=tr["t"]; ax=axs[0,j]
        ax.plot(t,tr["truth"],"o-",color="k",ms=4,lw=2.2,label="真値")
        for n,c,ls in show: ax.plot(t[1:],tr[n][0][1:],ls,color=c,lw=2,label=n)
        ax.set_title(LAB[case],fontsize=11.5); ax.set_ylabel("反り A−O [µm]"); ax.grid(alpha=.3)
        ax=axs[1,j]; tt=np.linspace(0,600,1201); ax.plot(tt,[QTRUE[case](x) for x in tt],color="k",lw=2,label="真の発熱")
        for n,c,ls in show: ax.plot(t,tr[n][1],ls,color=c,lw=2,marker="o",ms=3,label=n)
        ax.set_ylabel("推定した発熱 [W]"); ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3); ax.set_ylim(-1,32)
    axs[0,0].legend(fontsize=8.5,loc="upper right"); axs[1,0].legend(fontsize=8.5,loc="upper right")
    fig.suptitle("改善の試み：上＝反り A−O、下＝推定した発熱（ヒータが OFF と仮定した時刻は 0 W として表示）。温度P2＋P4＋変位A/O、5 seed 平均",fontsize=12)
    fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig(os.path.join(IMG,"improve_timeseries.png"),dpi=150); plt.close(fig)
    fig,axs=plt.subplots(1,3,figsize=(16,4.8),sharey=True)
    names=[n for n,_,_ in METHODS]+["回帰式（温度P2+P4）"]; cols=["#2E8B57","#C0392B","#85C1E9","#2E6FD8","#1B4F72","#9AA5B1"]
    for ax,case in zip(axs,["learned","q25","intermittent"]):
        v=[out["cases"][case][n]["AO_um"] for n in names]; ax.bar(range(len(v)),v,color=cols)
        for i,x in enumerate(v): ax.text(i,x,f"{x:.2f}",ha="center",va="bottom",fontsize=9.5)
        ax.set_xticks(range(len(v))); ax.set_xticklabels(["基準","方法1\nON/OFF時刻","方法2\nσ=0.1","方法2\nσ=0.2","方法2\nσ=0.4","回帰式"],fontsize=9.5)
        ax.set_title(LAB[case],fontsize=11); ax.grid(axis="y",alpha=.3)
    axs[0].set_ylabel("反り A−O の誤差 [µm]（30〜600 s 平均）")
    fig.suptitle("改善の試み：反りの誤差（温度P2＋P4＋変位A/O、5 seed 平均。回帰式は温度P2＋P4のみ）",fontsize=12)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"improve_bars.png"),dpi=150); plt.close(fig)


if __name__=="__main__": main()
