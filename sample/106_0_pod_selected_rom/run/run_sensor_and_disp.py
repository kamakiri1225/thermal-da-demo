"""① 熱感度(dT/dQ)による温度センサ選定の検証、② 変位推定が真値を捉えるかの検証.

① 5ノードそれぞれを「唯一の温度センサ」にしてDA(5seed)し、過渡期RMSEと dT/dQ の関係を見る
   → dT/dQ が高いノードほど推定が良い、を確認（熱感度でセンサを選ぶ根拠）。
② 温度2点+変位2点でDAし、変位QoI=Uz(A)-Uz(O)の推定が真値を捉えるかを時刻歴で見る。

出力(docs/img/): da_sensor_select.png, da_disp_track.png
再現: OMP_NUM_THREADS=4 python3 run/run_sensor_and_disp.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from run_da_compare import build_disp_operator
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
K=273.15; NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    xyz=d["xyz"]; heat_node=int(d["heat_node"]); C=d["C"]
    Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    pod_cells=kv["cell_idx"]; UP_pinv=np.linalg.pinv(U[pod_cells,:])
    uz_mean,Dmode=build_disp_operator(U,mean,Cc)
    def disp_of_T5(T5):
        a=(T5-mean[pod_cells])@UP_pinv.T; return uz_mean+a@Dmode.T
    def qoi_of_T5(T5):      # QoI = Uz(A)-Uz(O)
        u=disp_of_T5(T5); return u[...,0]-u[...,1]

    # dT/dQ
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens=(pert-base)/0.1

    # 真値トラジェクトリ
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tgrid=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(nodes,use_disp,seed):
        """nodes=None なら同化なし(free run)。戻り: rmse, qest, recT(各ノード推定[nt,5])."""
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        if nodes is not None:
            Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
        rmse=[np.sqrt(((Z[:,:NPT].mean(0)-Ttr[0])**2).mean())]
        qest=[qoi_of_T5(Z[:,:NPT].mean(0))]; recT=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Tn=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
            Z=Z.copy(); Z[:,:NPT]=Tn; tp=tb
            if nodes is not None:
                yv=list(Ttr[ci][nodes]); Yf=Z[:,nodes]
                if use_disp:
                    yv=yv+list(disp_of_T5(Ttr[ci])); Yf=np.column_stack([Yf,disp_of_T5(Z[:,:NPT])])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            rmse.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean()))
            qest.append(qoi_of_T5(Z[:,:NPT].mean(0))); recT.append(Z[:,:NPT].mean(0).copy())
        return np.array(rmse),np.array(qest),np.array(recT)

    # ① 各ノード単独センサの過渡期RMSE（5seed平均）
    heat=(tgrid>0)&(tgrid<=300)
    node_rmse=[]
    for i in range(NPT):
        rs=[run([i],False,s)[0] for s in SEEDS]
        node_rmse.append(np.mean([r[heat].mean() for r in rs]))
    node_rmse=np.array(node_rmse)
    print("[ss] dT/dQ =",np.round(sens,2))
    print("[ss] 各ノード単独センサの過渡期RMSE =",np.round(node_rmse,3))

    fig,ax=plt.subplots(figsize=(9,5.8))
    # 傾向線(最小二乗)
    cc=np.corrcoef(sens,node_rmse)[0,1]
    p=np.polyfit(sens,node_rmse,1); xs=np.linspace(sens.min(),sens.max(),20)
    ax.plot(xs,np.polyval(p,xs),"--",color="gray",lw=1.5,label=f"傾向線 (相関 r={cc:.2f})")
    for i in range(NPT):
        ax.scatter(sens[i],node_rmse[i],s=180,color=plt.cm.tab10(i),zorder=3)
        ax.annotate(f"P{i}",(sens[i],node_rmse[i]),textcoords="offset points",xytext=(8,4),fontsize=13)
    # 低感度域を薄く帯掛け
    ax.axvspan(sens.min()-0.2,3.5,color="tab:orange",alpha=0.08)
    ax.text(sens.min()+0.1,node_rmse.max()*0.97,"低感度域\n(明確に悪い→避ける)",fontsize=10,color="darkorange",va="top")
    ax.set_xlabel("温度センサ位置の熱感度 dT/dQ [K/(発熱倍率)]")
    ax.set_ylabel("加熱期(0-300s)の平均5点温度RMSE [K]")
    ax.set_title(f"熱感度と推定精度: 低感度は明確に悪い／高感度側は感度だけでは決まらない\n"
                 "（センサの価値 = 未知量への感度 × 場の代表性）")
    ax.grid(alpha=0.3); ax.legend(loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_sensor_select.png"),dpi=140); plt.close(fig)

    # --- per-node比較: 正解／同化なし／高感度センサ／低感度センサ ---
    hi=int(np.argmax(sens)); lo=int(np.argmin(sens))
    _,_,rT_free=run(None,False,SEEDS[0])
    _,_,rT_hi=run([hi],False,SEEDS[0])
    _,_,rT_lo=run([lo],False,SEEDS[0])
    fig,axes=plt.subplots(2,3,figsize=(15,8.2),sharex=True)
    for i in range(NPT):
        ax=axes[i//3][i%3]; ax.axvspan(0,300,color="orange",alpha=0.06)
        tag=" (高感度=観測点)" if i==hi else (" (低感度=観測点)" if i==lo else " (未観測)")
        ax.plot(tgrid,Ttr[:,i]-K,"-",color="k",lw=4,alpha=0.35)
        ax.plot(tgrid,rT_free[:,i]-K,"-",color="tab:gray",lw=1.6)
        ax.plot(tgrid,rT_hi[:,i]-K,"--",color="tab:blue",lw=1.8,dashes=(4,3))
        ax.plot(tgrid,rT_lo[:,i]-K,"--",color="tab:orange",lw=1.8,dashes=(2,2))
        ax.set_title(f"P{i}{tag}",fontsize=12); ax.grid(alpha=0.3)
        if i%3==0: ax.set_ylabel("温度 [degC]")
        if i//3==1: ax.set_xlabel("time [s]")
    style=[Line2D([0],[0],color="k",lw=4,alpha=0.35,label="正解(真値)"),
           Line2D([0],[0],color="tab:gray",lw=1.6,label="同化なし(free run)"),
           Line2D([0],[0],color="tab:blue",lw=1.8,ls="--",label="同化: 高感度センサで観測"),
           Line2D([0],[0],color="tab:orange",lw=1.8,ls="--",label="同化: 低感度センサで観測")]
    axes[1][2].axis("off"); axes[1][2].legend(handles=style,fontsize=12,loc="center")
    fig.suptitle("各ノードの温度: 正解 / 同化なし / 高感度センサ同化 / 低感度センサ同化",fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.96]); fig.savefig(os.path.join(IMG,"da_pernode_compare.png"),dpi=130); plt.close(fig)
    print(f"[ss] per-node比較図: 高感度=P{hi}, 低感度=P{lo}")

    # ② 変位QoIの推定 vs 真値（同化なし / 温度2点 / 温度2点+変位2点, 5seed平均）
    hi2=list(np.argsort(sens)[::-1][:2])
    qoi_true=qoi_of_T5(Ttr)
    qoi_free=np.mean([run(None,False,s)[1] for s in SEEDS],axis=0)   # 同化なし(free run)
    qoi_td  =np.mean([run(hi2,False,s)[1] for s in SEEDS],axis=0)    # 温度2点のみ
    qoi_tdu =np.mean([run(hi2,True ,s)[1] for s in SEEDS],axis=0)    # 温度2点+変位2点
    # 過渡期(0-300s)の変位QoI誤差[µm]で比較
    ht=(tgrid>0)&(tgrid<=300)
    e_free=np.abs(qoi_free-qoi_true)[ht].mean()
    e_td  =np.abs(qoi_td  -qoi_true)[ht].mean()
    e_tdu =np.abs(qoi_tdu -qoi_true)[ht].mean()
    print(f"[ss] 変位QoI 加熱期平均誤差[µm]: 同化なし={e_free:.2f} 温度2点={e_td:.2f} 温度2点+変位2点={e_tdu:.2f}")
    fig,ax=plt.subplots(figsize=(11,5.8)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(tgrid,qoi_true,"-",color="k",lw=4,alpha=0.35,label="真値")
    ax.plot(tgrid,qoi_free,":",color="tab:gray",lw=2.2,label=f"同化なし(free run)  誤差{e_free:.1f}µm")
    ax.plot(tgrid,qoi_td,"--",color="tab:green",lw=2,label=f"温度2点で同化→変位推定  誤差{e_td:.1f}µm")
    ax.plot(tgrid,qoi_tdu,"--",color="tab:red",lw=2.2,label=f"温度2点+変位2点で同化→変位推定  誤差{e_tdu:.1f}µm")
    ax.set_xlabel("time [s]"); ax.set_ylabel("変位QoI Uz(A)−Uz(O) [µm]")
    ax.set_title("変位の推定は真値を捉えているか（gappy-POD場→FrontISTR）\n"
                 "※誤差=加熱期(0-300s)の真値との平均絶対差")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_disp_track.png"),dpi=140); plt.close(fig)
    print("[ss] wrote da_sensor_select.png, da_disp_track.png")


if __name__=="__main__": main()
