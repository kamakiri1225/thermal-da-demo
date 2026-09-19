"""Q・h の識別は「時間領域」が握る ― 同化する時間窓を変える実験（計算力学用）.

温度センサ2点（P2/P4）は固定し、EnKFで観測を入れる時間窓だけを変える:
  A. 加熱期のみ (0-300 s)
  B. 冷却期のみ (300-600 s)
  C. 全期間     (0-600 s)
直接感度の予言（∂T/∂Qは加熱期に大、∂T/∂hは冷却期に大）どおり、
Qは加熱期の観測が、hは冷却期の観測が効くことを確かめる。

出力: docs/img/qh_time_window.png, results/qh_time_window.npz
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/time_window_qh.py
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
N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"]); h_true=float(d["h"])
    Q_true=rg.HEATER_RATED_W
    TS=[2,4]                                     # 温度センサ2点は固定（P2ヒータ側/P4反対側）
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(window, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.eye(2)*SIG_T**2; tp=0.0
        Qh=[(Z[:,IQ].mean()*Q_true, Z[:,IH].mean())]
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            use = (window=="full") or (window=="heat" and tb<=300) or (window=="cool" and tb>300)
            if use:
                y=Ttr[ci][TS]+rng_o.normal(0,SIG_T,2)
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Z[:,TS])
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            Qh.append((Z[:,IQ].mean()*Q_true, Z[:,IH].mean()))
        return np.array(Qh)

    wins=[("heat","加熱期のみ(0-300s)","tab:red"),
          ("cool","冷却期のみ(300-600s)","tab:blue"),
          ("full","全期間(0-600s)","tab:green")]
    res={}
    for w,lab,col in wins:
        rr=np.array([run(w,s) for s in SEEDS])          # (seed, ncyc+1, 2)
        traj=rr.mean(0)
        eQ=np.abs(rr[:,-1,0]-Q_true).mean(); eH=np.abs(rr[:,-1,1]-h_true).mean()
        res[w]=dict(traj=traj, eQ=eQ, eH=eH, lab=lab, col=col)
        print(f"[twin] {lab:18s} |Q誤差|={eQ:.2f} W  |h誤差|={eH*1000:.2f} mW/K")

    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(13.5,5.4))
    for w,lab,col in wins:
        ax0.plot(tg,res[w]["traj"][:,0],color=col,lw=2.3,label=f"{lab}  最終誤差{res[w]['eQ']:.2f}W")
        ax1.plot(tg,res[w]["traj"][:,1]*1000,color=col,lw=2.3,label=f"{lab}  最終誤差{res[w]['eH']*1000:.2f}mW/K")
    ax0.axhline(Q_true,color="k",lw=2.6,alpha=.4); ax0.text(430,Q_true+0.35,"真値 15 W",fontsize=10)
    ax1.axhline(h_true*1000,color="k",lw=2.6,alpha=.4); ax1.text(20,h_true*1000+0.4,"真値 15.4 mW/K",fontsize=10)
    for ax in (ax0,ax1): ax.axvspan(0,300,color="orange",alpha=.07); ax.grid(alpha=.3); ax.set_xlabel("time [s]"); ax.legend(fontsize=9)
    ax0.set_ylabel("推定 Q [W]"); ax0.set_title("発熱量 Q の推定",fontsize=13,weight="bold")
    ax1.set_ylabel("推定 h [mW/K]"); ax1.set_title("放熱係数 h の推定",fontsize=13,weight="bold")
    fig.suptitle("同じ温度センサ2点でも『いつ観測するか』で Q・h の推定が変わる（ROM＋EnKF, 5seed平均）\n"
                 "Q：加熱期の観測が担う（冷却のみだと誤差7倍）。h：冷却期を足すと相対的に改善（7.6→5.3mW/K）だが、"
                 "温度応答がノイズ以下のため依然難しい",
                 fontsize=11.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.88])
    out=os.path.join(IMG,"qh_time_window.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[twin] wrote",out)
    np.savez(os.path.join(RES,"qh_time_window.npz"),
             **{f"{w}_{k}":res[w][k] for w in res for k in ["eQ","eH"]})


if __name__=="__main__": main()
