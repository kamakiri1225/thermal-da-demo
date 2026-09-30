"""観測期間を延ばすと放熱 h は同定できるようになるのか.

blog_003 §4-0-b で「一様冷却の時定数は約15,800 s。600 s の観測はその 4% しか
見ていないので h の情報が蓄積されない。観測期間を時定数の1/3まで延ばせば
同定できるはず」と書いたが、これは推測だった。実際に回して確かめる。

観測窓を 600 s から 60,000 s まで延ばし、h の推定値とアンサンブルの
ばらつきがどう動くかを追う。加熱は 0-300 s のみ（以降は冷却）。

出力: docs/img/h_long_horizon.png, results/h_long_horizon.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/h_long_horizon_check.py
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
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0
N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915]
HORIZONS=[600,1800,5000,18000,60000]      # s


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    iAz,iOz=2,5
    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T
    wAO=(D[iAz]-D[iOz])@UPp; hi2=list(np.argsort(np.abs(wAO))[::-1][:2])

    # 一様冷却の時定数 τ = ΣC / (N*h)   （全ノードが同じ温度なら結合Kは効かない）
    tau=float(C.sum()/(NPT*h_true))
    print(f"ROM: ΣC={C.sum():.1f} J/K, h_true={h_true:.6f} W/K → 一様冷却の時定数 τ={tau:.0f} s "
          f"（{tau/3600:.1f} 時間）\n")

    # 信号（室温からの上昇）が観測ノイズを下回る時刻＝観測が情報を持たなくなる境界
    TMAX=max(HORIZONS)
    cyc_all=np.arange(OBS_DT,TMAX+1e-9,OBS_DT)
    # 真値トラジェクトリ（最長まで一度だけ）
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc_all[:-1]],cyc_all):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr)

    exc=(Ttr-rg.T_AIR_K).mean(1)
    tg_all=np.r_[0,cyc_all]
    after_peak=tg_all>300.0          # 加熱中は上昇0から始まるので、冷却期だけで探す
    below=np.where(after_peak & (exc<SIG_T))[0]
    t_noise=float(tg_all[below[0]]) if len(below) else float("inf")
    print(f"信号（室温からの上昇）が観測ノイズ {SIG_T} K を下回る時刻 = "
          f"{t_noise:.0f} s（{t_noise/tau:.2f}τ）以降は観測に情報がない\n")

    def run(T_END, seed):
        cyc=cyc_all[cyc_all<=T_END+1e-9]
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2)
        th=[0.0]; hm=[Z[:,IH].mean()]; hs=[Z[:,IH].std()]; qm=[Z[:,IQ].mean()]
        tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][hi2])+[Utr[ci][iAz],Utr[ci][iOz]]
            Yf=np.column_stack([Z[:,hi2],disp(Z[:,:NPT])[:,[iAz,iOz]]])
            y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            th.append(tb); hm.append(Z[:,IH].mean()); hs.append(Z[:,IH].std()); qm.append(Z[:,IQ].mean())
        return np.array(th),np.array(hm),np.array(hs),np.array(qm)

    res={}; rows=[]
    print(f"{'観測窓':>10}{'τ比':>8}{'h推定':>11}{'誤差':>9}{'ばらつき':>11}{'Q推定':>9}")
    for TE in HORIZONS:
        rs=[run(TE,s) for s in SEEDS]
        th=rs[0][0]
        hm=np.mean([r[1] for r in rs],axis=0); hs=np.mean([r[2] for r in rs],axis=0)
        qm=np.mean([r[3] for r in rs],axis=0)
        res[TE]=dict(t=th.tolist(),h_mean=hm.tolist(),h_std=hs.tolist())
        err=abs(hm[-1]-h_true)/h_true*100
        rows.append(dict(T_END=TE,tau_ratio=TE/tau,h_final=float(hm[-1]),
                         h_err_pct=float(err),h_std_final=float(hs[-1]),
                         h_std_init=float(hs[0]),q_final=float(qm[-1]*15)))
        print(f"{TE:8.0f}s{TE/tau:8.2f}{hm[-1]:11.5f}{err:8.0f}%{hs[-1]:11.5f}{qm[-1]*15:8.2f}W")
    print(f"\n真値 h={h_true:.5f} W/K,  Q=15.00 W")

    fig,(a0,a1,a2)=plt.subplots(1,3,figsize=(18.6,5.8))
    cols=plt.cm.viridis(np.linspace(0,.85,len(HORIZONS)))
    for (TE,c) in zip(HORIZONS,cols):
        t=np.array(res[TE]["t"])/3600; hm=np.array(res[TE]["h_mean"]); hs=np.array(res[TE]["h_std"])
        a0.plot(t,hm,lw=2.3,color=c,label=f"{TE/3600:.1f} h 観測")
        a0.fill_between(t,hm-hs,hm+hs,color=c,alpha=.10,linewidth=0)
    a0.axhline(h_true,color="k",ls="--",lw=2.4,label=f"真値 {h_true:.5f}")
    a0.axvline(tau/3600,color="#C0392B",ls=":",lw=2.2)
    a0.text(tau/3600,a0.get_ylim()[1]*0.96,f" 時定数 τ={tau/3600:.1f} h",
            color="#C0392B",fontsize=11,weight="bold",va="top")
    a0.set_xscale("log")
    a0.set_xlabel("時刻 [時間]（対数軸）",fontsize=12.5); a0.set_ylabel("放熱 $h$ [W/K]",fontsize=12.5)
    a0.set_title("① 観測窓を延ばすと $h$ はどう動くか\n帯＝±1標準偏差",fontsize=13,weight="bold")
    a0.legend(fontsize=9.5); a0.grid(alpha=.3)

    x=np.array([r["tau_ratio"] for r in rows])
    a1.plot(x,[r["h_err_pct"] for r in rows],"-o",lw=2.6,ms=9,color="#C0392B")
    for r in rows:
        a1.annotate(f"{r['T_END']/3600:.1f}h",(r["tau_ratio"],r["h_err_pct"]),
                    textcoords="offset points",xytext=(0,11),ha="center",fontsize=10,weight="bold")
    a1.axvline(1/3,color="#1F4E9C",ls="--",lw=2.2)
    a1.text(1/3,a1.get_ylim()[1]*0.55," 時定数の1/3\n（事前の推測）",color="#1F4E9C",
            fontsize=10.5,weight="bold",ha="center")
    a1.axvspan(t_noise/tau,a1.get_xlim()[1],color="#C0392B",alpha=.10)
    a1.text(t_noise/tau*1.15,a1.get_ylim()[1]*0.5,
            f" 信号 < ノイズ\n（{t_noise/tau:.1f}τ 以降）\n ここから発散",
            color="#C0392B",fontsize=10.5,weight="bold")
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlabel("観測窓 / 時定数 $\\tau$（対数軸）",fontsize=12.5)
    a1.set_ylabel("$h$ の推定誤差 [%]",fontsize=12.5)
    a1.set_title("② どこまで観測すれば当たるのか",fontsize=13,weight="bold")
    a1.grid(alpha=.3,which="both")

    a2.plot(x,[r["h_std_final"] for r in rows],"-o",lw=2.6,ms=9,color="#1F4E9C",label="最終ばらつき")
    a2.axhline(rows[0]["h_std_init"],color="0.5",ls="--",lw=2.2,
               label=f"初期ばらつき {rows[0]['h_std_init']:.4f}")
    a2.set_xscale("log"); a2.set_yscale("log")
    a2.set_xlabel("観測窓 / 時定数 $\\tau$（対数軸）",fontsize=12.5)
    a2.set_ylabel("$h$ アンサンブルの標準偏差 [W/K]",fontsize=12.5)
    a2.set_title("③ 情報は入っているか\n初期より下がれば観測から学べている",fontsize=13,weight="bold")
    a2.axvspan(t_noise/tau,a2.get_xlim()[1],color="#C0392B",alpha=.10)
    a2.legend(fontsize=10.5); a2.grid(alpha=.3,which="both")

    fig.suptitle("観測期間を延ばせば放熱 $h$ は同定できるようになるのか "
                 f"― ROMの一様冷却時定数 $\\tau$={tau:.0f} s（{tau/3600:.1f} 時間）\n"
                 "観測＝温度2点＋変位2点、加熱は 0–300 s のみ、60メンバー・3 seed平均",
                 fontsize=14.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.875])
    o=os.path.join(IMG,"h_long_horizon.png"); fig.savefig(o,dpi=130); plt.close(fig)
    json.dump(dict(h_true=h_true,tau_s=tau,t_signal_below_noise_s=t_noise,
                   sweet_spot="0.3tau - 1.1tau", rows=rows),
              open(os.path.join(RES,"h_long_horizon.json"),"w"),ensure_ascii=False,indent=2)
    print("wrote",o)


if __name__=="__main__":
    main()
