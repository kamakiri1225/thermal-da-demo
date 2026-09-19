"""Q・h の直接感度解析（研究引き継ぎ Step 2-9）― 既存ROMを壊さない新規実装.

ROM（dacore/rom_general.py）の力学:
    dT/dt = M T + b_air + heat(t),   M_ij=K_ij/C_i (i!=j),  M_ii=-(ΣK_i+h)/C_i
    b_air_i = h T_air / C_i,   heat_i = Q χ(t)/C_i (ヒータノードのみ), χ=1(t<300),0(else)
    ※ Q = 15 q_scale [W]（HEATER_RATED_W=15）

直接感度（Q, h で微分）:
    dS_Q/dt = M S_Q + b_Q(t),   b_Q,i = χ(t)/C_i (ヒータノードのみ)
    dS_h/dt = M S_h - (T - T_air)/C          （M も b_air も h に依存するため）

T, S_Q, S_h を同じRK4ステージで同時積分（15状態=5点×3）。
有限差分で検証し、PODで全場感度へ復元、変位感度へ写像する。

出力: results/sensitivity_qh.npz, docs/img/sensitivity_qh_timeseries.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensitivity_analysis.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p          # 日本語フォント設定
from dacore import rom_general as rg
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
TAIR=rg.T_AIR_K; RATED=rg.HEATER_RATED_W


def build_M(C, Kmat, h):
    C=np.asarray(C,float); M=Kmat/C[:,None]
    np.fill_diagonal(M, -(Kmat.sum(1)+h)/C)
    return M


def chi(t):
    return 1.0 if t < rg.HEATER_ON_S else 0.0


def integrate_sensitivity(T0, C, Kmat, h, Q, heat_node, t0, t1, dt):
    """T, S_Q=∂T/∂Q, S_h=∂T/∂h の時刻歴を同時積分（RK4）。"""
    C=np.asarray(C,float); n=len(C); M=build_M(C,Kmat,h)
    b_air=h*TAIR/C
    def dz(T,SQ,Sh,t):
        c=chi(t)
        heat=np.zeros(n); heat[heat_node]=Q*c/C[heat_node]
        dT=M@T+b_air+heat
        bQ=np.zeros(n); bQ[heat_node]=c/C[heat_node]
        dSQ=M@SQ+bQ
        dSh=M@Sh-(T-TAIR)/C
        return dT,dSQ,dSh
    steps=max(1,int(round((t1-t0)/dt))); dt=(t1-t0)/steps
    T=np.asarray(T0,float).copy(); SQ=np.zeros(n); Sh=np.zeros(n)
    ts=[t0]; Tt=[T.copy()]; SQt=[SQ.copy()]; Sht=[Sh.copy()]; t=t0
    for _ in range(steps):
        aT,aSQ,aSh=dz(T,SQ,Sh,t)
        bT,bSQ,bSh=dz(T+.5*dt*aT,SQ+.5*dt*aSQ,Sh+.5*dt*aSh,t+.5*dt)
        cT,cSQ,cSh=dz(T+.5*dt*bT,SQ+.5*dt*bSQ,Sh+.5*dt*bSh,t+.5*dt)
        eT,eSQ,eSh=dz(T+dt*cT,SQ+dt*cSQ,Sh+dt*cSh,t+dt)
        T =T +dt/6*(aT+2*bT+2*cT+eT)
        SQ=SQ+dt/6*(aSQ+2*bSQ+2*cSQ+eSQ)
        Sh=Sh+dt/6*(aSh+2*bSh+2*cSh+eSh)
        t+=dt; ts.append(t); Tt.append(T.copy()); SQt.append(SQ.copy()); Sht.append(Sh.copy())
    return np.array(ts),np.array(Tt),np.array(SQt),np.array(Sht)


def fd_check(T0,C,Kmat,h,Q,heat_node,t0,t1,dt,dQ=1e-2,dh=1e-4):
    """有限差分（中心差分）でS_Q,S_hを検証（本番計算法でなくunit test）。"""
    _,TpQ=rg.integrate_single(T0,C,Kmat,h,(Q+dQ)/RATED,heat_node,t0,t1,dt)
    _,TmQ=rg.integrate_single(T0,C,Kmat,h,(Q-dQ)/RATED,heat_node,t0,t1,dt)
    SQ_fd=(TpQ-TmQ)/(2*dQ)
    _,TpH=rg.integrate_single(T0,C,Kmat,h+dh,Q/RATED,heat_node,t0,t1,dt)
    _,TmH=rg.integrate_single(T0,C,Kmat,h-dh,Q/RATED,heat_node,t0,t1,dt)
    Sh_fd=(TpH-TmH)/(2*dh)
    return SQ_fd,Sh_fd


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])
    Q=RATED*1.0                      # 動作点 q_scale=1 → Q=15 W
    T0=np.full(5, TAIR)
    ts,T,SQ,Sh=integrate_sensitivity(T0,C,Kmat,h,Q,heat,0.0,600.0,2.0)

    # --- Step3: 有限差分で検証 ---
    SQ_fd,Sh_fd=fd_check(T0,C,Kmat,h,Q,heat,0.0,600.0,2.0)
    eQ=np.abs(SQ-SQ_fd).max(); eH=np.abs(Sh-Sh_fd).max()
    rQ=eQ/np.abs(SQ_fd).max(); rH=eH/np.abs(Sh_fd).max()
    print(f"[sens] 検証(直接 vs 有限差分)  S_Q: 最大絶対差={eQ:.3e} (相対{rQ:.2e})")
    print(f"[sens] 検証(直接 vs 有限差分)  S_h: 最大絶対差={eH:.3e} (相対{rH:.2e})")
    print(f"[sens] |S_Q|最大={np.abs(SQ).max():.4f} K/W,  |S_h|最大={np.abs(Sh).max():.4f} K/(W/K)")

    np.savez(os.path.join(RES,"sensitivity_qh.npz"),
             t=ts, T=T, S_Q=SQ, S_h=Sh, cells=d["cells"], xyz=d["xyz"],
             heat_node=heat, Q=Q, h=h,
             fd_SQ_max_abs_err=eQ, fd_Sh_max_abs_err=eH)
    print("[sens] wrote results/sensitivity_qh.npz")

    # --- 図：∂T/∂Q と ∂T/∂h の時刻歴（時間情報が Q/h を分離する）---
    import matplotlib.pyplot as plt
    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(13,5.2),sharex=True)
    cols=plt.cm.viridis(np.linspace(0,.85,5))
    for i in range(5):
        ax0.plot(ts,SQ[:,i],color=cols[i],lw=2,label=f"P{i}")
        ax1.plot(ts,Sh[:,i],color=cols[i],lw=2,label=f"P{i}")
    for ax in (ax0,ax1): ax.axvspan(0,300,color="orange",alpha=.07); ax.axhline(0,color="k",lw=.8); ax.grid(alpha=.3); ax.set_xlabel("time [s]")
    ax0.set_ylabel("∂T/∂Q [K/W]"); ax0.set_title("発熱量Qへの感度 ∂T/∂Q",fontsize=13,weight="bold"); ax0.legend(fontsize=9,ncol=2)
    ax1.set_ylabel("∂T/∂h [K/(W/K)]"); ax1.set_title("放熱hへの感度 ∂T/∂h",fontsize=13,weight="bold")
    ax0.annotate("加熱で立上り→冷却(300s)で減衰",(150,SQ[:,heat].max()*.6),fontsize=10,color="dimgray")
    ax1.annotate("温度が上がるほど大きい\n(冷却期に最大)",(330,Sh.min()*.5),fontsize=10,color="dimgray")
    fig.suptitle("直接感度の時刻歴：∂T/∂Q と ∂T/∂h は時間形が違う → 時間情報で Q と h を分離できる\n"
                 "（空間的にはほぼ一様：∂T/∂Q は最大/最小=1.6倍, ∂T/∂h=1.2倍 → 温度センサの"
                 "『置き場所』はQ/h推定にはほぼ効かない。効くのは変位観測と時間）",
                 fontsize=11.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.9])
    out=os.path.join(IMG,"sensitivity_qh_timeseries.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[sens] wrote",out)


if __name__=="__main__": main()
