"""低感度点 vs 高感度点に温度センサ → Q・h 推定精度（研究引き継ぎ Step 4-9）.

sensitivity_analysis.py で検証済みの直接感度 S_Q=∂T/∂Q を PODで全場へ復元し、
|∂T/∂Q| が最大のセル（高感度点）と小さいセル（低感度点）を観測点に選ぶ。
拡大状態 z=[T1..T5, Q, h] の固定ゲインOIで、その1点の温度観測から Q,h を推定し、
配置による推定誤差の差を出す（＝blog_002の変位版[低W/高W]の Q,h 版）。

出力: results/placement_qh.npz, docs/img/sensitivity_qh_placement.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/placement_qh_comparison.py
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
from sensitivity_analysis import integrate_sensitivity
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
TAIR=rg.T_AIR_K; RATED=rg.HEATER_RATED_W
NPT=5; IQ=5; IH=6; NAUG=7
DT=2.0; OBS_DT=30.0; T_END=600.0
SIG_T=0.30                          # 温度観測ノイズ [K]
SEEDS=[20260913,20260914,20260915,20260916,20260917]
N_B=200


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],5); h_true=float(d["h"]); heat=int(d["heat_node"])
    Q_true=RATED*1.0
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])       # (5,5): T5-mean5 -> モード係数

    # --- 全場 ∂T/∂Q をPOD復元し、高/低感度セルを選ぶ ---
    ts,T,SQ,Sh=integrate_sensitivity(np.full(5,TAIR),C,Kmat,h_true,Q_true,heat,0.0,T_END,DT)
    SQ_full=SQ@UPp.T@U.T                # (nt, N)  S_Q_full = U (UPp S_Q_rom)
    peak=np.abs(SQ_full).max(0)         # 各セルのピーク|∂T/∂Q|
    c_hi=int(np.argmax(peak))
    lo_thresh=np.percentile(peak,8)     # 低感度=下位8%くらいの実セル
    c_lo=int(np.argmin(np.abs(peak-lo_thresh)))
    print(f"[place] 高感度セル c={c_hi} |dT/dQ|peak={peak[c_hi]:.4f}K/W xyz={np.round(Cc[c_hi]*1000,0).astype(int)}")
    print(f"[place] 低感度セル c={c_lo} |dT/dQ|peak={peak[c_lo]:.4f}K/W xyz={np.round(Cc[c_lo]*1000,0).astype(int)}")

    def Hrow_cell(c):                   # セルc温度の T5 に対する線形係数（1x5）
        return U[c,:]@UPp
    def cell_temp(T5, c):               # T5(5,) -> セルc温度
        return mean[c]+(T5-mean[pod])@Hrow_cell(c)

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    # 真値トラジェクトリ
    Tt=np.full(NPT,TAIR); Ttr=[Tt.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(Tt,C,Kmat,h_true,Q_true/RATED,heat,a,b,DT); Tt=tr[-1]; Ttr.append(Tt.copy())
    Ttr=np.array(Ttr)

    N_ENS=120; R=np.array([[SIG_T**2]])
    def run(mode, seed):
        """mode: 'none'(DAなし) / 'lo' / 'hi'（観測セル）。拡大状態EnKF。"""
        c = None if mode=="none" else (c_lo if mode=="lo" else c_hi)
        rng=np.random.default_rng(seed)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=TAIR                                   # 初期は全メンバー室温→伝播でQ,hと相関
        Z[:,IQ]=rng.uniform(5,25,N_ENS)                  # Q事前（真値15を含む広め）
        Z[:,IH]=rng.uniform(0.005,0.040,N_ENS)           # h事前（真値0.0154を含む広め）
        Qh=[(Z[:,IQ].mean(),Z[:,IH].mean())]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ]/RATED,heat,tp,tb,DT); tp=tb
            if mode!="none":
                Yf=(mean[c]+(Z[:,:NPT]-mean[pod])@Hrow_cell(c))[:,None]      # (N,1) 予報観測
                y=np.array([cell_temp(Ttr[ci],c)+rng.normal(0,SIG_T)])       # (1,) 観測
                Z=enkf_update(Z,y,None,R,rng,inflation=1.03,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,40); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            Qh.append((Z[:,IQ].mean(),Z[:,IH].mean()))
        return np.array(Qh)

    res={}
    for m in ["none","lo","hi"]:
        rr=np.array([run(m,s) for s in SEEDS])      # (seed, ncyc+1, 2)
        Qf=rr[:,-1,0]; hf=rr[:,-1,1]
        res[m]=dict(Q_err=np.abs(Qf-Q_true).mean(), h_err=np.abs(hf-h_true).mean(),
                    Q_std=Qf.std(), traj=rr.mean(0))
    for m in ["none","lo","hi"]:
        print(f"[place] {m}: |Q-Q_true|={res[m]['Q_err']:.2f} W,  |h-h_true|={res[m]['h_err']:.4f} W/K")

    # --- 図：Q誤差・h誤差の棒グラフ（DAなし/低感度/高感度）---
    labels=["データ同化なし","低感度点に温度計","高感度点に温度計"]
    keys=["none","lo","hi"]; cols=["tab:gray","tab:orange","tab:blue"]
    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(12.5,5.2))
    ax0.bar(labels,[res[k]["Q_err"] for k in keys],color=cols)
    ax0.set_ylabel("|Q推定誤差| [W]（真値15W）"); ax0.set_title("発熱量Qの推定誤差",fontsize=13,weight="bold")
    ax0.tick_params(axis="x",labelsize=10)
    for i,k in enumerate(keys): ax0.text(i,res[k]["Q_err"],f"{res[k]['Q_err']:.1f}",ha="center",va="bottom",weight="bold")
    ax1.bar(labels,[res[k]["h_err"] for k in keys],color=cols)
    ax1.set_ylabel("|h推定誤差| [W/K]"); ax1.set_title("放熱係数hの推定誤差",fontsize=13,weight="bold")
    ax1.tick_params(axis="x",labelsize=10)
    for i,k in enumerate(keys): ax1.text(i,res[k]["h_err"],f"{res[k]['h_err']:.3f}",ha="center",va="bottom",weight="bold")
    fig.suptitle("温度センサを『低感度点 vs 高感度点』に置いたときの Q・h 推定誤差（5seed平均）\n"
                 f"高感度セル|∂T/∂Q|={peak[c_hi]:.3f} / 低感度セル={peak[c_lo]:.3f} K/W（直接感度から選定）",
                 fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.9])
    out=os.path.join(IMG,"sensitivity_qh_placement.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[place] wrote",out)
    np.savez(os.path.join(RES,"placement_qh.npz"),
             c_hi=c_hi,c_lo=c_lo,peak_hi=peak[c_hi],peak_lo=peak[c_lo],
             Q_true=Q_true,h_true=h_true,
             **{f"{m}_{k}":res[m][k] for m in res for k in ["Q_err","h_err","Q_std"]})


if __name__=="__main__": main()
