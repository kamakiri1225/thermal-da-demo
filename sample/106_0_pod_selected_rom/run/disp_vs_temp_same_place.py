"""同じ場所に「変位計」を置くか「温度計」を置くかで、推定精度がどう変わるかを比べる.

本番は 温度2点(P2+P0) ＋ 変位2点(上面A/OのUz)。
では、A/O と同じ場所に温度計を置いたらどうなるか。センサの本数（4本）と場所をそろえ、
変えるのは「そこで何を測るか」だけにした対照実験。

構成:
  (a) 温度2点のみ（P2+P0）                      … 参考
  (b) 温度2点 ＋ A/O に温度計2本                … 同じ場所で温度を測る
  (c) 温度2点 ＋ A/O に変位計2本（本番）        … 同じ場所で変位を測る
  (d) 温度4点（P2+P0＋A/O温度）＋ 変位2点       … 参考（全部入り）

A/O の温度は ROM 代表点ではないので、POD 復元を観測演算子にする（run/obs_points_not_limited_to_rom.py と同じ）。
評価: 全20,696セルの温度RMSE、上面の反り A−O の誤差、個別 Uz(A)/Uz(O) の誤差。
真値も同じ ROM の双子実験。EnKF・60メンバー・30秒ごと20回・5 seed・加熱期30〜300秒。

出力: results/disp_vs_temp_same_place.json, docs/img/disp_vs_temp_same_place.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/disp_vs_temp_same_place.py
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
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TBASE=[2,0]                       # 本番の温度2点（P2, P0）
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005)


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:]); A=U@UPp
    op=np.load(os.path.join(RES,"disp_operator.npz")); ub=op["uz_mean"]; Dm=op["Dmode"]
    Wd=Dm@UPp                                       # (2,5) 温度→A/Oの変位 [µm/K]
    # A/O に最も近い OpenFOAM セル（ここに温度計を置く想定）
    tcell=[int(np.linalg.norm(Cc-np.array(x),axis=1).argmin()) for x in (A_XYZ,O_XYZ)]
    Wt=np.array([A[c] for c in tcell])              # (2,5) 5点温度→そのセルの温度
    t0=np.array([mean[c] for c in tcell])
    print("A/Oに最も近いセル:",[(c,np.round(Cc[c]*1000,1).tolist()) for c in tcell])
    field=lambda T5: mean+A@(T5-mean[pod])
    disp =lambda T5: ub+Wd@(T5-mean[pod])
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(nodes,use_tAO,use_dAO,seed):
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_T**2]*2 if use_tAO else [])+([SIG_U**2]*2 if use_dAO else []))
        eF=[];eAO=[];eA=[];eO=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][nodes]); Yf=Z[:,nodes]
            if use_tAO:
                yv+=list(t0+Wt@(Ttr[ci]-mean[pod])); Yf=np.column_stack([Yf,t0+(Z[:,:NPT]-mean[pod])@Wt.T])
            if use_dAO:
                yv+=list(disp(Ttr[ci])); Yf=np.column_stack([Yf,ub+(Z[:,:NPT]-mean[pod])@Wd.T])
            y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                m=Z[:,:NPT].mean(0); u=disp(m); ut=disp(Ttr[ci])
                eF.append(np.sqrt(((field(m)-field(Ttr[ci]))**2).mean()))
                eAO.append(abs((u[0]-u[1])-(ut[0]-ut[1]))); eA.append(abs(u[0]-ut[0])); eO.append(abs(u[1]-ut[1]))
        return [np.mean(eF),np.mean(eAO),np.mean(eA),np.mean(eO),15*Z[:,IQ].mean()]

    cfgs=[("温度2点のみ（P2+P0）",TBASE,False,False),
          ("＋A/Oに温度計2本（計4本とも温度）",TBASE,True,False),
          ("＋A/Oに変位計2本（本番）",TBASE,False,True),
          ("＋A/Oに温度計2本と変位計2本（計6本）",TBASE,True,True)]
    out={"note":"誤差は加熱期30〜300秒・5 seed平均。温度場は全20,696セル","obs_cells_for_AO_temp":[int(c) for c in tcell],
         "cells_xyz_mm":[np.round(Cc[c]*1000,1).tolist() for c in tcell],"configs":{}}
    for name,nd,ta,da in cfgs:
        r=np.array([run(nd,ta,da,s) for s in SEEDS])
        m=r.mean(0); sd=r.std(0,ddof=1)
        out["configs"][name]=dict(field_K=float(m[0]),field_K_sd=float(sd[0]),AO_um=float(m[1]),AO_um_sd=float(sd[1]),
                                  A_um=float(m[2]),O_um=float(m[3]),Q_W=float(m[4]),per_seed=r.tolist())
        print(f"{name:34s} 温度場 {m[0]:.3f}±{sd[0]:.3f} K  反り {m[1]:.3f}±{sd[1]:.3f} µm  "
              f"Uz(A) {m[2]:.3f}  Uz(O) {m[3]:.3f} µm  Q {m[4]:.2f} W",flush=True)
    # seed ごとの対応をとった勝敗（温度計 vs 変位計）
    t_=np.array(out["configs"]["＋A/Oに温度計2本（計4本とも温度）"]["per_seed"])
    u_=np.array(out["configs"]["＋A/Oに変位計2本（本番）"]["per_seed"])
    win={k:int((u_[:,i]<t_[:,i]).sum()) for i,k in enumerate(["温度場","反りA−O","Uz(A)","Uz(O)"])}
    out["disp_wins_over_temp_per_seed"]=win
    print("seedごとに変位計が温度計より良かった回数（5中）:",win)
    json.dump(out,open(os.path.join(RES,"disp_vs_temp_same_place.json"),"w"),ensure_ascii=False,indent=1)
    names=list(out["configs"]); cols=["#9AA5B1","#E67E22","#C0392B","#2E6FD8"]
    fig,axs=plt.subplots(1,3,figsize=(16,4.8))
    for ax,key,ttl,un in [(axs[0],"field_K","温度場全体（20,696セル）","K"),
                          (axs[1],"AO_um","反り A−O","µm"),(axs[2],"A_um","上面 A の変位","µm")]:
        v=[out["configs"][n][key] for n in names]
        sd=[out["configs"][n].get(key+"_sd",0) for n in names]
        ax.bar(range(len(v)),v,yerr=sd,color=cols,capsize=3)
        for i,x in enumerate(v): ax.text(i,x,f"{x:.3f}",ha="center",va="bottom",fontsize=10)
        ax.set_xticks(range(len(v))); ax.set_xticklabels(["温度2点\nのみ","＋温度計\n2本","＋変位計\n2本\n(本番)","＋両方\n(6本)"],fontsize=10)
        ax.set_title(ttl,fontsize=12); ax.set_ylabel(f"誤差 [{un}]"); ax.grid(axis="y",alpha=.3)
    fig.suptitle("同じ場所（上面 A・O）で「温度」を測るか「変位」を測るか ― 温度2点 P2+P0 は共通、5 seed 平均",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.91)); fig.savefig(os.path.join(IMG,"disp_vs_temp_same_place.png"),dpi=150); plt.close(fig)
    print("wrote figure")


if __name__=="__main__": main()
