"""なぜ「熱感度Wが高い温度点」を選ぶと悪くなるのかを調べる.

直感：温度が1K変わると変位が大きく動く点(W大)を測るべき。
実測：温度1点だけなら W低(P3) のほうが W高(P1) より良い。

5点それぞれを唯一の温度センサにして同化し、
  ・温度RMSE ・変位差の誤差
を測る。そのうえで、結果を説明しうる指標を並べて相関を見る:
  (a) 熱感度 W_i = |d(Uz(A)-Uz(O))/dT_i|      「そこの温度誤差が変位にどれだけ響くか」
  (b) 観測の価値 Δ(u, T_i)                     「そこを測ると変位の分散がどれだけ減るか」
  (c) その点と他4点の温度の相関の平均           「その点は場全体を代表しているか」
  (d) その点の温度の事前分散                     「そこはそもそもばらつくのか」

出力: docs/img/why_lowW_wins.png, results/why_lowW_wins.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/why_lowW_wins_check.py
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
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
N_PRIOR=4000


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    Cc=kv["cell_centres"]; UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    iAz,iOz=2,5
    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T

    W=np.abs((D[iAz]-D[iOz])@UPp)        # (a) 熱感度
    print(f"発熱ノード = P{heat}")
    print("点の座標 [mm]:", {f"P{i}":list(np.round(Cc[pod[i]]*1000,1)) for i in range(NPT)},"\n")

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr); qtrue=Utr[:,iAz]-Utr[:,iOz]

    # ── 事前アンサンブル（指標の計算用）──
    rng=np.random.default_rng(20260930)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    tp=0.0
    for tb in np.arange(OBS_DT,150.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,tp,tb,DT); tp=tb
    upr=disp(Tp); Xq=upr[:,iAz]-upr[:,iOz]
    R=np.corrcoef(Tp.T)
    delta=np.array([np.cov(Xq,Tp[:,i])[0,1]**2/(Tp[:,i].var()+SIG_T**2) for i in range(NPT)])
    reprs=np.array([(R[i].sum()-1)/(NPT-1) for i in range(NPT)])   # (c) 代表性
    varT=Tp.var(axis=0)                                            # (d) 事前分散
    # (e) 発熱への温度感度 dT/dQ（有限差分, 加熱ピーク t=300s）
    b0=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.0,heat,0,300,DT)[1][-1]
    b1=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.1,heat,0,300,DT)[1][-1]
    dTdQ=(b1-b0)/0.1

    # ── 1点ずつ同化 ──
    def run(node,seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.array([[SIG_T**2]]); eT=[]; eU=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            y=np.array([Ttr[ci][node]])+rng_o.normal(0,SIG_T)
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Z[:,[node]])
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                m=Z[:,:NPT].mean(0)
                eT.append(np.sqrt(((m-Ttr[ci])**2).mean()))
                uu=disp(m); eU.append(abs((uu[iAz]-uu[iOz])-qtrue[ci]))
        return float(np.mean(eT)), float(np.mean(eU))

    rmseT=np.zeros(NPT); maeU=np.zeros(NPT)
    for i in range(NPT):
        r=[run(i,s) for s in SEEDS]
        rmseT[i]=np.mean([x[0] for x in r]); maeU[i]=np.mean([x[1] for x in r])

    print(f"{'点':>4}{'W [µm/K]':>11}{'dT/dQ':>9}{'Δ(u,T)':>10}{'代表性':>9}"
          f"{'温度RMSE':>11}{'変位MAE':>10}")
    for i in range(NPT):
        print(f"P{i:>3}{W[i]:11.4f}{dTdQ[i]:9.2f}{delta[i]:10.4f}{reprs[i]:9.3f}"
              f"{rmseT[i]:11.4f}{maeU[i]:10.4f}")

    names=["(a) 熱感度 $W_i$\n（温度→変位）","(b) 観測の価値 $\\Delta(u,T_i)$",
           "(c) 他点との相関の平均","(d) 事前分散 $\\mathrm{Var}(T_i)$",
           "(e) 発熱への温度感度\n$\\partial T_i/\\partial Q$"]
    crit=[W,delta,reprs,varT,dTdQ]
    print("\n各指標と実測の相関（温度RMSE・変位MAEは小さいほど良いので、良い指標なら負）")
    corrs=[]
    for nm,cv in zip(names,crit):
        cT=float(np.corrcoef(cv,rmseT)[0,1]); cU=float(np.corrcoef(cv,maeU)[0,1])
        corrs.append((nm,cT,cU))
        print(f"  {nm:<34} vs 温度RMSE {cT:+.3f}   vs 変位MAE {cU:+.3f}")

    fig,axes=plt.subplots(1,5,figsize=(22.5,5.2))
    for k,(ax,(nm,cv)) in enumerate(zip(axes,zip(names,crit))):
        # 温度センサの話なので、変位MAEと温度RMSEの両方を見る
        ax.scatter(cv,maeU,s=200,c="#C0392B",zorder=4,label="変位MAE")
        ax2=ax.twinx(); ax2.scatter(cv,rmseT,s=160,marker="s",
                                    c="#1F4E9C",zorder=4,alpha=.75,label="温度RMSE")
        ax2.set_ylabel("温度RMSE [K]",fontsize=10.5,color="#1F4E9C")
        ax2.tick_params(axis="y",labelcolor="#1F4E9C",labelsize=9)
        cT=float(np.corrcoef(cv,rmseT)[0,1])
        for i in range(NPT):
            ax.annotate(f"P{i}",(cv[i],maeU[i]),textcoords="offset points",
                        xytext=(9,7),fontsize=12,weight="bold")
        c=float(np.corrcoef(cv,maeU)[0,1])
        ax.set_xlabel(nm,fontsize=11.5)
        ax.set_ylabel("変位差の誤差 [µm]",fontsize=10.5,color="#C0392B")
        ax.tick_params(axis="y",labelcolor="#C0392B",labelsize=9)
        ax.set_title(f"変位 {c:+.3f} ／ 温度 {cT:+.3f}",fontsize=12.5,weight="bold",
                     color=("#2e7d32" if min(c,cT)<-0.5 else "#C0392B" if max(c,cT)>0.5 else "#556"))
        ax.grid(alpha=.3)
    fig.suptitle("なぜ「熱感度 $W$ が高い温度点」を選ぶと悪くなるのか ― 4つの指標と実測の突き合わせ\n"
                 "温度1点だけで同化（60メンバー・5 seed平均・加熱期）。"
                 "赤●＝変位差の誤差、青■＝温度RMSE。良い指標なら右下がり（相関が負）",
                 fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.86])
    o=os.path.join(IMG,"why_lowW_wins.png"); fig.savefig(o,dpi=130); plt.close(fig)
    json.dump(dict(heat_node=heat,
                   points={f"P{i}":[round(v*1000,1) for v in Cc[pod[i]]] for i in range(NPT)},
                   W=W.tolist(),dTdQ=dTdQ.tolist(),delta=delta.tolist(),repr_corr=reprs.tolist(),
                   prior_var=varT.tolist(),rmseT=rmseT.tolist(),maeU=maeU.tolist(),
                   correlations=[dict(criterion=n,vs_rmseT=a,vs_maeU=b) for n,a,b in corrs]),
              open(os.path.join(RES,"why_lowW_wins.json"),"w"),ensure_ascii=False,indent=2)
    print("wrote",o)


if __name__=="__main__":
    main()
