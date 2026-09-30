"""観測の価値 Δ は、実際のデータ同化の効き目を予測できているかを検証する.

Δ が予測するのは「注目量 X の分散がどれだけ減るか」なので、比較相手も分散にする:
  ・予測 : Δ/Var_f(X)   事前(=予報)共分散から計算した減少率
  ・実測 : 1 - Var_a(X)/Var_f(X)   60メンバーEnKFの各サイクルで、
           更新前後のアンサンブル分散から実際に測った減少率
参考として、真値に対するMAEも出す。ただし注目量が「差」なので、
同化なしは A と O が同方向に外れて相殺し、MAEでは有利に見える点に注意
（blog_004 §7-4 に記載の既知の罠）。

出力: docs/img/delta_vs_da.png, results/delta_vs_da.json

出力: docs/img/delta_vs_da.png, results/delta_vs_da.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/verify_delta_predicts_da.py
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
N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
N_PRIOR=4000; T_EVAL=150.0


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    iAz,iOz,iBz=2,5,8

    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T

    # ── 事前アンサンブル（Δ計算用）──
    # 実際のEnKFの「1サイクル目の予報」と同じ土俵にするため、t=OBS_DT まで前進させる
    rng=np.random.default_rng(20260930)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,0.0,OBS_DT,DT)
    up=disp(Tp)

    # 推定したい量 X（本研究の注目量：上面の反り）
    X=up[:,iAz]-up[:,iOz]
    varX=X.var()

    # 観測候補
    cands=[]
    for i in range(NPT):
        cands.append((f"温度P{i}", ("T",i), Tp[:,i], SIG_T**2))
    for lab,k in [("変位 Uz(A)",iAz),("変位 Uz(O)",iOz),("変位 Uz(B)",iBz)]:
        cands.append((lab, ("U",k), up[:,k], SIG_U**2))

    # ── 真値トラジェクトリ ──
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr); Qoi_true=Utr[:,iAz]-Utr[:,iOz]

    # ── 観測1つだけで同化し、更新前後のアンサンブル分散を記録 ──
    def run(kind, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        red=[]; mae=[]; tp=0.0
        u0=disp(Z[:,:NPT].mean(0)); mae.append(abs((u0[iAz]-u0[iOz])-Qoi_true[0]))
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            def qoi(Zs):
                uu=disp(Zs[:,:NPT]); return uu[:,iAz]-uu[:,iOz]
            vf=qoi(Z).var()                                  # 予報分散
            if ci==1: v1f=vf                                  # 1サイクル目の予報分散
            if kind is not None:
                t,i=kind
                if t=="T":
                    Rd=np.array([[SIG_T**2]]); yv=[Ttr[ci][i]]; Yf=Z[:,[i]]
                else:
                    Rd=np.array([[SIG_U**2]]); yv=[Utr[ci][i]]; Yf=disp(Z[:,:NPT])[:,[i]]
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                va=qoi(Z).var()                              # 解析分散
                if ci==1 and vf>0: red.append(1-va/vf)       # 1サイクル目だけ（Δと同じ土俵）
            uu=disp(Z[:,:NPT].mean(0)); mae.append(abs((uu[iAz]-uu[iOz])-Qoi_true[ci]))
        rmseT=None
        return (np.mean(red)*100 if red else 0.0), float(np.mean(mae[1:11]))

    e_free=np.mean([run(None,s)[1] for s in SEEDS])
    rows=[]
    for lab,kind,y,r in cands:
        dlt=np.cov(X,y)[0,1]**2/(y.var()+r)
        pred=dlt/varX*100
        res=[run(kind,s) for s in SEEDS]
        act=float(np.mean([x[0] for x in res])); e=float(np.mean([x[1] for x in res]))
        rows.append(dict(label=lab,delta=float(dlt),pred_pct=float(pred),
                         mae_um=e,act_pct=act))
        print(f"{lab:<12} Δ予測 {pred:6.1f}%   実測(分散減少) {act:6.1f}%   参考MAE {e:.3f} µm")
    print(f"\n参考: 同化なしのMAE = {e_free:.3f} µm"
          f"（注目量が『差』なので相殺で小さく見える。blog_004 §7-4 の罠）")
    P=np.array([r["pred_pct"] for r in rows]); A=np.array([r["act_pct"] for r in rows])
    rho=float(np.corrcoef(P,A)[0,1])
    from scipy.stats import spearmanr
    sp=float(spearmanr(P,A).correlation)
    print(f"\nピアソン相関 {rho:+.3f} / スピアマン順位相関 {sp:+.3f}")

    fig,(a0,a1)=plt.subplots(1,2,figsize=(15.0,6.2))
    a0.scatter(P,A,s=150,c=["#1F4E9C"]*NPT+["#C0392B"]*3,zorder=4)
    for r,x,y in zip(rows,P,A):
        a0.annotate(r["label"],(x,y),textcoords="offset points",xytext=(8,7),fontsize=11)
    lim=[min(P.min(),A.min())-6,max(P.max(),A.max())+8]
    a0.plot(lim,lim,"--",color="0.5",lw=2,label="完全一致の線")
    a0.set_xlim(lim); a0.set_ylim(lim)
    a0.set_xlabel("$\\Delta$ から予測した分散減少 [%]",fontsize=12.5)
    a0.set_ylabel("実測：EnKFでの分散減少 $1-V_a/V_f$ [%]",fontsize=12.5)
    a0.set_title(f"① 予測 vs 実測\nピアソン {rho:+.3f} ／ 順位相関 {sp:+.3f}",
                 fontsize=13.5,weight="bold")
    a0.legend(fontsize=11); a0.grid(alpha=.3)

    o=np.argsort(-P); x=np.arange(len(rows))
    a1.bar(x-0.2,P[o],width=.4,color="#1F4E9C",label="$\\Delta$ の予測")
    a1.bar(x+0.2,A[o],width=.4,color="#C0392B",label="実測")
    a1.set_xticks(x); a1.set_xticklabels([rows[i]["label"] for i in o],fontsize=10,rotation=20)
    a1.set_ylabel("注目量の分散減少 [%]",fontsize=12.5)
    a1.set_title("② 観測候補ごとの比較（$\\Delta$ の大きい順）",fontsize=13.5,weight="bold")
    a1.legend(fontsize=11); a1.grid(alpha=.3,axis="y")

    fig.suptitle("観測の価値 $\\Delta$ は、実際のデータ同化の効き目を予測できているか\n"
                 "推定したい量＝上面の反り $U_z(A)-U_z(O)$、観測は1つだけ、60メンバー・5 seed平均・加熱期",
                 fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.89])
    out=os.path.join(IMG,"delta_vs_da.png"); fig.savefig(out,dpi=130); plt.close(fig)
    json.dump(dict(target="Uz(A)-Uz(O) 加熱期MAE",free_mae_um=float(e_free),
                   pearson=rho,spearman=sp,rows=rows),
              open(os.path.join(RES,"delta_vs_da.json"),"w"),ensure_ascii=False,indent=2)
    print("wrote",out)


if __name__=="__main__":
    main()
