"""【配置設計の比較】既存の最適センサ配置法と、本研究の目的指向選定を同じ土俵で比べる.

blog_008 §5 は「選定2点」と「悪い例」しか比べていなかった。
計算力学・計算工学の聴衆が気にするのは「既存の配置法より良いのか」なので、
同じ候補集合（全5,040節点×3成分のうち A/O から 12 mm 以上離れた 14,754 通り）に対して
次の6つの選び方を適用し、同じ holdout 試験（A・O は一度も観測しない）で比べる。

  1. ランダム                 … 20 通りの無作為な2点（下限の目安）
  2. 熱感度の大きさ最大          … 「よく動く点に置く」という素朴な基準
  3. Q-DEIM（QR列ピボット）   … Manohar et al. のモードベース選定（DEIM/Q-DEIM）
  4. A最適（trace 最小）      … 5点の温度を一様に良く知る古典的な最適計画
  5. D最適（det 最小）        … 同上（情報行列の行列式）
  6. 目的指向（本研究）       … σ_S² = wᵀ(B − BHᵀ(HBHᵀ+R)⁻¹HB)w の最小化＋再採点

真値は OpenFOAM＋FrontISTR（results/truth_disp_all_<case>.npz のキャッシュを使う。
holdout_AO_realsolver.py が作る。FrontISTR の再実行は不要）。3条件 × 5 seed。

出力: results/placement_method_comparison.json, docs/img/placement_method_comparison.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/placement_method_comparison.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
from scipy.linalg import qr
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from improve_heating_schedule import integrate as integ_on
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TN=[2,0]; COMP=["Ux","Uy","Uz"]; N_PRIOR=4000; PRIOR_SEED=20260930
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005); EXCL_R=0.012
N_RANDOM=20; RANDOM_SEED=20261003
SCHED={"learned":lambda t:(t<300.0),"q25":lambda t:(t<300.0),
       "intermittent":lambda t:(t<150.0)|((t>=300.0)&(t<450.0))}
CASES=[("learned","15 W"),("q25","25 W"),("intermittent","間欠加熱")]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); mp=mean[pod]
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; Dall=op["D"]; coords=op["coords"]
    Wall=np.einsum("nck,kj->ncj",Dall,UPp)          # (節点, 成分, 5点温度) µm/K
    near=lambda x:int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    iA=near(A_XYZ); iO=near(O_XYZ)
    wAO=Wall[iA,2]-Wall[iO,2]
    far=(np.linalg.norm(coords-np.array(A_XYZ),axis=1)>EXCL_R)&(np.linalg.norm(coords-np.array(O_XYZ),axis=1)>EXCL_R)
    Wf=Wall.reshape(-1,NPT); cand=np.where(np.repeat(far,3))[0]
    print(f"候補 {cand.size} 通り（全 {Wf.shape[0]} から A・O 近傍を除外）  w(A−O)={np.round(wAO,3)}")

    # ---- 事前共分散 B（時刻ごと。測定値は使わない）----
    rng=np.random.default_rng(PRIOR_SEED)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Bs=[]; tp=0.0
    for tb in np.arange(OBS_DT,300.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,tp,tb,DT); tp=tb; Bs.append(np.cov(Tp.T))
    Ht=np.eye(NPT)[TN]

    def post(B,ws):
        H=np.vstack([Ht]+[w[None,:] for w in ws])
        R=np.diag([SIG_T**2]*len(TN)+[SIG_U**2]*len(ws))
        K=B@H.T@np.linalg.inv(H@B@H.T+R); return B-K@H@B

    def scores(ws):
        """(σ_AO [µm], trace [K²], logdet) を時刻平均で返す"""
        sg=[];tr=[];ld=[]
        for B in Bs:
            P=post(B,ws); sg.append(np.sqrt(max(wAO@P@wAO,0))); tr.append(np.trace(P))
            ld.append(np.linalg.slogdet(P)[1])
        return float(np.mean(sg)),float(np.mean(tr)),float(np.mean(ld))

    def sweep(base,key):
        """候補を全部採点して key（0:σ_AO, 1:trace, 2:logdet）が最小のものを返す"""
        best=None
        for k in cand:
            if k in base: continue
            v=scores([Wf[j] for j in base]+[Wf[k]])[key]
            if best is None or v<best[1]: best=(int(k),v)
        return best

    def greedy(key,label):
        picks=[]
        for st in range(2):
            k,v=sweep(picks,key); picks.append(k)
            i,c=divmod(k,3)
            print(f"  [{label}] {st+1}点目 {np.round(coords[i]*1000,1)} {COMP[c]}  指標={v:.5g}",flush=True)
        return picks

    # ---- 各手法の選定 ----
    methods={}
    # 2. 熱感度の大きさ最大（素朴）
    nw=np.linalg.norm(Wf,axis=1); order=cand[np.argsort(-nw[cand])]
    methods["熱感度の大きさ最大（素朴）"]=[int(order[0]),int(order[1])]
    # 3. Q-DEIM（候補行列の列ピボット付き QR）
    _,_,piv=qr(Wf[cand].T,pivoting=True)
    methods["Q-DEIM（QR列ピボット）"]=[int(cand[piv[0]]),int(cand[piv[1]])]
    # 4,5,6 貪欲
    methods["A最適（温度の trace 最小）"]=greedy(1,"A最適")
    methods["D最適（温度の det 最小）"]=greedy(2,"D最適")
    methods["目的指向（本研究）"]=greedy(0,"目的指向")
    # 1. ランダム（複数通り）
    rr=np.random.default_rng(RANDOM_SEED)
    randoms=[[int(x) for x in rr.choice(cand,2,replace=False)] for _ in range(N_RANDOM)]

    for nm,ks in methods.items():
        s=scores([Wf[k] for k in ks])
        print(f"{nm:26s} 予測σ(A−O)={s[0]:.4f} µm  trace={s[1]:.4f}  logdet={s[2]:.3f}  "
              f"{[f'{np.round(coords[k//3]*1000,1).tolist()} {COMP[k%3]}' for k in ks]}",flush=True)

    # ---- holdout 試験（真値＝OpenFOAM＋FrontISTR）----
    def run(case,obs,seed,t,T5,utrue,aot,on):
        ws=[Wall[i,c] for i,c in obs]; u0=[um[i,c] for i,c in obs]
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*len(TN)+[SIG_U**2]*len(ws))
        wA=Wall[iA,2]; wO=Wall[iO,2]; aA=um[iA,2]; aO=um[iO,2]
        eAO=[]; tp=0.0
        for ci,tb in enumerate(t[1:],1):
            Z=Z.copy(); Z[:,:NPT]=integ_on(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,on,tp,tb); tp=tb
            yv=list(T5[ci][TN]); Yf=Z[:,TN]
            for (i,c),w,b0 in zip(obs,ws,u0):
                yv.append(utrue[ci,i,c]); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
            y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                dm=Z[:,:NPT].mean(0)-mp
                eAO.append(abs(((aA+wA@dm)-(aO+wO@dm))-aot[ci]))
        return float(np.mean(eAO))

    out={"note":"同じ候補集合（14,754通り）に6つの配置法を適用し、A・Oを一度も観測しない holdout で比べた。"
                "真値＝OpenFOAM＋FrontISTR。誤差は加熱期30〜300秒・5 seed平均 [µm]",
         "n_candidates":int(cand.size),"methods":{},"random":{}}
    for nm,ks in methods.items():
        out["methods"][nm]=dict(xyz_mm=[np.round(coords[k//3]*1000,1).tolist() for k in ks],
                                comp=[COMP[k%3] for k in ks],
                                sigma_pred_um=scores([Wf[k] for k in ks])[0],AO_um={})
    out["methods"]["変位なし（温度2点のみ）"]=dict(xyz_mm=[],comp=[],sigma_pred_um=scores([])[0],AO_um={})
    out["methods"]["A・O 自身（循環・参考）"]=dict(xyz_mm=[np.round(coords[i]*1000,1).tolist() for i in (iA,iO)],
                                           comp=["Uz","Uz"],sigma_pred_um=scores([Wall[iA,2],Wall[iO,2]])[0],AO_um={})
    for case,lab in CASES:
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz"))
        t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]
        utrue=np.load(os.path.join(RES,f"truth_disp_all_{case}.npz"))["u"]
        aot=utrue[:,iA,2]-utrue[:,iO,2]; on=SCHED[case]
        for nm in out["methods"]:
            ks=methods.get(nm)
            obs=[divmod(k,3) for k in ks] if ks else ([(iA,2),(iO,2)] if "自身" in nm else [])
            v=[run(case,obs,s,t,T5,utrue,aot,on) for s in SEEDS]
            out["methods"][nm]["AO_um"][case]=float(np.mean(v))
            print(f"{case:13s} {nm:26s} 反りA−O {np.mean(v):.3f} µm",flush=True)
        rv=[]
        for ks in randoms:
            obs=[divmod(k,3) for k in ks]
            rv.append(float(np.mean([run(case,obs,s,t,T5,utrue,aot,on) for s in SEEDS])))
        out["random"][case]=dict(mean=float(np.mean(rv)),median=float(np.median(rv)),
                                 best=float(np.min(rv)),worst=float(np.max(rv)),n=N_RANDOM,values=rv)
        print(f"{case:13s} {'ランダム20通り':26s} 中央値 {np.median(rv):.3f}（最良 {np.min(rv):.3f} / 最悪 {np.max(rv):.3f}）µm",flush=True)
    json.dump(out,open(os.path.join(RES,"placement_method_comparison.json"),"w"),ensure_ascii=False,indent=1)

    # ---- 図 ----
    names=["変位なし（温度2点のみ）","ランダム（20通りの中央値）","熱感度の大きさ最大（素朴）",
           "Q-DEIM（QR列ピボット）","A最適（温度の trace 最小）","D最適（温度の det 最小）",
           "A・O 自身（循環・参考）","目的指向（本研究）"]
    cols=["#9AA5B1","#B0B7BF","#8E6FB0","#2E6FD8","#1F9D62","#E67E22","#BBBBBB","#C0392B"]
    fig,axs=plt.subplots(1,3,figsize=(16.5,5.2))
    for ax,(case,lab) in zip(axs,CASES):
        v=[out["random"][case]["median"] if nm.startswith("ランダム") else out["methods"][nm]["AO_um"][case] for nm in names]
        ax.barh(range(len(v)),v,color=cols)
        if True:
            r=out["random"][case]; i=names.index("ランダム（20通りの中央値）")
            ax.plot([r["best"],r["worst"]],[i,i],color="#44505c",lw=1.4)
            ax.plot([r["best"],r["worst"]],[i,i],"|",color="#44505c",ms=8)
        for i,x in enumerate(v): ax.text(x,i,f" {x:.3f}",va="center",fontsize=9.5)
        ax.set_yticks(range(len(v)))
        ax.set_yticklabels(names,fontsize=9) if ax is axs[0] else ax.set_yticklabels([])
        ax.invert_yaxis(); ax.set_xlabel("観測していない反り A−O の誤差 [µm]")
        ax.set_title(lab,fontsize=12); ax.grid(axis="x",alpha=.3)
    fig.suptitle("配置の選び方を比べる ― 同じ候補 14,754 通り・同じ holdout（A・O は一度も観測しない）"
                 "・真値＝OpenFOAM＋FrontISTR・5 seed 平均",fontsize=12.5)
    fig.tight_layout(rect=(0,0,1,0.91))
    fig.savefig(os.path.join(IMG,"placement_method_comparison.png"),dpi=150); plt.close(fig)
    print("wrote docs/img/placement_method_comparison.png")


if __name__=="__main__": main()
