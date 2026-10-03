"""A/O 以外の場所で変位を測り、観測していない A/O の変位・反りを当てられるかを調べる.

blog_004 §7-4 の構成は評価点 A・O をそのまま観測に使っており、循環だった。
ここでは観測を A/O 以外に限り、A・O は一度も観測しない（holdout）。

観測：温度2点（P2+P0、本番と同じ）＋ 変位2点（場所をいろいろ変える）
評価：観測していない Uz(A)、Uz(O)、反り A−O、および全20,696セルの温度場

変位の候補点は全5,040節点×3成分から選べる（results/dispop_allnodes.npz、FrontISTR 6回で作成）。
うち2点の組は、同化前に計算できる指標
    σ_S² = wᵀ( B − B Hᵀ(H B Hᵀ+R)⁻¹ H B ) w        （w＝A−O の熱感度）
が小さいものを貪欲に選ぶ。真値も同じ ROM の双子実験。
EnKF・60メンバー・30秒ごと20回・5 seed・加熱期30〜300秒。

出力: results/disp_elsewhere_predict_AO.json, docs/img/disp_elsewhere_predict_AO.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/disp_elsewhere_predict_AO.py
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
TNODES=[2,0]; COMP=["Ux","Uy","Uz"]; N_PRIOR=4000; PRIOR_SEED=20260930
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005)
EXCL_R=0.012      # A/O からこの距離[m]以内は「A/Oを測ったのと同じ」とみなして候補から除く


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); A=U@UPp; ncell=U.shape[0]; G=A.T@A
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; D=op["D"]; coords=op["coords"]; n=coords.shape[0]
    Wall=np.einsum("nck,kj->ncj",D,UPp)
    iA=int(np.linalg.norm(coords-np.array(A_XYZ),axis=1).argmin())
    iO=int(np.linalg.norm(coords-np.array(O_XYZ),axis=1).argmin())
    wAO=Wall[iA,2]-Wall[iO,2]; uAO0=um[iA,2]-um[iO,2]
    print(f"評価点 A={np.round(coords[iA]*1000,1)} O={np.round(coords[iO]*1000,1)}  W(A-O)={np.round(wAO,3)}")
    # 候補から A/O 近傍を除外
    far=(np.linalg.norm(coords-np.array(A_XYZ),axis=1)>EXCL_R)&(np.linalg.norm(coords-np.array(O_XYZ),axis=1)>EXCL_R)
    # 事前共分散
    rng=np.random.default_rng(PRIOR_SEED)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Bs=[]; tp=0.0
    for tb in np.arange(OBS_DT,300.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,tp,tb,DT); tp=tb; Bs.append(np.cov(Tp.T))
    Ht=np.eye(NPT)[TNODES]
    def post(B,ws):
        H=np.vstack([Ht]+[w[None,:] for w in ws])
        R=np.diag([SIG_T**2]*len(TNODES)+[SIG_U**2]*len(ws))
        K=B@H.T@np.linalg.inv(H@B@H.T+R); return B-K@H@B
    sig=lambda ws,w: float(np.mean([np.sqrt(max(w@post(B,ws)@w,0)) for B in Bs]))
    fld=lambda ws: float(np.mean([np.sqrt(float(np.einsum("ij,ji->",post(B,ws),G))/ncell) for B in Bs]))
    # 反り A−O を狙って変位2点を貪欲に選ぶ（A/O 近傍は除外）
    Wf=Wall.reshape(-1,NPT); idx=np.where(np.repeat(far,3))[0]
    pick=[]
    for step in range(2):
        vals=np.array([sig(pick+[Wf[k]],wAO) for k in idx])
        k=int(idx[int(vals.argmin())]); pick.append(Wf[k])
        i,c=divmod(k,3); print(f"  選定{step+1}点目: {np.round(coords[i]*1000,1)} の {COMP[c]}  σ(A−O)={vals.min():.4f} µm",flush=True)
        idx=idx[idx!=k]
    sel=[divmod(int(np.where((Wf==w).all(1))[0][0]),3) for w in pick]
    def near(x): return int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    hi=np.load(os.path.join(RES,"dispop_hiW.npz"))
    cfgs=[("変位なし（温度2点のみ）",[]),
          ("A・O 自身を観測（循環・参考）",[(iA,2),(iO,2)]),
          ("B・C：中段 z=75 mm の Uz",[(near((0.028,0,0.075)),2),(near((-0.028,0,0.075)),2)]),
          (f"選定2点：{np.round(coords[sel[0][0]]*1000,1).tolist()} {COMP[sel[0][1]]} と {np.round(coords[sel[1][0]]*1000,1).tolist()} {COMP[sel[1][1]]}",sel),
          ("側面中段 z=50 mm の Ux（半径方向）",[(near((0.0375,0,0.050)),0),(near((-0.0375,0,0.050)),0)]),
          ("底面近く z=5 mm の Uz（悪い例）",[(near((0.0375,0,0.005)),2),(near((-0.0375,0,0.005)),2)])]
    field=lambda T5: mean+A@(T5-mean[pod])
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    def run(obs,seed):
        ws=[Wall[i,c] for i,c in obs]; u0=[um[i,c] for i,c in obs]
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*len(TNODES)+[SIG_U**2]*len(ws)); eF=[];eAO=[];eA=[];eO=[]; tp=0.0
        aA=um[iA,2]; aO=um[iO,2]; wA=Wall[iA,2]; wO=Wall[iO,2]
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][TNODES]); Yf=Z[:,TNODES]
            for w,b0 in zip(ws,u0):
                yv.append(b0+w@(Ttr[ci]-mean[pod])); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mean[pod])@w])
            y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                m=Z[:,:NPT].mean(0); dm=m-mean[pod]; dt=Ttr[ci]-mean[pod]
                eF.append(np.sqrt(((field(m)-field(Ttr[ci]))**2).mean()))
                eA.append(abs((aA+wA@dm)-(aA+wA@dt))); eO.append(abs((aO+wO@dm)-(aO+wO@dt)))
                eAO.append(abs(wAO@dm-wAO@dt))
        return [np.mean(eF),np.mean(eAO),np.mean(eA),np.mean(eO)]
    out={"note":"観測＝温度2点(P2+P0)＋変位2点。A・O は（循環の構成を除き）一度も観測していない。誤差は加熱期30〜300秒・5 seed平均",
         "A_xyz_mm":np.round(coords[iA]*1000,1).tolist(),"O_xyz_mm":np.round(coords[iO]*1000,1).tolist(),"configs":{}}
    for name,obs in cfgs:
        r=np.array([run(obs,s) for s in SEEDS]); m=r.mean(0); sd=r.std(0,ddof=1)
        out["configs"][name]=dict(obs_xyz_mm=[np.round(coords[i]*1000,1).tolist() for i,_ in obs],
                                  obs_comp=[COMP[c] for _,c in obs],
                                  field_K=float(m[0]),AO_um=float(m[1]),AO_um_sd=float(sd[1]),
                                  A_um=float(m[2]),O_um=float(m[3]),
                                  sigma_pred_AO=(sig([Wall[i,c] for i,c in obs],wAO) if obs else sig([],wAO)))
        print(f"{name[:46]:48s} 温度場 {m[0]:.3f} K  反りA−O {m[1]:.3f}±{sd[1]:.3f}  Uz(A) {m[2]:.3f}  Uz(O) {m[3]:.3f} µm",flush=True)
    json.dump(out,open(os.path.join(RES,"disp_elsewhere_predict_AO.json"),"w"),ensure_ascii=False,indent=1)
    names=list(out["configs"]); cols=["#9AA5B1","#BBBBBB","#2E8B57","#C0392B","#2E6FD8","#E67E22"]
    fig,axs=plt.subplots(1,3,figsize=(16.5,5.0))
    for ax,key,ttl,un in [(axs[0],"AO_um","観測していない反り A−O","µm"),
                          (axs[1],"A_um","観測していない Uz(A)","µm"),(axs[2],"field_K","温度場全体","K")]:
        v=[out["configs"][nm][key] for nm in names]
        ax.barh(range(len(v)),v,color=cols[:len(v)])
        for i,x in enumerate(v): ax.text(x,i,f" {x:.3f}",va="center",fontsize=10)
        ax.set_yticks(range(len(v)))
        ax.set_yticklabels([nm.replace("：","\n")[:34] for nm in names],fontsize=8.5) if ax is axs[0] else ax.set_yticklabels([])
        ax.invert_yaxis(); ax.set_xlabel(f"誤差 [{un}]"); ax.set_title(ttl,fontsize=12); ax.grid(axis="x",alpha=.3)
    fig.suptitle("A・O 以外で変位を測り、観測していない A・O を当てられるか（温度2点 P2+P0 は共通、5 seed 平均）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.91)); fig.savefig(os.path.join(IMG,"disp_elsewhere_predict_AO.png"),dpi=150); plt.close(fig)
    print("wrote figure")


if __name__=="__main__": main()
