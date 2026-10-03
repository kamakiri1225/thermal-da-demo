"""blog_008 §6 用：変位2点をどう選んだかを、1点ずつ数値と図で追えるようにする.

手順（run/disp_elsewhere_predict_AO.py と同じ）
  1点目：15,120候補すべてについて σ（測った後に残る反りの分からなさ）を計算し、最小を選ぶ
  2点目：1点目を測った前提で、残りの候補について σ を計算し直し、最小を選ぶ
ここでは各段階の上位・下位の候補と σ を記録し、σ の分布図（1点目／2点目）も出す。

出力: results/explain_disp_selection.json, docs/img/disp_selection_steps.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/explain_disp_selection.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; DT=2.0; OBS_DT=30.0; SIG_T=0.30; SIG_U=0.30; N_PRIOR=4000; PRIOR_SEED=20260930
TNODES=[2,0]; COMP=["Ux","Uy","Uz"]
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005); EXCL_R=0.012


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    coords=op["coords"]; Wall=np.einsum("nck,kj->ncj",op["D"],UPp); n=coords.shape[0]
    iA=int(np.linalg.norm(coords-np.array(A_XYZ),axis=1).argmin())
    iO=int(np.linalg.norm(coords-np.array(O_XYZ),axis=1).argmin())
    w=Wall[iA,2]-Wall[iO,2]
    far=(np.linalg.norm(coords-np.array(A_XYZ),axis=1)>EXCL_R)&(np.linalg.norm(coords-np.array(O_XYZ),axis=1)>EXCL_R)
    rng=np.random.default_rng(PRIOR_SEED)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Bs=[]; tp=0.0
    for tb in np.arange(OBS_DT,300.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,tp,tb,DT); tp=tb; Bs.append(np.cov(Tp.T))
    Ht=np.eye(NPT)[TNODES]
    def post_after(B,ws):
        H=np.vstack([Ht]+[x[None,:] for x in ws])
        R=np.diag([SIG_T**2]*len(TNODES)+[SIG_U**2]*len(ws))
        K=B@H.T@np.linalg.inv(H@B@H.T+R); return B-K@H@B
    Wf=Wall.reshape(-1,NPT); mask=np.repeat(far,3)
    def scan(ws):
        """候補を1本足したときの σ を全候補について一括計算"""
        acc=np.zeros(Wf.shape[0])
        for B in Bs:
            P=post_after(B,ws)
            Pw=Wf@P; den=np.einsum("ij,ij->i",Pw,Wf)+SIG_U**2
            base=float(w@P@w); corr=(Pw@w)**2/den
            acc+=np.sqrt(np.maximum(base-corr,0.0))
        s=acc/len(Bs); s[~mask]=np.inf
        return s
    out={"note":"σ＝温度2点(P2+P0)と、選んだ変位を測った後に残る『反りA−Oの分からなさ』[µm]。加熱期10時刻の平均",
         "sigma_temp_only":float(np.mean([np.sqrt(w@post_after(B,[])@w) for B in Bs])),
         "W_AO":w.tolist(),"n_candidates":int(mask.sum()),"steps":[]}
    print(f"温度2点だけのとき σ = {out['sigma_temp_only']:.4f} µm   候補 {int(mask.sum())} 個")
    ws=[]; s1=None
    for step in range(2):
        s=scan(ws)
        if step==0: s1=s.copy()
        order=np.argsort(s)
        def info(k):
            i,c=divmod(int(k),3)
            return dict(xyz_mm=np.round(coords[i]*1000,1).tolist(),comp=COMP[c],sigma=float(s[k]),
                        W=float(np.abs(Wall[i,c]).sum()))
        fin=np.isfinite(s)
        rec=dict(step=step+1,best=[info(k) for k in order[:5]],
                 worst=[info(k) for k in order[fin.sum()-3:fin.sum()]],
                 median=float(np.median(s[fin])))
        out["steps"].append(rec)
        print(f"--- {step+1}点目 ---")
        for b in rec["best"][:3]: print(f"  良い: {b['xyz_mm']} {b['comp']}  σ={b['sigma']:.4f}")
        print(f"  中央値 σ={rec['median']:.4f}   悪い: {rec['worst'][-1]['xyz_mm']} {rec['worst'][-1]['comp']} σ={rec['worst'][-1]['sigma']:.4f}")
        k=int(order[0]); ws.append(Wf[k]); mask=mask.copy(); mask[k]=False
        i,c=divmod(k,3); rec["chosen"]=dict(xyz_mm=np.round(coords[i]*1000,1).tolist(),comp=COMP[c],sigma=float(s[k]))
    s2=scan(ws[:1])
    json.dump(out,open(os.path.join(RES,"explain_disp_selection.json"),"w"),ensure_ascii=False,indent=1)
    # 図：1点目と2点目の σ 分布（上面）
    xx,yy,zz=coords[:,0]*1000,coords[:,1]*1000,coords[:,2]*1000
    top=zz>100.0
    fig,axs=plt.subplots(1,3,figsize=(16.5,5.0),gridspec_kw=dict(width_ratios=[1,1,1.1],wspace=0.30))
    for ax,sc_,ttl in [(axs[0],s1,"① 1点目を選ぶとき"),(axs[1],s2,"② 1点目を測った前提で2点目を選ぶとき")]:
        v=sc_.reshape(n,3)[:,2]      # Uz 成分で地図を描く
        vv=np.where(np.isfinite(v),v,np.nan)
        m=ax.scatter(xx[top],yy[top],c=vv[top],s=46,cmap="viridis_r")
        ax.plot(28.7,0,"k*",ms=15); ax.plot(-28.7,0,"k*",ms=15)
        ax.annotate("評価点 A\n（測らない）",(28.7,0),textcoords="offset points",xytext=(-38,16),fontsize=9)
        ax.annotate("評価点 O\n（測らない）",(-28.7,0),textcoords="offset points",xytext=(-50,-30),fontsize=9)
        ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
        ax.set_title(f"{ttl}\n上面 z=100.5 mm の Uz を測る場合",fontsize=11.5)
        plt.colorbar(m,ax=ax,label="σ [µm]（小さいほど良い）")
    for k,st in enumerate(out["steps"]):
        p=st["chosen"]["xyz_mm"]; axs[k].plot(p[0],p[1],"^",ms=14,color="w",mec="k")
        axs[k].annotate(f"選んだ点\n{st['chosen']['comp']}",(p[0],p[1]),textcoords="offset points",xytext=(6,10),fontsize=9.5,weight="bold")
    ax=axs[2]
    lab=["温度2点だけ"]+[f"{s['chosen']['xyz_mm']}\n{s['chosen']['comp']} を足す" for s in out["steps"]]
    val=[out["sigma_temp_only"]]+[s["chosen"]["sigma"] for s in out["steps"]]
    ax.plot(range(3),val,"o-",color="#C0392B",lw=2.4,ms=10)
    for i,v_ in enumerate(val): ax.annotate(f"{v_:.3f} µm",(i,v_),textcoords="offset points",xytext=(8,8),fontsize=11)
    ax.set_xticks(range(3)); ax.set_xticklabels(lab,fontsize=9.5); ax.set_ylim(0,0.95)
    ax.set_ylabel("残る「反りの分からなさ」σ [µm]"); ax.grid(alpha=.3)
    ax.set_title("1点足すごとに σ が下がる",fontsize=11.5)
    fig.suptitle("変位2点の選び方：測った後に残る「反りの分からなさ」σ が最小の点を1本ずつ足す（同化前の計算のみ）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.91)); fig.savefig(os.path.join(IMG,"disp_selection_steps.png"),dpi=150); plt.close(fig)
    print("wrote disp_selection_steps.png")


if __name__=="__main__": main()
