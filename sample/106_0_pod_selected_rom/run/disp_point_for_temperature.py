"""どこに変位計を置くと「温度推定」が良くなるかを、全5,040節点から探す.

変位計は本来「変位を測る」ためのものだが、本研究では変位が5点温度の一次式なので、
変位計は温度場の推定にも効く（blog_006 §7-5b）。では、どこに置くのが良いか。

手順
  1. 全FEM節点×3成分（5,040×3＝15,120候補）の 温度→変位 写像を使う（results/dispop_allnodes.npz）
  2. 同化前に計算できる指標でふるいにかける
       σ_S² = tr( A (B - B Hᵀ(H B Hᵀ+R)⁻¹ H B) Aᵀ ) / n_cell     ← 全温度場の残る分からなさ
     A は 5点温度→全20,696セル の復元行列、B は温度のばらつき（ROMで作る）
  3. 上位・下位・代表的な候補について、実際に EnKF を回して温度場RMSEを確かめる

観測は「温度2点(P2+P0、本番と同じ) ＋ 変位1点」。変位計1本を足す効果だけを見る。
真値も同じ ROM（双子実験）。5 seed 平均、加熱期30〜300秒。

出力: results/disp_point_for_temperature.json, docs/img/disp_point_for_temperature.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/disp_point_for_temperature.py
"""
from __future__ import annotations
import os, sys, json, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TNODES=[2,0]                 # 温度センサ（本番と同じ P2, P0）
N_PRIOR=4000; PRIOR_SEED=20260930
COMP=["Ux","Uy","Uz"]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); ncell=U.shape[0]
    A=U@UPp                                        # 5点温度 → 全セル温度（復元行列、20,696×5）
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; D=op["D"]; coords=op["coords"]           # (n,3), (n,3,5), (n,3)
    n=coords.shape[0]
    Wall=np.einsum("nck,kj->ncj",D,UPp)            # (n,3,5)  温度1K→その成分の変位 [µm/K]
    # ── 事前共分散 B（時刻ごと。名目の前提）──
    rng=np.random.default_rng(PRIOR_SEED)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Bs=[]; tp=0.0
    for tb in np.arange(OBS_DT,300.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,C,Km,hp,qp,heat,tp,tb,DT); tp=tb; Bs.append(np.cov(Tp.T))
    Ht=np.eye(NPT)[TNODES]
    def post_field_sd(B,w=None):
        """温度2点（＋変位1点 w）を測った後に残る、全温度場の1セルあたり標準偏差 [K]"""
        H=Ht if w is None else np.vstack([Ht,w])
        R=np.diag([SIG_T**2]*len(TNODES)+([SIG_U**2] if w is not None else []))
        K=B@H.T@np.linalg.inv(H@B@H.T+R); P=B-K@H@B
        return np.sqrt(np.einsum("ij,jk,ik->",A,P,A)/ncell)
    base=float(np.mean([post_field_sd(B) for B in Bs]))
    # ── 全候補を指標でふるい分け ──
    score=np.zeros((n,3))
    for c in range(3):
        for i in range(n):
            score[i,c]=np.mean([post_field_sd(B,Wall[i,c]) for B in Bs])
    print(f"[scan] 温度2点のみの σ(全温度場) = {base:.4f} K   候補 {n*3} 個を評価",flush=True)
    flat=score.reshape(-1); order=np.argsort(flat)
    def info(k):
        i,c=divmod(int(k),3)
        return dict(node=i,comp=COMP[c],xyz_mm=np.round(coords[i]*1000,1).tolist(),
                    sigma=float(flat[k]),W=np.round(Wall[i,c],3).tolist())
    best=[info(k) for k in order[:6]]; worst=[info(k) for k in order[-3:]]
    for b in best: print("  良い候補",b["xyz_mm"],b["comp"],f"σ={b['sigma']:.4f}",flush=True)
    # ── 実際に EnKF で確かめる候補 ──
    def nearest(xyz): return int(np.linalg.norm(coords-np.array(xyz),axis=1).argmin())
    AO=[nearest((0.028,0,0.1005)),nearest((-0.028,0,0.1005))]
    cands=[("温度2点のみ（変位なし）",None,None)]
    for j,b in enumerate(best[:3]):
        cands.append((f"予測1〜3位：{b['xyz_mm']} の {b['comp']}",b["node"],COMP.index(b["comp"])))
    cands.append((f"予測最下位：{worst[-1]['xyz_mm']} の {worst[-1]['comp']}",worst[-1]["node"],COMP.index(worst[-1]["comp"])))
    cands.append((f"上面A の Uz（本番で使った点）",AO[0],2))
    cands.append((f"上面O の Uz（本番で使った点）",AO[1],2))
    cands.append(("底面近く(37.5,0,5)の Uz",nearest((0.0375,0,0.005)),2))
    field=lambda T5: mean+A@(T5-mean[pod])
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    def run(node,comp,seed):
        use=node is not None
        w=Wall[node,comp] if use else None; u0=um[node,comp] if use else 0.0
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*len(TNODES)+([SIG_U**2] if use else [])); e=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][TNODES]); Yf=Z[:,TNODES]
            if use:
                yv.append(u0+w@(Ttr[ci]-mean[pod])); Yf=np.column_stack([Yf,u0+(Z[:,:NPT]-mean[pod])@w])
            y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                m=Z[:,:NPT].mean(0); e.append(np.sqrt(((field(m)-field(Ttr[ci]))**2).mean()))
        return float(np.mean(e))
    res=[]
    for name,node,comp in cands:
        r=np.array([run(node,comp,s) for s in SEEDS])
        sg=float(np.mean([post_field_sd(B,Wall[node,comp]) for B in Bs])) if node is not None else base
        res.append(dict(name=name,xyz_mm=(np.round(coords[node]*1000,1).tolist() if node is not None else None),
                        comp=(COMP[comp] if node is not None else None),sigma_pred=sg,
                        field_K=float(r.mean()),field_K_sd=float(r.std(ddof=1))))
        print(f"  {name:42s} 予測σ {sg:.4f}  実測 温度場 {r.mean():.3f}±{r.std(ddof=1):.3f} K",flush=True)
    out=dict(note="観測＝温度2点(P2+P0)＋変位1点。温度場RMSEは全20,696セル・加熱期30〜300秒・5 seed平均",
             base_sigma=base,best=best,worst=worst,results=res)
    json.dump(out,open(os.path.join(RES,"disp_point_for_temperature.json"),"w"),ensure_ascii=False,indent=1)
    # ── 図 ──
    fig=plt.figure(figsize=(16.5,5.4)); gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.35],wspace=0.30)
    xx,yy,zz=coords[:,0]*1000,coords[:,1]*1000,coords[:,2]*1000
    vmin,vmax=float(score[:,2].min()),float(score[:,2].max())
    top=zz>100.0
    ax=fig.add_subplot(gs[0,0])
    sc=ax.scatter(xx[top],yy[top],c=score[top,2],s=46,cmap="viridis_r",vmin=vmin,vmax=vmax)
    ax.plot(28,0,"r*",ms=16,label="本番の変位計 A"); ax.plot(-28,0,"b*",ms=16,label="本番の変位計 O")
    bi=best[0]["node"]; ax.plot(coords[bi,0]*1000,coords[bi,1]*1000,"w^",ms=13,mec="k",label="予測1位")
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]"); ax.legend(fontsize=9,loc="upper right")
    ax.set_title("上面（z=100.5 mm）を真上から見た図\nヒータは +x 側",fontsize=11.5)
    plt.colorbar(sc,ax=ax,label="σ [K]（小さいほど良い）")
    ax=fig.add_subplot(gs[0,1])
    sl=np.abs(yy)<6
    sc=ax.scatter(xx[sl],zz[sl],c=score[sl,2],s=36,cmap="viridis_r",vmin=vmin,vmax=vmax)
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("z [mm]")
    ax.set_title("y≈0 の断面（横から見た図）",fontsize=11.5)
    plt.colorbar(sc,ax=ax,label="σ [K]")
    ax=fig.add_subplot(gs[0,2])
    nm=[r["name"].replace("：","\n") for r in res]; v=[r["field_K"] for r in res]; sd=[r["field_K_sd"] for r in res]
    cols=["#9AA5B1"]+["#2E8B57"]*3+["#C0392B"]+["#2E6FD8"]*2+["#E67E22"]
    ax.barh(range(len(v)),v,xerr=sd,color=cols[:len(v)],capsize=3)
    for i,(x,e_) in enumerate(zip(v,sd)): ax.text(x+e_+0.01,i,f"{x:.3f}",va="center",fontsize=10)
    ax.set_yticks(range(len(v))); ax.set_yticklabels(nm,fontsize=8.5); ax.invert_yaxis()
    ax.set_xlim(0,0.78); ax.set_xlabel("温度場RMSE [K]（加熱期・5 seed平均）"); ax.grid(axis="x",alpha=.3)
    ax.set_title("実際に同化して確かめた結果",fontsize=11.5)
    fig.suptitle("どこに変位計（上下方向 Uz）を置くと温度推定が良くなるか（温度2点 P2+P0 は固定、変位計1本を追加）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"disp_point_for_temperature.png"),dpi=150); plt.close(fig)
    print("wrote figure")


if __name__=="__main__": main()
