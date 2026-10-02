"""温度センサ2本の配置を「熱感度W＋温度のばらつき（共分散B）」で選べるかを調べる.

問い：A−Oの反り（相対変位）を当てたい。温度センサ2本をどこに置けばよいか。
  ・熱感度 W（FrontISTR, 温度→変位）だけでは決まらない
  ・反りの誤差 ＝ W ×「測った後に残る温度の分からなさ」
  ・そこで σ_S² = wᵀ(B − B H_Sᵀ(H_S B H_Sᵀ+R)⁻¹ H_S B) w を同化前に計算し、
    実際に EnKF を回した結果の順位と比べる
  ・さらに、真値の条件（発熱量・放熱・熱源位置・初期温度）を前提から外しても良いままかを調べる

真値も同化モデルも同じROM（双子実験）。OpenFOAM＋FrontISTRを真値にした試験ではない。

出力: results/sensor_pair_goal_oriented.json
      docs/img/sensor_pair_*.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensor_pair_goal_oriented.py
"""
from __future__ import annotations
import os, sys, json, itertools, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from scipy.stats import spearmanr
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
OUT=os.path.join(RES,"sensor_pair_goal_oriented.json")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
N_PRIOR=4000; PRIOR_SEED=20260930
iA,iO=2,5                                   # dispop_xyz_4pts の Uz(A), Uz(O)
PAIRS=[list(c) for c in itertools.combinations(range(NPT),2)]
lab=lambda S:"+".join(f"P{i}" for i in S)
SIDE={0:"A",1:"O",2:"A",3:"A",4:"O"}         # x座標の正負（A=+X, O=−X）
SPLIT={lab(S) for S in PAIRS if len({SIDE[i] for i in S})==2}


def load():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    m=dict(C=d["C"],Km=rg.tri_to_matrix(d["K_upper"],NPT),h0=float(d["h"]),heat=int(d["heat_node"]),xyz=d["xyz"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    # 変位を評価する全14点（Uz）。どれも POD モード空間の FrontISTR 応答
    ops=[(um[[2,5,8,11]],D[[2,5,8,11]])]
    for f,k1,k2 in [("dispop_hiW","uz_mean","D"),("dispop_loW","uz_mean","D"),
                    ("dispop_unobs_z75","uz_mean","Dmode"),("dispop_holdout_support_a4dadbb851dca4d7","uz_mean","D")]:
        z=np.load(os.path.join(RES,f+".npz")); ops.append((z[k1],z[k2]))
    UM=np.concatenate([o[0] for o in ops]); DM=np.vstack([o[1] for o in ops])
    m.update(U=U,mean=mean,pod=pod,UPp=UPp,um=um,D=D,UM=UM,DM=DM)
    m["w"]=(D[iA]-D[iO])@UPp                   # 熱感度W：反りA−Oに対する5点温度の重み [µm/K]
    return m


def prior_cov(m):
    """名目の前提（発熱0.3〜1.8倍、放熱N(0.02,0.01)、初期−3〜+12K）でROMを4000本回し、時刻ごとの共分散Bを作る."""
    rng=np.random.default_rng(PRIOR_SEED)
    Tp=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_PRIOR,NPT))
    qp=rng.uniform(0.3,1.8,N_PRIOR); hp=np.clip(rng.normal(0.02,0.01,N_PRIOR),1e-3,0.1)
    Bs={}; tp=0.0
    for tb in np.arange(OBS_DT,300.0+1e-9,OBS_DT):
        Tp=rg.integrate_ensemble(Tp,m["C"],m["Km"],hp,qp,m["heat"],tp,tb,DT); tp=tb; Bs[float(tb)]=np.cov(Tp.T)
    return Bs


def post_cov(B,S):
    if not S: return B
    H=np.eye(NPT)[S]; K=B@H.T@np.linalg.inv(H@B@H.T+SIG_T**2*np.eye(len(S)))
    return B-K@H@B


def da_run(m,Ttr,S,seed):
    disp=lambda T5: m["um"]+((T5-m["mean"][m["pod"]])@m["UPp"].T)@m["D"].T
    uz=lambda T5: m["UM"]+((T5-m["mean"][m["pod"]])@m["UPp"].T)@m["DM"].T
    field=lambda T5: m["mean"]+m["U"]@(m["UPp"]@(T5-m["mean"][m["pod"]]))
    rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
    Z=np.zeros((N_ENS,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.eye(len(S))*SIG_T**2; eF=[];eU=[];eAO=[]; tp=0.0
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    for ci,tb in enumerate(cyc,1):
        Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],m["C"],m["Km"],Z[:,IH],Z[:,IQ],m["heat"],tp,tb,DT); tp=tb
        y=Ttr[ci][S]+ro.normal(0,SIG_T,len(S))
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Z[:,S])
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        if tb<=300:                                   # 加熱期 30〜300 s
            mu=Z[:,:NPT].mean(0); uu=disp(mu); ut=disp(Ttr[ci])
            eF.append(np.sqrt(((field(mu)-field(Ttr[ci]))**2).mean()))
            eU.append(np.abs(uz(mu)-uz(Ttr[ci])).mean())
            eAO.append(abs((uu[iA]-uu[iO])-(ut[iA]-ut[iO])))
    return [float(np.mean(eF)),float(np.mean(eU)),float(np.mean(eAO))]


def truth(m,q,h,heat,dT0):
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K+dT0); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,m["C"],m["Km"],h,q,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    return np.array(Ttr)


def compute():
    m=load(); Bs=prior_cov(m); w=m["w"]
    pred={lab(S):float(np.mean([np.sqrt(w@post_cov(B,S)@w) for B in Bs.values()])) for S in PAIRS}
    pred_single={f"P{i}":float(np.mean([np.sqrt(w@post_cov(B,[i])@w) for B in Bs.values()])) for i in range(NPT)}
    B150=Bs[150.0]
    resid={}
    for S in [[],[0,2],[2,4]]:
        P=post_cov(B150,S)
        resid[lab(S) or "none"]=dict(sd=np.sqrt(np.diag(P)).tolist(),sigma=float(np.sqrt(w@P@w)))
    h0=m["h0"]; hn=m["heat"]
    SCEN=[("前提どおり",1.0,h0,hn,0.0),("発熱量 2.5倍",2.5,h0,hn,0.0),("発熱量 0.2倍",0.2,h0,hn,0.0),
          ("放熱 3倍",1.0,3*h0,hn,0.0),("放熱 0.3倍",1.0,0.3*h0,hn,0.0),
          ("初期温度 +15 K",1.0,h0,hn,15.0),("熱源がP0（A側）",1.0,h0,0,0.0),("熱源がP4（O側）",1.0,h0,4,0.0)]
    scen={}
    for name,q,h,heat,dT0 in SCEN:
        Ttr=truth(m,q,h,heat,dT0)
        scen[name]={lab(S):[da_run(m,Ttr,S,s) for s in SEEDS] for S in PAIRS}
        print(f"[pair] {name} done",flush=True)
    singles={}
    Ttr=truth(m,1.0,h0,hn,0.0)
    for i in range(NPT): singles[f"P{i}"]=[da_run(m,Ttr,[i],s) for s in SEEDS]
    out=dict(note="metrics per seed = [温度場全20,696セルRMSE K, 変位14点Uz平均絶対誤差 µm, 反りA−O絶対誤差 µm], 30〜300 s 平均",
             W_AO=w.tolist(),W_A=(m["D"][iA]@m["UPp"]).tolist(),W_O=(m["D"][iO]@m["UPp"]).tolist(),
             xyz=m["xyz"].tolist(),heat_node=hn,pred_sigma=pred,pred_sigma_single=pred_single,
             residual_150s=resid,scenarios=scen,singles=singles,split=sorted(SPLIT))
    json.dump(out,open(OUT,"w"),ensure_ascii=False,indent=1)
    return out


def figures(o):
    names=[lab(S) for S in PAIRS]
    nom={k:np.array(v) for k,v in o["scenarios"]["前提どおり"].items()}
    w=np.array(o["W_AO"])
    green="#2E8B57"; gray="#9AA5B1"
    # ── 図1：10通りの比較（温度場と反り、seedのばらつき付き）──
    order=sorted(names,key=lambda k:nom[k][:,2].mean())
    fig,axs=plt.subplots(1,2,figsize=(13,5.4),sharey=True)
    for ax,j,tt,un in [(axs[0],0,"温度場全体（20,696セル）の誤差","K"),(axs[1],2,"反り A−O の誤差","µm")]:
        mu=[nom[k][:,j].mean() for k in order]; sd=[nom[k][:,j].std(ddof=1) for k in order]
        cols=[green if k in o["split"] else gray for k in order]
        ax.barh(range(len(order)),mu,xerr=sd,color=cols,capsize=3)
        for y,(v,s) in enumerate(zip(mu,sd)): ax.text(v+s+0.01,y,f"{v:.2f}",va="center",fontsize=10)
        ax.set_title(tt,fontsize=13); ax.set_xlabel(f"0〜300 s の平均誤差 [{un}]（5 seed平均±標準偏差）")
        ax.grid(axis="x",alpha=.3)
    axs[0].set_yticks(range(len(order))); axs[0].set_yticklabels([f"{k}（{'両側' if k in o['split'] else '同じ側'}）" for k in order],fontsize=11)
    axs[0].invert_yaxis()
    fig.suptitle("温度センサ2本の全10通り：緑＝A側とO側に1本ずつ、灰＝同じ側に2本",fontsize=14)
    fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(os.path.join(IMG,"sensor_pair_compare.png"),dpi=150); plt.close(fig)
    # ── 図2：反りの誤差＝W×温度の分からなさ ──
    r=o["residual_150s"]
    fig,axs=plt.subplots(1,2,figsize=(13,4.8),gridspec_kw=dict(width_ratios=[1,1.5]))
    ax=axs[0]; cols=["#C0392B" if v>0 else "#2E6FD8" for v in w]
    ax.bar(range(NPT),w,color=cols)
    for i,v in enumerate(w): ax.text(i,v+(0.03 if v>=0 else -0.08),f"{v:+.2f}",ha="center",fontsize=11)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(range(NPT)); ax.set_xticklabels([f"P{i}\n({SIDE[i]}側)" for i in range(NPT)])
    ax.set_ylabel("熱感度 W [µm/K]"); ax.set_title("① 熱感度 W：その点が1 K上がると反りが何µm変わるか",fontsize=12)
    ax.set_ylim(-0.85,0.85); ax.grid(axis="y",alpha=.3)
    ax=axs[1]; x=np.arange(NPT); bw=0.27
    for k,(key,lb,c) in enumerate([("none","測らない","#BBBBBB"),("P0+P2","P0＋P2（A側だけ）","#E67E22"),("P2+P4","P2＋P4（両側）",green)]):
        sd=r[key]["sd"]; ax.bar(x+(k-1)*bw,sd,bw,color=c,label=f"{lb}：反りの分からなさ {r[key]['sigma']:.2f} µm")
    ax.set_xticks(x); ax.set_xticklabels([f"P{i}\nW={w[i]:+.2f}" for i in range(NPT)])
    ax.set_ylabel("測った後に残る温度の分からなさ [K]"); ax.set_ylim(0,3.6); ax.legend(fontsize=10,loc="upper center",ncol=1)
    ax.set_title("② W の大きい点の温度が決まっているか（150 s、1回観測）",fontsize=12); ax.grid(axis="y",alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"sensor_pair_residual.png"),dpi=150); plt.close(fig)
    # ── 図3：予測σと実際の誤差 ──
    p=np.array([o["pred_sigma"][k] for k in names]); a=np.array([nom[k][:,2].mean() for k in names])
    rho=spearmanr(p,a).correlation
    fig,ax=plt.subplots(figsize=(7.5,5.6))
    for k,x_,y_ in zip(names,p,a):
        ax.scatter(x_,y_,s=110,color=green if k in o["split"] else gray,edgecolor="k",zorder=3)
        ax.annotate(k,(x_,y_),textcoords="offset points",xytext=(7,4),fontsize=10)
    ax.set_xlabel("同化前に計算した「反りの分からなさ」σ [µm]（W＋共分散B）")
    ax.set_ylabel("実際に同化した反りの誤差 [µm]（EnKF・5 seed）")
    ax.set_title(f"同化を回す前に順位を当てられるか（順位相関 {rho:.2f}）",fontsize=13); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"sensor_pair_pred_vs_actual.png"),dpi=150); plt.close(fig)
    # ── 図4：前提が外れた真値（行ごとの順位で色付け）──
    sc=list(o["scenarios"].keys()); best=min(o["pred_sigma"],key=o["pred_sigma"].get)
    M=np.array([[np.mean(np.array(o["scenarios"][s][k])[:,2]) for k in names] for s in sc])
    R=np.argsort(np.argsort(M,axis=1),axis=1)+1
    fig,ax=plt.subplots(figsize=(13.5,5.2))
    ax.imshow(R,cmap="RdYlGn_r",aspect="auto",vmin=1,vmax=10)
    for i in range(len(sc)):
        for j in range(len(names)):
            ax.text(j,i,f"{M[i,j]:.2f}\n({R[i,j]}位)",ha="center",va="center",fontsize=9,
                    fontweight="bold" if names[j]==best else None)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([("★" if k==best else "")+k+("\n両側" if k in o["split"] else "\n同じ側") for k in names],fontsize=10)
    ax.set_yticks(range(len(sc))); ax.set_yticklabels(sc,fontsize=11)
    ax.set_title(f"前提から外れた真値での反りA−Oの誤差 [µm]と順位（★＝前提どおりの条件で選んだ {best}）",fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"sensor_pair_robustness.png"),dpi=150); plt.close(fig)
    # ── 図5：3D位置（Wの符号で色分け）──
    import pyvista as pv, vtk, cylinder_mesh
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([x for _n,x in mesh["nodes"]]); idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    xyz=np.array(o["xyz"]); Hh=0.1005; A=np.array([0.028,0,Hh]); O=np.array([-0.028,0,Hh])
    pl=pv.Plotter(off_screen=True,window_size=(1200,760))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.35)
    for i in range(NPT):
        col="red" if w[i]>0 else "royalblue"
        pl.add_mesh(pv.Sphere(radius=0.0024,center=xyz[i]),color=col)
        tag=" HEATER" if i==o["heat_node"] else ""
        pl.add_point_labels([xyz[i]+np.array([0,0,0.009 if i%2 else -0.009])],[f"P{i} W={w[i]:+.2f}{tag}"],
                            font_size=15,text_color="darkred" if w[i]>0 else "navy",shape=None,always_visible=True)
    for P,t,c in [(A,"A (+X top)","red"),(O,"O (-X top)","blue")]:
        pl.add_mesh(pv.Cube(center=P,x_length=0.004,y_length=0.004,z_length=0.004),color=c)
        pl.add_point_labels([P+np.array([0,0,0.012])],[t],font_size=16,text_color=c,shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.22),(0,0,0.05),(0,0,1)]; pl.set_background("white"); pl.camera.zoom(1.2)
    tmp=os.path.join(tempfile.mkdtemp(),"p.png"); pl.screenshot(tmp); pl.close()
    img=plt.imread(tmp); mk=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(mk)
    img=img[max(0,ys.min()-12):ys.max()+12,max(0,xs.min()-12):xs.max()+12]; h,wd=img.shape[:2]
    fig,ax=plt.subplots(figsize=(10,10*h/wd)); ax.imshow(img); ax.axis("off")
    ax.set_title("温度の候補点 P0〜P4 と反りの評価点 A・O\n赤＝W>0（上がると反りが増える・A側）、青＝W<0（O側）",fontsize=13)
    fig.tight_layout(pad=0.2); fig.savefig(os.path.join(IMG,"sensor_pair_positions.png"),dpi=150,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("[pair] figures written")


def main():
    o=json.load(open(OUT)) if (os.path.exists(OUT) and "--recompute" not in sys.argv) else compute()
    figures(o)


if __name__=="__main__": main()
