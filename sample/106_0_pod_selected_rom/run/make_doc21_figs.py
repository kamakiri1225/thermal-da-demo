"""docs/21_advantage_verification.md 用の図.

1) doc21_conditions.png   … 3つの真値条件のヒータ発熱と、真値の反り A−O の時刻歴
2) doc21_positions.png    … 温度センサ（P2, P0, P4）と変位計（A, O）の位置（3D）
3) doc21_timeseries.png   … 反り A−O の時刻歴：真値／回帰式／同化（温度2点）／同化（温度2点＋変位2点）（温度センサ P2+P4）
4) doc21_Q.png            … 推定した発熱量 Q の時刻歴（3条件、温度P2+P4＋変位A/O）
5) doc21_q_contribution.png … ③ 発熱量を真値へ近づけた量の内訳（温度センサ／変位計、5 seed）

同化の設定は run/limit_realsolver_truth.py と同じ（5 seed 平均）。
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_doc21_figs.py
"""
from __future__ import annotations
import os, sys, json, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02; TAIR=293.15
SEEDS=[20260913,20260914,20260915,20260916,20260917]
CASES=[("learned","15 W・0〜300 s（ROM の係数を決めた計算）",[(0,15),(300,15),(300,0),(600,0)]),
       ("q25","25 W・0〜300 s（新しく計算した条件）",[(0,25),(300,25),(300,0),(600,0)]),
       ("intermittent","15 W 間欠加熱（新しく計算した条件）",[(0,15),(150,15),(150,0),(300,0),(300,15),(450,15),(450,0),(600,0)])]
CK={"truth":"k","reg":"#7F8C8D","daT":"#2E8B57","daTU":"#C0392B"}


def model():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    kv=np.load(os.path.join(RES,"qdeim_points.npz")); U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); op=np.load(os.path.join(RES,"disp_operator.npz"))
    return dict(C=d["C"],Km=rg.tri_to_matrix(d["K_upper"],NPT),heat=int(d["heat_node"]),pod=pod,mean=mean,UPp=UPp,ub=op["uz_mean"],D=op["Dmode"])


def truth(case):
    d=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); pod=np.load(os.path.join(RES,"qdeim_points.npz"))["cell_idx"]
    return d["times"],d["Tfield"][:,pod],d["uz"]


def da(m,t,T5,uz,nodes,use_disp,seed):
    disp=lambda X: m["ub"]+((X-m["mean"][m["pod"]])@m["UPp"].T)@m["D"].T
    rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
    Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
    u0=disp(Z[:,:NPT].mean(0)); ao=[u0[0]-u0[1]]; Q=[15*Z[:,IQ].mean()]; tp=0.0
    for ci,tb in enumerate(t[1:],1):
        Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],m["C"],m["Km"],Z[:,IH],Z[:,IQ],m["heat"],tp,tb,DT); tp=tb
        yv=list(T5[ci][nodes]); Yf=Z[:,nodes]
        if use_disp: yv+=list(uz[ci]); Yf=np.column_stack([Yf,disp(Z[:,:NPT])])
        y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        u=disp(Z[:,:NPT].mean(0)); ao.append(u[0]-u[1]); Q.append(15*Z[:,IQ].mean())
    return np.array(ao),np.array(Q)


def main():
    m=model(); data={}
    tL,TL,uL=truth("learned"); nodes=[2,4]
    X0=np.column_stack([np.ones(len(tL)),TL[:,nodes]-TAIR]); c,*_=np.linalg.lstsq(X0,uL[:,0]-uL[:,1],rcond=None)
    for case,_,_ in CASES:
        t,T5,uz=truth(case)
        reg=np.mean([np.column_stack([np.ones(len(t)),T5[:,nodes]+np.random.default_rng(s+7).normal(0,SIG_T,(len(t),2))-TAIR])@c for s in SEEDS],axis=0)
        rT=[da(m,t,T5,uz,nodes,False,s) for s in SEEDS]; rU=[da(m,t,T5,uz,nodes,True,s) for s in SEEDS]
        data[case]=dict(t=t,truth=uz[:,0]-uz[:,1],reg=reg,daT=np.mean([r[0] for r in rT],0),daTU=np.mean([r[0] for r in rU],0),
                        QT=np.mean([r[1] for r in rT],0),QTU=np.mean([r[1] for r in rU],0)); print("[doc21]",case,flush=True)
    # 1) 条件
    fig,axs=plt.subplots(2,3,figsize=(15,6.4),sharex=True)
    for j,(case,lab,sch) in enumerate(CASES):
        x,y=zip(*sch); axs[0,j].plot(x,y,color="#E67E22",lw=3); axs[0,j].fill_between(x,y,color="#FDEBD0")
        axs[0,j].set_ylim(0,28); axs[0,j].set_title(lab,fontsize=12); axs[0,j].set_ylabel("ヒータ発熱 [W]"); axs[0,j].grid(alpha=.3)
        axs[0,j].axvspan(0,300,color="none",ec="#888",ls=":",lw=1)
        dd=data[case]; axs[1,j].plot(dd["t"],dd["truth"],"o-",color="k",ms=4)
        axs[1,j].set_ylim(-0.3,5); axs[1,j].set_ylabel("真値の反り A−O [µm]"); axs[1,j].set_xlabel("時刻 [s]"); axs[1,j].grid(alpha=.3)
    axs[0,0].text(5,25.5,"同化モデル(ROM)の仮定：\n0〜300 s に一定の発熱",fontsize=9.5,color="#555",va="top")
    fig.suptitle("3つの真値の条件（上：ヒータの発熱、下：OpenFOAM＋FrontISTR で計算した反り）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig(os.path.join(IMG,"doc21_conditions.png"),dpi=150); plt.close(fig)
    # 3) 反りの時刻歴
    fig,axs=plt.subplots(1,3,figsize=(15.5,4.6),sharey=True)
    for ax,(case,lab,_) in zip(axs,CASES):
        dd=data[case]; t=dd["t"]
        ax.plot(t,dd["truth"],"o-",color=CK["truth"],lw=2.2,ms=4,label="真値")
        ax.plot(t,dd["reg"],"s--",color=CK["reg"],lw=1.8,ms=4,label="回帰式（温度2点）")
        ax.plot(t[1:],dd["daT"][1:],"--",color=CK["daT"],lw=2,label="同化（温度2点）")
        ax.plot(t[1:],dd["daTU"][1:],"-",color=CK["daTU"],lw=2,label="同化（温度2点＋変位2点）")
        ax.set_title(lab,fontsize=12); ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3)
    axs[0].set_ylabel("反り A−O [µm]"); axs[0].legend(fontsize=9.5,loc="upper right")
    fig.suptitle("反り A−O の時刻歴（温度センサ P2＋P4、変位計 A・O、5 seed 平均。同化の最初の予測は 30 s から表示）",fontsize=12.5)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"doc21_timeseries.png"),dpi=150); plt.close(fig)
    # 4) Q
    fig,axs=plt.subplots(1,3,figsize=(15.5,4.3),sharey=True)
    for ax,(case,lab,sch) in zip(axs,CASES):
        dd=data[case]; truthQ={"learned":15,"q25":25,"intermittent":None}[case]
        ax.plot(dd["t"],dd["QTU"],"o-",color=CK["daTU"],ms=4,label="推定（温度2点＋変位2点）")
        ax.plot(dd["t"],dd["QT"],"s--",color=CK["daT"],ms=3,label="推定（温度2点）")
        if truthQ: ax.axhline(truthQ,color="k",lw=1.5,ls=":",label=f"真値 {truthQ} W")
        else: ax.text(300,24,"真値は一定ではない\n（15 W と 0 W の繰り返し）",ha="center",fontsize=10)
        ax.set_title(lab,fontsize=12); ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3); ax.set_ylim(0,30)
    axs[0].set_ylabel("推定した発熱量 Q [W]"); axs[0].legend(fontsize=9.5,loc="lower right")
    fig.suptitle("推定した発熱量の推移（温度センサ P2＋P4、5 seed 平均）",fontsize=12.5)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"doc21_Q.png"),dpi=150); plt.close(fig)
    # 5) ③
    q=json.load(open(os.path.join(RES,"q_contribution_seeds.json")))["seeds"]
    fig,ax=plt.subplots(figsize=(9,4.4)); x=np.arange(len(q)); w=0.38
    ax.bar(x-w/2,[r["toward_truth_temp_W"] for r in q],w,color="#E67E22",label="温度センサの分")
    ax.bar(x+w/2,[r["toward_truth_disp_W"] for r in q],w,color="#2E6FD8",label="変位計の分")
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(x); ax.set_xticklabels([str(r["seed"]) for r in q])
    ax.set_xlabel("seed"); ax.set_ylabel("真値 15 W へ近づけた量 [W]\n（＋＝近づけた、−＝遠ざけた）"); ax.grid(axis="y",alpha=.3); ax.legend()
    ax.set_title("③ 発熱量を真値へ近づけたのはどちらか（20 回の補正の合計）",fontsize=12.5)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"doc21_q_contribution.png"),dpi=150); plt.close(fig)
    # 2) 3D 位置
    import pyvista as pv, vtk, cylinder_mesh
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    xyz=np.load(os.path.join(RES,"qdeim_points.npz"))["xyz"]; H=0.1005
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    coords=np.array([p for _n,p in mesh["nodes"]]); idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    pl=pv.Plotter(off_screen=True,window_size=(1200,760)); pl.add_mesh(ug,color="lightsteelblue",opacity=0.35)
    lab={0:"P0 TEMP (A side)",2:"P2 TEMP / HEATER",4:"P4 TEMP (O side)"}
    for i in range(5):
        if i in lab:
            pl.add_mesh(pv.Sphere(radius=0.0026,center=xyz[i]),color="red" if i!=4 else "royalblue")
            pl.add_point_labels([xyz[i]+np.array([0,0,0.009 if i==0 else -0.009])],[lab[i]],font_size=15,text_color="darkred" if i!=4 else "navy",shape=None,always_visible=True)
        else: pl.add_mesh(pv.Sphere(radius=0.0015,center=xyz[i]),color="gray")
    for P,tx in [(np.array([0.028,0,H]),"A: DISP (+X top)"),(np.array([-0.028,0,H]),"O: DISP (-X top)")]:
        pl.add_mesh(pv.Cube(center=P,x_length=0.004,y_length=0.004,z_length=0.004),color="darkorange")
        pl.add_point_labels([P+np.array([0,0,0.014])],[tx],font_size=15,text_color="darkorange",shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.22),(0,0,0.05),(0,0,1)]; pl.set_background("white"); pl.camera.zoom(1.15)
    tmp=os.path.join(tempfile.mkdtemp(),"p.png"); pl.screenshot(tmp); pl.close()
    img=plt.imread(tmp); mk=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(mk)
    img=img[max(0,ys.min()-12):ys.max()+12,max(0,xs.min()-12):xs.max()+12]; h,w=img.shape[:2]
    fig,ax=plt.subplots(figsize=(9,9*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("使ったセンサの位置\n温度センサ：P2＋P0（2本ともA側）または P2＋P4（A側とO側）／変位計：上面 A・O",fontsize=12.5)
    fig.tight_layout(pad=0.2); fig.savefig(os.path.join(IMG,"doc21_positions.png"),dpi=150,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("[doc21] figures written")


if __name__=="__main__": main()
