"""未観測の温度3点・変位3点を並べて比較（凡例=上、説明=下）.

観測: 構成に応じ 温度P2 / P2+P0 / +変位A,O(上面)。ここで見るのはすべて未観測点:
  温度: gappy-POD復元で 3セル (mid+Y中央高さ / ヒータ側上端 / 反対側底面)
  変位: C,D (z=75mm, dispop_unobs_z75) と L1 (底面, dispop_loW)
出力: docs/img/unobs_multi_timeseries.png, docs/img/unobs_points_locations.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_unobs_multi_fig.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KC=273.15; NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])
    do=np.load(os.path.join(RES,"disp_operator.npz"))       # 観測に使うA/O
    z75=np.load(os.path.join(RES,"dispop_unobs_z75.npz"))   # 未観測C/D
    lo =np.load(os.path.join(RES,"dispop_loW.npz"))         # 未観測 底面
    def a_of(T5): return (T5-mean[pod])@UPp.T
    def uAO(T5): a=a_of(T5); return do["uz_mean"]+a@do["Dmode"].T
    # 未観測温度3セル（観測P2/P0から離す）
    probes={"T@mid +Y (0,30,50)mm":[0,0.030,0.050],
            "T@top heater (25,0,95)mm":[0.028,0,0.095],
            "T@bottom -X (-28,0,8)mm":[-0.028,0,0.008]}
    pcell=[int(np.linalg.norm(Cc-np.array(p),axis=1).argmin()) for p in probes.values()]
    def temp_at(T5,c): return mean[c]+U[c,:]@a_of(T5)
    # 未観測変位3点: C(z75,+X), D(z75,-X), L1(底面)
    def uz75(T5): a=a_of(T5); return z75["uz_mean"]+a@z75["Dmode"].T
    def ulo(T5):  a=a_of(T5); return lo["uz_mean"]+a@lo["D"].T

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]; ht=(tg>0)&(tg<=300)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(tsens, use_disp, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        recT=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        nobs=len(tsens)+(2 if use_disp else 0)
        Rd=np.diag([SIG_T**2]*len(tsens)+([SIG_U**2]*2 if use_disp else []))
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if nobs:
                yv=list(Ttr[ci][tsens]); Yf=Z[:,tsens]
                if use_disp:
                    yv+=list(uAO(Ttr[ci])); Yf=np.column_stack([Yf,uAO(Z[:,:NPT])])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            recT.append(Z[:,:NPT].mean(0).copy())
        return np.array(recT)

    cfgs=[("同化なし",[],False,"0.45"),
          ("温度1点:P2",[2],False,"tab:blue"),
          ("温度2点:P2+P0",[2,0],False,"tab:green"),
          ("温度2点+変位2点(A,O)",[2,0],True,"tab:red")]
    rec={}
    for name,ts_,ud,col in cfgs:
        rec[name]=np.mean([run(ts_,ud,s) for s in SEEDS],axis=0)
        print("[unobs]",name,"done")

    # --- 2×3 図 ---
    fig,axes=plt.subplots(2,3,figsize=(16.5,9.2))
    Tpan=[(k,pcell[i]) for i,k in enumerate(probes.keys())]
    for j,(lab,c) in enumerate(Tpan):
        ax=axes[0][j]; ax.axvspan(0,300,color="orange",alpha=.07)
        tru=np.array([temp_at(Ttr[k],c) for k in range(len(tg))])-KC
        ax.plot(tg,tru,color="black",lw=4.0,zorder=10)
        for name,_,_,col in cfgs:
            est=np.array([temp_at(rec[name][k],c) for k in range(len(tg))])-KC
            ax.plot(tg,est,lw=1.8,color=col,label=name)
        ax.set_title(f"未観測温度 {lab}",fontsize=11.5,weight="bold"); ax.grid(alpha=.3)
        ax.set_ylabel("温度 [degC]")
    Dpan=[("C: z=75mm ヒータ側",lambda T5: uz75(T5)[...,0]),
          ("D: z=75mm 反対側",  lambda T5: uz75(T5)[...,1]),
          ("L1: 底面(低W)",     lambda T5: ulo(T5)[...,0])]
    for j,(lab,fn) in enumerate(Dpan):
        ax=axes[1][j]; ax.axvspan(0,300,color="orange",alpha=.07)
        tru=np.array([fn(Ttr[k]) for k in range(len(tg))])
        ax.plot(tg,tru,color="black",lw=4.0,zorder=10)
        for name,_,_,col in cfgs:
            est=np.array([fn(rec[name][k]) for k in range(len(tg))])
            ax.plot(tg,est,lw=1.8,color=col)
        ax.set_title(f"未観測変位 {lab}",fontsize=11.5,weight="bold"); ax.grid(alpha=.3)
        ax.set_xlabel("時間 [s]"); ax.set_ylabel("Uz [µm]")
    fig.suptitle("観測に使っていない場所でも合うか ― 未観測の温度3点・変位3点（5seed平均）",
                 fontsize=13.5,weight="bold",y=0.985)
    hd,lb=axes[0][0].get_legend_handles_labels()
    hd=[mlines.Line2D([],[],color="black",lw=4)]+hd; lb=["真値"]+lb
    fig.legend(hd,lb,loc="upper center",bbox_to_anchor=(0.5,0.955),ncol=5,fontsize=12.5,frameon=False)
    fig.text(0.5,0.015,"同化なし(灰)は全点で大外れ。温度＋変位(赤)は未観測の温度・変位とも真値に最も近い"
             "＝少数観測が観測していない場所まで直すのがデータ同化の価値",
             ha="center",fontsize=11.5,weight="bold",color="#222")
    fig.subplots_adjust(left=0.05,right=0.99,top=0.865,bottom=0.115,hspace=0.32,wspace=0.24)
    out=os.path.join(IMG,"unobs_multi_timeseries.png"); fig.savefig(out,dpi=140,bbox_inches="tight"); plt.close(fig)
    print("[unobs] wrote",out)

    # --- 位置図（ParaView風）---
    import cylinder_mesh, vtk, pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([c2 for _n,c2 in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cl=[]
    for _e,conn in mesh["elements"]: cl.append(8); cl.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cl),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    pl=pv.Plotter(off_screen=True,window_size=(880,940))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.32,show_edges=False)
    xyzR=d["xyz"]
    def mark(c,col,txt,off):
        pl.add_mesh(pv.Sphere(radius=0.0026,center=np.asarray(c)),color=col)
        pl.add_point_labels([np.asarray(c)+np.asarray(off)],[txt],font_size=14,
                            text_color=col,shape=None,always_visible=True)
    # 観測（文脈用・小さめ色薄）
    mark(xyzR[2],"red","obs T:P2",(0.005,0,0.005))
    mark(xyzR[0],"salmon","obs T:P0",(0.005,0,-0.009))
    A=[0.028,0,0.1005]; O=[-0.028,0,0.1005]
    mark(A,"orange","obs disp A",(0.004,0,0.007)); mark(O,"orange","obs disp O",(0.004,0,0.012))
    # 未観測温度（緑）
    for (lab,p3),c in zip(probes.items(),pcell):
        short=lab.split(" ")[0]
        mark(Cc[c],"green",f"unobs {short}",(0.005,0,-0.008))
    # 未観測変位（紫）: C/D z75, L1底面
    mark([0.032,0,0.075],"purple","unobs disp C",(0.005,0,0.006))
    mark([-0.032,0,0.075],"purple","unobs disp D",(0.004,0,0.010))
    mark(coords[244],"purple","unobs disp L1",(0.005,0,-0.010))
    pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.2)
    p4=os.path.join(tempfile.mkdtemp(),"u.png"); pl.screenshot(p4); pl.close()
    img=plt.imread(p4); m=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(m); pad=12
    img=img[max(0,ys.min()-pad):ys.max()+pad, max(0,xs.min()-pad):xs.max()+pad]
    h,w=img.shape[:2]; fw=8.4
    fig,ax=plt.subplots(figsize=(fw,fw*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("観測点（赤/橙）と未観測の評価点（緑=温度3点、紫=変位3点）",fontsize=12.5)
    out2=os.path.join(IMG,"unobs_points_locations.png")
    fig.savefig(out2,dpi=130,bbox_inches="tight",pad_inches=0.04); plt.close(fig)
    print("[unobs] wrote",out2)


if __name__=="__main__": main()
