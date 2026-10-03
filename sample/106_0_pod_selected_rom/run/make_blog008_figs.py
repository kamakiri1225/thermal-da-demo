"""blog_008 用：観測点の位置（3D）と、観測していない A/O の変位の時刻歴アニメーション.

1) blog008_obs_points.png  … 温度2点・変位2点（構成ごと）・評価点A/Oを3Dで示す
2) blog008_AO_anim.gif     … 観測していない Uz(A)・Uz(O)・反りA−O が同化で合っていく様子

数値は run/disp_elsewhere_predict_AO.py と同じ設定を再現して計算する。
再現: OMP_NUM_THREADS=4 python3 run/make_blog008_figs.py
"""
from __future__ import annotations
import os, sys, json, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TNODES=[2,0]; COMP=["Ux","Uy","Uz"]
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005)
SEL=[((0.0097,-0.0362,0.1005),0),((-0.0346,0.0144,0.1005),2)]   # 選定2点（Ux, Uz）
BC =[((0.028,0.0,0.075),2),((-0.028,0.0,0.075),2)]              # B・C


def setup():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    m=dict(C=d["C"],Km=rg.tri_to_matrix(d["K_upper"],NPT),h=float(d["h"]),heat=int(d["heat_node"]))
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    m.update(U=U,mean=mean,pod=pod,UPp=UPp,um=op["u_mean"],D=op["D"],coords=op["coords"],
             Wall=np.einsum("nck,kj->ncj",op["D"],UPp),xyzP=kv["xyz"])
    return m


def near(m,x): return int(np.linalg.norm(m["coords"]-np.array(x),axis=1).argmin())


def run_da(m,obs,seed,Ttr,cyc,iA,iO):
    ws=[m["Wall"][i,c] for i,c in obs]; u0=[m["um"][i,c] for i,c in obs]
    rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
    Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.diag([SIG_T**2]*len(TNODES)+[SIG_U**2]*len(ws))
    wA=m["Wall"][iA,2]; wO=m["Wall"][iO,2]; aA=m["um"][iA,2]; aO=m["um"][iO,2]
    mp=m["mean"][m["pod"]]
    rec=[[aA+wA@(Z[:,:NPT].mean(0)-mp),aO+wO@(Z[:,:NPT].mean(0)-mp)]]; tp=0.0
    for ci,tb in enumerate(cyc,1):
        Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],m["C"],m["Km"],Z[:,IH],Z[:,IQ],m["heat"],tp,tb,DT); tp=tb
        yv=list(Ttr[ci][TNODES]); Yf=Z[:,TNODES]
        for w,b0 in zip(ws,u0):
            yv.append(b0+w@(Ttr[ci]-mp)); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
        y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        mu=Z[:,:NPT].mean(0); rec.append([aA+wA@(mu-mp),aO+wO@(mu-mp)])
    return np.array(rec)


def main():
    m=setup(); coords=m["coords"]
    iA=near(m,A_XYZ); iO=near(m,O_XYZ)
    sel=[(near(m,x),c) for x,c in SEL]; bc=[(near(m,x),c) for x,c in BC]
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,m["C"],m["Km"],m["h"],1.0,m["heat"],a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); mp=m["mean"][m["pod"]]
    truth=np.array([[m["um"][iA,2]+m["Wall"][iA,2]@(t-mp),m["um"][iO,2]+m["Wall"][iO,2]@(t-mp)] for t in Ttr])
    runs={}
    for nm,obs in [("変位なし（温度2点のみ）",[]),("B・C（中段 z=75 mm）",bc),("選定2点（上面）",sel)]:
        runs[nm]=np.mean([run_da(m,obs,s,Ttr,cyc,iA,iO) for s in SEEDS],axis=0)
        print("[fig]",nm,"done",flush=True)

    # ── 1) 3D 観測点図 ──
    import pyvista as pv, vtk, cylinder_mesh
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    H=0.1005
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    nc=np.array([x for _n,x in mesh["nodes"]]); idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),nc)
    pl=pv.Plotter(off_screen=True,window_size=(1300,820)); pl.add_mesh(ug,color="lightsteelblue",opacity=0.32)
    # 温度センサ P2,P0
    for i,lab in [(2,"P2 TEMP"),(0,"P0 TEMP")]:
        p0=m["xyzP"][i]
        pl.add_mesh(pv.Sphere(radius=0.0024,center=p0),color="red")
        pl.add_point_labels([p0+np.array([0,0,-0.010 if i==2 else 0.010])],[lab],font_size=15,text_color="darkred",shape=None,always_visible=True)
    # 評価点 A/O（測らない）
    for i,lab,col in [(iA,"A  EVAL (not measured)","black"),(iO,"O  EVAL (not measured)","dimgray")]:
        pl.add_mesh(pv.Sphere(radius=0.0030,center=coords[i]),color=col)
        pl.add_point_labels([coords[i]+np.array([0,0,0.020])],[lab],font_size=15,text_color=col,shape=None,always_visible=True)
    # 観測した変位点
    for (i,c),lab,col in zip(sel,["SEL-1  Ux","SEL-2  Uz"],["#1f9e4b","#1f9e4b"]):
        pl.add_mesh(pv.Cube(center=coords[i],x_length=0.004,y_length=0.004,z_length=0.004),color=col)
        dirv=[(1,0,0),(0,1,0),(0,0,1)][c]
        pl.add_mesh(pv.Arrow(start=coords[i],direction=dirv,scale=0.022),color=col)
        pl.add_point_labels([coords[i]+np.array([0,0,0.009])],[lab],font_size=15,text_color="#13702f",shape=None,always_visible=True)
    for (i,c),lab in zip(bc,["B  Uz","C  Uz"]):
        pl.add_mesh(pv.Cube(center=coords[i],x_length=0.0035,y_length=0.0035,z_length=0.0035),color="#E67E22")
        pl.add_point_labels([coords[i]+np.array([0,0,-0.012])],[lab],font_size=14,text_color="#A35200",shape=None,always_visible=True)
    pl.camera_position=[(0.27,-0.25,0.24),(0,0,0.05),(0,0,1)]; pl.set_background("white"); pl.camera.zoom(1.1)
    tmp=os.path.join(tempfile.mkdtemp(),"p.png"); pl.screenshot(tmp); pl.close()
    img=plt.imread(tmp); mk=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(mk)
    img=img[max(0,ys.min()-12):ys.max()+12,max(0,xs.min()-12):xs.max()+12]; hh,ww=img.shape[:2]
    fig=plt.figure(figsize=(15,15*0.60))
    ax=fig.add_axes([0.0,0.0,0.50,0.90]); ax.imshow(img); ax.axis("off")
    ax.set_title("斜めから見た図",fontsize=12.5)
    bx=fig.add_axes([0.54,0.06,0.44,0.80]); ang=np.linspace(0,2*np.pi,361)
    bx.fill(37.5*np.cos(ang),37.5*np.sin(ang),color="#DCE3EA"); bx.fill(20*np.cos(ang),20*np.sin(ang),color="white")
    bx.plot(37.5*np.cos(ang),37.5*np.sin(ang),color="#7F8C8D"); bx.plot(20*np.cos(ang),20*np.sin(ang),color="#7F8C8D")
    half=0.100/0.0375/2; th=np.linspace(-half,half,100)
    bx.plot(39.5*np.cos(th),39.5*np.sin(th),color="orangered",lw=7,solid_capstyle="butt")
    bx.text(44,-14,"ヒータ",color="orangered",fontsize=10.5)
    for i,lab,off in [(2,"P2 温度センサ",(10,8)),(0,"P0 温度センサ",(-98,-4))]:
        p0=m["xyzP"][i]*1000
        bx.plot(p0[0],p0[1],"o",ms=11,color="red",mec="k")
        bx.annotate(lab,(p0[0],p0[1]),textcoords="offset points",xytext=off,fontsize=10,color="darkred")
    for i,lab,off in [(iA,"評価点 A（測らない）",(-52,22)),(iO,"評価点 O（測らない）",(-66,20))]:
        p0=coords[i]*1000
        bx.plot(p0[0],p0[1],"o",ms=13,color="k")
        bx.annotate(lab,(p0[0],p0[1]),textcoords="offset points",xytext=off,fontsize=10.5,color="k",weight="bold")
    for (i,c),lab,off in zip(sel,["選定1点目\nUx（半径方向）","選定2点目\nUz（上下方向）"],[(-18,-46),(-94,6)]):
        p0=coords[i]*1000
        bx.plot(p0[0],p0[1],"s",ms=12,color="#1f9e4b",mec="k")
        bx.annotate(lab,(p0[0],p0[1]),textcoords="offset points",xytext=off,fontsize=10,color="#13702f",weight="bold")
    # B・C は A・O の真下（z=75 mm）なので、真上から見ると重なる。少し内側にずらして描く
    for (i,c),lab,off,dx in zip(bc,["B（Uz）\nA の真下","C（Uz）\nO の真下"],[(6,-34),(-52,-34)],[-6.5,6.5]):
        p0=coords[i]*1000
        bx.plot(p0[0]+dx,p0[1],"^",ms=10,color="#E67E22",mec="k")
        bx.annotate(lab,(p0[0]+dx,p0[1]),textcoords="offset points",xytext=off,fontsize=9,color="#A35200")
    bx.annotate("",xy=(58,-46),xytext=(45,-46),arrowprops=dict(arrowstyle="->")); bx.text(60,-47,"+x",fontsize=11,va="center")
    bx.set_xlim(-56,80); bx.set_ylim(-60,56); bx.set_aspect("equal"); bx.axis("off")
    bx.set_title("真上から見た図（寸法 mm）\nB・C は中段 z=75 mm、ほかは上面 z=100.5 mm",fontsize=12)
    fig.suptitle("観測点と評価点　赤＝温度センサ／緑＝選定した変位2点／橙＝B・C／黒＝評価点 A・O（一度も測らない）",fontsize=13)
    fig.savefig(os.path.join(IMG,"blog008_obs_points.png"),dpi=150,bbox_inches="tight",pad_inches=0.12); plt.close(fig)
    print("[fig] 3D done")

    # グラフのアニメーションは分かりにくいため廃止（3Dアニメは run/make_blog008_3d_anim.py）


if __name__=="__main__": main()
