"""blog_001用：①解析対象のセットアップ図（3D＋断面寸法） ②t=300s温度場スナップショット.

ヒータ実範囲はheaterSelection.stlのBBoxより: 外周+X側 θ≈±80°, z=24.75〜75.75mm。
出力: docs/img/blog001_setup.png, docs/img/blog001_temp_field300.png
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
from dacore import rom_general as rg
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
import cylinder_mesh, vtk
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results"); KC=273.15
R_IN,R_OUT,H=0.020,0.0375,0.1005
HZ0,HZ1=0.02475,0.07575; HTH=80  # ヒータ: z範囲, ±角度[deg]


def cyl_grid():
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,R_IN,R_OUT,H)
    import pyvista as pv
    coords=np.array([c for _n,c in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cl=[]
    for _e,conn in mesh["elements"]: cl.append(8); cl.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cl),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    return ug,coords


def heater_shell():
    import pyvista as pv
    th=np.radians(np.linspace(-HTH,HTH,40)); z=np.linspace(HZ0,HZ1,14)
    TH,Z=np.meshgrid(th,z)
    r=R_OUT*1.002
    pts=np.column_stack([r*np.cos(TH).ravel(), r*np.sin(TH).ravel(), Z.ravel()])
    g=pv.StructuredGrid(); g.points=pts; g.dimensions=(len(th),len(z),1)
    return g


def trim(imgpath,pad=12):
    img=plt.imread(imgpath); m=np.any(img[...,:3]<0.96,axis=-1)
    ys,xs=np.where(m)
    return img[max(0,ys.min()-pad):ys.max()+pad, max(0,xs.min()-pad):xs.max()+pad]


def setup_fig():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    ug,_=cyl_grid()
    pl=pv.Plotter(off_screen=True,window_size=(880,940))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.42,show_edges=False)
    pl.add_mesh(heater_shell(),color="red",opacity=0.95)
    pl.add_point_labels([[R_OUT+0.004,0,0.052]],["heater 15 W (0-300 s)"],font_size=16,
                        text_color="red",shape=None,always_visible=True)
    pl.add_point_labels([[0,-R_OUT-0.002,0.098]],["air (natural convection)"],font_size=14,
                        text_color="dimgray",shape=None,always_visible=True)
    pl.add_point_labels([[-R_OUT-0.002,0,0.002]],["bottom fixed (FEM)"],font_size=14,
                        text_color="#33475b",shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.15)
    p3=os.path.join(tempfile.mkdtemp(),"s.png"); pl.screenshot(p3); pl.close()
    img3=trim(p3)

    fig=plt.figure(figsize=(13.6,6.6))
    axA=fig.add_axes([0.015,0.03,0.46,0.86]); axA.imshow(img3); axA.axis("off")
    axA.set_title("① 3D：中空円筒＋ヒータ（赤）",fontsize=13,weight="bold")
    # --- 断面図 ---
    axB=fig.add_axes([0.52,0.10,0.30,0.80])
    for sgn in (+1,-1):
        axB.add_patch(Rectangle((sgn*R_IN*1000 if sgn>0 else -R_OUT*1000,0),
                     (R_OUT-R_IN)*1000,H*1000,fc="#c9d6e4",ec="#33475b",lw=1.5))
    axB.add_patch(Rectangle((R_OUT*1000,HZ0*1000),3.5,(HZ1-HZ0)*1000,fc="red",ec="k"))
    axB.text(R_OUT*1000+5,52,"ヒータ\n15 W\n(0-300 s)",color="red",fontsize=11,weight="bold",va="center")
    axB.add_patch(Rectangle((-R_OUT*1000,-4),2*R_OUT*1000,4,fc="none",ec="#33475b",hatch="///"))
    axB.text(0,-9,"底面固定（FEM）",ha="center",fontsize=10.5,color="#33475b")
    axB.annotate("",xy=(R_IN*1000,108),xytext=(0,108),arrowprops=dict(arrowstyle="<->"))
    axB.text(R_IN*500,111,"r_in=20",ha="center",fontsize=10)
    axB.annotate("",xy=(R_OUT*1000,118),xytext=(0,118),arrowprops=dict(arrowstyle="<->"))
    axB.text(R_OUT*500,121,"r_out=37.5",ha="center",fontsize=10)
    axB.annotate("",xy=(-46,0),xytext=(-46,H*1000),arrowprops=dict(arrowstyle="<->"))
    axB.text(-52,H*500,"H=100.5",rotation=90,va="center",fontsize=10)
    axB.text(0,55,"中空\n（空気）",ha="center",fontsize=10,color="dimgray")
    axB.text(-30,90,"周囲：空気\n（自然対流）",fontsize=9.5,color="dimgray")
    axB.set_xlim(-60,62); axB.set_ylim(-14,128); axB.set_aspect("equal"); axB.axis("off")
    axB.set_title("② 断面と寸法 [mm]",fontsize=13,weight="bold")
    # --- 条件表＋加熱プロファイル ---
    axC=fig.add_axes([0.845,0.52,0.14,0.36]); axC.axis("off")
    axC.text(0,1.0,"材料：鋼",fontsize=11.5,weight="bold",va="top")
    axC.text(0,0.86,"E = 205 GPa\nν = 0.3\nα = 1.2e-5 /K\nk = 50 W/(m·K)\nT_ref = 20 ℃",
             fontsize=10.5,va="top",linespacing=1.6)
    axD=fig.add_axes([0.85,0.12,0.135,0.30])
    t=[0,300,300,600]; q=[15,15,0,0]
    axD.plot(t,q,color="red",lw=2.5); axD.fill_between([0,300],[15,15],color="red",alpha=.15)
    axD.set_xlabel("t [s]",fontsize=9); axD.set_ylabel("Q [W]",fontsize=9)
    axD.set_xticks([0,300,600]); axD.set_yticks([0,15]); axD.tick_params(labelsize=8.5)
    axD.set_title("加熱プロファイル",fontsize=10.5)
    fig.suptitle("解析対象：片側ヒータ付き中空円筒（鋼）― OpenFOAM(CHT)＋FrontISTR(熱弾性)",
                 fontsize=14,weight="bold",y=0.985)
    out=os.path.join(IMG,"blog001_setup.png")
    fig.savefig(out,dpi=140,bbox_inches="tight"); plt.close(fig); print("wrote",out)


def field300_fig():
    import pyvista as pv
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; K=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])
    _,tr=rg.integrate_single(np.full(5,rg.T_AIR_K),C,K,h,1.0,heat,0,300,2.0)
    T5=tr[-1]
    a=np.linalg.pinv(U[pod,:])@(T5-mean[pod])
    field=mean+U@a-KC
    ug,coords=cyl_grid()
    _,near=cKDTree(Cc).query(coords)
    ug.point_data["T"]=field[near]
    pl=pv.Plotter(off_screen=True,window_size=(880,940))
    pl.add_mesh(ug,scalars="T",cmap="turbo",n_colors=16,
                scalar_bar_args={"title":"T [degC]","fmt":"%.1f"})
    pl.add_mesh(heater_shell(),color="white",opacity=0.25)
    pl.add_point_labels([[R_OUT+0.004,0,0.052]],["heater side"],font_size=15,
                        text_color="black",shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.15)
    p3=os.path.join(tempfile.mkdtemp(),"f.png"); pl.screenshot(p3); pl.close()
    img=trim(p3); hh,ww=img.shape[:2]; fw=7.4
    fig,ax=plt.subplots(figsize=(fw,fw*hh/ww)); ax.imshow(img); ax.axis("off")
    ax.set_title("計算結果の例：t=300 s（加熱終了直前）の固体温度場\n"
                 "ヒータ側（+X）が最も熱く、反対側・底面へ向かって温度が下がる",fontsize=12.5)
    out=os.path.join(IMG,"blog001_temp_field300.png")
    fig.savefig(out,dpi=140,bbox_inches="tight",pad_inches=0.04); plt.close(fig); print("wrote",out)


if __name__=="__main__":
    setup_fig(); field300_fig()
