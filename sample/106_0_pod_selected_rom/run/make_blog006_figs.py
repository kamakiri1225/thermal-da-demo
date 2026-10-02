"""blog_006 用の図：温度2点＋変位2点の位置（3D）と、センサが読む値の時系列.

センサの値は OpenFOAM(CHT, 15 W・0〜300 s) の温度場と、それを入れた FrontISTR の上面 Uz
（results/limit_truth_learned.npz、run/limit_realsolver_truth.py learned で作成）。

出力: docs/img/blog006_sensor_positions.png, docs/img/blog006_sensor_readings.png
再現: python3 run/make_blog006_figs.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
H=0.1005; A=np.array([0.028,0,H]); O=np.array([-0.028,0,H])


def positions():
    import pyvista as pv, vtk, cylinder_mesh
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    xyz=np.load(os.path.join(RES,"qdeim_points.npz"))["xyz"]
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    coords=np.array([x for _n,x in mesh["nodes"]]); idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    pl=pv.Plotter(off_screen=True,window_size=(1200,760)); pl.add_mesh(ug,color="lightsteelblue",opacity=0.35)
    for i in range(5):
        used=i in (0,2)
        pl.add_mesh(pv.Sphere(radius=0.0026 if used else 0.0016,center=xyz[i]),color="red" if used else "gray")
        lab=f"P{i} TEMP SENSOR" if used else f"P{i} (not measured)"
        if i==2: lab+=" / HEATER"
        pl.add_point_labels([xyz[i]+np.array([0,0,0.009 if i in (0,3) else -0.009])],[lab],font_size=15 if used else 12,
                            text_color="darkred" if used else "dimgray",shape=None,always_visible=True)
    for P,t in [(A,"A: DISP GAUGE (+X top)"),(O,"O: DISP GAUGE (-X top)")]:
        pl.add_mesh(pv.Cube(center=P,x_length=0.004,y_length=0.004,z_length=0.004),color="darkorange")
        pl.add_mesh(pv.Arrow(start=P+np.array([0,0,0.003]),direction=(0,0,1),scale=0.02),color="darkorange")
        pl.add_point_labels([P+np.array([0,0,0.026])],[t],font_size=15,text_color="darkorange",shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.22),(0,0,0.05),(0,0,1)]; pl.set_background("white"); pl.camera.zoom(1.15)
    tmp=os.path.join(tempfile.mkdtemp(),"p.png"); pl.screenshot(tmp); pl.close()
    img=plt.imread(tmp); mk=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(mk)
    img=img[max(0,ys.min()-12):ys.max()+12,max(0,xs.min()-12):xs.max()+12]; h,w=img.shape[:2]
    fig,ax=plt.subplots(figsize=(10,10*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("温度2点＋変位2点の位置\n赤＝温度センサ（P2・P0）、橙＝変位計（上面A・O、上下方向 Uz）、灰＝測っていない代表点",fontsize=13)
    fig.tight_layout(pad=0.2); fig.savefig(os.path.join(IMG,"blog006_sensor_positions.png"),dpi=150,bbox_inches="tight",pad_inches=0.04); plt.close(fig)


def readings():
    d=np.load(os.path.join(RES,"limit_truth_learned.npz")); pod=np.load(os.path.join(RES,"qdeim_points.npz"))["cell_idx"]
    t=d["times"]; T=d["Tfield"][:,pod]-273.15; uz=d["uz"]
    fig,axs=plt.subplots(1,2,figsize=(13,4.6))
    ax=axs[0]
    ax.plot(t,T[:,2],"o-",color="#C0392B",label="温度 P2（ヒータ横）")
    ax.plot(t,T[:,0],"s-",color="#E67E22",label="温度 P0")
    ax.axvspan(0,300,color="#FDEBD0",alpha=.5,label="加熱中（15 W）")
    ax.set_xlabel("時刻 [s]"); ax.set_ylabel("温度 [℃]"); ax.set_title("温度センサ2本が読む値",fontsize=13); ax.grid(alpha=.3); ax.legend(fontsize=10)
    ax=axs[1]
    ax.plot(t,uz[:,0],"o-",color="#C0392B",label="変位 A（ヒータ側）")
    ax.plot(t,uz[:,1],"s-",color="#2E6FD8",label="変位 O（反対側）")
    ax.plot(t,uz[:,0]-uz[:,1],"--",color="k",lw=2,label="反り A−O")
    ax.axvspan(0,300,color="#FDEBD0",alpha=.5)
    ax.set_xlabel("時刻 [s]"); ax.set_ylabel("上下方向の変位 Uz [µm]"); ax.set_title("変位計2本が読む値",fontsize=13); ax.grid(alpha=.3); ax.legend(fontsize=10)
    fig.suptitle("センサが読む値（ノイズなしの真値。OpenFOAM の温度場＋FrontISTR の変形）。同化ではここに温度 0.3 K・変位 0.3 µm のノイズを足して使う",fontsize=11.5)
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(IMG,"blog006_sensor_readings.png"),dpi=150); plt.close(fig)


if __name__=="__main__":
    positions(); readings(); print("wrote blog006 figs")
