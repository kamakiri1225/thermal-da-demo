"""実ソルバOIの観測点（hot/cold温度センサ・上面変位2点）の位置を図示する.

温度観測: hot=(32,0,50.25)mm ヒータ側(+X)中央高さ, cold=(-32,0,50.25)mm 反対側(-X)中央高さ
変位観測: 上面(z=100.5mm)のヒータ側(+X)平均Uz と 反対側(-X)平均Uz

出力: docs/img/oi_obs_points.png
再現: OMP_NUM_THREADS=4 python3 run/make_oi_obs_points_fig.py
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
import cylinder_mesh, vtk
IMG=os.path.join(ROOT,"docs","img")
# 観測点 [m]
HOT=(0.032,0,0.05025); COLD=(-0.032,0,0.05025)      # 温度センサ
DH=(0.032,0,0.1005);   DO=(-0.032,0,0.1005)         # 上面変位(ヒータ側/反対側)


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    pl=pv.Plotter(off_screen=True,window_size=(1000,900))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.35,show_edges=False)
    def mark(c,col,txt,dz):
        pl.add_mesh(pv.Sphere(radius=0.0022,center=c),color=col)
        pl.add_point_labels([np.array(c)+np.array([0,0,dz])],[txt],font_size=17,
                            text_color=col,shape=None,always_visible=True)
    mark(HOT ,"red",   "hot  (temp, heater side +X)", 0.010)
    mark(COLD,"blue",  "cold (temp, opposite -X)",   0.010)
    mark(DH  ,"darkorange","Disp A (heater side, top +X)", 0.012)
    mark(DO  ,"purple","Disp O (opposite, top -X)", 0.020)
    pl.add_axes(line_width=3)
    pl.camera_position=[(0.27,-0.24,0.24),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    p=os.path.join(tempfile.mkdtemp(),"o.png"); pl.screenshot(p); pl.close()
    # 白余白を自動トリミング
    img=plt.imread(p); m=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(m); pad=14
    img=img[max(0,ys.min()-pad):ys.max()+pad, max(0,xs.min()-pad):xs.max()+pad]
    h,w=img.shape[:2]
    fig,ax=plt.subplots(figsize=(9,9*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("実ソルバOIの観測点：温度センサ hot/cold（中央高さ）と 上面変位2点\n"
                 "ヒータ側=+X（赤・橙）／反対側=−X（青・紫）。中空円柱、上面自由・底面固定。",fontsize=13)
    fig.savefig(os.path.join(IMG,"oi_obs_points.png"),dpi=130,bbox_inches="tight",pad_inches=0.05); plt.close(fig)
    print("[obs] wrote docs/img/oi_obs_points.png")


if __name__=="__main__": main()
