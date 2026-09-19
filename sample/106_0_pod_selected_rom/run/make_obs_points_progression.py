"""段階比較（温度1点→2点→＋変位2点）で使った観測点の位置図（ParaView風）.

opencae_sensor_progression.py の観測点:
  温度: P2(ヒータ側・高感度) → 2点目に P4(反対側)
  変位: 高W=105のW行感度が高い上端2点 / 低W=底面2点
DA結果を示すときは必ずこの種の観測点位置図を併示する（プロジェクト方針）。

出力: docs/img/obs_points_progression.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_obs_points_progression.py
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
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
KINVH=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","kinvh_sensitivity.npz")


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz")); xyz=d["xyz"]
    kh=np.load(KINVH); hiN=[int(i) for i in kh["hi"]]; loN=[int(i) for i in kh["lo"]]

    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([c for _n,c in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)

    pl=pv.Plotter(off_screen=True,window_size=(1100,820))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.35,show_edges=False)
    def mark(c,col,txt,off):
        pl.add_mesh(pv.Sphere(radius=0.0026,center=np.asarray(c)),color=col)
        pl.add_point_labels([np.asarray(c)+np.asarray(off)],[txt],font_size=16,
                            text_color=col,shape=None,always_visible=True)
    # 温度センサ（ROM点 P2, P4）
    mark(xyz[2],"red",  "T1 = P2 (temp, heater side)",(0.006,0,0.006))
    mark(xyz[4],"blue", "T2 = P4 (temp, opposite)",  (0.004,0,-0.010))
    # 変位: 高W(上端2点) / 低W(底面2点)
    for k,i in enumerate(hiN):
        mark(coords[i],"darkorange",f"disp hi-W S{k+1} (top)",(0.005,0,0.008+0.006*k))
    for k,i in enumerate(loN):
        mark(coords[i],"gray",f"disp lo-W L{k+1} (bottom)",(0.005,0,-0.010-0.006*k))
    pl.add_axes(line_width=3)
    pl.camera_position=[(0.29,-0.26,0.26),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    p=os.path.join(tempfile.mkdtemp(),"o.png"); pl.screenshot(p); pl.close()

    img=plt.imread(p); m=np.any(img[...,:3]<0.96,axis=-1); ys,xs=np.where(m); pad=14
    img=img[max(0,ys.min()-pad):ys.max()+pad, max(0,xs.min()-pad):xs.max()+pad]
    h,w=img.shape[:2]
    fig_w=8.8; fig,ax=plt.subplots(figsize=(fig_w,fig_w*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("段階比較で使った観測点の位置\n"
                 "温度：T1=P2（ヒータ側・高感度）→2点目に T2=P4（反対側）。\n"
                 "変位：高W＝上端2点（橙）／低W＝底面2点（灰）。底面固定・中空円柱。",fontsize=12.5)
    out=os.path.join(IMG,"obs_points_progression.png")
    fig.savefig(out,dpi=130,bbox_inches="tight",pad_inches=0.05); plt.close(fig)
    print("[obs] wrote",out)


if __name__=="__main__": main()
