"""qdeim_points.png をコンパクト（横長）に作り直す.

select_points_qdeim.py はOpenFOAMスナップショット読み込みが要り重いので、
保存済み results/qdeim_points.npz（選点座標・エネルギー）とメッシュだけで3D図を再生成する。
縦長で「でかすぎる」問題を、俯瞰気味カメラ＋横長トリム上限で解消。

出力: docs/img/qdeim_points.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/remake_qdeim_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")


def main():
    d=np.load(os.path.join(RES,"qdeim_points.npz"))
    coords=d["xyz"]; energy=d["energy"]; cum=np.cumsum(energy)

    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    mco=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mco)
    # 円柱の縦横比が自然に見える低めのカメラ（俯瞰しすぎると寸詰まりに見える）
    pl=pv.Plotter(off_screen=True,window_size=(880,940))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.5,show_edges=False)
    for j,xyz in enumerate(coords):
        pl.add_mesh(pv.Sphere(radius=0.0022,center=xyz),color="crimson")
        p=np.round(xyz*1000,0).astype(int)
        pl.add_point_labels([xyz+np.array([0.004,0,0.004*(1 if j%2 else -1)])],
            [f"P{j} ({p[0]},{p[1]},{p[2]})mm"],font_size=15,text_color="crimson",
            shape=None,always_visible=True)
    # 低めのカメラで円柱の高さ/直径比が自然に見えるように
    pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.2)
    p4=os.path.join(IMG,"_tmp_qdeim.png"); pl.screenshot(p4); pl.close()

    img=plt.imread(p4)
    mask=np.any(img[...,:3]<0.96,axis=-1)
    ys,xs=np.where(mask); pad=12
    y0,y1=max(0,ys.min()-pad),min(img.shape[0],ys.max()+pad)
    x0,x1=max(0,xs.min()-pad),min(img.shape[1],xs.max()+pad)
    img=img[y0:y1,x0:x1]
    h,w=img.shape[:2]
    # 真のアスペクト比で描画（aspect autoやキャップは歪みの原因なので使わない）
    fig_w=7.6
    fig,ax=plt.subplots(figsize=(fig_w,fig_w*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("POD＋Q-DEIM で選んだROM代表点（勘でなくデータから選定）\n"
                 f"温度場は2モードで{cum[1]*100:.0f}%説明→少数点で表せる",fontsize=12)
    fig.savefig(os.path.join(IMG,"qdeim_points.png"),dpi=130,bbox_inches="tight",pad_inches=0.04)
    plt.close(fig); os.remove(p4)
    print("[qdeim] rewrote docs/img/qdeim_points.png (compact)")


if __name__=="__main__": main()
