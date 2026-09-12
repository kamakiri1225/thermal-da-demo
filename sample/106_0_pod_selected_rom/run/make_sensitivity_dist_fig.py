"""熱感度 dT/dQ の空間分布（どこが高感度/低感度か）を可視化する.

dT/dQ は「その場所の温度が発熱Qにどれだけ反応するか」。5つのROM代表点で有限差分で
求めた dT/dQ を、gappy-POD（PODモード）で全セルへ復元し、分布として描く。
高感度＝発熱に敏感（温度センサ向き）、低感度＝反応が鈍い。

出力: docs/img/sensitivity_dist.png
再現: OMP_NUM_THREADS=4 python3 run/make_sensitivity_dist_fig.py
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
from dacore import rom_general as rg
import cylinder_mesh, vtk
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img"); NPT=5; DT=2.0


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    xyz=d["xyz"]; C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT)
    heat_node=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); Cc=kv["cell_centres"]; pod_cells=kv["cell_idx"]
    # 5点の dT/dQ（有限差分, t=300s）
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens5=(pert-base)/0.1
    # gappy-POD で全セルの dT/dQ 分布へ（モード空間へ射影）
    a=np.linalg.lstsq(U[pod_cells,:], sens5, rcond=None)[0]
    field=U@a
    print(f"[sd] 5点 dT/dQ={np.round(sens5,2)}  復元場 範囲 {field.min():.2f}..{field.max():.2f} K")

    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz2 for _n,xyz2 in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    _,near=cKDTree(Cc).query(coords)
    g=ug.copy(); g.point_data["s"]=field[near]
    # 表と裏の2アングル
    TMP=tempfile.mkdtemp(); shots=[]
    for cam,lab in [([(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)],"ヒータ側(+X)から"),
                    ([(-0.24,0.22,0.20),(0,0,0.05),(0,0,1)],"反対側(−X)から")]:
        pl=pv.Plotter(off_screen=True,window_size=(620,660))
        pl.add_mesh(g,scalars="s",cmap="turbo",n_colors=12,
                    scalar_bar_args={"title":"dT/dQ [K/(発熱倍率)]","title_font_size":16,"label_font_size":13})
        for i in range(NPT):
            pl.add_mesh(pv.Sphere(radius=0.004,center=xyz[i]),color="white")
        pl.camera_position=cam; pl.set_background("white")
        p=os.path.join(TMP,lab[:3]+".png"); pl.screenshot(p); pl.close(); shots.append((p,lab))
    fig,axes=plt.subplots(1,2,figsize=(13,6.4))
    for ax,(p,lab) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(lab,fontsize=13)
    fig.suptitle("熱感度 dT/dQ の分布（発熱Qへの温度の反応。白球=ROM5点）\n"
                 "赤=高感度(発熱に敏感→温度センサ向き) / 青=低感度(反応が鈍い→避ける)",fontsize=14,weight="bold")
    fig.text(0.5,0.02,"※ 5点で有限差分したdT/dQをgappy-POD(PODモード)で全場に復元した分布。"
             "ヒータ側(+X)の柱が高感度、反対側・底面(P4付近)が低感度。",ha="center",fontsize=10,color="dimgray")
    fig.tight_layout(rect=[0,0.04,1,0.93]); fig.savefig(os.path.join(IMG,"sensitivity_dist.png"),dpi=130); plt.close(fig)
    print("[sd] wrote docs/img/sensitivity_dist.png")


if __name__=="__main__": main()
