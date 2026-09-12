"""温度センサ(熱感度dT/dQで選定)と変位観測点(上面A/O)の位置を可視化する.

- 5つのROM代表点を dT/dQ（発熱への感度）で色付け
- 高感度ノード=温度センサ候補、低感度ノードも明示
- 変位観測点 A(+X上面), O(-X上面) を別色で表示

dT/dQ の計算: 各ノード温度の「発熱倍率 q_scale への感度」を有限差分で求める。
   dT_i/dq ≈ [ T_i(q=1.1) − T_i(q=1.0) ] / 0.1   （加熱ピーク t=300s で評価）

出力: docs/img/sensor_positions.png
再現: OMP_NUM_THREADS=4 python3 run/make_position_fig.py
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
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; DT=2.0; H=0.1005
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    xyz=d["xyz"]; C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT)
    heat_node=int(d["heat_node"]); h_true=float(d["h"])
    # dT/dQ を有限差分で（加熱ピーク t=300s）
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens=(pert-base)/0.1
    hi=int(np.argmax(sens)); lo=int(np.argmin(sens))
    print("[pos] dT/dQ =",np.round(sens,3),f" 高感度=P{hi} 低感度=P{lo}")

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
    pl=pv.Plotter(off_screen=True,window_size=(1000,880))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.4,show_edges=False)
    # 5ノードを dT/dQ で色付け（球＋数値ラベル）
    smin,smax=sens.min(),sens.max()
    for i in range(NPT):
        frac=(sens[i]-smin)/(smax-smin+1e-9)
        col=plt.cm.viridis(frac)[:3]
        pl.add_mesh(pv.Sphere(radius=0.006,center=xyz[i]),color=col)
        role=" <- HIGH (temp sensor)" if i==hi else (" <- LOW (avoid)" if i==lo else "")
        dz=0.010 if i%2 else -0.010
        pl.add_point_labels([xyz[i]+np.array([0,0,dz])],
            [f"P{i} dT/dQ={sens[i]:.1f}{role}"],font_size=15,text_color="black",
            shape=None,always_visible=True)
    # 変位観測点 A/O（赤・青の別マーカー、ラベルはASCII）
    pl.add_mesh(pv.Sphere(radius=0.006,center=A_XYZ),color="red")
    pl.add_mesh(pv.Sphere(radius=0.006,center=O_XYZ),color="blue")
    pl.add_point_labels([A_XYZ+np.array([0,0,0.008])],["Disp A (+X top)"],font_size=16,text_color="red",shape=None,always_visible=True)
    pl.add_point_labels([O_XYZ+np.array([0,0,0.018])],["Disp O (-X top)"],font_size=16,text_color="blue",shape=None,always_visible=True)
    pl.camera_position=[(0.26,-0.24,0.22),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    p4=os.path.join(tempfile.mkdtemp(),"pos.png"); pl.screenshot(p4); pl.close()
    fig,ax=plt.subplots(figsize=(10,9)); ax.imshow(plt.imread(p4)); ax.axis("off")
    ax.set_title("観測点の位置: 温度センサ(ROM代表点をdT/dQで色付け)と変位観測点(上面A/O)\n"
                 f"温度センサは高感度ノードP{hi}(dT/dQ={sens[hi]:.1f})、避けるのは低感度P{lo}(dT/dQ={sens[lo]:.1f})",
                 fontsize=13)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"sensor_positions.png"),dpi=130); plt.close(fig)
    print("[pos] wrote docs/img/sensor_positions.png")


if __name__=="__main__": main()
