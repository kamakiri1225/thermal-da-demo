"""ROMの5ノードと結合ネットワークを3D表示（blog_004 §6-1用）.

集中定数ROMの式  C_i dT_i/dt = Σ_j K_ij (T_j-T_i) + q_i - h(T_i-T_air)
が「どのノードの式か」を示す。ノード＝POD+Q-DEIM選定5点(P0..P4)、
P2=ヒータ発熱ノード(q投入)、辺の太さ=コンダクタンス K_ij。

出力: docs/img/rom_network.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_rom_network_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
sys.path.insert(0, ROOT)
from dacore import plots as _p
from dacore import rom_general as rg
import matplotlib.pyplot as plt
import cylinder_mesh
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    xyz=d["xyz"]; heat=int(d["heat_node"]); K=rg.tri_to_matrix(d["K_upper"],5); C=d["C"]
    Kmax=K.max()

    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    mco=np.array([xyz2 for _n,xyz2 in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mco)

    pl=pv.Plotter(off_screen=True,window_size=(880,940))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.35,show_edges=False)

    # 辺：コンダクタンス K_ij を太さで表現（非ゼロのみ）
    for i in range(5):
        for j in range(i+1,5):
            if K[i,j]<=0.3: continue
            r=0.0004+0.0016*(K[i,j]/Kmax)
            tube=pv.Line(xyz[i],xyz[j]).tube(radius=r)
            pl.add_mesh(tube,color="#3b5b78",opacity=0.9)
            mid=(xyz[i]+xyz[j])/2
            pl.add_point_labels([mid],[f"K={K[i,j]:.1f}"],font_size=13,
                text_color="#1b3550",shape=None,always_visible=True)

    # ノード球：ヒータノードはオレンジ、他はクリムゾン
    for j,p in enumerate(xyz):
        col="#e8820c" if j==heat else "#c0392b"
        pl.add_mesh(pv.Sphere(radius=0.0030,center=p),color=col)
        tag=f"P{j} (heater, q)" if j==heat else f"P{j}"
        pl.add_point_labels([p+np.array([0.005,0,0.006*(1 if j%2 else -1)])],
            [tag],font_size=17,text_color=col,shape=None,always_visible=True)

    pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.2)
    p4=os.path.join(IMG,"_tmp_net.png"); pl.screenshot(p4); pl.close()

    img=plt.imread(p4)
    mask=np.any(img[...,:3]<0.96,axis=-1)
    ys,xs=np.where(mask); pad=12
    y0,y1=max(0,ys.min()-pad),min(img.shape[0],ys.max()+pad)
    x0,x1=max(0,xs.min()-pad),min(img.shape[1],xs.max()+pad)
    img=img[y0:y1,x0:x1]; h,w=img.shape[:2]
    fig_w=8.6; fig_h=fig_w*(h/w)          # 真のアスペクト比（つぶれ防止：aspect autoを使わない）
    fig,ax=plt.subplots(figsize=(fig_w,fig_h)); ax.imshow(img); ax.axis("off")
    ax.set_title("ROMのノード＝POD＋Q-DEIM選定5点。各ノードで集中定数の式を解く\n"
                 "辺の太さ＝コンダクタンス K_ij（P0-P2が最強）。オレンジ＝ヒータ発熱ノードP2にqを投入",
                 fontsize=12.5)
    fig.savefig(os.path.join(IMG,"rom_network.png"),dpi=130,bbox_inches="tight",pad_inches=0.04)
    plt.close(fig); os.remove(p4)
    print("[rom-net] wrote docs/img/rom_network.png")


if __name__=="__main__": main()
