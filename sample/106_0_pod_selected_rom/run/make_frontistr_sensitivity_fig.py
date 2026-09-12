"""FrontISTRの熱感度 W=K⁻¹H の分布を106に持ち込み、変位観測点A/Oを重ねる.

105でDUMPWパッチ版FrontISTRが計算した感度行列 W=K⁻¹H（温度→変位のヤコビアン）を
再利用する。W の「行ノルム」= 各点の変位が温度場にどれだけ敏感か = 変位観測点の情報量。
106の変位観測点 A(+X上面), O(-X上面) を重ねて、なぜ上面が変位観測に良いかを示す。

出力: docs/img/frontistr_sensitivity.png
再現: OMP_NUM_THREADS=4 python3 run/make_frontistr_sensitivity_fig.py
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
IMG=os.path.join(ROOT,"docs","img"); H=0.1005
KINVH=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","kinvh_sensitivity.npz")
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])


def main():
    d=np.load(KINVH)
    row=d["row_sens"]; valid=d["valid"]          # 5040節点(102_1 FEMメッシュ順)
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    print(f"[fs] W行感度 有効{int(valid.sum())}節点 範囲 {row[valid].min():.2f}..{row[valid].max():.2f} µm/K")

    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    g=ug.copy(); g.point_data["row"]=np.where(valid,row,np.nan)
    vmax=float(row[valid].max())
    TMP=tempfile.mkdtemp(); shots=[]
    for cam,lab in [([(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)],"ヒータ側(+X)から"),
                    ([(-0.24,0.22,0.20),(0,0,0.05),(0,0,1)],"反対側(−X)から")]:
        pl=pv.Plotter(off_screen=True,window_size=(620,660))
        pl.add_mesh(g,scalars="row",cmap="turbo",clim=[0,vmax],n_colors=10,nan_color="lightgray",
                    scalar_bar_args={"title":"W row-sens [um/K]","title_font_size":15,"label_font_size":12})
        pl.add_mesh(pv.Sphere(radius=0.005,center=A_XYZ),color="red")
        pl.add_mesh(pv.Sphere(radius=0.005,center=O_XYZ),color="blue")
        pl.camera_position=cam; pl.set_background("white")
        p=os.path.join(TMP,lab[:3]+".png"); pl.screenshot(p); pl.close(); shots.append((p,lab))
    fig,axes=plt.subplots(1,2,figsize=(13,6.6))
    for ax,(p,lab) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(lab,fontsize=13)
    fig.suptitle("FrontISTRの熱感度 $W=K^{-1}H$ の分布（温度→変位）と変位観測点\n"
                 "赤=W行感度が高い（変位観測に効く）／赤球=変位A(+X上面)・青球=変位O(-X上面)",fontsize=14,weight="bold")
    fig.text(0.5,0.015,"※ 105のDUMPWパッチ版FrontISTRで計算した $W=K^{-1}H$ の行ノルム（DUMPW直接出力と1.9e-7一致）。\n"
             "上端(自由端)ほど高感度で、変位観測点A/Oを上面に置くのが良い理由になる。灰色=固定部近傍の除外域。",
             ha="center",fontsize=10,color="dimgray")
    fig.tight_layout(rect=[0,0.06,1,0.92]); fig.savefig(os.path.join(IMG,"frontistr_sensitivity.png"),dpi=130); plt.close(fig)
    print("[fs] wrote docs/img/frontistr_sensitivity.png")


if __name__=="__main__": main()
