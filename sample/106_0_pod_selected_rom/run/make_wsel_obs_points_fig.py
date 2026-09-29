"""熱感度スライドの対照実験で「どこを測っているか」を示す図.

固定 : 温度センサ P2（ヒータ側+X・中央高さ）1点  ← どの構成でも共通
振る : 変位2点の置き場所だけ
        高W = 上面 z=100.5mm の2点（行感度 6.33 / 4.45 µm/K）
        低W = 底面近く z=5mm の2点（行感度 0.124 / 0.126 µm/K）
注目量: 上面A/Oの変位差（観測には使わない）

出力: docs/img/wsel_obs_points.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_wsel_obs_points_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KINVH=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","kinvh_sensitivity.npz")
H=0.1005


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    kv=np.load(KINVH, allow_pickle=True)
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    co=np.array([p for _n,p in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,cc in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in cc)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),co)
    q=np.load(os.path.join(RES,"qdeim_points.npz"))
    P2=q["cell_centres"][q["cell_idx"][2]]
    hi=[co[i] for i in kv["hi"]]; lo=[co[i] for i in kv["lo"]]
    rs=kv["row_sens"]
    print("P2",np.round(P2*1000,1),"hi",[np.round(x*1000,1) for x in hi],
          "lo",[np.round(x*1000,1) for x in lo])

    def trim(img):
        m=np.any(img[:,:,:3]<246,axis=-1); ys,xs=np.where(m); p=6
        return img[max(0,ys.min()-p):ys.max()+p, max(0,xs.min()-p):xs.max()+p]

    def shoot(cam, show):
        pl=pv.Plotter(off_screen=True, window_size=(700,800))
        pl.add_mesh(ug, color="#c3cedb", opacity=0.24)
        # 温度センサ P2（全構成で共通・固定）
        pl.add_mesh(pv.Sphere(radius=0.0040, center=P2), color="#1F4E9C")
        pl.add_point_labels([P2+np.array([0.004,0,0.006])],["P2 (temp.)"],
                            font_size=30, text_color="#1F4E9C", shape=None,
                            always_visible=True, bold=True)
        if show in ("hi","both"):
            for k,(x,i) in enumerate(zip(hi, kv["hi"])):
                pl.add_mesh(pv.Sphere(radius=0.0040, center=x), color="#C0392B")
                pl.add_mesh(pv.Arrow(start=x, direction=(0,0,1), scale=0.026), color="#C0392B")
            # ラベルは2点まとめて1つ（近接していて重なるため）
            c=(np.array(hi[0])+np.array(hi[1]))/2
            pl.add_point_labels([c+np.array([0,0,0.034])],
                                [f"W = {rs[kv[chr(39)+chr(39)] if False else kv['hi'][0]]:.2f} / {rs[kv['hi'][1]]:.2f} um/K"],
                                font_size=28, text_color="#C0392B", shape=None,
                                always_visible=True, bold=True)
        if show in ("lo","both"):
            for k,(x,i) in enumerate(zip(lo, kv["lo"])):
                pl.add_mesh(pv.Sphere(radius=0.0040, center=x), color="#e58f2a")
                pl.add_mesh(pv.Arrow(start=x, direction=(0,0,1), scale=0.026), color="#e58f2a")
            c=(np.array(lo[0])+np.array(lo[1]))/2
            pl.add_point_labels([c+np.array([0,0,0.034])],
                                [f"W = {rs[kv['lo'][0]]:.2f} / {rs[kv['lo'][1]]:.2f} um/K"],
                                font_size=28, text_color="#b3720f", shape=None,
                                always_visible=True, bold=True)
        pl.camera_position=cam; pl.set_background("white"); pl.camera.zoom(1.15)
        img=trim(pl.screenshot(return_img=True)); pl.close(); return img

    CAM_H=[(-0.24,-0.27,0.22),(0,0,0.06),(0,0,1)]     # 高Wは −X,−Y 側
    CAM_L=[(0.30,-0.26,0.16),(0,0,0.05),(0,0,1)]      # 低Wは +X 側の底
    a=shoot(CAM_H,"hi"); b=shoot(CAM_L,"lo")

    fig,(a0,a1)=plt.subplots(1,2,figsize=(12.4,6.8))
    a0.imshow(a); a0.axis("off")
    a0.set_title("高W：上面 z=100.5 mm の2点\n(−22.8, −29.8) と (−36.2, 9.7) mm",
                 fontsize=13, weight="bold", color="#C0392B")
    a1.imshow(b); a1.axis("off")
    a1.set_title("低W：底面近く z=5 mm の2点\n(37.5, 0.0) と (33.1, 0.0) mm",
                 fontsize=13, weight="bold", color="#b3720f")
    fig.suptitle("対照実験でどこを測っているか ― 温度は P2 に固定、変位2点だけ動かす",
                 fontsize=15, weight="bold")
    fig.text(0.5,0.020,
             "青 P2 (36.5, 6.9, 50.8) mm ＝ヒータ側(+X)・中央高さの温度センサ。3構成すべてで共通。"
             "　矢印＝測る変位の向き（上下方向 $U_z$）。\n"
             "底面は固定端に近く動けないので行感度が 0.12 µm/K しかなく、上面の 6.33 µm/K とは 51倍の差がある。",
             ha="center", fontsize=11.5, color="#444")
    fig.tight_layout(rect=[0,0.055,1,0.92])
    out=os.path.join(IMG,"wsel_obs_points.png")
    fig.savefig(out, dpi=135); plt.close(fig); print("wrote", out)


if __name__=="__main__":
    main()
