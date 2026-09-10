"""観測に使った点と、確認(未観測)に使った点を中空円筒上に示す図。

観測 (EnKFに渡す):
  温度2点  … hot(+X, ヒータ側), cold(-X, 反対側)  … 赤球
  変位2点  … 上面 heater側 / opposite側 の Uz     … 橙ダイヤ
確認 (未観測、DAの当たりを検証):
  mid(+Y 90°), top(上面近傍), core(肉厚中心)      … 青球

使い方: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_obs_vs_verify.py
出力: docs/img/obs_vs_verify.png（presentation/assets にもコピー）
"""
from __future__ import annotations
import os, sys, subprocess
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore.node_locations import (OUTER_RADIUS, INNER_RADIUS, HEIGHT, NODE_XYZ,
                                   HEATER_ARC_LENGTH, HEATER_AXIAL_HEIGHT, HEATER_CENTER_Z)
import numpy as np
IMG = os.path.join(ROOT, "docs", "img")

# 変位観測点（上面 z=HEIGHT の +X / -X 側）
DISP_OBS = {"uz_heater": (0.028, 0.0, HEIGHT), "uz_opp": (-0.028, 0.0, HEIGHT)}
TEMP_OBS = ["hot", "cold"]      # 温度観測
VERIFY = ["mid", "top", "core"] # 確認(未観測)


def hollow_cylinder(pv):
    disc = pv.Disc(center=(0,0,0), inner=INNER_RADIUS, outer=OUTER_RADIUS,
                   normal=(0,0,1), r_res=6, c_res=80)
    return disc.extrude((0,0,HEIGHT), capping=True)


def heater_patch(pv):
    ha = 0.5*HEATER_ARC_LENGTH/OUTER_RADIUS
    z0, z1 = HEATER_CENTER_Z-0.5*HEATER_AXIAL_HEIGHT, HEATER_CENTER_Z+0.5*HEATER_AXIAL_HEIGHT
    na, nz = 40, 8; angs=np.linspace(-ha,ha,na); zs=np.linspace(z0,z1,nz); r=OUTER_RADIUS*1.002
    pts=[(r*np.cos(a), r*np.sin(a), z) for z in zs for a in angs]; pts=np.array(pts)
    faces=[]
    for j in range(nz-1):
        for i in range(na-1):
            p=j*na+i; faces+=[4,p,p+1,p+na+1,p+na]
    return pv.PolyData(pts, np.array(faces))


def render(pv, pl, view_top=False):
    pl.set_background("white")
    pl.add_mesh(hollow_cylinder(pv), color="lightsteelblue", opacity=0.32, smooth_shading=True)
    pl.add_mesh(heater_patch(pv), color="orange", opacity=0.85)
    rs = 0.0045
    for n in TEMP_OBS:      # 温度観測=赤
        pl.add_mesh(pv.Sphere(radius=rs, center=NODE_XYZ[n]), color="red")
        pl.add_point_labels([NODE_XYZ[n]], [f"{n} [temp obs]"], font_size=20,
                            text_color="red", shape=None, always_visible=True)
    for n,p in DISP_OBS.items():   # 変位観測=橙
        pl.add_mesh(pv.Sphere(radius=rs, center=p), color="darkorange")
        pl.add_point_labels([p], [f"{n} [disp obs]"], font_size=20,
                            text_color="darkorange", shape=None, always_visible=True)
    for n in VERIFY:        # 確認=青
        pl.add_mesh(pv.Sphere(radius=rs, center=NODE_XYZ[n]), color="royalblue")
        pl.add_point_labels([NODE_XYZ[n]], [f"{n} [verify]"], font_size=20,
                            text_color="royalblue", shape=None, always_visible=True)
    pl.add_text("red=temperature obs   orange=displacement obs   blue=verify(unobserved)   +X=heater",
                position="upper_left", font_size=11, color="black")
    if view_top:
        pl.view_xy()
    else:
        pl.camera_position=[(0.22,-0.18,0.20),(0,0,0.05),(0,0,1)]; pl.camera.zoom(1.3)


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    os.makedirs(IMG, exist_ok=True)
    pl=pv.Plotter(off_screen=True, window_size=(1400,950)); render(pv, pl)
    out=os.path.join(IMG,"obs_vs_verify.png"); pl.screenshot(out); pl.close()
    subprocess.run(["cp", out, os.path.join(ROOT,"presentation","assets","obs_vs_verify.png")], check=False)
    print(f"wrote {os.path.relpath(out, ROOT)}")


if __name__ == "__main__":
    main()
