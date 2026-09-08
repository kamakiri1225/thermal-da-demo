"""観測点・未観測点を中空円筒ジオメトリ上に PyVista で可視化する.

config/da_config.yaml の observation.nodes を観測点(赤)、それ以外を
未観測点(青)として、102_0 と同じ中空円筒の上に球で表示する。
ヒータ領域(外周 +X 側)も色分けして示す。

使い方(このフォルダをカレントにして):
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_probe_locations.py

出力: docs/img/probe_locations_3d.png
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore.node_locations import (  # noqa: E402
    HEATER_ARC_LENGTH,
    HEATER_AXIAL_HEIGHT,
    HEATER_CENTER_Z,
    HEIGHT,
    INNER_RADIUS,
    NODE_XYZ,
    OUTER_RADIUS,
)

IMG = os.path.join(ROOT, "docs", "img")


def _pv():
    import pyvista as pv
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass
    return pv


def hollow_cylinder(pv):
    """中空円筒(肉厚のある管)を Disc の押し出しで作る."""
    disc = pv.Disc(center=(0, 0, 0), inner=INNER_RADIUS, outer=OUTER_RADIUS,
                   normal=(0, 0, 1), r_res=6, c_res=80)
    solid = disc.extrude((0, 0, HEIGHT), capping=True)
    return solid


def heater_patch(pv):
    """外周 +X 側のヒータ領域を円弧状の帯で表す(可視化用)."""
    half_ang = 0.5 * HEATER_ARC_LENGTH / OUTER_RADIUS   # 周方向 半角 [rad]
    z0 = HEATER_CENTER_Z - 0.5 * HEATER_AXIAL_HEIGHT
    z1 = HEATER_CENTER_Z + 0.5 * HEATER_AXIAL_HEIGHT
    na, nz = 40, 8
    angs = np.linspace(-half_ang, half_ang, na)
    zs = np.linspace(z0, z1, nz)
    r = OUTER_RADIUS * 1.002    # 表面のわずかに外側
    pts, faces = [], []
    for zi in zs:
        for a in angs:
            pts.append((r * np.cos(a), r * np.sin(a), zi))
    pts = np.array(pts)
    for j in range(nz - 1):
        for i in range(na - 1):
            p0 = j * na + i
            faces += [4, p0, p0 + 1, p0 + na + 1, p0 + na]
    return pv.PolyData(pts, np.array(faces))


def main():
    os.makedirs(IMG, exist_ok=True)
    pv = _pv()

    with open(os.path.join(ROOT, "config", "da_config.yaml")) as f:
        cfg = yaml.safe_load(f)
    obs_nodes = set(cfg["observation"]["nodes"])

    pl = pv.Plotter(off_screen=True, window_size=(1400, 900))
    pl.set_background("white")

    # 試験体(半透明)
    pl.add_mesh(hollow_cylinder(pv), color="lightsteelblue", opacity=0.35,
                smooth_shading=True, show_edges=False)
    # ヒータ領域
    pl.add_mesh(heater_patch(pv), color="orange", opacity=0.9)

    # ノード(観測=赤球, 未観測=青球)
    r_sphere = 0.004
    for name, xyz in NODE_XYZ.items():
        observed = name in obs_nodes
        color = "red" if observed else "royalblue"
        pl.add_mesh(pv.Sphere(radius=r_sphere, center=xyz), color=color)
        # VTK の点ラベルは日本語グリフを持たないため ASCII 表記にする
        label = f"{name} [obs]" if observed else f"{name} [unobs]"
        pl.add_point_labels([xyz], [label], font_size=22, point_size=1,
                            text_color=color, shape=None, always_visible=True,
                            font_family="arial")

    # 凡例代わりのテキスト
    pl.add_text("red = observed (hot, cold)   blue = unobserved (mid, top, core)"
                "   orange = heater",
                position="upper_left", font_size=12, color="black")

    pl.add_axes(line_width=3, labels_off=False)
    pl.camera_position = [(0.22, -0.18, 0.19), (0.0, 0.0, 0.05), (0, 0, 1)]
    pl.camera.zoom(1.35)
    out = os.path.join(IMG, "probe_locations_3d.png")
    pl.screenshot(out)
    pl.close()
    print(f"[probe-viz] wrote {os.path.relpath(out)}")

    # 上面図(X-Y)も出力
    pl2 = pv.Plotter(off_screen=True, window_size=(1000, 1000))
    pl2.set_background("white")
    pl2.add_mesh(hollow_cylinder(pv), color="lightsteelblue", opacity=0.35)
    pl2.add_mesh(heater_patch(pv), color="orange", opacity=0.9)
    for name, xyz in NODE_XYZ.items():
        observed = name in obs_nodes
        color = "red" if observed else "royalblue"
        pl2.add_mesh(pv.Sphere(radius=r_sphere, center=xyz), color=color)
        pl2.add_point_labels([xyz], [f"{name}"], font_size=24,
                             text_color=color, shape=None, always_visible=True)
    pl2.view_xy()
    pl2.add_text("top view (X-Y),  +X = heater side",
                 position="upper_left", font_size=12, color="black")
    out2 = os.path.join(IMG, "probe_locations_top.png")
    pl2.screenshot(out2)
    pl2.close()
    print(f"[probe-viz] wrote {os.path.relpath(out2)}")


if __name__ == "__main__":
    main()
