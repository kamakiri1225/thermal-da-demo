"""ROM のデータ同化後の温度「分布」を中空円筒上に可視化する.

ROM は5ノード(hot/mid/cold/top/core)なので、その5点の温度を円筒形状へ
逆距離補間(IDW)して連続的な温度分布として描く。指定時刻(既定 t=300s, 加熱ピーク)で
  でたらめ初期のまま(同化なし) / データ同化後 / 真値
の3面を比較する。

注意: これは5ノードからの再構成であり、細かな分布は OpenFOAM 版
(docs/img/openfoam_fem_enkf_field.png, 20696セル)が正確。ここでは
「ROMで同化した温度がどんな空間分布になるか」を掴むためのもの。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_rom_field.py [時刻s]
出力: docs/img/rom_fem_field.png
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore.calibrate import load_calibrated
from dacore.displacement import load_operator
from dacore.node_locations import (
    HEIGHT, INNER_RADIUS, NODE_XYZ, OUTER_RADIUS,
)
from dacore.twin_fem import run_rom_fem_twin

K = 273.15
IMG = os.path.join(ROOT, "docs", "img")


def cylinder_points(n_r=6, n_th=90, n_z=26):
    """中空円筒の体積を点群でサンプルして (M,3) を返す."""
    rs = np.linspace(INNER_RADIUS, OUTER_RADIUS, n_r)
    th = np.linspace(0, 2*np.pi, n_th, endpoint=False)
    zs = np.linspace(0, HEIGHT, n_z)
    pts = []
    for z in zs:
        for t in th:
            for r in rs:
                pts.append((r*np.cos(t), r*np.sin(t), z))
    return np.array(pts)


def idw_reconstruct(points, node_xyz, node_T, power=2.0):
    """5ノード温度を points へ逆距離補間する (M,)."""
    nodes = np.array(list(node_xyz.values()))
    T = np.asarray(node_T)
    out = np.empty(len(points))
    for m, p in enumerate(points):
        d = np.linalg.norm(nodes - p, axis=1)
        if d.min() < 1e-9:
            out[m] = T[d.argmin()]
        else:
            w = 1.0 / d**power
            out[m] = (w @ T) / w.sum()
    return out


def main():
    t_show = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    with open(os.path.join(ROOT, "config", "da_config.yaml")) as f:
        cfg = yaml.safe_load(f)
    calib = load_calibrated()
    op = load_operator()
    hist = run_rom_fem_twin(cfg, calib, op)

    ti = int(np.argmin(np.abs(hist["times"] - t_show)))
    t_actual = hist["times"][ti]
    panels = [
        ("initial guess (garbage, no DA)", hist["free_T"][ti]),
        ("after data assimilation", hist["da_T"][ti]),
        ("truth", hist["truth_T"][ti]),
    ]

    import pyvista as pv
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass

    pts = cylinder_points()
    clim = [min(hist["truth_T"][ti].min(), hist["da_T"][ti].min()) - K,
            max(panels[0][1].max(), hist["truth_T"][ti].max()) - K]

    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1650, 640))
    for j, (lab, nodeT) in enumerate(panels):
        pl.subplot(0, j)
        field = idw_reconstruct(pts, NODE_XYZ, nodeT) - K
        cloud = pv.PolyData(pts)
        cloud["T_degC"] = field
        pl.add_mesh(cloud, scalars="T_degC", cmap="turbo", clim=clim,
                    render_points_as_spheres=True, point_size=9,
                    scalar_bar_args={"title": "T [degC]"})
        pl.add_text(lab, font_size=20, color="black")
        pl.camera_position = [(0.22, -0.18, 0.19), (0.0, 0.0, 0.05), (0, 0, 1)]
    pl.set_background("white")
    out = os.path.join(IMG, "rom_fem_field.png")
    pl.screenshot(out)
    pl.close()
    print(f"[rom-field] t={t_actual:g}s  wrote {os.path.relpath(out, ROOT)}")
    print(f"[rom-field] node T (hot,mid,cold,top,core) degC:")
    for lab, nt in panels:
        print(f"    {lab:32s}: {np.round(nt-K,2)}")


if __name__ == "__main__":
    main()
