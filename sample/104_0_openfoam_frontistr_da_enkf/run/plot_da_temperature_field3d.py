"""データ同化後の温度分布を「本物の3Dソリッドメッシュ」で立体表示する.

点群ではなく、OpenFOAM の固体領域(20696セル)の実メッシュを温度で塗る。
保存済みの温度場 npy(でたらめ/同化後/真値)をメッシュのセル順に並べ替えて
セルデータとして与え、3面比較する。断面(clip)で内部の分布も見せる。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_da_temperature_field3d.py
出力: docs/img/da_temperature_field3d.png
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

K = 273.15
FIELDS = os.path.join(ROOT, "openfoam", "run_fem_enkf", "fields")
CASE = os.path.join(ROOT, "openfoam", "run_fem_enkf", "truth", "case.foam")
IMG = os.path.join(ROOT, "docs", "img")


def solid_internal_mesh():
    import pyvista as pv
    if not os.path.exists(CASE):
        open(CASE, "w").close()
    r = pv.POpenFOAMReader(CASE)
    r.set_active_time_value(r.time_values[-1])
    m = r.read()
    solid = m["solid"]
    grid = solid["internalMesh"] if "internalMesh" in solid.keys() else solid.combine()
    return grid


def main():
    import pyvista as pv
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass

    grid = solid_internal_mesh()
    print(f"[3d] solid internal mesh cells: {grid.n_cells}")

    centres = np.load(os.path.join(FIELDS, "cell_centres.npy"))
    tree = cKDTree(centres)                      # npy 側の順序
    mesh_c = np.asarray(grid.cell_centers().points)
    _, idx = tree.query(mesh_c)                  # メッシュセル → npy インデックス

    cases = [
        ("initial guess (garbage, no DA)", "ensmean_t0.npy"),
        ("after data assimilation", "ensmean_t20.npy"),
        ("truth", "truth_final.npy"),
    ]
    fields = {}
    for lab, fn in cases:
        f = np.load(os.path.join(FIELDS, fn))[idx] - K   # メッシュ順に並べ替え
        fields[lab] = f
    # 同化後と真値は共通スケール、でたらめは独自スケール
    tclim = [float(fields["truth"].min()), float(fields["truth"].max())]
    clims = {"initial guess (garbage, no DA)":
             [float(fields[cases[0][0]].min()), float(fields[cases[0][0]].max())],
             "after data assimilation": tclim, "truth": tclim}

    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1650, 680))
    for j, (lab, _fn) in enumerate(cases):
        pl.subplot(0, j)
        g = grid.copy()
        g.cell_data["T [degC]"] = fields[lab]
        # 半分を切り取って内部の温度分布も見せる(y<0 側を残す)
        clipped = g.clip(normal="y", origin=(0, 0, 0.05025))
        pl.add_mesh(clipped, scalars="T [degC]", cmap="turbo", clim=clims[lab],
                    show_edges=False, scalar_bar_args={
                        "title": f"T {lab.split('(')[0].strip()[:6]} [degC]",
                        "title_font_size": 22, "label_font_size": 18})
        pl.add_text(f"{lab}\nT {fields[lab].min():.1f}..{fields[lab].max():.1f}C",
                    font_size=20, color="black")
        pl.camera_position = [(0.24, -0.20, 0.20), (0.0, -0.01, 0.05), (0, 0, 1)]
    pl.set_background("white")
    out = os.path.join(IMG, "da_temperature_field3d.png")
    pl.screenshot(out)
    pl.close()
    print(f"[3d] wrote {os.path.relpath(out, ROOT)}")
    for lab, _ in cases:
        print(f"    {lab:32s}: T {fields[lab].min():.2f}..{fields[lab].max():.2f} C")


if __name__ == "__main__":
    main()
