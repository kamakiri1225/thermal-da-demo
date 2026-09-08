"""データ同化された変位「分布」を可視化する.

OpenFOAM 104 のデータ同化で得た固体温度場(全20696セル)を FrontISTR の
熱膨張線形静解析に通し、上面変位の3次元分布(変形形状)を求める。3つを比較:
    でたらめ初期のまま / データ同化後 / 真値

温度場は openfoam/run_fem_enkf/fields/*.npy(同化ドライバが保存)を使う:
    ensmean_t0.npy   … でたらめ初期(アンサンブル平均)
    ensmean_t20.npy  … データ同化後(アンサンブル平均)
    truth_final.npy  … 真値

各温度場を FrontISTR 節点へ逆距離補間 → 線形静解析 → 全節点変位、を実行する
(102_1 の連成コードを再利用)。FrontISTR は各16秒ほど。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_da_displacement_field.py
出力: docs/img/da_displacement_field.png
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
F102 = os.path.abspath(os.path.join(
    ROOT, "..", "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
sys.path.insert(0, F102)

import cylinder_mesh  # noqa: E402  (102_1)
import fistr_case     # noqa: E402  (102_1)
from openfoam_temperature import (  # noqa: E402
    align_cell_centers_to_node_mesh, interpolate_to_nodes)
from fem.fem_obs import MATERIAL  # noqa: E402

NR, NTH, NZ = 4, 48, 20
R_IN, R_OUT, H = 0.020, 0.0375, 0.1005
K = 273.15
FIELDS = os.path.join(ROOT, "openfoam", "run_fem_enkf", "fields")
IMG = os.path.join(ROOT, "docs", "img")


def fem_displacement_field(cell_centres, cell_T, work_dir):
    """セル温度場 → FrontISTR → (node_coords[N,3], disp[N,3][m], node_ids)."""
    os.makedirs(work_dir, exist_ok=True)
    mesh = fistr_case.write_mesh(
        Path(work_dir), NR, NTH, NZ, R_IN, R_OUT, H,
        young_modulus=MATERIAL["young_modulus_Pa"],
        poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],
        thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work_dir))
    node_ids = [nid for nid, _ in mesh["nodes"]]
    node_coords = np.array([xyz for _, xyz in mesh["nodes"]], dtype=float)
    aligned, _ = align_cell_centers_to_node_mesh(cell_centres, node_coords)
    node_T = interpolate_to_nodes(node_coords, aligned, cell_T, k=8)
    fistr_case.write_cnt(
        Path(work_dir), node_ids, node_T,
        reference_temperature=MATERIAL["reference_temperature_K"],
        young_modulus=MATERIAL["young_modulus_Pa"],
        poisson_ratio=MATERIAL["poisson_ratio"],
        thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.run_fistr(Path(work_dir))
    disp = fistr_case.read_displacement(Path(work_dir))
    U = np.array([disp[n] for n in node_ids], dtype=float)   # (N,3) [m]
    return node_coords, U, node_T


def main():
    import pyvista as pv
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass

    centres = np.load(os.path.join(FIELDS, "cell_centres.npy"))
    cases = [
        ("initial guess (garbage, no DA)", "ensmean_t0.npy"),
        ("after data assimilation", "ensmean_t20.npy"),
        ("truth", "truth_final.npy"),
    ]
    results = []
    for lab, fn in cases:
        cell_T = np.load(os.path.join(FIELDS, fn))
        work = os.path.join(ROOT, "openfoam", "fem_disp_" + fn.replace(".npy", ""))
        coords, U, nodeT = fem_displacement_field(centres, cell_T, work)
        uz_um = U[:, 2] * 1e6
        results.append((lab, coords, U, uz_um))
        print(f"[da-disp] {lab:32s}: Uz range [{uz_um.min():.2f}, {uz_um.max():.2f}] um  "
              f"nodeTmax={nodeT.max()-K:.2f}C", flush=True)

    # でたらめ(t0)は桁違いに大きいので独自スケール、
    # データ同化後と真値は共通スケール(真値の範囲)にして細部を比較できるようにする
    truth_uz = results[2][3]
    tclim = [float(truth_uz.min()), float(truth_uz.max())]
    clims = [[float(results[0][3].min()), float(results[0][3].max())], tclim, tclim]
    # 変形の誇張倍率は全パネル共通(相対比較のため)。最大9umを約2cmに。
    gmax = max(abs(results[0][3].min()), abs(results[0][3].max()), 1e-9)
    scale = 0.02 / gmax * 1e6

    bar_titles = ["Uz garbage [um]", "Uz DA [um]", "Uz truth [um]"]
    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1650, 660))
    for j, (lab, coords, U, uz_um) in enumerate(results):
        pl.subplot(0, j)
        warped = coords + U * scale
        cloud = pv.PolyData(warped)
        cloud[bar_titles[j]] = uz_um
        pl.add_mesh(cloud, scalars=bar_titles[j], cmap="coolwarm", clim=clims[j],
                    render_points_as_spheres=True, point_size=6,
                    scalar_bar_args={"title": bar_titles[j],
                                     "title_font_size": 16, "label_font_size": 13})
        rng = f"Uz {uz_um.min():.2f}..{uz_um.max():.2f} um"
        pl.add_text(f"{lab}\n{rng}", font_size=9, color="black")
        pl.camera_position = [(0.22, -0.18, 0.19), (0.0, 0.0, 0.05), (0, 0, 1)]
    pl.set_background("white")
    out = os.path.join(IMG, "da_displacement_field.png")
    pl.screenshot(out)
    pl.close()
    print(f"[da-disp] wrote {os.path.relpath(out, ROOT)} "
          f"(変形 {scale:.0f}倍に誇張)")


if __name__ == "__main__":
    main()
