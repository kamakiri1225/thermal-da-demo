"""データ同化後 vs 正解 の時刻歴を GIF アニメーションで比較する.

ROM の 0→600s 時刻歴(31ステップ)を円筒メッシュに再構成し、各時刻を2面
(左=データ同化後 / 右=真値)で描いて GIF にする。

  docs/img/da_vs_truth_temperature.gif … 温度分布(断面)。IDW補間のみで高速
  docs/img/da_vs_truth_displacement.gif … 変位(変形形状)。各フレームで FrontISTR

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_da_gifs.py
"""

from __future__ import annotations

import glob
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
F102 = os.path.abspath(os.path.join(
    ROOT, "..", "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
sys.path.insert(0, F102)

import cylinder_mesh   # noqa: E402  (102_1)
import fistr_case      # noqa: E402  (102_1)
import vtk             # noqa: E402
from dacore.calibrate import load_calibrated       # noqa: E402
from dacore.displacement import load_operator      # noqa: E402
from dacore.node_locations import NODE_XYZ         # noqa: E402
from dacore.twin_fem import run_rom_fem_twin        # noqa: E402
from fem.fem_obs import MATERIAL                    # noqa: E402

K = 273.15
IMG = os.path.join(ROOT, "docs", "img")
NR, NTH, NZ = 4, 48, 20
R_IN, R_OUT, H = 0.020, 0.0375, 0.1005


def cylinder_grid():
    import pyvista as pv
    mesh = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, R_IN, R_OUT, H)
    coords = np.array([xyz for _n, xyz in mesh["nodes"]], dtype=float)
    id_row = {nid: i for i, (nid, _x) in enumerate(mesh["nodes"])}
    cells = []
    for _e, conn in mesh["elements"]:
        cells.append(8)
        cells.extend(id_row[n] for n in conn)
    ctypes = np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, np.uint8)
    grid = pv.UnstructuredGrid(np.array(cells), ctypes, coords)
    node_ids = [nid for nid, _x in mesh["nodes"]]
    return grid, coords, node_ids


def idw_weights(coords, rom_nodes):
    W = np.zeros((len(coords), len(rom_nodes)))
    for i, p in enumerate(coords):
        d = np.linalg.norm(rom_nodes - p, axis=1)
        if d.min() < 1e-9:
            W[i, d.argmin()] = 1.0
        else:
            w = 1.0 / d**2
            W[i] = w / w.sum()
    return W


def _make_gif(frame_glob, out_gif, delay=18):
    frames = sorted(glob.glob(frame_glob))
    subprocess.run(["convert", "-delay", str(delay), "-loop", "0",
                    *frames, "-layers", "Optimize", "-colors", "128", out_gif],
                   check=True)
    for f in frames:
        os.remove(f)
    print(f"[gif] wrote {os.path.relpath(out_gif, ROOT)} ({len(frames)} frames)")


def temperature_gif(hist, grid, W, tmp):
    import pyvista as pv
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass
    clim = [float(hist["truth_T"].min()-K), float(hist["truth_T"].max()-K)]
    for k, t in enumerate(hist["times"]):
        pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1100, 620))
        for j, (lab, nodeT) in enumerate([("data assimilation", hist["da_T"][k]),
                                          ("truth", hist["truth_T"][k])]):
            pl.subplot(0, j)
            g = grid.copy()
            g.point_data["T [degC]"] = W @ (nodeT - K)
            clipped = g.clip(normal="y", origin=(0, 0, 0.05025))
            pl.add_mesh(clipped, scalars="T [degC]", cmap="turbo", clim=clim,
                        scalar_bar_args={"title": "T [degC]",
                                         "title_font_size": 20, "label_font_size": 16})
            pl.add_text(f"{lab}", font_size=20, color="black")
            pl.camera_position = [(0.24, -0.20, 0.20), (0, -0.01, 0.05), (0, 0, 1)]
        pl.add_text(f"t = {t:g} s", position="lower_edge", font_size=18, color="black")
        pl.set_background("white")
        pl.screenshot(os.path.join(tmp, f"T_{k:03d}.png"))
        pl.close()
    _make_gif(os.path.join(tmp, "T_*.png"),
              os.path.join(IMG, "da_vs_truth_temperature.gif"))


def displacement_gif(hist, grid, coords, node_ids, W, tmp):
    import pyvista as pv
    work = os.path.join(ROOT, "openfoam", "fem_gif")
    os.makedirs(work, exist_ok=True)
    # メッシュは1回だけ書く。以降は cnt+run のみで高速化。
    fistr_case.write_mesh(
        Path(work), NR, NTH, NZ, R_IN, R_OUT, H,
        young_modulus=MATERIAL["young_modulus_Pa"], poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],
        thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))

    def fistr_U(node_T_K):
        fistr_case.write_cnt(
            Path(work), node_ids, node_T_K,
            reference_temperature=MATERIAL["reference_temperature_K"],
            young_modulus=MATERIAL["young_modulus_Pa"], poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work))
        disp = fistr_case.read_displacement(Path(work))
        return np.array([disp[n] for n in node_ids], dtype=float)

    # 全フレームの変位を先に計算(共通スケール決定のため)
    U_da, U_tr = [], []
    for k in range(len(hist["times"])):
        U_da.append(fistr_U(W @ hist["da_T"][k]))
        U_tr.append(fistr_U(W @ hist["truth_T"][k]))
        print(f"[gif-disp] FrontISTR {k+1}/{len(hist['times'])}", flush=True)
    U_da = np.array(U_da); U_tr = np.array(U_tr)
    uz_all = np.concatenate([U_tr[..., 2].ravel()]) * 1e6
    clim = [float(uz_all.min()), float(uz_all.max())]
    gmax = max(abs(clim[0]), abs(clim[1]), 1e-9)
    scale = 0.02 / gmax * 1e6

    pv.OFF_SCREEN = True
    for k, t in enumerate(hist["times"]):
        pl = pv.Plotter(off_screen=True, shape=(1, 2), window_size=(1100, 640))
        for j, (lab, U) in enumerate([("data assimilation", U_da[k]),
                                      ("truth", U_tr[k])]):
            pl.subplot(0, j)
            g = grid.copy()
            g.points = coords + U * scale
            g["Uz [um]"] = U[:, 2] * 1e6
            pl.add_mesh(g, scalars="Uz [um]", cmap="coolwarm", clim=clim,
                        scalar_bar_args={"title": "Uz [um]",
                                         "title_font_size": 20, "label_font_size": 16})
            pl.add_text(f"{lab}", font_size=20, color="black")
            pl.camera_position = [(0.24, -0.20, 0.20), (0, -0.01, 0.05), (0, 0, 1)]
        pl.add_text(f"t = {t:g} s  (warp x{scale:.0f})", position="lower_edge",
                    font_size=18, color="black")
        pl.set_background("white")
        pl.screenshot(os.path.join(tmp, f"U_{k:03d}.png"))
        pl.close()
    _make_gif(os.path.join(tmp, "U_*.png"),
              os.path.join(IMG, "da_vs_truth_displacement.gif"))


def main():
    os.makedirs(IMG, exist_ok=True)
    tmp = os.path.join(ROOT, "openfoam", "gif_frames")
    os.makedirs(tmp, exist_ok=True)

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "da_config.yaml")))
    hist = run_rom_fem_twin(cfg, load_calibrated(), load_operator())
    grid, coords, node_ids = cylinder_grid()
    W = idw_weights(coords, np.array(list(NODE_XYZ.values())))

    temperature_gif(hist, grid, W, tmp)
    displacement_gif(hist, grid, coords, node_ids, W, tmp)
    print("[gif] done")


if __name__ == "__main__":
    main()
