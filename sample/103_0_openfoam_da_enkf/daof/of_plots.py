"""OpenFOAM field-space データ同化の可視化(RMSE・Q推定・温度場スナップショット)."""

from __future__ import annotations

import glob
import os
import sys

import numpy as np

# 日本語フォント設定を dacore.plots から流用
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt

K = 273.15


def plot_field_rmse(hist, title, out_path):
    t = np.array(hist["times"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t, hist["rmse_field"], "-o", ms=4, label="固体温度場 全セル")
    ax.plot(t, hist["rmse_obs"], "-s", ms=4, label="観測セル")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("温度 RMSE [K]  (アンサンブル平均 vs 真値場)")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_Q(hist, title, out_path):
    t = np.array(hist["times"])
    m = np.array(hist["q_mean"])
    s = np.array(hist["q_std"])
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.fill_between(t, m - s, m + s, alpha=0.2, color="tab:green")
    ax.plot(t, m, "-o", color="tab:green", label="推定 Q 平均 ±1σ")
    ax.axhline(hist["Q_true"], color="k", ls="--", label=f"真値 {hist['Q_true']} W")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("ヒータ発熱 Q [W]")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_history_csv(hist, out_path):
    t = hist["times"]
    header = "time_s,rmse_field_K,rmse_obs_K,q_mean_W,q_std_W,ess"
    rows = np.column_stack([t, hist["rmse_field"], hist["rmse_obs"],
                            hist["q_mean"], hist["q_std"], hist["ess"]])
    np.savetxt(out_path, rows, delimiter=",", header=header, comments="", fmt="%.6g")


def field_snapshots(fields_dir, title, out_path):
    """でたらめ初期 / DA後 / 真値 の固体温度場をコンタ面で3面比較(PyVista).

    OpenFOAMのセル中心値を、同一形状の六面体FEMメッシュ節点へ最近傍で写して
    サーフェス表示する（点群表示は分かりにくいため廃止）。
    """
    import pyvista as pv
    import vtk
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(__file__), "..", "..",
        "102_1_frontistr_hollow_cylinder_thermal_expansion", "python")))
    import cylinder_mesh

    centres = np.load(os.path.join(fields_dir, "cell_centres.npy"))
    t0 = np.load(os.path.join(fields_dir, "ensmean_t0.npy"))
    truth = np.load(os.path.join(fields_dir, "truth_final.npy"))
    finals = sorted(glob.glob(os.path.join(fields_dir, "ensmean_t*.npy")),
                    key=lambda p: float(p.split("ensmean_t")[1][:-4]))
    da_final = np.load(finals[-1])

    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    coords = np.array([xyz for _n, xyz in mesh["nodes"]])
    idr = {nid: i for i, (nid, _x) in enumerate(mesh["nodes"])}
    cells = []
    for _e, conn in mesh["elements"]:
        cells.append(8)
        cells.extend(idr[n] for n in conn)
    grid = pv.UnstructuredGrid(
        np.array(cells),
        np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, np.uint8), coords)
    _, nearest = cKDTree(centres).query(coords)   # 節点→最寄りセル中心

    panels = [("initial guess (garbage, ensemble mean)", t0),
              ("after DA (ensemble mean)", da_final),
              ("truth", truth)]
    clim = [min(truth.min(), da_final.min()) - K,
            max(truth.max(), da_final.max()) - K]

    pl = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1650, 620))
    for j, (lab, field) in enumerate(panels):
        pl.subplot(0, j)
        g = grid.copy()
        g.point_data["T_degC"] = field[nearest] - K
        pl.add_mesh(g, scalars="T_degC", cmap="turbo", clim=clim, n_colors=14,
                    scalar_bar_args={"title": "T [degC]"})
        pl.add_text(lab, font_size=9, color="black")
        pl.camera_position = [(0.22, -0.18, 0.19), (0.0, 0.0, 0.05), (0, 0, 1)]
    pl.set_background("white")
    pl.screenshot(out_path)
    pl.close()
