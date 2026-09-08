"""OpenFOAM×FrontISTR データ同化(温度＋変位観測 EnKF)を実行する.

103 の EnKF に変位観測(FrontISTR 上面 Uz 2点)を追加した双子実験。

使い方(このフォルダをカレントにして):
    cd openfoam && sh setup_base_case.sh && cd ..
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
      nohup python3 run/run_openfoam_fem_enkf.py > openfoam/run_fem_enkf.log 2>&1 &
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from daof import of_plots
from daof.of_case import BASE_CASE
from daof.of_fem_twin import run_fem_twin

RESULTS = os.path.join(ROOT, "results")
IMG = os.path.join(ROOT, "docs", "img")
WORKDIR = os.path.join(ROOT, "openfoam", "run_fem_enkf")


def plot_displacement(hist, out_path):
    import matplotlib.pyplot as plt
    t = np.array(hist["times"])
    labels = ["ヒータ側上面 Uz", "反ヒータ側上面 Uz"]
    fig, ax = plt.subplots(figsize=(8.5, 5))
    colors = ["tab:orange", "tab:purple"]
    for j, lab in enumerate(labels):
        m = hist["uz_mean"][:, j] * 1000.0   # mm -> um
        s = hist["uz_std"][:, j] * 1000.0
        tr = hist["uz_truth"][:, j] * 1000.0
        ax.fill_between(t, m - s, m + s, alpha=0.15, color=colors[j])
        ax.plot(t, m, "-o", ms=4, color=colors[j], label=f"{lab} (アンサンブル)")
        ax.plot(t, tr, "--", color=colors[j], lw=2, label=f"{lab} (真値)")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("上面変位 Uz [um]")
    ax.set_title("OpenFOAM×FrontISTR×EnKF: 変位観測の予報値と真値")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def main():
    if not os.path.isdir(os.path.join(BASE_CASE, "constant", "solid", "polyMesh")):
        raise SystemExit("base_case のメッシュがありません。先に "
                         "`sh openfoam/setup_base_case.sh` を実行してください。")
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(IMG, exist_ok=True)

    with open(os.path.join(ROOT, "openfoam", "da_openfoam_config.yaml")) as f:
        cfg = yaml.safe_load(f)

    hist = run_fem_twin(cfg, WORKDIR)

    # RMSE / Q 図(103 と同じ形式)
    hist_c = dict(hist)
    hist_c["rmse_obs"] = hist["rmse_field"]  # of_plots 互換(観測セル別は省略)
    hist_c["ess"] = [np.nan] * len(hist["times"])
    of_plots.save_history_csv(
        {**hist_c, "times": hist["times"]},
        os.path.join(RESULTS, "openfoam_fem_enkf_history.csv"))
    of_plots.plot_field_rmse(
        hist_c, "OpenFOAM×FrontISTR×EnKF: 固体温度場 RMSE(温度2点+変位2点観測)",
        os.path.join(IMG, "openfoam_fem_enkf_rmse.png"))
    of_plots.plot_Q(
        hist, "OpenFOAM×FrontISTR×EnKF: ヒータ発熱 Q の推定",
        os.path.join(IMG, "openfoam_fem_enkf_Q.png"))
    of_plots.field_snapshots(
        os.path.join(WORKDIR, "fields"), "",
        os.path.join(IMG, "openfoam_fem_enkf_field.png"))
    plot_displacement(hist, os.path.join(IMG, "openfoam_fem_enkf_disp.png"))

    summary = {
        "filter": "EnKF (field-space, OpenFOAM + FrontISTR displacement obs)",
        "n_members": cfg["ensemble"]["n_members"],
        "solid_cells": hist["Nc"],
        "obs": ["T_hot", "T_cold", "uz_heater(FrontISTR)", "uz_opposite(FrontISTR)"],
        "t_end_s": cfg["experiment"]["t_end_s"],
        "rmse_field_initial_K": float(hist["rmse_field"][0]),
        "rmse_field_final_K": float(hist["rmse_field"][-1]),
        "Q_final_W": float(hist["q_mean"][-1]),
        "Q_true_W": float(hist["Q_true"]),
    }
    with open(os.path.join(RESULTS, "openfoam_fem_enkf_summary.yaml"), "w") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
    print("[run] summary:", summary)


if __name__ == "__main__":
    main()
