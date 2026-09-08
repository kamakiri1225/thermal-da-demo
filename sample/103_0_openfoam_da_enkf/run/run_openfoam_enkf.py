"""OpenFOAM-in-the-loop の EnKF データ同化を実行する(field-space).

各メンバーが実際に chtMultiRegionFoam(既存メッシュ・輻射あり)を回し、
固体温度場全セル＋ヒータ発熱 Q を EnKF で更新する。既存メッシュ・直列前提で
1メンバーの短ウィンドウでも数分かかるため、まず openfoam/setup_base_case.sh で
ベースケースを用意してから実行すること。

使い方(このフォルダをカレントにして):
    sh openfoam/setup_base_case.sh
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_openfoam_enkf.py

計算が長いのでバックグラウンド実行推奨:
    nohup python3 run/run_openfoam_enkf.py > openfoam/run_enkf.log 2>&1 &
"""

from __future__ import annotations

import os
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from daof import of_plots
from daof.of_case import BASE_CASE
from daof.of_twin import run_openfoam_twin

RESULTS = os.path.join(ROOT, "results")
IMG = os.path.join(ROOT, "docs", "img")
WORKDIR = os.path.join(ROOT, "openfoam", "run_enkf")


def main():
    if not os.path.isdir(os.path.join(BASE_CASE, "constant", "solid", "polyMesh")):
        raise SystemExit("base_case のメッシュがありません。先に "
                         "`sh openfoam/setup_base_case.sh` を実行してください。")
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(IMG, exist_ok=True)

    with open(os.path.join(ROOT, "openfoam", "da_openfoam_config.yaml")) as f:
        cfg = yaml.safe_load(f)

    hist = run_openfoam_twin(cfg, "enkf", WORKDIR)

    of_plots.save_history_csv(hist, os.path.join(RESULTS, "openfoam_enkf_history.csv"))
    of_plots.plot_field_rmse(
        hist, "OpenFOAM×EnKF: 固体温度場 RMSE の収束(でたらめ初期→真値)",
        os.path.join(IMG, "openfoam_enkf_rmse.png"))
    of_plots.plot_Q(
        hist, "OpenFOAM×EnKF: ヒータ発熱 Q の推定",
        os.path.join(IMG, "openfoam_enkf_Q.png"))
    of_plots.field_snapshots(
        os.path.join(WORKDIR, "fields"),
        "OpenFOAM×EnKF: 固体温度場のデータ同化",
        os.path.join(IMG, "openfoam_enkf_field.png"))

    summary = {
        "filter": "EnKF (field-space, OpenFOAM-in-the-loop)",
        "n_members": cfg["ensemble"]["n_members"],
        "solid_cells": hist["Nc"],
        "obs_probes": list(cfg["observation"]["probes"]),
        "t_end_s": cfg["experiment"]["t_end_s"],
        "obs_interval_s": cfg["experiment"]["obs_interval_s"],
        "rmse_field_initial_K": float(hist["rmse_field"][0]),
        "rmse_field_final_K": float(hist["rmse_field"][-1]),
        "Q_final_W": float(hist["q_mean"][-1]),
        "Q_true_W": float(hist["Q_true"]),
    }
    with open(os.path.join(RESULTS, "openfoam_enkf_summary.yaml"), "w") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)
    print("[run] summary:", summary)


if __name__ == "__main__":
    main()
