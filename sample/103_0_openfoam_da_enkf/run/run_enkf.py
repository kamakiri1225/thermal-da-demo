"""アンサンブルカルマンフィルタ(EnKF)による双子実験を実行する.

手順:
    1. ROM が 102_0 に未校正なら calibrate.py を実行して config/rom_calibrated.yaml を作る
    2. 真値ラン + 合成観測を生成し、でたらめな初期アンサンブルから EnKF で同化
    3. 温度時刻歴の補正・RMSE・パラメータ推定を results/ と docs/img/ に出力

使い方(このフォルダをカレントにして):
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_enkf.py
"""

from __future__ import annotations

import os
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore import plots
from dacore.calibrate import calibrate
from dacore.twin import run_twin

RESULTS = os.path.join(ROOT, "results")
IMG = os.path.join(ROOT, "docs", "img")
FILTER = "enkf"
LABEL = "EnKF (アンサンブルカルマンフィルタ)"


def main():
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(IMG, exist_ok=True)

    with open(os.path.join(ROOT, "config", "da_config.yaml")) as f:
        cfg = yaml.safe_load(f)

    calib_path = os.path.join(ROOT, "config", "rom_calibrated.yaml")
    if not os.path.exists(calib_path):
        print("[run] 校正ファイルが無いので calibrate を実行します")
        calibrate()
    with open(calib_path) as f:
        calib = yaml.safe_load(f)

    print(f"[run] {LABEL}: n_members={cfg['ensemble']['n_members']}, "
          f"obs_nodes={cfg['observation']['nodes']}, "
          f"obs_interval={cfg['experiment']['obs_interval_s']} s")

    hist = run_twin(cfg, calib, FILTER)

    plots.save_history_csv(hist, os.path.join(RESULTS, "enkf_history.csv"))
    plots.plot_state_convergence(
        hist, f"{LABEL}: 温度の補正時刻歴(でたらめな初期状態→真値)",
        os.path.join(IMG, "enkf_state_convergence.png"))
    plots.plot_rmse(hist, f"{LABEL}: 温度 RMSE の収束",
                    os.path.join(IMG, "enkf_rmse.png"))
    plots.plot_params(hist, f"{LABEL}: パラメータ推定",
                      os.path.join(IMG, "enkf_params.png"))

    # 要約
    r0, r1 = hist["rmse_all"][0], hist["rmse_all"][-1]
    summary = {
        "filter": "EnKF",
        "n_members": cfg["ensemble"]["n_members"],
        "obs_nodes": cfg["observation"]["nodes"],
        "obs_interval_s": cfg["experiment"]["obs_interval_s"],
        "rmse_all_initial_K": float(r0),
        "rmse_all_final_K": float(r1),
        "rmse_unobs_final_K": float(hist["rmse_unobs"][-1]),
        "q_scale_final": float(hist["q_mean"][-1]),
        "q_scale_true": float(hist["q_scale_true"]),
        "h_final": float(hist["h_mean"][-1]),
        "h_true": float(hist["h_true"]),
    }
    with open(os.path.join(RESULTS, "enkf_summary.yaml"), "w") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    print(f"[run] RMSE(全ノード): 初期 {r0:.2f} K -> 最終 {r1:.3f} K")
    print(f"[run] 未観測ノード最終 RMSE: {hist['rmse_unobs'][-1]:.3f} K")
    print(f"[run] q_scale: 推定 {summary['q_scale_final']:.3f} (真値 1.0)")
    print(f"[run] h: 推定 {summary['h_final']:.4f} (真値 {summary['h_true']:.4f})")
    print(f"[run] 図を {os.path.relpath(IMG)} に、CSV/要約を {os.path.relpath(RESULTS)} に出力")


if __name__ == "__main__":
    main()
