"""粒子フィルタ(PF, ブートストラップ/SIR)による双子実験を実行する.

手順は 103_0 の EnKF 版と同じで、解析ステップだけ粒子フィルタに差し替える。
真値ラン・合成観測は seed が同じなので EnKF とまったく同じ観測を同化し、
両フィルタを公平に比較できる。

使い方(このフォルダをカレントにして):
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_pf.py
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
FILTER = "pf"
LABEL = "PF (粒子フィルタ)"


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

    print(f"[run] {LABEL}: n_particles={cfg['ensemble']['n_members']}, "
          f"obs_nodes={cfg['observation']['nodes']}, "
          f"obs_interval={cfg['experiment']['obs_interval_s']} s")

    hist = run_twin(cfg, calib, FILTER)

    plots.save_history_csv(hist, os.path.join(RESULTS, "pf_history.csv"))
    plots.plot_state_convergence(
        hist, f"{LABEL}: 温度の補正時刻歴(でたらめな初期状態→真値)",
        os.path.join(IMG, "pf_state_convergence.png"))
    plots.plot_rmse(hist, f"{LABEL}: 温度 RMSE の収束",
                    os.path.join(IMG, "pf_rmse.png"))
    plots.plot_params(hist, f"{LABEL}: パラメータ推定",
                      os.path.join(IMG, "pf_params.png"))

    r0, r1 = hist["rmse_all"][0], hist["rmse_all"][-1]
    ess_min = float(min(x for x in hist["ess"][1:] if x == x))  # NaN除外
    summary = {
        "filter": "PF",
        "n_particles": cfg["ensemble"]["n_members"],
        "obs_nodes": cfg["observation"]["nodes"],
        "obs_interval_s": cfg["experiment"]["obs_interval_s"],
        "rmse_all_initial_K": float(r0),
        "rmse_all_final_K": float(r1),
        "rmse_unobs_final_K": float(hist["rmse_unobs"][-1]),
        "ess_min": ess_min,
        "q_scale_final": float(hist["q_mean"][-1]),
        "q_scale_true": float(hist["q_scale_true"]),
        "h_final": float(hist["h_mean"][-1]),
        "h_true": float(hist["h_true"]),
    }
    with open(os.path.join(RESULTS, "pf_summary.yaml"), "w") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    print(f"[run] RMSE(全ノード): 初期 {r0:.2f} K -> 最終 {r1:.3f} K")
    print(f"[run] 未観測ノード最終 RMSE: {hist['rmse_unobs'][-1]:.3f} K")
    print(f"[run] 最小 ESS: {ess_min:.1f} / {cfg['ensemble']['n_members']}")
    print(f"[run] q_scale: 推定 {summary['q_scale_final']:.3f} (真値 1.0)")
    print(f"[run] h: 推定 {summary['h_final']:.4f} (真値 {summary['h_true']:.4f})")
    print(f"[run] 図を {os.path.relpath(IMG)} に、CSV/要約を {os.path.relpath(RESULTS)} に出力")


if __name__ == "__main__":
    main()
