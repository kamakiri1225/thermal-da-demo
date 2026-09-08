"""EnKF と PF を同一の真値・観測で走らせ、収束性を比較する.

EnKF の設定は隣の 103_0_openfoam_da_enkf/config/da_config.yaml を、
PF の設定は本フォルダの config/da_config.yaml を用いる。
真値・観測は seed が同じなので、両フィルタはまったく同じ観測を同化する。

出力: docs/img/enkf_vs_pf_rmse.png と results/compare_summary.yaml
"""

from __future__ import annotations

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore import plots  # noqa: F401  (日本語フォント設定を有効化)
from dacore.twin import run_twin

ENKF_DIR = os.path.join(ROOT, "..", "103_0_openfoam_da_enkf")
IMG = os.path.join(ROOT, "docs", "img")
RESULTS = os.path.join(ROOT, "results")


def _load(path):
    with open(path) as f:
        return yaml.safe_load(f)


def main():
    os.makedirs(IMG, exist_ok=True)
    os.makedirs(RESULTS, exist_ok=True)

    calib = _load(os.path.join(ROOT, "config", "rom_calibrated.yaml"))
    cfg_pf = _load(os.path.join(ROOT, "config", "da_config.yaml"))
    cfg_enkf = _load(os.path.join(ENKF_DIR, "config", "da_config.yaml"))

    h_enkf = run_twin(cfg_enkf, calib, "enkf")
    h_pf = run_twin(cfg_pf, calib, "pf")

    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    ax.plot(h_enkf["times"], h_enkf["rmse_all"], "-o", ms=3,
            label=f"EnKF (N={cfg_enkf['ensemble']['n_members']})")
    ax.plot(h_pf["times"], h_pf["rmse_all"], "-s", ms=3,
            label=f"PF (N={cfg_pf['ensemble']['n_members']})")
    ax.axhline(cfg_pf["observation"]["noise_C"], color="gray", ls=":",
               label=f"観測ノイズ {cfg_pf['observation']['noise_C']} K")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("温度 RMSE(全ノード)[K]")
    ax.set_yscale("log")
    ax.set_title("EnKF vs PF: でたらめな初期状態からの収束比較")
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(IMG, "enkf_vs_pf_rmse.png"), dpi=130)
    plt.close(fig)

    summary = {
        "enkf": {
            "n": cfg_enkf["ensemble"]["n_members"],
            "rmse_all_final_K": float(h_enkf["rmse_all"][-1]),
            "rmse_unobs_final_K": float(h_enkf["rmse_unobs"][-1]),
        },
        "pf": {
            "n": cfg_pf["ensemble"]["n_members"],
            "rmse_all_final_K": float(h_pf["rmse_all"][-1]),
            "rmse_unobs_final_K": float(h_pf["rmse_unobs"][-1]),
        },
    }
    with open(os.path.join(RESULTS, "compare_summary.yaml"), "w") as f:
        yaml.safe_dump(summary, f, sort_keys=False, allow_unicode=True)

    print("[compare] EnKF final RMSE(all) = "
          f"{summary['enkf']['rmse_all_final_K']:.3f} K")
    print("[compare] PF   final RMSE(all) = "
          f"{summary['pf']['rmse_all_final_K']:.3f} K")
    print(f"[compare] 図を {os.path.relpath(os.path.join(IMG,'enkf_vs_pf_rmse.png'))} に出力")


if __name__ == "__main__":
    main()
