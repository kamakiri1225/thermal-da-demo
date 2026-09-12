"""集中定数ROMを 102_0 の temperature_history.csv に合わせて同定する.

出力: config/rom_calibrated.yaml
    - C, k_circ, k_core, k_axial : 既知として固定する校正済みパラメータ
    - h_true, q_scale_true       : 双子実験(OSSE)の「真値」パラメータ
    - fit_rmse_C                 : ROM と OpenFOAM の当てはめ残差 RMSE

この校正により、ROM の「真値ラン」が本物の chtMultiRegionFoam の
温度履歴とよく一致するようになり、データ同化の観測(真値)が
102_0 の実挙動に根ざしたものになる。
"""

from __future__ import annotations

import os

import numpy as np
import yaml
from scipy.optimize import least_squares

from .cht_rom import (
    COLD,
    HOT,
    MID,
    TOP,
    ROMParams,
    integrate_single,
)

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG_DIR = os.path.join(HERE, "..", "config")
CSV_PATH = os.path.join(CONFIG_DIR, "openfoam_102_0_temperature_history.csv")
OUT_YAML = os.path.join(CONFIG_DIR, "rom_calibrated.yaml")

OBS_NODES = [HOT, MID, COLD, TOP]      # 校正で使う4観測ノード
OBS_COLS = ["T_hot_C", "T_mid_C", "T_cold_C", "T_top_C"]


def load_openfoam_curves():
    """102_0 の CSV を読み、(times[s], T_obs[nt,4] in Kelvin) を返す."""
    raw = np.genfromtxt(CSV_PATH, delimiter=",", names=True)
    times = raw["time_s"].astype(float)
    cols = np.column_stack([raw[c].astype(float) for c in OBS_COLS])
    return times, cols + 273.15


def _forward(vec, times):
    """校正ベクトル vec=[C0..C4,k_circ,k_core,k_axial,h] を前進積分し、
    観測4ノードの温度 (nt,4) [K] を返す(q_scale=1 の真値ヒータ)。"""
    h = vec[8]
    params = ROMParams.calib_vector_to_params(vec[:8], h, q_scale=1.0)
    T0 = np.full(5, 293.15)
    t1 = float(times[-1])
    # 5 s 刻みの観測時刻に必ず乗るよう dt=1 s で積分
    ts, traj = integrate_single(T0, params, 0.0, t1, dt=1.0)
    # 観測時刻へ最近傍サンプリング
    idx = np.searchsorted(ts, times)
    idx = np.clip(idx, 0, len(ts) - 1)
    return traj[idx][:, OBS_NODES]


def calibrate(verbose=True):
    times, T_obs = load_openfoam_curves()

    def residual(vec):
        model = _forward(vec, times)
        return (model - T_obs).ravel()

    # 初期値: 鋼材2.49kg, c~500 J/kgK -> C_tot~1245 J/K を core 中心に配分
    x0 = np.array([120.0, 120.0, 120.0, 120.0, 700.0,  # C
                   2.0, 5.0, 1.0,                        # k_circ, k_core, k_axial
                   0.3])                                 # h
    lb = np.array([10, 10, 10, 10, 100, 1e-3, 1e-3, 1e-3, 1e-4])
    ub = np.array([2000, 2000, 2000, 2000, 5000, 100, 100, 100, 50])

    sol = least_squares(residual, x0, bounds=(lb, ub), xtol=1e-12, ftol=1e-12)
    vec = sol.x
    rmse = float(np.sqrt(np.mean(sol.fun ** 2)))

    C = vec[:5].tolist()
    result = {
        "C": [float(c) for c in C],
        "k_circ": float(vec[5]),
        "k_core": float(vec[6]),
        "k_axial": float(vec[7]),
        "h_true": float(vec[8]),
        "q_scale_true": 1.0,
        "fit_rmse_C": rmse,
        "note": (
            "集中定数ROMを 102_0 の temperature_history.csv に同定した結果。"
            "C, k_* は既知として固定。h_true, q_scale_true は双子実験の真値。"
        ),
    }
    os.makedirs(CONFIG_DIR, exist_ok=True)
    with open(OUT_YAML, "w") as f:
        yaml.safe_dump(result, f, sort_keys=False, allow_unicode=True)

    if verbose:
        print(f"[calibrate] fit RMSE = {rmse:.4f} K  ({rmse:.4f} degC)")
        print(f"[calibrate] C       = {[round(c,1) for c in C]}")
        print(f"[calibrate] k_circ={vec[5]:.4f} k_core={vec[6]:.4f} k_axial={vec[7]:.4f}")
        print(f"[calibrate] h_true  = {vec[8]:.5f} W/K")
        print(f"[calibrate] wrote {os.path.relpath(OUT_YAML)}")
    return result


def load_calibrated():
    with open(OUT_YAML) as f:
        return yaml.safe_load(f)


if __name__ == "__main__":
    calibrate()
