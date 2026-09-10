"""ROM版の「メモリ上のアンサンブル」を実際に CSV へ書き出して可視化する.

普段 ROM は途中経過をディスクに残さない(メモリ内で計算)。この見えない中身を
「実態」として確認できるよう、各同化サイクルの 60メンバー × [T5点, q_scale, h] を
CSV に保存する。OpenFOAM が member_00/ 等のフォルダを作るのに相当するものを、
ROM でも目に見える形にする。

出力: results/rom_ensemble/cycle_00_before.csv, cycle_00_after.csv, ...
      各ファイル = 60行(メンバー) × 7列 + 末尾に mean / std 行
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore.calibrate import load_calibrated
from dacore.cht_rom import N_NODES, NODE_NAMES
from dacore.displacement import load_operator
from dacore.enkf import enkf_update
from dacore.ensemble import I_H, I_QSCALE, clip_params, forecast, init_ensemble
from dacore.observations import generate_truth, obs_node_indices
from dacore.twin_fem import build_H_R

K = 273.15
OUT = os.path.join(ROOT, "results", "rom_ensemble")
COLS = [f"T_{n}_degC" for n in NODE_NAMES] + ["q_scale", "h_W/K"]


def _save(Z, path, note):
    """Z (60,7) を CSV に。温度は degC、末尾に mean/std 行を付ける。"""
    A = Z.copy()
    A[:, :N_NODES] -= K
    mean = A.mean(axis=0)
    std = A.std(axis=0)
    with open(path, "w") as f:
        f.write(f"# {note}\n")
        f.write("member," + ",".join(COLS) + "\n")
        for i, row in enumerate(A):
            f.write(f"{i}," + ",".join(f"{x:.4g}" for x in row) + "\n")
        f.write("mean," + ",".join(f"{x:.4g}" for x in mean) + "\n")
        f.write("std," + ",".join(f"{x:.4g}" for x in std) + "\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "da_config.yaml")))
    calib = load_calibrated()
    op = load_operator()

    H, R, n_t = build_H_R(cfg)
    tnodes = obs_node_indices(cfg)
    fb = cfg["filter"]
    rng_obs = np.random.default_rng(cfg["experiment"]["seed"])
    rng = np.random.default_rng(cfg["experiment"]["seed"] + 1)

    times, truth = generate_truth(calib, cfg)
    dt_obs = cfg["experiment"]["obs_interval_s"]
    cyc_t = np.arange(dt_obs, cfg["experiment"]["t_end_s"] + 1e-9, dt_obs)
    cyc_idx = np.clip(np.searchsorted(times, cyc_t), 0, len(times) - 1)

    Z = init_ensemble(cfg, rng)
    _save(Z, os.path.join(OUT, "cycle_00_initial.csv"),
          "t=0s でたらめ初期アンサンブル(まだ同化していない)")
    print(f"[dump] {Z.shape[0]}メンバー × {Z.shape[1]}列 (= 温度5点 + q_scale + h)")

    t_prev = 0.0
    for k, (t1, ci) in enumerate(zip(cyc_t, cyc_idx)):
        Z = forecast(Z, t_prev, t1, calib, cfg, rng)
        _save(Z, os.path.join(OUT, f"cycle_{k+1:02d}_before.csv"),
              f"t={t1:g}s 予報後(観測で補正する前)")
        y = truth[ci][tnodes] + rng_obs.normal(0, cfg["observation"]["noise_C"], n_t)
        Z = enkf_update(Z, y, H, R, rng, inflation=fb["inflation"])
        clip_params(Z, cfg)
        _save(Z, os.path.join(OUT, f"cycle_{k+1:02d}_after.csv"),
              f"t={t1:g}s 解析後(EnKFで補正済み) 観測y={np.round(y-K,3)}")
        t_prev = t1

    n_files = len(os.listdir(OUT))
    print(f"[dump] {n_files}個のCSVを {os.path.relpath(OUT, ROOT)} に書き出した")
    print(f"[dump] 例: cycle_01_before.csv(補正前) と cycle_01_after.csv(補正後)を見比べる")


if __name__ == "__main__":
    main()
