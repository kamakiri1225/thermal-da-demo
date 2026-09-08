"""双子実験(OSSE)ドライバ: EnKF / PF を同一条件で走らせて履歴を記録する.

真値ランと合成観測は seed から決定的に生成するため、EnKF と PF は
「まったく同じ観測」を同化する。フィルタだけを differ させて比較できる。

記録する履歴(各同化サイクル時刻ごと):
    - 解析アンサンブル平均・標準偏差(5ノード温度)
    - 予報(更新前)アンサンブル平均(補正の効き具合を見るため)
    - パラメータ推定(q_scale, h)の平均・標準偏差
    - 真値に対する温度 RMSE(全ノード / 観測ノード / 未観測ノード)
"""

from __future__ import annotations

import numpy as np

from . import enkf as enkf_mod
from . import pf as pf_mod
from .cht_rom import N_NODES, NODE_NAMES
from .ensemble import I_H, I_QSCALE, N_AUG, clip_params, forecast, init_ensemble
from .observations import (
    generate_truth,
    make_observations,
    obs_cov,
    obs_matrix,
    obs_node_indices,
)


def _rmse(mean_T, truth_T, node_subset):
    d = mean_T[node_subset] - truth_T[node_subset]
    return float(np.sqrt(np.mean(d ** 2)))


def run_twin(cfg, calib, filter_kind):
    """filter_kind: 'enkf' または 'pf'. history(dict) を返す."""
    seed = cfg["experiment"]["seed"]
    rng_obs = np.random.default_rng(seed)          # 真値・観測用(両フィルタで共通)
    rng_filt = np.random.default_rng(seed + 1)     # アンサンブル・フィルタ用

    # 1. 真値ランと合成観測
    times, truth = generate_truth(calib, cfg)
    obs = make_observations(times, truth, calib, cfg, rng_obs)
    H = obs_matrix(cfg)
    R = obs_cov(cfg)
    obs_nodes = obs_node_indices(cfg)
    all_nodes = np.arange(N_NODES)
    unobs_nodes = np.array([i for i in all_nodes if i not in set(obs_nodes.tolist())])

    fb = cfg["filter"]

    # 2. でたらめな初期アンサンブル
    Z = init_ensemble(cfg, rng_filt)

    # 履歴コンテナ
    hist = {
        "node_names": NODE_NAMES,
        "obs_nodes": obs_nodes.tolist(),
        "unobs_nodes": unobs_nodes.tolist(),
        "q_scale_true": calib["q_scale_true"],
        "h_true": calib["h_true"],
        "times": [0.0],
        "truth_T": [truth[0].copy()],
        "ana_mean_T": [Z[:, :N_NODES].mean(axis=0)],
        "ana_std_T": [Z[:, :N_NODES].std(axis=0)],
        "fore_mean_T": [Z[:, :N_NODES].mean(axis=0)],
        "q_mean": [Z[:, I_QSCALE].mean()],
        "q_std": [Z[:, I_QSCALE].std()],
        "h_mean": [Z[:, I_H].mean()],
        "h_std": [Z[:, I_H].std()],
        "rmse_all": [_rmse(Z[:, :N_NODES].mean(axis=0), truth[0], all_nodes)],
        "rmse_obs": [_rmse(Z[:, :N_NODES].mean(axis=0), truth[0], obs_nodes)],
        "rmse_unobs": [_rmse(Z[:, :N_NODES].mean(axis=0), truth[0], unobs_nodes)],
        "ess": [float(cfg["ensemble"]["n_members"])],
        "obs_times": obs["step_times"].tolist(),
        "obs_values_C": (obs["values"] - 273.15).tolist(),
    }

    t_prev = 0.0
    for c, t_c in enumerate(obs["step_times"]):
        # --- 予報(forecast) ---
        Z = forecast(Z, t_prev, t_c, calib, cfg, rng_filt)
        fore_mean = Z[:, :N_NODES].mean(axis=0)

        # --- 解析(analysis) ---
        y = obs["values"][c]
        if filter_kind == "enkf":
            Z = enkf_mod.enkf_update(Z, y, H, R, rng_filt, inflation=fb["inflation"])
            ess = np.nan
        elif filter_kind == "pf":
            Z, _w, ess, _res = pf_mod.pf_update(
                Z, y, H, R, rng_filt,
                ess_frac=fb["pf_ess_frac"],
                jitter_T=fb["pf_jitter_T_C"],
                n_nodes=N_NODES,
            )
        else:
            raise ValueError(filter_kind)
        clip_params(Z, cfg)

        ana_mean = Z[:, :N_NODES].mean(axis=0)
        truth_c = truth[obs["step_idx"][c]]

        hist["times"].append(float(t_c))
        hist["truth_T"].append(truth_c.copy())
        hist["ana_mean_T"].append(ana_mean)
        hist["ana_std_T"].append(Z[:, :N_NODES].std(axis=0))
        hist["fore_mean_T"].append(fore_mean)
        hist["q_mean"].append(float(Z[:, I_QSCALE].mean()))
        hist["q_std"].append(float(Z[:, I_QSCALE].std()))
        hist["h_mean"].append(float(Z[:, I_H].mean()))
        hist["h_std"].append(float(Z[:, I_H].std()))
        hist["rmse_all"].append(_rmse(ana_mean, truth_c, all_nodes))
        hist["rmse_obs"].append(_rmse(ana_mean, truth_c, obs_nodes))
        hist["rmse_unobs"].append(_rmse(ana_mean, truth_c, unobs_nodes))
        hist["ess"].append(float(ess))
        t_prev = t_c

    # numpy 化
    for k in ["truth_T", "ana_mean_T", "ana_std_T", "fore_mean_T"]:
        hist[k] = np.array(hist[k])
    for k in ["times", "q_mean", "q_std", "h_mean", "h_std",
              "rmse_all", "rmse_obs", "rmse_unobs", "ess"]:
        hist[k] = np.array(hist[k])
    return hist
