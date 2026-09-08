"""観測(数点の温度)の定義と、双子実験の真値・合成観測の生成.

双子実験(OSSE)の流れ:
    1. 校正済み ROM に真パラメータ(q_scale=1, h=h_true)と真の初期温度を与え、
       「真値ラン」の温度時刻歴を作る。
    2. 真値ランのうち観測ノード(例: hot, cold)だけを obs_interval ごとに抜き出し、
       正規ノイズを乗せて「合成観測」とする。
    3. アンサンブルはこの数点の観測だけを頼りに、でたらめな初期状態から
       真値の温度場へ収束できるか(=データ同化の効果)を見る。
"""

from __future__ import annotations

import numpy as np

from .cht_rom import NODE_NAMES, N_NODES, integrate_single, ROMParams

NAME_TO_IDX = {name: i for i, name in enumerate(NODE_NAMES)}


def obs_node_indices(cfg):
    """設定の観測ノード名リストをインデックス配列に変換する."""
    names = cfg["observation"]["nodes"]
    return np.array([NAME_TO_IDX[n] for n in names], dtype=int)


def obs_matrix(cfg):
    """観測演算子 H (n_obs, 7): 拡大状態から観測ノード温度を抜き出す."""
    from .ensemble import N_AUG

    idx = obs_node_indices(cfg)
    H = np.zeros((len(idx), N_AUG))
    for r, c in enumerate(idx):
        H[r, c] = 1.0
    return H


def obs_cov(cfg):
    """観測誤差共分散 R (n_obs, n_obs), 対角 = noise_C^2 [K^2]."""
    idx = obs_node_indices(cfg)
    var = cfg["observation"]["noise_C"] ** 2
    return np.eye(len(idx)) * var


def true_params(calib):
    """校正結果から真値パラメータの ROMParams を組む."""
    return ROMParams(
        C=calib["C"],
        k_circ=calib["k_circ"],
        k_core=calib["k_core"],
        k_axial=calib["k_axial"],
        h=calib["h_true"],
        q_scale=calib["q_scale_true"],
    )


def generate_truth(calib, cfg):
    """真値ランを積分し、(times, truth_traj[nt,5]) を返す [K]."""
    p = true_params(calib)
    T0 = np.full(N_NODES, cfg["truth"]["T0_C"] + 273.15)
    t_end = cfg["experiment"]["t_end_s"]
    dt = cfg["experiment"]["model_dt_s"]
    times, traj = integrate_single(T0, p, 0.0, t_end, dt=dt)
    return times, traj


def make_observations(times, truth_traj, calib, cfg, rng):
    """真値時刻歴から合成観測を作る.

    戻り値 obs: dict
        step_times : (n_cyc,) 観測時刻 [s]
        step_idx   : (n_cyc,) truth_traj 上のインデックス
        values     : (n_cyc, n_obs) 観測値 [K] (ノイズ入り)
        nodes      : 観測ノードのインデックス
    """
    idx_nodes = obs_node_indices(cfg)
    t_end = cfg["experiment"]["t_end_s"]
    dt_obs = cfg["experiment"]["obs_interval_s"]
    noise = cfg["observation"]["noise_C"]

    step_times = np.arange(dt_obs, t_end + 1e-9, dt_obs)
    step_idx = np.searchsorted(times, step_times)
    step_idx = np.clip(step_idx, 0, len(times) - 1)

    clean = truth_traj[step_idx][:, idx_nodes]
    values = clean + rng.normal(0.0, noise, size=clean.shape)
    return {
        "step_times": step_times,
        "step_idx": step_idx,
        "values": values,
        "clean": clean,
        "nodes": idx_nodes,
    }
