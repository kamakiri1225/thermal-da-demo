"""拡大状態(augmented state)アンサンブルの生成と前進(forecast).

拡大状態ベクトル z (長さ 7):
    z[0:5] = 5ノードの温度 [K]  (hot, mid, cold, top, core)
    z[5]   = q_scale  (ヒータ発熱倍率, 真値 1.0)
    z[6]   = h        (放熱係数 [W/K], 真値 ~0.024)

温度(状態)と物理パラメータを同じベクトルに詰めることで、
EnKF/PF は観測とパラメータの相関を通じて未知パラメータも同時推定できる
(パラメータ推定 = 拡大状態フィルタリング)。
"""

from __future__ import annotations

import numpy as np

from .cht_rom import N_NODES, integrate_ensemble

I_QSCALE = N_NODES        # =5
I_H = N_NODES + 1         # =6
N_AUG = N_NODES + 2       # =7


def init_ensemble(cfg, rng):
    """「でたらめな状態」からアンサンブルを初期化して Z (n_ens, 7) を返す."""
    ec = cfg["ensemble"]
    n = ec["n_members"]
    Z = np.empty((n, N_AUG))
    # 温度: 各ノード・各メンバーを一様乱数で散らす(真の 20 degC を知らない)
    Tmin = ec["init_T_min_C"] + 273.15
    Tmax = ec["init_T_max_C"] + 273.15
    Z[:, :N_NODES] = rng.uniform(Tmin, Tmax, size=(n, N_NODES))
    # パラメータ: 真値からずらした正規分布
    Z[:, I_QSCALE] = rng.normal(ec["init_q_scale_mean"], ec["init_q_scale_std"], n)
    Z[:, I_H] = rng.normal(ec["init_h_mean"], ec["init_h_std"], n)
    clip_params(Z, cfg)
    return Z


def clip_params(Z, cfg):
    """パラメータを物理的な範囲へクリップ(in-place)."""
    fb = cfg["filter"]
    Z[:, I_QSCALE] = np.clip(Z[:, I_QSCALE], *fb["q_scale_bounds"])
    Z[:, I_H] = np.clip(Z[:, I_H], *fb["h_bounds"])
    return Z


def forecast(Z, t0, t1, calib, cfg, rng):
    """アンサンブルを t0->t1 まで ROM で前進積分して返す.

    パラメータ(q_scale, h)には微小なランダムウォーク(jitter)を加える。
    これは静的パラメータでアンサンブルが1点に潰れる(サンプル貧困化)のを防ぐ、
    パラメータ推定における定石的な処理。
    """
    Z = Z.copy()
    fb = cfg["filter"]
    # パラメータのランダムウォーク
    Z[:, I_QSCALE] += rng.normal(0.0, fb["param_jitter_q"], Z.shape[0])
    Z[:, I_H] += rng.normal(0.0, fb["param_jitter_h"], Z.shape[0])
    clip_params(Z, cfg)

    T = Z[:, :N_NODES]
    dt = cfg["experiment"]["model_dt_s"]
    Z[:, :N_NODES] = integrate_ensemble(
        T,
        C=calib["C"],
        k_circ=calib["k_circ"],
        k_core=calib["k_core"],
        k_axial=calib["k_axial"],
        h=Z[:, I_H],
        q_scale=Z[:, I_QSCALE],
        t0=t0,
        t1=t1,
        dt=dt,
    )
    return Z
