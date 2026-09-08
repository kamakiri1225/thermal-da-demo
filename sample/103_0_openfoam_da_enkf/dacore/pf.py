"""ブートストラップ粒子フィルタ(SIR, 系統リサンプリング + 正則化ジッタ).

理論:
    各粒子 z^{(i)} に尤度で重みを付ける:
        w^{(i)} ∝ exp(-1/2 (y - H z^{(i)})^T R^{-1} (y - H z^{(i)}))
    重みを正規化し、有効サンプル数 ESS = 1/Σ w^2 が閾値を下回ったら
    重みに比例して粒子を系統リサンプリングする。
    リサンプリング後は同一粒子が重複するため、温度に微小ジッタを加える
    (正則化粒子フィルタ)。静的パラメータの貧困化は forecast 側の
    ランダムウォークで緩和する。

EnKF と違いガウス性を仮定しないため、非ガウス・多峰な事後分布も表現できるが、
状態次元が上がると必要粒子数が急増する(次元の呪い)。本ケースは低次元なので
粒子フィルタが安定に動く好例になっている。
"""

from __future__ import annotations

import numpy as np


def _log_likelihood(Zf, y, H, R):
    Yf = Zf @ H.T                        # (n_ens, n_obs)
    innov = y[None, :] - Yf              # (n_ens, n_obs)
    Rinv = np.linalg.inv(R)
    d2 = np.einsum("ei,ij,ej->e", innov, Rinv, innov)
    return -0.5 * d2


def systematic_resample(weights, rng):
    """系統リサンプリング: 重みに比例した親インデックスを返す."""
    n = len(weights)
    positions = (rng.random() + np.arange(n)) / n
    cumsum = np.cumsum(weights)
    cumsum[-1] = 1.0
    idx = np.searchsorted(cumsum, positions)
    return np.clip(idx, 0, n - 1)


def pf_update(Zf, y, H, R, rng, ess_frac=0.5, jitter_T=0.0, n_nodes=5):
    """粒子フィルタの解析ステップ.

    戻り値:
        Za        : (n_ens, n_aug) 解析(リサンプリング後)アンサンブル
        weights   : (n_ens,) 更新後の重み(リサンプリング時は一様)
        ess       : 有効サンプル数(リサンプリング判定に用いた値)
        resampled : bool
    """
    n_ens = Zf.shape[0]
    logw = _log_likelihood(Zf, y, H, R)
    logw -= logw.max()                   # オーバーフロー防止
    w = np.exp(logw)
    w /= w.sum()

    ess = 1.0 / np.sum(w ** 2)
    resampled = False
    Za = Zf.copy()
    if ess < ess_frac * n_ens:
        idx = systematic_resample(w, rng)
        Za = Zf[idx].copy()
        if jitter_T > 0.0:
            Za[:, :n_nodes] += rng.normal(0.0, jitter_T, size=(n_ens, n_nodes))
        w = np.full(n_ens, 1.0 / n_ens)
        resampled = True
    return Za, w, ess, resampled
