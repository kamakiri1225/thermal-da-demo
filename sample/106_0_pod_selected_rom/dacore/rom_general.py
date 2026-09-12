"""任意配置の代表点に対する一般化ROM（集中定数熱モデル）.

105までのROMは hot/mid/cold/top/core の固定トポロジだったが、
106では POD＋Q-DEIM で選んだ任意の5点をノードにするため、
ノード間コンダクタンスを「全点対の一般グラフ」として持つ.

各ノードのエネルギー収支（線形）:
    C_i dT_i/dt = Σ_j K_ij (T_j - T_i) + q_i(t) - h (T_i - T_air)
    - K_ij : ノード間コンダクタンス [W/K]（対称・全点対）
    - q_i  : ヒータ発熱（ヒータ最近傍ノードに投入）
    - h    : 放熱係数 [W/K]（全ノード共通）

未知として同化で推定するのは q_scale, h。C と K は calibrate で固定する.
"""
from __future__ import annotations
import numpy as np

T_AIR_K = 293.15
HEATER_RATED_W = 15.0
HEATER_ON_S = 300.0


def heater_power_W(t, q_scale=1.0):
    return np.where(np.asarray(t) < HEATER_ON_S, q_scale * HEATER_RATED_W, 0.0)


def tri_to_matrix(k_upper, n):
    """上三角(i<j)のコンダクタンス列 -> 対称行列 K(n,n)."""
    K = np.zeros((n, n)); idx = 0
    for i in range(n):
        for j in range(i + 1, n):
            K[i, j] = K[j, i] = k_upper[idx]; idx += 1
    return K


def n_edges(n):
    return n * (n - 1) // 2


def integrate_ensemble(T, C, Kmat, h, q_scale, heat_node, t0, t1, dt):
    """アンサンブル一括前進（RK4）. T:(n_ens,n), h/q_scale:(n_ens,).

    Kmat:(n,n) 固定コンダクタンス, heat_node:ヒータ投入ノード index.
    """
    T = np.asarray(T, float).copy(); n_ens, n = T.shape
    C = np.asarray(C, float); h = np.asarray(h, float); q_scale = np.asarray(q_scale, float)
    Ksum = Kmat.sum(axis=1)
    M = np.broadcast_to(Kmat, (n_ens, n, n)) / C[None, :, None]
    M = M.copy()
    diag = -(Ksum[None, :] + h[:, None]) / C[None, :]
    idx = np.arange(n); M[:, idx, idx] = diag
    b_air = h[:, None] * T_AIR_K / C[None, :]
    steps = max(1, int(round((t1 - t0) / dt))); dt = (t1 - t0) / steps

    def deriv(Tv, t):
        heat = np.zeros_like(Tv)
        heat[:, heat_node] = heater_power_W(t, q_scale) / C[heat_node]
        return np.einsum("eij,ej->ei", M, Tv) + b_air + heat

    t = t0
    for _ in range(steps):
        k1 = deriv(T, t); k2 = deriv(T + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = deriv(T + 0.5 * dt * k2, t + 0.5 * dt); k4 = deriv(T + dt * k3, t + dt)
        T = T + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
    return T


def integrate_single(T0, C, Kmat, h, q_scale, heat_node, t0, t1, dt):
    """1本の時刻歴 (times, T[nt,n]) を返す（校正・検証用）."""
    T = np.asarray(T0, float).reshape(1, -1)
    steps = max(1, int(round((t1 - t0) / dt)))
    times = [t0]; traj = [T[0].copy()]
    tt = t0; sub = (t1 - t0) / steps
    for _ in range(steps):
        T = integrate_ensemble(T, C, Kmat, np.array([h]), np.array([q_scale]),
                               heat_node, tt, tt + sub, sub)
        tt += sub; times.append(tt); traj.append(T[0].copy())
    return np.array(times), np.array(traj)
