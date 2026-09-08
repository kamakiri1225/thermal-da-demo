"""中空円筒CHTケース(102_0)の集中定数(ROM)熱モデル.

102_0 の chtMultiRegionFoam は 67k セル・1メンバー約1時間かかるため、
アンサンブルを何十メンバーも回すデータ同化には重すぎる。そこで、
固体の温度場を少数の代表ノードに縮約した線形の集中定数モデル(ROM)を用意する。

ノード(いずれも 102_0 の固体プローブ位置に対応):
    0: hot   (+X, ヒータ直下側)      -> 102_0 の T_hot
    1: mid   (+Y, 周方向90度)        -> 102_0 の T_mid
    2: cold  (-X, 反ヒータ側)        -> 102_0 の T_cold
    3: top   (上面近傍)              -> 102_0 の T_top
    4: core  (肉厚中心の潜在ノード)   -> 直接は観測しない

各ノードのエネルギー収支(集中定数):
    C_i dT_i/dt = Σ_j K_ij (T_j - T_i) + q_i(t) - h_i (T_i - T_air)

    - K_ij : ノード間の熱コンダクタンス [W/K]
    - q_i  : 発熱入力 [W] (ヒータは hot ノードへ入れる)
    - h_i  : 周囲空気への実効放熱係数 [W/K]
    - T_air: 周囲空気温度 [K] (一定と仮定)

これは dT/dt = M T + f(t) の線形システムなので、
アンサンブル全体をベクトル化して高速に前進積分できる。

データ同化で未知数として推定する物理パラメータ:
    - q_scale : ヒータ発熱の倍率 (真値 1.0, 定格 15 W に掛かる)
    - h       : 放熱係数 [W/K]
それ以外(コンダクタンス K・熱容量 C)は calibrate.py で
102_0 の temperature_history.csv に合わせて同定し、既知として固定する。
"""

from __future__ import annotations

import numpy as np

NODE_NAMES = ["hot", "mid", "cold", "top", "core"]
N_NODES = len(NODE_NAMES)
HOT, MID, COLD, TOP, CORE = range(N_NODES)

# 102_0 の設定と一致させる基準値
T_AIR_K = 293.15          # 周囲空気温度(=初期温度 20 degC)
HEATER_RATED_W = 15.0     # ヒータ定格 15 W
HEATER_ON_S = 300.0       # 0-300 s 通電、以降 0 W


def heater_power_W(t, q_scale=1.0):
    """時刻 t [s] のヒータ発熱 [W]. 0-300 s は q_scale*15 W, 以降 0 W."""
    on = np.asarray(t) < HEATER_ON_S
    return np.where(on, q_scale * HEATER_RATED_W, 0.0)


class ROMParams:
    """ROM の物理パラメータ.

    calibrated (既知・固定): 熱容量 C[5], コンダクタンス k_circ/k_core/k_axial
    uncertain (同化で推定): q_scale, h
    """

    def __init__(self, C, k_circ, k_core, k_axial, h, q_scale=1.0):
        self.C = np.asarray(C, dtype=float)      # [W*s/K] 各ノード熱容量 (5,)
        self.k_circ = float(k_circ)              # hot-mid, mid-cold の周方向コンダクタンス
        self.k_core = float(k_core)              # 各ノード <-> core の径方向コンダクタンス
        self.k_axial = float(k_axial)            # hot-top, top-core の軸方向コンダクタンス
        self.h = float(h)                        # 放熱係数 [W/K] (全ノード共通)
        self.q_scale = float(q_scale)

    # --- 校正パラメータのベクトル化(least_squares 用) -----------------
    @staticmethod
    def calib_vector_to_params(v, h, q_scale=1.0):
        """校正ベクトル v=[C0..C4, k_circ, k_core, k_axial] から ROMParams を作る."""
        C = v[0:5]
        return ROMParams(C, v[5], v[6], v[7], h, q_scale)

    def as_calib_vector(self):
        return np.concatenate([self.C, [self.k_circ, self.k_core, self.k_axial]])


def conductance_matrix(k_circ, k_core, k_axial):
    """ノード間コンダクタンス K_ij [W/K] (対称, 5x5) を組む."""
    K = np.zeros((N_NODES, N_NODES))

    def link(a, b, k):
        K[a, b] = k
        K[b, a] = k

    # 周方向: hot -> mid -> cold
    link(HOT, MID, k_circ)
    link(MID, COLD, k_circ)
    # 径方向: 各表面ノード -> core(肉厚中心)
    link(HOT, CORE, k_core)
    link(MID, CORE, k_core)
    link(COLD, CORE, k_core)
    link(TOP, CORE, k_core)
    # 軸方向: hot -> top, top -> core
    link(HOT, TOP, k_axial)
    return K


def system_matrix(params: ROMParams):
    """線形システム dT/dt = M T + b_const + b_heater(t)/C の M と定数項を返す.

    戻り値:
        M       : (5,5) システム行列 (熱容量で正規化済み)
        b_air   : (5,) 放熱による定数項 (h*T_air/C_i)
    ヒータ項は時刻依存なので別途 heater_input() で加える。
    """
    K = conductance_matrix(params.k_circ, params.k_core, params.k_axial)
    C = params.C
    M = np.zeros((N_NODES, N_NODES))
    for i in range(N_NODES):
        # 隣接ノードとの伝導
        M[i, :] = K[i, :] / C[i]
        M[i, i] = -(K[i, :].sum() + params.h) / C[i]
    b_air = params.h * T_AIR_K / C
    return M, b_air


def heater_input(t, params: ROMParams):
    """ヒータ発熱による項 q/C を hot ノードにだけ与える (5,)."""
    b = np.zeros(N_NODES)
    b[HOT] = heater_power_W(t, params.q_scale) / params.C[HOT]
    return b


def integrate_single(T0, params: ROMParams, t0, t1, dt):
    """1メンバーを RK4 で t0->t1 積分し、時刻歴 (times, T[nt,5]) を返す."""
    n = max(1, int(round((t1 - t0) / dt)))
    dt = (t1 - t0) / n
    M, b_air = system_matrix(params)
    T = np.asarray(T0, dtype=float).copy()
    times = [t0]
    traj = [T.copy()]

    def deriv(Tv, t):
        return M @ Tv + b_air + heater_input(t, params)

    t = t0
    for _ in range(n):
        k1 = deriv(T, t)
        k2 = deriv(T + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = deriv(T + 0.5 * dt * k2, t + 0.5 * dt)
        k4 = deriv(T + dt * k3, t + dt)
        T = T + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
        times.append(t)
        traj.append(T.copy())
    return np.array(times), np.array(traj)


def integrate_ensemble(T, C, k_circ, k_core, k_axial, h, q_scale, t0, t1, dt):
    """アンサンブルをベクトル化して t0->t1 積分する.

    引数:
        T        : (n_ens, 5) 各メンバーのノード温度 [K]
        C        : (5,) 共通熱容量(校正済み・固定)
        k_*      : スカラー(校正済み・固定コンダクタンス)
        h        : (n_ens,) メンバーごとの放熱係数 [W/K]
        q_scale  : (n_ens,) メンバーごとのヒータ倍率
    戻り値:
        T        : (n_ens, 5) t1 での温度
    メンバーごとに h と q_scale が異なるため、M をメンバー軸でスタックして
    einsum で一括計算する。
    """
    T = np.asarray(T, dtype=float).copy()
    n_ens = T.shape[0]
    C = np.asarray(C, dtype=float)
    h = np.asarray(h, dtype=float)
    q_scale = np.asarray(q_scale, dtype=float)

    K = conductance_matrix(k_circ, k_core, k_axial)     # (5,5) 固定
    Ksum = K.sum(axis=1)                                # (5,)

    # メンバーごとのシステム行列 M_e (n_ens,5,5)
    M = np.broadcast_to(K, (n_ens, N_NODES, N_NODES)) / C[None, :, None]
    M = M.copy()
    diag = -(Ksum[None, :] + h[:, None]) / C[None, :]   # (n_ens,5)
    idx = np.arange(N_NODES)
    M[:, idx, idx] = diag
    b_air = h[:, None] * T_AIR_K / C[None, :]           # (n_ens,5)

    n = max(1, int(round((t1 - t0) / dt)))
    dt = (t1 - t0) / n

    def deriv(Tv, t):
        heat = np.zeros_like(Tv)
        heat[:, HOT] = heater_power_W(t, q_scale) / C[HOT]
        return np.einsum("eij,ej->ei", M, Tv) + b_air + heat

    t = t0
    for _ in range(n):
        k1 = deriv(T, t)
        k2 = deriv(T + 0.5 * dt * k1, t + 0.5 * dt)
        k3 = deriv(T + 0.5 * dt * k2, t + 0.5 * dt)
        k4 = deriv(T + dt * k3, t + dt)
        T = T + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt
    return T
