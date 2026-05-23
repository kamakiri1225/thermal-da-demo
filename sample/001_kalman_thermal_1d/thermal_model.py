"""
thermal_model.py
================
1次元熱伝導棒の集中定数モデル。

10個の節点で棒の温度分布を近似し、節点0に熱入力を与える。
両端では周囲温度基準の対流放熱を簡易的に表現する。
状態量は周囲温度からの温度上昇ではなく、温度そのものとして扱う。
"""

import numpy as np


class ThermalBarModel:
    """
    1次元熱伝導棒の離散時間モデル。

    Parameters
    ----------
    n_nodes : int
        節点数。
    dt : float
        時間刻み [s]。
    alpha : float
        熱拡散率 [m^2/s]。鋼材の代表値として 1.2e-5 程度を使う。
    dx : float
        節点間距離 [m]。
    h_conv : float
        対流熱伝達率 [W/(m^2 K)]。両端境界の放熱に使う。
    rho_cp : float
        体積熱容量 [J/(m^3 K)]。鋼材の代表値として 3.5e6 程度を使う。
    """

    def __init__(
        self,
        n_nodes: int = 10,
        dt: float = 10.0,
        alpha: float = 1.2e-5,
        dx: float = 0.1,
        h_conv: float = 8.0,
        rho_cp: float = 3.5e6,
    ) -> None:
        self.n = n_nodes
        self.dt = dt
        self.alpha = alpha
        self.dx = dx
        self.h_conv = h_conv
        self.rho_cp = rho_cp

        self.A_c, self.B_c = self._build_continuous()
        self.A_d, self.B_d = self._discretize()

    def _build_continuous(self) -> tuple[np.ndarray, np.ndarray]:
        """連続時間系の状態行列 A_c と入力行列 B_c を作る。"""
        n = self.n
        r = self.alpha / self.dx**2
        h = self.h_conv / (self.rho_cp * self.dx)

        # 内部節点は中央差分、端点は片側差分と対流放熱で近似する。
        A = np.zeros((n, n))
        for i in range(1, n - 1):
            A[i, i - 1] = r
            A[i, i] = -2.0 * r
            A[i, i + 1] = r

        A[0, 0] = -r - h
        A[0, 1] = r
        A[-1, -1] = -r - h
        A[-1, -2] = r

        # 熱入力は節点0に集中して入るものとして扱う。
        B = np.zeros((n, 1))
        B[0, 0] = 1.0 / (self.rho_cp * self.dx)

        return A, B

    def _discretize(self) -> tuple[np.ndarray, np.ndarray]:
        """
        前進オイラー法で連続時間系を離散化する。

        x_{k+1} = (I + A_c dt) x_k + (B_c dt) u_k
        """
        A_d = np.eye(self.n) + self.A_c * self.dt
        B_d = self.B_c * self.dt
        return A_d, B_d

    def step(
        self,
        x: np.ndarray,
        u: np.ndarray,
        noise_std: float = 0.0,
    ) -> np.ndarray:
        """
        1ステップだけ温度状態を進める。

        Parameters
        ----------
        x : np.ndarray
            現在の温度状態。shape は (n,)。
        u : np.ndarray
            熱入力 [W/m]。shape は (1,)。
        noise_std : float
            プロセスノイズの標準偏差 [degC]。
        """
        x_next = self.A_d @ x + self.B_d @ u
        if noise_std > 0.0:
            x_next = x_next + np.random.randn(self.n) * noise_std
        return x_next
