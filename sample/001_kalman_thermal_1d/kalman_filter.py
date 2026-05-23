"""
kalman_filter.py
================
線形離散時間カルマンフィルタ。

状態方程式:
    x_{k+1} = A x_k + B u_k + w_k,  w ~ N(0, Q)

観測方程式:
    y_k = H x_k + v_k,  v ~ N(0, R)
"""

from typing import Tuple

import numpy as np


class KalmanFilter:
    """
    線形離散時間カルマンフィルタ。

    Parameters
    ----------
    A : np.ndarray
        状態遷移行列。shape は (n, n)。
    B : np.ndarray
        入力行列。shape は (n, m)。
    H : np.ndarray
        観測行列。shape は (p, n)。
    Q : np.ndarray
        プロセスノイズ共分散。shape は (n, n)。
    R : np.ndarray
        観測ノイズ共分散。shape は (p, p)。
    x0 : np.ndarray
        初期状態推定値。shape は (n,)。
    P0 : np.ndarray
        初期誤差共分散。shape は (n, n)。
    """

    def __init__(
        self,
        A: np.ndarray,
        B: np.ndarray,
        H: np.ndarray,
        Q: np.ndarray,
        R: np.ndarray,
        x0: np.ndarray,
        P0: np.ndarray,
    ) -> None:
        self.A = A
        self.B = B
        self.H = H
        self.Q = Q
        self.R = R
        self.x = x0.copy()
        self.P = P0.copy()

    def predict(self, u: np.ndarray) -> None:
        """
        予測ステップ。

        x^- = A x + B u
        P^- = A P A^T + Q
        """
        self.x = self.A @ self.x + self.B @ u
        self.P = self.A @ self.P @ self.A.T + self.Q

    def update(self, y: np.ndarray) -> None:
        """
        更新ステップ。観測値 y を使って状態推定値を補正する。

        K = P^- H^T (H P^- H^T + R)^-1
        x = x^- + K (y - H x^-)
        P = (I - K H) P^-
        """
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)
        innov = y - self.H @ self.x
        self.x = self.x + K @ innov
        self.P = (np.eye(self.x.shape[0]) - K @ self.H) @ self.P

    def step(
        self,
        u: np.ndarray,
        y: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        予測と更新を1ステップ分まとめて実行する。

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            更新後の状態推定値と誤差共分散。
        """
        self.predict(u)
        self.update(y)
        return self.x.copy(), self.P.copy()

    @property
    def uncertainty(self) -> np.ndarray:
        """各状態の推定標準偏差を返す。"""
        return np.sqrt(np.diag(self.P))
