"""最適内挿（Optimal Interpolation, OI）の解析更新.

EnKF との違い：背景誤差共分散 B を**事前に決め打ちで固定**する（アンサンブルで
毎サイクル作り直さない）。だから予報は「1本の背景トラジェクトリ」だけでよく、
アンサンブル N 本を回す必要がない＝計算が軽い。

    z_a = z_b + K ( y - h(z_b) ),   K = B Hᵀ ( H B Hᵀ + R )⁻¹

B は固定なのでゲイン K も（H,R が同じなら）毎サイクル同じ。
本実装では観測演算子 H は線形（温度=状態抜き出し、変位=温度場のPODモード応答）とし、
予報観測の残差 innov = y - h(z_b) は非線形評価 h(z_b) を渡してもよい。
"""
from __future__ import annotations
import numpy as np


def oi_gain(H, B, R):
    """固定ゲイン K = B Hᵀ (H B Hᵀ + R)⁻¹ を返す（B,H,R が固定なら1回だけ計算）."""
    HB = H @ B                      # (n_obs, n)
    S = HB @ H.T + R                # (n_obs, n_obs)
    K = np.linalg.solve(S, HB).T    # (n, n_obs) = B Hᵀ S⁻¹
    return K


def oi_update(z_b, y, y_pred, K):
    """OI解析：z_a = z_b + K ( y - y_pred ). y_pred=h(z_b)（線形なら H z_b）."""
    return z_b + K @ (y - y_pred)
