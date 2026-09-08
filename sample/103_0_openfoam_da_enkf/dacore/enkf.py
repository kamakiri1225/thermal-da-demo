"""確率的アンサンブルカルマンフィルタ(stochastic EnKF, 摂動観測版)の解析更新.

理論:
    予報アンサンブル Z_f (n_ens, n_aug) から標本共分散を作り、
    カルマンゲイン K = C_zy (C_yy + R)^{-1} で観測へ引き寄せる。
        C_zy = Cov(Z_f, Y_f),  Y_f = H Z_f   (予報観測)
        C_yy = Cov(Y_f, Y_f)
    各メンバーに独立な観測摂動 eps~N(0,R) を与えて更新する(摂動観測法):
        Z_a^{(i)} = Z_f^{(i)} + K ( y + eps^{(i)} - Y_f^{(i)} )
    これによりアンサンブルの解析共分散が理論値へ一致する。

線形観測(温度の抜き出し)なので H は行列で表せるが、非線形観測にも
そのまま拡張できるよう Y_f = H(Z_f) を明示的に受け取る形にしている。
"""

from __future__ import annotations

import numpy as np


def enkf_update(Zf, y, H, R, rng, inflation=1.0):
    """確率的 EnKF の解析ステップ.

    引数:
        Zf        : (n_ens, n_aug) 予報アンサンブル
        y         : (n_obs,) 観測ベクトル
        H         : (n_obs, n_aug) 観測演算子
        R         : (n_obs, n_obs) 観測誤差共分散
        inflation : 乗算的共分散インフレーション係数(>=1)
    戻り値:
        Za        : (n_ens, n_aug) 解析アンサンブル
    """
    n_ens = Zf.shape[0]

    # 乗算インフレーション: 平均まわりに広げてサンプル貧困化を緩和
    if inflation and inflation != 1.0:
        zbar = Zf.mean(axis=0, keepdims=True)
        Zf = zbar + inflation * (Zf - zbar)

    Yf = Zf @ H.T                       # (n_ens, n_obs) 予報観測
    zbar = Zf.mean(axis=0)
    ybar = Yf.mean(axis=0)
    dZ = Zf - zbar                      # (n_ens, n_aug) 状態アノマリ
    dY = Yf - ybar                      # (n_ens, n_obs) 観測アノマリ

    C_zy = dZ.T @ dY / (n_ens - 1)      # (n_aug, n_obs)
    C_yy = dY.T @ dY / (n_ens - 1)      # (n_obs, n_obs)

    S = C_yy + R
    # 摂動観測: 各メンバーに独立ノイズ
    eps = rng.multivariate_normal(np.zeros(len(y)), R, size=n_ens)
    innov = (y[None, :] + eps) - Yf     # (n_ens, n_obs)

    # K = C_zy S^{-1}; 更新 = innov @ (S^{-1} C_zy^T) を解いて求める
    gain_T = np.linalg.solve(S, C_zy.T)  # (n_obs, n_aug) = S^{-1} C_zy^T
    Za = Zf + innov @ gain_T
    return Za
