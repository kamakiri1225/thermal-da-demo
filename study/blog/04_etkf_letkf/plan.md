# 記事04 計画：ETKF・LETKF（拡張アンサンブルカルマンフィルタ）

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | ETKF・LETKFとは？大規模データ同化を支えるアンサンブル変換をPythonで解説 |
| **メタディスクリプション** | ETKF（Ensemble Transform KF）とLETKF（Local ETKF）は、観測摂動なしで安定したアンサンブル更新を実現する拡張EnKF手法です。Lorenz 96モデルへの適用をPythonで実装します。 |
| **メインキーワード** | ETKF LETKF アンサンブルカルマンフィルタ |
| **関連キーワード** | 拡張アンサンブルカルマンフィルタ, LETKF Python, アンサンブル変換, 局所化 データ同化, Lorenz 96 |
| **想定文字数** | 5,000〜6,500字 |
| **Python可視化** | あり（Lorenz 96・EnKF vs ETKF比較・局所化効果） |

---

## EnKF（Stochastic）の問題点 → ETKFの動機

### EnKF（記事03）のおさらいと課題

| 課題 | 内容 |
|---|---|
| **サンプリング誤差** | 観測摂動 $\boldsymbol{\epsilon}^{(i)} \sim \mathcal{N}(0,\mathbf{R})$ を加えるため、$N$ が小さいと解析スプレッドが不正確 |
| **フィルター発散** | スプレッドが過小評価されると観測を無視し始める |
| **行列逆演算コスト** | $(H\mathbf{P}H^T + \mathbf{R})^{-1}$ は観測次元 $m$ が大きいと重い |

**ETKFのアイデア**：観測摂動を使わず、**アンサンブル変換行列** $\mathbf{T}$ を使って解析アンサンブルを直接構成する。

---

## 記事構成

### リード文
```
記事03で解説したEnKF（Stochastic）は観測に確率的な摂動を加えるため、
アンサンブルサイズが小さいと誤差が入ります。
ETKF（Ensemble Transform Kalman Filter）はこの問題を
「アンサンブル変換行列」という発想で解決した、
より安定した拡張アンサンブルカルマンフィルタです。
```

### この記事でわかること
- Stochastic EnKFとETKFの違いと、なぜETKFが安定かがわかる
- アンサンブル変換行列 $\mathbf{T}$ の導出と意味が理解できる
- LETKFの局所化が大規模問題でなぜ有効かがわかり、Pythonで実装できる

---

### H2-1：ETKF（Ensemble Transform Kalman Filter）

**キーワード**：ETKF 仕組み, アンサンブル変換行列

#### H3：アンサンブル空間での更新

アンサンブル偏差行列 $\mathbf{X}^f \in \mathbb{R}^{n \times N}$：
$$\mathbf{X}^f = \frac{1}{\sqrt{N-1}}[\mathbf{x}^{f,(1)} - \bar{\mathbf{x}}^f, \ldots, \mathbf{x}^{f,(N)} - \bar{\mathbf{x}}^f]$$

解析アンサンブル偏差行列：
$$\mathbf{X}^a = \mathbf{X}^f \mathbf{T}$$

変換行列：
$$\mathbf{T} = \left[(N-1)\mathbf{I} + (\mathbf{Y}^f)^T \mathbf{R}^{-1} \mathbf{Y}^f\right]^{-\frac{1}{2}}$$

ここで $\mathbf{Y}^f = H\mathbf{X}^f$（観測空間のアンサンブル偏差）。

| 量 | サイズ | 意味 |
|---|---|---|
| $\mathbf{X}^f$ | $n \times N$ | 状態空間のアンサンブル偏差 |
| $\mathbf{Y}^f$ | $m \times N$ | 観測空間のアンサンブル偏差 |
| $\mathbf{T}$ | $N \times N$ | アンサンブル変換行列（$N \ll n$ なので安い） |

**ETKFの核心**：$n \times n$ 行列を扱わず、$N \times N$ 行列だけで全更新が完結する。

#### H3：解析平均の更新

$$\bar{\mathbf{x}}^a = \bar{\mathbf{x}}^f + \mathbf{X}^f \mathbf{w}^a$$

$$\mathbf{w}^a = \mathbf{\tilde{P}}^a (\mathbf{Y}^f)^T \mathbf{R}^{-1} (\mathbf{y} - H\bar{\mathbf{x}}^f)$$

$$\mathbf{\tilde{P}}^a = \left[(N-1)\mathbf{I} + (\mathbf{Y}^f)^T \mathbf{R}^{-1} \mathbf{Y}^f\right]^{-1}$$

---

### H2-2：LETKF（Local Ensemble Transform Kalman Filter）

**キーワード**：LETKF 局所化, LETKF Python, データ同化 大規模

#### H3：なぜ局所化が必要か

- $N = 50$ のアンサンブルで $n = 10^6$ の共分散を推定 → 偽の遠距離相関が発生
- 物理的に離れた点の間に相関はないはず → **虚偽相関（spurious correlation）**

#### H3：局所化の仕組み

各グリッド点 $i$ ごとに、**周囲の観測だけ**を使って独立にETKFを実行：

$$\mathbf{R}^{loc} = \mathbf{R} \oslash \rho(\mathbf{r})$$

- $\rho(\mathbf{r})$：距離 $\mathbf{r}$ に応じた減衰関数（GC関数など）
- 遠い観測の影響をゼロに近づける

**LETKFの計算の流れ**：
1. 各グリッド点 $i$ ごとにローカルウィンドウ内の観測を選択
2. ローカルETKFを実行（$N \times N$ 行列演算）
3. 各点の解析値を集めて全体の解析場を構成

→ 並列化が容易（各グリッド点は独立に計算できる）

---

### H2-3：テストモデル — Lorenz 96

**キーワード**：Lorenz 96 データ同化, LETKF Lorenz 96

$$\frac{dx_j}{dt} = (x_{j+1} - x_{j-2})x_{j-1} - x_j + F, \quad j = 1, \ldots, J$$

- $J = 40$、$F = 8$（カオス的挙動）
- 周期境界条件：$x_{J+1} = x_1$
- すべての点で5点おきに観測（$m = 8$）

**なぜLorenz 96か**：
- 40次元で「大規模」の挙動をシミュレートできる
- 局所的な伝播構造を持つ → 局所化の効果が明確に現れる

---

### H2-4：PythonによるLETKF実装

**キーワード**：LETKF Python 実装, アンサンブルカルマンフィルタ 拡張

```python
# letkf_lorenz96.py
import numpy as np
import matplotlib.pyplot as plt

# ===== Lorenz 96 モデル =====
def lorenz96(x, F=8.0):
    J = len(x)
    dx = np.zeros(J)
    for j in range(J):
        dx[j] = (x[(j+1) % J] - x[(j-2) % J]) * x[(j-1) % J] - x[j] + F
    return dx

def rk4_step(x, dt, model):
    k1 = model(x)
    k2 = model(x + 0.5*dt*k1)
    k3 = model(x + 0.5*dt*k2)
    k4 = model(x + dt*k3)
    return x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

# ===== LETKF =====
def letkf_update(ensemble, y_obs, obs_idx, R_diag, rho_loc, inflation=1.05):
    """
    Parameters
    ----------
    ensemble  : (n, N) アンサンブル行列
    y_obs     : (m,) 観測値
    obs_idx   : (m,) 観測点のグリッドインデックス
    R_diag    : (m,) 観測誤差分散（対角成分）
    rho_loc   : 局所化半径（グリッド数）
    inflation : アンサンブルインフレーション係数
    """
    n, N = ensemble.shape
    m = len(y_obs)
    xa_mean = np.zeros(n)
    Xa = np.zeros((n, N))

    # アンサンブル平均と偏差
    xf_mean = ensemble.mean(axis=1)
    Xf = (ensemble - xf_mean[:, None]) / np.sqrt(N - 1)

    for i in range(n):
        # --- 局所化：グリッド点 i に近い観測を選択 ---
        dist = np.array([min(abs(i - j), n - abs(i - j)) for j in obs_idx])
        loc_weight = np.maximum(0, 1 - (dist / rho_loc) ** 2) ** 2  # 簡易GC関数
        active = loc_weight > 0
        if active.sum() == 0:
            xa_mean[i] = xf_mean[i]
            Xa[i, :] = Xf[i, :]
            continue

        # ローカル観測量
        y_loc = y_obs[active]
        R_loc = np.diag(R_diag[active] / loc_weight[active])
        Yf_loc = Xf[i:i+1, :] * 0  # ダミー（本来はH_locを適用）
        # H は単純に obs_idx の位置から値を取るとして実装
        Yf_loc = Xf[obs_idx[active], :]  # (m_loc, N)

        # ETKFの中心計算（N×N の逆行列）
        C = Yf_loc.T @ np.linalg.solve(R_loc, Yf_loc)  # (N, N)
        Pa_tilde = np.linalg.inv((N - 1) * np.eye(N) + C)

        # 平均更新
        innov = y_loc - (xf_mean[obs_idx[active]] + Yf_loc @ np.zeros(N))
        w = Pa_tilde @ Yf_loc.T @ np.linalg.solve(R_loc, innov)
        xa_mean[i] = xf_mean[i] + inflation * Xf[i, :] @ w

        # 偏差更新（変換行列 T = sqrt(Pa_tilde) * sqrt(N-1)）
        T = np.real(np.linalg.matrix_power(
            np.linalg.cholesky((N-1) * Pa_tilde), 1))
        Xa[i, :] = inflation * Xf[i, :] @ T

    return xa_mean[:, None] + np.sqrt(N - 1) * Xa
```

**可視化パターン**
```
1. 真値 vs EnKF推定 vs LETKF推定の時系列（x_1〜x_10）
2. RMSE時系列比較（EnKF / ETKF / LETKF）
3. 局所化半径を変えたときのRMSE感度
```

---

### H2-5：EnKF vs ETKF vs LETKF の比較まとめ

| 項目 | EnKF（Stochastic） | ETKF | LETKF |
|---|---|---|---|
| 観測摂動 | あり | なし | なし |
| 局所化 | なし | なし | あり |
| 行列サイズ | $m \times m$ | $N \times N$ | $N \times N$（ローカル） |
| サンプリング誤差 | 大きい | 小さい | 小さい |
| 大規模対応 | ✗ | △ | ○ |
| 並列化 | △ | △ | ○（グリッド点独立） |

---

### FAQ セクション

```
Q1: ETKFとEnKFはどちらを使うべきですか？
A: 精度・安定性の観点ではETKF/LETKFが優れています。
   EnKF（Stochastic）はアルゴリズムが直感的でデバッグが容易なため、
   勉強・試作段階に向いています。
   実用的な大規模システムではLETKFが標準的に使われています。

Q2: 局所化半径はどうやって決めるのですか？
A: 観測密度・相関長さ・アンサンブルサイズに依存します。
   交差検証やRMSE最小化でチューニングします。
   気象モデルでは数百〜数千 kmが典型的です。

Q3: IEnKF（反復EnKF）とETKFの違いは？
A: ETKFは1回の線形更新で解析アンサンブルを求めます。
   IEnKF（Iterative EnKF）は観測演算子Hが強く非線形な場合に、
   更新を反復して非線形効果を考慮します。
   Hが弱非線形ならETKFで十分です。

Q4: LETKFはOpenFOAMに使えますか？
A: 可能です。CFDメッシュの各セルをグリッド点として扱い、
   センサーとの距離に基づいて局所化できます。
   並列化とも相性がよく、大規模CFD解析への適用が期待されています。
```

---

## Pythonファイル

- `python/letkf_lorenz96.py`：LETKF実装 + Lorenz 96デモ
- `python/enkf_vs_letkf.py`：EnKF / ETKF / LETKFのRMSE比較
- `python/localization_sensitivity.py`：局所化半径の感度分析
