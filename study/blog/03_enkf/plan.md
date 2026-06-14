# 記事03 計画：アンサンブルカルマンフィルタ（EnKF）

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | アンサンブルカルマンフィルタ（EnKF）とは？非線形・大規模問題への拡張をPythonで実装 |
| **メタディスクリプション** | アンサンブルカルマンフィルタ（EnKF）は、アンサンブルで誤差共分散を近似し高次元・非線形問題に対応するデータ同化手法です。Lorenz 63モデルへの適用例をPythonでゼロから実装します。 |
| **メインキーワード** | アンサンブルカルマンフィルタ EnKF |
| **関連キーワード** | EnKF Python 実装, アンサンブルカルマンフィルタ 仕組み, Lorenz 63 データ同化, EnKF 非線形, 誤差共分散 アンサンブル |
| **想定文字数** | 5,500〜7,000字 |
| **Python可視化** | あり（Lorenz 63 アンサンブル軌跡・RMSE vs スプレッド） |

---

## 記事構成

### リード文
```
カルマンフィルタは最適な状態推定アルゴリズムですが、
気象予報のような超高次元問題（状態変数 $n \sim 10^8$）では
$n \times n$ の共分散行列を扱うことは不可能です。
アンサンブルカルマンフィルタ（EnKF）は、
この問題を「複数のモデルランを同時実行する」
モンテカルロ的アプローチで解決します。
```

### この記事でわかること
- なぜKFが大規模問題に使えないかと、EnKFがどう解決するかがわかる
- アンサンブルによる誤差共分散の近似の数式と意味が理解できる
- PythonでLorenz 63モデルにEnKFを適用し、カオスのトラッキングを体験できる

---

### H2-1：カルマンフィルタの計算上の限界

**キーワード**：カルマンフィルタ 高次元 問題, EnKF 動機

#### H3：次元の呪い

| 問題規模 | 状態次元 $n$ | $\mathbf{P}$ のサイズ | メモリ（64bit float） |
|---|---|---|---|
| 教科書例題 | 10 | 100要素 | 0.8 KB |
| 地域気象モデル | $10^5$ | $10^{10}$ 要素 | 80 GB |
| 全球気象モデル | $10^8$ | $10^{16}$ 要素 | **80,000 TB** |

**結論**：$n > 10^3$ 程度でKFは実用不可能。

#### H3：EnKFのアイデア

$N$ 個のアンサンブルメンバー $\{\mathbf{x}^{(i)}\}_{i=1}^{N}$ で統計を近似：

$$\bar{\mathbf{x}} = \frac{1}{N}\sum_{i=1}^{N}\mathbf{x}^{(i)}$$

$$\mathbf{P} \approx \mathbf{P}^{ens} = \frac{1}{N-1}\sum_{i=1}^{N}(\mathbf{x}^{(i)} - \bar{\mathbf{x}})(\mathbf{x}^{(i)} - \bar{\mathbf{x}})^T$$

- $N = 20 \sim 100$ 程度で実用的な精度が得られる
- $n \times N$ のアンサンブル行列だけ保持すればよい（$n \times n$ 行列は不要）

---

### H2-2：EnKFのアルゴリズム（Stochastic EnKF）

**キーワード**：EnKF アルゴリズム, Stochastic EnKF

#### H3：ステップ1 — アンサンブル予測

各メンバーを個別にシステムモデルで発展させる：
$$\mathbf{x}^{f,(i)}_k = M(\mathbf{x}^{a,(i)}_{k-1}) + \boldsymbol{\eta}^{(i)}_k, \quad \boldsymbol{\eta}^{(i)}_k \sim \mathcal{N}(0, \mathbf{Q})$$

- $M$ が非線形でも適用可能（差分化・ヤコビアン不要）

#### H3：ステップ2 — アンサンブル更新

観測摂動（perturbed observations）：
$$\mathbf{y}^{(i)}_k = \mathbf{y}_k + \boldsymbol{\epsilon}^{(i)}_k, \quad \boldsymbol{\epsilon}^{(i)}_k \sim \mathcal{N}(0, \mathbf{R})$$

アンサンブルカルマンゲイン：
$$\mathbf{K}^{ens}_k = \mathbf{P}^{f,ens}_k H^T(H\mathbf{P}^{f,ens}_k H^T + \mathbf{R})^{-1}$$

各メンバーの更新：
$$\mathbf{x}^{a,(i)}_k = \mathbf{x}^{f,(i)}_k + \mathbf{K}^{ens}_k(\mathbf{y}^{(i)}_k - H\mathbf{x}^{f,(i)}_k)$$

#### H3：KFとEnKFの対応表

| 量 | カルマンフィルタ | EnKF |
|---|---|---|
| 誤差共分散 | $\mathbf{P}$（行列演算で厳密計算） | $\mathbf{P}^{ens}$（アンサンブルで近似） |
| 非線形対応 | 不可（線形化が必要） | 可能（各メンバーを非線形モデルで発展） |
| 計算コスト | $O(n^3)$ | $O(N \cdot n)$ |
| 必要な条件 | ガウス・線形 | 近似的にガウス（緩い条件） |

---

### H2-3：Lorenz 63モデル — カオスへの挑戦

**キーワード**：Lorenz 63 データ同化, カオス 状態推定, EnKF Lorenz

#### H3：Lorenz 63モデルとは

カオスを示す3変数の常微分方程式：
$$\frac{dx}{dt} = \sigma(y - x)$$
$$\frac{dy}{dt} = x(\rho - z) - y$$
$$\frac{dz}{dt} = xy - \beta z$$

標準パラメータ：$\sigma = 10$, $\rho = 28$, $\beta = 8/3$

**EnKFのテスト問題として最適な理由**：
- 3次元で可視化が容易
- カオス的な感度（わずかな誤差が急速に拡大）
- 非線形の時間発展

#### H3：観測設定

- $x$ 成分のみを $\Delta t_{obs} = 0.08$ ごとに観測
- 観測ノイズ：$\sigma_r = 2.0$
- アンサンブルサイズ：$N = 20$

---

### H2-4：PythonによるEnKF実装（Lorenz 63）

**キーワード**：EnKF Python コード, アンサンブルカルマンフィルタ 実装

```python
# enkf_lorenz63.py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# ===== Lorenz 63 モデル =====
def lorenz63(state, sigma=10, rho=28, beta=8/3):
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return np.array([dx, dy, dz])

def rk4_step(state, dt, model):
    k1 = model(state)
    k2 = model(state + 0.5 * dt * k1)
    k3 = model(state + 0.5 * dt * k2)
    k4 = model(state + dt * k3)
    return state + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

# ===== アンサンブルカルマンフィルタ =====
class EnKF:
    def __init__(self, N, n, H, R, Q_std):
        self.N = N          # アンサンブルサイズ
        self.n = n          # 状態次元
        self.H = H          # 観測演算子
        self.R = R          # 観測誤差共分散
        self.Q_std = Q_std  # システムノイズの標準偏差

    def predict(self, ensemble, dt):
        """各メンバーを非線形モデルで発展"""
        new_ens = np.zeros_like(ensemble)
        for i in range(self.N):
            new_ens[:, i] = rk4_step(ensemble[:, i], dt, lorenz63)
            new_ens[:, i] += np.random.normal(0, self.Q_std, self.n)
        return new_ens

    def update(self, ensemble, y_obs):
        """観測でアンサンブルを更新"""
        m = len(y_obs)

        # アンサンブル平均と偏差
        x_mean = ensemble.mean(axis=1, keepdims=True)
        A = ensemble - x_mean  # アンサンブル偏差行列 (n×N)

        # アンサンブルによる共分散近似
        Pf = A @ A.T / (self.N - 1)

        # カルマンゲイン
        S = self.H @ Pf @ self.H.T + self.R
        K = Pf @ self.H.T @ np.linalg.inv(S)

        # 観測摂動と更新
        new_ens = np.zeros_like(ensemble)
        for i in range(self.N):
            eps = np.random.multivariate_normal(np.zeros(m), self.R)
            y_perturbed = y_obs + eps
            innovation = y_perturbed - self.H @ ensemble[:, i]
            new_ens[:, i] = ensemble[:, i] + K @ innovation

        return new_ens

# ===== シミュレーション =====
def run_enkf_lorenz63():
    np.random.seed(42)
    n, N = 3, 20
    dt_model = 0.01
    dt_obs = 0.08
    steps_per_obs = int(dt_obs / dt_model)
    n_obs = 200

    H = np.array([[1, 0, 0]])       # x成分のみ観測
    R = np.array([[4.0]])           # 観測ノイズ分散
    Q_std = 0.1

    enkf = EnKF(N, n, H, R, Q_std * np.eye(n))

    # 真値の生成（スピンアップ付き）
    x_true = np.array([1.0, 1.0, 1.0])
    for _ in range(2000):
        x_true = rk4_step(x_true, dt_model, lorenz63)

    # アンサンブル初期化（真値付近にランダム散布）
    ensemble = x_true[:, None] + np.random.normal(0, 2, (n, N))

    # 時間発展
    true_states, obs_list, mean_estimates = [], [], []

    for k in range(n_obs):
        # モデル積分
        for _ in range(steps_per_obs):
            x_true = rk4_step(x_true, dt_model, lorenz63)
        ensemble = enkf.predict(ensemble, dt_model * steps_per_obs)

        # 観測
        y = H @ x_true + np.random.normal(0, 2, 1)

        # 更新
        ensemble = enkf.update(ensemble, y)

        true_states.append(x_true.copy())
        obs_list.append(y[0])
        mean_estimates.append(ensemble.mean(axis=1))

    return (np.array(true_states), np.array(obs_list),
            np.array(mean_estimates), ensemble)

# ===== 可視化 =====
true_states, obs_list, estimates, final_ens = run_enkf_lorenz63()
t = np.arange(len(true_states)) * 0.08
rmse = np.sqrt(((true_states - estimates)**2).mean(axis=1))

fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
labels = ['x', 'y', 'z']
for i, ax in enumerate(axes):
    ax.plot(t, true_states[:, i], 'k-', lw=1.5, label='真値')
    ax.plot(t, estimates[:, i], 'r--', lw=1.5, label='EnKF推定値')
    if i == 0:
        ax.scatter(t, obs_list, s=5, color='gray', alpha=0.5, label='観測値')
    ax.set_ylabel(labels[i])
    ax.legend(loc='upper right', fontsize=8)
    ax.grid(True, alpha=0.3)
axes[-1].set_xlabel('時刻 t')
plt.suptitle('EnKF による Lorenz 63 モデルの状態推定', y=1.01)
plt.tight_layout()
plt.savefig('enkf_lorenz63_timeseries.png', dpi=150)
plt.show()
```

**可視化パターン2：3Dアトラクター**
```python
# 真値・EnKF推定値をLorentzアトラクター上で3Dプロット
# → カオス軌道を追跡できていることを視覚化
```

**可視化パターン3：RMSE vs アンサンブルスプレッド**
```python
# スプレッド = sqrt(trace(P^ens) / n)
# RMSEとスプレッドが一致 → アンサンブルが正しく誤差を表現している指標
```

---

### H2-5：EnKFの実装上の工夫

**キーワード**：EnKF 局所化, EnKF インフレーション, アンサンブルカルマンフィルタ 実装

#### H3：アンサンブルの劣化問題

- **フィルター発散**：アンサンブルスプレッドが小さくなりすぎると観測を無視する
- **対策1 — インフレーション**：$\mathbf{P}^{ens} \to \alpha \mathbf{P}^{ens}$（$\alpha > 1$）
  
#### H3：局所化（Localization）

- 遠く離れた観測点の影響をゼロに近づける
- $\mathbf{P}^{ens}$ にSchur積で距離減衰関数をかける
- 大規模問題で必須のテクニック

---

### H2-6：シリーズまとめ — OI・KF・EnKFの統一的理解

**キーワード**：OI カルマンフィルタ EnKF 比較, データ同化 まとめ

| 手法 | 時間発展 | 非線形 | $n$ の規模 | 共分散 |
|---|---|---|---|---|
| OI | なし | 不可 | 小〜中 | 固定 $\mathbf{B}$ |
| KF | あり | 不可 | 小 | 厳密計算 |
| EnKF | あり | 可能 | 大規模 | アンサンブル近似 |

**今後の発展**：
- 変分法（3D-Var / 4D-Var）との比較
- OpenFOAMへの組み込み（[thermal-da-demo](https://github.com/kamakiri1225/thermal-da-demo)リポジトリで公開予定）

**内部リンク**：← [カルマンフィルタ](../02_kalman_filter/plan.md)

---

### FAQ セクション

```
Q1: アンサンブルサイズNはいくつにすればよいですか？
A: 問題の次元や非線形性によりますが、気象予報では20〜100程度が一般的です。
   Nが小さすぎるとサンプリング誤差が大きくなります。
   局所化と組み合わせることで少ないNでも実用的な精度が得られます。

Q2: EnKFとパーティクルフィルタの違いは何ですか？
A: パーティクルフィルタはガウス分布を仮定せず、
   任意の確率分布を扱えます。ただし高次元では重み縮退が起きやすく、
   EnKFより多くのメンバーが必要です。
   EnKFは近似的にガウス性を仮定しますが、高次元問題に実用的です。

Q3: EnKFはOpenFOAMに適用できますか？
A: 可能です。OpenFOAMの各タイムステップを「予測ステップ」として扱い、
   センサーデータをもとに各アンサンブルメンバーの場を更新します。
   PDAfやOpenDAなどのフレームワークを使う方法もあります。

Q4: Deterministic EnKF（ETKF）との違いは？
A: Stochastic EnKFは観測に確率的な摂動を加えます。
   ETKF（Ensemble Transform KF）は摂動なしで
   アンサンブル変換行列を使って更新し、サンプリング誤差を低減します。
   実用的にはETKFやLESTKFが多く使われます。
```

---

## Pythonファイル

- `python/enkf_lorenz63.py`：EnKFクラス + Lorenz 63デモ（時系列・3Dアトラクター）
- `python/enkf_diagnostics.py`：RMSE vs スプレッド・インフレーション効果の分析
