# 記事05 計画：粒子フィルタ（Particle Filter）

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | 粒子フィルタとは？非ガウス・非線形の状態推定をPythonでゼロから実装【データ同化入門】 |
| **メタディスクリプション** | 粒子フィルタ（パーティクルフィルタ）はガウス分布を仮定しない完全非線形な状態推定手法です。重要度サンプリングとリサンプリングの仕組みを丁寧に解説し、PythonでLorenz 63への適用を実装します。 |
| **メインキーワード** | 粒子フィルタ とは |
| **関連キーワード** | パーティクルフィルタ Python, 粒子フィルタ 実装, 非ガウス 状態推定, 重要度サンプリング, リサンプリング |
| **想定文字数** | 5,000〜6,000字 |
| **Python可視化** | あり（重み分布・リサンプリング前後・EnKFとの比較） |

---

## 記事構成

### リード文
```
カルマンフィルタ系の手法は誤差がガウス分布に従うことを前提としています。
しかし現実の非線形システムでは、事後分布が多峰性になったり、
非対称な形をとることがあります。
粒子フィルタはそうした「ガウス分布の仮定が使えない問題」への
究極の答えです。
```

### この記事でわかること
- 重要度サンプリングとリサンプリングの仕組みが理解できる
- EnKF系との違い（ガウス性の仮定なし）とトレードオフがわかる
- PythonでSIR粒子フィルタをスクラッチ実装できる

---

### H2-1：粒子フィルタの基本アイデア

**キーワード**：粒子フィルタ 仕組み, 重要度サンプリング

#### H3：確率分布を「粒子（サンプル）」で表す

EnKF系の限界：
$$p(\mathbf{x}|\mathbf{y}) \approx \mathcal{N}(\bar{\mathbf{x}}^a, \mathbf{P}^a) \quad \text{（ガウス近似）}$$

粒子フィルタのアプローチ：
$$p(\mathbf{x}|\mathbf{y}) \approx \sum_{i=1}^{N} w^{(i)} \delta(\mathbf{x} - \mathbf{x}^{(i)})$$

- $N$ 個の粒子 $\mathbf{x}^{(i)}$ と重み $w^{(i)}$ で任意の分布を表現
- $\sum_i w^{(i)} = 1$
- ガウス性・線形性の仮定が一切不要

---

### H2-2：SIR粒子フィルタのアルゴリズム

**キーワード**：SIR粒子フィルタ, Sequential Importance Resampling

#### H3：ステップ1 — 予測（状態の伝播）

各粒子をシステムモデルで発展：
$$\mathbf{x}^{(i)}_k = M(\mathbf{x}^{(i)}_{k-1}) + \boldsymbol{\eta}^{(i)}_k$$

#### H3：ステップ2 — 重み更新（尤度評価）

観測 $\mathbf{y}_k$ に対する各粒子の尤度で重みを更新：
$$\tilde{w}^{(i)}_k = w^{(i)}_{k-1} \cdot p(\mathbf{y}_k | \mathbf{x}^{(i)}_k)$$

ガウス観測モデルの場合：
$$p(\mathbf{y}_k | \mathbf{x}^{(i)}_k) = \mathcal{N}(\mathbf{y}_k; H\mathbf{x}^{(i)}_k, \mathbf{R})$$

正規化：
$$w^{(i)}_k = \frac{\tilde{w}^{(i)}_k}{\sum_j \tilde{w}^{(j)}_k}$$

#### H3：ステップ3 — リサンプリング（重み縮退の防止）

**重み縮退（weight degeneracy）問題**：
- 時間が経つと1つの粒子に重みが集中し、他は $w \approx 0$ になる
- 有効サンプルサイズ：$N_{eff} = \frac{1}{\sum_i (w^{(i)})^2}$
- $N_{eff} < N/2$ になったらリサンプリング実行

**多項リサンプリング**：
- 重み $\{w^{(i)}\}$ を確率とみなして $N$ 個を復元抽出
- リサンプリング後は全重みを $1/N$ にリセット

---

### H2-3：EnKFとの比較

**キーワード**：粒子フィルタ EnKF 違い, データ同化 手法比較

| 観点 | EnKF | 粒子フィルタ |
|---|---|---|
| 分布の表現 | ガウス（平均＋共分散） | 任意（粒子＋重み） |
| 非線形 | 近似的に対応 | 完全対応 |
| 非ガウス | 不可 | 可能 |
| 高次元 | 比較的得意（局所化で対応） | 苦手（次元の呪い） |
| 必要粒子数 | $N = 20 \sim 100$ | $N$ は次元に指数的に増加 |

**実用上の選択基準**：
- 低次元（$n < 10$）かつ非ガウス → 粒子フィルタ
- 高次元かつ近似的にガウス → LETKF
- 中間 → ハイブリッド手法（粒子EnKFなど）

---

### H2-4：「次元の呪い」を直感的に理解する

**キーワード**：次元の呪い 粒子フィルタ, 粒子フィルタ 高次元 問題

- $n$ 次元のガウス分布を $N$ 粒子でカバーするには？
- 1次元：$N = 100$ で十分
- 10次元：各次元10点 → $10^{10}$ 粒子が必要（理論的）
- 気象モデル（$n \sim 10^8$）には事実上不可能

**現実的な対策**：
- 局所化粒子フィルタ（Localized PF）
- 粒子フィルタ × EnKFのハイブリッド（Ensemble Particle Filter）

---

### H2-5：PythonによるSIR粒子フィルタ実装（Lorenz 63）

**キーワード**：粒子フィルタ Python コード, SIR フィルタ 実装

```python
# particle_filter_lorenz63.py
import numpy as np
import matplotlib.pyplot as plt

def lorenz63(state, sigma=10, rho=28, beta=8/3):
    x, y, z = state
    return np.array([
        sigma * (y - x),
        x * (rho - z) - y,
        x * y - beta * z
    ])

def rk4_step(state, dt, model):
    k1 = model(state)
    k2 = model(state + 0.5*dt*k1)
    k3 = model(state + 0.5*dt*k2)
    k4 = model(state + dt*k3)
    return state + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)

class ParticleFilter:
    def __init__(self, N, n, H, R, Q_std):
        self.N = N
        self.n = n
        self.H = H
        self.R = R
        self.Q_std = Q_std
        self.weights = np.ones(N) / N

    def predict(self, particles, dt):
        new_particles = np.zeros_like(particles)
        for i in range(self.N):
            new_particles[:, i] = rk4_step(particles[:, i], dt, lorenz63)
            new_particles[:, i] += np.random.normal(0, self.Q_std, self.n)
        return new_particles

    def update(self, particles, y_obs):
        """尤度で重みを更新"""
        log_weights = np.zeros(self.N)
        for i in range(self.N):
            residual = y_obs - self.H @ particles[:, i]
            # ガウス対数尤度
            log_weights[i] = -0.5 * residual @ np.linalg.solve(self.R, residual)

        # オーバーフロー対策（log-sum-exp）
        log_weights -= log_weights.max()
        self.weights = np.exp(log_weights)
        self.weights /= self.weights.sum()
        return self.weights

    def resample(self, particles):
        """多項リサンプリング"""
        N_eff = 1.0 / (self.weights ** 2).sum()
        if N_eff < self.N / 2:
            indices = np.random.choice(self.N, self.N, p=self.weights)
            particles = particles[:, indices]
            self.weights = np.ones(self.N) / self.N
        return particles

    def estimate(self, particles):
        """加重平均で状態推定"""
        return particles @ self.weights

# --- シミュレーション ---
np.random.seed(0)
n, N = 3, 500
dt = 0.01
dt_obs = 0.08
steps_per_obs = int(dt_obs / dt)
n_obs = 200

H = np.array([[1, 0, 0]])
R = np.array([[4.0]])
Q_std = 0.1

pf = ParticleFilter(N, n, H, R, Q_std)

# 真値スピンアップ
x_true = np.array([1.0, 1.0, 1.0])
for _ in range(2000):
    x_true = rk4_step(x_true, dt, lorenz63)

# 粒子初期化
particles = x_true[:, None] + np.random.normal(0, 3, (n, N))

true_states, obs_list, pf_estimates = [], [], []

for k in range(n_obs):
    for _ in range(steps_per_obs):
        x_true = rk4_step(x_true, dt, lorenz63)
    particles = pf.predict(particles, dt * steps_per_obs)

    y = H @ x_true + np.random.normal(0, 2, 1)
    pf.update(particles, y)
    particles = pf.resample(particles)

    true_states.append(x_true.copy())
    obs_list.append(y[0])
    pf_estimates.append(pf.estimate(particles))

# --- プロット ---
true_states = np.array(true_states)
pf_estimates = np.array(pf_estimates)
t = np.arange(n_obs) * dt_obs

fig, axes = plt.subplots(3, 1, figsize=(12, 8), sharex=True)
for i, (ax, label) in enumerate(zip(axes, ['x', 'y', 'z'])):
    ax.plot(t, true_states[:, i], 'k-', lw=1.5, label='真値')
    ax.plot(t, pf_estimates[:, i], 'b--', lw=1.5, label='粒子フィルタ推定')
    if i == 0:
        ax.scatter(t, obs_list, s=5, color='gray', alpha=0.5, label='観測値')
    ax.set_ylabel(label)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
axes[-1].set_xlabel('時刻 t')
plt.suptitle(f'粒子フィルタ（N={N}）による Lorenz 63 状態推定')
plt.tight_layout()
plt.savefig('particle_filter_lorenz63.png', dpi=150)
plt.show()
```

**可視化パターン2：粒子の重み分布の可視化**
```
→ 観測直後に重みが1点に集中する様子 → リサンプリング後に均等化
→ 「重み縮退」を視覚的に体験させる
```

**可視化パターン3：粒子数 N の影響**
```
N = 50, 100, 500, 1000 で RMSE を比較
→ N が多いほど精度が上がることを確認
→ EnKF（N=20）と粒子フィルタ（N=500）のRMSE比較
```

---

### FAQ セクション

```
Q1: 粒子フィルタは実用的に使われていますか？
A: 低次元問題（ロボット自己位置推定・金融など）では広く使われています。
   高次元の気象・海洋分野では次元の呪いのため直接適用は難しく、
   局所化粒子フィルタや粒子EnKFなどの改良版が研究されています。

Q2: リサンプリング方法にはどんな種類がありますか？
A: 多項リサンプリング（最もシンプル）、系統的リサンプリング
   （サンプリング分散が小さい）、残差リサンプリングなどがあります。
   実用的には系統的リサンプリングがよく使われます。

Q3: 粒子フィルタと変分型（4D-Var）の違いは？
A: 粒子フィルタは状態の確率分布全体を推定します。
   4D-Varはコスト関数を最小化して最良推定値を1点求めます。
   確率分布が必要かどうかで使い分けます。

Q4: N_eff の閾値 N/2 は絶対的なものですか？
A: いいえ。問題に応じてN/4〜N/2の範囲でチューニングします。
   リサンプリングが多すぎると粒子の多様性が失われ（サンプリング貧困化）、
   少なすぎると重み縮退が進みます。
```

---

## Pythonファイル

- `python/particle_filter.py`：ParticleFilterクラス + Lorenz 63デモ
- `python/pf_weight_visualization.py`：重み縮退とリサンプリングの可視化
- `python/pf_vs_enkf.py`：粒子フィルタ vs EnKFのRMSE比較
