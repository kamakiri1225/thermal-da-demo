# 記事02 計画：カルマンフィルタ

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | カルマンフィルタとは？予測と更新の仕組みを数式・Python・図で完全解説【データ同化入門】 |
| **メタディスクリプション** | カルマンフィルタは予測ステップと更新ステップを繰り返す最適な状態推定アルゴリズムです。数式の意味をゼロから解説し、PythonでトラッキングとLorenz系への応用を実装します。 |
| **メインキーワード** | カルマンフィルタ とは |
| **関連キーワード** | カルマンフィルタ 仕組み, カルマンフィルタ Python, 状態空間モデル, 予測ステップ 更新ステップ, カルマンフィルタ データ同化 |
| **想定文字数** | 5,000〜6,000字 |
| **Python可視化** | あり（1次元追跡・誤差共分散の時間発展） |

---

## 記事構成

### リード文
```
カルマンフィルタは1960年にRudolf E. Kálmánが提案した、
動的システムの状態を最適に推定するアルゴリズムです。
GPS・航空機の誘導制御・気象予報まで、
あらゆる「動くものの状態推定」に使われてきました。
前回解説したOI法を時間方向に拡張したものとして理解できます。
```

### この記事でわかること
- 予測ステップ・更新ステップの数式と、各変数が持つ意味がわかる
- OI法との違い（誤差共分散が時間とともに最適化される点）が理解できる
- PythonでカルマンフィルタをスクラッチImplementしてトラッキングを体験できる

---

### H2-1：状態空間モデルとは

**キーワード**：状態空間モデル, システムモデル 観測モデル

#### H3：2つの方程式

**システムモデル（時間発展）**：
$$\mathbf{x}_k = M \mathbf{x}_{k-1} + \boldsymbol{\eta}_k, \quad \boldsymbol{\eta}_k \sim \mathcal{N}(0, \mathbf{Q})$$

**観測モデル**：
$$\mathbf{y}_k = H \mathbf{x}_k + \boldsymbol{\epsilon}_k, \quad \boldsymbol{\epsilon}_k \sim \mathcal{N}(0, \mathbf{R})$$

| 記号 | 名称 | 意味 |
|---|---|---|
| $M$ | 状態遷移行列 | 時刻 $k-1$ から $k$ へのシステムの発展 |
| $\mathbf{Q}$ | システムノイズ共分散 | モデル化されない擾乱の大きさ |
| $H$ | 観測演算子 | 状態空間→観測空間への写像 |
| $\mathbf{R}$ | 観測ノイズ共分散 | センサー誤差の大きさ |

#### H3：OI法との違い
OI法は時刻を固定した「静的」な問題でした。
カルマンフィルタは $k = 1, 2, 3, \ldots$ と時刻を進め、
**毎ステップ最良推定を更新し続けます**。

---

### H2-2：カルマンフィルタの2ステップ

**キーワード**：カルマンフィルタ 予測 更新, カルマンフィルタ アルゴリズム

#### H3：ステップ1 — 予測（Forecast）

状態の予測：
$$\mathbf{x}^f_k = M \mathbf{x}^a_{k-1}$$

誤差共分散の予測：
$$\mathbf{P}^f_k = M \mathbf{P}^a_{k-1} M^T + \mathbf{Q}$$

- $\mathbf{P}^f_k$：時間発展によって誤差が「膨らむ」
- $\mathbf{Q}$ を加えることでシステムノイズ分の不確かさが増加

#### H3：ステップ2 — 更新（Analysis）

カルマンゲイン：
$$\mathbf{K}_k = \mathbf{P}^f_k H^T (H \mathbf{P}^f_k H^T + \mathbf{R})^{-1}$$

状態の更新（OI法と全く同じ形！）：
$$\mathbf{x}^a_k = \mathbf{x}^f_k + \mathbf{K}_k(\mathbf{y}_k - H\mathbf{x}^f_k)$$

誤差共分散の更新：
$$\mathbf{P}^a_k = (I - \mathbf{K}_k H)\mathbf{P}^f_k$$

#### H3：OI法との統一的理解

```
OI法  ：B は固定（手動チューニング）
KF    ：P^f は予測式で自動更新 → 時間とともに最適な B が得られる
```

---

### H2-3：直感的な理解 — 誤差楕円で考える

**キーワード**：カルマンフィルタ 直感, 誤差共分散 意味

**図の説明**（本文中に模式図を挿入）：
1. 初期状態：大きな誤差楕円（不確かさが大きい）
2. 予測後：楕円が少し膨らむ（システムノイズが加わる）
3. 観測後：楕円が潰れる（情報が増え不確かさが減る）
4. 次のステップへ：また少し膨らむ → 繰り返し

---

### H2-4：数値例 — 等速運動する物体の追跡

**キーワード**：カルマンフィルタ 例題, カルマンフィルタ トラッキング

#### H3：問題設定

- 1次元空間で等速運動する物体（位置 $x$、速度 $v$）
- 状態ベクトル：$\mathbf{x} = [x, v]^T$
- 状態遷移行列（$\Delta t = 1$ s）：

$$M = \begin{bmatrix} 1 & \Delta t \\ 0 & 1 \end{bmatrix}$$

- 位置のみ観測：$H = [1, 0]$
- 真値：$x_0 = 0$, $v = 2$ m/s、観測ノイズ $\sigma_r = 3$ m

#### H3：期待される結果
- 初期は観測に大きく引きずられる（不確かさが大きいため）
- 時間が経つとKFの推定値が滑らかになる（誤差共分散が収束）
- カルマンゲイン $K$ の時間発展グラフ：定常値に収束することを確認

---

### H2-5：PythonによるカルマンフィルタのFull実装

**キーワード**：カルマンフィルタ Python 実装, カルマンフィルタ コード

```python
# kalman_filter.py
import numpy as np
import matplotlib.pyplot as plt

class KalmanFilter:
    def __init__(self, M, H, Q, R, x0, P0):
        self.M = M
        self.H = H
        self.Q = Q
        self.R = R
        self.x = x0.copy()
        self.P = P0.copy()

    def predict(self):
        self.x = self.M @ self.x
        self.P = self.M @ self.P @ self.M.T + self.Q

    def update(self, y):
        S = self.H @ self.P @ self.H.T + self.R
        K = self.P @ self.H.T @ np.linalg.inv(S)
        innovation = y - self.H @ self.x
        self.x = self.x + K @ innovation
        self.P = (np.eye(len(self.x)) - K @ self.H) @ self.P
        return K

def simulate_tracking():
    dt = 1.0
    n_steps = 50
    np.random.seed(0)

    # システム設定
    M = np.array([[1, dt], [0, 1]])
    H = np.array([[1, 0]])
    Q = np.diag([0.1, 0.01])   # システムノイズ
    R = np.array([[9.0]])       # 観測ノイズ分散 (3m)^2

    # 真値生成
    x_true = np.zeros((n_steps, 2))
    x_true[0] = [0, 2]
    for k in range(1, n_steps):
        x_true[k] = M @ x_true[k-1] + np.random.multivariate_normal([0,0], Q)

    # 観測値生成（位置のみ）
    obs = x_true[:, 0] + np.random.normal(0, 3, n_steps)

    # カルマンフィルタ適用
    kf = KalmanFilter(M, H, Q, R,
                      x0=np.array([0.0, 0.0]),
                      P0=np.diag([100.0, 10.0]))

    estimates = np.zeros((n_steps, 2))
    gains = []
    stds = np.zeros(n_steps)

    for k in range(n_steps):
        kf.predict()
        K = kf.update(obs[k:k+1])
        estimates[k] = kf.x
        gains.append(K[0, 0])
        stds[k] = np.sqrt(kf.P[0, 0])

    return x_true, obs, estimates, gains, stds

# プロット
x_true, obs, est, gains, stds = simulate_tracking()
t = np.arange(len(x_true))

fig, axes = plt.subplots(2, 1, figsize=(10, 8))

# 上段：位置の推定
ax = axes[0]
ax.plot(t, x_true[:, 0], 'k--', label='真値')
ax.scatter(t, obs, s=20, alpha=0.5, color='gray', label='観測値')
ax.plot(t, est[:, 0], 'r-', linewidth=2, label='KF推定値')
ax.fill_between(t, est[:, 0] - 2*stds, est[:, 0] + 2*stds,
                alpha=0.2, color='red', label='±2σ')
ax.set_ylabel('位置 [m]')
ax.legend()
ax.grid(True, alpha=0.3)

# 下段：カルマンゲインの収束
ax = axes[1]
ax.plot(t, gains, 'b-o', markersize=3)
ax.set_xlabel('時刻ステップ')
ax.set_ylabel('カルマンゲイン K')
ax.set_title('カルマンゲインの時間発展（定常値に収束）')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('kalman_filter_tracking.png', dpi=150)
plt.show()
```

---

### H2-6：カルマンフィルタの限界と次のステップ

**キーワード**：カルマンフィルタ 限界, 非線形 カルマンフィルタ, EnKF

| 問題 | 内容 |
|---|---|
| 線形性の仮定 | $M$, $H$ が非線形の場合は厳密解ではない |
| 計算コスト | $\mathbf{P} \in \mathbb{R}^{n \times n}$：$n=10^6$ では $10^{12}$ 要素 |
| 正規分布の仮定 | 非ガウス誤差には対応不可 |

**次の記事へ**：高次元・非線形問題に対応するため、
$\mathbf{P}$ をアンサンブルで近似するEnKFを解説します。

**内部リンク**：← [最適補間法（OI）](../01_optimal_interpolation/plan.md) | → [EnKF](../03_enkf/plan.md)

---

### FAQ セクション

```
Q1: カルマンフィルタと移動平均の違いは何ですか？
A: 移動平均は過去の観測に等しい重みをかけるシンプルな平滑化です。
   カルマンフィルタは物理モデル（状態遷移行列M）と誤差統計（P, R）を
   使って時変の最適ゲインを計算します。
   システムの動特性を知っている場合はKFが大幅に優れています。

Q2: 定常カルマンゲインとは何ですか？
A: 定常系（M, H, Q, R が時不変）ではKが一定値に収束します。
   これを定常カルマンゲインといい、事前に計算しておけます。
   代数リカッチ方程式（DARE）で求められます。

Q3: 拡張カルマンフィルタ（EKF）とは？
A: MやHが非線形の場合、ヤコビアンで線形化してKFを適用する手法です。
   線形化誤差があるため、大きな非線形性では精度が低下します。
   EnKFはこの問題をサンプリングで回避します。

Q4: カルマンフィルタはリアルタイム処理に使えますか？
A: はい。カルマンフィルタは逐次更新型なので、
   観測が届くたびに1ステップだけ計算すれば済み、
   リアルタイム処理に適しています。
```

---

## Pythonファイル

- `python/kalman_filter.py`：KalmanFilterクラスの実装と追跡デモ
- `python/kalman_convergence.py`：カルマンゲイン・誤差共分散の収束分析
