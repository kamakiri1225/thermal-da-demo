# 記事01 計画：最適補間法（OI法）

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | 最適補間法（OI法）とは？データ同化の基礎を数式とPythonで徹底解説 |
| **メタディスクリプション** | 最適補間法（OI法）は、背景場と観測データを最小分散推定で統合するデータ同化の基礎手法です。カルマンゲインの直感的な意味から、Pythonによる1次元実装まで丁寧に解説します。 |
| **メインキーワード** | 最適補間法 OI法 |
| **関連キーワード** | 最適補間法 数式, OI法 データ同化, カルマンゲイン, 最小分散推定, 背景誤差共分散 |
| **想定文字数** | 4,000〜5,000字 |
| **Python可視化** | あり（1次元補間・B/Rパラメータ感度） |

---

## 記事構成

### リード文
```
数値シミュレーションの予測値（背景場）と、
センサーから得た観測値——どちらも誤差を含んでいます。
最適補間法（OI法）は、両者の誤差を考慮して
「最も確からしい状態」を導く、データ同化の出発点となる手法です。
```

### この記事でわかること
- 最小分散推定の考え方と、なぜOI法が「最適」なのかがわかる
- カルマンゲイン行列の意味と、B・Rが果たす役割が理解できる
- PythonでOI法を実装して、補間結果を可視化できる

---

### H2-1：最適補間法の基本アイデア

**キーワード**：最適補間法 とは, 最小分散推定

**内容の流れ**
1. 問題設定：背景場 $\mathbf{x}^b$ と観測 $\mathbf{y}$ をどう統合するか
2. 直感：「信頼度に応じた加重平均」として捉える
3. 最小分散推定：解析場の誤差分散を最小化する $\mathbf{K}$ を求める

**解析方程式（OI の核心）**：
$$\mathbf{x}^a = \mathbf{x}^b + \mathbf{K}(\mathbf{y} - H\mathbf{x}^b)$$

| 項 | 名称 | 直感的な意味 |
|---|---|---|
| $\mathbf{x}^a$ | 解析場（Analysis） | 融合後の最良推定値 |
| $\mathbf{x}^b$ | 背景場（Background） | モデルの事前予測 |
| $\mathbf{y} - H\mathbf{x}^b$ | イノベーション | 観測と予測のズレ（残差） |
| $\mathbf{K}$ | カルマンゲイン | イノベーションをどれだけ信頼するか |

---

### H2-2：カルマンゲインの導出

**キーワード**：カルマンゲイン 導出, 最小分散推定 データ同化

#### H3：コスト関数による導出（変分的視点）

コスト関数（背景場と観測からの二乗誤差の加重和）：
$$J(\mathbf{x}) = \frac{1}{2}(\mathbf{x}-\mathbf{x}^b)^T \mathbf{B}^{-1}(\mathbf{x}-\mathbf{x}^b) + \frac{1}{2}(\mathbf{y}-H\mathbf{x})^T \mathbf{R}^{-1}(\mathbf{y}-H\mathbf{x})$$

$\nabla J = 0$ を解くと：
$$\mathbf{x}^a = \mathbf{x}^b + \mathbf{B}H^T(H\mathbf{B}H^T + \mathbf{R})^{-1}(\mathbf{y} - H\mathbf{x}^b)$$

よってカルマンゲイン：
$$\mathbf{K} = \mathbf{B}H^T(H\mathbf{B}H^T + \mathbf{R})^{-1}$$

#### H3：カルマンゲインの直感的意味

| 条件 | $\mathbf{K}$ の振る舞い | 解釈 |
|---|---|---|
| $\mathbf{B} \gg \mathbf{R}$（モデル誤差大） | $\mathbf{K} \to H^{-1}$（大きい） | 観測を強く信頼する |
| $\mathbf{B} \ll \mathbf{R}$（観測誤差大） | $\mathbf{K} \to 0$（小さい） | モデル予測を維持する |
| $\mathbf{B} = \mathbf{R}$ | 中間の値 | 半々に混ぜる |

#### H3：解析場の誤差共分散

$$\mathbf{P}^a = (I - \mathbf{K}H)\mathbf{B}$$

- データ同化後は必ず $\mathbf{P}^a \leq \mathbf{B}$（誤差は減少する）

---

### H2-3：1次元の具体例で理解する

**キーワード**：最適補間法 例題, OI法 1次元

#### H3：問題設定

- 1次元空間（x = 0〜100 km）での温度場を推定
- 背景場：$x^b = 20 + 0.05x$（傾きのある線形分布）
- 観測点：$x = 20, 50, 80$ km で温度を観測（ノイズあり）
- 背景誤差分散：$\sigma_b^2 = 4.0$、観測誤差分散：$\sigma_r^2 = 1.0$

#### H3：結果の解釈
- 観測点付近では解析場が観測値に引き寄せられる
- 観測点から離れると背景場に戻っていく
- 引き寄せの強さは $\sigma_b/\sigma_r$ の比で決まる

---

### H2-4：PythonによるOI法の実装

**キーワード**：最適補間法 Python, OI法 実装

```python
# oi_1d.py
import numpy as np
import matplotlib.pyplot as plt

def gaussian_covariance(x1, x2, sigma2, L):
    """ガウス型背景誤差共分散"""
    return sigma2 * np.exp(-0.5 * ((x1 - x2) / L) ** 2)

def optimal_interpolation(x_grid, xb, obs_locs, obs_vals, sigma_b2, L, sigma_r2):
    """
    Parameters
    ----------
    x_grid   : グリッド点座標 (n,)
    xb       : 背景場 (n,)
    obs_locs : 観測点座標 (m,)
    obs_vals : 観測値 (m,)
    sigma_b2 : 背景誤差分散
    L        : 背景誤差の相関長さ
    sigma_r2 : 観測誤差分散
    """
    n = len(x_grid)
    m = len(obs_locs)

    # 背景誤差共分散行列 B (n×n) と BH^T (n×m)
    B = np.array([[gaussian_covariance(x_grid[i], x_grid[j], sigma_b2, L)
                   for j in range(n)] for i in range(n)])

    BHT = np.array([[gaussian_covariance(x_grid[i], obs_locs[j], sigma_b2, L)
                     for j in range(m)] for i in range(n)])

    # HBH^T (m×m)
    HBHT = np.array([[gaussian_covariance(obs_locs[i], obs_locs[j], sigma_b2, L)
                      for j in range(m)] for i in range(m)])

    # 観測誤差共分散 R
    R = sigma_r2 * np.eye(m)

    # カルマンゲイン K = BH^T (HBH^T + R)^{-1}
    K = BHT @ np.linalg.inv(HBHT + R)

    # 背景場を観測点で評価（H xb）
    xb_at_obs = np.interp(obs_locs, x_grid, xb)

    # 解析場
    innovation = obs_vals - xb_at_obs
    xa = xb + K @ innovation

    return xa, K

# メイン
np.random.seed(42)
x_grid = np.linspace(0, 100, 200)
xb = 20 + 0.05 * x_grid                        # 背景場（真値に近い）
x_true = 20 + 0.05 * x_grid + 2 * np.sin(x_grid / 15)  # 真の状態

obs_locs = np.array([20.0, 50.0, 80.0])
obs_vals = np.interp(obs_locs, x_grid, x_true) + np.random.normal(0, 1, 3)

xa, K = optimal_interpolation(
    x_grid, xb, obs_locs, obs_vals,
    sigma_b2=4.0, L=15.0, sigma_r2=1.0
)

# プロット
fig, ax = plt.subplots(figsize=(10, 5))
ax.plot(x_grid, x_true, 'k--', label='真の状態')
ax.plot(x_grid, xb, 'b-', alpha=0.6, label='背景場 $x^b$')
ax.plot(x_grid, xa, 'r-', linewidth=2, label='解析場 $x^a$（OI）')
ax.scatter(obs_locs, obs_vals, s=80, color='green', zorder=5,
           label='観測値', marker='o')
ax.set_xlabel('x [km]')
ax.set_ylabel('温度 [°C]')
ax.set_title('最適補間法（OI）による1次元データ同化')
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('oi_1d_result.png', dpi=150)
plt.show()
```

**可視化パターン2：B/Rパラメータ感度**
```python
# sigma_b2 を変えたとき解析場がどう変わるか比較プロット
# → 「Bが大きいほど観測に引き寄せられる」を視覚的に確認
```

---

### H2-5：OI法の限界と次のステップ

**キーワード**：OI法 限界, データ同化 カルマンフィルタ

| 特性 | OI法 | カルマンフィルタ（次回） |
|---|---|---|
| 時間発展 | なし（静的） | あり（動的） |
| $\mathbf{B}$ の更新 | 固定（チューニング必要） | 時間とともに最適化 |
| 計算コスト | 低い | 中程度 |
| 適用問題 | 空間補間・再解析 | 時系列状態推定 |

**次の記事へ**：OI法は時間が止まった状態での最良推定でした。
次回は時間発展するシステムで繰り返しデータ同化を行う「カルマンフィルタ」を解説します。

**内部リンク**：← [データ同化とは？](../00_introduction/plan.md) | → [カルマンフィルタ](../02_kalman_filter/plan.md)

---

### FAQ セクション

```
Q1: OI法と補間（スプライン補間など）の違いは何ですか？
A: スプライン補間は観測点を通る曲線を求めますが、
   OI法は背景場と観測誤差の統計的な重みを考慮します。
   観測値に誤差がある場合、OI法は観測点を「通り過ぎず」
   適切に平滑化します。

Q2: 観測演算子 H が非線形の場合はどうなりますか？
A: 非線形の場合は H をヤコビアン行列で線形近似します（拡張KFなど）。
   大きく非線形な場合はEnKFが有効です。

Q3: 背景誤差共分散 B はどうやって決めるのですか？
A: 実務では NMC法（モデル予測差の統計）や
   アンサンブル法で推定します。
   チューニングが必要な重要パラメータです。

Q4: OI法はOpenFOAMに適用できますか？
A: 可能です。OpenFOAMのフィールドデータを状態ベクトルとして扱い、
   センサーデータをOI法で同化する研究事例があります。
```

---

## Pythonファイル

- `python/oi_1d.py`：1次元OI法の実装と可視化
- `python/oi_sensitivity.py`：B/Rパラメータ感度分析
