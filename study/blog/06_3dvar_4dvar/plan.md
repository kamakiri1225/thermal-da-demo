# 記事06 計画：3D-Var・4D-Var（変分型データ同化）

## メタ情報（SEO設定）

| 項目 | 内容 |
|---|---|
| **タイトル** | 3D-Var・4D-Varとは？変分型データ同化の理論と気象予報への応用をPythonで解説 |
| **メタディスクリプション** | 3D-Var・4D-Varは気象数値予報で広く使われる変分型データ同化手法です。コスト関数の最小化と随伴法の理論をわかりやすく解説し、Pythonで3D-Varを実装します。 |
| **メインキーワード** | 3D-Var 4D-Var データ同化 |
| **関連キーワード** | 変分型データ同化, 3D-Var 数式, 4D-Var 随伴法, 気象予報 データ同化, コスト関数 最小化 |
| **想定文字数** | 5,500〜7,000字 |
| **Python可視化** | あり（コスト関数の収束・3D-Var解析場） |

---

## 記事構成

### リード文
```
気象庁・ECMWFなどの数値天気予報では、
毎日世界中の観測データを「変分型データ同化」で取り込んでいます。
3D-VarはOI法と同じゴール（コスト関数の最小化）を別の角度から解き、
4D-Varは時間次元も含めてシステム全体を最適化します。
逐次型（EnKF系）とは根本的に異なるアプローチを理解しましょう。
```

### この記事でわかること
- 3D-Varのコスト関数がOI法と数学的に等価であることが理解できる
- 4D-Varが随伴法を使って時間方向にも最適化する仕組みがわかる
- Pythonでシンプルなコスト関数最小化（3D-Var）を実装できる

---

### H2-1：変分型データ同化のアイデア

**キーワード**：変分型データ同化 とは, コスト関数 最小化

#### H3：OI法との関係

実はOI法のカルマンゲイン式は、次のコスト関数の最小化点と完全に一致します：

$$J(\mathbf{x}) = \underbrace{\frac{1}{2}(\mathbf{x}-\mathbf{x}^b)^T \mathbf{B}^{-1}(\mathbf{x}-\mathbf{x}^b)}_{J_b：背景場からのズレ} + \underbrace{\frac{1}{2}(\mathbf{y}-H\mathbf{x})^T \mathbf{R}^{-1}(\mathbf{y}-H\mathbf{x})}_{J_o：観測からのズレ}$$

- $J_b$：モデル予測を大きく変えるほどペナルティが増える
- $J_o$：観測との残差が大きいほどペナルティが増える
- $\mathbf{B}^{-1}$, $\mathbf{R}^{-1}$ が「どちらをどれだけ信頼するか」の重み

**「変分型」と呼ぶ理由**：$J(\mathbf{x})$ を最小化する $\mathbf{x}^a$ を変分法で求めるから。

#### H3：3D-VarとOI法の違いは実装にある

| 観点 | OI法 | 3D-Var |
|---|---|---|
| 数学的等価性 | ✓（同じ解） | ✓（同じ解） |
| Hの非線形性 | 線形のみ | 非線形OK（反復最小化） |
| $\mathbf{B}$の扱い | 明示的に行列 | 前処理子として暗黙的に適用 |
| 計算スケール | 小規模向け | 超大規模（$n \sim 10^8$）向け |

---

### H2-2：3D-Varの実装

**キーワード**：3D-Var 実装, コスト関数 勾配

#### H3：勾配の計算

コスト関数の勾配：
$$\nabla_{\mathbf{x}} J = \mathbf{B}^{-1}(\mathbf{x} - \mathbf{x}^b) - H^T \mathbf{R}^{-1}(\mathbf{y} - H\mathbf{x})$$

反復最小化（共役勾配法、L-BFGS）で解く：
1. 初期値 $\mathbf{x}^{(0)} = \mathbf{x}^b$ からスタート
2. 勾配 $\nabla J$ に沿って $\mathbf{x}$ を更新
3. 収束まで繰り返す

#### H3：変数変換（プレコンディショニング）

$\mathbf{x} = \mathbf{x}^b + \mathbf{B}^{1/2}\boldsymbol{\xi}$ と変換すると：
$$J(\boldsymbol{\xi}) = \frac{1}{2}\boldsymbol{\xi}^T\boldsymbol{\xi} + \frac{1}{2}(\mathbf{y} - H(\mathbf{x}^b + \mathbf{B}^{1/2}\boldsymbol{\xi}))^T \mathbf{R}^{-1}(\cdots)$$

→ 第1項が単位行列になり収束が速い（実用的な3D-Varの標準形）

---

### H2-3：4D-Varへの拡張

**キーワード**：4D-Var 随伴法, 4D-Var 気象予報

#### H3：3D-Varの限界と4D-Varのアイデア

3D-Varは「ある1時刻」の状態を推定する。
4D-Varは「時間窓 $[t_0, t_N]$ 内の全観測」を使って、初期状態 $\mathbf{x}_0$ を最適化する。

$$J(\mathbf{x}_0) = \frac{1}{2}(\mathbf{x}_0-\mathbf{x}^b_0)^T \mathbf{B}^{-1}(\mathbf{x}_0-\mathbf{x}^b_0) + \frac{1}{2}\sum_{k=0}^{N}(\mathbf{y}_k - H_k\mathbf{x}_k)^T \mathbf{R}_k^{-1}(\mathbf{y}_k - H_k\mathbf{x}_k)$$

ここで $\mathbf{x}_k = M_{k-1}(\mathbf{x}_{k-1})$（モデル積分で推進）。

#### H3：随伴法（Adjoint Method）

$J$ を $\mathbf{x}_0$ で微分するには、**時間を逆向きに積分する「随伴モデル」**が必要：

1. **順方向積分**：$\mathbf{x}_0 \to \mathbf{x}_1 \to \cdots \to \mathbf{x}_N$（観測残差を蓄積）
2. **逆方向積分**：随伴変数 $\boldsymbol{\lambda}$ を $t_N \to t_0$ に伝播させて $\nabla J$ を計算

**随伴モデルのコスト**：順方向モデルの3〜4倍の開発工数が必要
→ 大規模システムへの適用が難しい理由のひとつ

---

### H2-4：Python実装（3D-Var）

**キーワード**：3D-Var Python, 変分データ同化 実装

```python
# 3dvar_simple.py
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def cost_function(x, xb, y, H, B_inv, R_inv):
    """3D-Varコスト関数"""
    dx = x - xb
    innov = y - H @ x
    Jb = 0.5 * dx @ B_inv @ dx
    Jo = 0.5 * innov @ R_inv @ innov
    return Jb + Jo

def cost_gradient(x, xb, y, H, B_inv, R_inv):
    """コスト関数の勾配"""
    dx = x - xb
    innov = y - H @ x
    return B_inv @ dx - H.T @ R_inv @ innov

def run_3dvar(x_true, xb, obs_locs, obs_vals, sigma_b, L, sigma_r):
    n = len(xb)
    m = len(obs_vals)

    # ガウス型背景誤差共分散
    B = np.array([[sigma_b**2 * np.exp(-0.5*((i-j)/L)**2)
                   for j in range(n)] for i in range(n)])
    B_inv = np.linalg.inv(B + 1e-8 * np.eye(n))

    R = sigma_r**2 * np.eye(m)
    R_inv = np.linalg.inv(R)

    # 観測演算子（線形補間）
    H = np.zeros((m, n))
    x_grid = np.linspace(0, 1, n)
    for i, loc in enumerate(obs_locs):
        idx = np.searchsorted(x_grid, loc)
        idx = np.clip(idx, 1, n-1)
        alpha = (loc - x_grid[idx-1]) / (x_grid[idx] - x_grid[idx-1])
        H[i, idx-1] = 1 - alpha
        H[i, idx] = alpha

    # コスト関数をscipy.optimizeで最小化
    J_vals = []
    def callback(x):
        J_vals.append(cost_function(x, xb, obs_vals, H, B_inv, R_inv))

    result = minimize(
        cost_function, xb,
        jac=cost_gradient,
        args=(xb, obs_vals, H, B_inv, R_inv),
        method='L-BFGS-B',
        callback=callback,
        options={'maxiter': 100}
    )
    return result.x, J_vals

# --- 実行 ---
np.random.seed(42)
n = 100
x_grid = np.linspace(0, 1, n)
x_true = np.sin(2*np.pi*x_grid) + 0.5*np.sin(4*np.pi*x_grid)
xb = np.zeros(n)  # 背景場は0

obs_locs = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
obs_vals = np.interp(obs_locs, x_grid, x_true) + np.random.normal(0, 0.1, 5)

xa, J_vals = run_3dvar(x_true, xb, obs_locs, obs_vals,
                        sigma_b=1.0, L=0.15, sigma_r=0.1)

# プロット
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.plot(x_grid, x_true, 'k-', label='真値')
ax1.plot(x_grid, xb, 'b--', alpha=0.5, label='背景場 $x^b$')
ax1.plot(x_grid, xa, 'r-', lw=2, label='3D-Var解析場 $x^a$')
ax1.scatter(obs_locs, obs_vals, s=80, color='green', zorder=5, label='観測')
ax1.legend(); ax1.grid(alpha=0.3)
ax1.set_title('3D-Var 解析結果')

ax2.semilogy(J_vals, 'b-o', markersize=4)
ax2.set_xlabel('反復回数'); ax2.set_ylabel('コスト関数 J')
ax2.set_title('コスト関数の収束')
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('3dvar_result.png', dpi=150)
plt.show()
```

---

### H2-5：EnKFと4D-Varの比較・使い分け

**キーワード**：EnKF 4D-Var 比較, En-4DVar ハイブリッド

| 観点 | EnKF/LETKF | 4D-Var |
|---|---|---|
| 誤差共分散 | アンサンブルで流動的に推定 | 固定 $\mathbf{B}$（4D-Varの弱点） |
| 非線形 | 近似対応 | 随伴モデルで対応 |
| 開発コスト | 低（モデルをブラックボックス扱い） | 高（随伴モデルが必要） |
| 時間窓の最適化 | なし | あり（観測を過去〜未来から活用） |
| 現在のトレンド | ECMWF等でEnKF+4D-VarのHybridが主流 | |

**ハイブリッド手法（En-4DVar）**：
- 4D-VarのBをアンサンブルで推定（流動的な背景誤差共分散）
- ECMWF・気象庁などで採用

---

### FAQ セクション

```
Q1: 3D-VarはOI法と結果が同じですか？
A: Hが線形かつBが同じ行列の場合は数学的に同じ解になります。
   3D-Varの利点は非線形Hへの対応と、
   大規模問題でBを陽に持たない実装が可能な点です。

Q2: 随伴モデルとは何ですか？
A: 順方向モデルの転置（随伴）を使って勾配を効率的に計算する技術です。
   数学的にはモデル方程式を時間逆方向に解くことに相当します。
   開発が複雑で、4D-Varが普及しにくい原因のひとつです。

Q3: 変分型とアンサンブル型のどちらを選ぶべきですか？
A: 現在の実務ではハイブリッド（En-4DVar）が最も精度が高いとされています。
   コンピュータリソースと随伴モデル開発コストが許すなら4D-Var系、
   そうでなければLETKFが現実的な選択肢です。

Q4: CFD（OpenFOAM）への変分型データ同化の適用は？
A: 可能ですが随伴CFDソルバーの開発が必要です。
   OpenFOAMにはAdjoint最適化機能があり、
   これをデータ同化に応用した研究が進んでいます。
```

---

## Pythonファイル

- `python/3dvar_simple.py`：3D-Varの実装（scipy.optimizeによる最小化）
- `python/cost_function_visualization.py`：2次元コスト関数の等高線プロット
