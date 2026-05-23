# 理論解説：状態空間モデルの導出とカルマンフィルタの位置付け

---

## 目次

1. [熱方程式から状態空間モデルへ（A・B行列の導出）](#1-熱方程式から状態空間モデルへ)
   - 1.1 出発点：熱伝導方程式（PDE）
   - 1.2 空間離散化：有限差分法（FDM）
   - 1.3 内部節点の方程式
   - 1.4 境界節点の方程式（対流境界条件）
   - 1.5 行列形式への整理（連続時間状態空間モデル）
   - 1.6 観測方程式と観測行列 H
   - 1.7 前進オイラー法による離散化：A_d・B_d の計算（①なぜ必要か → ②前進差分の考え方 → ③離散化式の導出 → ④安定条件の確認 → ⑤実装）
2. [カルマンフィルタの理論的位置付け](#2-カルマンフィルタの理論的位置付け)
   - 2.1 データ同化とは何か
   - 2.2 ベイズ推定としての解釈
   - 2.3 最小分散推定（MMSE）としての導出
   - 2.4 カルマンゲインの物理的意味
   - 2.5 他の手法との位置関係
   - 2.6 工作機械熱誤差補償への応用

---

## 1. 熱方程式から状態空間モデルへ

### 1.1 出発点：熱伝導方程式（PDE）

1次元の熱棒における温度 $T(x,\,t)$ の時間発展は、**熱伝導方程式**（放物型偏微分方程式）で支配される：

$$
\rho c_p \frac{\partial T}{\partial t} = k \frac{\partial^2 T}{\partial x^2} + \dot{q}(x,\,t)
$$

| 記号 | 単位 | 意味 |
|---|---|---|
| $T(x,t)$ | °C | 温度（位置・時刻の関数）|
| $\rho$ | kg/m³ | 密度 |
| $c_p$ | J/(kg·K) | 比熱容量 |
| $k$ | W/(m·K) | 熱伝導率 |
| $\dot{q}$ | W/m³ | 体積発熱率（熱源）|

熱拡散率 $\alpha$ を導入すると：

$$
\alpha = \frac{k}{\rho c_p} \quad [\text{m}^2/\text{s}]
$$

PDE は次のように書き直せる：

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2} + \frac{\dot{q}}{\rho c_p}
$$

> **物理的意味**：左辺 $\partial T/\partial t$ は「温度の時間的変化率」。右辺第1項は「空間的な温度の曲率に比例した熱拡散」、第2項は「熱源による温度上昇」を表す。

---

### 1.2 空間離散化：有限差分法（FDM）

PDE をコンピュータで解くために、空間方向 $x$ を有限個の**節点**（node）に離散化する。

棒の全長 $L$ を $n-1$ 等分し、節点間距離を $\Delta x = L/(n-1)$ とする：

$$
x_i = i \cdot \Delta x, \quad i = 0,\, 1,\, \ldots,\, n-1
$$

```
x=0                               x=L
 |     Δx    |     Δx    |     Δx    |
[0] ─────── [1] ─────── [2] ─────── [3] ── ··· ── [n-1]
 ↑                                                    ↑
左端（熱源・対流）                               右端（対流）
```

各節点の温度を $T_i(t) \equiv T(x_i,\, t)$ と表記する。

**2階微分の中心差分近似**（テイラー展開から導出）：

$$
\frac{\partial^2 T}{\partial x^2}\bigg|_{x=x_i} \approx \frac{T_{i-1}(t) - 2T_i(t) + T_{i+1}(t)}{(\Delta x)^2} + O((\Delta x)^2)
$$

> **導出**：テイラー展開より
> $T_{i+1} = T_i + \Delta x \frac{\partial T}{\partial x} + \frac{(\Delta x)^2}{2}\frac{\partial^2 T}{\partial x^2} + \cdots$
> $T_{i-1} = T_i - \Delta x \frac{\partial T}{\partial x} + \frac{(\Delta x)^2}{2}\frac{\partial^2 T}{\partial x^2} + \cdots$
> 辺々加算して $\frac{\partial^2 T}{\partial x^2} \approx (T_{i-1} - 2T_i + T_{i+1})/(\Delta x)^2$

---

### 1.3 内部節点の方程式（$i = 1, 2, \ldots, n-2$）

内部節点では両隣の節点からの熱伝導のみが働く：

$$
\frac{dT_i}{dt} = \frac{\alpha}{(\Delta x)^2}\left(T_{i-1} - 2T_i + T_{i+1}\right)
$$

ここで $r = \alpha / (\Delta x)^2$ と置くと：

$$
\boxed{\frac{dT_i}{dt} = r\,T_{i-1} - 2r\,T_i + r\,T_{i+1}}
$$

$r$ は**拡散係数**（数値安定性の指標 $r \cdot \Delta t \leq 0.5$ が陽的差分の安定条件）。

本プログラムの値：$r = 1.2\times10^{-5} / (0.1)^2 = 1.2\times10^{-3}$ s⁻¹

---

### 1.4 境界節点の方程式（対流境界条件）

両端節点では、隣接節点からの**熱伝導**に加えて、周囲空気との**対流熱伝達**が働く。

**Newton の冷却則**（対流）：

$$
q_\text{conv} = h_\text{conv}(T_i - T_\text{amb}) \quad [\text{W/m}^2]
$$

#### 節点 $0$（左端：熱源 + 対流）

エネルギーバランス（制御体積 $\Delta x$ に対して）：

$$
\rho c_p \Delta x \frac{dT_0}{dt} = \underbrace{k \cdot \frac{T_1 - T_0}{\Delta x}}_{\text{右隣からの熱伝導}} - \underbrace{h_\text{conv}(T_0 - T_\text{amb})}_{\text{対流冷却}} + \underbrace{\dot{q}_\text{in}}_{\text{熱源}}
$$

$T_\text{amb} = 0$（周囲温度を基準としてゼロ設定）として整理：

$$
\frac{dT_0}{dt} = \underbrace{\frac{\alpha}{(\Delta x)^2}}_{r}(T_1 - T_0) - \underbrace{\frac{h_\text{conv}}{\rho c_p \Delta x}}_{h}\,T_0 + \underbrace{\frac{1}{\rho c_p \Delta x}}_{b}\,u
$$

$$
\boxed{\frac{dT_0}{dt} = -(r + h)\,T_0 + r\,T_1 + b\,u}
$$

| 係数 | 定義 | 本プログラムの値 |
|---|---|---|
| $r = \alpha/(\Delta x)^2$ | 拡散係数 [1/s] | $1.2\times10^{-3}$ |
| $h = h_\text{conv}/(\rho c_p \Delta x)$ | 対流係数 [1/s] | $2.3\times10^{-5}$ |
| $b = 1/(\rho c_p \Delta x)$ | 入力係数 [m³·K/J] | $2.86\times10^{-7}$ |

**節点 $n-1$（右端：対流のみ）** も同様：

$$
\frac{dT_{n-1}}{dt} = r\,T_{n-2} - (r + h)\,T_{n-1}
$$

---

### 1.5 行列形式への整理（連続時間状態空間モデル）

全節点の温度をベクトルにまとめる：

$$
\boldsymbol{\theta}(t) = \begin{bmatrix} T_0(t) \\ T_1(t) \\ \cdots \\ T_9(t) \end{bmatrix} \in \mathbb{R}^{10}
$$

全節点の微分方程式を行列形式で書くと：

$$
\dot{\boldsymbol{\theta}} = \mathbf{A}_c\,\boldsymbol{\theta} + \mathbf{B}_c\,u
$$

#### A_c 行列（$n=10$ の場合）

$$
\mathbf{A}_c = \begin{bmatrix}
-(r+h) & r      & 0      & 0      & 0      \\
r      & -2r    & r      & 0      & 0      \\
0      & r      & -2r    & r      & 0      \\
0      & 0      & r      & -2r    & r      \\
0      & 0      & 0      & r      & -(r+h)
\end{bmatrix}
$$

**構造の特徴**：
- **三重対角行列**（tridiagonal）：隣接節点のみが結合
- **対角要素が負**：各節点は熱を「失う」（散逸系）
- **対角外要素が正**：隣接節点から熱が「来る」
- **両端の対角要素が内部節点より小さい**：対流冷却の追加損失 $-h$

> **安定性**：$\mathbf{A}_c$ は負定値行列（全固有値 $< 0$）であり、入力なし・断熱境界なら温度は必ず周囲温度に収束する。

#### B_c 行列

熱源は節点 $0$ のみに集中：

$$
\mathbf{B}_c = \begin{bmatrix} b \\ 0 \\ 0 \\ 0 \\ 0 \end{bmatrix} = \begin{bmatrix} 1/(\rho c_p \Delta x) \\ 0 \\ 0 \\ 0 \\ 0 \end{bmatrix}
$$

入力 $u$ は熱入力 $Q$ [W/m]（単位奥行き当たりの線熱流束）。

---

### 1.6 観測方程式と観測行列 H

センサは節点 $0$ と $9$ のみに配置されているとする。観測方程式：

$$
\mathbf{y}_k = \mathbf{H}\,\boldsymbol{\theta}_k + \mathbf{v}_k, \quad \mathbf{v}_k \sim \mathcal{N}(\mathbf{0}, \mathbf{R})
$$

$$
\mathbf{H} = \begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
\end{bmatrix}
$$

$\mathbf{H}$ は「どの節点を観測するか」を選択する行列。$\mathbf{H} \in \mathbb{R}^{p \times n}$（$p$: センサ数、$n$: 節点数）。

> **可観測性条件**：$\mathbf{H}$ と $\mathbf{A}_c$（または $\mathbf{A}_d$）の組が可観測であること。
> 可観測行列 $\mathcal{O}$ がフルランクであれば、有限個のセンサで全状態が推定可能：
> $$\text{rank}(\mathcal{O}) = \text{rank}\begin{pmatrix}\mathbf{H} \\ \mathbf{H}\mathbf{A}_d \\ \vdots \\ \mathbf{H}\mathbf{A}_d^{n-1}\end{pmatrix} = n$$

---

### 1.7 前進オイラー法による離散化：A_d・B_d の計算

この節では「なぜ離散化が必要か」という動機から始め、前進差分の考え方を使って離散化行列を導出し、最後に安定条件の確認と実装を説明する。

---

#### ステップ① なぜ離散化が必要か

熱系の物理モデルは連続時間の微分方程式（ODE）で書かれている：

$$
\dot{\boldsymbol{\theta}}(t) = \mathbf{A}_c\,\boldsymbol{\theta}(t) + \mathbf{B}_c\,u(t)
$$

しかし現実のセンサとコンピュータは**離散的な時刻**でしか動作しない。

```
物理（連続時間）          センサ・カルマンフィルタ（離散時間）
──────────────           ────────────────────────────────────
θ(t) が滑らかに変化        θ_0, θ_1, θ_2, … の列として扱う
                              ↑        ↑        ↑
                           t=0    t=Δt   t=2Δt  ← Δt ごとにサンプル
```

目標：「時刻 $k\cdot\Delta t$ の状態 $\boldsymbol{\theta}_k$ から、次の時刻 $(k+1)\cdot\Delta t$ の状態 $\boldsymbol{\theta}_{k+1}$ を直接計算する漸化式」を求めること。

$$
\boldsymbol{\theta}_{k+1} = \mathbf{A}_d\,\boldsymbol{\theta}_k + \mathbf{B}_d\,u_k
$$

---

#### ステップ② 前進差分の考え方

微分 $\dot{\boldsymbol{\theta}} = d\boldsymbol{\theta}/dt$ を、**前進差分**で近似する：

$$
\frac{d\boldsymbol{\theta}}{dt}\bigg|_{t=k\Delta t}
\approx \frac{\boldsymbol{\theta}_{k+1} - \boldsymbol{\theta}_k}{\Delta t}
$$

```
θ(t)
 │        ●  θ_{k+1}
 │       ╱
 │      ╱  ← この傾き ≈ dθ/dt を使って θ_{k+1} を予測
 │     ●
 │   θ_k
 │
 ├───┼───┼──→ t
    k   k+1
    ←Δt→
```

これを状態方程式 $\dot{\boldsymbol{\theta}} = \mathbf{A}_c\boldsymbol{\theta} + \mathbf{B}_c u$ に代入すると：

$$
\frac{\boldsymbol{\theta}_{k+1} - \boldsymbol{\theta}_k}{\Delta t}
= \mathbf{A}_c\,\boldsymbol{\theta}_k + \mathbf{B}_c\,u_k
$$

---

#### ステップ③ 離散化式の導出

両辺に $\Delta t$ をかけて $\boldsymbol{\theta}_{k+1}$ を左辺に移項する：

$$
\boldsymbol{\theta}_{k+1} = \boldsymbol{\theta}_k + \Delta t\,(\mathbf{A}_c\,\boldsymbol{\theta}_k + \mathbf{B}_c\,u_k)
$$

$$
= (\mathbf{I} + \mathbf{A}_c\,\Delta t)\,\boldsymbol{\theta}_k + \mathbf{B}_c\,\Delta t\,u_k
$$

これにより離散化行列が得られる：

$$
\boxed{
\mathbf{A}_d = \mathbf{I} + \mathbf{A}_c\,\Delta t,
\qquad
\mathbf{B}_d = \mathbf{B}_c\,\Delta t
}
$$

**意味**：$\mathbf{I}$ は「現在の温度をそのまま持ち越す」部分、$\mathbf{A}_c\,\Delta t$ は「$\Delta t$ 秒間の熱拡散による温度変化」を加える部分。

スカラー例（$\dot{\theta} = -a\theta + bu$）で確認：

| 連続系 | 離散化後 |
|---|---|
| $a$ | $A_d = 1 - a\Delta t$ |
| $b$ | $B_d = b\,\Delta t$ |

---

#### ステップ④ 安定条件の確認

前進オイラー法には安定条件がある。スカラー版 $A_d = 1 - a\Delta t$ が安定には $|A_d| < 1$ が必要：

$$
|1 - a\Delta t| < 1 \quad \Longrightarrow \quad 0 < a\Delta t < 2
$$

行列版では、$\mathbf{A}_d = \mathbf{I} + \mathbf{A}_c\,\Delta t$ の全固有値の絶対値が1未満であること：

$$
|\lambda_i(\mathbf{A}_d)| = |1 + \lambda_i(\mathbf{A}_c)\,\Delta t| < 1
$$

| 離散化方法 | 安定条件 | 備考 |
|---|---|---|
| **前進オイラー法（陽的）** | $r\cdot\Delta t \leq 0.5$ | **本プログラムで使用** |
| 後退オイラー法（陰的） | 無条件安定 | 計算が複雑（逆行列が必要） |
| ZOH（行列指数関数） | 無条件安定 | `scipy.linalg.expm` が必要 |

本プログラムの安定性確認：

$$
r \cdot \Delta t = \frac{\alpha}{(\Delta x)^2}\,\Delta t
= \frac{1.2\times10^{-5}}{(0.1)^2} \times 10
= 1.2\times10^{-3} \times 10 = 0.012 \ll 0.5 \quad \checkmark
$$

Δt=10s・dx=0.1m の条件では十分安定であることが確認できる。

---

#### ステップ⑤ 実装（Python コード）

前進オイラー法の実装は行列の加算・乗算だけで済む：

```python
import numpy as np

n = 10  # 状態数（節点数）

# A_d = I + A_c·Δt
A_d = np.eye(n) + A_c * dt

# B_d = B_c·Δt
B_d = B_c * dt
```

`scipy` などの外部ライブラリは不要。`numpy` だけで計算できる。

---

## 2. カルマンフィルタの理論的位置付け

### 2.1 データ同化とは何か

**データ同化（Data Assimilation）** とは、**物理モデルの予測**と**実際の観測データ**を統計的に最適な方法で融合し、システムの状態をより正確に推定する手法である。

```
物理モデル（第一原理）             観測データ（センサ計測）
  ・原理的に正確                    ・直接的な情報
  ・モデル誤差・パラメータ誤差あり   ・ノイズ・欠損・空間的疎らさ
           ↓                                  ↓
           └─────────── データ同化 ────────────┘
                              ↓
                     最良状態推定（Best Estimate）
```

データ同化は元々**気象予報**（大気・海洋の数値シミュレーション）で発展した。現在は工学・医学・経済学など広範囲に応用されている。

**本研究での位置付け**：FEM デジタルツイン（物理モデル）＋温度センサ計測値（観測）→ 工作機械の全温度場推定

---

### 2.2 ベイズ推定としての解釈

カルマンフィルタは **ベイズ推定（Bayesian Estimation）** の特殊ケースとして厳密に導出できる。

#### ベイズの定理

$$
\underbrace{p(\boldsymbol{\theta}_k \mid \mathbf{y}_{0:k})}_{\text{事後分布}}
\propto
\underbrace{p(\mathbf{y}_k \mid \boldsymbol{\theta}_k)}_{\text{尤度}} \times
\underbrace{p(\boldsymbol{\theta}_k \mid \mathbf{y}_{0:k-1})}_{\text{事前分布}}
$$

| 項 | 意味 | カルマンフィルタでの表現 |
|---|---|---|
| 事前分布 | モデル予測から得られる信念 | $\mathcal{N}(\hat{\boldsymbol{\theta}}_k^-, \mathbf{P}_k^-)$ |
| 尤度 | センサ計測がその状態から生成される確率 | $\mathcal{N}(\mathbf{H}\hat{\boldsymbol{\theta}}_k^-, \mathbf{R})$ |
| 事後分布 | 計測後に更新された信念 | $\mathcal{N}(\hat{\boldsymbol{\theta}}_k, \mathbf{P}_k)$ |

**ガウス分布の積はガウス分布**であることを利用すると、更新式は解析的に解ける：

$$
p(\boldsymbol{\theta}_k \mid \mathbf{y}_{0:k})
= \mathcal{N}\!\left(\hat{\boldsymbol{\theta}}_k,\, \mathbf{P}_k\right)
$$

$$
\hat{\boldsymbol{\theta}}_k = \hat{\boldsymbol{\theta}}_k^- + \mathbf{K}_k\left(\mathbf{y}_k - \mathbf{H}\hat{\boldsymbol{\theta}}_k^-\right)
$$

$$
\mathbf{P}_k = \left(\mathbf{I} - \mathbf{K}_k\mathbf{H}\right)\mathbf{P}_k^-
$$

これがカルマンフィルタの更新式そのものである。つまり、**カルマンフィルタは線形ガウス系に対する「ベイズ推定の閉形式解」**である。

---

### 2.3 最小分散推定（MMSE）としての導出

カルマンゲイン $\mathbf{K}_k$ を「推定誤差の分散を最小化する」条件から導出する。

**目標**：推定誤差の共分散行列のトレース（= 平均二乗誤差の和）を最小化する：

$$
\min_{\mathbf{K}_k} \;\mathrm{tr}\!\left(\mathbf{P}_k\right)
= \min_{\mathbf{K}_k} \;\mathrm{tr}\!\left[\left(\mathbf{I} - \mathbf{K}_k\mathbf{H}\right)\mathbf{P}_k^-\left(\mathbf{I} - \mathbf{K}_k\mathbf{H}\right)^\top + \mathbf{K}_k\mathbf{R}\mathbf{K}_k^\top\right]
$$

$\mathbf{K}_k$ で微分してゼロとおく（$\partial\,\mathrm{tr}(\mathbf{P}_k)/\partial\mathbf{K}_k = \mathbf{0}$）：

$$
-2\left(\mathbf{I} - \mathbf{K}_k\mathbf{H}\right)\mathbf{P}_k^-\mathbf{H}^\top + 2\mathbf{K}_k\mathbf{R} = \mathbf{0}
$$

$$
\mathbf{K}_k\left(\mathbf{H}\mathbf{P}_k^-\mathbf{H}^\top + \mathbf{R}\right) = \mathbf{P}_k^-\mathbf{H}^\top
$$

$$
\boxed{
\mathbf{K}_k = \mathbf{P}_k^-\mathbf{H}^\top\underbrace{\left(\mathbf{H}\mathbf{P}_k^-\mathbf{H}^\top + \mathbf{R}\right)^{-1}}_{\text{イノベーション共分散}^{-1}}
}
$$

これがカルマンゲインの式である。この $\mathbf{K}_k$ は**線形不偏推定量の中で分散最小（MVUE: Minimum Variance Unbiased Estimator）**であることが証明されている。

---

### 2.4 カルマンゲインの物理的意味

カルマンゲイン $\mathbf{K}_k$ は「モデル予測とセンサ計測の信頼度の比」で自動的に決まる重みである。

$$
\hat{\boldsymbol{\theta}}_k
= \underbrace{\hat{\boldsymbol{\theta}}_k^-}_{\text{モデル予測}} + \mathbf{K}_k\underbrace{\left(\mathbf{y}_k - \mathbf{H}\hat{\boldsymbol{\theta}}_k^-\right)}_{\text{イノベーション（予測残差）}}
$$

#### 2つの極端なケース

**ケース1：センサが非常に正確（$\mathbf{R} \to \mathbf{0}$）**

$$
\mathbf{K}_k \to \mathbf{H}^+ \quad (\mathbf{H} \text{ の擬似逆行列})
\implies \hat{\boldsymbol{\theta}}_k \approx \mathbf{H}^+\mathbf{y}_k \quad \text{（センサを全信頼）}
$$

**ケース2：モデルが非常に正確（$\mathbf{P}_k^- \to \mathbf{0}$）**

$$
\mathbf{K}_k \to \mathbf{0}
\implies \hat{\boldsymbol{\theta}}_k \approx \hat{\boldsymbol{\theta}}_k^- \quad \text{（モデル予測を全信頼）}
$$

**一般的なケース（最適バランス）**：

```
モデル不確かさ P が大きい ──→ K が大きい ──→ センサを重視
センサノイズ R が大きい   ──→ K が小さい ──→ モデルを重視
```

> **直感的な理解**：GPS と地図を使ったナビゲーションに例えると、
> - GPS 信号が弱いとき（$\mathbf{R}$ 大）→ 地図（モデル）を重視
> - 地図が古いとき（$\mathbf{P}^-$ 大）→ GPS（センサ）を重視
> - カルマンフィルタはこの「どちらをどれだけ信頼するか」を最適に計算する

---

### 2.5 他の手法との位置関係

#### データ同化手法の分類

```
データ同化
├── 逐次型（Sequential）
│   ├── カルマンフィルタ（KF）           ← 本プログラム
│   │   ├── 拡張カルマンフィルタ（EKF）
│   │   ├── アンサンテッドKF（UKF）
│   │   └── アンサンブルKF（EnKF）
│   └── 粒子フィルタ（Particle Filter）
└── 変分型（Variational）
    ├── 3D-Var
    └── 4D-Var
```

#### 各手法の比較

| 手法 | 対象系 | 計算コスト | 必要条件 | 工作機械への適合性 |
|---|---|---|---|---|
| **標準 KF** | 線形・ガウス | 低 | 線形性 | ○（小規模・線形なら最適）|
| EKF | 弱非線形 | 中 | ヤコビアン計算可 | ◎（温度依存材料特性に対応）|
| UKF | 中程度非線形 | 中 | — | ◎（高精度・コード簡潔）|
| EnKF | 強非線形・高次元 | 高 | アンサンブル数 | △（FEMの全状態にそのまま適用困難）|
| 粒子フィルタ | 任意 | 非常に高 | — | △（リアルタイムに不向き）|

#### ルエンベルガーオブザーバとの対比

制御工学では**ルエンベルガーオブザーバ**（状態推定器）が知られている：

$$
\hat{\boldsymbol{\theta}}_{k+1}
= \mathbf{A}_d\hat{\boldsymbol{\theta}}_k + \mathbf{B}_d u_k
+ \underbrace{\mathbf{L}}_{\text{ゲイン（設計者が決める）}}\!\left(\mathbf{y}_k - \mathbf{H}\hat{\boldsymbol{\theta}}_k\right)
$$

カルマンフィルタはルエンベルガーオブザーバの特殊ケースであり、ゲイン $\mathbf{L} = \mathbf{K}_k$ が**ノイズ統計量 $\mathbf{Q}, \mathbf{R}$ から最適に計算される点が異なる**。

#### 変分型データ同化（4D-Var）との等価性

時間窓が無限に長い場合、カルマンフィルタの解は4D-Var（最適制御問題として定式化するデータ同化）の解と等価になることが知られている。

---

### 2.6 工作機械熱誤差補償への応用

本プログラムが模擬する論文（Lang et al., 2024）との対応を整理する。

#### 役割分担

```
FEM デジタルツイン                   カルマンフィルタ
─────────────────────               ────────────────────────────
「熱的にこう動くはずだ」              「モデルとセンサ両方から
  という物理的な予測                   最良の状態推定を計算する」
      ↓                                       ↑
  A_d, B_d を提供                      y_k を受け取る
      └────────────────────────────────┘
              センサ計測 y_k（実機から）
```

#### 本手法の優位性

| 課題 | 従来のデータ駆動手法 | カルマンフィルタ + FEM |
|---|---|---|
| 訓練データ | 大量の実験データが必要 | **不要**（FEM が代替）|
| センサ数 | 多数必要（16本など）| **少数可能**（4本で83%補償）|
| 汎化性 | 学習条件外で劣化 | FEM の物理的制約により頑健 |
| リアルタイム更新 | 困難（再学習必要）| 自然に対応（逐次更新）|

#### 今後の拡張方向

| テーマ | 内容 |
|---|---|
| **適応 KF** | $\mathbf{Q}$ をオンラインで更新し、FEM モデル誤差を吸収 |
| **EKF への拡張** | 温度依存の材料特性（非線形 $\mathbf{A}_c$）に対応 |
| **センサ配置最適化** | 可観測性グラミアン $\mathcal{W}_o = \sum_{k=0}^{\infty} (\mathbf{A}_d^k)^\top \mathbf{H}^\top \mathbf{H} \mathbf{A}_d^k$ を最大化する $\mathbf{H}$ を設計 |
| **OpenFOAM 連成** | CFD（流体・熱連成）解析から $\mathbf{A}_c, \mathbf{B}_c$ を抽出 |

---

## まとめ：A・B 行列とカルマンフィルタの関係

$$
\underbrace{
  \text{熱方程式 (PDE)}
  \xrightarrow{\text{FDM 空間離散化}}
  \dot{\boldsymbol{\theta}} = \mathbf{A}_c\boldsymbol{\theta} + \mathbf{B}_c u
  \xrightarrow{\text{前進オイラー離散化}}
  \boldsymbol{\theta}_{k+1} = \mathbf{A}_d\boldsymbol{\theta}_k + \mathbf{B}_d u_k
}_{\text{FEM デジタルツイン（物理モデル）}}
$$

$$
\downarrow \quad \text{カルマンフィルタに入力}
$$

$$
\underbrace{
  \hat{\boldsymbol{\theta}}_k = \hat{\boldsymbol{\theta}}_k^- + \mathbf{K}_k(\mathbf{y}_k - \mathbf{H}\hat{\boldsymbol{\theta}}_k^-)
  \xrightarrow{\delta = \alpha_\text{exp}\sum_i\hat{\theta}_i\Delta x}
  \hat{\delta}_k
}_{\text{データ同化 → 熱変位の推定（補償値の計算）}}
$$

A・B 行列はデジタルツインが保持する「物理知識の結晶」であり、カルマンフィルタはその知識とセンサ情報を最適に融合する「推論エンジン」である。
