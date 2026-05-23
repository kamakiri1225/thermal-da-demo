# 離散化ガイド：前進オイラー法

連続時間の微分方程式をコンピュータで扱える差分方程式に変換する方法を説明する。
本プログラムでは**前進オイラー法**を使用している。

---

## 1. 「離散化」とは何か

物理の方程式（熱拡散など）は**連続時間**で書かれる：

$$\frac{d\boldsymbol{\theta}}{dt} = \mathbf{A}_c\,\boldsymbol{\theta} + \mathbf{B}_c\,u$$

コンピュータは**離散時間**（= 一定時間 Δt ごとのステップ）でしか動けない。

```
連続時間の世界:
  θ(t) ─────────────────────────────→ 時間 t
         なめらかに変化（微分方程式）

コンピュータの世界:
  θ₀, θ₁, θ₂, θ₃, ...  ← 飛び飛びの時刻の値（数列）
   ↑Δt↑Δt↑Δt

「離散化」= 微分方程式 → 差分方程式 への変換
```

---

## 2. 前進オイラー法の考え方

### 微分を「前進差分」で近似する

$$\frac{d\boldsymbol{\theta}}{dt}\bigg|_{t=k\Delta t} \approx \frac{\boldsymbol{\theta}_{k+1} - \boldsymbol{\theta}_k}{\Delta t}$$

「$\Delta t$ 秒前の傾き（= 現在の変化率）を使って、次の時刻の値を予測する」という考え方。

```
θ(t)
 │        ● θ_{k+1}（← 予測したい）
 │       ╱
傾き╱  ← dθ/dt ≈ A_c θ_k + B_c u_k を使う
 │     ●
 │   θ_k（← 現在の値）
 ├───┼───┼──→ t
    k   k+1
    ←Δt→
```

---

## 3. 離散化行列 A_d, B_d の導出

状態方程式 $\dot{\boldsymbol{\theta}} = \mathbf{A}_c\boldsymbol{\theta} + \mathbf{B}_c u$ に前進差分を代入：

$$\frac{\boldsymbol{\theta}_{k+1} - \boldsymbol{\theta}_k}{\Delta t} = \mathbf{A}_c\boldsymbol{\theta}_k + \mathbf{B}_c u_k$$

両辺に $\Delta t$ をかけて整理：

$$\boldsymbol{\theta}_{k+1} = \underbrace{(\mathbf{I} + \mathbf{A}_c\Delta t)}_{\mathbf{A}_d}\,\boldsymbol{\theta}_k + \underbrace{\mathbf{B}_c\Delta t}_{\mathbf{B}_d}\,u_k$$

$$\boxed{\mathbf{A}_d = \mathbf{I} + \mathbf{A}_c\,\Delta t, \qquad \mathbf{B}_d = \mathbf{B}_c\,\Delta t}$$

---

## 4. コードとの対応（`thermal_model.py`）

```python
def _discretize(self):
    # A_d = I + A_c·Δt
    A_d = np.eye(self.n) + self.A_c * self.dt

    # B_d = B_c·Δt
    B_d = self.B_c * self.dt

    return A_d, B_d
```

`numpy` だけで実装できる。外部ライブラリは不要。

---

## 5. 安定条件

前進オイラー法には安定条件がある。
スカラー版（$\dot{\theta} = -a\theta$、$a > 0$）の場合：

$$A_d = 1 - a\Delta t \quad \text{が安定} \iff |A_d| < 1 \iff 0 < a\Delta t < 2$$

行列版では、$\mathbf{A}_d$ の全固有値の絶対値が1未満であることが必要。

**本プログラムでの確認：**

$$r \cdot \Delta t = \frac{\alpha}{(\Delta x)^2} \cdot \Delta t = \frac{1.2\times10^{-5}}{0.1^2} \times 10 = 0.012 \ll 0.5 \quad \checkmark$$

Δt=10s・dx=0.1m の設定では十分安定。

---

## 6. 主な離散化手法の比較

| 手法 | $A_d$ の計算 | 安定条件 | ライブラリ |
|---|---|---|---|
| **前進オイラー（本プログラム）** | $I + A_c\Delta t$ | $r\Delta t < 0.5$ | numpy のみ |
| 後退オイラー | $(I - A_c\Delta t)^{-1}$ | 無条件安定 | numpy（逆行列） |
| ZOH（行列指数関数） | $e^{A_c\Delta t}$ | 無条件安定 | scipy.linalg.expm |

本プログラムの条件（$r\Delta t = 0.012$）ではどの手法でも結果はほぼ同じになる。
前進オイラーはコードが最もシンプルなため採用している。
