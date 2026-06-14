# データ同化手法の選定

## 問題設定

| 項目 | 値 |
|------|-----|
| 状態次元 | $n = 2400$（40×60 グリッド） |
| 観測数 | $m \approx 20$〜$50$（スパースな点観測） |
| モデル | 静的な正解場（時間発展なし） |
| 非線形性 | 線形（$H$ = 最近傍サンプリング） |
| 目標 | スパース観測から2次元ドット絵を再構築 |

---

## 手法の比較

### OI（最適内挿法）

**原理**：一回限りの線形更新

$$\mathbf{x}^a = \mathbf{x}^b + \mathbf{K}(\mathbf{y} - H\mathbf{x}^b)$$

$$\mathbf{K} = \mathbf{B}H^\top(H\mathbf{B}H^\top + \mathbf{R})^{-1}$$

| 長所 | 短所 |
|------|------|
| シンプル、解析解が存在する | $\mathbf{B}$ を手動で指定（例：ガウス相関）する必要がある |
| 1回の更新で完結 | $\mathbf{B}$ は更新されない → サイクルで改善できない |
| 一発補正に適している | 反復的な画像再構築には不向き |

**このデモへの適用可否**：1ステップの補正には使えるが、
「おっさんが徐々に浮かび上がる」多サイクル演出ができない。

---

### カルマンフィルタ（KF）

**原理**：モデル $\mathbf{x}_{k+1} = F\mathbf{x}_k + \text{noise}$ を持つ最適逐次推定法

**予測ステップ：**

$$\mathbf{x}^f_k = F\mathbf{x}^a_{k-1}$$

$$P^f_k = FP^a_{k-1}F^\top + Q$$

**更新ステップ：**

$$\mathbf{K}_k = P^f_k H^\top(HP^f_k H^\top + R)^{-1}$$

$$\mathbf{x}^a_k = \mathbf{x}^f_k + \mathbf{K}_k(\mathbf{y}_k - H\mathbf{x}^f_k)$$

$$P^a_k = (I - \mathbf{K}_k H) P^f_k$$

| 長所 | 短所 |
|------|------|
| 線形ガウス系に対して最適 | 陽なモデル行列 $F$ が必要 |
| 誤差共分散を解析的に伝播させる | $P$ の格納コストが $n \times n = 2400 \times 2400$（float64 で約 46 MB） |

**このデモへの適用可否**：静的な正解場では $F = I$（恒等行列）とすることでKFはOIの逐次版に帰着するが、
$P$ の格納・逆行列演算コストが $n=2400$ では非現実的。

---

### EnKF（アンサンブルカルマンフィルタ）← **推奨**

**原理**：$N$ 個のアンサンブルメンバーで $\mathbf{B}$ を近似

$$\text{アンサンブル：} \{\mathbf{x}^{(1)}, \ldots, \mathbf{x}^{(N)}\}$$

$$\mathbf{B} \approx \frac{1}{N-1} \mathbf{X}' (\mathbf{X}')^\top, \quad \mathbf{X}' = \mathbf{X} - \bar{\mathbf{x}}\mathbf{1}^\top$$

ゲイン行列：

$$\mathbf{K} = \mathbf{B}H^\top(H\mathbf{B}H^\top + \mathbf{R})^{-1}$$

各メンバーの更新（確率的 EnKF）：

$$\mathbf{x}^{a,(i)} = \mathbf{x}^{f,(i)} + \mathbf{K}\!\left(\mathbf{y} + \boldsymbol{\varepsilon}^{(i)} - H\mathbf{x}^{f,(i)}\right), \quad \boldsymbol{\varepsilon}^{(i)} \sim \mathcal{N}(\mathbf{0}, \mathbf{R})$$

| 長所 | 短所 |
|------|------|
| $\mathbf{B}$ の陽な格納が不要（$N \times n$ のアンサンブルのみ） | $N$ が小さいとサンプリング誤差が生じる |
| 非線形モデルにも適用可能 | 1サイクルごとに $N$ 回の前向き計算が必要 |
| 反復サイクルで画像が徐々に浮かび上がる演出ができる ✓ | 大次元系ではローカライゼーションが必要 |

**このデモへの適用可否**：最適。

- $N = 50$ メンバー × $n = 2400$ セル = 120,000 個の値（メモリは無視できる）
- 各 DA サイクルでスパース観測を同化 → 正解ドット絵に収束
- 収束プロセスをアニメーション化 →「おっさんが浮かび上がる」可視化が実現できる

---

## 推奨実装計画

```
Cycle  0 : 背景場 = ランダムノイズ（$\sigma_b = 20$）
Cycle  1 : 20 点の観測を同化 → うっすら輪郭が現れる
Cycle  5 : メガネ・顔の構造が見えてくる
Cycle 10 : おっさんの姿が認識できる
Cycle 20 : 正解場に収束
```

### 推奨パラメータ

| パラメータ | 値 | 備考 |
|-----------|-----|------|
| アンサンブルサイズ $N$ | 50 | $n=2400$ 次元に十分 |
| 観測数 $m$ | 20〜30 | 正解場からランダムサンプリング |
| 観測誤差 $\sigma_r$ | 5.0 | 温度レンジ $[0, 100]$ に対して相対的な値 |
| 背景誤差 $\sigma_b$ | 20.0 | 初期ランダムノイズの振幅 |
| ローカライゼーション半径 | 10 cells | 遠距離の擬似相関を除去 |
| DA サイクル数 | 20〜50 | 収束するまで反復 |

### RMSE の監視

各サイクル $k$ で以下を計算し収束を確認する：

$$\text{RMSE}_k = \frac{1}{\sqrt{n}} \left\|\bar{\mathbf{x}}^a_k - \mathbf{x}^{\text{true}}\right\|_2$$

---

## フォルダ構成

```
practice/A00_plan/
  plan.md                  ← プロジェクト概要
  da_method_selection.md   ← このファイル

practice/A01_ossan/        ← 正解温度場（おっさんのドット絵）
practice/A01_miffy/        ← 正解温度場（ミッフィーのドット絵）
practice/A02_enkf_ossan/   ← おっさんに対する EnKF データ同化
practice/A02_miffy_enkf_T/ ← ミッフィーに対する EnKF データ同化
practice/A03_analysis/     ← アニメーション・RMSE 収束プロット
```

---

## 参考文献

- Evensen (2003): *The Ensemble Kalman Filter: theoretical formulation and practical implementation*
- Hunt et al. (2007): LETKF — 大規模系のための高効率 EnKF
- 本プロジェクト：`blog/03_enkf/article.md`
