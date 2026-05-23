# 論文2：Sensor Placement Utilizing a Digital Twin for Thermal Error Compensation of Machine Tools

**著者**: Sebastian Lang, Mario Zavator, Stephan Schäfer, Stefan Blaser, Josef Mayr, Markus Bambach  
**掲載誌**: Journal of Manufacturing Systems  
**年**: 2024

---

## 表紙・アブストラクト

![p1](img/sensor_p01.png)

---

## 1. 背景と問題設定

### 熱誤差補償におけるセンサ配置の重要性

データ駆動型の熱誤差補償モデルでは、**どのセンサで温度を計測するか** が補償精度を大きく左右する。

- センサ数が多い → コスト増大・過学習リスク・設置困難
- センサ数が少ない → 情報不足・精度低下
- **配置が悪い → センサ数に関わらず精度が出ない**

### 従来の問題点

| アプローチ | 問題 |
|---|---|
| エンジニアの経験則による手動選択 | 属人的・最適解の保証なし |
| 試行錯誤による実機実験 | 時間・コストが大 |
| 静的な単一条件での最適化 | 多様な運転条件に汎化しない |

→ **実機実験前に、多様な運転条件に対して最適なセンサ配置を決定する手法** が必要。

![p2](img/sensor_p02.png)

---

## 2. 提案手法：デジタルツインベースのセンサ配置最適化

### システム全体のアーキテクチャ

![p3](img/sensor_p03.png)

デジタルツイン（FEM）を活用した **仮想的なセンサ選択** パイプライン：

```
① FEM シミュレーション
   └─ 多様な運転条件（スピンドル回転数・切削条件等）での
      温度場・熱変位を計算

② 仮想センサデータの生成
   └─ FEM 節点全体から「全候補センサ位置」のデータを抽出

③ センサ選択最適化（Group LASSO + k-means）
   └─ 補償モデルの精度を最大化するセンサ組み合わせを探索

④ 実機検証
   └─ 選ばれたセンサ位置に実センサを配置して精度確認
```

---

## 3. デジタルツインの構成

### FEM モデル

![p4](img/sensor_p04.png)

- Ansys / MORe ベースの熱・機械連成 FEM モデル
- 熱源（スピンドル・軸受・モータ）・熱伝導・対流境界条件を含む
- **多種の運転サイクルを仮想実験** として生成可能

### 代表的な温度場

![p5](img/sensor_p05.png)

様々な熱負荷条件での温度分布をシミュレーション。センサ候補は機械全面の格子点全体（数百〜数千点）。

---

## 4. センサ選択アルゴリズム

### 4-1. k-means クラスタリング

![p6](img/sensor_p06.png)

全候補センサ位置 $\mathcal{S} = \{s_1, s_2, \ldots, s_N\}$ における温度時系列データを行列 $\mathbf{T} \in \mathbb{R}^{N \times L}$（$L$: 時刻数）で表す。$k$-means により $k$ 個のクラスタに分割：

$$
\min_{\mathcal{C}_1,\ldots,\mathcal{C}_k} \sum_{j=1}^{k} \sum_{s_i \in \mathcal{C}_j} \left\| \mathbf{t}_i - \boldsymbol{\mu}_j \right\|_2^2
$$

| 記号 | 意味 |
|---|---|
| $\mathbf{t}_i \in \mathbb{R}^L$ | センサ候補 $s_i$ の温度時系列 |
| $\boldsymbol{\mu}_j$ | クラスタ $\mathcal{C}_j$ の重心 |
| $k$ | センサ数の上限（ハイパーパラメータ） |

各クラスタの重心 $\boldsymbol{\mu}_j$ に最も近い節点をセンサ候補として選択し、候補数を $N \to k$ に削減する。

---

### 4-2. Group LASSO による適応的センサ選択

補償モデルを ARX/回帰モデルとして構築し、**不必要なセンサグループの係数を正則化によってゼロに縮退** させる：

$$
\min_{\boldsymbol{\beta}} \left\| \mathbf{y} - \sum_{i=1}^{k} \mathbf{X}_i \boldsymbol{\beta}_i \right\|_2^2 + \lambda \sum_{i=1}^{k} \left\| \boldsymbol{\beta}_i \right\|_2
$$

| 記号 | 意味 |
|---|---|
| $\mathbf{y} \in \mathbb{R}^{T}$ | 目標値：熱変位の時系列（$T$: サンプル数） |
| $\mathbf{X}_i \in \mathbb{R}^{T \times d}$ | $i$ 番目センサグループの温度特徴量行列 |
| $\boldsymbol{\beta}_i \in \mathbb{R}^{d}$ | 対応する回帰係数ベクトル |
| $\lambda \geq 0$ | スパース正則化パラメータ |

$\lambda$ を大きくすると $\boldsymbol{\beta}_i = \mathbf{0}$ となるグループが増え、**選択センサ数が減少**する。Group LASSO は通常の LASSO と異なり、センサ $i$ に対応する係数グループを一括でゼロに縮退させる（グループスパース性）。

$$
\hat{\boldsymbol{\beta}}_i = \mathbf{0} \implies \text{センサ } i \text{ を除外}
$$

正則化パスを走査して、RMSE が最小となる $\lambda^*$ と対応するセンサセット $\mathcal{S}^* \subseteq \mathcal{S}$ を決定する：

$$
\mathcal{S}^* = \arg\min_{\mathcal{S}' \subseteq \mathcal{S}} \mathrm{RMSE}\!\left(\mathbf{y},\, \hat{\mathbf{y}}_{\mathcal{S}'}\right)
$$

---

## 5. 実験設定

### 工作機械とセンサ候補配置

![p7](img/sensor_p07.png)

- 対象機械: 5軸工作機械
- センサ候補位置: 機械全面（Front / Back / Top / Bottom / Left / Right の各面）
- 候補総数: 数十〜数百箇所

### FEM 由来のセンサ最適配置

![p8](img/sensor_p08.png)

FEM シミュレーション上で $\mathcal{S}^*$ として選ばれた配置（4〜8センサ）。熱源近傍・熱的に敏感な領域が優先的に選ばれる。

### 物理実験での検証

![p9](img/sensor_p09.png)

- $\mathcal{S}^*$ の位置に実センサを配置して実機運転サイクルで検証
- 比較対象: ランダム配置・エンジニア手動選択・全センサ使用

---

## 6. 結果

### 6-1. センサ数と補償精度のトレードオフ

![p10](img/sensor_p10.png)

| 配置方法 | センサ数 | RMSE |
|---|---|---|
| ランダム配置 | 4〜8 | 高い（ばらつき大） |
| エンジニア手動選択 | 8〜12 | 中程度 |
| **デジタルツイン最適化** ($\mathcal{S}^*$) | **4〜8** | **最小（安定）** |
| 全センサ使用 | 16+ | 過学習により劣化することも |

**少数センサでもデジタルツイン由来の配置が最善の精度** を達成。

### 6-2. 補償結果の時系列比較

![p11](img/sensor_p11.png)

$X$・$Y$・$Z$ 方向の熱変位について補償前後を比較：
- $\mathcal{S}^*$ センサによる補償後残差が最小
- 複数の運転サイクルにわたって安定

### 6-3. 詳細な補償性能評価

![p12](img/sensor_p12.png)

各センサ数・各配置戦略の RMSE を網羅的に比較。デジタルツインアプローチが一貫して優位。

---

## 7. さらなる詳細解析

![p13](img/sensor_p13.png)

![p14](img/sensor_p14.png)

- センサを増やしても精度が必ずしも向上しない（過学習）ことを確認
- 重要なのは **センサの位置** であり、数ではない
- $\mathcal{S}^*$ は熱的に重要な箇所を適切にカバーしている

---

## 8. 結論と参考文献

![p15](img/sensor_p15.png)

---

## 9. 論文の貢献まとめ

| 貢献 | 詳細 |
|---|---|
| 新規フレームワーク | FEM デジタルツイン上での仮想センサ選択最適化パイプライン |
| 実験工数削減 | 実機実験前にセンサ配置 $\mathcal{S}^*$ を決定できる |
| 精度向上 | 手動選択・ランダム配置を上回る補償精度 |
| センサ削減 | 4〜8センサで十分な精度（16+センサと同等以上） |
| 汎化性 | FEM で学んだ $\mathcal{S}^*$ が実機に転移する |

---

## 10. 残された課題・限界

- FEM モデルの精度が $\mathcal{S}^*$ の選択結果に直結（モデル誤差の影響が未検証）
- センサ配置（$\mathbf{C}$ の構造）とデータ同化（カルマンフィルタ）は独立に設計——統合最適化はなし
- 対象は固体熱伝導のみ（クーラントや冷却風の流体効果は未考慮）
- 異なる機械クラスへの転移可能性は未検討

---

## 11. 論文1（Kalman Filter）との対応関係

$$
\underbrace{\mathcal{S}^* \;\text{の決定}}_{\text{論文2}} \longrightarrow \underbrace{\mathbf{C} = \mathbf{C}(\mathcal{S}^*)}_{\text{観測行列の構築}} \longrightarrow \underbrace{\mathbf{K}_k = \mathbf{P}_k^-\mathbf{C}^\top(\mathbf{C}\mathbf{P}_k^-\mathbf{C}^\top + \mathbf{R})^{-1}}_{\text{論文1：カルマンゲインの計算}}
$$

センサ配置 $\mathcal{S}^*$ が観測行列 $\mathbf{C}$ の構造を決め、それがカルマンゲイン $\mathbf{K}_k$ の性質（推定誤差の大きさ）を直接支配する。両論文を組み合わせることで「**最適センサ × 効率的データ同化**」が実現する。
