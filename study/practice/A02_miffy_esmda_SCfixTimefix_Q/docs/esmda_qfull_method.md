# 完全逆問題版: 状態変数を熱源場 Q にした ESMDA

`python/esmda_qfull.py` の手法解説です。
温度場を同化する `enkf_heatid.py` と対になるケースで、
HANDOFF の「次の技術課題」（状態変数を $\mathbf{Q}_{\rm full}$ に変更し、
$\mathbf{L}^{-1}$ の順問題を組み込んだ完全逆問題へ進める）を実装したものです。

---

## 1. 何が問題だったか

`enkf_heatid.py` は温度場 $\mathbf{x}$ を同化し、熱源は最後に

$$
\mathbf{q}_{\rm derived}
= -\alpha \mathbf{L}\mathbf{x}_{\rm mean}
$$

で逆算していました。ラプラシアン $\mathbf{L}$ は差分オペレータなので、
同化後の温度場に残るセルスケールのノイズを強く増幅します。
その結果 `fig03_qest_final.png` の熱源図はノイズが支配的で、
熱源分布として直接読める図にはなりませんでした。

## 2. 完全逆問題の定式化

状態変数を熱源場 $\mathbf{q}$ そのものに変更します。

| 項目 | 温度場版 (`enkf_heatid.py`) | 完全逆問題版 (`esmda_qfull.py`) |
|---|---|---|
| 状態変数 | $\mathbf{x}$（温度、900セル） | $\mathbf{q}$（熱源、900セル） |
| 観測演算子 | $\mathbf{H} = \mathbf{I}_{\rm obs}$（恒等の抜き出し） | $\mathcal{H}(\mathbf{q}) = [\mathbf{L}^{-1}(-\mathbf{q}/\alpha)]_{\rm obs}$（順問題） |
| 事前分布 | 20 °C 一様 + 相関ノイズ | 0 中心 + 相関ノイズ |
| 熱源の出し方 | 事後に $\mathbf{q} = -\alpha\mathbf{L}\mathbf{x}$ で逆算 | $\mathbf{q}$ が直接の推定対象 |

順問題は定常熱伝導（ポアソン方程式）です。

$$
\mathbf{L}\mathbf{x}
= -\frac{\mathbf{q}}{\alpha}
\quad\Longrightarrow\quad
\mathbf{x}(\mathbf{q})
= \mathbf{L}^{-1}\left(-\frac{\mathbf{q}}{\alpha}\right)
$$

$\mathbf{L}$ は一度だけ LU 分解（`scipy.sparse.linalg.splu`）しておき、
毎反復はアンサンブル行列 (900×300) をまとめて後退代入するだけなので
1 ランは数秒で終わります。

ESMDA の更新式は温度場版と同一で、アンサンブル摂動の取り方だけが
異なります。温度場版は $\mathbf{HX}' = \mathbf{X}'_{\rm obs}$ でしたが、完全逆問題版は
順問題を通した温度摂動 $\mathbf{HT}' = \mathbf{T}'_{\rm obs}$ を使います。

$$
\mathbf{S}
=
\frac{\mathbf{HT}'{\mathbf{HT}'}^{\!\top}}{N-1}
+ N_{\rm iter}\sigma_r^2\mathbf{I}
$$

$$
\mathbf{K}
=
\left(
\frac{\mathbf{Q}'{\mathbf{HT}'}^{\!\top}}{N-1}
\right)
\mathbf{S}^{-1}
$$

$$
\mathbf{Q}
\leftarrow
\mathbf{Q}
+
(\mathbf{K}\circ\boldsymbol{\rho})
\left(
\mathbf{y}
+
\boldsymbol{\varepsilon}
-
\mathbf{T}_{\rm obs}
\right),
\qquad
\boldsymbol{\varepsilon}
\sim
\mathcal{N}
\left(
\mathbf{0},
N_{\rm iter}\sigma_r^2\mathbf{I}
\right)
$$

## 3. これはなぜデータ同化なのか

「熱源 $\mathbf{q}$ の候補をたくさん作り、それぞれが予測するセンサ温度と
実測センサ温度のズレを使って $\mathbf{q}$ を更新する」と書くと、単なる
試行錯誤に見えます。しかし実際には、データ同化の標準形

$$
\text{解析値}
=
\text{背景値}
+
\text{ゲイン}
\times
\text{観測残差}
$$

を、状態変数を熱源 $\mathbf{q}$ にして実行しています。

データ同化に必要な要素は次の 4 つです。

| 要素 | 一般のデータ同化 | 今回の完全逆問題版 |
|---|---|---|
| 状態 | 推定したい変数 $\mathbf{x}$ | 熱源場 $\mathbf{q}$ |
| 背景 | 観測前の推定値・候補 | 300 本の熱源候補 $\mathbf{Q}$ |
| 観測 | 実測データ $\mathbf{y}$ | 固定センサ温度 $\mathbf{y}_{\rm obs}$ |
| 観測演算子 | 状態から観測値を予測する $\mathbf{H}$ | $\mathcal{H}(\mathbf{q})=\mathbf{H}\mathbf{L}^{-1}(-\mathbf{q}/\alpha)$ |

ポイントは、観測が熱源 $\mathbf{q}$ を直接測っていないことです。
観測しているのは温度なので、候補 $\mathbf{q}$ をいったん順問題に通して
「その熱源ならセンサ温度は何度になるか」を計算します。

順問題は

$$
\mathbf{L}\mathbf{x}
=
-\frac{\mathbf{q}}{\alpha}
$$

なので、温度場は

$$
\mathbf{x}(\mathbf{q})
=
\mathbf{L}^{-1}
\left(
-\frac{\mathbf{q}}{\alpha}
\right)
$$

です。センサ位置だけを取り出す行列を $\mathbf{H}$ とすると、熱源 $\mathbf{q}$ が予測する
センサ温度は

$$
\mathbf{y}_{\rm pred}(\mathbf{q})
=
\mathbf{H}\mathbf{x}(\mathbf{q})
=
\mathbf{H}\mathbf{L}^{-1}
\left(
-\frac{\mathbf{q}}{\alpha}
\right)
$$

です。この写像を観測演算子として

$$
\mathcal{H}(\mathbf{q})
=
\mathbf{H}\mathbf{L}^{-1}
\left(
-\frac{\mathbf{q}}{\alpha}
\right)
$$

と置きます。すると、今回の逆問題は

$$
\mathbf{y}_{\rm obs}
\approx
\mathcal{H}(\mathbf{q})
$$

を満たす $\mathbf{q}$ を探す問題になります。

ESMDA では $\mathbf{q}$ を 1 本だけ持たず、$N$ 本の候補を持ちます。

$$
\mathbf{Q}
=
\begin{bmatrix}
\mathbf{q}_1 & \mathbf{q}_2 & \cdots & \mathbf{q}_N
\end{bmatrix}
$$

各候補を観測演算子に通すと、

$$
\mathbf{Y}
=
\mathcal{H}(\mathbf{Q})
=
\begin{bmatrix}
\mathcal{H}(\mathbf{q}_1) &
\mathcal{H}(\mathbf{q}_2) &
\cdots &
\mathcal{H}(\mathbf{q}_N)
\end{bmatrix}
$$

です。ここで $\mathbf{Q}$ は「熱源候補の集団」、$\mathbf{Y}$ は「各候補が予測した
センサ温度の集団」です。

まず平均を引いた偏差を作ります。

$$
\mathbf{Q}'
=
\mathbf{Q}
-
\bar{\mathbf{q}}\mathbf{1}^{\top},
\qquad
\mathbf{Y}'
=
\mathbf{Y}
-
\bar{\mathbf{y}}\mathbf{1}^{\top}
$$

この 2 つから、アンサンブル標本共分散を作ります。

$$
\mathbf{C}^{qy}
=
\frac{1}{N-1}
\mathbf{Q}'{\mathbf{Y}'}^{\!\top}
$$

$$
\mathbf{C}^{yy}
=
\frac{1}{N-1}
\mathbf{Y}'{\mathbf{Y}'}^{\!\top}
$$

$\mathbf{C}^{qy}$ は「どのセルの熱源 $\mathbf{q}$ が、どのセンサ温度 $\mathbf{y}$ と一緒に動くか」
を表します。$\mathbf{C}^{yy}$ は「センサ温度同士がどう一緒に動くか」を表します。

カルマンゲインは

$$
\mathbf{K}
=
\mathbf{C}^{qy}
\left(
\mathbf{C}^{yy}
+
\alpha_k\mathbf{R}
\right)^{-1}
$$

です。観測誤差共分散 $\mathbf{R}$ は、センサノイズの標準偏差が
3 °C、つまり $\sigma_r=3$ なら

$$
\mathbf{R}
=
\sigma_r^2\mathbf{I}
=
9\mathbf{I}
$$

です。ESMDA では同じ観測を $N_{\rm iter}$ 回に分けて使うため、
各反復では観測誤差を $\alpha_k$ 倍に膨らませます。今回の一様分割では

$$
\alpha_k
=
N_{\rm iter}
=
30
$$

です。

各メンバーの更新式は

$$
\mathbf{q}_i^{a}
=
\mathbf{q}_i^{b}
+
\mathbf{K}
\left(
\mathbf{y}_{\rm obs}
+
\boldsymbol{\varepsilon}_i
-
\mathcal{H}(\mathbf{q}_i^{b})
\right)
$$

です。

ここで、

- `b` は before（更新前）
- `a` は after（更新後）
- $\mathbf{y}_{\rm obs}-\mathcal{H}(\mathbf{q}_i^b)$ は観測残差、または innovation
- $\boldsymbol{\varepsilon}_i$ は摂動観測ノイズで、$\boldsymbol{\varepsilon}_i \sim \mathcal{N}(\mathbf{0},\alpha_k\mathbf{R})$

です。摂動観測を入れるのは、更新後のアンサンブルの広がりを
正しい事後分布の広がりに保つためです。

さらに実装では局所化を入れています。セル i とセンサ j の距離を
`d_ij` とすると、

$$
\rho_{ij}
=
\exp
\left(
-\frac{d_{ij}^2}{2R_{\rm loc}^2}
\right)
$$

を作り、ゲインに要素ごとに掛けます。

$$
\mathbf{K}_{\rm loc}
=
\mathbf{K}
\circ
\boldsymbol{\rho}
$$

したがって実装上の更新は

$$
\mathbf{q}_i^{a}
=
\mathbf{q}_i^{b}
+
\mathbf{K}_{\rm loc}
\left(
\mathbf{y}_{\rm obs}
+
\boldsymbol{\varepsilon}_i
-
\mathcal{H}(\mathbf{q}_i^{b})
\right)
$$

です。

以上をまとめると、今回の計算は

1. 観測前の熱源候補 $\mathbf{Q}$ を作る
2. 物理モデル $\mathbf{x}(\mathbf{q})=\mathbf{L}^{-1}(-\mathbf{q}/\alpha)$ で各候補の温度場を計算する
3. センサ位置だけを抜き出して観測予測 $\mathcal{H}(\mathbf{q})$ を作る
4. 実測センサ温度との差 $\mathbf{y}_{\rm obs}-\mathcal{H}(\mathbf{q})$ を計算する
5. アンサンブル共分散から作ったゲイン $\mathbf{K}$ で $\mathbf{q}$ を更新する
6. これを 30 回繰り返す

という流れです。これは、観測と物理モデルを組み合わせて状態を更新する
データ同化そのものです。違いは、通常の天気予報のように時刻を進めながら
逐次同化するのではなく、静的な逆問題に対して同じ観測を複数回に分けて
使うスムーザー型のデータ同化、つまり ESMDA である点です。

## 4. 実行

```bash
cd python
python3 esmda_qfull.py
```

必要パッケージは `numpy`、`scipy`、`matplotlib`、`Pillow`。
BLAS スレッドはスクリプト内で 4 に制限しています
（このサイズの行列では全コアスレッド化はかえって大幅に遅くなるため）。

## 5. 設定

| パラメータ | 値 | 備考 |
|---|---|---|
| グリッド | 30×30 = 900 | 温度場版と共通 |
| アンサンブル数 N | 300 | 〃 |
| 反復数 N_ITER | 30 | 〃 |
| 観測ノイズ $\sigma_r$ | 3 °C | 〃 |
| 事前分布の広がり `SIGMA_QB` | 40 | $\mathbf{q}_{\rm true}$ の std ≈ 40 に一致 |
| 相関スケール CORR_L | 1.5 cells | 〃 |
| 局所化半径 R_LOC | **10 cells** | 温度場版（5）より広く取る必要がある |

チューニングの結果:

- `R_LOC = 3` では RMSE が 3 桁発散する。順問題 $\mathbf{L}^{-1}$ は非局所なので、
  実在する長距離の $\mathbf{q}$–$\mathbf{T}$ 共分散を強く切り落とすと更新が壊れる。
- `R_LOC = 5` でも最終値は収束するが、反復 4〜7 で温度 RMSE が
  最大 5000 °C 超まで暴走する過渡が出る。局所化された更新が
  センサ空白域の温度を崩し、後の反復で修正される挙動。
- `R_LOC = 10` では過渡が完全に消え、最終精度も同等以上。
  少数センサ時の悪化（下記 5.5）も解消するため 10 を採用。
- `SIGMA_QB` を 80, 120 に増やしても、$N=600$ にしても、$N_{\rm iter}=60$ にしても
  q RMSE は改善しない（≈ 39 で飽和）。

## 6. 結果と解釈

### 6.1 温度場は同等の精度で復元できる

$m=300$ で $\mathbf{x}(\mathbf{q}_{\rm mean})$ の温度 RMSE は **20.06 °C**。
温度場を直接同化した場合（19.82 °C）とほぼ同じです。
収束も速く、実質 5 反復で最終値に達します。

### 6.2 熱源のセルスケールのエッジは原理的に復元できない

$\mathbf{q}_{\rm true}=-\alpha\mathbf{L}\mathbf{x}_{\rm ss}$ はミッフィーの「輪郭線」状のエッジ場で、
std は 40.0。推定 $\mathbf{q}$ の RMSE は $m=400$ でも 38.8 までしか下がりません。
順問題 $\mathbf{L}^{-1}$ が強い平滑化作用を持つため、点温度観測には
熱源の高周波成分の情報がほとんど残っていない
（逆問題の不適切性 ill-posedness）ためです。

### 6.3 復元できるのは「平滑化されたスケール」まで

$\sigma=1$ セルでガウス平滑化した $\mathbf{q}_{\rm true}$ と比較すると、
推定 $\mathbf{q}$ との相関はセンサ数とともに単調に向上します。

| センサ数 $m$ | corr($\mathbf{q}_{\rm est}$, $\mathbf{q}_{\rm true,smoothed}$) |
|---:|---:|
| 30 | 0.33 |
| 100 | 0.59 |
| 200 | 0.65 |
| 300 | **0.76** |

$m=300$ の推定 $\mathbf{q}$ には目・鼻・輪郭の構造が視認できます
（`fig08_qfull_vs_derived.png`）。

### 6.4 旧方式との違いは「情報量」ではなく「使いやすさ」

意外なことに、温度場版の $\mathbf{q}_{\rm derived}$ も平滑化真値との相関は
大差ありません（$m=300$ で 0.74、Q-state は 0.76）。線形問題なので、
どちらの定式化でも観測から取り出せる低周波情報はほぼ同じ、
というのが理論通りの結果です。

違いは推定場の性格に出ます。

- $\mathbf{q}_{\rm derived}$: 信号がセルスケールの増幅ノイズ（±70 程度）に埋もれ、
  そのままでは熱源図として読めない
- Q-state: 事前分布が $\mathbf{q}$ 自体を正則化するため、滑らかで
  そのまま表示・利用できる熱源場が得られる

### 6.5 局所化半径が狭いと同化が温度場を壊す（R_LOC=5 の教訓）

開発中、温度場版と同じ R_LOC=5 で走らせたところ、

- 反復途中で温度 RMSE が最大 5000 °C 超まで暴走する過渡が出る
- m=30〜50 では最終値でも 47〜62 °C となり、
  同化しない場合（≈38 °C）より悪化する

という現象が出ました。局所化はセンサ近傍の $\mathbf{q}$ だけを更新しますが、
更新された熱源は $\mathbf{L}^{-1}$ を通じて遠方の温度まで変えます。センサが疎だと
この遠方への副作用を観測で拘束できず、センサ空白域の温度場が
崩れるためです。R_LOC=10 に広げると両方とも解消しました。
順問題が非局所な逆問題に局所化を使うときの代表的な落とし穴です。

## 7. 出力ファイル

| ファイル | 内容 |
|---|---|
| `img/fig05_qfull_rmse.png` | $\mathbf{q}$ / 温度 RMSE の収束（旧方式の水準線付き） |
| `img/fig06_qfull_qest_final.png` | 最終推定熱源場（m 別） |
| `img/fig07_qfull_xrecon_final.png` | 推定熱源からの順解析温度場（m 別） |
| `img/fig08_qfull_vs_derived.png` | 真値・平滑化真値・Q-state・旧方式の比較 |
| `img/fig09_qfull_sensor_count_sweep.png` | センサ数スイープ（q と温度） |
| `img/anim_esmda_qfull.gif` | 同化過程アニメーション |
| `docs/sensor_count_sweep_qfull.csv` | スイープ数値データ |

## 8. 関連

- [`enkf_heatid_method.md`](enkf_heatid_method.md) — 温度場版の手法
- [`esmda_explanation.md`](esmda_explanation.md) — ESMDA 一般論
- [`sensor_count_study.md`](sensor_count_study.md) — 温度場版のセンサ数検討
