# sample

このフォルダには、熱データ同化の検証サンプルを段階別に置いています。

ルートの `README.md` はプロジェクト全体の説明です。各手法の目的、計算条件、実行方法はこの `sample/README.md` と各サンプル配下の `docs/` を参照してください。

## サンプル一覧

| サンプル | 目的 | 使う手法 | OpenFOAM | スライド |
|---|---|---|---|---|
| `001_kalman_thermal_1d` | Python だけでデータ同化の基本動作を確認する | Kalman filter | 不要 | [slides.html](001_kalman_thermal_1d/slides.html) |
| `002-0_laplacian_da_1d` | OpenFOAM の温度場予測をセンサ観測で補正する流れを確認する | optimal interpolation | 必要 | [slides.html](002-0_laplacian_da_1d/slides.html) |
| `002-1_laplacian_da_round_bar` | 丸棒の軸方向温度場を OpenFOAM / FrontISTR の2ソルバーで OI データ同化し比較する | optimal interpolation | 必要 | [slides.html](002-1_laplacian_da_round_bar/slides.html) |

## 001: Kalman filter による 1D 熱データ同化

場所:

```text
sample/001_kalman_thermal_1d/
```

### 目的

1 次元熱伝導棒を Python でモデル化し、Kalman filter によって温度場を推定します。

OpenFOAM は使いません。熱モデル、状態遷移行列、Kalman filter の予測・更新をすべて Python 側で扱います。

### 手法

Kalman filter は、モデル予測と観測値を誤差共分散に応じて混ぜる逐次推定手法です。

$$
\begin{aligned}
\text{予測:}\quad & \mathbf{x}^- = \mathbf{A}_d \mathbf{x} + \mathbf{B}_d \mathbf{u} \\
& \mathbf{P}^- = \mathbf{A}_d \mathbf{P}\mathbf{A}_d^\top + \mathbf{Q} \\
\text{更新:}\quad & \mathbf{K} = \mathbf{P}^- \mathbf{H}^\top \left(\mathbf{H}\mathbf{P}^- \mathbf{H}^\top + \mathbf{R}\right)^{-1} \\
& \mathbf{x}^+ = \mathbf{x}^- + \mathbf{K}(\mathbf{y} - \mathbf{H}\mathbf{x}^-) \\
& \mathbf{P}^+ = (\mathbf{I} - \mathbf{K}\mathbf{H})\mathbf{P}^-
\end{aligned}
$$

001 では、Python 側で `A_d` を作ります。そのため、手法の理解や小規模モデルの検証に向いています。

### 主なファイル

| ファイル | 役割 |
|---|---|
| `thermal_model.py` | 1 次元熱伝導モデルを作る |
| `kalman_filter.py` | Kalman filter の予測・更新を実装する |
| `main.py` | 基本デモを実行する |
| `main_da_demo.py` | 初期誤差・モデル誤差がある場合のデータ同化効果を見る |
| `main_disp_obs.py` | 変位観測を追加した場合を比較する |
| `slides.html` | 解説スライド |
| `docs/00_index.md` | 詳細ドキュメントの入口 |

### スライド

[001_kalman_thermal_1d/slides.html](001_kalman_thermal_1d/slides.html)

### 実行例

```bash
cd sample/001_kalman_thermal_1d
python main.py
python main_da_demo.py
python main_disp_obs.py
```

## 002: laplacianFoam + optimal interpolation

場所:

```text
sample/002-0_laplacian_da_1d/
```

### 目的

OpenFOAM の `laplacianFoam` を温度場の予測器として使い、Python 側でセンサ観測による補正を行います。

001 と違い、OpenFOAM の内部の状態遷移行列 `A_d` は作りません。OpenFOAM を 1 ステップ進めて出てきたセル温度を読み、optimal interpolation で補正し、補正後の温度場を次ステップの OpenFOAM 入力として書き戻します。

### 比較する 3 条件

| 名前 | 意味 |
|---|---|
| `truth` | 左端の加熱を正しい強さにした、正解として扱う計算結果 |
| `OpenFOAM only` | 左端の加熱を半分に間違えた、補正なしの計算結果 |
| `OpenFOAM DA (OI)` | 同じく加熱を半分に間違えた計算を、センサ観測で補正した結果 |

### 手法

002 では optimal interpolation を使います。

$$
\begin{aligned}
\mathbf{x}_a &= \mathbf{x}_f + \mathbf{K}(\mathbf{y} - \mathbf{H}\mathbf{x}_f) \\
\mathbf{K} &= \mathbf{B}\mathbf{H}^\top \left(\mathbf{H}\mathbf{B}\mathbf{H}^\top + \mathbf{R}\right)^{-1}
\end{aligned}
$$

ここで、$\mathbf{x}_f$ は OpenFOAM が 1 ステップ進めた予測温度場、$\mathbf{y}$ はセンサ観測、$\mathbf{H}$ はセンサ位置を取り出す行列、$\mathbf{B}$ は温度場の誤差が空間的にどう広がるかを仮定した行列、$\mathbf{R}$ はセンサノイズです。

### 主なファイル

| ファイル・フォルダ | 役割 |
|---|---|
| `of_interface.py` | Python から OpenFOAM の `T` ファイルを書き、`laplacianFoam` を実行し、結果を読む |
| `case_base/` | OpenFOAM ケースのひな形 |
| `oi/da_main.py` | optimal interpolation によるデータ同化の本体 |
| `oi/run.sh` | 実行スクリプト |
| `kf/README.md` | Kalman filter 版へ拡張する場合のメモ。今回の計算では未使用 |
| `slides.html` | 解説スライド |
| `docs/00_index.md` | 詳細ドキュメントの入口 |

### スライド

[002-0_laplacian_da_1d/slides.html](002-0_laplacian_da_1d/slides.html)

### 実行例

OpenFOAM 2512 が使える環境で実行します。

```bash
cd sample/002-0_laplacian_da_1d/oi
./run.sh
```

### 注意

`oi/results/` や `oi/paraview_cases/` には計算結果や ParaView 可視化用ファイルが生成されます。これらは重くなりやすいため、GitHub には基本的に push しません。

## 002-1: 丸棒 laplacianFoam + FrontISTR + optimal interpolation

場所:

```text
sample/002-1_laplacian_da_round_bar/       ← OpenFOAM 側
sample/002-1-FrontISTR_round_bar_da/       ← FrontISTR 側
```

長さ 0.3 m、直径 0.01 m の丸棒を 2 種類のソルバー (OpenFOAM FVM / FrontISTR FEM) で解き、同じ OI ループを接続して温度場を同化します。

- 同化センサ: N4, N16（2 点）
- 検証センサ: N28, N36（2 点）
- 熱入力: G=120 (ON 450s) → OFF (150s) のサイクル × 2
- 結果: OpenFOAM +93.8% / FrontISTR +94.0% の RMSE 改善

**スライド**: [▶ GitHub Pages で開く](https://kamakiri1225.github.io/thermal-da-demo/sample/002-1_laplacian_da_round_bar/slides.html)

```bash
# OpenFOAM 側
cd sample/002-1_laplacian_da_round_bar/oi
python da_main.py

# FrontISTR 側
cd sample/002-1-FrontISTR_round_bar_da/oi
python da_main.py
```

詳細は [002-1_laplacian_da_round_bar/README.md](002-1_laplacian_da_round_bar/README.md) を参照してください。

## 次に進めること

1. 002 の `truth` を実測温度データに置き換える
2. 同化用センサと検証用センサを分けて評価する
3. OI のセンサ配置・相関長・背景誤差標準偏差を検討する
4. `A_d` を作らない方針のまま、EnKF / local EnKF を検討する
