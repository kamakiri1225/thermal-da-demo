# 002_laplacian_da_1d: laplacianFoam + OI データ同化テスト

## このサンプルの目的

このサンプルは、OpenFOAM の `laplacianFoam` を使って温度場を予測し、Python側でセンサ観測を使って温度場を補正するデータ同化テストです。

今回の同化手法は **Optimal Interpolation** です。略して **OI** と呼びます。

重要な点は、以前のように Python側で `A_d` という時間発展行列を作らないことです。つまり、

```text
OpenFOAM がどんな離散化行列で温度場を進めているか
```

を Python側で再現しません。

OpenFOAM はブラックボックスの予測器として扱い、Python側は

```text
OpenFOAM予測値 + センサ観測による補正
```

だけを行います。

## 何をしたか

1. `G_TRUE = 100 degC/m` で `laplacianFoam` を走らせ、合成真値を作る
2. `G_MODEL = 50 degC/m` で `laplacianFoam` を独立に走らせ、DAなしケースを作る
3. `G_MODEL = 50 degC/m` の OpenFOAM 予測を、Node 0 と Node 4 の観測で OI 補正する
4. 補正後の温度場を次ステップの OpenFOAM 初期温度場として書き戻す
5. DAなしとDAありのRMSEを比較する

## 現在の設定

| 項目 | 設定 |
|---|---|
| ソルバ | `laplacianFoam` |
| 同化手法 | Optimal Interpolation (OI) |
| セル数 | 10 |
| 棒の長さ | 0.3 m |
| セル幅 | 0.03 m |
| 時間刻み | 10 s |
| ステップ数 | 90 |
| 真値の境界勾配 | `G_TRUE = 100 degC/m` |
| 予測モデルの境界勾配 | `G_MODEL = 50 degC/m` |
| センサ位置 | Node 0, Node 4 |
| 観測ノイズ | `OBS_NOISE_STD = 0.05 degC` |
| 背景誤差標準偏差 | `BACKGROUND_ERROR_STD = 3.0 degC` |
| 空間相関長 | `CORRELATION_LENGTH_M = 0.09 m` |

`G_MODEL` を真値の50%にしているのは、データ同化を使わない場合と使った場合の差を見やすくするためです。

## ドキュメント

目次は [`00_index.md`](00_index.md) を参照。

| ファイル | 内容 |
|---|---|
| [`02_theory.md`](02_theory.md) | OIとは何か、データ同化手法の中での位置づけ、理論式、ゲインの意味、OpenFOAMとの接続 |
| [`03_program_guide.md`](03_program_guide.md) | プログラムの流れと、理論式がコードのどこに対応するか |
| [`04_results_guide.md`](04_results_guide.md) | 結果の見方、DAなしとDAありの違い |
| [`07_blog.md`](07_blog.md) | 今回の手法、利点・欠点、設定手順、結果、応用先、研究の問いをブログ形式で整理 |
| [`06_nextstep.md`](06_nextstep.md) | 丸棒実験で実測温度データを使うための具体的な手順と次のアクション |

## 実行方法

OpenFOAM 2512 が WSL 側にある前提です。

```bash
cd sample/002_laplacian_da_1d/oi
./run.sh
```

結果画像は次に保存されます。

```text
results/results_of_da.png
```

DAなし、DAあり、真値の温度時系列とOpenFOAM時刻フォルダは、`results/` 以下に分けて保存されます。

この実装は `A_d` を使わないため、OpenFOAM の離散化行列をPython側で推定・再現する必要がありません。これは、将来 `laplacianFoam` 以外のOpenFOAMソルバへ広げるときにも扱いやすい構成です。
