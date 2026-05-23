# 結果の見方

## 生成される図

実行すると次の画像が作られます。

```text
results/results_of_da.png
```

互換用に、同じ画像を次にもコピーしています。

```text
results_of_da.png
```

## 結果フォルダの構成

実行後の結果は `results/` フォルダに分けて保存します。

```text
results/
  truth/
    temperature_history.csv
    openfoam_time_dirs/
  model_only/
    temperature_history.csv
    openfoam_time_dirs/
  with_da/
    temperature_history.csv
    openfoam_forecast_before_oi.csv
    openfoam_time_dirs/
  summary_rmse.csv
  results_of_da.png
```

各フォルダの意味は次です。

| フォルダ | 意味 |
|---|---|
| `truth` | `G_TRUE = 100` でOpenFOAMを走らせた真値 |
| `model_only` | `G_MODEL = 50` でOpenFOAMを独立に走らせたDAなし結果 |
| `with_da` | `G_MODEL = 50` のOpenFOAM予測をOIで補正したDAあり結果 |

`openfoam_time_dirs/` には、OpenFOAMの時刻ディレクトリをコピーしています。たとえば `openfoam_time_dirs/100/T` は、時刻100秒のOpenFOAMフィールドです。

## temperature_history.csv とは

`temperature_history.csv` は、各ケースの温度時系列を表形式で保存したものです。

OpenFOAMの `T` ファイルは時刻フォルダごとに分かれているため、そのままだと全時刻の温度変化を見にくいです。そこで、Python側で読み取った温度を1つのCSVにまとめています。

列の意味は次です。

```text
time_s,N0,N1,N2,N3,N4,N5,N6,N7,N8,N9
```

| 列 | 意味 |
|---|---|
| `time_s` | 時刻 [s] |
| `N0` ... `N9` | 各節点の温度 [degC] |

たとえば `results/truth/temperature_history.csv` は、真値ケースの全節点温度時系列です。

`results/model_only/temperature_history.csv` は、データ同化なしのOpenFOAM結果です。

`results/with_da/temperature_history.csv` は、OIで補正したデータ同化ありの結果です。

`results/with_da/openfoam_forecast_before_oi.csv` は、DAありループの中で、OI補正をかける直前のOpenFOAM予測値です。

## 計算回数について

今回のデータ同化ありケース `with_da` では、OpenFOAMを1回だけ最後まで走らせているわけではありません。

各時刻で次を繰り返しています。

```text
OpenFOAMを1ステップ進める
温度場を読む
OIで補正する
補正後の温度場を書き戻す
次ステップへ進む
```

現在は `N_STEPS = 90` なので、`with_da` だけで `laplacianFoam` を90回起動しています。

今回のツイン実験では、比較のために `truth`, `model_only`, `with_da` の3ケースを回しているため、合計では概念的に270回の `laplacianFoam` 実行になります。

実測データを使う本番では `truth` は不要です。最低限は `with_da` だけでよく、DAなし比較もするなら `model_only` と `with_da` の2ケースになります。

図では、次の3つを比較しています。

| 名前 | 意味 |
|---|---|
| truth | `G_TRUE = 100` で OpenFOAM を走らせた真値 |
| OF only | `G_MODEL = 50` で OpenFOAM を独立に走らせたDAなし予測 |
| OI (DA) | `G_MODEL = 50` のOpenFOAM予測を、OIで補正した推定値 |

今回の真値は実験値ではなく、OpenFOAMで作った合成真値です。これはツイン実験です。

## 今回の結果

現在の設定では、`G_MODEL` を真値の50%にしているため、DAなしケースは真値から外れます。

一方で、DAありケースは Node 0 と Node 4 の観測に引き戻されます。

つまり今回の結果は、次のことを確認するためのものです。

```text
かなり間違った境界条件で計算すると、温度場は真値から大きく外れる。
しかし、一部の節点温度を観測として与えると、
OIにより観測していない節点も含めて真値に近づく傾向がある。
```

ただし、すべての点が必ず正解になるわけではありません。補正の広がり方は、センサ位置と背景誤差共分散 `B` に依存します。

相関長を調整した最新のOI版では、後半1/3の全節点平均RMSEは次のようになりました。

```text
DAなし : 1.186 degC
DAあり : 0.136 degC
改善率 : 90.3 %
```

## 図の各パネル

| パネル | 内容 |
|---|---|
| (1) Sensor nodes | センサ節点 Node 0, Node 4 の真値、DAなし、OIあり |
| (2) Hidden nodes | センサがない節点の真値、DAなし、OIあり |
| (3) RMSE | 後半1/3のRMSE |
| (4) Error | 未観測節点の誤差時系列 |

## OIで補正が広がる仕組み

Node 0 と Node 4 は観測で直接補正されます。

その他の節点は直接観測されませんが、背景誤差共分散 `B` によって間接的に補正されます。

今回の `B` は距離で作っています。

```text
近い節点  : 強く補正される
遠い節点  : 弱く補正される
```

この補正の広がり方は、主に次のパラメータで変わります。

```text
BACKGROUND_ERROR_STD
CORRELATION_LENGTH_M
OBS_NOISE_STD
```

## A_d がないことの意味

この結果は、OpenFOAMの離散化行列 `A_d` を使っていません。

つまり、Python側は、

```text
OpenFOAMがどのような行列で時間発展しているか
```

を知りません。

知っているのは、

```text
OpenFOAMの予測温度場
センサ観測値
センサ位置
背景誤差の空間的な広がり方
```

だけです。

このため、OpenFOAMをブラックボックスとして扱うデータ同化の第一段階として使いやすい構成です。

## 注意点

OIはシンプルですが、万能ではありません。

- `B` の相関長をどう決めるかで結果が変わる
- 境界勾配 `G` 自体は推定していない
- センサが少ない場所では補正が弱くなる

より本格的に進めるなら、次は次の候補があります。

- `G` を状態変数に含める
- センサ配置を変える
- アンサンブルカルマンフィルタに拡張する
