# 002-1: laplacianFoam + OI round bar data assimilation

`002_laplacian_da_1d` を丸棒実験向けにした派生サンプルです。

OpenFOAM の `laplacianFoam` を温度場予測器として使い、Python 側で optimal interpolation (OI) による補正を行います。丸棒は軸方向の温度分布を推定する目的に合わせ、長さ 0.3 m、直径 0.01 m の擬似 1 次元モデルとして扱います。

## OI 更新式

$$
\mathbf{x}_a = \mathbf{x}_f + \mathbf{K}(\mathbf{y} - \mathbf{H}\mathbf{x}_f)
$$

$$
\mathbf{K} = \mathbf{B}\mathbf{H}^\top\left(\mathbf{H}\mathbf{B}\mathbf{H}^\top + \mathbf{R}\right)^{-1}
$$

ここで、\(\mathbf{x}_f\) はソルバー予測、\(\mathbf{x}_a\) は OI 補正後の状態、\(\mathbf{y}\) は観測、\(\mathbf{H}\) は観測行列、\(\mathbf{B}\) は背景誤差共分散、\(\mathbf{R}\) は観測誤差共分散です。

## モデル条件

| 項目 | 値 |
|---|---:|
| 丸棒長さ | 0.30 m |
| 丸棒直径 | 0.01 m |
| 軸方向セル数 | 40 |
| 断面セル数 | 12 |
| 総セル数 | 480 |
| 時間刻み | 5 s |
| 計算時間 | 900 s |
| 熱拡散率 | 6.6e-5 m2/s |
| 真値の左端温度勾配 | 120 degC/m |
| 予測モデルの左端温度勾配 | 80 degC/m |

## センサ配置

| 用途 | セル | 位置 |
|---|---:|---:|
| 同化 | N4 | x = 0.034 m |
| 同化 | N16 | x = 0.124 m |
| 検証 | N28 | x = 0.214 m |
| 検証 | N36 | x = 0.274 m |

同化センサだけを OI の観測 `y` に使います。検証センサは、同化に使っていない位置でも推定が改善しているかを見るために残しています。

## 実行方法

OpenFOAM 2512 が使える環境で実行します。

```bash
cd sample/002-1_laplacian_da_round_bar/oi
chmod +x run.sh
./run.sh
```

結果は `oi/results/` に保存されます。

| 出力 | 内容 |
|---|---|
| `results/results_of_da.png` | 真値、DAなし、OIありの比較図 |
| `results/summary_rmse.csv` | 節点ごとのRMSE |
| `results/truth/temperature_history.csv` | 真値温度履歴 |
| `results/model_only/temperature_history.csv` | DAなし温度履歴 |
| `results/with_da/temperature_history.csv` | OI補正後の温度履歴 |

## 実測データへ置き換える場合

現在はツイン実験です。`oi/da_main.py` の `run_da()` 内で、

```python
y = H @ truth_T[k] + np.random.randn(n_obs) * OBS_NOISE_STD
```

として合成観測を作っています。実測丸棒データを使う場合は、この `y` を CSV から読んだ同化センサ温度に置き換えます。
