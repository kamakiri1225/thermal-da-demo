# プログラム解説: OIとコードの対応

## 全体の流れ

`da_main.py` は、OpenFOAMの `laplacianFoam` を1ステップずつ実行し、そのたびに OI で温度場を補正します。

```text
1. G_TRUE で laplacianFoam を走らせ、真値 T_true を作る
2. G_MODEL だけで laplacianFoam を独立に走らせ、DAなしケースを作る
3. G_MODEL で laplacianFoam を1ステップ進め、DAありケースの予測値 x_pred を得る
4. T_true のセンサ位置にノイズを加えて観測値 y を作る
5. OIで x_pred を補正し、x_da を得る
6. x_da を次ステップの OpenFOAM 初期温度場として書き戻す
```

## 今回のデータ同化手法

今回使っているのは **Optimal Interpolation (OI)** です。

OIの更新式は次です。

$$
\mathbf{x}^a
=
\mathbf{x}^f
+
\mathbf{K}
\left(
\mathbf{y}
-
\mathbf{H}\mathbf{x}^f
\right)
$$

ゲインは、

$$
\mathbf{K}
=
\mathbf{B}\mathbf{H}^T
\left(
\mathbf{H}\mathbf{B}\mathbf{H}^T
+
\mathbf{R}
\right)^{-1}
$$

です。

ここで、`B` は背景誤差共分散、`R` は観測誤差共分散です。

## A_d を使わない理由

現在のコードには、OpenFOAMの時間発展を近似する `A_d` はありません。

これは意図的です。

`A_d` を使う方法では、OpenFOAMが内部でどのような時間発展をしているかをPython側である程度再現する必要があります。しかし、OpenFOAMソルバが複雑になるほど、その行列を作るのは難しくなります。

OIでは、OpenFOAMの予測値をそのまま受け取り、固定した背景誤差共分散 `B` で補正します。

```text
OpenFOAM予測 x_pred : OpenFOAMが計算
補正ゲイン K        : Python側の B, H, R から計算
同化後 x_da         : x_pred + K(y - H x_pred)
```

## コード対応

### 観測行列 H、背景誤差 B、観測誤差 R

対応する関数は `build_oi_matrices()` です。

```python
H, B_mat, R_mat = build_oi_matrices()
```

`H` はセンサ位置だけを取り出す行列です。

```python
for i, node in enumerate(SENSOR_NODES):
    H[i, node] = 1.0
```

`B` はセル間距離から作ります。

```python
distance = np.abs(x[:, None] - x[None, :])
B = BACKGROUND_ERROR_STD**2 * np.exp(
    -(distance**2) / (2.0 * CORRELATION_LENGTH_M**2)
)
```

これは「近いセルほど同じ誤差を持ちやすい」という仮定です。

### OIゲイン

`run_da()` の冒頭で、OIゲインを計算します。

```python
gain = B @ H.T @ np.linalg.inv(H @ B @ H.T + R)
```

これは理論式の、

$$
\mathbf{K}
=
\mathbf{B}\mathbf{H}^T
\left(
\mathbf{H}\mathbf{B}\mathbf{H}^T
+
\mathbf{R}
\right)^{-1}
$$

に対応します。

### 観測値

ツイン実験なので、真値から観測値を作ります。

```python
y = H @ truth_T[k] + np.random.randn(n_obs) * OBS_NOISE_STD
```

実験データを使う場合は、この `y` を実測センサ値に置き換えます。

### OpenFOAM予測

温度場を進めるのは Python ではなく OpenFOAM です。

```python
ofi.write_field(k * DT, x_current, gradient=g_model_series[k])
ofi.run_step(k * DT, (k + 1) * DT)
x_pred = ofi.read_field((k + 1) * DT)
```

### OI更新

OpenFOAM予測と観測との差を使って補正します。

```python
innovation = y - H @ x_pred
x_da = x_pred + gain @ innovation
```

この `x_da` がデータ同化後の温度場です。

### 次ステップへの書き戻し

補正後の温度場を次の OpenFOAM 初期条件として使います。

```python
x_current = x_da
```

これにより、センサ情報が次ステップ以降のOpenFOAM計算に反映されます。

## 計算の回し方と計算時間

この実装では、DAありケース `with_da` の中で OpenFOAM を時間ステップごとに止めています。

```text
1. OpenFOAMを1ステップ進める
2. Pythonで温度場を読む
3. OIで補正する
4. 補正後の温度場をOpenFOAMへ書き戻す
5. 次のステップへ進む
```

つまり、OpenFOAMを一度だけ最後まで走らせてから後処理しているわけではありません。

今回の設定では `N_STEPS = 90` なので、DAありケースだけで `laplacianFoam` を90回起動します。

さらにツイン実験では、比較のために次の3ケースを実行しています。

| ケース | 目的 |
|---|---|
| `truth` | 正解データを作る |
| `model_only` | DAなし結果を作る |
| `with_da` | DAあり結果を作る |

そのため、今回のツイン実験全体では `laplacianFoam` を概念的に270回実行します。

実測データを使う場合は `truth` は不要です。最低限は `with_da` だけで実行できます。DAなし比較も欲しい場合は、`model_only` と `with_da` の2ケースを実行します。

計算時間が問題になる場合は、同化間隔を長くします。

```text
10秒ごとに同化 -> 90回
30秒ごとに同化 -> 約30回
60秒ごとに同化 -> 約15回
```

同化間隔を長くすると計算は軽くなりますが、観測を反映する頻度は下がります。

## この仕様で何がうれしいか

この仕様では、OpenFOAMソルバの方程式や離散化行列をPython側で作りません。

そのため、

- `laplacianFoam` 以外のソルバにも考え方を持ち込みやすい
- ソルバを改造しなくてもよい
- Python側の実装がシンプル
- 何の同化手法を使っているかが明確

という利点があります。

一方で、`B` の作り方が結果に効くため、相関長や背景誤差標準偏差の設定は検討が必要です。
