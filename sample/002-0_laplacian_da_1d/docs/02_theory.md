# 理論説明: laplacianFoam と Optimal Interpolation

## 1. 今回使っているデータ同化手法

今回の手法は **Optimal Interpolation (OI)** です。

日本語では「最適内挿」と呼ばれます。考え方は単純で、

```text
OpenFOAMの予測値を、センサ観測との差で補正する
```

というものです。

以前のカルマンフィルタ型の実装では、Python側に `A_d` という時間発展行列を作っていました。しかし、それには「OpenFOAM がどのような方程式・離散化で温度場を進めているかを、Python側でもある程度再現できる」という前提がありました。

今回のOIでは、その前提を置きません。

```text
OpenFOAMの時間発展 : OpenFOAMに任せる
観測による補正     : PythonのOIで行う
```

このため、Python側に `A_d` はありません。

OIは、次のような場面で使いやすい手法です。

- シミュレーション結果はある
- センサ観測もある
- ただし、シミュレーションモデルの時間発展行列は作りたくない
- どの場所の誤差がどの場所へ似た形で広がるかは、ある程度仮定できる

今回で言えば、OpenFOAMの内部離散化は追いかけずに、

```text
近いセル同士は似た誤差を持ちやすい
```

という空間相関だけを仮定しています。

## 2. OpenFOAM が解く問題

`laplacianFoam` は熱拡散型の方程式を解きます。

$$
\frac{\partial T}{\partial t}
= \nabla \cdot (\alpha \nabla T)
$$

このサンプルでは1次元棒として扱っています。

左端には温度勾配を与え、右端は20 degCに固定します。

```text
左端: fixedGradient
右端: fixedValue 20 degC
```

真値と予測モデルでは、左端境界勾配だけを変えています。

```text
真値   : G_TRUE  = 100 degC/m
モデル : G_MODEL =  50 degC/m
```

`G_MODEL` を真値の50%にしているのは、DAなしとDAありの差を見やすくするためです。

## 3. 状態ベクトル

10セルの温度を状態ベクトルとします。

$$
\mathbf{x}
=
\begin{bmatrix}
T_0 & T_1 & \cdots & T_9
\end{bmatrix}^T
$$

OpenFOAMによる1ステップ予測を次のように書きます。

$$
\mathbf{x}_k^f
=
\mathcal{M}_{OF}
\left(
\mathbf{x}_{k-1}^a,\,
G_{\mathrm{model}}
\right)
$$

ここで、

| 記号 | 意味 |
|---|---|
| `f` | forecast、予測値 |
| `a` | analysis、同化後の解析値 |
| `M_OF` | OpenFOAMの時間発展 |

`M_OF` の中身を Python側で行列化しないことがポイントです。

## 4. 観測式

センサは Node 0 と Node 4 にあります。

```text
SENSOR_NODES = [0, 4]
```

観測式は、

$$
\mathbf{y}_k
=
\mathbf{H}\mathbf{x}^{true}_k
+
\mathbf{v}_k
$$

です。

観測行列 `H` は、状態ベクトルからセンサ位置だけを取り出す行列です。

$$
\mathbf{H}
=
\begin{bmatrix}
1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 & 0
\end{bmatrix}
$$

観測ノイズ共分散は、

$$
\mathbf{R}
=
\sigma_{obs}^2 \mathbf{I}
$$

です。

現在は、

```text
OBS_NOISE_STD = 0.05 degC
```

です。

## 5. OIの基本的な考え方

OIでは、まず OpenFOAM の予測値を「背景値」と呼びます。

$$
\mathbf{x}^f
$$

ここで `f` は forecast の意味です。

一方、センサ観測は、

$$
\mathbf{y}
$$

です。

OpenFOAM予測が完全に正しいなら、そのセンサ位置の値 `H x_f` は観測 `y` と一致するはずです。しかし実際には差があります。

$$
\mathbf{d}
=
\mathbf{y}
-
\mathbf{H}\mathbf{x}^f
$$

この `d` を **イノベーション** と呼びます。

日本語で言えば、

```text
観測値 - 予測値
```

です。

OIは、この差 `d` を使って、予測値全体を補正します。

$$
\mathbf{x}^a
=
\mathbf{x}^f
+
\mathbf{K}\mathbf{d}
$$

ここで `a` は analysis の意味で、同化後の解析値です。

重要なのは、センサがある節点だけを補正するのではなく、行列 `K` によってセンサがない節点にも補正を広げることです。

## 6. OIの更新式

OIでは、OpenFOAMの予測値 `x_f` を、観測との差で補正します。

$$
\mathbf{x}_k^a
=
\mathbf{x}_k^f
+
\mathbf{K}
\left(
\mathbf{y}_k
-
\mathbf{H}\mathbf{x}_k^f
\right)
$$

括弧の中はイノベーションです。

$$
\mathbf{d}_k
=
\mathbf{y}_k
-
\mathbf{H}\mathbf{x}_k^f
$$

これは「センサ値と予測値の差」です。

OIゲインは次で計算します。

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

ここで `B` は背景誤差共分散です。

この式は、次の2つの誤差情報のバランスで補正量を決めています。

| 行列 | 意味 |
|---|---|
| `B` | OpenFOAM予測がどれくらい、どの場所で、どのように間違いやすいか |
| `R` | センサ観測がどれくらい不確かか |

`B` が大きく `R` が小さい場合、モデルより観測を強く信じます。

`B` が小さく `R` が大きい場合、観測よりモデルを強く信じます。

## 7. OIゲインの意味

OIゲイン `K` は、観測との差をどの節点へどれだけ配るかを決める行列です。

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

この `K` は、広い意味では **カルマンゲインと同じ形のゲイン行列** です。

ただし、今回のOIでは次の点が通常のカルマンフィルタと異なります。

```text
カルマンフィルタ:
  時間発展させた予測誤差共分散 P^- を使う
  K = P^- H^T (H P^- H^T + R)^-1

OI:
  固定した背景誤差共分散 B を使う
  K = B H^T (H B H^T + R)^-1
```

したがって、OIの `K` は「カルマンゲインそのもの」と呼んでも式の構造としては近いですが、より正確には **OIゲイン**、または **固定した背景誤差共分散 `B` を使ったカルマンゲイン型の重み** と呼ぶのが分かりやすいです。

この式は一見難しく見えますが、1つのセンサだけで考えると直感的です。

ある1点だけを観測している場合、補正はおおよそ次の形になります。

$$
x_i^a
=
x_i^f
+
\frac{B_{i,s}}{B_{s,s} + R}
\left(
y_s - x_s^f
\right)
$$

ここで `s` はセンサ節点です。

| 記号 | 意味 |
|---|---|
| `x_i` | 補正したい節点 i の温度 |
| `x_s` | センサ節点 s の温度 |
| `B_{i,s}` | 節点 i とセンサ節点 s の誤差の似ている度合い |
| `B_{s,s}` | センサ節点自身の背景誤差分散 |
| `R` | 観測誤差分散 |

この式は、次の3段階に分けると読みやすくなります。

### 1. センサ位置で予測のズレを求める

$$
y_s - x_s^f
$$

これは、センサ位置での

```text
実測値 - OpenFOAM予測値
```

です。

たとえば、

```text
y_s   = 30 degC
x_s^f = 26 degC
```

なら、

```text
y_s - x_s^f = 4 degC
```

です。

### 2. そのズレを節点 i にどれくらい伝えるかを決める

$$
\frac{B_{i,s}}{B_{s,s} + R}
$$

これは重みです。

```text
1に近い -> センサのズレを強く反映する
0に近い -> センサのズレをあまり反映しない
```

`B_{i,s}` が大きいほど、節点 `i` とセンサ `s` の誤差は似ていると考えます。

`R` が大きいほど、センサ値は不確かなので、重みは小さくなります。

この重みは、現在のプログラムでは次の順番で計算しています。

1. 節点間の距離を求める
2. 距離から背景誤差共分散 `B` を作る
3. 観測ノイズから `R` を作る
4. `B` と `R` から重み、つまりOIゲイン `K` を作る

コードではここに対応します。

```python
x = np.arange(N_NODES) * DX
distance = np.abs(x[:, None] - x[None, :])
B = BACKGROUND_ERROR_STD**2 * np.exp(
    -(distance**2) / (2.0 * CORRELATION_LENGTH_M**2)
)
R = np.eye(n_obs) * OBS_NOISE_STD**2
gain = B @ H.T @ np.linalg.inv(H @ B @ H.T + R)
```

現在の設定では、

```text
BACKGROUND_ERROR_STD = 3.0 degC
CORRELATION_LENGTH_M = 0.06 m
OBS_NOISE_STD = 0.05 degC
DX = 0.03 m
```

です。

たとえば1つのセンサだけで考えると、距離が近い節点ほど重みが大きくなります。

| センサからの距離 | 重みの目安 |
|---:|---:|
| 0.00 m | 約 1.00 |
| 0.03 m | 約 0.88 |
| 0.06 m | 約 0.61 |
| 0.09 m | 約 0.32 |
| 0.12 m | 約 0.14 |

これは、相関長 `0.06 m` に対して、センサから遠ざかるほど「同じ誤差を持っている可能性が低い」と仮定しているためです。

複数センサの場合は、各センサからの補正が重なります。そのため、単純に各センサの重みを足すのではなく、

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

で、センサ同士の重複や観測ノイズも含めて一度に計算します。

### 3. 補正前の温度に補正量を足す

$$
x_i^a
=
x_i^f
+
\text{重み}
\times
\text{センサ位置でのズレ}
$$

たとえば、

```text
x_i^f = 24 degC
センサ位置でのズレ = 4 degC
重み = 0.8
```

なら、

```text
x_i^a = 24 + 0.8 * 4 = 27.2 degC
```

です。

重みが `0.1` なら、

```text
x_i^a = 24 + 0.1 * 4 = 24.4 degC
```

です。

つまり、この式は

```text
センサで分かったズレを、節点ごとの重みに応じて配る
```

という意味です。

この式から分かることは次です。

```text
B_{i,s} が大きい  -> センサとの差が節点 i に強く伝わる
B_{i,s} が小さい  -> センサとの差は節点 i にあまり伝わらない
R が大きい        -> 観測をあまり信用しない
R が小さい        -> 観測を強く信用する
```

つまりOIは、

```text
観測との差を、背景誤差共分散 B に従って空間的に配る方法
```

と見ることができます。

### BからKを作る直感

`B` は、そのまま補正量ではありません。

`B` はあくまで、

```text
どの節点の誤差が、どのセンサの誤差と似ていそうか
```

を表す情報です。

たとえば、センサが1点だけ、Node 4 にあるとします。

このとき、`B` の中で重要なのは、各節点とNode 4の関係を表す列です。

```text
B[:, 4]
```

これは、

```text
Node 4で分かったズレを、
Node 0, Node 1, ..., Node 9 にどれくらい配るか
```

の候補になります。

ただし、そのまま `B[:, 4]` を使うと、観測誤差の大きさが考慮されません。センサがノイズだらけなら、観測値をあまり信用してはいけません。

そこで、センサ位置で予測される誤差の大きさと観測誤差を足したもので割ります。

$$
K_i
=
\frac{B_{i,s}}{B_{s,s} + R}
$$

ここで `s` はセンサ位置です。

この式は、

```text
分子   : センサのズレを節点 i へ伝える強さ
分母   : センサ位置で見えるズレ全体の不確かさ
```

と読めます。

センサが信頼できる、つまり `R` が小さい場合は、分母が小さくなり、補正が強くなります。

センサが信頼できない、つまり `R` が大きい場合は、分母が大きくなり、補正が弱くなります。

複数センサになると、これを行列で同時に行います。

```text
1センサ:
  K_i = B_{i,s} / (B_{s,s} + R)

複数センサ:
  K = B H^T (H B H^T + R)^-1
```

ここで、

| 部分 | 役割 |
|---|---|
| `B H^T` | 各節点と各センサ位置の誤差のつながり |
| `H B H^T` | センサ位置同士で予測される誤差の共分散 |
| `H B H^T + R` | センサ空間で見た予測誤差と観測誤差の合計 |
| `(...)^{-1}` | 複数センサの重複や信頼度を考慮して重みを調整する |

つまり `K` は、`B` を観測空間に写し、観測誤差 `R` も考慮して作った「観測ズレの配分表」です。

## 8. OIゲインの導出

ここでは、なぜ

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

になるのかを整理します。

まず、真の状態を `x^true`、OpenFOAMの予測値を `x^f` とします。

予測誤差を、

$$
\mathbf{e}^f
=
\mathbf{x}^f
-
\mathbf{x}^{true}
$$

とします。

観測値は、

$$
\mathbf{y}
=
\mathbf{H}\mathbf{x}^{true}
+
\mathbf{v}
$$

です。ここで `v` は観測誤差です。

OIでは、解析値を次の線形な形で作ると仮定します。

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

このとき、イノベーションは、

$$
\mathbf{y}
-
\mathbf{H}\mathbf{x}^f
=
\mathbf{H}\mathbf{x}^{true}
+
\mathbf{v}
-
\mathbf{H}\mathbf{x}^f
$$

なので、

$$
\mathbf{y}
-
\mathbf{H}\mathbf{x}^f
=
-
\mathbf{H}\mathbf{e}^f
+
\mathbf{v}
$$

となります。

解析誤差は、

$$
\mathbf{e}^a
=
\mathbf{x}^a
-
\mathbf{x}^{true}
$$

です。上の更新式を代入すると、

$$
\mathbf{e}^a
=
\mathbf{e}^f
+
\mathbf{K}
\left(
-
\mathbf{H}\mathbf{e}^f
+
\mathbf{v}
\right)
$$

つまり、

$$
\mathbf{e}^a
=
\left(
\mathbf{I}
-
\mathbf{K}\mathbf{H}
\right)
\mathbf{e}^f
+
\mathbf{K}\mathbf{v}
$$

です。

ここで、予測誤差共分散と観測誤差共分散を、

$$
\mathbf{B}
=
E
\left[
\mathbf{e}^f
(\mathbf{e}^f)^T
\right]
$$

$$
\mathbf{R}
=
E
\left[
\mathbf{v}
\mathbf{v}^T
\right]
$$

とします。

また、予測誤差と観測誤差は相関しないと仮定します。

$$
E
\left[
\mathbf{e}^f
\mathbf{v}^T
\right]
=
\mathbf{0}
$$

このとき、解析誤差共分散は、

$$
\mathbf{A}
=
E
\left[
\mathbf{e}^a
(\mathbf{e}^a)^T
\right]
$$

であり、計算すると、

$$
\mathbf{A}
=
\left(
\mathbf{I}
-
\mathbf{K}\mathbf{H}
\right)
\mathbf{B}
\left(
\mathbf{I}
-
\mathbf{K}\mathbf{H}
\right)^T
+
\mathbf{K}
\mathbf{R}
\mathbf{K}^T
$$

となります。

OIでは、この解析誤差をできるだけ小さくする `K` を選びます。具体的には、解析誤差分散の総和である `trace(A)` を最小にします。

$$
\frac{\partial}{\partial \mathbf{K}}
\mathrm{tr}(\mathbf{A})
=
\mathbf{0}
$$

この条件から、

$$
\mathbf{K}
\left(
\mathbf{H}\mathbf{B}\mathbf{H}^T
+
\mathbf{R}
\right)
=
\mathbf{B}\mathbf{H}^T
$$

が得られます。

右から逆行列を掛けると、

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

になります。

これがOIゲインです。

直感的には、

```text
解析誤差が最も小さくなるように、
モデル予測と観測値の信頼度を B と R で重み付けしている
```

ということです。

## 9. OIはどの意味で「Optimal」なのか

OIの「Optimal」は、何でも完全に最適という意味ではありません。

次の仮定を置いたときに、解析誤差分散を小さくする線形推定として最適、という意味です。

- 予測誤差の平均は0
- 観測誤差の平均は0
- 予測誤差共分散 `B` が分かっている
- 観測誤差共分散 `R` が分かっている
- 補正は線形な形 `x_a = x_f + K(y - Hx_f)` で行う

この条件のもとで、

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

という形の推定を考え、解析誤差

$$
\mathbf{e}^a
=
\mathbf{x}^a - \mathbf{x}^{true}
$$

の分散が小さくなるように `K` を選ぶと、

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

が得られます。

つまりOIは、「背景値と観測値のどちらをどれだけ信じるか」を、誤差共分散 `B` と `R` から計算する方法です。

## 10. 背景誤差共分散 B

`B` は、OpenFOAM予測にどのような誤差の広がり方を仮定するかを表します。

今回の実装では、支配方程式から `B` を作るのではなく、距離に基づいて作っています。

$$
B_{ij}
=
\sigma_b^2
\exp
\left(
-
\frac{d_{ij}^2}{2L_c^2}
\right)
$$

| 記号 | 意味 |
|---|---|
| `sigma_b` | 背景誤差標準偏差 |
| `d_ij` | セル i とセル j の距離 |
| `L_c` | 空間相関長 |

現在の値は、

```text
BACKGROUND_ERROR_STD = 3.0 degC
CORRELATION_LENGTH_M = 0.09 m
```

です。

この仮定により、センサに近い節点ほど強く補正され、遠い節点ほど弱く補正されます。

### Bの各成分の意味

`B` の対角成分は、各節点の予測誤差の大きさを表します。

$$
B_{i,i}
=
\sigma_b^2
$$

非対角成分は、異なる節点同士の誤差がどれくらい一緒に動きやすいかを表します。

$$
B_{i,j}
$$

が大きいほど、節点 `j` の観測情報が節点 `i` にも伝わりやすくなります。

今回のように距離で `B` を作ると、

```text
近い節点  : 誤差が似ていると仮定する
遠い節点  : 誤差の関係は弱いと仮定する
```

という意味になります。

## 11. 観測誤差共分散 R

`R` はセンサ観測の不確かさを表します。

$$
\mathbf{R}
=
\sigma_{obs}^2 \mathbf{I}
$$

現在は2つのセンサを使っており、それぞれ独立な観測ノイズを仮定しています。

```text
OBS_NOISE_STD = 0.05 degC
```

`R` を大きくすると、観測をあまり信じなくなります。

`R` を小さくすると、観測を強く信じます。

ただし、観測ノイズを小さくしすぎると、ノイズまで過剰に追いかける可能性があります。

## 12. なぜ A_d が不要なのか

カルマンフィルタでは、予測誤差共分散を時間発展させるために、

$$
\mathbf{P}_k^-
=
\mathbf{F}_k
\mathbf{P}_{k-1}^+
\mathbf{F}_k^T
+
\mathbf{Q}
$$

のような式を使います。

ここで `F_k` は、OpenFOAMの時間発展を線形化したものです。これを厳密に作るには、OpenFOAMが内部でどのような離散化をしているかを知る必要があります。

今回のOIでは、この時間発展共分散 `P` を使いません。代わりに、固定した背景誤差共分散 `B` を使います。

```text
カルマンフィルタ: OpenFOAMの時間発展に対応する F や A_d が必要になりやすい
OI              : 固定した B を使うため A_d は不要
```

そのため、OpenFOAM側の方程式や離散化をPython側で再現しなくても使えます。

なお、`A_d` を使う方法が間違いというわけではありません。

`A_d` を使うカルマンフィルタ型の方法は、次のような場合には有効です。

- 支配方程式が分かっている
- 離散化方法も分かっている
- 状態数が比較的小さい
- 時間発展行列、またはその近似を作れる
- 線形化しても大きな破綻がない

たとえば、1次元熱伝導方程式を自分で差分法で解いている場合は、`A_d` を作りやすいため、カルマンフィルタを素直に適用できます。

一方でOpenFOAMのような汎用CFDソルバでは、メッシュ、境界条件、離散化スキーム、非線形項、ソルバ設定が複雑になります。その場合、Python側で `A_d` を作ると、

```text
Python側の A_d が、本当に OpenFOAM の時間発展を表しているのか
```

という問題が出ます。

今回のOIは、この問題を避けるために `A_d` を使わない仕様にしています。

## 13. 今回のプログラムとの対応

理論式とコードの対応は次の通りです。

| 理論 | コード |
|---|---|
| `H` | `build_oi_matrices()` 内で `H[i, node] = 1.0` |
| `B` | `BACKGROUND_ERROR_STD` と `CORRELATION_LENGTH_M` から作る距離相関行列 |
| `R` | `OBS_NOISE_STD**2` の対角行列 |
| `K` | `gain = B @ H.T @ inv(H @ B @ H.T + R)` |
| `y - Hx_f` | `innovation = y - H @ x_pred` |
| `x_a = x_f + Kd` | `x_da = x_pred + gain @ innovation` |

OpenFOAMが計算するのは `x_f` です。

```python
x_pred = ofi.read_field((k + 1) * DT)
```

Python側のOIが計算するのは `x_a` です。

```python
x_da = x_pred + gain @ innovation
```

そして、次ステップでは `x_da` をOpenFOAMの初期温度場として使います。

```python
x_current = x_da
```

## 14. データ同化手法の中でのOIの位置づけ

データ同化にはいくつかの流派があります。大きく分けると、次のように整理できます。

| 分類 | 代表的な手法 | 特徴 |
|---|---|
| 逐次同化 | OI, 3D-Var, Kalman Filter, EnKF | 時刻ごとに予測を観測で補正する |
| 時間窓同化 | 4D-Var | ある時間範囲全体で観測に合うように初期値やパラメータを調整する |
| パラメータ同化 | 最小二乗、Bayesian calibration, EnKFなど | 温度場ではなく、境界条件や物性値などを推定する |

今回のOIは、この中の **逐次同化** に入ります。

```text
OpenFOAMで1ステップ予測
        ↓
その時刻のセンサ観測で補正
        ↓
補正後の場を次ステップへ渡す
```

という流れです。

### OIは「固定Bを使う逐次同化」

OIの特徴は、背景誤差共分散 `B` を時間発展させず、あらかじめ決めた `B` を使うことです。

```text
OI:
  B は固定
  K = B H^T (H B H^T + R)^-1
```

これに対して、カルマンフィルタでは誤差共分散 `P` を時間発展させます。

```text
Kalman Filter:
  P をモデルで時間発展させる
  P_k^- = F P_{k-1}^+ F^T + Q
```

この `F` はモデルの時間発展を線形化したものです。OpenFOAMでこれを作るには、ソルバの中身や離散化にかなり踏み込む必要があります。

つまり、

```text
OI              : モデル内部をあまり知らなくても使いやすい
Kalman Filter   : モデルの線形化・時間発展行列が必要になりやすい
EnKF            : 行列Fは不要だが、OpenFOAMを多数ケース実行する必要がある
4D-Var          : 強力だが、随伴モデルや勾配計算が必要になりやすい
```

という位置づけです。

### OIと3D-Varの関係

OIは、3D-Varとも近い手法です。

3D-Varは次の評価関数を最小にする方法です。

$$
J(\mathbf{x})
=
\frac{1}{2}
(\mathbf{x}-\mathbf{x}^f)^T
\mathbf{B}^{-1}
(\mathbf{x}-\mathbf{x}^f)
+
\frac{1}{2}
(\mathbf{y}-\mathbf{H}\mathbf{x})^T
\mathbf{R}^{-1}
(\mathbf{y}-\mathbf{H}\mathbf{x})
$$

第1項は「OpenFOAM予測から離れすぎない」ことを表します。

第2項は「観測値から離れすぎない」ことを表します。

この式を線形観測 `H` のもとで解くと、OIと同じ形の更新式が得られます。

$$
\mathbf{x}^a
=
\mathbf{x}^f
+
\mathbf{B}\mathbf{H}^T
\left(
\mathbf{H}\mathbf{B}\mathbf{H}^T
+
\mathbf{R}
\right)^{-1}
\left(
\mathbf{y}
-
\mathbf{H}\mathbf{x}^f
\right)
$$

したがって、OIは

```text
固定した背景誤差共分散Bを使う、軽量な3D-Var的な逐次同化
```

と見ることもできます。

### 今回OIを選ぶ理由

今回の目的は、OpenFOAMソルバを改造せず、まず外部Pythonからデータ同化の流れを作ることです。

その目的に対して、OIはちょうどよい位置にあります。

| 観点 | OIが合う理由 |
|---|---|
| OpenFOAMの扱い | ブラックボックスの予測器として扱える |
| 実装量 | 少ない |
| 計算コスト | OpenFOAMを1本ずつ逐次実行すればよい |
| 説明しやすさ | `B`, `R`, `H` の意味が比較的分かりやすい |
| 発展性 | 次に EnKF やパラメータ同化へ進みやすい |

一方で、OIは最終形というより、OpenFOAMデータ同化を始めるための基礎的な方法です。

### 発展先

OIの次に考えられる発展は、目的によって変わります。

| 目的 | 発展先 |
|---|---|
| 温度場をより自然に補正したい | `B` の作り方を改善する |
| 熱入力や境界条件を推定したい | パラメータ同化にする |
| OpenFOAMの非線形性や複雑な誤差を扱いたい | EnKFにする |
| 時間履歴全体に合う初期値やパラメータを探したい | 4D-Varや最適化にする |

今回のサンプルでは、OpenFOAMソルバを改造せず、かつ `A_d` を作らないことを優先して OI を採用しています。

## 15. 限界

OIでは `B` を固定で与えます。そのため、補正の広がり方は `BACKGROUND_ERROR_STD` と `CORRELATION_LENGTH_M` に依存します。

また、現在は境界勾配 `G` 自体は推定していません。毎ステップのOpenFOAM予測では `G_MODEL = 50` の誤差が入り続けます。

今後さらに本格化するなら、次の方向があります。

- `G` も状態変数に含めて同化する
- `B` の相関長を検討する
- 複数ケースを使う Ensemble Kalman Filter に拡張する
