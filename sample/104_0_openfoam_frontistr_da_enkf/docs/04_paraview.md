# ParaView で見る(温度分布・変位分布・時刻歴アニメーション)

`run/export_paraview.py` が、データ同化の結果を ParaView で開ける VTU / PVD に書き出す。

```
cd sample/104_0_openfoam_frontistr_da_enkf
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/export_paraview.py
```

出力(`paraview/`、再生成可能なので Git 管理外):

| ファイル | 中身 |
|----------|------|
| `temperature_fields.vtu` | 固体実メッシュ(20696セル)。セルデータ `T_garbage_degC / T_da_degC / T_truth_degC` |
| `displacement_fields.vtu` | 円筒メッシュ(5040節点)。点データ `U_garbage_m / U_da_m / U_truth_m`, `Uz_*_um` |
| `temperature_timehistory.pvd` | **0→600s の温度時刻歴アニメーション**(31ステップ)。`T_da / T_free / T_truth_degC` |

---

## 1. データ同化後の3D温度分布(静止画)

1. ParaView で `temperature_fields.vtu` を開く(Apply)
2. Coloring を **`T_da_degC`**(データ同化後)に
3. Filters → **Clip**(Plane, Y方向)で断面にすると内部分布が見える
4. `T_truth_degC` に切り替えると真値と一致することが確認できる。
   `T_garbage_degC` は一様(でたらめ)

## 2. データ同化された変位分布(変形形状)

1. `displacement_fields.vtu` を開く(Apply)
2. Filters → **Warp By Vector** → Vectors に **`U_da_m`**、Scale Factor を大きく(例 2000)
3. Coloring を **`Uz_da_um`** にするとヒータ側が持ち上がる傾きが見える
4. Vectors を `U_truth_m` に変えると真値の変形と一致。`U_garbage_m` は一様膨張で過大

## 3. 時刻歴アニメーション(0→600s)★

1. **`temperature_timehistory.pvd`** を開く(Apply)
2. Coloring を **`T_da_degC`** に、Clip で断面表示
3. 上部の **再生ボタン(▶)** で 0→600s を再生
   - ヒータ ON(0-300s)で温まり、OFF(300-600s)で冷える
   - `T_free_degC`(同化なし)に切り替えると、でたらめなまま真値からズレ続ける
   - `T_da_degC`(同化あり)は序盤で真値に吸い付き、以降ずっと一致

> 時刻歴の温度場は ROM(5ノード)を円筒メッシュへ逆距離補間して作っている。
> OpenFOAM 版のフル解像度時刻歴が要る場合は、`da_openfoam_config.yaml` の
> `t_end_s` を伸ばして `run/run_openfoam_fem_enkf.py` を回す(計算は重い)。

---

## 補足A: 色スケールを揃える

3つの配列(garbage/da/truth)を見比べるときは、ParaView の
**Rescale to Custom Range** で全ケース同じ範囲にすると比較しやすい。
でたらめは桁違いに大きいので、da と truth を比較するときは真値の範囲に合わせるとよい。

---

## 補足B: ROM(5点)なのに、なぜ3D分布/時刻歴を出せるのか

ROM は固体温度を **5点(hot/mid/cold/top/core)** しか持たない。そのままでは
「分布」にならない。ParaView 用の連続場は、**5点を円筒メッシュ全体へ逆距離加重
補間(IDW: Inverse Distance Weighting)して作り直している**。仕組みは3段階
(`run/export_paraview.py` の `export_timehistory()`)。

### 1. 表示用の器(メッシュ)を借りる

中身が空の、節点 5040 個・六面体要素の中空円筒メッシュを用意する
(`cylinder_mesh.build_cylinder_mesh`)。これは ROM 自身の格子ではなく、
温度を塗るための**入れ物**。

### 2. 5点 → 全節点へ逆距離加重補間(IDW)

各メッシュ節点 $p$ の温度を、ROM の5ノード温度 $T_i$(座標 $x_i$ は既知)から
距離の重みで合成する。近いノードほど強く効く、という単純な内挿:

$$T(p) = \frac{\sum_{i=1}^{5} w_i\, T_i}{\sum_{i=1}^{5} w_i},
\qquad w_i = \frac{1}{\lvert p - x_i \rvert^{2}}$$

節点が5ノードのどれかに一致するとき($\lvert p-x_i\rvert \to 0$)はその値をそのまま使う。

重みは各節点について時間に依らず一定なので、**正規化した重みを行列
$W$(5040行 × 5列)に1回だけ組んでおく**。各行 $W_{p,:}$ が「節点 $p$ を5ノードから
作るための重み」で、行和は 1:

$$W_{p,i} = \frac{w_i}{\sum_{j} w_j},
\qquad \sum_{i=1}^{5} W_{p,i} = 1$$

すると、ある時刻の5ノード温度ベクトル $\mathbf{T}_\mathrm{node}\in\mathbb{R}^5$ から
全節点温度 $\mathbf{T}_\mathrm{mesh}\in\mathbb{R}^{5040}$ は**行列積1発**で得られる:

$$\boxed{\ \mathbf{T}_\mathrm{mesh} = W\, \mathbf{T}_\mathrm{node}\ }$$

```python
# 重み行列 W (5040 x 5) を最初に1回だけ作る
W = np.zeros((len(node_coords), 5))
for i, p in enumerate(node_coords):
    d = np.linalg.norm(rom_nodes - p, axis=1)
    w = 1.0 / d**2
    W[i] = w / w.sum()

# 各時刻はこれだけ(行列積1発なので31ステップでも一瞬)
T_mesh = W @ T_node          # (5040,) = (5040,5) @ (5,)
```

これを `da / free / truth` それぞれに適用し、`T_da_degC` などの点データとして与える。

### 3. 各時刻を VTU に、束ねて .pvd 時系列に

0→600s の各サイクル(31ステップ)で $W\mathbf{T}$ を計算し `.vtu` を書き、
`.pvd` で時系列に束ねる。ParaView はこれをアニメーションとして再生する。

### 正直な限界

これは **5点からの内挿**なので、分布は「なめらかだが粗い」再構成。ヒータ直下の
急な勾配など細かい構造は5点では表現しきれない。**本物の高解像度分布**は、
20696セルを直接持つ **OpenFOAM 版**(`temperature_fields.vtu` /
`da_temperature_field3d.png`)が正確。

| | 温度の持ち方 | ParaView 分布の作り方 |
|---|---|---|
| **ROM 版** | 5点 | 5点を IDW で全メッシュに補間($\mathbf{T}_\mathrm{mesh}=W\mathbf{T}_\mathrm{node}$) |
| **OpenFOAM 版** | 20696セル | 実メッシュのセルにそのまま貼る |

ROM は「0→600s の時刻歴アニメが数秒で作れる」、OpenFOAM は「分布が正確」という住み分け。
