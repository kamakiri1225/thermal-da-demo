# 手法の詳細ウォークスルー ― 104 は具体的に何をどうやっているか

温度＋変位観測のデータ同化(EnKF)が、1サイクルで**どのデータがどのファイルの
どの関数を通って流れるか**を、実際の数値とともに追う。理論の数式は
`00_temp_displacement_da.md` と `../103_0_openfoam_da_enkf/docs/05_theory_to_code.md`。

---

## 全体像(1枚で)

```
                       ┌─────────────── 1 同化サイクル [t_k → t_{k+1}] ───────────────┐
 でたらめ初期          │                                                              │
 温度場+Q      ─────►  │  ① 前進(予報)        ② 観測演算子       ③ EnKF 解析         │ ──► 解析
 (5メンバー)           │  chtMultiRegionFoam   温度抜出+FrontISTR  標本共分散で更新     │     温度場+Q
                       │  で各メンバー前進     で予報観測 y_f      → 書き戻し           │     (真値へ近づく)
                       └──────────────────────────────────────────────────────────────┘
                            daof/of_case.py      daof/of_fem_twin.py   dacore/enkf.py
                                                 fem/fem_obs.py
```

---

## ステップ0: 準備(サイクル前)

**状態ベクトル**(1メンバー分、長さ 20697):

```
z = [ T_0, T_1, ..., T_20695 ,  Q ]
     └─── 固体全セル温度[K] ───┘  └ヒータ発熱[W]
```

**観測ベクトル**(長さ 4):

```
y = [ T_hot, T_cold ,  uz_heater, uz_opposite ]
     └ 温度2点[K] ┘    └ 上面変位2点[mm](FrontISTR) ┘
```

**観測誤差**(対角、単位混在は σ で正規化):
`R = diag(0.3², 0.3², 1e-4², 1e-4²)` → `daof/of_fem_twin.py` の `R = np.diag(...)`

**でたらめ初期状態**(実際の値、seed=20260908):
- メンバー温度 `T0 = [46.8, 25.4, 11.0, 17.5, 37.3] ℃`(真の 20 ℃ を全く知らない)
- ヒータ発熱 `Q = [8.4, 11.2, 6.3, 19.2, 6.2] W`(真値 15 W からずれ)

→ `of_fem_twin.py` の `rng.uniform(...)`。各メンバーの `0/solid/T` に書き込み
(`of_case.set_solid_state`)。

---

## ステップ1: 真値ランと合成観測(サイクルの前に1回)

真パラメータ(Q=15 W, 初期 20 ℃)で OpenFOAM を回し、真の温度場を得る。
その温度場を **FrontISTR に通して**真の変位も得る。両方にノイズを乗せて観測とする。

実測値(この実行):

| 時刻 | T_hot | T_cold | uz_heater | uz_opposite |
|------|-------|--------|-----------|-------------|
| t=10 s | 20.6 ℃ | 20.1 ℃ | +0.19 µm | −0.21 µm |
| t=20 s | 21.0 ℃ | 19.9 ℃ | +0.63 µm | −0.41 µm |

コード: `of_fem_twin.py` の truth ループ → `member_obs()` →
`of_case.read_solid_T()`(温度) + `fem_obs.displacement_obs()`(変位)。

---

## ステップ2: 予報(forecast)― 各メンバーを実際に OpenFOAM で回す

```python
for member in members:
    of_case.run_window(member, t_k, t_{k+1})   # chtMultiRegionFoam を直列実行
```

- `daof/of_case.py run_window()` → `subprocess` で `chtMultiRegionFoam`
- 起動直後は adaptive dt が小さく、1メンバー数分(輻射あり67kセル)
- でたらめに熱い/冷たい初期場から、ヒータで少し温まった温度場になる

---

## ステップ3: 観測演算子 ― 予報観測 y_f をサンプリング

各メンバーの予報温度場について、観測に対応する量を評価:

```python
T = of_case.read_solid_T(member, t_{k+1})          # 固体全セル温度
u = fem_obs.displacement_obs(member, t, fem_work)  # ← FrontISTR を回す
y_f = [ T[hot_cell], T[cold_cell], u[0], u[1] ]
```

- 温度2点は線形(セルの抜き出し)
- **変位2点は非線形**: 温度場 → IDW補間で FrontISTR 節点温度 → 線形静熱膨張解析 →
  上面 Uz。これが `fem/fem_obs.py displacement_obs()`(102_1 の `run_one_time` 再利用)。
  1メンバー約16秒。
- でたらめに +27 K 熱いメンバーは一様膨張で uz ~ +32 µm となり、真値 +0.6 µm から
  大きく外れる → 変位1点で温度場全体のオフセットが強く拘束される

予報観測を全メンバー分並べた行列 `Yf` (5×4) を作る。

---

## ステップ4: EnKF 解析 ― 標本共分散で状態を観測へ引き寄せる

```python
Za = enkf_update(Zf, y, H=None, R, rng, inflation=1.05, Yf=Yf)
```

`dacore/enkf.py enkf_update()` の中身(H 行列ではなく Yf を使う経路):

1. 状態アノマリ `dZ = Zf - mean` と観測アノマリ `dY = Yf - mean`
2. 標本共分散 `C_zy = dZ.T @ dY /(N-1)`、`C_yy = dY.T @ dY /(N-1)`
   ← ここに「温度場の各セル」と「変位観測」の相関が入る(これが未観測セルを動かす源)
3. ゲイン `K = C_zy (C_yy + R)^{-1}`(逆行列は 4×4 だけ)
4. 摂動観測付き更新 `Za = Zf + (y + ε - Yf) @ Kᵀ`

**キモ**: 変位は行列で書けない写像だが、`Yf` を各メンバーで**実際に計算**して渡すだけで、
温度観測とまったく同じ枠組みで同時同化できる。

---

## ステップ5: 書き戻し ― 次サイクルの再開点にする

```python
for member in members:
    of_case.set_solid_state(member, t_{k+1}, Za[member, :Nc], Za[member, Q])
```

解析した温度場と Q を `t_{k+1}/solid/T` に上書き。次サイクルは `startTime=t_{k+1}` で
ここから再開するので、補正が次の物理発展に反映される。温度場と Q だけ書き戻し、
流体場は各メンバー自身のものを使う(安定)。

---

## 1サイクルの計算コスト

```
① 前進:   chtMultiRegionFoam × 5メンバー   … 数分 × 5
③ 観測:   FrontISTR × 5メンバー(+真値1)   … 約16秒 × 6
```

OpenFOAM が支配的。FrontISTR は軽いので、変位観測を足しても総コストはほぼ変わらない。

---

## 出力(完走後)

| 図 | 内容 |
|----|------|
| `docs/img/openfoam_fem_enkf_disp.png` | 変位観測の予報 vs 真値(µm) |
| `docs/img/openfoam_fem_enkf_rmse.png` | 固体温度場 RMSE の収束 |
| `docs/img/openfoam_fem_enkf_Q.png` | ヒータ発熱 Q の推定 |
| `docs/img/openfoam_fem_enkf_field.png` | 温度場スナップショット(でたらめ→DA後→真値) |
| `results/openfoam_fem_enkf_summary.yaml` | 数値要約 |

103(温度のみ)と比べ、変位観測を足すことで温度場・Q の推定がどう変わるかを見る。
