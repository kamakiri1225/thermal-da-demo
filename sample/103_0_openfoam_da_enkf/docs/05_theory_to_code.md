# 数式で学ぶ EnKF / PF ― 理論と実装箇所の完全対応表

データ同化の手法を**数式で**理解し、その各項が**このリポジトリのどのファイル・どの関数**で
実装されているかを1対1で追えるようにした勉強用ノート。式番号を付けてコードと突き合わせる。

対象コード(103_0 / 103_1 で共通):

```
dacore/enkf.py          EnKF 解析ステップ(モデル非依存)
dacore/pf.py            粒子フィルタ解析ステップ(モデル非依存)
dacore/ensemble.py      拡大状態アンサンブルの生成・前進(ROM版)
dacore/observations.py  観測演算子 H・観測誤差 R・双子実験の観測生成(ROM版)
dacore/twin.py          ROM版 双子実験ドライバ
dacore/cht_rom.py       集中定数ROM(前進モデル M)
daof/of_twin.py         OpenFOAM版 双子実験ドライバ(field-space)
daof/of_case.py         OpenFOAM メンバーの実行・場の読み書き
daof/of_io.py           OpenFOAM フィールドファイル入出力
```

---

## 1. 状態空間モデルとベイズフィルタリング

データ同化はすべて次の状態空間モデルの上に立つ。

**システムモデル(前進モデル)** …(1):

$$x_k = M(x_{k-1}) + \eta_k$$

**観測モデル** …(2):

$$y_k = H x_k + \varepsilon_k, \quad \varepsilon_k \sim N(0, R)$$

- $x_k$: 時刻 $k$ の状態。本ケースでは**拡大状態**(§4)
- $M$: 前進モデル。ROM版は `dacore/cht_rom.py` の熱回路方程式、
  OpenFOAM版は **chtMultiRegionFoam そのもの**
- $H$: 観測演算子(状態→観測点の温度の抜き出し。線形)
- $R$: 観測誤差共分散(対角、$\sigma = 0.3\,\mathrm{K}$)

**実装箇所:**

| 式の要素 | ROM版 | OpenFOAM版 |
|---------|-------|-----------|
| 前進 $M$ | `dacore/cht_rom.py:150 integrate_ensemble()` (RK4) | `daof/of_case.py:56 run_window()` → chtMultiRegionFoam 実行 |
| $H$ の構築 | `dacore/observations.py:27 obs_matrix()` | `daof/of_twin.py:44 run_openfoam_twin()` 内 `H[np.arange(n_obs), obs_cells] = 1` |
| $R$ の構築 | `dacore/observations.py:38 obs_cov()` | 同上 `R = np.eye(n_obs) * noise**2` |

ベイズフィルタは「予報(時間発展)→解析(ベイズ更新)」を繰り返す …(3):

$$p(x_k \mid y_{1:k}) \propto p(y_k \mid x_k)\, p(x_k \mid y_{1:k-1})$$

この事後分布の近似の仕方が手法の分かれ目:
**EnKF はガウス近似**、**PF は粒子(重み付きサンプル)近似**。

---

## 2. アンサンブルカルマンフィルタ(EnKF)

### 2.1 カルマンフィルタの解析式(出発点)

線形・ガウスなら厳密解がカルマンフィルタ。

カルマンゲイン …(4):

$$K = P_f H^\top \left( H P_f H^\top + R \right)^{-1}$$

解析(平均の更新) …(5):

$$x_a = x_f + K \left( y - H x_f \right)$$

解析共分散 …(6):

$$P_a = (I - KH)\, P_f$$

問題: 状態次元 $n$ が大きいと $P_f$($n \times n$)を持てない
(本ケース $n \approx 20697$ → $P$ は約 $4 \times 10^8$ 要素)。

### 2.2 アンサンブル近似

$N$ 個のメンバー $\lbrace x_f^{(i)} \rbrace$ で予報分布を代表させ、**標本統計**で置き換える。

標本平均 …(7):

$$\bar{x}_f = \frac{1}{N} \sum_{i=1}^{N} x_f^{(i)}$$

状態アノマリ行列($n \times N$) …(8):

$$\delta X = \left[ x_f^{(1)} - \bar{x}_f, \; \dots, \; x_f^{(N)} - \bar{x}_f \right]$$

観測空間アノマリ($m \times N$) …(9):

$$\delta Y = H\, \delta X$$

標本共分散($P_f H^\top$ と $H P_f H^\top$ の近似) …(10),(11):

$$C_{zy} = \frac{\delta X\, \delta Y^\top}{N-1}, \qquad C_{yy} = \frac{\delta Y\, \delta Y^\top}{N-1}$$

ゲイン …(12):

$$K = C_{zy} \left( C_{yy} + R \right)^{-1}$$

ポイント: $P_f$ 本体は一切作らない。逆行列は**観測数** $m \times m$(本ケース $2 \times 2$)だけ。
これが「低ランク更新」であり、$n = 2$万でも $N = 5$ で動く理由。

### 2.3 摂動観測法(stochastic EnKF)

式(5)を全メンバーに同じ $y$ で適用すると解析アンサンブルの分散が式(6)より小さく
なりすぎる。そこで各メンバーに独立な**観測摂動**を与える …(13):

$$x_a^{(i)} = x_f^{(i)} + K \left( y + \varepsilon^{(i)} - H x_f^{(i)} \right), \quad \varepsilon^{(i)} \sim N(0, R)$$

こうすると解析アンサンブルの標本共分散の期待値が式(6)に一致する(Burgers et al. 1998)。

### 2.4 共分散インフレーション

$N$ が小さいと標本共分散は真の共分散を過小評価し、フィルタが観測を無視するようになる
(フィルタ発散)。対策として平均まわりに広げる …(14):

$$x_f^{(i)} \leftarrow \bar{x}_f + \lambda \left( x_f^{(i)} - \bar{x}_f \right), \quad \lambda \gtrsim 1$$

本ケースは $\lambda = 1.02 \sim 1.05$。

### 2.5 実装対応表(`dacore/enkf.py:21 enkf_update()`)

| 式 | コード(enkf.py 内) |
|----|--------------------|
| (14) インフレーション | `Zf = zbar + inflation * (Zf - zbar)` |
| (9) $\delta Y$ | `Yf = Zf @ H.T` → `dY = Yf - ybar` |
| (8) $\delta X$ | `dZ = Zf - zbar` |
| (10) $C_{zy}$ | `C_zy = dZ.T @ dY / (n_ens - 1)` |
| (11) $C_{yy}$ | `C_yy = dY.T @ dY / (n_ens - 1)` |
| (13) 観測摂動 | `eps = rng.multivariate_normal(zeros, R, size=n_ens)` |
| (12)(13) ゲインと更新 | `gain_T = np.linalg.solve(S, C_zy.T)` → `Za = Zf + innov @ gain_T` |

`np.linalg.solve(S, ...)` は $S^{-1}$ を陽に作らない数値的定石(逆行列より安定・高速)。

---

## 3. 粒子フィルタ(PF, ブートストラップ/SIR)

### 3.1 重要度サンプリング

事後分布を重み付き粒子で表す …(15):

$$p(x_k \mid y_{1:k}) \approx \sum_{i=1}^{N} w_k^{(i)}\, \delta\!\left(x - x_k^{(i)}\right), \qquad \sum_i w_k^{(i)} = 1$$

提案分布に前進モデルをそのまま使う(ブートストラップ)と、重み更新は尤度のみ …(16):

$$w_k^{(i)} \propto w_{k-1}^{(i)} \exp\!\left( -\tfrac{1}{2} \left(y - H x^{(i)}\right)^\top R^{-1} \left(y - H x^{(i)}\right) \right)$$

ガウス近似を一切していない点が EnKF との本質的な違い。

### 3.2 退化と有効サンプル数(ESS)

逐次更新すると重みが少数粒子に集中する(退化)。指標が有効サンプル数 …(17):

$$\mathrm{ESS} = \frac{1}{\sum_i \left(w^{(i)}\right)^2}, \qquad 1 \le \mathrm{ESS} \le N$$

$\mathrm{ESS} < 0.5N$ になったら**リサンプリング**(重みに比例して粒子を複製・淘汰)する。

### 3.3 系統リサンプリング

一様乱数 $u \sim U(0, 1/N)$ から等間隔の位置を作り、重みの累積分布に射影 …(18):

$$\mathrm{positions}_j = \frac{u + j}{N}, \quad j = 0, \dots, N-1$$

親のインデックスは累積重みが $\mathrm{positions}_j$ を超える最小の $i$。
層化されているため単純多項リサンプリングより分散が小さい。

### 3.4 正則化(ジッタ)

リサンプリング後は同一粒子のコピーだらけになる(粒子貧困化)。微小ノイズで多様性を回復 …(19):

$$x^{(i)} \leftarrow x^{(i)} + \xi^{(i)}, \quad \xi \sim N(0, h^2)$$

静的パラメータ($Q$ など)は前進で分布が広がらないため、ジッタが特に重要。

### 3.5 実装対応表(`dacore/pf.py`)

| 式 | コード |
|----|--------|
| (16) 対数尤度 | `pf.py:22 _log_likelihood()` — `d2 = einsum("ei,ij,ej->e", innov, Rinv, innov)` |
| 数値安定化 | `pf.py:40 pf_update()` — `logw -= logw.max()` してから `exp`(オーバーフロー防止の定石) |
| (17) ESS | `ess = 1.0 / np.sum(w ** 2)` |
| (18) 系統リサンプリング | `pf.py:30 systematic_resample()` — `positions = (rng.random() + arange(n)) / n` |
| (19) ジッタ | `pf_update()` 内 `Za[:, :n_nodes] += rng.normal(0, jitter_T, ...)`; OpenFOAM版の $Q$ ジッタは `daof/of_twin.py` の `resampled` 分岐 |

### 3.6 次元の呪い(本ケースで実際に起きたこと)

観測・状態の次元が上がると、式(16)の尤度は最良粒子とそれ以外で桁違いになり、
1粒子に重みが集中する。必要粒子数は状態の実効次元に対して**指数的**に増える
(Snyder et al. 2008)。

OpenFOAM版 PF($N=5$, 状態 $\approx 2$万次元)の実測: **1回目の解析で ESS = 1.0**
(`results/openfoam_pf_history.csv`)。全粒子が最良1体のコピーになり、
場 RMSE は 2.38 K で停滞、$Q$ は 20.2 W に固着した。
同条件の EnKF は 0.085 K / $Q = 15.00$ W。式(12)の低ランク**連続**更新(全メンバーを
観測方向へ引き寄せる)と、式(16)の**選択**のみの更新(良い粒子が既に存在しないと
補正できない)の差がそのまま数字に出ている。

---

## 4. 拡大状態によるパラメータ同時推定

未知パラメータ $\theta$(ヒータ発熱 $Q$ など)は状態に連結して同時推定する …(20):

$$z = \begin{bmatrix} x \\ \theta \end{bmatrix}, \qquad \theta_k = \theta_{k-1} + \zeta_k$$

観測は温度のみ($H$ の $\theta$ 列はゼロ)だが、式(10)の標本共分散に
**温度と $\theta$ の相関**が現れるため、$\theta$ も観測から更新される。
これがパラメータ推定の仕組み。

ランダムウォーク項 $\zeta$(パラメータジッタ)はアンサンブルの $\theta$ が
1点に潰れるのを防ぐ。

| 実装 | 場所 |
|------|------|
| ROM版 $z = [T(5), q_\mathrm{scale}, h]$ | `dacore/ensemble.py:24 init_ensemble()` / `I_QSCALE, I_H` |
| ROM版 $\zeta$ | `dacore/ensemble.py:48 forecast()` — `Z[:, I_QSCALE] += rng.normal(0, param_jitter_q, ...)` |
| OpenFOAM版 $z = [T(20696), Q]$ | `daof/of_twin.py:44` — `Z[i, :Nc] = read_solid_T(...)`, `Z[i, iQ] = Q[i]` |
| 物理範囲クリップ | `dacore/ensemble.py:40 clip_params()` / `of_twin.py` の `np.clip` |

---

## 5. 双子実験(OSSE)の数式

真値 $x^t$ を自分で作り、そこから観測を合成して手法を検証する。

真値ラン …(21):

$$x^t_k = M\!\left(x^t_{k-1};\, \theta_\mathrm{true}\right)$$

合成観測 …(22):

$$y_k = H x^t_k + \varepsilon_k, \quad \varepsilon_k \sim N(0, R)$$

評価指標 …(23):

$$\mathrm{RMSE}_k = \frac{1}{\sqrt{n}} \left\| \bar{x}_{a,k} - x^t_k \right\|$$

| 実装 | ROM版 | OpenFOAM版 |
|------|-------|-----------|
| (21) | `dacore/observations.py:57 generate_truth()` | `of_twin.py` 冒頭の truth ループ |
| (22) | `dacore/observations.py:67 make_observations()` | `obs_values[t1] = Tf[obs_cells] + rng_obs.normal(...)` |
| (23) | `dacore/twin.py:30 _rmse()` | `of_twin.py` 内 `rmse()` |
| サイクル制御 | `dacore/twin.py:35 run_twin()` | `daof/of_twin.py:44 run_openfoam_twin()` |

重要: 乱数を `rng_obs`(真値・観測用)と `rng`(フィルタ用)に分離し、同じ seed なら
EnKF と PF が**同一の観測**を同化するようにしている(公平な比較の要)。

---

## 6. 前進モデルの数式

### 6.1 ROM(集中定数熱回路) — `dacore/cht_rom.py`

各ノード $i$ のエネルギー収支 …(24):

$$C_i \frac{dT_i}{dt} = \sum_j K_{ij} \left(T_j - T_i\right) + q_i(t) - h \left(T_i - T_\mathrm{air}\right)$$

行列形 $\dot{T} = MT + b$ にして RK4 積分。

| 式の要素 | 実装 |
|---------|------|
| $K_{ij}$ ネットワーク | `cht_rom.py:78 conductance_matrix()` |
| $M, b$ | `cht_rom.py:99 system_matrix()` |
| $q_i(t)$ (0-300s ON) | `cht_rom.py:46 heater_power_W()` / `:118 heater_input()` |
| RK4(1本) | `cht_rom.py:125 integrate_single()` |
| RK4(アンサンブル一括, einsum) | `cht_rom.py:150 integrate_ensemble()` |
| パラメータ同定(least_squares) | `dacore/calibrate.py` — 102_0 実履歴への当てはめ、残差 RMSE 0.012 K |

### 6.2 OpenFOAM(本物のCHT) — `daof/`

$M$ = chtMultiRegionFoam(輻射 fvDOM・自然対流・固体伝導)。DA ループに必要な操作:

| 操作 | 実装 |
|------|------|
| 固体場→状態ベクトル | `of_io.py:62 read_internal_scalar()`(nonuniform List の正規表現パース) |
| 状態ベクトル→固体場 | `of_io.py:98 write_solid_T()`(internalField 書き出し＋ヒータ $Q$ を BC に埋め込み) |
| 観測セルの決定 | `of_case.py:27 nearest_cells()`(プローブ座標→最近傍セル) |
| 時間ウィンドウ実行 | `of_io.py:120 set_control_window()` + `of_case.py:56 run_window()` |
| 解析場の書き戻し→再開 | `of_twin.py` — `set_solid_state(m, t1, Za[i,:Nc], Q)` で t1 から次窓を再開 |

実装上の落とし穴(実際に踏んだ): `adjustTimeStep yes` のとき `writeControl runTime` では
終了時刻を正確に踏めず時刻ディレクトリが `2.001591` になる。
`writeControl adjustableRunTime` が正解(`of_io.py` 参照)。

---

## 7. さらに学ぶための文献

- Evensen (1994, 2003): EnKF の原典・レビュー
- Burgers, van Leeuwen & Evensen (1998): 摂動観測法の理論的正当化
- Anderson & Anderson (1999): 共分散インフレーション
- Gordon, Salmond & Smith (1993): ブートストラップ粒子フィルタの原典
- Snyder et al. (2008): PF の次元の呪いの定量的解析
- 淡路ほか『データ同化 観測・実験とモデルを融合するイノベーション』(京都大学学術出版会)
- 樋口編『データ同化入門』(朝倉書店): 粒子フィルタ中心の和書
