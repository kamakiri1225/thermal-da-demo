# アンサンブルは何ケース・どこで計算し、平均と共分散をどう求めたか

実ソルバ版では、**初期温度と発熱量を変えた5個のOpenFOAMケースを計算し、その5個の結果から平均・共分散を求めてEnKFで更新**している。真値を作るケースは別に1個あり、平均・共分散の計算には入れない。

以下のパスは `sample/104_0_openfoam_frontistr_da_enkf/` を基準とする。設定・プログラムに加え、保存済みのケースフォルダと `openfoam/run_fem_enkf.log` を確認した内容である。

## 1. 何ケース計算したか

| 項目 | OpenFOAM＋FrontISTR版 | 軽量なROM版 |
|---|---|---|
| 同化用メンバー数 | **5** | **60** |
| 真値 | 別の1ケース | 別の1軌道 |
| 同化なしとの比較 | このドライバには別の同化なしアンサンブルはない | 同じ初期値の60メンバーを別途前進 |
| 時間範囲 | 0〜20秒 | 0〜600秒 |
| 同化時刻 | 10秒、20秒の2回 | 20秒ごと、30回 |
| 1メンバーの推定量 | 固体20,696セル温度＋発熱量Q | 5ノード温度＋発熱倍率q＋放熱h |
| 同化に使う観測 | 温度2成分＋変位2成分 | **温度2成分のみ**。変位は温度から予測 |
| 設定 | [openfoam/da_openfoam_config.yaml](../openfoam/da_openfoam_config.yaml) | [config/da_config.yaml](../config/da_config.yaml) |

「5メンバー×2サイクル」は、毎回新しい10ケースを作る意味ではない。同じ5ケースを0→10秒まで計算し、同化してから10→20秒へ継続する。

この実ソルバ版の同化実行でのソルバ呼び出し回数は次のとおり。準備・別の可視化スクリプトによる計算は含めない。

| 対象 | OpenFOAMの区間計算 | FrontISTRの静解析 |
|---|---:|---:|
| 推定用5メンバー | 5×2＝10回 | 5×2＝10回 |
| 真値1ケース | 1×2＝2回 | 1×2＝2回 |
| 合計 | **12回** | **12回** |

## 2. どこで、どの順に計算したか

入口は [run/run_openfoam_fem_enkf.py](../run/run_openfoam_fem_enkf.py) の `main()`。ここで設定を読み、`WORKDIR` を `openfoam/run_fem_enkf/` に定め、`run_fem_twin(cfg, WORKDIR)` を呼ぶ。

```text
openfoam/
├── base_case/                  共通の元ケース
├── run_fem_enkf.log            全体の進行・各メンバーの値
└── run_fem_enkf/
    ├── truth/                  真値専用（アンサンブルに含めない）
    ├── member_00/              推定用メンバー0
    ├── member_01/              推定用メンバー1
    ├── member_02/              推定用メンバー2
    ├── member_03/              推定用メンバー3
    ├── member_04/              推定用メンバー4
    └── fields/                 同化後の平均温度場など
```

各 `member_XX/` の中に以下がある。

```text
0/solid/T                      そのメンバーの初期温度・発熱量
10/solid/T                     10秒の温度・発熱量（同化後に上書き）
20/solid/T                     20秒の温度・発熱量（同化後に上書き）
log.run_0_10                   OpenFOAMの0→10秒ログ
log.run_10_20                  OpenFOAMの10→20秒ログ
fem_t10/                       10秒の予報温度から計算したFrontISTRケース
fem_t20/                       20秒の予報温度から計算したFrontISTRケース
```

実際のループは [daof/of_fem_twin.py](../daof/of_fem_twin.py) の `run_fem_twin()` にある。

1. `truth/` を0→10→20秒と前進させ、合成観測を作る。
2. 5メンバーを用意する。
3. 0→10秒について、`member_00` のOpenFOAM→FrontISTR、次に `member_01` … `member_04` と計算する。
4. **5個そろってから**平均・共分散を計算し、5メンバーそれぞれを同化する。
5. 同化した各メンバーを再開点として、10→20秒でも同じ処理を行う。

[daof/of_case.py](../daof/of_case.py) の `run_window()` は、次のように各ケースを作業ディレクトリにして起動し、終了を待つ。

```python
subprocess.run(["chtMultiRegionFoam"], cwd=member_dir,
               stdout=lf, stderr=subprocess.STDOUT)
```

したがって、**5メンバーは順番に実行される**。この処理は `mpirun` で5ケースを並列実行するものではない。スレッド数の環境変数とアンサンブルのケース数も別である。

## 3. 5メンバーは何を変えたのか

`run_fem_twin()` で `rng = np.random.default_rng(seed + 1)` とし、設定の `seed=20260908` に対して初期化用シードは **20260909**。温度を280〜320 K、Qを5〜30 Wの一様乱数で各5個生成する。観測生成には別の乱数生成器 `rng_obs`（シード20260908）を使う。

| メンバー | 初期温度［℃］ | 初期Q［W］ |
|---|---:|---:|
| member_00 | 46.8108 | 8.4225 |
| member_01 | 25.4481 | 11.1967 |
| member_02 | 10.9867 | 6.2731 |
| member_03 | 17.5002 | 19.1881 |
| member_04 | 37.2997 | 6.2407 |
| 真値（平均から除外） | 20.0000 | 15.0000 |

初期温度は、**メンバーごとに異なるが、そのメンバー内では固体全セル一様**。20,696セルを独立な乱数で初期化しているわけではない。共通メッシュを `prepare_member()` で複製し、`set_solid_state()` で温度とQを書き込む。前進すると、各ケース内に温度分布が生じる。

## 4. 計算結果をどの配列に集めたか

`run_fem_twin()` は各サイクルで以下を確保する。**行がメンバー、列が物理量**である。

| 配列 | サイズ | 1行の中身 |
|---|---|---|
| `Z`（EnKF内では `Zf`） | 5×20,697 | 全セル温度20,696個［K］とQ［W］ |
| `Yf` | 5×4 | hot温度、cold温度［K］、ヒータ側Uz、反対側Uz［mm］ |
| `y` | 4 | その時刻の合成観測。全メンバー共通 |
| `Za` | 5×20,697 | 同化後の各メンバーの状態 |

コードの該当箇所は以下。

```python
Z = np.empty((N, n_aug))
Yf = np.empty((N, n_obs))
for i, m in enumerate(members):
    ok, log = of_case.run_window(m, t0, t1)
    # 実コードではここで失敗判定を行う
    yv, Tm = member_obs(m, t1, os.path.join(m, f"fem_t{t1:g}"))
    Z[i, :Nc] = Tm
    Z[i, iQ] = Q[i]
    Yf[i] = yv
```

`member_obs()` は同じファイル内の関数。温度は最近傍セルから取り出し、変位は [fem/fem_obs.py](../fem/fem_obs.py) の `displacement_obs()` を通じてFrontISTRを実行する。

変位の2成分は厳密には各側の上面節点群の平均Uzである。この**1ケース内の節点平均**は、これから説明する**5ケース間のアンサンブル平均**とは別の処理。節点平均とm→mm変換は、隣接ケースの [run_thermal_expansion.py](../../102_1_frontistr_hollow_cylinder_thermal_expansion/python/run_thermal_expansion.py) の `run_one_time()` が担当する。

## 5. 平均・分散・共分散はどう計算するか

実装は [dacore/enkf.py](../dacore/enkf.py) の **`enkf_update()`**。OpenFOAMやFrontISTRの内部ではなく、Python／NumPy側で5ケースの結果をまとめて計算する。

### 5.1 最初に、平均まわりのばらつきを1.05倍する

設定の `inflation: 1.05` に従い、次の処理を行う。

```python
zbar = Zf.mean(axis=0, keepdims=True)
Zf = zbar + inflation * (Zf - zbar)
ybar_i = Yf.mean(axis=0, keepdims=True)
Yf = ybar_i + inflation * (Yf - ybar_i)
```

平均は維持し、各メンバーと平均の差を1.05倍する。したがって、後で計算する共分散は、拡大前に比べて **1.05²＝1.1025倍**になる。実装は渡された `Yf` の偏差も同じ倍率で拡大し、この時点でFrontISTRを再実行はしない。

### 5.2 各セル・各観測成分について、5ケースの平均を取る

以下の式では、拡大処理後の状態・予報観測をそれぞれ $z^{(i)},y_f^{(i)}$ と書く。

$$
\bar z_j=\frac{1}{5}\sum_{i=1}^{5}z_j^{(i)},\qquad
\bar y_a=\frac{1}{5}\sum_{i=1}^{5}y_{f,a}^{(i)}
$$

```python
zbar = Zf.mean(axis=0)  # 長さ20,697：同じセル・同じQ列を5ケースで平均
ybar = Yf.mean(axis=0)  # 長さ4：同じ観測成分を5ケースで平均
dZ = Zf - zbar
dY = Yf - ybar
```

例えば `zbar[100]` は「セル100の温度を5ケースで平均した値」。円筒全体の空間平均でも、時間平均でもない。`dZ` と `dY` は各メンバーがその平均からどれだけ外れたかを表す。

### 5.3 平均からの差の積を足し、N−1＝4で割る

$$
(C_{zy})_{ja}=\frac{1}{4}\sum_{i=1}^{5}
(z_j^{(i)}-\bar z_j)(y_{f,a}^{(i)}-\bar y_a)
$$

$$
(C_{yy})_{ab}=\frac{1}{4}\sum_{i=1}^{5}
(y_{f,a}^{(i)}-\bar y_a)(y_{f,b}^{(i)}-\bar y_b)
$$

```python
C_zy = dZ.T @ dY / (n_ens - 1)  # 20,697×4
C_yy = dY.T @ dY / (n_ens - 1)  # 4×4
```

- `C_zy`：各セル温度・Qと、各観測成分が一緒にどう変わるか。未観測セルやQを補正する手がかりになる。
- `C_yy` の対角：温度2成分と変位2成分、それぞれの**標本分散**。
- `C_yy` の非対角：異なる観測成分同士の**標本共分散**。

セルjの標本分散を単独で求めるなら `sum(dZ[:, j]**2)/4` だが、この更新コードは状態の全共分散 `dZ.T @ dZ / 4`（20,697×20,697）を作らない。必要な `C_zy` と `C_yy` だけ計算する。

5メンバーから作る偏差行列のランクは最大4である。したがって、20,696セル分の独立した不確かさを網羅しているわけではなく、この5ケースが表す変動を使う最小規模の検証である。

### 5.4 観測誤差を加え、5メンバーをそれぞれ補正する

観測誤差共分散は `run_fem_twin()` で設定値から作る。

$$
R=\operatorname{diag}(0.30^2,\ 0.30^2,\ (10^{-4})^2,\ (10^{-4})^2)
$$

温度成分の分散はK²、変位成分はmm²。これは5ケースから推定する `C_yy` と異なり、想定した測定誤差である。

$$
K=C_{zy}(C_{yy}+R)^{-1},\qquad
z_a^{(i)}=z_f^{(i)}+K\bigl(y+\varepsilon^{(i)}-y_f^{(i)}\bigr),
\quad\varepsilon^{(i)}\sim\mathcal N(0,R)
$$

実際には逆行列を明示的に作らず、4×4の連立方程式を解く。

```python
S = C_yy + R
eps = rng.multivariate_normal(np.zeros(len(y)), R, size=n_ens)
innov = (y[None, :] + eps) - Yf
gain_T = np.linalg.solve(S, C_zy.T)
Za = Zf + innov @ gain_T
```

共通観測 `y` 自体に加えた模擬測定ノイズと、ここで各メンバーに加える `eps` は別である。後者は確率的EnKFの観測摂動であり、5行×4成分を生成する。

## 6. 同化後の平均と、表示する標準偏差

`run_fem_twin()` に戻った後、温度を250〜400 K、Qを0〜60 Wにクリップし、**各メンバーの `Za[i]` をそのメンバーへ書き戻す**。5ケース全部を平均値に置き換える処理ではない。

記録用の平均・標準偏差は次の計算である。

```python
Q = Za[:, iQ].copy()
ens_mean = Za[:, :Nc].mean(axis=0)
q_mean = Q.mean()
q_std = Q.std()  # NumPy既定のddof=0
uz_mean = Yf[:, n_t:].mean(axis=0)
uz_std = Yf[:, n_t:].std(axis=0)
```

**EnKFの共分散はN−1で割るが、表示用の `.std()` はNで割った分散の平方根**である。

$$
\bar Q=\frac{1}{5}\sum_i Q_i,\qquad
\sigma_{Q,\mathrm{display}}=\sqrt{\frac{1}{5}\sum_i(Q_i-\bar Q)^2}
$$

初期Qの具体例では、平均 **10.264224 W**、表示用分散 **23.202949 W²**、表示用標準偏差 **4.816944 W**。同じ初期QをN−1で割った標本分散は **29.003686 W²**になる。これは初期Qによる計算例であり、インフレーション後の共分散そのものではない。

| 時刻［s］ | Q平均［W］ | Q表示用標準偏差［W］ | 平均温度場のRMSE［K］ |
|---|---:|---:|---:|
| 0 | 10.2642 | 4.81694 | 7.60909 |
| 10（同化後） | 16.3594 | 2.94875 | 0.0308066 |
| 20（同化後） | 14.4442 | 0.878855 | 0.0202975 |

出典：[results/openfoam_fem_enkf_history.csv](../results/openfoam_fem_enkf_history.csv)。RMSEはアンサンブルのばらつきではなく、**同化後の平均温度場と真値との差**を20,696セルで二乗平均して平方根を取ったもの。

変位図の `uz_mean` / `uz_std` は、同化前にFrontISTRで求めた**予報**の `Yf` から計算する。`enkf_update()` 内のインフレーションは新しい配列への代入なので、呼び出し元で変位表示に使う `Yf` は拡大前のまま。温度場・Qの同化後統計とは段階が異なる。0秒の変位統計はコードでゼロを入れた初期表示値であり、5ケースの初期変位を計算した値ではない。

## 7. 何が保存され、何が保存されないか

| 場所 | 確認できる内容 |
|---|---|
| `openfoam/run_fem_enkf.log` | 初期5メンバーの値、各区間の予報温度・変位の一部、同化後RMSE・Q統計（丸め値） |
| `openfoam/run_fem_enkf/member_XX/` | 各ケースの設定・場・OpenFOAMログ |
| 各ケースの `fem_t10/`, `fem_t20/` | 各予報温度から作ったFrontISTR入力・結果 |
| `openfoam/run_fem_enkf/fields/ensmean_t0.npy`, `ensmean_t10.npy`, `ensmean_t20.npy` | 初期および同化後のセルごとの平均温度 |
| `results/openfoam_fem_enkf_history.csv` | Q平均・標準偏差と温度RMSEなど |
| `results/openfoam_fem_enkf_summary.yaml` | メンバー数・最終RMSE・最終Qなど |

`C_zy`、`C_yy`、ゲイン、観測摂動は現在のコードではファイル保存していない。また、各時刻の `solid/T` は同化後に上書きするため、そのファイルだけから同化前の全 `Zf` と共分散を厳密には再現できない。共分散計算そのものは `enkf_update()` の上記コードで確認する。

なおCSVの `rmse_obs_K` は、実行スクリプトが描画コードとの互換用に `rmse_field` を複写した列であり、観測点だけで別途計算したRMSEではない。

## 8. 関連プログラムを追う順番

| 順 | ファイル | 関数・設定 | 担当 |
|---|---|---|---|
| 1 | [openfoam/da_openfoam_config.yaml](../openfoam/da_openfoam_config.yaml) | `ensemble`, `experiment`, `filter` | 5メンバー、時間、初期範囲、誤差、拡大倍率 |
| 2 | [run/run_openfoam_fem_enkf.py](../run/run_openfoam_fem_enkf.py) | `main()`, `WORKDIR` | 設定読み込み・開始・結果保存・描画 |
| 3 | [daof/of_fem_twin.py](../daof/of_fem_twin.py) | `run_fem_twin()`, 内部の `member_obs()` | 真値生成、5ケースの生成・前進ループ、配列集約、EnKF呼び出し、統計 |
| 4 | [daof/of_case.py](../daof/of_case.py) | `prepare_member()`, `run_window()`, `read_solid_T()`, `set_solid_state()` | ケース複製、OpenFOAM起動、温度・Qの読み書き |
| 5 | [daof/of_io.py](../daof/of_io.py) | `read_internal_scalar()`, `write_solid_T()`, `set_control_window()` | OpenFOAMファイル形式での入出力、開始・終了時刻設定 |
| 6 | [fem/fem_obs.py](../fem/fem_obs.py) | `displacement_obs()` | 各メンバーの温度場から変位2成分を取得 |
| 7 | [102_1のrun_thermal_expansion.py](../../102_1_frontistr_hollow_cylinder_thermal_expansion/python/run_thermal_expansion.py) | `run_one_time()` | 温度補間、FEM入力生成、構造解析、上面節点平均 |
| 8 | [102_1のfistr_case.py](../../102_1_frontistr_hollow_cylinder_thermal_expansion/python/fistr_case.py) | `run_fistr()` | FrontISTR実行 |
| 9 | [dacore/enkf.py](../dacore/enkf.py) | **`enkf_update()`** | **平均・偏差・共分散・ゲイン・5メンバーの更新** |

## 9. ROM版の60メンバーはどこにあるか

ROM版はケースフォルダを60個作らず、**Pythonのメモリ上に60行×7列の配列**として持つ。

- [run/run_rom_fem.py](../run/run_rom_fem.py) の `main()` が入口。
- [dacore/twin_fem.py](../dacore/twin_fem.py) の `run_rom_fem_twin()` が、同化用 `Zda` と比較用 `Zfree` を作り、30サイクル回す。真値は別の軌道で、平均に含めない。
- [dacore/ensemble.py](../dacore/ensemble.py) の `init_ensemble()` が60メンバーの初期値を生成し、`forecast()` が前進する。ROM版では各メンバーの5ノード温度も個別に乱数化する。
- [dacore/cht_rom.py](../dacore/cht_rom.py) の `integrate_ensemble()` が60メンバーを配列として時間積分する。毎サイクルOpenFOAMやFrontISTRを起動する処理ではない。
- 更新は共通の [dacore/enkf.py](../dacore/enkf.py) の `enkf_update()`。`Zf` は60×7、`Yf=Zf @ H.T` は60×2、`C_zy` は7×2、`C_yy` は2×2。共分散の分母は **60−1＝59**、偏差拡大倍率は **1.02**。
- 表示する温度は `Zda[:, :N_NODES].mean(0)` で60メンバー平均を取り、変位はその平均温度から校正済み線形写像で計算する。変位の模擬観測点は図示するが、EnKF更新には使わない。

実ソルバ版の5ケースとROM版の60メンバーは、別の設定・別の実験である。
