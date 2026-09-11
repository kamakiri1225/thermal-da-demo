# 何を・どこで計算し、いま104に何があるのか(作業の流れ)

「結局どこで何を計算して、この結果ができたのか」を最初から順に、正確にたどる。
成果物ファイルと、それを作った計算を1対1で対応させる。

---

## 全体の流れ(4段階)

```
[段階A] 102_0 (別フォルダ)          [段階B] 102_1 (別フォルダ)
  OpenFOAM CHT を実計算               FrontISTR 熱膨張を実計算
  → 固体の温度履歴                     → 温度→変位の関係
        │                                   │
        │ 温度履歴CSVを借用                  │ 変位データを借用
        ▼                                   ▼
[段階C] 104 の準備(校正)          config/ に校正結果を保存
        │
        ▼
[段階D] 104 の本番(データ同化)   2つのやり方で実行
        ├─ ROM版   (数秒, NumPyだけ)      → rom_fem_*.png
        └─ OpenFOAM版 (数十分, 実ソルバ)  → openfoam_fem_*.png ほか
```

段階Aと段階Bは **104ではなく隣のフォルダ**(102_0 / 102_1)で先に済ませた計算。
104はその結果を「借りて」使っている。

---

## 段階A: 102_0 で OpenFOAM CHT を実計算(過去に実施済み)

- 場所: `../102_0_openfoam_hollow_cylinder_heat_transfer`
- 何を: 中空円筒を外周ヒータで温める chtMultiRegionFoam(輻射あり)を 0→600s 計算
- 結果のうち104が使うもの:
  **温度履歴 CSV** → 104 の `config/openfoam_102_0_temperature_history.csv` にコピー済み

## 段階B: 102_1 で FrontISTR 熱膨張を実計算(過去に実施済み)

- 場所: `../102_1_frontistr_hollow_cylinder_thermal_expansion`
- 何を: 102_0 の温度場を入力に、線形静解析で熱膨張変位を 0→600s 計算
- 結果のうち104が使うもの:
  **変位の時刻歴**(`data/timehistory.csv`)を、104の校正で参照

> 104の `fem/fem_obs.py` と `run/*displacement*` は、102_1 の連成コード
> (`cylinder_mesh.py` / `fistr_case.py` など)を **import して再利用**している。
> つまり FrontISTR 本体の計算ロジックは 102_1 のものをそのまま使う。

---

## 段階C: 104 の校正(段階A・Bの結果を数式モデルに焼き込む)

104 で最初に1回だけ走らせる軽い計算。結果は `config/` に保存され、以降の入力になる。

| 実行コマンド | 何を計算 | 使った入力 | 生成物(104内) |
|--------------|----------|-----------|----------------|
| `python3 dacore/calibrate.py` | 5ノードROMを102_0の温度履歴に最小二乗フィット | `config/openfoam_102_0_temperature_history.csv`(段階A) | `config/rom_calibrated.yaml`(C,K,h、残差0.012K) |
| `python3 -m dacore.displacement` | 温度→上面変位の線形写像Dを102_1の変位に同定 | `../102_1/data/timehistory.csv`(段階B) | `config/displacement_operator.yaml`(残差0.003µm) |

**この段階の成果 = `config/rom_calibrated.yaml` と `config/displacement_operator.yaml`。**
「本物のOpenFOAM/FrontISTRの挙動を、軽い数式に写し取った校正データ」。

---

## 段階D: 104 の本番 = データ同化(EnKF)

でたらめな初期状態から、数点観測で真値に補正する双子実験。**2つのやり方で実施した。**

### D-1. ROM版(軽い・数秒・NumPyだけ) ― 実行済み

- 実行: `python3 run/run_rom_fem.py`
- 場所: **このマシンのPython**(OpenFOAM/FrontISTRは呼ばない)
- 中身: 60メンバー×30サイクル。前進=`dacore/cht_rom.py`、解析=`dacore/enkf.py`
- 段階Cの校正データを使うので、挙動は本物そっくり
- **いま104にある成果物**:
  - `docs/img/rom_fem_temperature.png`(5点温度の時刻歴 でたらめ/同化/真値)
  - `docs/img/rom_fem_displacement.png`(上面変位の時刻歴)
  - `docs/img/rom_fem_rmse.png`(誤差の収束)
  - `docs/img/da_vs_truth_temperature.gif` / `da_vs_truth_displacement.gif`(0→600s動画)
  - `docs/img/da_temperature_field3d.png` は下のD-2の場を使うが、時刻歴系はROM

### D-2. OpenFOAM版(重い・数十分・実ソルバ) ― 実行済み

- 実行: `sh openfoam/setup_base_case.sh` → `python3 run/run_openfoam_fem_enkf.py`
- 場所: **このマシンで chtMultiRegionFoam と fistr1 を実際に起動**
- 中身: **5メンバー×10サイクル、0→600s**(ヒータON 0-300s→OFF)。状態=固体全20696セル温度＋Q
  (初回は0→20s・2サイクルの最小実行で検証し、その後600sフルで再実行)
  - 前進=`daof/of_case.py`(chtMultiRegionFoam を subprocess 実行)
  - 変位観測=`fem/fem_obs.py`(FrontISTR を実行)
  - 解析=`dacore/enkf.py`(ROMと共通)
- **いま104にある成果物**:
  - `results/openfoam_fem_enkf_summary.yaml`(600s: **RMSE 7.61→0.011K, Q=14.93W**。20s版のバックアップは `*_20s.*`)
  - `results/openfoam_fem_enkf_history.csv`(サイクルごとの記録)
  - `docs/img/openfoam_fem_enkf_rmse.png` / `_Q.png` / `_disp.png` / `_field.png`
  - `openfoam/run_fem_enkf/fields/*.npy`(でたらめ/同化後/真値の温度場。Git管理外)

### D-3. D-2の温度場を使った追加の可視化 ― 実行済み

D-2が保存した温度場npyを入力に、後処理として実行した:

| 実行コマンド | 何を計算 | 生成物 |
|--------------|----------|--------|
| `python3 run/plot_da_temperature_field3d.py` | DA温度場を実メッシュに貼り断面表示 | `docs/img/da_temperature_field3d.png` |
| `python3 run/plot_da_displacement_field.py` | DA温度場→FrontISTR→変位分布 | `docs/img/da_displacement_field.png` |
| `python3 run/export_paraview.py` | ParaView用に書き出し | `paraview/*.vtu`, `*.pvd`(Git管理外) |

---

## いま104フォルダにある成果物の由来(早見表)

| ファイル | どの段階/実行で作られたか | 元になった計算 |
|----------|--------------------------|----------------|
| `config/openfoam_102_0_temperature_history.csv` | 段階A のコピー | 102_0 OpenFOAM 実計算 |
| `config/rom_calibrated.yaml` | 段階C `calibrate.py` | 102_0 履歴へのフィット |
| `config/displacement_operator.yaml` | 段階C `dacore.displacement` | 102_1 変位へのフィット |
| `docs/img/rom_fem_*.png`, `da_vs_truth_*.gif` | 段階D-1 `run_rom_fem.py` / `make_da_gifs.py` | ROM(NumPy)のEnKF |
| `results/openfoam_fem_enkf_*` | 段階D-2 `run_openfoam_fem_enkf.py` | chtMultiRegionFoam+FrontISTR の実EnKF |
| `docs/img/openfoam_fem_enkf_*.png` | 段階D-2 | 同上 |
| `openfoam/run_fem_enkf/fields/ensmean_t*.npy` | 段階D-2(Git管理外) | 解析後(補正済み)の平均温度場 |
| `openfoam/run_fem_enkf/fields/foremean_t*.npy` | 段階D-2(Git管理外) | **予報(補正前)の平均温度場**。解析でsolid/Tを上書きする前に保存し、補正量(イノベーション)の空間分布を後から可視化できる |
| `docs/img/da_temperature_field3d.png` | 段階D-3 | D-2の温度場を後処理 |
| `docs/img/da_displacement_field.png` | 段階D-3 | D-2の温度場→FrontISTR |
| `paraview/*`(Git管理外) | 段階D-3 `export_paraview.py` | D-2の温度場＋ROM時刻歴 |

---

## 一言でまとめると

1. **102_0(OpenFOAM)** と **102_1(FrontISTR)** で本物の物理計算を先に済ませた
2. その結果を **104のconfig/** に校正データとして焼き込んだ(段階C)
3. 104本番のデータ同化を、**軽いROM版(数秒)** と **重いOpenFOAM実機版(数十分)** の
   両方で実行した(段階D)
4. いま104にある図・CSV・yamlは、上のどれかの実行が出した結果(早見表のとおり)

再現手順そのものは [`05_folder_structure.md`](05_folder_structure.md) の「目的別」と
各 `run/` スクリプト冒頭のコメントに書いてある。
