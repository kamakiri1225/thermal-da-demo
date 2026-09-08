# 103_1 OpenFOAM × 粒子フィルタ(PF)によるデータ同化

`102_0_openfoam_hollow_cylinder_heat_transfer` をベースに、**数点の温度観測だけで、
でたらめな初期状態から固体温度場を真値へ補正していく時刻歴**を、粒子フィルタ(PF)で
確認する双子実験(OSSE)。姉妹フォルダ `103_0_openfoam_da_enkf`(EnKF)と同一の
真値・観測を用いるので、両手法を直接比較できる。

## 何を推定するか

- **状態**: 固体温度場の全セル値(約2万セル)＋ ヒータ発熱 Q
- **観測**: 固体表面近傍の数点温度(既定は hot・cold の2点)＋ 観測ノイズ 0.3 K
- **狙い**: でたらめな初期温度場＋誤った Q から、2点観測だけで温度場と Q が真値へ収束するか、
  そして PF が高次元(場)でどう振る舞うか(次元の呪い)

## 2つの実装

| 実装 | 場所 | 速度 | 用途 |
|------|------|------|------|
| **OpenFOAM 連成版(本命)** | `daof/`, `run/run_openfoam_pf.py` | 数十分〜 | 各粒子が実際に chtMultiRegionFoam を回す field-space PF |
| **ROM 高速版(プロトタイプ)** | `dacore/`, `run/run_pf.py` | 数秒 | 集中定数モデルで PF を高速確認・EnKF と比較 |

解析コアは共通(`dacore/pf.py`)。

## クイックスタート

```
pip install --user -r requirements.txt

# 1) ROM 版(数秒)＋ EnKF との比較図
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 dacore/calibrate.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_pf.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/compare_enkf_pf.py

# 2) OpenFOAM 連成版(要 OpenFOAM v2512、既存 102_0 メッシュ)
cd openfoam && sh setup_base_case.sh && cd ..
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  nohup python3 run/run_openfoam_pf.py > openfoam/run_pf.log 2>&1 &
```

## ドキュメント(`docs/`)

- `00_data_assimilation_methods.md` ― データ同化手法の俯瞰(103_0 と共通)
- `01_twin_experiment_and_setup.md` ― 双子実験の問題設定(103_0 と共通)
- `02_particle_filter_method.md` ― PF の理論と OpenFOAM 連成の実装・結果
- `03_run_procedure.md` ― 再現手順
- `04_enkf_vs_pf.md` ― EnKF と PF の比較

計算結果の重いファイルは `.gitignore` 済み。
