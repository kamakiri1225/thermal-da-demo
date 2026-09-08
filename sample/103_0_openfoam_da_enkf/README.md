# 103_0 OpenFOAM × アンサンブルカルマンフィルタ(EnKF)によるデータ同化

`102_0_openfoam_hollow_cylinder_heat_transfer`(中空円筒を外周ヒータで温める
chtMultiRegionFoam・輻射あり CHT)をベースに、**数点の温度観測だけで、でたらめな
初期状態から固体温度場を真値へ補正していく時刻歴**を、EnKF(アンサンブルカルマン
フィルタ)で確認する双子実験(OSSE)。

姉妹フォルダ `103_1_openfoam_da_pf` は同じ問題を粒子フィルタ(PF)で解く。

## 何を推定するか

- **状態**: 固体温度場の全セル値(約2万セル)＋ ヒータ発熱 Q(拡大状態=パラメータ同時推定)
- **観測**: 固体表面近傍の数点温度(既定は hot・cold の2点)＋ 観測ノイズ 0.3 K
- **狙い**: でたらめな初期温度場＋誤った Q から、2点観測だけで温度場全体と Q が真値へ収束するか

## 2つの実装

| 実装 | 場所 | 速度 | 用途 |
|------|------|------|------|
| **OpenFOAM 連成版(本命)** | `daof/`, `run/run_openfoam_enkf.py` | 数十分〜 | 各メンバーが実際に chtMultiRegionFoam を回す field-space EnKF |
| **ROM 高速版(プロトタイプ)** | `dacore/`, `run/run_enkf.py` | 数秒 | 102_0 に校正した集中定数モデルでアルゴリズムを高速確認 |

どちらも解析コアは共通(`dacore/enkf.py`)。

## クイックスタート

```
pip install --user -r requirements.txt

# 1) まず ROM 版で挙動確認(数秒)
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 dacore/calibrate.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_enkf.py

# 2) OpenFOAM 連成版(要 OpenFOAM v2512、既存 102_0 メッシュ)
cd openfoam && sh setup_base_case.sh && cd ..
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  nohup python3 run/run_openfoam_enkf.py > openfoam/run_enkf.log 2>&1 &
```

## ドキュメント(`docs/`)

- `00_data_assimilation_methods.md` ― どんなデータ同化手法があるか、なぜ EnKF/PF か
- `01_twin_experiment_and_setup.md` ― 双子実験の問題設定・観測点・状態ベクトル
- `02_enkf_method.md` ― EnKF の理論と OpenFOAM 連成の実装・結果
- `03_run_procedure.md` ― 再現手順(ROM 版・OpenFOAM 版)

## ディレクトリ

```
dacore/     ROM(集中定数モデル)+ 解析コア(enkf.py/pf.py)+ 可視化
daof/       OpenFOAM 連成(場の入出力・メンバー実行・双子実験ドライバ)
openfoam/   base_case 生成スクリプト・DA設定・(実行時に生成されるケース群=Git管理外)
run/        実行スクリプト
config/     ROM 版の設定・校正結果・102_0 温度履歴
results/    RMSE 等の CSV・要約(軽量なので追跡)
docs/       解説と図
```

計算結果の重いファイル(メッシュ入りケース・時間ディレクトリ・中間場)は `.gitignore` 済み。
