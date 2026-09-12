# 105_0 センサ配置と熱感度 ― FrontISTRの感度行列 W=K⁻¹H で観測点を選ぶ

104（温度＋変位2点のデータ同化）の発展。**熱感度 = FrontISTR(KinvH)の感度行列
W=K⁻¹H**（DUMPWパッチ版fistr1がダンプするK・Hから構築）を本ケース円筒に適用し、
「感度が高い場所を観測に選ぶとデータ同化の精度が出る／低い場所では出ない」を
実際の同化計算で検証する。

## 問い（仮説）

- W=K⁻¹H は「節点温度→変位」のヤコビアン。**行**で読めば変位観測点の価値、
  **列**で読めば温度位置の影響が事前に分かるはず
- ならば「W行感度が高い2点」の変位観測は同化精度が出る／低い2点では出ない ― のでは？

## 方法

1. **熱感度（本命）**: `run/dumpw_cylinder.py`
   六面体メッシュを6分割テトラ化(23,040四面体)し、**DUMPWパッチ版fistr1**で
   `sensitivity_Wdiff.vtk`＋K・Hダンプを出力 → `run/build_W_from_dumps.py` で
   W=K⁻¹H全体(5040×5040)を構築（DUMPW出力と最大相対差1.9e-7で一致検証）
   → `run/kinvh_sensitivity.py` が行/列感度・観測点選定を確定
2. **センサ配置実験**（`run/sensor_placement_study.py`）
   温度センサ1点を hot/mid/cold/top/core の各候補に置き、
   (a) 温度1点のみ、(b) 温度1点＋変位2点（M演算子、104の手法）で
   ROMデータ同化を5seed実行。最終温度RMSEと、その位置のW列感度を突き合わせる。
3. **変位観測点の選定実験**（`run/displacement_point_selection.py` ほか）
   W行感度で高感度2点/低感度2点/現行/ランダムを選び分けてDA精度を比較。
   GIF・QoI時刻歴は `run/make_selection_gif.py` / `run/plot_qoi_points.py`
4. 交差検証: `run/sensitivity_map.py`（120パッチ摂動、W射影と最大差5.8%）

## 実行

```
# DUMPWパッチ版fistr1 (~/src/FrontISTR-dumpw) が必要
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/dumpw_cylinder.py          # 15秒
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/build_W_from_dumps.py      # 約2分
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/kinvh_sensitivity.py       # 数十秒
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensor_placement_study.py  # 数分
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/displacement_point_selection.py
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_selection_gif.py      # 長い(FrontISTR多数)
```

## 出力

- `docs/img/sensitivity_frontistr.png` … W行ノルム/W列の面表示（熱感度の本体）
- `docs/img/disp_point_sensitivity_3d.png` … W行感度と選定点（全点ラベル付き）
- `docs/img/sensor_placement_rmse.png` / `sensitivity_vs_rmse.png` … 実験1
- `docs/img/selection_*.png|gif` … 実験2（QoI時刻歴・温度/変形アニメ）
- `paraview/sensitivity_Wdiff_cylinder.vtk` … DUMPW直接出力（ParaView用）
- `results/kinvh_sensitivity.npz` / `sensor_placement.csv` / `disp_point_selection.csv`
  （`results/Wz_full.npy` ≈100MB は git 管理外・`build_W_from_dumps.py`で再生成）

考察は `docs/00_discussion.md`、証拠計算の詳細は `docs/02_evidence_sensitivity_selection.md`、
**実務ガイド（実験前のセンサ配置設計・資料作り用）は `docs/03_sensor_design_guide.md`**、
作業ログ（成功も失敗も）は `docs/90_worklog_103-105.md`。

## 前提

- 隣に `102_1_frontistr_hollow_cylinder_thermal_expansion`（FrontISTR連成コード）
- `config/` は104の校正結果（ROM・M演算子）を持ち込み済み。`fistr1` がPATHにあること
- DUMPWパッチ版: `~/src/FrontISTR-dumpw/build-dumpw/fistr1/fistr1`
  （パッチ: `20260810_KinvH/frontistr/patch/frontistr_dumpw_tet.patch`、四面体341専用）
