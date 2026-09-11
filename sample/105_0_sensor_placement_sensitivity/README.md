# 105_0 センサ配置と熱感度 ― どこに温度センサを置くとデータ同化の精度が出るか

104（温度＋変位2点のデータ同化）の発展。**変位観測の熱感度分布**を FrontISTR で
実計算し、「熱感度が高い場所に温度センサを置くほどデータ同化の精度が出やすいか」を
実際の同化計算で検証する。

## 問い（仮説）

- 変位 uz は温度分布の重み付き積分。**どの場所の温度に敏感か**は一様ではないはず
- ならば「感度が高い場所」の温度センサは変位観測と強く連携し、同化精度が出る／
  「感度が低い場所」のセンサは情報が薄く、精度が出にくい ― のでは？

## 方法

1. **感度分布の実計算**（`run/sensitivity_map.py`）
   円筒を周方向12×軸方向10=120パッチに分割し、各パッチ+1KでFrontISTRを実行。
   2観測点（uz_heater / uz_opp）の応答 ∂uz/∂T(x) を分布として取得（線形なので厳密）。
2. **センサ配置実験**（`run/sensor_placement_study.py`）
   温度センサ1点を hot/mid/cold/top/core の各候補に置き、
   (a) 温度1点のみ、(b) 温度1点＋変位2点（M演算子、104の手法）で
   ROMデータ同化を5seed実行。最終温度RMSEと、その位置の熱感度を突き合わせる。

## 実行

```
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensitivity_map.py        # 数分(FrontISTR×120)
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/sensor_placement_study.py # 数十秒
```

## 出力

- `docs/img/sensitivity_map.png` … 変位2点の熱感度分布（円筒展開図）
- `docs/img/sensor_placement_rmse.png` … センサ位置ごとの同化精度（棒グラフ）
- `docs/img/sensitivity_vs_rmse.png` … 感度 vs 精度の散布図（仮説の検証）
- `results/sensitivity_uz.npz` / `sensor_placement.csv`

考察は `docs/00_discussion.md`（結果確定後に記載）。

## 前提

- 隣に `102_1_frontistr_hollow_cylinder_thermal_expansion`（FrontISTR連成コード）
- `config/` は104の校正結果（ROM・M演算子）を持ち込み済み。`fistr1` がPATHにあること
