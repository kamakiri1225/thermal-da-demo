# 002-1 FrontISTR round bar data assimilation

`002-1_laplacian_da_round_bar` と同じ丸棒データ同化を、FrontISTR の熱伝導解析で行うための準備フォルダです。

FrontISTR は `/home/kamakiri/local/frontistr/bin/fistr1` にビルド済みです。
このフォルダでは、その `fistr1` を使って丸棒の熱伝導ツイン実験と OI データ同化を実行します。

## 目的

OpenFOAM 版で行っている次の流れを FrontISTR に置き換えます。

```text
丸棒温度場を FrontISTR で1ステップ予測
  ↓
同化用センサ温度で OI 補正
  ↓
補正後温度場を次ステップ初期温度として書き戻す
  ↓
DAなし / OIあり / truth を比較してグラフ化
```

## 丸棒条件

OpenFOAM 版 `sample/002-1_laplacian_da_round_bar` と同じ条件を使う予定です。

| 項目 | 値 |
|---|---:|
| 丸棒長さ | 0.30 m |
| 丸棒直径 | 0.01 m |
| 軸方向分割 | 40 |
| 時間刻み | 1 s |
| 計算時間 | 900 s |
| 熱拡散率 | 6.6e-5 m2/s |
| 同化センサ | N1, N4 |
| 検証センサ | N28, N36 |

## 実行方法

```bash
cd /mnt/d/work/002_CAE/openfoam/20260505_datadoka/sample/002-1-FrontISTR_round_bar_da/oi
./run.sh
```

結果は次に出力されます。画像は `oi/results/img/` に集約します。

| 出力 | 内容 |
|---|---|
| `oi/results/img/results_frontistr_da.png` | truth / FrontISTR-only / OI の比較グラフ |
| `oi/results/img/measurement_points_frontistr_da.png` | 代表測定点の比較グラフ |
| `oi/results/img/all_axial_nodes_timeseries_frontistr_da.png` | 全 axial node の時刻歴 |
| `oi/results/img/axial_nodes_heatmap_frontistr_da.png` | 全 axial node の温度ヒートマップ |
| `oi/results/img/axial_nodes_error_heatmap_frontistr_da.png` | 全 axial node の誤差ヒートマップ |
| `oi/results/summary_rmse.csv` | 軸方向 RMSE |
| `oi/results/*/temperature_history.csv` | 節点温度履歴 |
| `oi/results/*/axial_temperature_history.csv` | 軸方向平均温度履歴 |
| `oi/results/vtk/*/*.vtk` | ParaView 用の 1 秒刻み VTK 連番 |

## 関連リンク

- [FrontISTR の比較画像](oi/results/img/results_frontistr_da.png)
- [FrontISTR の測定点比較](oi/results/img/measurement_points_frontistr_da.png)
- [OpenFOAM 版スライド](../002-1_laplacian_da_round_bar/slides.html)

## インストール

FrontISTR の導入方法は [INSTALL.md](INSTALL.md) を参照してください。

## 作業ログ

今回行ったことは [WORK_LOG.md](WORK_LOG.md) に記録しています。
