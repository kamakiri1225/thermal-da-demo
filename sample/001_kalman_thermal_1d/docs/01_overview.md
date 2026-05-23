# 001 カルマンフィルタによる熱データ同化デモ

このサンプルは、1次元熱伝導モデルとカルマンフィルタを使って、2つの温度センサから未観測節点の温度を推定するデータ同化のデモです。

現在の設定では、熱棒を **10節点** に分割し、温度センサは両端の **Node 0** と **Node 9** に置いています。センサのない Node 1 から Node 8 の温度は、熱モデルとカルマンフィルタで推定します。

## フォルダ構成

```text
001_kalman_thermal_1d/
├── main.py                 基本デモ。2点温度センサで全節点温度と熱変位を推定
├── main_compare_da.py      データ同化なし/ありの比較グラフを作成
├── main_disp_obs.py        変位観測を追加した場合の比較
├── thermal_model.py        1次元熱伝導棒モデル
├── kalman_filter.py        線形離散時間カルマンフィルタ
├── results.png             main.py の出力
├── results_compare_da.png  main_compare_da.py の出力
├── results_disp_obs.png    main_disp_obs.py の出力
└── docs/                   ドキュメント類
```

## ドキュメント

目次は [`00_index.md`](00_index.md) を参照。

- [`08_blog.md`](08_blog.md) — 2つの温度センサから熱棒全体の温度を推定する
- [`07_nextstep.md`](07_nextstep.md) — 次のステップ: laplacianFoamで行うソース改造なしデータ同化

## 実行方法

```bash
cd sample/001_kalman_thermal_1d
python main.py
python main_compare_da.py
python main_disp_obs.py
```

## 出力

- `results.png`: 10節点モデルに対するデータ同化の基本結果
- `results_compare_da.png`: データ同化なし/ありの比較
- `results_disp_obs.png`: 温度観測に変位観測を加えた場合の比較

グラフ内の凡例、軸ラベル、タイトルは、環境による文字化けを避けるため英語表記にしています。
