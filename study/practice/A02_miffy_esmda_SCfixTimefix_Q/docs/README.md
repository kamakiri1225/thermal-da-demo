# ミッフィー温度場の固定センサ ESMDA 復元

形状を知らない一様な初期温度分布から、固定温度センサの観測だけを
ESMDA（Ensemble Smoother with Multiple Data Assimilation）で同化し、
ミッフィー温度場を復元するデモです。

初期アンサンブルの平均は全セル 20 degC であり、ミッフィー形状を含みません。
センサ位置も真値を参照せず、格子座標だけを使った空間充填配置です。

## ドキュメント

| ファイル | 内容 |
|---------|------|
| [`enkf_heatid_method.md`](enkf_heatid_method.md) | 手法の詳細説明・精度改善策 |
| [`esmda_explanation.md`](esmda_explanation.md) | ESMDAの仕組みを図と数式で解説 |
| [`sensor_count_study.md`](sensor_count_study.md) | 必要なセンサ数をシミュレーションで評価 |

## 出力ファイル

| ファイル | 内容 |
|---------|------|
| `img/anim_enkf_heatid.gif` | センサ数別の温度・導出熱源の同化過程 |
| `img/fig00_sensor_layout.png` | 固定センサ配置 |
| `img/fig01_rmse.png` | RMSE 収束 |
| `img/fig02_xpred_final.png` | 最終温度場 |
| `img/fig03_qest_final.png` | 最終温度場から導出した熱源場 |
| `img/fig04_sensor_count_sweep.png` | センサ数と精度・追加効果 |

## 実行

```bash
cd python
python3 enkf_heatid.py
```

## 主な設定

| パラメータ | 値 |
|-----------|-----|
| グリッド | 30×30 (n=900) |
| 初期平均温度 | 全セル 20 degC |
| 同化手法 | 確率的ESMDA |
| アンサンブル数 | 300 |
| 初期相関スケール | 1.5 cells |
| 局所化半径 | 5 cells |
| 代表固定センサ数 | 30, 100, 200, 300 |
| センサ数スイープ | 10～400点、各3回 |
| ESMDA反復数 | 30 |

固定点温度だけから900セルの任意の温度場を一意に決めることはできません。
センサ数スイープの結果、概形には75点、実用的なバランスには150点、
高精度側には300点が目安です。300点から400点へ増やしてもRMSEは
19.82から19.78 degCにしか改善せず、現設定ではほぼ飽和しています。

熱量図は、外れ値で色が潰れないよう全最終ケースの絶対値99%分位を
共通表示範囲に使います。各パネルのタイトルにはクリップ前の実際の
最小値・最大値も表示します。

## 関連デモ

- [`A003_miffy_enkfT_fixsensor`](../../A003_miffy_enkfT_fixsensor/docs/README.md) — 固定センサ数比較（ケース1）
- [`A02_miffy_enkf_T`](../../A02_miffy_enkf_T/docs/enkf_method.md) — 基本 EnKF（時間発展）
- [`A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q`](../../A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q/docs/README.md) — 決定論的ESMDA/OI版
