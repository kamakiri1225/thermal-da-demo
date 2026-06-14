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

## 出力ファイル

| ファイル | 内容 |
|---------|------|
| `img/anim_enkf_heatid.gif` | センサ数別の温度・導出熱源の同化過程 |
| `img/fig00_sensor_layout.png` | 固定センサ配置 |
| `img/fig01_rmse.png` | RMSE 収束 |
| `img/fig02_xpred_final.png` | 最終温度場 |
| `img/fig03_qest_final.png` | 最終温度場から導出した熱源場 |

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
| 同化手法 | 決定論的ESMDA / 最適内挿 |
| 事前共分散 | Matérn 3/2、相関長1.4 cells |
| 固定センサ数 m | 30, 100, 300, 600 |
| ESMDA表示反復数 | 30 |

固定点温度だけから900セルの任意の温度場を一意に決めることはできません。
この例ではMatérn共分散を仮定してセンサ間を補間します。最終RMSEは
30点で45.1 degC、100点で28.3 degC、300点で18.2 degC、
600点で11.9 degCです。輪郭を得るには十分なセンサ密度、時間発展する
物理モデル、または妥当な追加事前知識が必要です。

熱量図は、外れ値で色が潰れないよう全最終ケースの絶対値99%分位を
共通表示範囲に使います。各パネルのタイトルにはクリップ前の実際の
最小値・最大値も表示します。

## 関連デモ

- [`A003_miffy_enkfT_fixsensor`](../../A003_miffy_enkfT_fixsensor/docs/README.md) — 固定センサ数比較（ケース1）
- [`A02_miffy_enkf_T`](../../A02_miffy_enkf_T/docs/enkf_method.md) — 基本 EnKF（時間発展）
