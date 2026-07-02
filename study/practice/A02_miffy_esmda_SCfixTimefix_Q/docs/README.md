# ミッフィー温度場の固定センサ ESMDA 復元

形状を知らない一様な初期温度分布から、固定温度センサの観測だけを
ESMDA（Ensemble Smoother with Multiple Data Assimilation）で同化し、
ミッフィー温度場を復元するデモです。

初期アンサンブルの平均は全セル 20 degC であり、ミッフィー形状を含みません。
センサ位置も真値を参照せず、格子座標だけを使った空間充填配置です。

スクリプトは 2 本あります。

| スクリプト | 状態変数 | 熱源の扱い |
|---|---|---|
| `python/enkf_heatid.py` | 温度場 x | 事後に `q = -α·L·x` で逆算（ノイズ増幅で図として読めない） |
| `python/esmda_qfull.py` | **熱源場 q** | 順モデル `x = L⁻¹(−q/α)` を組み込んだ完全逆問題として直接推定 |

## ドキュメント

| ファイル | 内容 |
|---------|------|
| [`enkf_heatid_method.md`](enkf_heatid_method.md) | 温度場版の手法詳細・精度改善策 |
| [`esmda_qfull_method.md`](esmda_qfull_method.md) | 完全逆問題版の手法・チューニング記録・結果解釈 |
| [`esmda_explanation.md`](esmda_explanation.md) | ESMDAの仕組みを図と数式で解説 |
| [`sensor_count_study.md`](sensor_count_study.md) | 必要なセンサ数をシミュレーションで評価（温度場版） |

## 出力ファイル

温度場版（`enkf_heatid.py`）:

| ファイル | 内容 |
|---------|------|
| `img/anim_enkf_heatid.gif` | センサ数別の温度・導出熱源の同化過程 |
| `img/fig00_sensor_layout.png` | 固定センサ配置 |
| `img/fig01_rmse.png` | RMSE 収束 |
| `img/fig02_xpred_final.png` | 最終温度場 |
| `img/fig03_qest_final.png` | 最終温度場から導出した熱源場 |
| `img/fig04_sensor_count_sweep.png` | センサ数と精度・追加効果 |

完全逆問題版（`esmda_qfull.py`）:

| ファイル | 内容 |
|---------|------|
| `img/anim_esmda_qfull.gif` | 温度（上段）と熱源（下段）の同化過程 |
| `img/fig05_qfull_rmse.png` | q / 温度 RMSE 収束（旧方式の水準線付き） |
| `img/fig06_qfull_qest_final.png` | 最終推定熱源場 |
| `img/fig07_qfull_xrecon_final.png` | 推定熱源からの順解析温度場 |
| `img/fig08_qfull_vs_derived.png` | 真値・平滑化真値・Q-state・旧方式の比較 |
| `img/fig09_qfull_sensor_count_sweep.png` | センサ数スイープ（q と温度） |

## 実行

必要なPythonパッケージは `numpy`、`scipy`、`matplotlib`、`Pillow` です。

```bash
cd python
python3 enkf_heatid.py    # 温度場版
python3 esmda_qfull.py    # 完全逆問題版
```

## 主な設定

| パラメータ | 値 |
|-----------|-----|
| グリッド | 30×30 (n=900) |
| 初期平均温度 | 全セル 20 degC |
| 同化手法 | 確率的ESMDA |
| アンサンブル数 | 300 |
| 初期相関スケール | 1.5 cells |
| 局所化半径 | 5 cells（温度場版）/ 10 cells（完全逆問題版） |
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

## 完全逆問題版の要点

- 温度場の復元精度は温度場版と同等（m=300 で 20.06 degC）。
- 熱源のセルスケールのエッジは点温度観測から原理的に復元できず、
  q RMSE は 38.8 で飽和（q_true の std は 40.0）。
- ただし σ=1 で平滑化したスケールでは相関 0.76 まで復元でき、
  目・鼻・輪郭が視認できる熱源場が直接得られる。
- 局所化半径を温度場版と同じ 5 にすると反復途中に温度 RMSE が暴走する
  過渡が出る（順モデルの非局所性との干渉）。10 に広げると解消。
- 詳細は [`esmda_qfull_method.md`](esmda_qfull_method.md) と
  [`../../WORKLOG_20260702_esmda_qfull.md`](../../WORKLOG_20260702_esmda_qfull.md)。

## 関連デモ

- [`A02_miffy_enkf_SCrandomTimefix_T`](../../A02_miffy_enkf_SCrandomTimefix_T/) — EnKF固定センサ数比較
- [`A02_miffy_enkf_SCrandomTime_T`](../../A02_miffy_enkf_SCrandomTime_T/docs/enkf_method.md) — 基本EnKF
- [`A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q`](../../A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q/docs/README.md) — 決定論的ESMDA/OI版
