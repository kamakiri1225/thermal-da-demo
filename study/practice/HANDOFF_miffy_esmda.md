# 引き継ぎメモ: ミッフィー固定センサ ESMDA

更新日: 2026-07-02

この文書は、他のPCへそのまま引き継ぐための作業メモです。  
現時点の主対象は [A02_miffy_esmda_SCfixTimefix_Q](./A02_miffy_esmda_SCfixTimefix_Q) です。

---

## 1. 何を作っているか

- 一様な初期温度場から開始
- 固定センサの温度観測だけを使ってミッフィー温度場を復元
- 同化手法は確率的 ESMDA
- 熱量は最終温度場から `q = -alpha * L * x` で導出

決定論的 ESMDA/OI 版は別フォルダに分離しています。

- [A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q](./A02_miffy_deterministic_esmda_oi_SCfixTimefix_Q)

---

## 2. 主要ファイル

### 現行ケース

- [python/enkf_heatid.py](./A02_miffy_esmda_SCfixTimefix_Q/python/enkf_heatid.py) — 温度場版
- [python/esmda_qfull.py](./A02_miffy_esmda_SCfixTimefix_Q/python/esmda_qfull.py) — **完全逆問題版（2026-07-02 追加）**
- [docs/README.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/README.md)
- [docs/enkf_heatid_method.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/enkf_heatid_method.md)
- [docs/esmda_qfull_method.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/esmda_qfull_method.md)
- [docs/esmda_explanation.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/esmda_explanation.md)
- [docs/sensor_count_study.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/sensor_count_study.md)
- [docs/sensor_count_sweep.csv](./A02_miffy_esmda_SCfixTimefix_Q/docs/sensor_count_sweep.csv)
- [docs/sensor_count_sweep_qfull.csv](./A02_miffy_esmda_SCfixTimefix_Q/docs/sensor_count_sweep_qfull.csv)

### 生成物（温度場版）

- [img/fig00_sensor_layout.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig00_sensor_layout.png)
- [img/fig01_rmse.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig01_rmse.png)
- [img/fig02_xpred_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig02_xpred_final.png)
- [img/fig03_qest_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig03_qest_final.png)
- [img/fig04_sensor_count_sweep.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig04_sensor_count_sweep.png)
- [img/anim_enkf_heatid.gif](./A02_miffy_esmda_SCfixTimefix_Q/img/anim_enkf_heatid.gif)

### 生成物（完全逆問題版）

- [img/fig05_qfull_rmse.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig05_qfull_rmse.png)
- [img/fig06_qfull_qest_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig06_qfull_qest_final.png)
- [img/fig07_qfull_xrecon_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig07_qfull_xrecon_final.png)
- [img/fig08_qfull_vs_derived.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig08_qfull_vs_derived.png)
- [img/fig09_qfull_sensor_count_sweep.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig09_qfull_sensor_count_sweep.png)
- [img/anim_esmda_qfull.gif](./A02_miffy_esmda_SCfixTimefix_Q/img/anim_esmda_qfull.gif)

---

## 3. 再実行手順

```bash
cd practice/A02_miffy_esmda_SCfixTimefix_Q/python
python3 enkf_heatid.py
```

このスクリプトは実行開始時に `img/` 内の既知の生成物を初期化してから再生成します。  
古い画像が残って混ざることは避けています。

---

## 4. 現在の設定

| 項目 | 値 |
|---|---|
| グリッド | 30×30 = 900セル |
| 初期平均温度 | 20 degC 一様 |
| 同化手法 | 確率的 ESMDA |
| アンサンブル数 | 300 |
| 初期相関スケール | 1.5 cells |
| 局所化半径 | 5 cells |
| 反復数 | 30 |
| 代表センサ数 | 30, 100, 200, 300 |
| センサ数スイープ | 10, 20, 30, 50, 75, 100, 150, 200, 250, 300, 400 |
| スイープ回数 | 各3シード |

---

## 5. 最新の検証結果

### 代表ケース

- m=30 -> RMSE 40.50 degC
- m=100 -> RMSE 26.30 degC
- m=200 -> RMSE 22.87 degC
- m=300 -> RMSE 19.82 degC

### センサ数スイープ平均

| センサ数 | 平均RMSE [degC] |
|---:|---:|
| 10 | 50.09 |
| 20 | 45.67 |
| 30 | 40.50 |
| 50 | 30.92 |
| 75 | 27.19 |
| 100 | 26.30 |
| 150 | 23.91 |
| 200 | 22.87 |
| 250 | 20.96 |
| 300 | 19.82 |
| 400 | 19.78 |

### 実務上の目安

- 75点以上: 概形が見え始める
- 150点程度: 精度とコストのバランスが良い
- 300点程度: 高精度側
- 400点: 300点からほぼ改善なし

---

## 6. 図の読み方

- [fig00_sensor_layout.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig00_sensor_layout.png)
  - センサ配置。真値を見ずに格子座標だけで決めた空間充填配置。
- [fig01_rmse.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig01_rmse.png)
  - センサ数別のRMSE収束。
- [fig02_xpred_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig02_xpred_final.png)
  - 最終温度場。センサが増えるほど輪郭が明瞭になる。
- [fig03_qest_final.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig03_qest_final.png)
  - 最終温度場から導出した熱量分布。
  - 色は絶対値99%分位で共通スケール化している。
- [fig04_sensor_count_sweep.png](./A02_miffy_esmda_SCfixTimefix_Q/img/fig04_sensor_count_sweep.png)
  - 必要センサ数の検討図。

---

## 7. 引き継ぎ時の注意点

1. `practice/README_demo_summary.md` は全体概要として残っていますが、このケースの最新結果は古いです。
2. このケースの正式な説明は、まず [A02_miffy_esmda_SCfixTimefix_Q/docs/README.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/README.md) を見てください。
3. センサ数の根拠は [sensor_count_study.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/sensor_count_study.md) が最新です。
4. 決定論的版は別フォルダなので、手法を混ぜないでください。

---

## 8. 引き継ぎ後に完了した作業

- `practice/README.md` のESMDA節を現行設定・最新結果へ更新
- `README_demo_summary.md` を `README.md` への互換リンクに変更
- 実行依存（`numpy`, `scipy`, `matplotlib`, `Pillow`）をケースREADMEへ追記
- ケースREADMEの存在しない関連デモリンクを修正

## 9. 次の技術課題 → 完了（2026-07-02）

「状態変数を `Q_full` に変更し、`L^-1` の順問題を組み込んだ完全逆問題へ進める」
は [python/esmda_qfull.py](./A02_miffy_esmda_SCfixTimefix_Q/python/esmda_qfull.py)
として実装・検証済みです。

- 温度場は m=300 で RMSE 20.06 degC（温度場版と同等）
- 熱源はセルスケールのエッジは復元不可（q RMSE 38.8 で飽和、真値 std 40.0）だが、
  σ=1 平滑化スケールで相関 0.76、目・鼻・輪郭が視認できる
- 局所化半径は温度場版の 5 では過渡発散するため 10 に拡大

詳細な経緯と数値は [WORKLOG_20260702_esmda_qfull.md](./WORKLOG_20260702_esmda_qfull.md)、
手法は [docs/esmda_qfull_method.md](./A02_miffy_esmda_SCfixTimefix_Q/docs/esmda_qfull_method.md) を参照。

## 10. さらに次の課題案

- 非定常問題化（時系列観測で熱源の時間変化を追う）
- エッジを表現できる事前分布（非ガウス・全変動正則化系）の検討
- OpenFOAM / FrontISTR ケース（sample/002-1 系）への Q-state 同化の展開
