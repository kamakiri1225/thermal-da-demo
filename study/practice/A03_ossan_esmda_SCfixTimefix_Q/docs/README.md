# おっさん熱源場の固定センサ ESMDA 同定（完全逆問題）

ミッフィーで確立した完全逆問題版 ESMDA
（[A02_miffy_esmda_SCfixTimefix_Q](../../A02_miffy_esmda_SCfixTimefix_Q/docs/README.md)）を、
60×60 のおっさん温度場（`A01_ossan` の PNG ピクセルアート由来）へ展開したケースです。

状態変数は熱源場 q（3600 セル）。順モデル `x = L⁻¹(−q/α)` で温度を計算し、
固定温度センサの観測だけを ESMDA で同化します。
センサ配置は真値を見ない空間充填配置、事前アンサンブルはおっさん形状を含みません。

## 実行

必要な Python パッケージは `numpy`、`scipy`、`matplotlib`、`Pillow` です。

```bash
cd python
python3 esmda_qfull_ossan.py
```

`../A01_ossan/img/001_ossan.png` を読み込むため、リポジトリ構成のまま実行してください。

## 主な設定

| パラメータ | 値 | 備考 |
|-----------|-----|------|
| グリッド | 60×60 (n=3600) | ミッフィー（30×30）の 4 倍の未知数 |
| 同化手法 | 確率的 ESMDA（状態 = 熱源場 q） | |
| アンサンブル数 | 300 | |
| 反復数 | 30 | |
| 観測ノイズ σ_r | 3 °C | |
| 事前分布の広がり | q_true の std（≈28.6）に一致 | |
| 初期相関スケール | 1.5 cells | |
| 局所化半径 R_LOC | **20 cells** | 下記参照 |
| 代表センサ数 | 100, 300, 600, 1200 | |
| センサ数スイープ | 50〜1200 点、各 3 シード | |

**局所化半径について**: ミッフィー（30×30）では R_LOC=10 で安定でしたが、
60×60 ではドメインがセル数で 2 倍になり、順モデル L⁻¹ の非局所的な感度も
セル単位でより遠くまで届きます。R_LOC=10 では反復途中の過渡暴走
（温度 RMSE が一時 1600 °C 級）が再発し、20 に広げると解消しました。
「局所化半径の下限は観測演算子の感度の物理的到達距離で決まる」という
ミッフィーの教訓の再確認になっています。

## 結果

| 指標 | 値 |
|---|---|
| 温度 RMSE（m=1200, x(q_mean)） | **17.1 °C** |
| 熱源 q RMSE | 27.7 で飽和（q_true の std は 28.6） |
| corr(q_est, σ=1 平滑化 q_true), m=100/300/600/1200 | 0.18 / 0.29 / 0.55 / **0.72** |
| 旧方式（温度同化 → q 逆算）の同相関（m=1200） | 0.71 |

ミッフィーと同じ構図です。

- セルスケールのエッジ状熱源はポイントワイズには復元不可（q RMSE ≈ 事前 std）
- 平滑化スケールでは相関 0.72 まで復元でき、輪郭・帽子・口が視認できる
- 3600 セルに対してセンサ 600 点（密度 1/6）で相関 0.55、1200 点（1/3）で 0.72。
  ミッフィー（900 セル、300 点 = 1/3 で 0.76）とセンサ**密度**でほぼ対応する

## 出力ファイル

| ファイル | 内容 |
|---------|------|
| `img/anim_esmda_qfull.gif` | 温度（上段）と熱源（下段）の同化過程 |
| `img/fig00_sensor_layout.png` | 固定センサ配置 |
| `img/fig05_qfull_rmse.png` | q / 温度 RMSE 収束（旧方式の水準線付き） |
| `img/fig06_qfull_qest_final.png` | 最終推定熱源場 |
| `img/fig07_qfull_xrecon_final.png` | 推定熱源からの順解析温度場 |
| `img/fig08_qfull_vs_derived.png` | 真値・平滑化真値・Q-state・旧方式の比較 |
| `img/fig09_qfull_sensor_count_sweep.png` | センサ数スイープ |
| `docs/sensor_count_sweep_qfull.csv` | スイープ数値データ |

## 関連

- [`A02_miffy_esmda_SCfixTimefix_Q`](../../A02_miffy_esmda_SCfixTimefix_Q/docs/README.md) — 手法確立元（ミッフィー、30×30）
- [`esmda_qfull_method.md`](../../A02_miffy_esmda_SCfixTimefix_Q/docs/esmda_qfull_method.md) — 手法の詳細解説
- [`A03_ossan_enkf_SCrandomTime_Q`](../../A03_ossan_enkf_SCrandomTime_Q/) — 同じおっさん場の EnKF 版（ランダムセンサ）
- [ESMDA 理論解説スライド（ブラウザで表示）](https://kamakiri1225.github.io/thermal-da-demo/study/practice/A02_miffy_esmda_SCfixTimefix_Q/slides.html)
