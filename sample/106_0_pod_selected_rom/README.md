# 106_0 POD選定ROM × データ同化

研究の出発点は104のOpenFOAM＋FrontISTRによる実ソルバ同化。
600秒の現象・5メンバー・10回更新で、温度場RMSEは7.609→0.01064 K、Qは14.926 W（真値15 W）となった。
一方、当時の作業記録では直列実行に約13時間を要したため、ROMで多数の観測配置・乱数条件を比較する。
計算条件・所要時間の出典と限界は [本編冒頭](docs/02_full_story.md) に記載。

105までのROMは「勘で置いた測定プローブ」を代表点に流用していた。
106では **POD（固有直交分解）で温度場の型を抽出し、Q-DEIM で代表点をデータから系統的に選ぶ**。
その5点でROMを組んで校正し、データ同化で温度・発熱量Q・放熱hを推定、全温度場を復元する。

## ストーリー（絵・数式・結果）

- **`docs/01_pod_qdeim_algorithm.md`** … POD/Q-DEIMの数式展開とプログラム手続き（学習用）
- **`docs/02_full_story.md`** … 全体を絵・数式・結果グラフで追う本編

## 流れと再現

```
OMP_NUM_THREADS=4 python3 run/select_points_qdeim.py    # ① POD+Q-DEIMで代表5点を選定
OMP_NUM_THREADS=4 python3 run/build_calibrate_rom.py    # ② 5点ROMを102_0に校正(残差0.014K)
OMP_NUM_THREADS=4 python3 run/make_pod_modes_fig.py     # PODモードの絵
OMP_NUM_THREADS=4 python3 run/run_da.py                 # ③ 拡大状態EnKFで同化(T,Q,h)＋温度追従図
OMP_NUM_THREADS=4 python3 run/run_da_compare.py         # ④ 観測構成の比較(有用性, 5seed)
OMP_NUM_THREADS=4 python3 run/run_sensor_and_disp.py    # ⑤ 熱感度センサ選定/per-node比較/変位追従の検証
OMP_NUM_THREADS=4 python3 run/run_make_gifs.py          # ⑥ 温度・変形GIF(高感度/低感度/真値)
```

観測構成の比較は **アンサンブル N=60 メンバー、5seed 平均**。RMSE は全5点温度の
$\sqrt{\frac15\sum_i(\hat T_i-T_i^{true})^2}$ を過渡期(0-300s)平均したもの。

## 主な結果

- 温度場は **平均場＋2モードで99.9%** 説明できる → 少数点で表せる
- Q-DEIMが選んだ5点は流用点と違い**底面も含めて場を広くカバー**
- POD選定5点ROMは **OpenFOAMを残差0.014Kで再現**（ΣC=1211 J/K で物理的にも妥当）
- **有用性（過渡期の平均5点温度RMSE, 5seed平均）**:

| 観測構成 | 温度RMSE | 変位差RMSE |
|---|---|---|
| 同化なし(free run) | 4.585 K | 0.276 µm |
| 温度1点（熱感度 $W$ 低：P3 底面） | 0.534 K | 0.794 µm |
| 温度1点（熱感度 $W$ 高：P1 中央高さ） | 0.841 K | 1.132 µm |
| 温度2点（熱感度 $W$ 上位2：P1+P2） | 0.272 K | 0.335 µm |
| 温度2点+変位2点（上面A/O） | **0.147 K** | **0.162 µm** |

高低の基準は **熱感度 $W=K_s^{-1}H_T$**（その点の温度が1 K上がると注目変位が何µm動くか）であって、
発熱への温度感度 $\partial T/\partial Q$ ではない。両者の相関はわずか +0.045 で順位が一致しない
（比較は `run/sensitivity_criterion_compare.py`、解説は `docs/blog_005_optimal_sensor_placement.md` §4-2b/§4-2c）。

→ ①データ同化は必須、②**温度1点では $W$ の高低で逆転もする＝1点では足りない**、③変位観測が決定的。

### 5点から復元した温度分布は、測っていないセルでも合うか

`run/recon_unobserved_check.py` で、観測に使った5セルを**除外して**評価した。

| 評価対象 | 復元RMSE |
|---|---|
| 観測点（5セル） | 0 K（定義上ぴったり一致） |
| **未観測 20,691セル（訓練内）** | **0.00081 K**（最大 0.0168 K＝温度変動幅5.78 Kの0.29 %） |
| **未観測 20,691セル（訓練外）** | **0.00234 K**（前半60枚でPOD・選点→後半61枚で評価） |

観測ノイズ 0.3 K より2〜3桁小さい。図は `docs/img/recon_unobserved.png`。

## 発表資料

- **[オープンCAE学会シンポジウム B-16 発表スライド（20枚, reveal.js）](https://kamakiri1225.github.io/thermal-da-demo/opencae2026.html)**
  … ローカルでは `presentation/opencae_2026_slides.html` を直接ブラウザで開く
- 原稿・調査メモ: `presentation/opencae_2026_slide02_intro.md`（導入の文献調査・メーカー比較）、
  `presentation/opencae_2026_qa.md`（想定質問Q1〜Q20）
- 学会ポスター（4学会分）: `presentation/conference_posters.html`

## 前提

- 隣に `102_0`（OpenFOAM CHT）, `102_1`（FrontISTR熱膨張連成コード）
- `dacore/rom_general.py` … 任意配置ノードの一般化ROM（105の固定トポロジを一般化）
- WSL2では `OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4` 必須
