# 104_0 OpenFOAM × FrontISTR × EnKF ― 温度＋変位観測のデータ同化

103 は温度だけの観測だったが、**104 は FrontISTR 熱膨張解析の変位も観測に加える**。
実験(熱電対＋変位センサ)と同じ観測構成で、でたらめな初期状態から固体温度場と
ヒータ発熱 Q を補正する双子実験(OSSE)。

## 📊 発表スライド（ブラウザで閲覧）

**[reveal.js プレゼンをブラウザで開く（GitHub Pages）](https://kamakiri1225.github.io/thermal-da-demo/sample/104_0_openfoam_frontistr_da_enkf/presentation/conference_reveal.html)**
― 章＝横スクロール(←→)、章内＝縦スクロール(↑↓)、数式はTeX、温度/流速/変形のGIFアニメつき。
配布用は `presentation/conference_slides.pdf`。構成: 第1部 OpenFOAM＋FrontISTR →
第2部 ROM → 全体まとめ。

## 解説書（`docs/`, データ同化を理解するために）

- `09_what_104_does.md` … まず全体像 / `02_beginner_guide.md` … 初学者向け超ていねい
- `08_what_was_computed.md` … 何をどこで計算したか（作業の流れ・成果物の由来）
- `05_folder_structure.md` … フォルダ地図 / `07_cycle_and_timing_faq.md` … サイクル・分散のFAQ
- `00_temp_displacement_da.md`・`10_q_estimation_explained.md` … 手法の数式（Q推定の全記号定義）
- `03_rom_derivation.md`・`11_rom_concept.md` … ROMの導出・発想・変位同化の成否
- `01_method_walkthrough.md`・`06_ensemble_calculation.md` … 実装の流れ・アンサンブル計算
- `04_paraview.md` … ParaViewでの分布・時刻歴の見方

## 観測構成

| 観測 | 由来 | ノイズ |
|------|------|--------|
| T_hot / T_cold(温度2点) | 固体表面近傍セル | 0.3 K |
| uz_heater / uz_opposite(上面変位2点) | **メンバーの温度場→FrontISTR線形静解析→上面Uz** | 0.1 µm |

変位の観測演算子は行列ではなく **物理チェーンそのもの**
(OpenFOAM 固体温度場 → IDW補間 → FrontISTR → 上面変位)。
EnKF はメンバーごとに FrontISTR を回して予報観測 Yf をサンプリングし、
標本共分散で更新する(非線形観測演算子の標準的な扱い)。

## 状態・モデル(103 と同じ)

- 状態: 固体温度場 全20,696セル ＋ ヒータ発熱 Q(拡大状態)
- 前進モデル: chtMultiRegionFoam(輻射fvDOM、102_0の既存メッシュ、直列)
- 更新して書き戻すのは温度場と Q のみ(安定)。変位は観測にだけ使う

## 実行方法

前提: `../102_0`(メッシュ)と `../102_1`(FrontISTR連成コード)が隣にあり、
`fistr1` が PATH にあること。

```
pip install --user -r requirements.txt
cd openfoam && sh setup_base_case.sh && cd ..
OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
  nohup python3 run/run_openfoam_fem_enkf.py > openfoam/run_fem_enkf.log 2>&1 &
tail -f openfoam/run_fem_enkf.log
```

出力:
- `docs/img/openfoam_fem_enkf_rmse.png` … 温度場 RMSE の収束
- `docs/img/openfoam_fem_enkf_Q.png` … Q 推定
- `docs/img/openfoam_fem_enkf_disp.png` … 変位観測の予報 vs 真値
- `docs/img/openfoam_fem_enkf_field.png` … 温度場スナップショット
- `results/openfoam_fem_enkf_summary.yaml` / `_history.csv`

## 構成

```
dacore/enkf.py     EnKF(Yf サンプリング対応版)
daof/of_fem_twin.py  双子実験ドライバ(温度＋変位観測)
fem/fem_obs.py     FrontISTR 変位観測演算子(102_1 の連成コードを再利用)
openfoam/          base_case 生成・設定(重い生成物は Git 管理外)
```

理論(EnKF の数式・実装対応)は `../103_0_openfoam_da_enkf/docs/05_theory_to_code.md`、
FrontISTR 連成の仕組みは `../102_1_frontistr_hollow_cylinder_thermal_expansion/docs/` を参照。
