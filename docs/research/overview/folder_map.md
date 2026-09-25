# 本研究で使用したフォルダ一覧

熱データ同化研究（OpenFOAM＋FrontISTR＋データ同化／ROM）で実際に使ったフォルダの地図です。
リポジトリのルートは `20260505_datadoka/`、計算ケースはすべて `sample/` 配下にあります。

作成日: 2026-09-25

---

## 1. 研究の本体：中空円筒のケース群（101〜106）

番号は研究の進行順です。**下に行くほど新しく、上の成果を材料にしています**。

| フォルダ | 容量 | 役割 | この研究での位置づけ |
|---|---|---|---|
| `sample/101_0_openfoam_cht_radiation_box` | 1.0 GB | 箱形状のCHT＋輻射＋ヒートマット発熱 | **予備検討**。連成解析の練習台 |
| `sample/101_1_frontistr_cht_box_thermal_expansion` | 35 MB | 101_0の温度場でFrontISTR熱膨張 | 予備検討。温度→変位の受け渡し手順を確立 |
| `sample/102_0_openfoam_hollow_cylinder_heat_transfer` | 1.1 GB | **中空円筒の過渡CHT解析** | **本研究の土台**。ヒータ15 W・0〜300 s加熱の温度場を生成（blog_001） |
| `sample/102_1_frontistr_hollow_cylinder_thermal_expansion` | 635 MB | **中空円筒のFrontISTR熱膨張** | 温度→変位の写像。`cylinder_mesh.py` は他ケースからも参照（blog_001） |
| `sample/103_0_openfoam_da_enkf` | 478 MB | OpenFOAM × EnKF | 温度のみのデータ同化。EnKF実装の検証 |
| `sample/103_1_openfoam_da_pf` | 291 MB | OpenFOAM × 粒子フィルタ | 手法比較（EnKFとの対比） |
| `sample/104_0_openfoam_frontistr_da_enkf` | 4.3 GB | **実ソルバでの温度＋変位同化** | **最大のケース**。OI/EnKFを実ソルバで回した比較の出どころ（blog_002 §6、ポスターP1） |
| `sample/105_0_sensor_placement_sensitivity` | 148 MB | **センサ配置と熱感度** | 感度行列 $W=K_s^{-1}H_T$ で観測点を選ぶ。高W/低W節点の根拠（blog_002 §4、blog_005） |
| `sample/106_0_pod_selected_rom` | 128 MB | **POD選定ROM × データ同化（主フォルダ）** | **現在の作業場所**。POD＋Q-DEIMで代表点を選び、軽量ROMで同化。ブログ・ポスターもここ |

### 初期の1次元・丸棒ケース（手法の学習段階）

| フォルダ | 容量 | 内容 |
|---|---|---|
| `sample/001_python_kalman_thermal_1d` | 30 MB | Pythonのみの1次元カルマンフィルタ |
| `sample/002_0_openfoam_laplacian_da_1d` | 39 MB | laplacianFoam＋1次元同化 |
| `sample/002_1_openfoam_round_bar_da` | 76 MB | 丸棒のOI同化（OpenFOAM） |
| `sample/002_2_frontistr_round_bar_da` | 51 MB | 丸棒のFrontISTR版 |
| `sample/002_3_calculix_round_bar_da` | 8 KB | 丸棒のCalculiX版（スクリプトのみ） |

---

## 2. 主フォルダ `sample/106_0_pod_selected_rom/` の内部

現在の研究・ブログ・学会ポスターはすべてここで作っています。

| サブフォルダ | 容量 | 中身 |
|---|---|---|
| `dacore/` | 132 KB | **共通ライブラリ**。`rom_general.py`（集中定数ROMとRK4積分）、`enkf.py`（EnKF更新）、`plots.py`（日本語フォント設定） |
| `run/` | 428 KB | **実行スクリプト42本**。ROM構築・校正、Q-DEIM選点、感度解析、EnKF比較、各種作図 |
| `results/` | 1.7 MB | 計算結果の `.npz`（校正済みROM、POD モード、変位オペレータ、感度）と集計表 |
| `docs/` | 36 MB | **ブログ本文5本**（`blog_001`〜`blog_005`）、`img/`（図・アニメ）、`pdf/`（二段組TeX版PDFとビルド系） |
| `presentation/` | 3.3 MB | 学会発表ポスター（`conference_posters.html` / `.pdf`、4学会分） |
| `openfoam/` | 87 MB | 実ソルバ連携用の作業ディレクトリ（大きいため非公開） |
| `config/`, `fem/` | 28 KB | 解析条件とFEM設定 |

### 他ケースへの依存（`run/` のコードが参照している先）

- `102_1_...` … メッシュ生成 `cylinder_mesh.py`（20か所で参照。最多）
- `105_0_...` … 熱感度 `kinvh_sensitivity.npz`（高W/低W節点の座標）
- `104_0_...` … 実ソルバの真値・変位場（`displacement_fields.vtu` など）
- `102_0_...` … CHTスナップショット（PODの学習データ）

---

## 3. 文書・公開まわり

| 場所 | 内容 |
|---|---|
| `docs/research/overview/` | 研究ロードマップ・サマリ（`research_roadmap.md`, `research_summary.md`）、本ファイル |
| `sample/106_0_pod_selected_rom/docs/blog_00N_*.md` | ブログ本文（原本。ここを直すとPDF・Webに反映） |
| `sample/106_0_pod_selected_rom/docs/pdf/` | 論文形式PDF（A4二段組）とビルドスクリプト `build_pdf.py` |
| `sample/105_0_sensor_placement_sensitivity/docs/` | センサ配置の実務ガイド（`03_sensor_design_guide.md` など） |
| gh-pages ブランチ | 公開サイト <https://kamakiri1225.github.io/thermal-da-demo/> （ブログHTML・ポスター・PDF） |

---

## 4. 研究の流れとフォルダの対応

```
102_0（CHTで温度場）
   ↓ 温度
102_1（FrontISTRで熱膨張）          ← blog_001：土台づくり
   ↓ 温度場・変位
103_0 / 103_1（EnKF・PFの検証）
   ↓ 手法
104_0（実ソルバで温度＋変位を同化）  ← blog_002 §6・blog_003：重いが本物
   ↓ 「実ソルバでは遅い」という課題
105_0（感度Wで観測点を選ぶ）        ← blog_002 §4・blog_005：どこに置くか
   ↓ 観測設計
106_0（PODで代表点→軽量ROMで同化）  ← blog_004・blog_005：現在地
```

---

## 5. 補足：この一覧に入れていないもの

- `IoT/` … 熱膨張実験装置（ESP32・ダイヤルゲージ）。CAE解析とは別系統
- `program/`, `study/` … 学習用のコードと練習。研究本体では使用していない
- `docs/career_building/` … 研究以外の個人文書

容量の大きいケース（`101_0`, `102_*`, `104_0` など）は計算結果を含むため、
GitHubには `.gitignore` で除外し、スクリプトと図・文書のみを公開しています。
