# ドキュメント目次

`002-0_laplacian_da_1d` — laplacianFoam + Optimal Interpolation データ同化サンプル

---

## ドキュメント一覧

| # | ファイル | 種別 | 内容概要 |
|---|---|---|---|
| 01 | [01_overview.md](01_overview.md) | 概要 | プロジェクト全体の目的・設定・実行方法 |
| 02 | [02_theory.md](02_theory.md) | 理論 | OI の理論式・ゲインの導出・OpenFOAM との接続 |
| 03 | [03_program_guide.md](03_program_guide.md) | 実装解説 | コードと理論式の対応・実行フロー |
| 04 | [04_results_guide.md](04_results_guide.md) | 結果解説 | グラフの読み方・CSVの構造・RMSEの見方 |
| 05 | [05_taxonomy_and_refs.md](05_taxonomy_and_refs.md) | 分類・参考文献 | DA 手法の全体マップ・位置づけ・参考文献リスト |
| 06 | [06_nextstep.md](06_nextstep.md) | 研究展望 | 実測データへの展開・次のアクション |
| 07 | [07_blog.md](07_blog.md) | ブログ | 手法・結果・応用先・研究の問いをブログ形式で整理 |

PDF 資料（参考文献）：
- `MTI_Handbook_Nakano_Assimilation_Ver_1_2.pdf` — 中野 慎也ほか「データ同化入門」ハンドブック版
- `20170827_ito_06.pdf` — データ同化関連講義資料

---

## 各ドキュメントの詳細

---

### 01 — プロジェクト概要

**ファイル**: [`01_overview.md`](01_overview.md)
**種別**: 概要

**何が書いてあるか**

このサンプルの目的・何を行ったか・現在の設定値・実行方法を1枚でまとめたドキュメント。
はじめて読む人は、まずここから読む。

**主な内容**

- ツイン実験の構成（G_TRUE / G_MODEL / センサ配置）
- シミュレーション設定（セル数・時間刻み・ステップ数・各パラメータ）
- 実行コマンド（`./run.sh`）
- 他ドキュメントへのリンク

**読むべき人**: 全員（最初に読む）

---

### 02 — 理論説明

**ファイル**: [`02_theory.md`](02_theory.md)
**種別**: 理論

**何が書いてあるか**

Optimal Interpolation (OI) の理論を詳細に解説する。

OI ゲインがなぜ

$$
\mathbf{K} = \mathbf{B}\mathbf{H}^T \left(\mathbf{H}\mathbf{B}\mathbf{H}^T + \mathbf{R}\right)^{-1}
$$

という形になるのかを、解析誤差分散の最小化から導出する。

**主な内容**

- `laplacianFoam` が解く方程式
- 状態ベクトル・観測式の定義
- OI 更新式の導出（解析誤差の最小化）
- OI ゲインの直感的な意味（センサのズレを各節点へ配る重み）
- 背景誤差共分散 $\mathbf{B}$（ガウシアン相関カーネル）と観測誤差共分散 $\mathbf{R}$ の意味
- なぜ $A_d$ が不要なのか
- OI と 3D-Var・カルマンフィルタ・EnKF との位置関係
- OI の限界（固定 $\mathbf{B}$・境界条件 $G$ を推定しない）

**読むべき人**: 理論を理解したい人・数式を確認したい人

---

### 03 — プログラム解説

**ファイル**: [`03_program_guide.md`](03_program_guide.md)
**種別**: 実装解説

**何が書いてあるか**

`da_main.py` のコードと `02_theory.md` の理論式が、どの行でどのように対応するかを示す。

**主な内容**

- 全体フロー（真値生成 → DA なし → DA あり）
- `build_oi_matrices()` の実装（$\mathbf{H}$, $\mathbf{B}$, $\mathbf{R}$）
- `run_da()` の OI 更新ステップ
- OpenFOAM 1 ステップ起動の仕組み（`write_field` → `run_step` → `read_field`）
- 計算時間・同化間隔の考え方
- $A_d$ を使わない設計方針

**読むべき人**: コードを読む人・パラメータを変更したい人

---

### 04 — 結果の見方

**ファイル**: [`04_results_guide.md`](04_results_guide.md)
**種別**: 結果解説

**何が書いてあるか**

実行後に生成されるファイル・グラフ・CSV の意味を解説する。

**主な内容**

- `results/` フォルダの構成（`truth/`, `model_only/`, `with_da/`）
- `temperature_history.csv` の列定義
- `results_of_da.png` の 4 パネルの読み方
  - (1) センサ節点の比較
  - (2) 隠れ節点の比較
  - (3) RMSE 棒グラフ
  - (4) 誤差時系列
- 最新の RMSE 結果（DA なし 1.186°C → DA あり 0.136°C、改善率 90.3%）
- OI で補正が空間的に広がる仕組み

**読むべき人**: 実行後に結果を確認したい人

---

### 05 — DA 手法の分類と参考文献

**ファイル**: [`05_taxonomy_and_refs.md`](05_taxonomy_and_refs.md)
**種別**: 分類・参考文献

**何が書いてあるか**

データ同化手法の全体マップを示し、本プログラムの OI がその中のどこに位置するかを整理する。また、理論書・原著論文・CFD + DA 論文・工作機械熱誤差補償論文を参考文献リストとして収録する。

**主な内容**

- DA 手法の全体マップ（変分法 / 逐次法）
  - 3DVar, 4DVar, OI, KF, EKF, UKF, EnKF, Particle Filter
- OI と KF の比較表
- KF と EnKF の比較（次ステップとしての EnKF の利点）
- 参考文献リスト
  - 基礎理論書（日本語・英語）
  - Kalman (1960) 原著論文
  - CFD + DA 論文（Meldi & Poux 2017 など）
  - 工作機械熱誤差補償 + DA 論文（Lang et al. 2024, 2025 など）
  - OpenFOAM + DA ツール（DAPPER, OpenDA）
- 本プロジェクトと工作機械熱問題との関係
- 研究ステップのロードマップ

**読むべき人**: 手法の位置づけを知りたい人・論文を探したい人

---

### 06 — 次のステップ

**ファイル**: [`06_nextstep.md`](06_nextstep.md)
**種別**: 研究展望

**何が書いてあるか**

現在のツイン実験を実測データに置き換えるための具体的な手順と、研究の次のアクションを示す。

**主な内容**

- 実測温度データを使うための変更箇所（`y` の生成部分）
- センサ配置・サンプリング間隔の検討
- `G` をパラメータ推定する方向への発展
- OI → EnKF への拡張ステップ
- 丸棒実験での検証計画

**読むべき人**: 実験データへの適用を検討している人

---

### 07 — ブログ形式解説

**ファイル**: [`07_blog.md`](07_blog.md)
**種別**: ブログ

**何が書いてあるか**

今回の手法・実装・結果・応用先・研究の問いを、技術ブログの形式でまとめたドキュメント。理論の数式より「なぜ・どのように」を中心に書いており、`02_theory.md` より読みやすい。

**主な内容**

- なぜ OI を使ったのか（カルマンフィルタの $A_d$ 問題を避けるため）
- $\mathbf{B}$ の直感的な意味（「一緒に間違っていそう」の地図）
- 1 センサで見た補正の直感例（重みが 0.8 なら 80% を反映）
- 解析設定の手順と OpenFOAM 1 ステップ起動コード
- 今回の結果（RMSE 90.3% 改善）
- 図の各パネルの読み方
- 手法の利点・デメリット・応用先
- 6 つの研究の問い

**読むべき人**: 全体像を直感的に理解したい人・発表資料・外部共有用

---

## 推奨読書順

```
はじめて読む場合:
  01_overview.md  →  07_blog.md  →  04_results_guide.md

理論を深く理解したい場合:
  01_overview.md  →  02_theory.md  →  03_program_guide.md

手法の位置づけを把握したい場合:
  05_taxonomy_and_refs.md

次のステップを検討したい場合:
  06_nextstep.md  →  05_taxonomy_and_refs.md（発展先の論文）
```

---

## フォルダ構成（参考）

```
002-0_laplacian_da_1d/
├── case_base/              ← OpenFOAM テンプレート（実行しない）
├── of_interface.py         ← OpenFOAM-Python ブリッジ（モデル共有）
├── docs/                   ← このフォルダ
│   ├── 00_index.md         ← 目次（このファイル）
│   ├── 01_overview.md      ← 概要・実行方法
│   ├── 02_theory.md        ← OI 理論
│   ├── 03_program_guide.md ← コード解説
│   ├── 04_results_guide.md ← 結果の見方
│   ├── 05_taxonomy_and_refs.md ← 手法分類・参考文献
│   ├── 06_nextstep.md      ← 次のステップ
│   ├── 07_blog.md          ← ブログ形式解説
│   └── *.pdf               ← 参考 PDF 資料
├── oi/                     ← Optimal Interpolation 実験
│   ├── case/               ← OI 用 OpenFOAM 作業ディレクトリ
│   ├── da_main.py
│   ├── run.sh
│   └── results/
└── kf/                     ← Kalman Filter（次ステップ）
    └── README.md
```
