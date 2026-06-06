# ドキュメント目次

`001_kalman_thermal_1d` — 線形カルマンフィルタによる1次元熱データ同化デモ

---

## ドキュメント一覧

| # | ファイル | 種別 | 内容概要 |
|---|---|---|---|
| 01 | [01_overview.md](01_overview.md) | 概要 | プロジェクト全体の目的・フォルダ構成・実行方法 |
| 02 | [02_theory.md](02_theory.md) | 理論 | 状態空間モデルの導出・KF の理論的位置づけ |
| 03 | [03_theory_to_code.md](03_theory_to_code.md) | 理論↔実装対応 | 数式がコードのどこに対応するかを1対1で示す |
| 04 | [04_discretization_guide.md](04_discretization_guide.md) | 離散化解説 | 連続時間 ODE を前進オイラー法で離散化する手順（**本コードで使用**）。ZOH は未使用の代替として比較表に記載 |
| 05 | [05_program_guide.md](05_program_guide.md) | 実装解説 | 工作機械の問題設定からプログラムの全体フローを解説 |
| 06 | [06_results_guide.md](06_results_guide.md) | 結果解説 | 出力グラフ・RMSE の見方と各ケースの意味 |
| 07 | [07_nextstep.md](07_nextstep.md) | 研究展望 | OpenFOAM + OI への橋渡し（002 サンプルへの接続） |
| 08 | [08_blog.md](08_blog.md) | ブログ | 全体を物語形式で解説した読み物 |

---

## 各ドキュメントの詳細

---

### 01 — プロジェクト概要

**ファイル**: [`01_overview.md`](01_overview.md)
**種別**: 概要

**何が書いてあるか**

このサンプルの目的・フォルダ構成・実行方法を1枚にまとめたドキュメント。
Python だけで完結するデモ（OpenFOAM 不要）。

**主な内容**

- 10 節点モデル・センサ配置（Node 0 と Node 9）の概要
- フォルダ構成（`main.py`, `kalman_filter.py`, `thermal_model.py`）
- 実行コマンドと出力ファイル一覧

**読むべき人**: 全員（最初に読む）

---

### 02 — 理論解説

**ファイル**: [`02_theory.md`](02_theory.md)
**種別**: 理論

**何が書いてあるか**

熱伝導 PDE から離散時間状態空間モデルを導出し、カルマンフィルタの理論的位置づけを解説する。
最も詳細な理論ドキュメント。

**主な内容**

- PDE → 有限差分 → 連続時間 ODE（$\dot{\boldsymbol{\theta}} = \mathbf{A}_c \boldsymbol{\theta} + \mathbf{B}_c u$）
- 境界条件（対流境界）の行列への組み込み方
- 前進オイラー法による離散化（$\mathbf{A}_d$, $\mathbf{B}_d$ の導出）
- 安定条件（$r \Delta t < 0.5$）の確認
- カルマンフィルタのベイズ推定・最小分散推定（MMSE）としての解釈
- カルマンゲインの物理的意味
- KF と OI・EnKF・4DVar との位置関係

**読むべき人**: 理論を深く理解したい人

---

### 03 — 理論↔コード対応ガイド

**ファイル**: [`03_theory_to_code.md`](03_theory_to_code.md)
**種別**: 理論↔実装対応

**何が書いてあるか**

`02_theory.md` の数式が `thermal_model.py`・`kalman_filter.py`・`main.py` のどの行に対応するかを1対1で示す早見表。

**主な内容**

- 物理（熱拡散 PDE）→ 連続時間 ODE → 離散時間方程式 → KF → 熱変位の全フロー図
- 各数式とコード行の対応表
- $\mathbf{A}_c$, $\mathbf{B}_c$, $\mathbf{A}_d$, $\mathbf{B}_d$ の実装箇所
- 予測ステップ・更新ステップのコード対応

**読むべき人**: コードと理論を照合したい人

---

### 04 — 離散化ガイド（前進オイラー法）

**ファイル**: [`04_discretization_guide.md`](04_discretization_guide.md)
**種別**: 離散化解説

> **注意**: 旧ファイル名は `zoh_guide.md` だったが、内容は ZOH ではなく**前進オイラー法**（コードで実際に使用）の解説。
> ZOH は「代替手法の比較表」として Section 6 に1箇所出てくるだけで、本プログラムでは使用していない。

**何が書いてあるか**

連続時間の微分方程式（$\dot{\boldsymbol{\theta}} = \mathbf{A}_c \boldsymbol{\theta} + \mathbf{B}_c u$）を、コンピュータで扱える差分方程式に変換する方法を解説する。

**主な内容**

- 「離散化」の概念（連続時間 vs 離散時間）
- 前進オイラー法の導出（前進差分近似から）
- $\mathbf{A}_d = \mathbf{I} + \mathbf{A}_c \Delta t$、$\mathbf{B}_d = \mathbf{B}_c \Delta t$ の計算（**thermal_model.py で実装済み**）
- 安定条件（$r \Delta t = 0.012 \ll 0.5$）の確認
- 前進オイラー / 後退オイラー / ZOH の比較表（ZOH はここのみ登場・未使用）

**読むべき人**: 離散化の手順を確認したい人・数値安定性を理解したい人

---

### 05 — プログラム解説

**ファイル**: [`05_program_guide.md`](05_program_guide.md)
**種別**: 実装解説

**何が書いてあるか**

「工作機械の熱変位補償」という現実の問題設定から出発し、Python プログラム全体の構成と動作を解説する。

**主な内容**

- 問題設定（工作機械スピンドルの熱変位・少数センサの制約）
- `thermal_model.py`: 1次元熱伝導棒モデルの実装
- `kalman_filter.py`: 線形離散時間 KF の実装（予測・更新ステップ）
- `main.py`: メインループ（DA なし / DA あり の比較）
- プロセスノイズ $\mathbf{Q}$ と観測ノイズ $\mathbf{R}$ の設定指針
- センサ位置・初期値・ステップ数のパラメータ説明

**読むべき人**: コードを変更・拡張したい人

---

### 06 — 結果の見方

**ファイル**: [`06_results_guide.md`](06_results_guide.md)
**種別**: 結果解説

**何が書いてあるか**

`main.py`, `main_compare_da.py`, `main_disp_obs.py` が出力する 3 つのグラフの読み方を解説する。

**主な内容**

- `results.png`: 10 節点すべての真値 vs KF 推定値
- `results_compare_da.png`: DA なし（モデルのみ）vs DA あり（KF）の RMSE 比較
- `results_disp_obs.png`: 変位観測を追加した場合の改善効果
- RMSE の計算区間（定常域後半）と数値の読み方
- センサ節点と非センサ節点での挙動の違い

**読むべき人**: 実行後に結果を確認したい人

---

### 07 — 次のステップ（OpenFOAM 連成へ）

**ファイル**: [`07_nextstep.md`](07_nextstep.md)
**種別**: 研究展望

**何が書いてあるか**

Python 版（001）から OpenFOAM 版（002）への橋渡しとなるドキュメント。
`laplacianFoam` をソース改造なしで使うデータ同化の方針を提案する。

**主な内容**

- Python 版 KF（001）の限界（$\mathbf{A}_d$ が Python 近似）
- OpenFOAM をブラックボックスとして使う考え方
- OI への切り替え方針（$\mathbf{A}_d$ 不要 → 固定 $\mathbf{B}$ に置き換え）
- 002 サンプルへの接続（`laplacianFoam + OI`）

**読むべき人**: 001 から 002 への移行を理解したい人

---

### 08 — ブログ形式解説

**ファイル**: [`08_blog.md`](08_blog.md)
**種別**: ブログ

**何が書いてあるか**

「2つの温度センサから熱棒全体の温度を推定する」というテーマで、手法・結果・意義を物語形式でまとめた読み物。理論の数式より「なぜ・どのように」を中心に書いている。

**主な内容**

- 問題の設定（少数センサ・未観測節点の推定）
- カルマンフィルタが「モデル予測」と「観測」を混ぜる仕組みの直感的説明
- カルマンゲインの意味（$\mathbf{K}$ が小さい = モデルを信頼 / 大きい = 観測を信頼）
- 「真値」「モデルのみ」「KF あり」の比較グラフの見方
- 工作機械熱変位補償への応用イメージ
- Python コードとの対応

**読むべき人**: 全体像を直感的に理解したい人・外部共有用

---

## 推奨読書順

```
はじめて読む場合:
  01_overview.md  →  08_blog.md  →  06_results_guide.md

理論を深く理解したい場合:
  02_theory.md  →  04_zoh_guide.md  →  03_theory_to_code.md

コードを変更・拡張したい場合:
  05_program_guide.md  →  03_theory_to_code.md

002 サンプルへ進む場合:
  07_nextstep.md  →  002-0_laplacian_da_1d/docs/00_index.md
```

---

## フォルダ構成（参考）

```
001_kalman_thermal_1d/
├── thermal_model.py       ← 1次元熱伝導棒モデル（A_c, A_d, B_c, B_d）
├── kalman_filter.py       ← 線形離散時間カルマンフィルタ
├── main.py                ← 基本デモ（温度推定）
├── main_compare_da.py     ← DA なし vs DA あり の比較
├── main_disp_obs.py       ← 変位観測追加版
├── results*.png           ← 出力グラフ
└── docs/
    ├── 00_index.md        ← 目次（このファイル）
    ├── 01_overview.md     ← 概要・実行方法
    ├── 02_theory.md       ← 状態空間モデル導出・KF 理論
    ├── 03_theory_to_code.md ← 理論↔コード対応
    ├── 04_discretization_guide.md ← 前進オイラー離散化（使用中）・ZOH比較（未使用）
    ├── 05_program_guide.md ← プログラム解説
    ├── 06_results_guide.md ← 結果の見方
    ├── 07_nextstep.md     ← 002 サンプルへの橋渡し
    ├── 08_blog.md         ← ブログ形式解説
    └── pdf/               ← PDF 出力（自動生成）
```
