# オープンCAE 2026 発表スライド案 ― 導入「はじめに」

作成日: 2026-09-29 ／ 文献調査に基づく導入（1枚目）の草案

---

## 本番用テキスト（そのままスライドに貼れる形）

### 課題（箇条書き3行）

- 工作機械の加工誤差のうち、**熱変位が占める割合は最大75 %** であり、最大の誤差要因である。
- 対策は「**熱的に素性の良い設計**」「**温度制御**」「**モデルによる補償**」の3系統だが、
  実用の主役は補償であり、その前提は**熱変形を正しく予測できること**である。
- しかし予測が難しい。**TCP変位は温度の絶対値ではなく温度勾配で決まる**ため、
  **温度の傾向が合っても変形の傾向が合うとは限らない**（時定数が複数あると変位の向きすら途中で変わる）。
  さらに**熱伝達率などの境界条件は直接測れず**、その不確かさがTCP予測精度を直接左右する。

### ゆえに求められること（1行）

> **限られた実測から境界条件ごと補正し、熱変形そのものを推定する枠組みが求められる。**

---

## 各行の根拠（引用元）

### 1行目：最大75 %（2026-09-29 訂正）

**訂正**：当初「40〜70 %［Bryan 1990］」としましたが、Bryan (1990) は有料で**原典を確認できません**でした。
代わりに**全文無料で確認できた「最大75 %」**に差し替えます。

- **Li, Z., Vogl, G. W., Kinzel, E. C., Santa, B., Landers, R. G. (2024)**
  "Machine Tool Thermal Error Measurement and Prediction via Wireless Microscope",
  *Manufacturing Letters* 41, 1440–1451（オープンアクセス）
  📄 <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076>
  > "Thermal errors can contribute **up to 75 percent** of the overall machining errors of a machined part."
- **Bünger, A. et al. (2023)** arXiv:2306.12736（全文無料）
  📄 <https://arxiv.org/abs/2306.12736>
  > "According to **Mayr et al., 2012, the thermal error accounts for 75 %** of the total manufacturing error."

→ どちらも **Mayr et al. (2012)** を典拠としており、**孫引きの連鎖が確認できる**。

**スライドでの書き方の例**：「熱変位は加工誤差の**最大75 %**を占める［Mayr et al. 2012］」

### 2行目：対策の3系統

- **Mayr, J. et al. (2012)** "Thermal issues in machine tools",
  *CIRP Annals* 61(2), 771–791.（CIRPキーノート、この分野の標準レビュー）
  → 研究領域を「**測定**」「**計算（モデリング）**」「**熱影響の低減・補償**」の3つに整理している。
  設計による低減（thermally symmetric design など）と、モデルによる補償の両方を扱う。

**補足**：設計で消せるのは一部で、残りを補償で取るのが実務の構図。
これを1行目と3行目の橋渡しに使う。

### 3行目：温度が合っても変形は合わない

これが**本発表の動機そのもの**なので、根拠を3つ重ねると強い:

1. **勾配が効く（てこの原理）**
   柱の上下にわずかな温度勾配があるだけで曲げ変形が生じ、
   主軸先端・TCPでは**てこで拡大されて大きな変位**になる。
   → だから「平均温度が合っている」ことは何の保証にもならない。

2. **変位の向きが途中で変わる**
   Mayr et al. (2012) では、**異なる時定数が2つ関与する場合や、熱源からの急な勾配が
   時間とともに均一化していく場合に、TCP変位が運転中に向きを変えうる**ことが指摘されている。
   → 温度の単調な上昇トレンドが合っていても、変形の符号を外す可能性がある。

3. **境界条件が測れない**
   **Bünger, A., Herzog, R., Naumann, A., Stoll, M. (2023)**
   "Uncertainty Propagation of Initial Conditions in Thermal Models", arXiv:2306.12736
   → 「**the accuracy of the TCP prediction depends highly on the accuracy of the model parameters,
   such as heat exchange parameters, and the initial temperature**」と明記。
   熱交換パラメータ（熱伝達率）と初期温度の不確かさが、TCP予測精度を直接左右する。
   さらに **CHTC（対流熱伝達率）を直接測る計測器が存在しない**ことも複数の論文で指摘されている。

**スライドでの書き方の例**：
「TCP変位は温度勾配で決まる。温度の傾向が合っても変形の傾向が合う保証はない［Mayr 2012］。
しかも熱伝達率は直接測れない［Bünger 2023］」

---

## 最も近い先行研究（質疑への備え）

**Teshima, Y., Tanaka, S., Kizaki, T., Sugita, N. (2024)**
"Sensor placement strategy based on reduced-order models for thermal error estimation in machine tools",
*CIRP JMST* 55, 403–410. DOI: 10.1016/j.cirpj.2024.10.015

- **やっていること**：温度入力とTCP変位を結ぶ**縮約モデル（ROM）**を作り、
  回帰係数を各測定点の**温度感度**とみなして、**センサの最適な数と配置の指針**を与える。
- **彼らが挙げるギャップ**：「入力点を増やせば精度は上がるが、
  **最適な点数と配置を決める方法がこれまで無かった**」
- 続報（2025）もある：
  "Strategy for Sensor Placement to Estimate Thermal Errors Using Temperature-Sensitivity
  Distribution Based on a Reduced-Order Model of Machine Tools"

**本研究との違い（言えること）**

| | Teshima et al. 2024 | 本研究 |
|---|---|---|
| ROMの作り方 | 温度→TCP変位の**回帰** | **POD＋Q-DEIM**で代表点を選び、集中定数ROMを校正 |
| 配置の基準 | 温度感度（回帰係数） | **観測の価値** $\Delta=\mathrm{Cov}(X,y)^2/(\mathrm{Var}(y)+r)$ 、推定対象ごとに変わる |
| 未知パラメータ | 扱わない（温度は入力） | **発熱量 $Q$ ・放熱 $h$ を状態に入れて同時推定**（EnKF） |
| 実装 | ― | **OSSのみ**（OpenFOAM＋FrontISTR＋Python）で公開 |

→ **「センサ配置をROMで決める」という流れは既にあるが、
データ同化で境界条件ごと推定する点と、OSS一気通貫で公開する点が本研究の位置づけ。**

---

## 使い方のメモ

- オープンCAEの聴衆向けなので、**1行目で数字を出して「大きい問題だ」と示し**、
  **3行目で「だからシミュレーションだけでは足りない」へ落とす**のが効く。
- 図を入れるなら、1行目の横に工作機械の熱変形の模式図、
  3行目の横に「温度は合っているのに変位がずれている」グラフがあると理想。
  → 本研究の `oi_fullsolver_compare.png`（温度RMSEは0.09 Kまで下がるのに $Q$ は8.9 W止まり）は
  **まさにこの構図**なので、導入の伏線として使える。

## 詳しい実務調査

熱変位補正が実際に何をしているか、メーカー各社（オークマ・マザック・DMG MORI・FANUC）が
どう実装しているかは、**`thermal_compensation_survey.md`** に無料リンク付きでまとめています。

## 確認しておくとよい点

- Bryan (1990) と Mayr et al. (2012) は**有料論文**のため、本メモの記述は
  アブストラクト・レビュー論文での引用・検索結果に基づく。
  **スライドに数字を載せる前に、原典または入手可能なレビューで一度確認**しておくと安全。
- 「40〜70 %」は Bryan (1990) 由来、「up to 70 %」は Teshima et al. (2024) の冒頭。
  **どちらを引くかで表現を揃える**こと。
