# オープンCAE 2026 発表スライド案 ― 導入「はじめに」

作成日: 2026-09-26 ／ 最終更新: 2026-09-29（文献調査・数字の訂正・実務調査を統合）

**方針**：主張には**無料で読める一次情報のリンク**を必ず付ける。
有料論文しか根拠が無い場合は、**そのことを明記**して使用可否を判断できるようにする。

---

## 1. 本番用テキスト（そのままスライドに貼れる形）

### 課題（箇条書き3行）

- 工作機械の加工誤差のうち、**熱変位が占める割合は最大75 %** であり、最大の誤差要因である。
- 対策は「**熱的に素性の良い設計（回避）**」と「**モデルによる補償＝熱変位補正**」の2系統だが、
  **回避はコストが高く**、実用の主役は**補償**である。その前提は熱変形を正しく予測できることである。
- しかし予測が難しい。**TCP変位は温度の絶対値ではなく温度勾配で決まる**ため、
  **温度の傾向が合っても変形の傾向が合うとは限らない**。
  さらに**熱伝達率を直接測る計測器が存在せず**、境界条件の不確かさがTCP予測精度を直接左右する。

### ゆえに求められること（1行）

> **限られた実測から境界条件ごと補正し、熱変形そのものを推定する枠組みが求められる。**

### 図（2枚並べる）

| 位置 | ファイル | 内容 |
|---|---|---|
| 左 | `docs/img/intro_temp_history.png` | 測温2点の温度の時間変化（温度差2.9 K） |
| 右 | `docs/img/intro_deform_anim.gif` | 熱変形の3D（誇張）＋変位の時間変化（36フレームのアニメ） |

レイアウト見本：`opencae_2026_intro_slide.html`（ブラウザで開くと確認できます）

---

## 2. 各行の根拠（すべて無料で読める出典）

### 2-1. 1行目：最大75 %

> ⚠️ **訂正（2026-09-29）**：当初「40〜70 %［Bryan 1990］」としていましたが、
> Bryan (1990) は有料で**原典を確認できませんでした**。
> **全文無料で確認できた「最大75 %」に差し替えます。**

**根拠①（オープンアクセス・原文確認済み）**

**Li, Z., Vogl, G. W., Kinzel, E. C., Santa, B., Landers, R. G. (2024)**
"Machine Tool Thermal Error Measurement and Prediction via Wireless Microscope",
*Manufacturing Letters* **41**, 1440–1451（NAMRC 52, CC BY-NC-ND）
📄 全文無料: <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076>

> "Thermal errors can contribute **up to 75 percent** of the overall machining errors of a machined part."

**根拠②（arXiv・原文確認済み）**

**Bünger, A., Herzog, R., Naumann, A., Stoll, M. (2023)**
"Uncertainty Propagation of Initial Conditions in Thermal Models", arXiv:2306.12736
📄 全文無料: <https://arxiv.org/abs/2306.12736>

> "According to **Mayr et al., 2012, the thermal error accounts for 75 %** of the total
> manufacturing error in the final product."

→ **どちらも Mayr et al. (2012) を典拠**としており、孫引きの連鎖が確認できる。

**スライドでの書き方**：「熱変位は加工誤差の**最大75 %**を占める［Mayr et al. 2012］」

### 2-2. 2行目：対策は「回避」と「補償」

**Li et al. (2024)**（上記、全文無料）の Introduction より原文引用:

> "In general, there are two methods for thermal error reduction:
> **thermal error avoidance** and **thermal error compensation**."

> "Thermal error avoidance techniques make the machine tool less sensitive to temperature variations.
> These techniques include the use of materials with reduced friction and low coefficients of
> thermal expansion, optimization of lubrication and cooling systems, etc."

> "**Thermal error avoidance is typically a more costly solution than thermal error compensation**
> and is more sensitive to modeling errors and unknown disturbances."

→ **「素性の良い機械を作る（回避）」だけでは高くつき、補償が実務の主役**、という構図が
文献で明言されている。これが1行目と3行目の橋渡しになる。

（参考：**Mayr, J. et al. (2012)** "Thermal issues in machine tools", *CIRP Annals* **61**(2), 771–791。
この分野の標準レビュー。研究領域を「測定」「計算」「低減・補償」に整理。
🔗 <https://www.nist.gov/publications/thermal-issues-machine-tools> ※**有料**、無料PDFは見つからず）

### 2-3. 3行目：温度が合っても変形は合わない

**根拠①：境界条件と初期温度が精度を直接左右する**（原文確認済み・無料）

**Bünger et al. (2023)** arXiv:2306.12736 📄 <https://arxiv.org/abs/2306.12736>

> "the accuracy of the **TCP prediction depends highly on the accuracy of the model parameters,
> such as heat exchange parameters, and the initial temperature**"

**根拠②：熱伝達率（CHTC）を直接測る計測器が存在しない**

> ⚠️ 前回「複数の論文で指摘」と書いて具体名を挙げていませんでした。以下が該当します。
> **いずれも無料公開**ですが、私の環境からは本文取得がブロックされたため
> **検索結果に基づく記述**です。スライドに載せる前に原文をご確認ください。

1. "Intelligent Soft Sensor for Spindle Convective Heat Transfer Coefficient Under Varying
   Operating Conditions Using Improved Grey Wolf Optimization Algorithm" (2025)
   📄 無料全文（PMC）: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12473924/>
   → 「**CHTCを直接測定する専用計測器が存在しない**ため、主軸の熱解析は大きな困難に直面する」
   「従来の熱流束センサや非接触赤外サーモグラフィは温度分布は得られるが、
   **CHTCを直接測ることはできず**、熱モデルの精度を制限している」の趣旨

2. "The Thermal Error Estimation of the Machine Tool Spindle Based on Machine Learning",
   *Machines* **9**(9), 184 (2021)。MDPI **オープンアクセス**
   📄 <https://www.mdpi.com/2075-1702/9/9/184>
   → 「FEMの精度は熱源・熱伝達係数・境界条件が明確に定義されているかに依存する。
   しかし主軸は材質の異なる多数の部品からなり、熱源と境界条件は組立や加工条件に強く依存するため、
   **汎用的に通用する定義は難しい**」の趣旨

**根拠③：勾配とてこ、そして変位の向きの反転**

- 柱の上下にわずかな温度勾配があるだけで曲げが生じ、主軸先端・TCPでは
  **てこで拡大されて大きな変位**になる → 「平均温度が合っている」ことは保証にならない
- **Mayr et al. (2012)** では、**異なる時定数が2つ関与する場合や、熱源からの急な勾配が
  時間とともに均一化していく過程で、TCP変位が運転中に向きを変えうる**ことが指摘されている
  ※こちらは**有料論文の記述**であり、検索結果に基づく。スライドに載せるなら要確認

**スライドでの書き方（安全な版）**：
「TCP変位は温度勾配で決まる。**熱伝達率は直接測る計測器が存在せず**［PMC 2025］、
その不確かさがTCP予測精度を直接左右する［Bünger et al. 2023］」

---

## 3. 熱変位補正は実際に何をしているのか

### 3-1. 結論（先に答え）

| 質問 | 答え |
|---|---|
| 「補償」＝熱変位補正か？ | **はい**。文献は対策を avoidance（回避）と compensation（補償＝熱変位補正）に二分 |
| 何をしている？ | 機械各所の**温度を測り**、事前に取った**温度→刃先変位の関係モデル**で変位を予測し、**その逆符号を軸の指令位置に足し込む** |
| 「多点温度＋変位実測→回帰でモデル化」で合っているか | **合っています**（下記の原文引用のとおり） |
| いつ補正する？ | **運転中に連続して**。CNCが温度を取り込み**オフセットを常時更新**する |
| 機械に埋め込まれているか | **はい**。オークマOSP、マザックAI Thermal Shield、DMG MORI SGS など**CNC側の標準／オプション機能** |

### 3-2. 文献での記述（Li et al. 2024、全文無料）

> "In contrast, **thermal error compensation is normally based on a predictive model established by
> the thermal error measurement of a machine tool**."

> "Thermal error compensation strategies typically employ a **predictive error model, which are
> inverted to determine compensation amounts**. Common thermal error models include
> **least-square regression, finite element, neural network, gray system**, etc."

→ 「いろいろな温度条件で多点温度と刃先変位を測り、その関係をモデル化し、
**モデルを逆に使って補正量を決める**」という理解は**文献の記述と一致**。
線形回帰・FEM・ニューラルネット・グレーシステムが代表例として挙がっている
（**ガウス過程回帰もこの系列の一つ**）。

### 3-3. メーカー各社の実装

> ⚠️ 以下は**各社の公開情報（製品ページ）に基づく**もので、学術論文による検証ではありません。
> 数値は各社の公称値です。

#### オークマ ― サーモフレンドリーコンセプト（TFC）
🔗 <https://www.okuma.co.jp/onlyone/thermo/>

考え方が特徴的で、温度変化を抑え込むのではなく **"受け入れる"**。

| 要素 | 内容 |
|---|---|
| 設計側 | 「熱変形の単純化構造」「温度分布均一化」により**機械を素直に変形させ**、ねじれ・傾きを抑えて**熱変位を予測可能な状態にする** |
| **TAS-C**（環境熱変位制御） | **適切に配置されたセンサの温度情報＋送り軸の位置情報**から、環境温度変化による構造体の熱変位を**推定して制御** |
| **TAS-S**（主軸熱変位制御） | 主軸の温度変化を**リアルタイム管理**し、回転・停止時の熱変位を制御 |
| 実績 | 2019年10月時点で全127機種中83機種に搭載、出荷5万台超 |

**本研究との関係**：TAS-Cの「センサ温度＋軸位置から熱変位を推定」は、
本研究の「**少数観測から状態を推定**」と同じ構図。
ただしオークマは**設計で変形を単純化してから**推定する点が巧み。

#### ヤマザキマザック ― AI Thermal Shield / Intelligent Thermal Shield
🔗 <https://www.mazak.com/jp-ja/technology/accuracy/>

| 要素 | 内容 |
|---|---|
| 入力 | **主軸回転速度**＋**温度センサ情報** |
| 考慮する条件 | 温度変化、**機械位置**、**クーラントON/OFF** |
| 学習 | 加工後の計測データを**蓄積・学習**し、顧客の加工環境に適した補正へ最適化 |
| 公称性能 | **室温が8 ℃変化しても連続加工精度6 µmを維持** |
| ハード側 | 主軸軸受・モータ外筒に温調冷却油、**X/Y/Z全軸のボールねじ軸冷却を標準装備** |

**注目点**：**クーラントON/OFFを入力に入れている**（冷却液は熱境界条件そのもの）。

#### DMG MORI ― Spindle Growth Sensor（SGS）
🔗 <https://en.dmgmori.com/products/machines/milling/vertical-milling/nvx/nvx-5100>

| 要素 | 内容 |
|---|---|
| **SGS** | 主軸の**軸方向伸びを直接センシング**して補正（モデル予測ではなく**実測フィードバック**） |
| 構造側 | 冷却されたリニアガイド、主軸成長のアクティブ制御 |
| 制御 | 機械各所の温度センサで**リアルタイム監視**し、**工具オフセットを自動調整** |

**注目点**：**変位を直接測る**アプローチ。本研究の「変位観測を足すと精度が上がる」
（温度2点 0.197 K → ＋変位 0.159 K）と**同じ思想が既に商用化されている**。

#### FANUC ― CNC側の機能
🔗 <https://www.fanuc.co.jp/ja/product/cnc/index.html>

- **AI熱変位補正**（2022年6月発表）
- 多点の温度センサ入力に対応する小型・省配線のI/Oユニットを提供
→ **補正は最終的にCNCの中で座標オフセットとして効く**という実装の裏付け。

### 3-4. 「いつ補正するのか」

1. **事前（オフライン）**：さまざまな温度条件で運転し、**多点温度**と**刃先変位**を
   同時計測 → 回帰などで**予測モデル**を作る
2. **運転中（オンライン）**：CNCが温度センサを**常時読み込み**、モデルで変位を予測 →
   **符号を逆にして軸の指令位置に加算**
   （オークマTAS-Sは「リアルタイム管理」、DMG MORIは「real-time／自動調整」と明記）
3. **随時（学習）**：マザックは加工後の計測データを蓄積・学習して補正を更新

→ **「加工の合間に測って直す」のではなく、運転中ずっと補正がかかり続けている**のが標準。

---

## 4. 最も近い先行研究（質疑への備え）

### 4-1. Teshima, Y., Tanaka, S., Kizaki, T., Sugita, N. (2024)

"Sensor placement strategy based on reduced-order models for thermal error estimation in machine tools",
*CIRP Journal of Manufacturing Science and Technology* **55**, 403–410
🔗 <https://doi.org/10.1016/j.cirpj.2024.10.015> ※**有料**（アブストラクトは無料）

- 温度入力とTCP変位を結ぶ**縮約モデル（ROM）**を作り、**回帰係数を各点の温度感度**とみなす
- 冒頭に「thermal errors in machine tools account for **up to 70 %** of machining errors」
- 彼ら自身が挙げるギャップ：
  「入力点を増やせば精度は上がるが、**最適な点数と配置を決める方法がこれまで無かった**」

### 4-2. Ando, S., Tanaka, S., Teshima, Y., Morishita, J., Kizaki, T.（続報）

"Strategy for Sensor Placement to Estimate Thermal Errors Using **Temperature-Sensitivity
Distribution** Based on a Reduced-Order Model of Machine Tools"
**ICTIMT2025**（4th International Conference on Thermal Issues in Machine Tools）論文集
🔗 <https://doi.org/10.1007/978-3-032-01194-7_31> ※**有料**
研究室: 東京大学 先端加工学（木崎研） <https://mfg.t.u-tokyo.ac.jp/>

- **伝達関数行列**を定義し、TCP-ワーク間の相対変位をセンサ位置の温度変化と結ぶ
- ROMに基づく**センサ感度関数**を定義して配置を決める

### 4-3. 本研究との違い

| | Teshima / Ando ら | 本研究 |
|---|---|---|
| ROMの作り方 | 温度→TCP変位の**伝達関数・回帰** | **POD＋Q-DEIM**で代表点を選び集中定数ROMを校正 |
| 配置の基準 | **温度感度分布**（伝達関数から） | **観測の価値** $\Delta=\mathrm{Cov}(X,y)^2/(\mathrm{Var}(y)+r)$ ＝相互情報量。**推定対象ごとに変わる** |
| 未知パラメータ | 扱わない（温度は入力） | **発熱量 $Q$ ・放熱 $h$ を状態に入れて同時推定**（EnKF） |
| 観測の種類 | 温度センサ | **温度＋変位**を同じ枠組みで扱う |
| 実装 | ― | **OSSのみ**（OpenFOAM＋FrontISTR＋Python）で全公開 |

→ **「ROMでセンサ配置を決める」流れは既にある。本研究の差分は
①データ同化で境界条件ごと推定する ②温度と変位を同じ価値尺度で比べる
③OSS一気通貫で再現可能にした、の3点。**

---

## 5. 発表での使い方

- オープンCAEの聴衆向けなので、**1行目で数字を出して「大きい問題だ」と示し**、
  **3行目で「だからシミュレーションだけでは足りない」へ落とす**のが効く。
- 2行目で「補償が実務の主役」「各社が実装済み」と言っておくと、
  **本研究が実務の延長線上にある**ことが伝わる。
- 本研究の `oi_fullsolver_compare.png`（温度RMSEは0.09 Kまで下がるのに $Q$ は8.9 W止まり）は
  **まさに3行目の構図**なので、導入の伏線として使える。

**→ そして本研究：「境界条件そのものを観測から推定してしまえばよい」** につながる。

---

## 付録：主張ごとの検証状況

| 主張 | 出典 | 無料で読めるか | 検証 |
|---|---|---|---|
| 熱誤差は最大75 % | Li et al. 2024 / Bünger et al. 2023 | ✅ 両方とも全文無料 | ✅ **原文確認済み** |
| avoidance と compensation の2分類 | Li et al. 2024 | ✅ 無料 | ✅ **原文確認済み** |
| 回避の方がコストが高い | Li et al. 2024 | ✅ 無料 | ✅ **原文確認済み** |
| 補正は予測モデルを逆に使う | Li et al. 2024 | ✅ 無料 | ✅ **原文確認済み** |
| モデルは回帰・FEM・NN・グレー系 | Li et al. 2024 | ✅ 無料 | ✅ **原文確認済み** |
| 熱交換パラメータ・初期温度が精度を左右 | Bünger et al. 2023 | ✅ 無料 | ✅ **原文確認済み** |
| CHTCを直接測る計測器が無い | PMC12473924 (2025) | ✅ 無料 | ⚠️ 検索結果のみ（**要原文確認**） |
| FEMの境界条件定義が難しい | *Machines* 9(9) 184 (2021) | ✅ 無料 | ⚠️ 検索結果のみ（**要原文確認**） |
| 変位の向きが途中で変わる | Mayr et al. 2012 | ❌ 有料 | ⚠️ 検索結果のみ（**要原文確認**） |
| 研究領域の3分類 | Mayr et al. 2012 | ❌ 有料 | ⚠️ 検索結果のみ |
| 40〜70 % | Bryan 1990 | ❌ 有料 | ❌ **未確認。使用非推奨** |
| 各社の実装内容 | 各社公開ページ | ✅ 無料 | ⚠️ メーカー公称値 |
