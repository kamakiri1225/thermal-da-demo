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
| **TAS-C**（Thermo Active Stabilizer - Construction／環境熱変位制御） | **適切に配置されたセンサの温度情報＋送り軸の位置情報**から、環境温度変化による構造体の熱変位を**推定して制御** |
| **TAS-S**（Thermo Active Stabilizer - Spindle／主軸熱変位制御） | 主軸の温度情報に加え、**主軸回転・回転速度変更・停止**といった状況変化も考慮。回転速度が頻繁に変わっても熱変位を制御 |
| 実績 | 2019年10月時点で全127機種中83機種に搭載、出荷5万台超 |

**オークマは何が他社と違うのか（3点）**

1. **「熱を抑える」のではなく「予測しやすい変形にする」**
   他社が熱を**入れない／追従する／実測する**方向なのに対し、オークマだけが
   **設計段階で「モデル化しやすさ」を作り込んでいる**。
   「熱変形の単純化構造」「温度分布均一化」により、**複雑なねじれ・傾きを排して
   単純な伸び縮みに落とす**。結果、少ないセンサ・単純なモデルで当たるようになる。
2. **機械とCNCを両方自社で作っている**
   制御装置OSPが自社製なので、**構造設計と補正アルゴリズムを一体で最適化**できる。
   「設計で変形を単純化 → その単純な変形を自社CNCで補正」という往復が社内で閉じる。
3. **送り軸の位置情報を入力に使っている**
   熱変形の効き方は**機械が今どこにいるか**で変わる（構造の伸びが軸位置で拡大される）。
   公開情報の範囲では、**軸位置を明示的に入力に入れているのはオークマの特徴**。

**本研究の言葉で言い直すと**

オークマがやっているのは、**「場を低次元にする設計」**です。
本研究でPODをかけたら**2モードで99.9 %**を説明できましたが、これは
**対象がそもそも低次元だったから**当たったのであって、条件が複雑なら当たりません。

> **オークマは、その低次元性を"設計変数"として作り込んでいる。**
> つまり **モデル化しやすさそのものを設計している**のが、他社との決定的な違い。

一方、本研究が扱うのは**「与えられた（設計を変えられない）機械で、
どこを測れば何が分かるか」**という問題。
**オークマ＝設計で解く／本研究＝観測設計で解く**、という補完関係にあります。

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

#### 牧野フライス製作所 ― 冷却と断熱で「熱を入れない」
🔗 <https://www.makino.co.jp/ja-jp/>

4社の中で**最も「回避（avoidance）」寄り**のアプローチ。

| 要素 | 内容 |
|---|---|
| 主軸 | 高速主軸に**軸芯冷却＋ジャケット冷却**の両方を採用 |
| 送り系 | **ボールねじ内部**とX軸案内面を集中冷却し、機械本体温度に同調させる |
| 機体全体 | **断熱カバー**で温度変化速度を均一化、**機内空気の循環・攪拌** |
| 熱源 | 温調油による冷却 |
| 補正 | **熱伝導遅れを考慮した熱変位補正**を併用 |

**注目点**：「**熱伝導遅れを考慮した**補正」＝温度が変わってから変形が現れるまでの
**時間遅れを明示的にモデル化**している。これは本研究の集中定数ROM
（ $C_i\,dT_i/dt=\sum_j K_{ij}(T_j-T_i)+\dots$ が時定数を持つ）と同じ発想。

#### FANUC ― CNC側の機能
🔗 <https://www.fanuc.co.jp/ja/product/cnc/index.html>

- **AI熱変位補正**（2022年6月発表）
- 多点の温度センサ入力に対応する小型・省配線のI/Oユニットを提供
→ **補正は最終的にCNCの中で座標オフセットとして効く**という実装の裏付け。

### 3-4. メーカー4社の実装まとめ（一覧）

> ⚠️ 各社の**公開情報（製品ページ・技術紹介）**に基づく。数値は各社の公称値であり、
> 第三者による検証ではありません。非公開の実装詳細は当然含まれません。

#### (a) 基本方針の違い

| | **オークマ** | **ヤマザキマザック** | **DMG MORI** | **牧野フライス** |
|---|---|---|---|---|
| 技術名 | サーモフレンドリーコンセプト（TFC）<br>TAS-C / TAS-S ＝ **Thermo Active Stabilizer**（読み：ティー・エー・エス） | AI Thermal Shield／Intelligent Thermal Shield | Spindle Growth Sensor（SGS）ほか | （個別技術の集合） |
| **一言でいうと** | **変形を予測可能な形にしてから推定する** | **学習で補正を育てる** | **変位を直接測って直す** | **そもそも熱を入れない** |
| 回避 vs 補償 | 回避（設計）→補償 の**二段構え** | 補償が主、ハード冷却も併用 | **実測フィードバック**が主 | **回避が主**、補償で仕上げ |
| 思想 | 温度変化を **"受け入れる"** | 現場ごとに**最適化していく** | **モデルを信じず測る** | **熱源を断つ・均一化する** |

#### (b) 何を測り、何を入力にしているか

| | オークマ | マザック | DMG MORI | 牧野 |
|---|---|---|---|---|
| 温度センサ | **適切に配置された**センサ群 | 主軸まわり＋環境 | 機械各所 | 機体各所（同調冷却用） |
| 変位センサ | （公開情報になし） | （公開情報になし） | **SGS＝主軸伸びを直接測定** | （公開情報になし） |
| **特徴的な入力** | **送り軸の位置情報** | **クーラントON/OFF**、機械位置、主軸回転速度 | **主軸伸びの実測値** | **熱伝導の時間遅れ** |
| 補正の作り方 | 機械の熱変位特性に基づく推定 | 加工後の計測データを**蓄積・学習** | 実測に基づくオフセット自動調整 | 遅れを考慮した補正式 |

#### (c) ハード側の熱対策

| | オークマ | マザック | DMG MORI | 牧野 |
|---|---|---|---|---|
| 主軸冷却 | ○ | **軸受・モータ外筒に温調冷却油** | ○ | **軸芯冷却＋ジャケット冷却** |
| 送り系冷却 | ○ | **X/Y/Z全軸のボールねじ軸冷却を標準装備** | **冷却されたリニアガイド** | **ボールねじ内部＋X軸案内面を集中冷却** |
| 構造・環境 | **熱変形の単純化構造**、温度分布均一化 | ― | ― | **機体全体の断熱カバー**、**機内空気の循環・攪拌** |

#### (d) 公称性能・実績

| | 内容 |
|---|---|
| オークマ | 2019年10月時点で**全127機種中83機種に搭載、出荷5万台超** |
| マザック | **室温が8 ℃変化しても連続加工精度6 µmを維持** |
| DMG MORI | 冷却リニアガイド＋主軸成長のアクティブ制御で「**長期安定精度**」 |
| 牧野 | V99等で「**ミクロンオーダーの中・大物金型加工に最適**」 |

#### (e) 本研究との対応

| 各社がやっていること | 本研究での対応物 |
|---|---|
| オークマ：センサ温度＋軸位置から**構造体の熱変位を推定** | **少数観測から全場を推定**（gappy-POD＋EnKF） |
| オークマ：設計で**変形を単純化して予測可能に** | PODで**場が2モードで99.9 %**＝低次元であることを利用 |
| マザック：**クーラントON/OFF**を入力に | クーラント＝**熱境界条件**。本研究は $h$ を**状態として推定** |
| マザック：加工後データで**学習** | EnKFの**逐次更新**（毎サイクル共分散を作り直す） |
| DMG MORI：**変位を直接測る** | 「**変位観測を足すと精度が上がる**」（温度2点0.197 K → ＋変位0.159 K）を定量化 |
| 牧野：**熱伝導遅れ**を考慮 | 集中定数ROM $C_i\,dT_i/dt=\sum_j K_{ij}(T_j-T_i)+\dots$ が**時定数を持つ** |

**→ 各社が経験的に到達している設計を、本研究は
「観測の価値 $\Delta$ 」という一つの尺度で説明・設計できる形にしている。**

### 3-5. 「いつ補正するのか」

1. **事前（オフライン）**：さまざまな温度条件で運転し、**多点温度**と**刃先変位**を
   同時計測 → 回帰などで**予測モデル**を作る
2. **運転中（オンライン）**：CNCが温度センサを**常時読み込み**、モデルで変位を予測 →
   **符号を逆にして軸の指令位置に加算**
   （オークマTAS-Sは「リアルタイム管理」、DMG MORIは「real-time／自動調整」と明記）
3. **随時（学習）**：マザックは加工後の計測データを蓄積・学習して補正を更新

→ **「加工の合間に測って直す」のではなく、運転中ずっと補正がかかり続けている**のが標準。

---

## 4. 計測の実際 ― 何を、どこで、どう測るか

### 4-1. 測定方法の2分類（文献の整理）

**Li et al. (2024)**（全文無料 📄 <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076>）
の Introduction が、測り方をきれいに2つに分けています（原文引用）:

> "**Indirect measurement methods** estimate the thermal deformation of a machine tool based on a
> temperature-deformation simulation model (e.g., finite element, neural network) and
> **measure temperatures via thermocouples or thermal cameras**."

> "Alternatively, **direct measurement methods** measure the changes of the tool position relative to
> the workpiece via **displacement sensors such as laser interferometers, capacitive sensors,
> and touch probes**."

> "For direct measurements, **multi-position measurements are necessary** to determine the thermal
> errors effects on the entire machine tool workspace, which is usually a **very time-consuming**
> measurement process."

| | 間接法（indirect） | 直接法（direct） |
|---|---|---|
| 測るもの | **温度** | **刃先とワークの相対変位** |
| 道具 | 熱電対、**サーマルカメラ** | レーザ干渉計、静電容量センサ、**タッチプローブ** |
| 変位の求め方 | 温度→変形モデル（FEM・NN）で予測 | そのまま実測 |
| 長所 | 運転中も常時測れる、安い | モデル誤差が無い |
| 短所 | **モデルの精度に依存**（境界条件の不確かさ） | **多点を測るので時間がかかる**、加工を止める必要 |

**→ 実機の補正は「間接法（温度）」が主役。直接法は事前のモデル作成と検証に使う**、という役割分担。

### 4-2. 温度はどう測るか・どこで測るか

**測る道具**

| 種類 | 特徴 | 用途 |
|---|---|---|
| **熱電対**（K型など） | 安い、応答が速い、多点化しやすい | **最も一般的**。機械各部に多数貼る |
| 測温抵抗体（Pt100） | 精度・安定性が高い、やや高価 | 基準温度・精密測定 |
| サーミスタ | 感度が高い、狭い範囲 | 組み込み用 |
| **サーマルカメラ** | **非接触で面分布**が一度に取れる | 熱源の特定、モデル検証 |

**どこに置くか** ― ここが本研究のテーマそのものです。各社の公開情報では:

- オークマ：「**適切に配置されたセンサ**の温度情報」（配置が効くことを示唆）
- マザック：主軸まわり＋環境温度
- 一般には：**主軸軸受近傍、主軸モータ、ボールねじナット、案内面、コラム、ベッド、機外環境**

**→ 「どこに置けば効くか」に定説が無いのが、先行研究（Teshima 2024ほか）と本研究の出発点。**

### 4-3. 温度条件はどう取るか ― ISO 230-3

**ISO 230-3:2020** "Test code for machine tools — Part 3: Determination of thermal effects"
🔗 <https://www.iso.org/standard/73291.html> ※本体は有料
📄 **無料プレビュー（Scope全文が読めます）**:
<https://cdn.standards.iteh.ai/samples/73291/b10e76761d1945c6b7648d2fea86b6a8/ISO-230-3-2020.pdf>

規定される試験は**4種類**（2020年版。2007年版は3種類）:

1. **ETVE**（Environmental Temperature Variation Error）＝**環境温度変動誤差**試験
   … 機械を動かさず、**周囲温度の変動だけ**でどれだけ動くかを測る
2. **主軸回転**による熱変形試験（主軸を回し続けて変位を測る）
3. **直線軸の運動**による熱変形試験
4. **回転軸の運動**による熱変形試験（2020年版で追加）

**ここが重要** ― 規格自身が次のように断っています（プレビューから原文引用）:

> "It is a recognized fact that **the ultimate thermo-elastic deformation of a machine tool is
> closely linked to the operating conditions**. **The test conditions described in this document are
> not intended to simulate the normal operating conditions** but are to facilitate performance
> estimation..."

> "For example, **use of coolants can significantly affect the actual thermal behaviour of the
> machine tool**. Therefore, these tests are considered only as the **preliminary tests** towards
> the determination of actual thermo-elastic behaviour of the machine tool..."

> "The tests are designed to measure the **relative displacements between the component that holds
> the tool and the component that holds the workpiece**."

**→ 発表で効く一文**：
**「標準試験は実運転を模擬していない、と規格自身が明記している。
とくにクーラントの影響は大きい」**
つまり**実運転の熱挙動は、標準試験だけでは決まらない**。
だからこそ**運転中の観測からその場で推定する**枠組みに意味がある。

### 4-4. 変位はどう測るか

#### タッチプローブとは

**主軸に装着する接触式の測定器**です。先端の球（スタイラス）がワークや基準球に
**触れた瞬間に信号を出す**ので、そのときの機械座標を読めば「その面・その点がどこにあるか」が分かります。

- 用途：ワーク原点出し、加工後の寸法測定、そして**熱変位の測定**
- 熱変位測定での使い方：**基準球（マスターボール）を定盤に固定**しておき、
  一定時間ごとにプローブで測る → **球の見かけの位置が動いた量＝その間に生じた熱変位**
- 長所：**機械に元から付いていることが多い**（追加投資が小さい）、3方向まとめて測れる
- 短所：**接触するたびに加工を止める**必要がある、測定に時間がかかる

#### そのほかの変位測定

| 方法 | 原理 | 特徴 |
|---|---|---|
| **静電容量変位センサ** | 対象との距離で静電容量が変わる | 非接触、**nm〜サブµm分解能**。主軸の伸びを連続測定できる |
| **渦電流センサ** | 金属中の渦電流で距離を検出 | 非接触、油・切粉に強い。実機に組み込みやすい |
| **レーザ干渉計** | 光の干渉で長さを測る | 最高精度。**軸の位置決め誤差**の評価に使う。設置が大がかり |
| **テストバー（マンドレル）＋変位センサ** | 主軸に基準の丸棒を付け、周囲から複数センサで測る | **ISO 230-3の主軸試験の標準構成**。軸方向の伸びと傾きを分離できる |
| **DMG MORI の SGS** | 主軸伸びの専用センサ | 実機に**常設**して運転中ずっと測る |

**補足**：ISO 230-3では、角度偏差も評価する場合、
「**変位センサの間隔**は十分な測定範囲・分解能・不確かさが得られるように選ぶ」とされ、
また「**軸方向センサは主軸端面に直接当てて**、テストバー自身の熱膨張の影響を除いてもよい」
とされています（プレビュー記載の趣旨）。

#### カメラで測れるか ― **測れます**（2024年のオープンアクセス論文）

**Li, Z., Vogl, G. W., Kinzel, E. C., Santa, B., Landers, R. G. (2024)**
"Machine Tool **Thermal Error Measurement and Prediction via Wireless Microscope**",
*Manufacturing Letters* **41**, 1440–1451
📄 **全文無料**: <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076>

やっていること（原文の要約）:

> "A novel method is proposed to measure the thermal errors of a three-axis machine tool by
> **taking images of unique custom-designed fiducials attached to a worktable using a wireless
> microscope mounted to the spindle**."

- **主軸にワイヤレス顕微鏡**を装着し、**テーブルに貼った基準マーク（フィデューシャル）**を撮影
- マークはフェムト秒レーザで**黒アルミ板に加工**（楕円の格子、間隔350 µm、加工分解能75 nm）
- 楕円1つ1つに**行・列の番号が振ってある**ので、**格子間隔より大きくずれても追跡できる**
- 画像解析で**3方向の変位**を算出

**精度**（原文引用）:

> "The results show the method can **measure thermal errors within four times the positioning
> resolution of the machine tool**, with most of the errors being **smaller in magnitude than twice
> the positioning resolution**."

**長所**：測定が**速く・安く・実用的**（著者らの主張）。ワイヤレスなので配線不要。
複数マークを置けば**作業空間の複数位置**を測れる。
**短所**：加工中は測れない（主軸に顕微鏡を付ける必要がある）、視野内の平面変位が主。

**→ 「カメラでも測れるか？」の答え：YES。2024年のNAMRC論文で実証済み**で、
しかも**全文無料で読めます**。

---

## 5. 最新研究の動向（2023〜2025）

| 論文 | 何をしているか | 無料か |
|---|---|---|
| **Teshima et al. (2024)** *CIRP JMST* 55, 403–410<br>🔗 <https://doi.org/10.1016/j.cirpj.2024.10.015> | 温度→TCP変位の**ROM**を作り、回帰係数＝温度感度として**センサ配置の指針**を与える | ❌ 有料 |
| **Ando et al. (ICTIMT2025)**<br>🔗 <https://doi.org/10.1007/978-3-032-01194-7_31> | 上の続報。**伝達関数行列**と**温度感度分布**で配置を決める | ❌ 有料 |
| **Bünger et al. (2023)** arXiv:2306.12736<br>📄 <https://arxiv.org/abs/2306.12736> | **初期条件の不確かさ伝播**を低ランク近似で高速計算し、**センサ配置の評価**に使う | ✅ 無料 |
| **Li et al. (2024)** *Manufacturing Letters* 41<br>📄 <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076> | **ワイヤレス顕微鏡＋基準マーク**で熱変位を直接測る新手法 | ✅ 無料 |
| **Coelho et al. (2025)** arXiv:2510.03261<br>📄 <https://arxiv.org/abs/2510.03261> | 温度・熱流束**場そのもの**をニューラルネットで予測。RNN/GRU/LSTM/Transformer等**6種をベンチマーク**。**相関ベースで測定点を選ぶ**戦略も提案 | ✅ 無料 |
| **Intelligent Soft Sensor for CHTC** (2025)<br>📄 <https://pmc.ncbi.nlm.nih.gov/articles/PMC12473924/> | **測れないCHTCを最適化アルゴリズムで推定**する（ソフトセンサ） | ✅ 無料 |

### 熱の「データ同化」そのものを扱った研究

本研究に**手法として最も近い**2本。どちらも「物理モデル＋データ同化で温度場を推定する」構図です。

#### Bünger, A., Herzog, R., Naumann, A., Stoll, M. (2023) ― 工作機械・センサ配置の評価

"Uncertainty Propagation of Initial Conditions in Thermal Models", arXiv:2306.12736
📄 **全文無料**: <https://arxiv.org/abs/2306.12736>

**何を主張しているか**（本文から要約）

1. **問題設定**：工作機械のTCP位置を知るには、幾何的に結合した熱伝導方程式群（FEM）が要る。
   その予測精度は**モデルパラメータ（熱交換パラメータ）と初期温度の精度に強く依存**する。
2. **やったこと**：**初期温度場を未知**として、温度センサの時系列から推定する**逆問題**を、
   **ベイズの枠組み**で解く。関心は推定値そのものではなく、
   **事後分散（posterior variance）＝推定がどれだけ絞り込めたか**の評価。
3. **技術的な壁**：状態空間が巨大（FEM節点数）なのに出力（センサ）はごく少数。
   事後共分散行列は**密で巨大**になり、そのまま扱えない。
4. **解決策**：データ誤差ヘッシアンを**低ランク近似**（主要固有対のみ）。
   さらに**テンソルトレイン（TT）分解**による方法を提案し、直接法と比較。
   → **TT法はCPU時間は増えるがメモリは劇的に少ない**（FEM行列と同程度）。
5. **目的**：「与えられたセンサ配置の良さを評価すること」は
   **最適センサ配置の必須の前提**である、という位置づけ。

**本研究との関係**

| | Bünger et al. 2023 | 本研究 |
|---|---|---|
| 推定対象 | **初期温度場** | **温度場＋発熱量 $Q$ ＋放熱 $h$** |
| 手法 | ベイズ逆問題＋低ランク／TT近似 | **EnKF**（アンサンブルで共分散を近似） |
| 低次元化 | 事後共分散を低ランク近似 | **POD＋Q-DEIM**で状態自体を5点に縮約 |
| 配置の扱い | 「**評価**する」道具を提供 | **評価＋対象ごとの比較**（ $\Delta$ ＝相互情報量） |
| 観測の種類 | 温度のみ | **温度＋変位** |

→ **「配置を評価する前提を整える」のがBünger。本研究は「対象が変われば最適配置も変わる」
を示す側**。相補的で、質疑で引き合いに出せる。

#### Peet, B. J. A. (2019) ― 赤外線シグネチャ監視・屋外物体の表面温度推定

"Accurate estimation of temperature distributions for IR signature monitoring with a dynamic
thermal model and data assimilation", *Proc. SPIE* **11158**, Target and Background Signatures V,
111580C. DOI: 10.1117/12.2532755
🔗 SPIE: <https://doi.org/10.1117/12.2532755> ※**有料**
🔗 TNOリポジトリ（書誌・キーワード。本文はメール請求制）:
<https://repository.tno.nl/islandora/object/uuid:3886b677-1d88-4198-9925-39d68360fb77>

**何を主張しているか**（公開されている記述に基づく）

- **目的**：赤外線シグネチャ（物体がIRセンサにどう見えるか）の監視には、
  **物体表面の温度分布**をリアルタイムで知る必要がある。
- **やったこと**：TNOが、**標準的な気象観測データだけを入力**として、
  物体表面の温度分布の**実時間変化を推定するモデル**を開発した。
- **方法**：気象データから表面要素への**熱流束を推定**し、**動的な熱収支を有限要素法で数値的に解く**。
  そこに**データ同化（カルマンフィルタ）**を組み合わせる
  （TNO登録キーワード：IR signatures / Thermal modelling / **Data assimilation** /
  Finite element method / **Kalman filter**）。

**本研究との関係**

| | Peet 2019 | 本研究 |
|---|---|---|
| 分野 | **赤外線シグネチャ**（艦船・車両の被発見性） | **工作機械の熱変形** |
| 物理モデル | FEMによる動的熱収支 | CHT（OpenFOAM）→縮約ROM |
| 同化手法 | **カルマンフィルタ** | **EnKF**（非線形・パラメータ拡大に対応） |
| 入力の不確かさ | **気象データから熱流束を推定** | **発熱量 $Q$ ・放熱 $h$ を状態として推定** |
| 出力 | 表面温度分布 | 温度場**＋熱変形** |
| 観測設計 | （主題でない） | **観測の価値 $\Delta$ で配置を設計** |

**→ 重要な示唆**：「**熱の問題にデータ同化を持ち込む**」という発想は、
工作機械分野の外（**IRシグネチャ監視**）では2019年時点で既に実用検討されている。
同じ構図（**測れない境界条件＝環境からの熱流束を、観測で補正しながら温度場を推定**）であり、
**本研究のアプローチが分野を超えて筋の良いものである傍証**として使える。
一方、**熱変形（変位）まで出す点**と**観測設計を定量化する点**は本研究の独自性。

### 動向の読み取り

1. **「誤差を直接予測」から「場を予測」へ**
   Coelho et al. (2025) は、誤差を直接出すのではなく**温度場・熱流束場を予測**し、
   下流で誤差に変換するモジュール構成を提案。**本研究の「温度場を推定してから変位へ」と同じ思想**。
2. **センサ配置が独立した研究テーマになってきた**
   Teshima (2024)、Ando (2025)、Bünger (2023) がいずれも配置問題を扱う。
   **本研究もこの流れに乗っている**。
3. **測れない境界条件を「推定する」方向**
   CHTCソフトセンサ (2025) は、測れないCHTCを最適化で推定。
   **本研究がEnKFで $Q,h$ を推定するのと動機が同じ**。
4. **測定側の革新**
   カメラ（Li 2024）のように、**安く速く測る**方法が出てきている。
   測定が安くなれば、**データ同化に使える観測が増える**。

**→ 本研究の位置づけ**：①場を推定する ②配置を設計する ③境界条件を推定する、
という**3つの潮流の交点**にあり、さらに**OSSで全部公開**している点が独自。

---

## 6. 最も近い先行研究（質疑への備え）

### 6-1. Teshima, Y., Tanaka, S., Kizaki, T., Sugita, N. (2024)

"Sensor placement strategy based on reduced-order models for thermal error estimation in machine tools",
*CIRP Journal of Manufacturing Science and Technology* **55**, 403–410
🔗 <https://doi.org/10.1016/j.cirpj.2024.10.015> ※**有料**（アブストラクトは無料）

- 温度入力とTCP変位を結ぶ**縮約モデル（ROM）**を作り、**回帰係数を各点の温度感度**とみなす
- 冒頭に「thermal errors in machine tools account for **up to 70 %** of machining errors」
- 彼ら自身が挙げるギャップ：
  「入力点を増やせば精度は上がるが、**最適な点数と配置を決める方法がこれまで無かった**」

### 6-2. Ando, S., Tanaka, S., Teshima, Y., Morishita, J., Kizaki, T.（続報）

"Strategy for Sensor Placement to Estimate Thermal Errors Using **Temperature-Sensitivity
Distribution** Based on a Reduced-Order Model of Machine Tools"
**ICTIMT2025**（4th International Conference on Thermal Issues in Machine Tools）論文集
🔗 <https://doi.org/10.1007/978-3-032-01194-7_31> ※**有料**
研究室: 東京大学 先端加工学（木崎研） <https://mfg.t.u-tokyo.ac.jp/>

- **伝達関数行列**を定義し、TCP-ワーク間の相対変位をセンサ位置の温度変化と結ぶ
- ROMに基づく**センサ感度関数**を定義して配置を決める

### 6-3. 本研究との違い

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

## 7. 発表での使い方

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
| 各社の実装内容 | 各社公開ページ | ✅ 無料 | ⚠️ メーカー公称値 |
