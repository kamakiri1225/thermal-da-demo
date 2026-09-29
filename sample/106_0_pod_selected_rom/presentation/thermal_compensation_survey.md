# 熱変位補正（thermal error compensation）の実務調査

作成日: 2026-09-29 ／ オープンCAE2026 発表の背景調査

**本メモの方針**：主張には**無料で読める一次情報のリンク**を必ず付ける。
有料論文しか根拠が無い場合は、そのことを明記する。

---

## 0. 結論（先に答え）

| 質問 | 答え |
|---|---|
| 「実用の主役は補償」＝熱変位補正か？ | **はい**。文献では対策を **avoidance（回避＝設計・冷却）** と **compensation（補償＝熱変位補正）** の2つに分け、**compensationの方が安価で実用的**とされる |
| 熱変位補正は何をしている？ | 機械各所の**温度を測り**、事前に取った**温度→刃先変位の関係モデル**で変位を予測し、**その逆符号を軸の指令位置に足し込む** |
| ご認識（多点温度＋変位実測→回帰でモデル化）は正しいか | **正しい**。文献でも「予測モデルを作り、それを逆に使って補正量を決める」と明記 |
| いつ補正する？ | **運転中に連続して**。CNCが温度を取り込み、**座標系オフセットを常時更新**する（加工の合間ではなく常時） |
| 機械に埋め込まれているか | **はい**。オークマOSP、マザックのAI Thermal Shield、DMG MORIのSGSなど、**CNC側の標準／オプション機能**として実装 |

---

## 1. 一番の根拠（オープンアクセス・全文無料）

**Li, Z., Vogl, G. W., Kinzel, E. C., Santa, B., Landers, R. G. (2024)**
"Machine Tool Thermal Error Measurement and Prediction via Wireless Microscope",
*Manufacturing Letters* **41**, 1440–1451（NAMRC 52）。**CC BY-NC-ND のオープンアクセス**
📄 無料PDF: <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=957076>
（NIST共著のため NIST サイトで全文公開）

この論文の Introduction に、知りたいことがほぼ全部書いてあります（原文引用）:

> "Thermal errors can contribute **up to 75 percent** of the overall machining errors of a machined part."

> "In general, there are two methods for thermal error reduction:
> **thermal error avoidance** and **thermal error compensation**."

> "Thermal error avoidance techniques make the machine tool less sensitive to temperature variations.
> These techniques include the use of materials with reduced friction and low coefficients of
> thermal expansion, optimization of lubrication and cooling systems, etc."

> "In contrast, **thermal error compensation is normally based on a predictive model established by
> the thermal error measurement of a machine tool**."

> "**Thermal error avoidance is typically a more costly solution than thermal error compensation**
> and is more sensitive to modeling errors and unknown disturbances."

> "Thermal error compensation strategies typically employ a **predictive error model, which are
> inverted to determine compensation amounts**. Common thermal error models include
> **least-square regression, finite element, neural network, gray system**, etc."

**→ ご認識どおり**です。「多点で温度を測り、刃先変位も測り、その関係を回帰などでモデル化し、
そのモデルを逆に使って補正量を出す」が、そのまま文献の記述と一致します。
線形回帰・FEM・ニューラルネット・グレーシステムが代表的なモデルとして挙がっています
（ガウス過程回帰もこの系列の一つ）。

---

## 2. 「熱誤差は何 % か」の数字（訂正あり）

**結論：引用するなら「最大75 %」が安全**です。

| 数字 | 出典 | 検証状況 |
|---|---|---|
| **75 %** | Mayr et al. (2012) CIRP Annals を、**無料で読める2本が引用**している | ✅ **裏取り済み** |
| 40〜70 % | Bryan (1990) CIRP Annals 39(2) 645–656 | ⚠️ **原典未確認**（有料。検索結果と二次引用のみ） |

**無料で確認できる引用元（2本とも全文無料）**

1. **Bünger, A., Herzog, R., Naumann, A., Stoll, M. (2023)**
   "Uncertainty Propagation of Initial Conditions in Thermal Models", arXiv:2306.12736
   📄 <https://arxiv.org/abs/2306.12736>
   > "According to **Mayr et al., 2012, the thermal error accounts for 75 % of the total
   > manufacturing error in the final product**."

2. **Li et al. (2024)**（上記）
   > "Thermal errors can contribute **up to 75 percent** of the overall machining errors"

**スライドでの書き方（推奨）**：
「熱変位は加工誤差の**最大75 %**を占める［Mayr et al. 2012］」
※前回お渡しした「40〜70 %［Bryan 1990］」は**原典未確認**でしたので、
上記に差し替えることをお勧めします。

---

## 3. 「温度が合っても変形が合うとは限らない」の根拠

こちらも**無料で読める**根拠に差し替えます。

**Bünger et al. (2023) arXiv:2306.12736**（全文無料）
📄 <https://arxiv.org/abs/2306.12736>
> "the accuracy of the **TCP prediction depends highly on the accuracy of the model parameters,
> such as heat exchange parameters, and the initial temperature**"

→ **熱交換パラメータ（＝熱伝達率）と初期温度**の不確かさが、TCP予測精度を直接左右する、と明記。

### 「熱伝達率が直接測れない」の出典（ご指摘の"複数の論文"）

前回「複数の論文で指摘」と書きましたが、**具体名を挙げていませんでした**。確認できたものを挙げます。

1. **Intelligent Soft Sensor for Spindle Convective Heat Transfer Coefficient Under Varying
   Operating Conditions Using Improved Grey Wolf Optimization Algorithm** (2025)
   📄 無料全文（PMC）: <https://pmc.ncbi.nlm.nih.gov/articles/PMC12473924/>
   > 「CHTCを直接測定する専用計測器が存在しないため、主軸の熱解析は大きな困難に直面する」
   > 「従来の熱流束センサや非接触赤外サーモグラフィは温度分布は得られるが、**CHTCを直接測ることはできず**、
   > 熱モデルの精度を制限している」の趣旨の記述あり
   ※私の環境からは本文取得がブロックされたため、**検索結果の引用に基づく**。
   リンクは無料公開なので**ご自身で原文確認を推奨**。

2. **The Thermal Error Estimation of the Machine Tool Spindle Based on Machine Learning**,
   *Machines* **9**(9), 184 (2021)。MDPI **オープンアクセス**
   📄 <https://www.mdpi.com/2075-1702/9/9/184>
   > 「FEMの精度は、熱源・熱伝達係数・境界条件が明確に定義されているかに依存する。
   > しかし主軸は材質の異なる多数の部品からなり、熱源と境界条件は組立や加工条件に強く依存するため、
   > **汎用的に通用する熱源・熱伝達係数・境界条件を定義するのは難しい**」の趣旨

**→ 発表では②（MDPI、確実にオープンアクセス）を主に引くのが安全**です。

---

## 4. メーカー各社の実装

⚠️ 以下は**各社の公開情報（製品ページ・技術紹介）に基づく**もので、
学術論文による検証ではありません。数値は各社の公称値です。

### オークマ ― サーモフレンドリーコンセプト（TFC）

🔗 <https://www.okuma.co.jp/onlyone/thermo/>

**考え方が特徴的**：温度変化を「抑え込む」のではなく **"受け入れる"**。

| 要素 | 内容 |
|---|---|
| **設計側** | 「熱変形の単純化構造」「温度分布均一化」により、**機械を素直に変形させ**、ねじれ・傾きを抑えて**熱変位を予測可能な状態にする** |
| **TAS-C**（環境熱変位制御） | 機械の熱変位特性を踏まえ、**適切に配置されたセンサの温度情報と送り軸の位置情報**から、環境温度変化による構造体の熱変位を**推定して制御** |
| **TAS-S**（主軸熱変位制御） | 主軸の温度変化を**リアルタイム管理**し、回転・停止時の熱変位を制御 |
| 実績 | 2019年10月時点で全127機種中83機種に搭載、出荷5万台超 |

**本研究との関係**：TAS-Cの「センサ温度＋軸位置から熱変位を推定」は、
まさに本研究の「**少数観測から状態を推定**」と同じ構図。
ただしオークマは**設計で変形を単純化してから**推定する点が巧み。

### ヤマザキマザック ― AI Thermal Shield / Intelligent Thermal Shield

🔗 <https://www.mazak.com/jp-ja/technology/accuracy/>

| 要素 | 内容 |
|---|---|
| 入力 | **主軸回転速度**＋**温度センサ情報** |
| 考慮する条件 | 温度変化、**機械位置**、**クーラントON/OFF** |
| 学習 | 加工後の計測データを**蓄積・学習**し、顧客の加工環境に適した補正に最適化 |
| 公称性能 | **室温が8 ℃変化しても連続加工精度6 µmを維持** |
| ハード側 | 主軸軸受・モータ外筒に温調冷却油、**X/Y/Z全軸のボールねじ軸冷却を標準装備** |

**注目点**：**クーラントON/OFFを入力に入れている**。
（前スライドのクーラント温度の話と直結。冷却液は熱境界条件そのもの）

### DMG MORI ― Spindle Growth Sensor（SGS）ほか

🔗 <https://en.dmgmori.com/products/machines/milling/vertical-milling/nvx/nvx-5100>

| 要素 | 内容 |
|---|---|
| **SGS** | 主軸の**軸方向伸びを直接センシング**して補正（モデル予測ではなく実測フィードバック） |
| 構造側 | 冷却されたリニアガイド、主軸成長のアクティブ制御 |
| 制御 | 機械各所の温度センサで**リアルタイム監視**し、**工具オフセットを自動調整** |

**注目点**：**変位を直接測る**アプローチ。
本研究の「変位観測を足すと精度が上がる」（温度2点0.197 K → ＋変位0.159 K）と
**同じ思想**で、実機で既に商用化されている。

### FANUC ― CNC側の機能

🔗 <https://www.fanuc.co.jp/ja/product/cnc/index.html>

- **AI熱変位補正**（2022年6月発表）
- 多点の温度センサ入力に対応する小型・省配線のI/Oユニットを提供

**→ 補正は最終的にCNCの中で座標オフセットとして効く**、という実装の裏付け。

---

## 5. 「いつ補正するのか」への答え

文献とメーカー情報を総合すると:

1. **事前（オフライン）**：さまざまな温度条件で運転し、
   **多点温度**と**刃先変位**を同時計測 → 回帰などで**予測モデル**を作る
2. **運転中（オンライン）**：CNCが温度センサを**常時読み込み**、
   モデルで変位を予測 → **符号を逆にして軸の指令位置に加算**
   （オークマTAS-Sは「リアルタイム管理」、DMG MORIは「real-time／自動調整」と明記）
3. **随時（学習）**：マザックは加工後の計測データを蓄積・学習して補正を更新

**つまり「加工の合間に測って直す」のではなく、運転中ずっと補正がかかり続けている**のが標準です。

---

## 6. 最近接の先行研究（東京大学・木崎研）

### (a) Teshima, Y., Tanaka, S., Kizaki, T., Sugita, N. (2024)

"Sensor placement strategy based on reduced-order models for thermal error estimation in machine tools",
*CIRP Journal of Manufacturing Science and Technology* **55**, 403–410
🔗 <https://doi.org/10.1016/j.cirpj.2024.10.015>（有料。アブストラクトは無料で読めます）

- 温度入力とTCP変位を結ぶ**縮約モデル（ROM）**を作り、**回帰係数を各点の温度感度**とみなす
- 彼ら自身が挙げるギャップ：
  「入力点を増やせば精度は上がるが、**最適な点数と配置を決める方法がこれまで無かった**」

### (b) Ando, S., Tanaka, S., Teshima, Y., Morishita, J., Kizaki, T.（続報）

"Strategy for Sensor Placement to Estimate Thermal Errors Using **Temperature-Sensitivity
Distribution** Based on a Reduced-Order Model of Machine Tools"
**ICTIMT2025**（4th International Conference on Thermal Issues in Machine Tools）論文集
🔗 <https://doi.org/10.1007/978-3-032-01194-7_31>（有料）

- **伝達関数行列**を定義し、TCP-ワーク間の相対変位をセンサ位置の温度変化と結ぶ
- ROMに基づく**センサ感度関数**を定義して配置を決める
- 研究室ページ: <https://mfg.t.u-tokyo.ac.jp/>

### 本研究との違い（発表・質疑用）

| | Teshima / Ando ら | 本研究 |
|---|---|---|
| ROMの作り方 | 温度→TCP変位の**伝達関数・回帰** | **POD＋Q-DEIM**で代表点を選び集中定数ROMを校正 |
| 配置の基準 | **温度感度分布**（伝達関数から） | **観測の価値** $\Delta=\mathrm{Cov}(X,y)^2/(\mathrm{Var}(y)+r)$ ＝相互情報量。**推定対象ごとに変わる** |
| 未知パラメータ | 扱わない（温度は入力） | **発熱量 $Q$ ・放熱 $h$ を状態に入れて同時推定**（EnKF） |
| 観測の種類 | 温度センサ | **温度＋変位**を同じ枠組みで扱う |
| 実装 | ― | **OSSのみ**（OpenFOAM＋FrontISTR＋Python）で全公開 |

**→ 「ROMでセンサ配置を決める」流れは既にある。
本研究の差分は ①データ同化で境界条件ごと推定する ②温度と変位を同じ価値尺度で比べる
③OSS一気通貫で再現可能にした、の3点。**

---

## 7. 発表での使い方（提案）

導入スライドの3行目を、**この調査を踏まえて強化**できます:

> - 熱変位は加工誤差の**最大75 %**を占める［Mayr et al. 2012］
> - 対策は**回避（設計・冷却）**と**補償（熱変位補正）**の2つ。
>   **補償の方が安価**で、実機（オークマ・マザック・DMG MORI）に広く実装されている
> - しかし補償の前提となる予測が難しい。**熱伝達率は直接測る計測器が存在せず**、
>   境界条件の不確かさがTCP予測精度を直接左右する［Bünger et al. 2023］

**→ そして本研究：「境界条件そのものを観測から推定してしまえばよい」**
という流れに自然につながります。

---

## 付録：検証状況の一覧

| 主張 | 出典 | 無料で読めるか | 検証 |
|---|---|---|---|
| 熱誤差は最大75 % | Li et al. 2024 / Bünger et al. 2023 | ✅ 両方とも全文無料 | ✅ 原文確認済み |
| avoidance と compensation の2分類 | Li et al. 2024 | ✅ 無料 | ✅ 原文確認済み |
| 補正は予測モデルを逆に使う | Li et al. 2024 | ✅ 無料 | ✅ 原文確認済み |
| モデルは回帰・FEM・NN・グレー系 | Li et al. 2024 | ✅ 無料 | ✅ 原文確認済み |
| 熱交換パラメータ・初期温度が精度を左右 | Bünger et al. 2023 | ✅ 無料 | ✅ 原文確認済み |
| CHTCを直接測る計測器が無い | PMC12473924 (2025) | ✅ 無料 | ⚠️ 検索結果のみ（要原文確認） |
| FEMの境界条件定義が難しい | Machines 9(9) 184 (2021) | ✅ 無料 | ⚠️ 検索結果のみ（要原文確認） |
| 40〜70 % | Bryan 1990 | ❌ 有料 | ❌ **未確認。使用非推奨** |
| 各社の実装内容 | 各社公開ページ | ✅ 無料 | ⚠️ メーカー公称値 |
