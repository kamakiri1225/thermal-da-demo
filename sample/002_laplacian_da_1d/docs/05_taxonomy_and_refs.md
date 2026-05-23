# データ同化手法の分類と参考文献

## 1. 本プログラムの手法の位置づけ

### 1.1 DA 手法の全体マップ

```
データ同化（Data Assimilation）
│
├── 変分法（Variational DA）
│   ├── 3DVar（3次元変分法）
│   │   全タイムウィンドウの観測を一度に使ってコスト関数を最小化する
│   └── 4DVar（4次元変分法）
│       時間方向にも最適化する（気象予報の標準手法）
│
└── 逐次法（Sequential DA） ← 本プログラムの分類
    観測が来るたびに逐次的に状態を更新する
    │
    ├── Optimal Interpolation（OI）← 最もシンプル
    │   固定の背景誤差共分散 B を使う（P を時間発展させない）
    │   ソルバの状態遷移行列 A_d が不要
    │
    ├── Kalman Filter（KF）
    │   誤差共分散 P を時間発展させる（P = A P A^T + Q）
    │   線形系に対して最適（MMSE 推定量）
    │   │
    │   ├── Extended KF（EKF）
    │   │   非線形 f(x) のヤコビアン F = ∂f/∂x で線形化して KF を適用
    │   │
    │   ├── Unscented KF（UKF）
    │   │   シグマ点でサンプリングして非線形を扱う（ヤコビアン不要）
    │   │
    │   └── Ensemble KF（EnKF）
    │       アンサンブルのばらつきから P を推定（A_d 不要・大規模に強い）
    │       OpenFOAM 連成の実用手法として多く使われる
    │
    └── Particle Filter（PF）
        任意の非線形・非ガウス系を扱える（計算コスト大）
```

---

### 1.2 本プログラム（da_main.py）の手法

**Optimal Interpolation（OI, 最適内挿）**

| 特徴 | 説明 |
|---|---|
| 更新方式 | 逐次（観測が来るたびに更新） |
| 背景誤差共分散 B | 固定（距離相関から事前に設定） |
| ゲイン K | `K = B H^T (H B H^T + R)^-1` |
| 最適性 | 固定した B, R と線形観測 H の仮定のもとで、解析誤差分散を小さくする線形推定 |
| A_d の要否 | 不要 |
| 実装難度 | 低（OpenFOAMをブラックボックス予測器として扱いやすい） |

---

### 1.3 OI との違い

OI（Optimal Interpolation）は KF を「P の時間発展なし」に単純化したものである。

| 項目 | OI（本プログラム） | KF |
|---|---|---|
| 背景誤差共分散 | 固定 B（事前設定） | 時変 P（毎ステップ伝播） |
| 状態遷移行列 A_d | 不要 | 必要 |
| ゲイン K | 固定（または B から毎回計算） | P から毎ステップ計算 |
| 実装の難しさ | 簡単 | やや複雑 |
| 理論的最適性 | 近似的 | 線形系で厳密に最適 |
| 大規模適用 | OpenFOAMをブラックボックス扱いしやすい | 小〜中規模または低次元化モデルに向く |

現在の `da_main.py` は、Python側で `A_d` を作らず、固定した `B` を使うため **OI の実装**である。

---

### 1.4 EnKF との違いと次のステップ

EnKF（Ensemble Kalman Filter）は OpenFOAM などの大規模 CFD ソルバと  
組み合わせる際の実用的な標準手法である。

**KF（本プログラム）の限界**:
- A_d を Python 側で近似（OpenFOAM の厳密な線形化ではない）
- 状態次元が大きくなると P 行列が巨大になる（セル数 N なら N×N 行列）

**EnKF の利点**:
- A_d が不要（アンサンブルのばらつきから P を推定）
- 大規模メッシュにスケールしやすい
- 非線形モデルにも（近似的に）対応できる

---

## 2. データ同化の関連文献

### 2.1 基礎理論・入門書

**日本語**

| # | 書籍 | 著者 | 出版 | 特徴 |
|---|---|---|---|---|
| J1 | 『データ同化入門』 | 中野慎也ほか | 朝倉書店, 2015 | 日本語唯一の体系的入門書。気象分野中心だが KF・EnKF を丁寧に解説 |
| J2 | 『カルマンフィルタの基礎』 | 片山徹 | 東京電機大学出版局, 2000 | 制御工学視点の KF 入門。証明付きで厳密 |
| J3 | 『カルマンフィルタ：Pythonによる実践入門』 | Roger Labbe | 和訳 Web 版あり | Python コード付き。直感的な説明が多い |

**英語**

| # | 書籍 | 著者 | 出版 | 特徴 |
|---|---|---|---|---|
| E1 | *Data Assimilation: The Ensemble Kalman Filter* | Evensen, G. | Springer, 2009 (2nd) | EnKF の原著者による標準テキスト |
| E2 | *Data Assimilation: Methods, Algorithms, and Applications* | Asch, Bocquet, Nodet | SIAM, 2016 | 変分法・KF・EnKF を網羅。無料 PDF あり |
| E3 | *Applied Optimal Estimation* | Gelb (Ed.) | MIT Press, 1974 | KF の古典的テキスト。航空宇宙分野向け |
| E4 | *Introduction to the Kalman Filter* | Welch & Bishop | UNC TR-95-041, 2006 | 無料 PDF。10 ページで KF を完結に説明 |

---

### 2.2 カルマンフィルタ原著論文

| # | 論文 | 発表年 | 収録誌 |
|---|---|---|---|
| P1 | Kalman, R.E., "A New Approach to Linear Filtering and Prediction Problems" | 1960 | *Journal of Basic Engineering*, 82(D), 35–45 |
| P2 | Kalman, R.E. and Bucy, R., "New Results in Linear Filtering and Prediction Theory" | 1961 | *Journal of Basic Engineering*, 83(D), 95–108 |

---

### 2.3 CFD・熱解析 + データ同化

| # | 論文 | 発表年 | 収録誌 | 内容 |
|---|---|---|---|---|
| C1 | Meldi & Poux, "A reduced order model based on Kalman filtering for sequential data assimilation of turbulent flows" | 2017 | *Journal of Computational Physics*, 347, 207–234 | CFD（流体）と KF の組み合わせ。低次元モデルで P を近似する本プログラムと同様の発想 |
| C2 | Foures et al., "A data-assimilation method for Reynolds-averaged Navier–Stokes-driven mean flow reconstruction" | 2014 | *Journal of Fluid Mechanics*, 759, 404–431 | RANS + 変分 DA。CFD 流れ場再構成 |
| C3 | Mons et al., "Reconstruction of unsteady viscous flows using data assimilation schemes" | 2016 | *Journal of Computational Physics*, 316, 255–280 | 非定常流れへの DA 適用 |

---

### 2.4 工作機械 熱誤差補償 + データ同化（本プロジェクト直接関連）

工作機械の温度問題には、データ同化・状態推定を適用できる可能性が高い。

理由は、工作機械の熱変位問題が次の構造を持つためである。

```text
熱源・境界条件・切削液・環境温度などにより温度場が変化する
        ↓
機械構造が熱変形する
        ↓
工具中心点（TCP）や主軸位置に熱変位が生じる
        ↓
加工誤差になる
```

一方で、全ての場所の温度や変位を実測することは難しい。

そこで、

```text
少数の温度センサ
物理モデルまたは低次元モデル
データ同化・状態推定
```

を組み合わせ、未観測の温度場や熱変位を推定する研究が行われている。

工作機械分野では、必ずしも「データ同化」という名前で呼ばれないことが多い。
文献検索では、次のキーワードの方が見つけやすい。

```text
machine tool thermal error compensation
machine tool thermal state observer
Kalman filter thermal error machine tool
digital twin thermal compensation machine tool
sensor placement thermal error compensation
```

| # | 論文 | 発表年 | 収録誌 | 内容 |
|---|---|---|---|---|
| M1 | [Lang et al., "Kalman Filter-Driven State Observer for Thermal Error Compensation in Machine Tool Digital Twins"](https://www.sciencedirect.com/science/article/pii/S2213846324000877) | 2024 | *Manufacturing Letters*, 41, 208–218 | **本プログラムの直接的な参考論文**。KF を使った工作機械熱誤差のオブザーバ設計 |
| M2 | [Lang et al., "Sensor placement utilizing a digital twin for thermal error compensation of machine tools"](https://www.sciencedirect.com/science/article/pii/S0278612525000640) | 2025 | *Journal of Manufacturing Systems*, 80, 243–257 | デジタルツインを使った温度センサ配置最適化。少数センサで熱誤差低減を狙う |
| M3 | [Mayr et al., "Thermal Issues in Machine Tools"](https://www.nist.gov/publications/thermal-issues-machine-tools) | 2012 | *CIRP Annals*, 61(2), 771–791 | 工作機械熱変位のサーベイ論文（背景知識） |
| M4 | [Gomez-Acedo et al., "Methodology for the design of a thermal distortion compensation for large machine tools based in state-space representation with Kalman filter"](https://www.sciencedirect.com/science/article/pii/S0890695513001454) | 2013 | *International Journal of Machine Tools and Manufacture*, 75, 100–108 | 状態空間モデルとKFによる大型工作機械の熱変位補償 |
| M5 | [Brecher et al., "Application of an Unscented Kalman Filter for Modeling Multiple Types of Machine Tool Errors"](https://www.sciencedirect.com/science/article/pii/S2212827117306546) | 2017 | *Procedia CIRP*, 63, 449–454 | UKFを用いて工作機械の幾何誤差・熱誤差モデルを扱う |
| M6 | 安藤ほか, "工作機械の低次元モデルに基づく熱変位推定のための温度センサ配置戦略" | 2025 | 精密工学会春季大会 | 多点温度計測と低次元モデルによる熱変位推定・センサ配置 |
| M7 | 田中ほか, "切削液の影響を考慮した工作機械の熱変位補償に関する研究" | 2020 | 精密工学会 | 多点温度測定による熱変位補償。切削液外乱を含む実機課題 |

#### 本プロジェクトとの関係

本プロジェクトの `laplacianFoam + OI` は、工作機械全体をすぐに扱うものではない。

しかし、考え方はかなり近い。

| 本プロジェクト | 工作機械熱問題 |
|---|---|
| 熱棒の温度場 | 工作機械構造の温度場 |
| 2点の温度センサ | 主軸・ベッド・コラム等の温度センサ |
| 未観測節点の温度推定 | 測っていない部位の温度場・熱変位推定 |
| OIで温度場を補正 | OI/KF/EnKFで熱状態を補正 |
| DAなしとDAありの比較 | 補償なし加工誤差と補償後誤差の比較 |

次の研究ステップとしては、

```text
丸棒・単純構造で実測温度同化
        ↓
簡単な機械構造部品で温度場推定
        ↓
熱変位センサも加えて温度場と変位を比較
        ↓
工作機械の熱変位補償へ展開
```

が自然である。

---

### 2.5 OpenFOAM + データ同化

| # | 資料 | 内容 |
|---|---|---|
| O1 | Villanueva & Reguly, "Flow Estimation Using EnKF and OpenFOAM" | OpenFOAM EnKF の実装事例 |
| O2 | DAPPER (Python library) | https://github.com/nansencenter/DAPPER  Python による DA ライブラリ（EnKF 含む） |
| O3 | OpenDA | https://www.openda.org/  汎用 DA フレームワーク。OpenFOAM との連成実績あり |

---

## 3. 発展方向のまとめ

```
現在（本プログラム）
  線形 KF + laplacianFoam
  Python 側で A_d を近似
     │
     ▼ より自然な拡張
  OI（Optimal Interpolation）
  A_d を使わず固定 B で補正
  気象・海洋分野で広く実績
     │
     ▼ より大規模・非線形
  EnKF（Ensemble Kalman Filter）
  OpenFOAM を N_ensemble 本走らせて P を推定
  A_d 不要・大規模に対応可能
  → OpenFOAM + EnKF は次の実装ステップとして有力
     │
     ▼ 最も本格的
  4DVar + adjoint
  OpenFOAM の随伴コードが必要（連成最適化）
  最も精度が高いが実装が複雑
```

**短期的に最も自然な次ステップは OI**：  
A_d を作らず、固定 B だけで補正する。現在の P 更新ロジックを
削除するだけで OI に変換できる。実験データへの対応が容易。

---

*本ドキュメントは [`02_theory.md`](02_theory.md)・[`03_program_guide.md`](03_program_guide.md) と合わせて読むこと。*
