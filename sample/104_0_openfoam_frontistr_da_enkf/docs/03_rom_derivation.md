# ROM(縮約モデル)の数学的導出 ― なぜ5点で本物に一致するのか

`02_beginner_guide.md` の §3 を、雰囲気ではなく**数式で**きちんと展開する。
「なぜ2万セルの偏微分方程式が、たった5個の常微分方程式で近似できるのか」
「校正とは何を最小化しているのか」を順に示す。

---

## 1. 出発点: 本物が解いている偏微分方程式(PDE)

OpenFOAM(chtMultiRegionFoam)は、固体内の温度 $T(\mathbf{x},t)$ について
**熱伝導方程式**を解いている:

$$\rho c_p \frac{\partial T}{\partial t} = \nabla\!\cdot\!\left(k\,\nabla T\right) + \dot{q}$$

- $\rho$: 密度 [kg/m³]、$c_p$: 比熱 [J/(kg·K)]、$k$: 熱伝導率 [W/(m·K)]
- $\dot q$: 単位体積あたり発熱 [W/m³]

境界条件は、ヒータ面で熱流束 $q_\mathrm{heater}$、その他の面で周囲空気への放熱

$$-k\,\nabla T\cdot \mathbf{n} = h_s\left(T - T_\mathrm{air}\right) \quad(\text{放熱面})$$

これを約2万個の有限体積セルで離散化するので、**2万本の連立ODE**になる。正確だが重い。

---

## 2. 縮約の原理: 有限体積法を「粗くする」極限

有限体積法は、各セル $i$ の体積 $V_i$ で PDE を積分し、発散定理を使う:

$$\int_{V_i} \rho c_p \frac{\partial T}{\partial t}\,dV
= \int_{V_i} \nabla\!\cdot\!(k\nabla T)\,dV + \int_{V_i}\dot q\,dV
= \underbrace{\sum_{f\in\partial V_i}\!\! \big(k\nabla T\cdot\mathbf{n}\big)_f A_f}_{\text{面を通る熱流}} + \dot Q_i$$

セル平均温度 $T_i \equiv \frac{1}{V_i}\int_{V_i}T\,dV$ を使うと、左辺は $\rho c_p V_i \dfrac{dT_i}{dt}$。
隣接セル間の熱流を**温度差に比例**すると近似(これが離散化の本質):

$$\big(k\nabla T\cdot \mathbf n\big)_f A_f \;\approx\; \frac{k A_f}{d_{ij}}\,(T_j - T_i) \;\equiv\; K_{ij}\,(T_j - T_i)$$

ここで $K_{ij} = \dfrac{k A_f}{d_{ij}}$ は面 $f$ の**熱コンダクタンス** [W/K]
($A_f$: 面積、$d_{ij}$: セル中心間距離)。まとめると、**任意の粒度**で同じ形の式が立つ:

$$\boxed{\,C_i \frac{dT_i}{dt} = \sum_{j} K_{ij}\,(T_j - T_i) - h_i\,(T_i - T_\mathrm{air}) + \dot Q_i\,}\tag{★}$$

- $C_i = \rho c_p V_i$: ノード $i$ の**熱容量** [J/K]
- $K_{ij}$: ノード間コンダクタンス [W/K]、$h_i$: 放熱コンダクタンス [W/K]
- $\dot Q_i$: ノードへの発熱 [W]

**縮約モデル(ROM)とは、この式(★)をセル2万個ではなく、代表5ノードで書いたもの**。
つまり ROM は「有限体積離散化を極端に粗くした同じ物理式」であり、別物の近似ではない。
これが「5点で本物に一致しうる」数学的な理由。

---

## 3. 5ノードの選び方とコンダクタンス構造

102_0 のプローブ位置に対応する5ノードを取る:

| index | 名前 | 位置 | 役割 |
|-------|------|------|------|
| 0 | hot | +X(ヒータ直下) | 発熱の入口 |
| 1 | mid | +Y(周方向90°) | 側面中間 |
| 2 | cold | −X(反ヒータ) | 一番冷たい |
| 3 | top | 上面近傍 | 軸方向 |
| 4 | core | 肉厚中心(潜在) | 蓄熱の芯 |

熱の通り道(コンダクタンス $K_{ij}$)を物理的なつながりで与える(`conductance_matrix`):

- 周方向: hot–mid, mid–cold(係数 $k_\mathrm{circ}$)
- 径方向: hot/mid/cold/top ↔ core(係数 $k_\mathrm{core}$)
- 軸方向: hot–top(係数 $k_\mathrm{axial}$)

対称行列 $K$($K_{ij}=K_{ji}$)として書くと、校正で決めるのは
$\{C_0,\dots,C_4,\ k_\mathrm{circ},\ k_\mathrm{core},\ k_\mathrm{axial},\ h\}$ の**9個の実数**だけ。

---

## 4. 行列形と厳密解の構造

式(★)を全ノードまとめてベクトル $\mathbf{T}=[T_0,\dots,T_4]^\top$ で書く:

$$\frac{d\mathbf T}{dt} = A\,\mathbf T + \mathbf b(t)$$

システム行列 $A$ の各成分は(熱容量で正規化):

$$A_{ij} = \frac{K_{ij}}{C_i}\ (i\neq j),\qquad
A_{ii} = -\frac{1}{C_i}\Big(\sum_{j} K_{ij} + h\Big)$$

定数項は放熱とヒータ:

$$\mathbf b(t) = \underbrace{\frac{h\,T_\mathrm{air}}{C_i}}_{\text{放熱}} + \underbrace{\frac{\dot Q_i(t)}{C_i}}_{\text{ヒータは hot ノードのみ}}$$

$A$ は(対角優位で固有値が負の)安定な行列。ヒータ一定区間では $\mathbf b$ が定数なので
厳密解は

$$\mathbf T(t) = e^{A(t-t_0)}\mathbf T(t_0) + A^{-1}\!\left(e^{A(t-t_0)}-I\right)\mathbf b$$

実装では一般性のため RK4 で数値積分する(`integrate_single` / `integrate_ensemble`)。
コード対応: `system_matrix()` が $A,\ \mathbf b_\mathrm{air}$、`heater_input()` が $\dot Q_i/C_i$。

---

## 5. 校正 = 逆問題(最小二乗)

未知の9パラメータ $\boldsymbol\theta=\{C_i, k_\mathrm{circ}, k_\mathrm{core}, k_\mathrm{axial}, h\}$ を、
**本物(102_0)の温度履歴に一致するように**決める。これは最小二乗の逆問題:

$$\hat{\boldsymbol\theta}
= \arg\min_{\boldsymbol\theta}\ \sum_{n=1}^{N_t}\ \sum_{p\in\{\mathrm{hot,mid,cold,top}\}}
\Big( T_p^\mathrm{ROM}(t_n;\boldsymbol\theta) - T_p^\mathrm{OF}(t_n) \Big)^2$$

- $T_p^\mathrm{OF}(t_n)$: OpenFOAM の観測4点の温度(`temperature_history.csv`、$N_t=121$ 時刻)
- $T_p^\mathrm{ROM}(t_n;\boldsymbol\theta)$: 式(★)を真のヒータ($q\!=\!1$)で積分した予測
- core は測っていない潜在変数なので目的関数に入れない(残り4点で拘束)

制約 $\boldsymbol\theta>0$(すべて正)のもとで
`scipy.optimize.least_squares`(信頼領域反射法)で解く(`calibrate.py`)。
残差ベクトルは $484 = 121\times4$ 本、パラメータ9個の**過剰決定**問題。

### 校正結果(実際の値)

$$C = [\,206.8,\ 10.0,\ 473.4,\ 131.5,\ 369.3\,]\ \mathrm{J/K}$$
$$k_\mathrm{circ}=1.925,\quad k_\mathrm{core}=3.288,\quad k_\mathrm{axial}=2.324\ \mathrm{W/K}$$
$$h = 0.0236\ \mathrm{W/K}$$

当てはめ残差:

$$\mathrm{RMSE} = \sqrt{\frac{1}{484}\sum_{n,p}\big(T^\mathrm{ROM}-T^\mathrm{OF}\big)^2} = 0.0119\ \mathrm{K}$$

**物理的な妥当性チェック**: 全熱容量 $\sum_i C_i \approx 1191\ \mathrm{J/K}$。
試験体は鋼 2.49 kg、$c_p\approx480\ \mathrm{J/(kg\,K)}$ なので
$\rho c_p V \approx 2.49\times480 \approx 1195\ \mathrm{J/K}$ ― ほぼ一致する。
つまり校正は数字合わせでなく、**物理的に正しい熱容量**を復元している。
また $h$ が非常に小さい($0.024\ \mathrm{W/K}$)ことは、600秒では放熱が遅く、
熱が主に内部再分配される(300秒後も約23.6℃で平衡)という実挙動と整合する。

---

## 6. なぜ0.01℃で一致するのか(次元削減の理論的背景)

一般に、拡散方程式の解は固有モード展開できる:

$$T(\mathbf x,t) = T_\infty(\mathbf x) + \sum_{m=1}^{\infty} a_m\,\phi_m(\mathbf x)\,e^{-\lambda_m t}$$

固有値 $\lambda_m$ は $m$ とともに急速に大きくなり($\lambda_m \sim m^2$)、
高次モードは**すぐ減衰する**。したがって、ゆっくりした緩やかな加熱
(ヒータのオンオフは秒〜分スケール)では、**低次の数モードだけ**で解がほぼ再現できる。
5ノードROMは、この「効いている低次モード」を張るのに十分な自由度を持つ。
これが「2万自由度→5自由度」でも 0.01℃ で一致する理論的な理由
(縮約基底法・モード打ち切りの考え方)。

---

## 7. データ同化での使われ方

校正で $\{C, K, h\}$ は既知として固定し、**データ同化で推定するのは
状態 $\mathbf T$ と、未知パラメータ $q$(ヒータ倍率)・$h$(放熱)** だけにする
(拡大状態、`02_beginner_guide.md` §8)。ROM は式(★)を数千ステップ積分しても
一瞬なので、アンサンブル60本×30サイクルでも数秒で終わる。

要するに ROM は「**同じ物理式(★)を、実データで係数校正した粗い離散化**」であり、
- 前進が速い(データ同化を何度でも試せる)
- 挙動が本物と一致(校正 RMSE 0.012 K)
という2条件を満たす、理論的に正当な代用モデルである。

---

### コード対応

| 数式 | 実装(`dacore/`) |
|------|------------------|
| 式(★)/ 行列 $A,\mathbf b$ | `cht_rom.py: system_matrix(), heater_input()` |
| $K_{ij}$ 構造 | `cht_rom.py: conductance_matrix()` |
| 時間積分(RK4) | `cht_rom.py: integrate_single(), integrate_ensemble()` |
| 最小二乗校正(§5) | `calibrate.py: calibrate()` |
| 校正結果 | `config/rom_calibrated.yaml` |
