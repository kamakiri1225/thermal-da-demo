# blog_005：推定したい量で「最適なセンサ」は変わる ― 観測の価値を1本の式で設計する

センサ配置の素朴な疑問から始めます。

> **温度センサや変位センサは、どこに置くのが一番いいのか？**

**「何を推定したいか」で有効な観測の位置・種類・時間が異なる可能性があります。** 手計算例でその理由を確認し、実モデルでは候補・ノイズ・時間条件を揃えて検証します。異なる推定対象でも同じ配置が最良になる場合はあり、対象の数だけ必ず別の最適点があるという意味ではありません。
この記事では、それを **1本の式 → 小さな行列の手計算 → 実データ** の順で示します。

> **この記事の流れ**
> 1. 観測の価値を測る式 $\Delta$ を導く（理論）
> 2. 3変数モデルの行列手計算で「対象ごとに最適観測が変わる」を見る（★が動く表）
> 3. Q/h直接感度を微分から導出 → RK4 → 全温度場・変位へ復元 → Fisher情報で識別性を評価
> 4. 実データ（ROM＋EnKF/OI）で検証
>
> 関連記事：OIの仕組みは [blog_002](blog_002_oi_data_assimilation.md)、ROM/PODは [blog_004](blog_004_pod_selected_rom.md)。

---

## 1. 観測の価値を測る「1本の式」

ここでは線形観測・ガウス誤差モデルでの分散減少を考えます。 $x$ と $v,c$ は列ベクトル、 $B$ は状態の事前共分散です。
式中のyはノイズを加える前の観測予報 $y=v^{\mathsf T}x$ 、実際の測定は $y_{\rm obs}=y+\epsilon$ 、 $\mathrm{Var}(\epsilon)=r$ とします。これにより分母でノイズを二重加算しません。
対象Xの更新前後の分散は $c^{\mathsf T}Bc$ と $c^{\mathsf T}B^ac$ なので、差を取ると次のΔが得られます。

データ同化の更新式（blog_002）から出発します。観測 $y$ を1つ入れると、状態の共分散は

$$B^a=B-\frac{(Bh)(Bh)^{\mathsf T}}{v^{\mathsf T}Bh+r}$$

と縮みます（ $v$ は観測行列の1行を転置した列ベクトル、 $r$ は観測ノイズ分散）。いま知りたい量（推定対象）を
$X=c^{\mathsf T}x$ と書くと（ $c$ は対象を取り出すベクトル）、**$X$ の分散がどれだけ減るか**は

$$\boxed{\ \Delta(X,\ y)=\frac{\mathrm{Cov}(X,\ y)^2}{\mathrm{Var}(y)+r}
=\frac{\bigl(c^{\mathsf T}Bh\bigr)^2}{v^{\mathsf T}Bh+r}\ }$$

これが**観測の価値**です。読み方は簡単で、

- 分子 $\mathrm{Cov}(X,y)^2$ ：**その観測が、知りたい量とどれだけ連動しているか**
- 分母 $\mathrm{Var}(y)+r$ ：その観測自身のばらつき＋ノイズ

つまり「**知りたい量とよく連動し、ノイズの小さい観測ほど価値が高い**」。
ここで大事なのは、分子に **対象 $c$ と観測 $v$ の両方**が入っていること。
**観測の価値は対象に依存する**＝全対象に共通する最適センサ位置は一般には保証されない、が式の形から既に見えています。

---

## 2. 小さな行列で手計算 ―「★」が対象ごとに動く

### 2-1. 舞台（blog_002 §3-5と同じ3変数）

状態 $x=(T_1,\ T_2,\ Q)$ 。温度は $Q$ で決まり（感度 $s=(0.3,\ 0.1)$ ）、変位は温度で決まる
（感度 $w=(0.5,\ 0.2)$ ）。背景共分散は $\sigma_T=1,\ \sigma_Q=2$ として

$$B=\begin{pmatrix}
1.36 & 0.12 & 1.2\\
0.12 & 1.04 & 0.4\\
1.2 & 0.4 & 4
\end{pmatrix},\qquad
B=\mathrm{diag}(\sigma_T^2,\ \sigma_T^2,\ 0)+\sigma_Q^2\,a\,a^{\mathsf T},\quad a=(s_1,\ s_2,\ 1)^{\mathsf T}.$$

（成分の出どころ：対角 $\sigma_T^2+\sigma_Q^2 s_i^2$ 、温度と $Q$ の相関 $\sigma_Q^2 s_i$ 。導出は blog_002 §3-5）

この手計算では温度差をK、QをW、変位をµmで表します。 $s_i$ はK/W、 $w_i$ はµm/K、温度観測分散はK²、変位観測分散はµm²です。相関の効果を示す説明用の値であり、実機の物性値ではありません。

**観測の候補は3つ**（それぞれ $v$ と $r$ ）:

| 候補 | $v$ | ノイズ $r$ |
|---|---|---|
| 温度点1 | $(1,\ 0,\ 0)$ | $0.3^2=0.09$ |
| 温度点2 | $(0,\ 1,\ 0)$ | $0.09$ |
| 変位 $u$ | $(w_1,\ w_2,\ 0)=(0.5,\ 0.2,\ 0)$ | $0.1^2=0.01$ |

### 2-2. 部品をそろえる（行列×ベクトルを成分で）

まず各候補の $Bh$ （＝「その観測と各変数の共分散」の列）と $\mathrm{Var}(y)+r$ :

$$B\begin{pmatrix}1\\0\\0\end{pmatrix}=\begin{pmatrix}1.36\\ 0.12\\ 1.2\end{pmatrix},\qquad
B\begin{pmatrix}0\\1\\0\end{pmatrix}=\begin{pmatrix}0.12\\ 1.04\\ 0.4\end{pmatrix},\qquad
B\begin{pmatrix}0.5\\0.2\\0\end{pmatrix}
=\begin{pmatrix}0.5\cdot1.36+0.2\cdot0.12\\ 0.5\cdot0.12+0.2\cdot1.04\\ 0.5\cdot1.2+0.2\cdot0.4\end{pmatrix}
=\begin{pmatrix}0.704\\ 0.268\\ 0.68\end{pmatrix}$$

$$\mathrm{Var}(y_1)+r=1.36+0.09=1.45,\qquad
\mathrm{Var}(y_2)+r=1.04+0.09=1.13,$$

$$\mathrm{Var}(u)+r=w^{\mathsf T}Bw+r=0.5\cdot0.704+0.2\cdot0.268+0.01=0.4156.$$

### 2-3. 価値 $\Delta$ を対象ごとに計算（例を3つ）

**例1：対象＝発熱量 $Q$ （ $c=(0,0,1)$ ）を温度点1で測る**

$$\Delta(Q,\ y_1)=\frac{\bigl(c^{\mathsf T}Bh\bigr)^2}{1.45}=\frac{1.2^2}{1.45}=0.993$$

**例2：対象＝ $Q$ を変位 $u$ で測る**（分子は $Bw$ の第3成分）

$$\Delta(Q,\ u)=\frac{0.68^2}{0.4156}=1.113\ \gt\ 0.993$$

→ **$Q$ には温度点1より変位の方が価値が高い**。

**例3：対象＝ $T_2$ を変位 $u$ で測る**（分子は $Bw$ の第2成分）

$$\Delta(T_2,\ u)=\frac{0.268^2}{0.4156}=0.173\ \ll\ \Delta(T_2,\ y_2)=\frac{1.04^2}{1.13}=0.957$$

→ **$T_2$ には温度点2が最適**（変位は $w_2=0.2$ としか連動しない）。

### 2-4. 全部並べると「★」が行ごとに動く

全温度場の行は $\Delta(T_1,y)+\Delta(T_2,y)$ 、つまり2点の分散減少の和です。温度和 $T_1+T_2$ の分散ではありません。温度の行はK²、Qの行はW²なので、大小は同じ行の中で比較します。

同じ計算を4対象×3候補で並べたのが下の表です（**すべて上の $B$ から電卓で出せます**）:

![観測の価値Δの表（対象×候補）](img/blog005_delta_table.png)

| 推定対象＼観測候補 | 温度点1 | 温度点2 | 変位 $u$ |
|---|---|---|---|
| $T_1$ | **1.276 ★** | 0.013 | 1.193 |
| $T_2$ | 0.010 | **0.957 ★** | 0.173 |
| 発熱量 $Q$ | 0.993 | 0.142 | **1.113 ★** |
| 全温度場（2点の分散の和） | 1.286 | 0.970 | **1.365 ★** |

**★（最適な観測）が行ごとに違います。** これがこの記事の主張の核心で、しかも理由は明快：

- $T_1$ を知りたい → $T_1$ 自身との連動が最大の「温度点1」
- $T_2$ を知りたい → 「温度点2」（変位は $w_2=0.2$ としか連動しない）
- $Q$ を知りたい → **変位**（複数温度の和として $Q$ の効果をまとめて拾い、ノイズも小さい）
- 全温度場 → **変位**（ $T_1$ とも $T_2$ とも連動するので合計で勝つ）

> **一言でいうと**：観測の価値は「対象との共分散」で決まり、共分散は対象ごとに違う。
> だから**対象によって最良の観測候補が変わり得る**。

### 2-5. 図解：なぜ変位が $Q$ ・全場に強いのか

![変位観測が温度を直す仕組み](img/disp_to_temp_mechanism.png)

*変位は複数点の温度誤差の**射影**（左・実線）を1つの数で読み、ハズレが共分散に比例して
**逆流**する（左・破線）。誤差が「共通の型」で起きるとき、型と重み $w$ が揃う高W点は
信号がノイズを超える（右）。詳しい導出は blog_002 §3-7。*

---

## 3. 第3の軸「時間」― いつ観測するかも対象で変わる

### 3-1. ここからの問いと、記号の定義

§1〜2は「共分散Bが与えられたとき、どの観測が有効か」を考えました。ここからは一段戻り、**Qやhを変えると、どの場所の温度・変位がどれだけ変わるかを、熱方程式から計算**します。その応答が直接感度です。

```mermaid
flowchart LR
 A["熱方程式をQ・hで微分"] --> B["5点温度と感度を同時積分"]
 B --> C["PODで全セルの感度へ復元"]
 C --> D["FrontISTR応答で変位感度へ"]
 D --> E["ノイズで規格化して観測点を比較"]
 E --> F["選んだ配置をEnKFで検証"]
```

| 記号 | 大きさ・単位 | 意味 |
|---|---|---|
| $m,N,N_t,r$ | 5、20,696、121、5 | ROM代表点数、全セル数、保存時刻数、採用モード数 |
| $\mathbf T(t)$ | $m\times1$ 、K | ROM代表点の温度。添字はコードに合わせ0〜4 |
| $C$ | $m\times m$ 、J/K | 対角成分が熱容量 $C_i$ の対角行列 |
| $K_{ij},L$ | W/K | ノード間の結合係数、結合をまとめた行列 |
| $h$ | W/K | ROMの各点に共通の放熱コンダクタンス |
| $Q_0,\chi(t)$ | W、無次元 | 加熱中の発熱量、ヒータON/OFF関数 |
| $\mathbf b,\mathbf1,I$ | $m\times1,m\times1,m\times m$ | 発熱配分、全成分1のベクトル、単位行列 |
| $\mathbf S_Q,\mathbf S_h$ | $m\times1$ | Qとhに対する温度感度 |

$T_\infty=293.15$ Kは外気温。**hは壁面熱伝達率[W/(m² K)]ではありません。** 論文で区別するなら $G_{\rm loss}$ と書き、現行コードの`h`に対応させます。

この節の $C$ は熱容量で、EnKFの共分散ではありません。また、§1〜2で観測演算子を表す記号と、ここでの放熱hを区別します。

### 3-2. 1点の熱収支を、行列の熱方程式にする

各点について、

$$
C_i\dot T_i=\sum_{j\ne i}K_{ij}(T_j-T_i)-h(T_i-T_\infty)+b_iQ_0\chi(t).
$$

右辺の結合項を展開すると、

$$
\sum_{j\ne i}K_{ij}T_j-\left(\sum_{j\ne i}K_{ij}\right)T_i.
$$

したがって、非対角を $L_{ij}=K_{ij}$ 、対角を $L_{ii}=-\sum_{j\ne i}K_{ij}$ とすれば、結合項は $L\mathbf T$ です。例えば2点なら、

$$
L=\begin{pmatrix}-k&k\\k&-k\end{pmatrix},\qquad
L\begin{pmatrix}T_0\\T_1\end{pmatrix}
=\begin{pmatrix}k(T_1-T_0)\\k(T_0-T_1)\end{pmatrix}.
$$

高温側から低温側へ流れる熱が、2点の和で打ち消し合います。全点をまとめると、

$$
C\dot{\mathbf T}=L\mathbf T-h(\mathbf T-T_\infty\mathbf1)+\mathbf bQ_0\chi(t).
$$

106ではヒータ点P2だけに入熱するので $\mathbf b=(0,0,1,0,0)^{\mathsf T}$ 。 $\chi(t)=1$ は $0\le t<300$ 秒、 $\chi(t)=0$ は300秒以降です。
**微分する未知量は、時間とともにON/OFFするQ(t)ではなく、その振幅Q₀**とします。

### 3-3. Q感度：Q₀で1項ずつ微分する

$C,L,h,T_\infty,\mathbf b$ 、初期温度をQ₀に依存しないものとします。定義は

$$
\mathbf S_Q(t)=\frac{\partial\mathbf T(t)}{\partial Q_0}\quad[\mathrm{K/W}].
$$

熱方程式の各項を微分すると、

$$
\begin{aligned}
\frac{\partial}{\partial Q_0}(C\dot{\mathbf T})&=C\dot{\mathbf S}_Q,\\
\frac{\partial}{\partial Q_0}(L\mathbf T)&=L\mathbf S_Q,\\
\frac{\partial}{\partial Q_0}[-h(\mathbf T-T_\infty\mathbf1)]&=-h\mathbf S_Q,\\
\frac{\partial}{\partial Q_0}[\mathbf bQ_0\chi(t)]&=\mathbf b\chi(t).
\end{aligned}
$$

これらをまとめて、

$$
\boxed{C\dot{\mathbf S}_Q=(L-hI)\mathbf S_Q+\mathbf b\chi(t)},\qquad \mathbf S_Q(0)=\mathbf0.
$$

**Qを少し増やした影響が、ヒータ点から入り、結合Lで周囲へ伝わり、hで失われる**式です。
冷却に入ると右辺の入力だけが消えます。300秒時点の感度は引き継ぎ、ゼロに初期化しません。

現行コードは $Q_0=15q_{\rm scale}$ Wなので、連鎖律より

$$
\frac{\partial\mathbf T}{\partial q_{\rm scale}}
=15\frac{\partial\mathbf T}{\partial Q_0},\qquad
\mathbf S_Q=\frac1{15}\mathbf S_{q_{\rm scale}}.
$$

q_scaleで微分した場合の駆動項は $15\mathbf b\chi(t)$ です。感度のグラフにK/Wと書くなら、この換算が必要です。

### 3-4. h感度：積の微分で項を落とさない

次に $\mathbf S_h=\partial\mathbf T/\partial h$ とします。単位はK/(W/K)＝K²/Wです。
放熱項にはhが明示的に入り、温度Tもhに依存します。

$$
\frac{\partial}{\partial h}[-h(\mathbf T-T_\infty\mathbf1)]
=-(\mathbf T-T_\infty\mathbf1)-h\mathbf S_h.
$$

左の積のうち、hを微分した分が第1項、温度を微分した分が第2項です。Q₀を固定すればヒータ項のh微分は0なので、

$$
\boxed{C\dot{\mathbf S}_h=(L-hI)\mathbf S_h-(\mathbf T-T_\infty\mathbf1)},\qquad \mathbf S_h(0)=\mathbf0.
$$

初期に $\mathbf T=T_\infty\mathbf1$ なら、h感度の駆動項は0です。温度が上がって初めて放熱係数の影響が現れます。
通常の加熱条件では、hを増やすと温度が下がるので感度は負になります。比較では符号と絶対値を区別します。

**1点のモデルでも確認できます。** 初期温度を外気温とし、加熱中の温度差を $\vartheta=T-T_\infty$ とすると、

$$
C\dot\vartheta=Q_0-h\vartheta,\qquad
\vartheta(t)=\frac{Q_0}{h}(1-e^{-ht/C}).
$$

これを直接微分して、

$$
S_Q=\frac{1-e^{-ht/C}}h,\qquad
S_h=-\frac{Q_0}{h^2}(1-e^{-ht/C})+\frac{Q_0t}{hC}e^{-ht/C}.
$$

加熱開始直後の展開は $S_Q\approx t/C$ 、 $S_h\approx-Q_0t^2/(2C^2)$ 。
Q感度は時間の1次、h感度は2次で立ち上がるため、**加熱初期の短い観測だけではhが特に見えにくい**と分かります。
h感度が最大となる時刻は条件によります。「必ず冷却期が最適」とはこの式だけでは言えません。

### 3-5. 15変数を同じRK4ステージで進める

$$
\boldsymbol\zeta=\begin{bmatrix}\mathbf T\\\mathbf S_Q\\\mathbf S_h\end{bmatrix}\in\mathbb R^{15},\qquad
\dot{\boldsymbol\zeta}=f(t,\boldsymbol\zeta).
$$

この15変数は**感度解析用**の拡大変数であり、EnKFの7状態とは別です。 $M=C^{-1}(L-hI)$ と置くと右辺は

$$
f=\begin{bmatrix}
M\mathbf T+C^{-1}hT_\infty\mathbf1+C^{-1}\mathbf bQ_0\chi(t)\\
M\mathbf S_Q+C^{-1}\mathbf b\chi(t)\\
M\mathbf S_h-C^{-1}(\mathbf T-T_\infty\mathbf1)
\end{bmatrix}.
$$

RK4の各段を $k_1,\ldots,k_4$ 、時間刻みを $\Delta t$ とすると、

$$
\begin{aligned}
k_1&=f(t_n,\zeta_n),\\
k_2&=f(t_n+\Delta t/2,\zeta_n+\Delta t k_1/2),\\
k_3&=f(t_n+\Delta t/2,\zeta_n+\Delta t k_2/2),\\
k_4&=f(t_n+\Delta t,\zeta_n+\Delta t k_3),\\
\zeta_{n+1}&=\zeta_n+\frac{\Delta t}{6}(k_1+2k_2+2k_3+k_4).
\end{aligned}
$$

h感度の右辺に入る温度も、各段の途中温度を使います。最初の温度だけを4段で使い回すと、この連立系に対するRK4ではなくなります。

実装は [`integrate_sensitivity()`](../run/sensitivity_analysis.py)、予報側は [`rom_general.py`](../dacore/rom_general.py) です。前者の`dz()`は次の構造です。

```python
# Cは対角成分だけを保存した長さ5の配列。配列で割る操作がC^{-1}に対応。
dT  = M @ T  + h * T_air / C + b * Q0 * chi(t) / C
dSQ = M @ SQ + b * chi(t) / C
dSh = M @ Sh - (T - T_air) / C
```

**切替時刻の注意**：現行実装は既存予報と同じ`chi(t)`を各段で評価し、例えば298→300秒の最終段はOFFです。以下の有限差分一致は、この同じ離散化に対する検証です。連続方程式に対する精度は、刻み幅を減らす検証と区別します。切替を厳密に扱う拡張では区間を分け、加熱区間の終端段は左側入力、冷却区間の始端は右側入力を使います。

### 3-6. 有限差分は「検算」に使う

直接感度とは別に、Q₀またはhを少しだけ上下させて温度を計算します。

$$
S_Q^{\rm FD}=\frac{T(Q_0+\delta Q)-T(Q_0-\delta Q)}{2\delta Q},\qquad
S_h^{\rm FD}=\frac{T(h+\delta h)-T(h-\delta h)}{2\delta h}.
$$

同じ初期値・時刻・刻みで比べます。保存済み [`sensitivity_qh.npz`](../results/sensitivity_qh.npz) の検算値は以下です。

| 条件・比較 | 値 |
|---|---|
| 動作点 | Q₀=15 W、h=0.01535104094 W/K |
| 時間 | 0〜600秒、2秒刻み |
| 差分幅（既存コード） | δQ=0.01 W、δh=0.0001 W/K |
| Q感度：全時刻・5点の最大絶対差 | 約2.34×10⁻¹¹ K/W |
| h感度：同最大絶対差 | 約8.14×10⁻⁸ K/(W/K) |

これは既存出力の確認であり、今回新しく再計算した結果ではありません。数値的一致はROMの微分実装の検証で、実機やOpenFOAM全場への妥当性の証明ではありません。

### 3-7. 5点の感度から20,696セルの感度へ

第4回の復元式を出発点にします。各セルの平均温度は

$$
\bar T_i=\frac1{121}\sum_{k=0}^{120}T_i(5k),\qquad
\bar{\mathbf T}\in\mathbb R^{20696}.
$$

全空間を平均した1個の数ではありません。0〜300秒の61枚平均と305〜600秒の60枚平均に分けて書くなら、

$$
\bar{\mathbf T}_{0:600}=
\frac{61\bar{\mathbf T}_{\rm heat}+60\bar{\mathbf T}_{\rm cool}}{121}.
$$

今回はこの平均場と基底 $\Phi\in\mathbb R^{20696\times5}$ を全時刻で固定します。300秒で平均場だけを変えません。
$P\in\mathbb R^{5\times20696}$ は代表セルを取り出す行列、 $^+$ は疑似逆行列です。

$$
\mathbf a=(P\Phi)^+(\mathbf T_r-P\bar{\mathbf T}),\qquad
\widehat{\mathbf T}_{\rm full}=\bar{\mathbf T}+\Phi\mathbf a.
$$

平均場・基底・選点をパラメータpに依存しないものとして微分すると、

$$
\frac{\partial\mathbf a}{\partial p}=(P\Phi)^+\mathbf S_p^r,
\qquad
\boxed{\mathbf S_p^{\rm full}=A_T\mathbf S_p^r},\quad A_T:=\Phi(P\Phi)^+.
$$

$A_T$ は20,696×5です。観測ノイズ共分散と混同しないため、ここでは復元行列をRではなくAと表記します。pごとにPODを学習し直す場合は、平均場・基底の微分項も必要になり、この式だけでは済みません。

```mermaid
flowchart LR
 A["5点の感度 S：5×1"] --> B["係数感度：(PΦ)⁺S：5×1"]
 B --> C["全場感度：Φ(PΦ)⁺S：20,696×1"]
 C --> D["セル座標と対応づけて分布を表示"]
```

各セルの座標 $(x_i,y_i,z_i)$ に $\partial T_i/\partial p$ を配置すれば感度分布になります。**空間微分 $\partial T/\partial x$ を求める操作ではありません。** CSVなら`cell_id,x_m,y_m,z_m,dT_dQ_K_per_W,dT_dh_K_per_WperK`という列で保存できます。

この全場は、固定POD部分空間で復元したROM感度です。OpenFOAMの全方程式を直接微分した感度とは区別します。

### 3-8. 温度感度を変位感度へ写す

線形熱弾性で、基準温度からの差を $\Delta T$ とすると、

$$
K_u\mathbf u=H_T\Delta\mathbf T,\qquad W=K_u^{-1}H_T.
$$

$K_u$ は構造剛性、 $H_T$ は温度差から熱荷重を作る行列です。剛性・拘束・膨張係数がpで変わらないなら、

$$
\frac{\partial\mathbf u}{\partial p}=W\frac{\partial\mathbf T_{\rm FEM}}{\partial p}.
$$

セル温度をFEM節点へ写す固定行列をJとすれば、 $\mathbf S_p^{u}=WJ\mathbf S_p^{\rm full}$ 。CFDセルとFEM節点の数・位置が異なるため、写像Jを省略して両者を同じ行列として掛けてはいけません。

106のPODモード応答では $\mathbf u=\mathbf u_0+D\mathbf a$ なので、

$$
\boxed{\mathbf S_p^u=A_u\mathbf S_p^r},\qquad A_u:=D(P\Phi)^+.
$$

$D$ にはFrontISTRで事前計算したモード変位応答を使います。現在の`disp_operator.npz`の`Dmode`は**上面A/OのUz 2成分×5モード**で、全FEM節点の3方向変位ではありません。全節点版を出力するにはDを全節点・全成分へ拡張する必要があります。

変位差 $d=u_{z,A}-u_{z,O}$ なら、

$$
\frac{\partial d}{\partial p}
=\frac{\partial u_{z,A}}{\partial p}-\frac{\partial u_{z,O}}{\partial p}.
$$

AとOが同じだけ動けば差の感度は0です。したがって「大きく動く点」と「目的の変位差を識別できる点」は必ずしも同じではありません。

### 3-9. 大きな感度でもQ・hが区別できない場合

候補温度点iで、観測時刻 $t_1,\ldots,t_M$ の感度を縦に並べます。

$$
G_i=\begin{pmatrix}
\partial T_i(t_1)/\partial Q_0&\partial T_i(t_1)/\partial h\\
\vdots&\vdots\\
\partial T_i(t_M)/\partial Q_0&\partial T_i(t_M)/\partial h
\end{pmatrix}\in\mathbb R^{M\times2}.
$$

小さなパラメータ変化に対して $\delta\mathbf y_i\approx G_i\delta\boldsymbol\theta$ 、 $\boldsymbol\theta=(Q_0,h)^{\mathsf T}$ です。
2列が比例していれば、Qの変化とhの変化が同じ観測パターンを作り、どちらが原因か区別できません。

温度と変位の単位を揃えるため、観測誤差共分散を $R_{\rm obs}$ として白色化します。さらにQとhの尺度も異なるので、比較全体で固定した尺度 $s_Q,s_h$ を用います。

$$
S_\theta=\mathrm{diag}(s_Q,s_h),\quad
\delta\theta=S_\theta\delta\eta,\quad
\widetilde G=R_{\rm obs}^{-1/2}G S_\theta,\quad F=\widetilde G^{\mathsf T}\widetilde G.
$$

$\eta$ は無次元パラメータ。例えば独立ノイズなら温度行をσT、変位行をσuで割ります。変位をµmで扱うならσuもµmにします。sQ・shは事前の不確かさ等で設定し、候補点ごとに変えません。
**ノイズだけ規格化してもQとhの単位差は残る**ため、最小固有値による順位付けにはこのパラメータ尺度の明示が必要です。

次は実データではなく、規格化後の2時刻×2パラメータの手計算例です。

$$
\widetilde G_A=\begin{pmatrix}1&1\\2&2\end{pmatrix},\quad
F_A=\begin{pmatrix}5&5\\5&5\end{pmatrix},\quad
\lambda(F_A)=(10,0).
$$

Aでは信号は大きくても1方向しか区別できません。一方、

$$
\widetilde G_B=\begin{pmatrix}1&0\\0&1\end{pmatrix},\quad F_B=I,\quad\lambda(F_B)=(1,1).
$$

Bは両パラメータを別々に区別できます。最小固有値を最大化する基準なら、AよりBを選びます。

| 指標 | 計算・読み方 |
|---|---|
| 列ノルム | $\lVert\widetilde g_Q\rVert,\lVert\widetilde g_h\rVert$ ：各パラメータがどれほど観測へ現れるか |
| 余弦類似度 | $\widetilde g_Q^{\mathsf T}\widetilde g_h/(\lVert\widetilde g_Q\rVert\lVert\widetilde g_h\rVert)$ ：絶対値1なら比例。ゼロ列では未定義 |
| 最小特異値 | $\sigma_{\min}(\widetilde G)$ ：最も見えにくいパラメータ方向の情報 |
| 最小固有値 | $\lambda_{\min}(F)=\sigma_{\min}(\widetilde G)^2$ |
| 条件数 | $\sigma_{\max}/\sigma_{\min}$ ：0割なら無限大。小さくても信号全体が微小なら十分ではない |

独立な温度観測iと変位観測jを組み合わせると、行を縦に積むため $F_{ij}=F_i^T+F_j^u$ 。この上付きTは温度の種類を表し、転置は $\mathsf T$ で書き分けます。
ノイズに相関がある場合は単純に加算せず、結合した $R_{\rm obs}$ で計算します。

$$
(i^*,j^*)=\mathop{\rm arg\,max}_{i,j}\lambda_{\min}(F_{ij}).
$$

これは指定した候補・時刻・ノイズ・尺度の中での選定基準です。Fisher情報は局所的な情報量であり、その順位が有限メンバーEnKFの誤差順位と一致するかは次に検証します。
§1の分散減少Δは事前共分散を含む評価、Fは観測からの情報の評価です。線形ガウス近似なら無次元パラメータの事後共分散は $(B_\eta^{-1}+F)^{-1}$ となり、両者がつながります。

### 3-10. フルFEMで直接求める場合と実装の範囲

ROMを使わない熱FEMの半離散式 $C\dot T+K(p)T=f(p)$ で、Cがpに依存しないなら、積の微分から

$$
C\dot S_p+K S_p=f_p-K_pT.
$$

$f_p=\partial f/\partial p$ 、 $K_p=\partial K/\partial p$ です。Cも依存する場合は右辺に $-C_p\dot T$ が加わります。
対流境界で $K=K_c+h_sK_h$ 、 $f=f_Q+h_sf_h$ と置くと、

$$
C\dot S_{h_s}+(K_c+h_sK_h)S_{h_s}=f_h-K_hT.
$$

ここで $h_s$ は表面熱伝達率です。ROMのhと同じ数値をそのまま入れてはいけません。
定常なら時間微分を除いて $KS_p=f_p-K_pT$ を解きます。これらは熱FEMの一般式で、OpenFOAMのCHT全体では流体の速度・圧力・界面も含む連成系の微分が必要です。

| 段階 | 現行ファイルと確認できた範囲 |
|---|---|
| 熱ROMとQ換算 | `dacore/rom_general.py`：Q=15 q_scale、hはW/K |
| 15状態RK4・中心差分 | `run/sensitivity_analysis.py`：実装済み。既存NPZに検算差を保存 |
| POD基底・平均・選点 | `results/qdeim_points.npz`：mean、pod_modes、cell_idx、cell_centres、times |
| 変位写像 | `results/disp_operator.npz`：現行Dmodeは2×5。全節点の3成分出力ではない |
| 感度保存 | `results/sensitivity_qh.npz`：t、T、S_Q、S_h、検算値。保存される感度は5代表点 |
| 全場CSV/VTK・全節点感度・Fisherランキング | 上記の式に沿う拡張範囲。`sensitivity_analysis.py`単体の既存出力には含まれない |

**次の比較設計**：ROMノードと真値・EnKF・観測ノイズ・seed群・更新時刻・センサ数を固定し、Q-DEIM、単純Q/h感度、Q/h識別性の3配置を比べます。温度場RMSE、Q/h誤差、観測に使わない変位差RMSEを、平均とseed間のばらつきで示します。
学習と同じROMでの双子実験に加え、未学習条件・独立した高忠実度真値を使う検証が必要です。以下の既存図は出発点となる実験であり、全場Fisher最適化まで実証した図ではありません。

### 3-11. 既存の感度時刻歴と、観測時間窓の結果



$\Delta$ の分子 $\mathrm{Cov}(X,y)$ は**時刻の関数**でもあります。直接感度（熱方程式を
$Q,h$ で微分したもの）を計算すると:

![∂T/∂Qと∂T/∂hの時刻歴](img/sensitivity_qh_timeseries.png)

- $\partial T/\partial Q$ は**加熱期に大きく**、冷却で減衰 → $Q$ と温度観測の連動は**加熱期**にある
- $\partial T/\partial h$ は温度が上がった**冷却期に最大** → $h$ の情報は**後半**にある

実験でも、同じ温度2点で**同化する時間窓だけ**変えると:

![同化時間窓を変えたQ/h推定](img/qh_time_window.png)

*$Q$ ：加熱期のみ0.57 W ≈ 全期間0.58 W、**冷却期のみだと4.23 W（7倍悪化）**＝ $Q$ は加熱期の観測が担う。
$h$ ：冷却期を足すと7.6→5.3 mW/Kと相対改善（ただし温度応答＜ノイズで依然難）。*

つまり最適観測は **どこで（空間）× 何を（種類）× いつ（時間）** の3軸で、どれも対象依存です。

---

## 4. 実データで検証（ROM＋EnKF/OI、すべて本リポジトリの計算）

### 4-1. まず観測点の位置

![観測点の位置](img/obs_points_progression.png)

### 4-2. 対象＝全温度場：今回の候補比較では「散らす」配置が改善

![配置方法×推定対象の誤差表](img/optimal_placement_table.png)

*温度2点のEnKF比較。**感度最大＝近接2点は0.91 Kと2.6倍悪化**、散らす配置は0.35 K。
今回の比較は観測の非冗長性の重要さを示唆しますが、Q-DEIMが全候補中で最適であることを証明するものではありません。*

### 4-3. 対象＝変位：変位センサを「高W点」へ（種類と場所が変わる）

![高W/低W変位点の比較（OI）](img/blog_oi_disp_selection.png)

*同じ変位2点でも**低感度W点はほぼ効かず（2.43 K）、高W点は0.86 K**。
§2の言葉では「対象＝変位・温度場」に対し $g=w^{\mathsf T}\varphi$ が大きい観測を選んだ、ということ。*

さらに個別の変位で見ると、**温度センサだけでは「センサに近い側」しか合いません**:

![個別のUz(A)/Uz(O)と差](img/blog_disp_timeseries_truth_vs_da_points.png)

### 4-4. 観測を増やす順番も式の通り

![温度1点→2点→＋変位の段階比較](img/opencae_sensor_progression.png)

*温度RMSE 0.617→0.197→0.159 K、未観測変位差 0.869→0.151→**0.077 µm**。
「まず独立な温度2点で場を張り、変位で $Q$ ・変位系の対象を締める」＝ $\Delta$ 表の順序そのもの。*

### 4-5. 仕上げ：未観測の6点でも（動く絵つき）

![未観測の温度3点・変位3点](img/unobs_multi_timeseries.png)

同化が全場を直せる理由（少数点→PODで全場復元）はアニメで一目です:

![5点から全温度場を復元](img/blog_temp_recon_anim.gif)

---

## 5. まとめ ― 対象別・最適観測の設計表

![推定対象×最適なセンサ種類・領域のまとめ](img/placement_summary_table.png)

| 推定対象 | 既存の限定条件で有効だった観測 | 理屈（ $\Delta$ の言葉） |
|---|---|---|
| 全温度場 | 温度を**空間に散らす** | 独立な方向を張ると合計の $\mathrm{Cov}$ が最大 |
| 変位（未観測含む） | **変位を高W点**へ | $g=w^{\mathsf T}\varphi$ が大きい観測 |
| 発熱量 $Q$ | 温度で可・今回の候補間では場所の効果が小さい、**時間＝加熱期** | $\mathrm{Cov}(Q,y(t))$ が加熱期に集中 |
| 放熱 $h$ | **冷却期**の観測で相対改善（依然難） | $\partial T/\partial h$ が冷却期に最大だが信号＜ノイズ |

1. **観測の価値は1本の式** $\Delta=\mathrm{Cov}(X,y)^2/(\mathrm{Var}(y)+r)$ で測れる。
2. 分子に**対象**が入るので、**対象によって最良候補が異なる場合がある**（3変数の手計算で★が行ごとに動いた）。
3. 軸は**場所だけでなく種類と時間**。実データ（ROM＋EnKF/OI）でも表の通りになった。

> 学会発表との対応：この論理は
> [発表ポスター](https://kamakiri1225.github.io/thermal-da-demo/posters.html) の
> **計算工学講演会（直接感度と観測設計）**・**計算力学講演会（推定対象依存の観測配置）**の骨格です。

---

### 参考（元になった実装・詳細）

- $\Delta$ の元になる更新式・変位→温度の導出：`blog_002_oi_data_assimilation.md` §2, §3-7
- 直接感度の実装と検証：`run/sensitivity_analysis.py`（有限差分と一致）
- 時間窓実験：`run/time_window_qh.py` ／ 空間配置比較：`run/optimal_placement_table.py`
- 段階比較・未観測評価：`run/opencae_sensor_progression.py`, `run/make_unobs_multi_fig.py`
