# POD と Q-DEIM ― 数式展開とプログラム手続き（学習用）

「代表点を勘で置く」のをやめ、**POD（固有直交分解）でモードを抽出 → Q-DEIM で代表点を選定 →
その点でROMを組んで校正**するまでを、数式とコードの両方でたどる。
本フォルダ 106 の実装（`run/select_points_qdeim.py`, `run/build_calibrate_rom.py`,
`dacore/rom_general.py`）に対応。

---

## 前提：104 の実ソルバ・データ同化について

106 のROMは、**104 の実ソルバ・データ同化を高速化するために作った**。
「どのプログラムで・どんなアルゴリズム（確率的アンサンブルEnKF）で同化していたか」は、
独立ドキュメント **[`00_data_assimilation_algorithm.md`](00_data_assimilation_algorithm.md)** に
丁寧にまとめた（プログラム構成 A → 状態・観測 B → 双子実験 C → 解析の数式 D → 重さと106との関係 E）。

要点だけ:

- 実行入口 `run/run_openfoam_fem_enkf.py` → ドライバ `daof/of_fem_twin.py::run_fem_twin()`
  → 解析本体 `dacore/enkf.py::enkf_update()`。
- アルゴリズムは**確率的EnKF（摂動観測法）**。ゲイン $K=C_{zy}(C_{yy}+R)^{-1}$ で観測へ引き寄せる。
- **106 は予報モデルだけを ROM に置き換えた同じEnKF**。以下 §0 からその ROM の作り方を追う。

---

## 0. 全体の流れ

```
OpenFOAMのCHT温度場(121枚 × 20696セル)
        │  ① 平均を引いてSVD
        ▼
   PODモード U（場の「型」）と エネルギー σ²
        │  ② Q-DEIM（枢軸QR）で r 点を選ぶ
        ▼
   代表5点（データから最適に選定）
        │  ③ その5点をノードに一般化ROMを組み、校正
        ▼
   軽いROM（0→600sが数秒、残差0.014K）
        │  ④ 拡大状態EnKFでデータ同化（別ドキュメント）
        ▼
   温度場・Q・h を少数観測から推定
```

---

## 1. POD の数式展開 ― なぜ SVD で「型」が出るのか

### 1.1 スナップショット行列

各時刻の全セル温度を縦ベクトルにし、時間方向に並べる（$N$=セル数, $m$=時刻数）:

$$X=\big[\,u(t_1)-\bar u,\ \dots,\ u(t_m)-\bar u\,\big]\in\mathbb{R}^{N\times m},\qquad
\bar u=\frac1m\sum_k u(t_k)$$

- $u(t_k)$ … **その時刻の温度場まるごと**（1列＝全セル分の縦ベクトル、長さ $N$）
- $\bar u$ … **平均場**。ただし全体で1つの数ではなく、**セルごとに、そのセル自身の
  時間平均**を並べた縦ベクトル（＝1枚の空間分布, 長さ $N$）
- $u(t_k)-\bar u$ … **同じセルどうしの引き算**（各セルから、そのセルの時間平均を引く）

#### 行と列を展開して確認する

**はい、 $\bar u$ は各セルごとの時間平均である。** 本節の小文字 $u$ は温度場を表す記号で、
熱変形の変位ではない。 $u_i(t_k)=T_i(t_k)$ と読み替えてよい。
まず平均を引く前の温度行列を $S_T$ と書き、平均を引いた行列 $X$ と区別する：

$$u(t_k)=\begin{bmatrix}u_1(t_k)\\u_2(t_k)\\\vdots\\u_N(t_k)\end{bmatrix},\qquad
S_T=\begin{bmatrix}
u_1(t_1)&u_1(t_2)&\cdots&u_1(t_m)\\
u_2(t_1)&u_2(t_2)&\cdots&u_2(t_m)\\
\vdots&\vdots&\ddots&\vdots\\
u_N(t_1)&u_N(t_2)&\cdots&u_N(t_m)
\end{bmatrix}\quad(N\times m).$$

**1行は同じセルの温度時刻歴、1列は同じ時刻の全セル温度場。**
今回の $N=20,696$、 $m=121$ で、 $t_1=0$ s、 $t_2=5$ s、…、 $t_{121}=600$ s。
各行を横方向に平均するので、

$$\bar u=\begin{bmatrix}
\bar u_1\\\bar u_2\\\vdots\\\bar u_N
\end{bmatrix}
=\begin{bmatrix}
\{u_1(t_1)+u_1(t_2)+\cdots+u_1(t_m)\}/m\\
\{u_2(t_1)+u_2(t_2)+\cdots+u_2(t_m)\}/m\\
\vdots\\
\{u_N(t_1)+u_N(t_2)+\cdots+u_N(t_m)\}/m
\end{bmatrix}\quad(N\times1).$$

セル1にはセル1の平均 $\bar u_1$、セル2にはセル2の平均 $\bar u_2$ を使う。
この平均ベクトルを121列分並べてから引くことに相当する：

$$X=S_T-\bar u\,\mathbf1_m^\top
=S_T-\begin{bmatrix}
\bar u_1&\bar u_1&\cdots&\bar u_1\\
\bar u_2&\bar u_2&\cdots&\bar u_2\\
\vdots&\vdots&\ddots&\vdots\\
\bar u_N&\bar u_N&\cdots&\bar u_N
\end{bmatrix},\qquad \mathbf1_m^\top=[1,1,\ldots,1].$$

従って、実際にPODへ渡す行列は

$$X=\begin{bmatrix}
u_1(t_1)-\bar u_1&u_1(t_2)-\bar u_1&\cdots&u_1(t_m)-\bar u_1\\
u_2(t_1)-\bar u_2&u_2(t_2)-\bar u_2&\cdots&u_2(t_m)-\bar u_2\\
\vdots&\vdots&\ddots&\vdots\\
u_N(t_1)-\bar u_N&u_N(t_2)-\bar u_N&\cdots&u_N(t_m)-\bar u_N
\end{bmatrix}.$$

各行の平均は0になる： $\frac1m\sum_k X_{ik}=0$。各列の空間平均が0になるとは限らない。
これは学習温度履歴の時間平均であり、EnKFのメンバー平均とは別の操作である。

| 本節の記号 | 実装 `run/select_points_qdeim.py` | 配列の形 |
|---|---|---|
| 生の温度行列 $S_T$ | `X` | `(20696, 121)` |
| 各セルの時間平均 $\bar u$ | `mean = X.mean(axis=1)` | `(20696,)` |
| 縦ベクトルとしての $\bar u$ | `mean[:, None]` | `(20696, 1)` |
| 平均を引いた行列 $X$ | `Xc = X - mean[:, None]` | `(20696, 121)` |

`axis=1`は列の方向、つまり121時刻を平均する指定。
`mean[:, None]`は平均を縦に並べ直し、NumPyが各列に同じ平均を適用する。
**本文の $X$ とコードの`X`は意味が違い、本文の $X$ はコードの`Xc`に対応する。**

### 1.1b 具体例で見る「平均場を場所ごとに引く」

言葉だと分かりにくいので、**3セル(A,B,C)×4時刻**の小さな例で示す。

**まず生の温度**（行＝場所、列＝時刻。右端が各行の時間平均 $\bar u_i$）:

| セル | $t_1$ | $t_2$ | $t_3$ | $t_4$ | 時間平均 $\bar u_i$ |
|---|---|---|---|---|---|
| A（ヒータ側） | 20 | 22 | 24 | 26 | **23** |
| B（中間） | 20 | 21 | 22 | 23 | **21.5** |
| C（反対側） | 20 | 20.5 | 21 | 21.5 | **20.75** |

- $u(t_1)$ ＝ 1列まるごと ＝ `[20, 20, 20]`（時刻 $t_1$ の温度場）
- $\bar u$ ＝ `[23, 21.5, 20.75]`（**場所ごとに違う**時間平均。A は自分の23、B は21.5…）

**各セルから「そのセル自身の平均」を引く** → これが $X=[u(t_k)-\bar u]$:

上の表を実際の行列演算として書くと（温度は℃、差はK）、

$$S_T=\begin{bmatrix}20&22&24&26\\20&21&22&23\\20&20.5&21&21.5\end{bmatrix},\qquad
\bar u=\frac14\begin{bmatrix}20+22+24+26\\20+21+22+23\\20+20.5+21+21.5\end{bmatrix}
=\begin{bmatrix}23\\21.5\\20.75\end{bmatrix},$$

$$X=\begin{bmatrix}20&22&24&26\\20&21&22&23\\20&20.5&21&21.5\end{bmatrix}
-\begin{bmatrix}23&23&23&23\\21.5&21.5&21.5&21.5\\20.75&20.75&20.75&20.75\end{bmatrix}
=\begin{bmatrix}-3&-1&1&3\\-1.5&-0.5&0.5&1.5\\-0.75&-0.25&0.25&0.75\end{bmatrix}.$$

例えばセルB・時刻 $t_3$ の成分は $X_{B,3}=22-21.5=0.5$ K。
復元するときは同じ平均を戻し、 $u_B(t_3)=21.5+0.5=22$ ℃となる。

| セル | $t_1$ | $t_2$ | $t_3$ | $t_4$ |
|---|---|---|---|---|
| A（−23） | −3 | −1 | +1 | +3 |
| B（−21.5） | −1.5 | −0.5 | +0.5 | +1.5 |
| C（−20.75） | −0.75 | −0.25 | +0.25 | +0.75 |

各列は「その時刻に、各場所が**自分の平均からどれだけ上下したか（変動）**」。
PODはこの変動表から**共通の型**を取り出す（この例なら「Aが強く・Bが中・Cが弱く、
時間とともに一様に増える」という型1つでほぼ説明できる）。
平均場＝定常的な下地、モード＝そこからの変動の型、という役割分担。

コードでは `X - mean[:, None]`（`mean = X.mean(axis=1)` ＝**行ごと＝場所ごとの時間平均**）
がこの「各行から自分の平均を引く」操作にあたる。

> **記号の注意**：`load_snapshots()` の返り値は `(ts, cell_xyz, X)`。
> `ts`=時刻配列 [s]（長さ $m$）、`cell_xyz`=各セル中心の座標 $(N\times3)$、
> `X`=温度のスナップショット行列 $(N\times m)$。
> （ROM側の熱容量 $C$ とは別物。座標を $C$ と書くと紛らわしいので `cell_xyz` と呼ぶ。）

### 1.2 最適な型を求める変分問題

「場を1つの型 $\varphi$（$\|\varphi\|=1$）で最もよく表す」とは、
全スナップショットの $\varphi$ 方向成分の二乗和を最大化すること:

$$\max_{\|\varphi\|=1}\ \sum_{k=1}^{m}\big|\varphi^\top (u(t_k)-\bar u)\big|^2
=\max_{\|\varphi\|=1}\ \varphi^\top \big(XX^\top\big)\varphi$$

ラグランジュ未定乗数 $\lambda$ で $\varphi^\top XX^\top\varphi-\lambda(\varphi^\top\varphi-1)$ を
$\varphi$ で微分してゼロと置くと、**固有値問題**になる:

$$\boxed{\ (XX^\top)\,\varphi=\lambda\,\varphi\ }$$

- 固有ベクトル $\varphi_k$ ＝ **PODモード（型）**
- 固有値 $\lambda_k=\sigma_k^2$ ＝ その型の**エネルギー**（大きいほど支配的）
- $C=XX^\top$ は空間相関行列（共分散に相当, $N\times N$）

### 1.3 SVD との同値（実際の計算法）

$X$ を特異値分解 $X=U\Sigma V^\top$ すると

$$XX^\top=U\Sigma V^\top V\Sigma U^\top=U\,\Sigma^2\,U^\top$$

なので、**$U$ の各列がそのまま $C=XX^\top$ の固有ベクトル＝PODモード**、
$\sigma_k^2$ が固有値。つまり巨大な $C$（$20696^2$）を作らず、SVD一発でPODが得られる。

### 1.4 なぜ「最良」か（Eckart–Young の定理）

上位 $r$ モードで打ち切った近似
$X_r=\sum_{k=1}^r\sigma_k u_k v_k^\top$ は、
**ランク $r$ の全行列の中でフロベニウス誤差 $\|X-X_r\|_F$ を最小にする**。
これが「POD＝最もエネルギー効率の良いモード分解」の数学的な裏づけ。

### 1.5 スナップショット法（なぜ速いか）

空間 $N$≫時間 $m$（本ケース 20696≫121）なので、 $N\times N$ でなく
**小さい $m\times m$ の $X^\top X$ の固有問題**を解けばよい。
$X^\top X\,v_k=\sigma_k^2 v_k$ を解き、モードは
$$\varphi_k=\frac{1}{\sigma_k}X v_k\quad(\text{＝スナップショットの重み付き和})$$
で得る。economy SVD（`full_matrices=False`）が実質これを行う。

---

## 2. POD の実装アルゴリズム（プログラム手続き）

`run/select_points_qdeim.py` の中核。numpy標準の数行で済む:

```python
# ① スナップショット行列 X (N=20696セル × m=121時刻) を作る
ts, C, X = load_snapshots()            # OpenFOAMの各時刻 solid/T を読む

# ② 平均を引く（変動成分だけ見る）
mean = X.mean(axis=1)
Xc   = X - mean[:, None]

# ③ SVD（＝PODモード抽出）
U, S, Vt = np.linalg.svd(Xc, full_matrices=False)   # U:(N,m) S:(m,) Vt:(m,m)

# ④ 各モードのエネルギー比と累積
energy = S**2 / (S**2).sum()
cum    = np.cumsum(energy)             # 何モードで何%説明できるか
```

- `U[:,k]` … 第k PODモード（空間の型）
- `S[k]**2/ΣS**2` … そのモードのエネルギー比
- 本ケースの結果：**1モードで91.8%、2モードで99.9%**（＝場は実質2〜3自由度）

---

## 3. Q-DEIM ― 代表点をどう選ぶか（アルゴリズム）

### 3.1 狙い：少数点の値から場を復元できる点を選ぶ

$r$ 個の点 $P=\{p_1,\dots,p_r\}$ の温度だけから、PODモードを使って全場を復元したい
（**gappy 再構成 / DEIM 補間**）:

$$\hat u=\bar u+U_r\,a,\qquad a=\big(U_r[P,:]\big)^{+}\big(u[P]-\bar u[P]\big)$$

ここで $U_r[P,:]$ は「モード行列の、選んだ点の行だけ」を抜いた $r\times r$ 行列。
この復元が安定なのは $U_r[P,:]$ が**良条件（可逆で条件数が小さい）**なとき。
だから「$U_r[P,:]$ の条件数を小さくする点集合 $P$」を選ぶ＝Q-DEIM。

### 3.2 アルゴリズム：列枢軸QR（Drmač–Gugercin, 2016）

Q-DEIM は「モード行列の転置 $U_r^\top$ に**列枢軸付きQR分解**をかけ、
最初の $r$ 個の枢軸列（＝行）を補間点にする」だけ:

```python
from scipy.linalg import qr
# 先頭 r モード U[:, :r] の転置に列枢軸QR
_, _, piv = qr(U[:, :r].T, pivoting=True)   # piv: 重要な行(=点)の順
points = list(piv[:r])                       # 最初の r 個が代表点(セル番号)
coords = C[points]                           # その点の座標
```

枢軸QRは「互いに最も独立な（張る空間が大きい）行を順に選ぶ」ので、
**モードを最もよく分離できる点＝場を最もよく覆う点**が出る。
これが「勘」ではなく「データから系統的に選ぶ」の中身。

### 3.3 本ケースで選ばれた5点

| 点 | 座標(mm) | 特徴 |
|---|---|---|
| P0 | (21, −2, 52) | 内周・ヒータ側・中央 |
| P1 | (−14, 33, 52) | +Y側・中央 |
| P2 | (37, 7, 51) | 外周・ヒータ側・中央（←ヒータ最近傍） |
| P3 | (35, 1, 1) | 外周・ヒータ側・**底面** |
| P4 | (−35, 1, 1) | 外周・反対側・**底面** |

勘の流用点（すべて中央高さ）と違い、**底面(P3,P4)も選んで場を広くカバー**している。

---

## 4. 選んだ点で ROM を組み、校正する

### 4.1 一般化した集中定数モデル（数式）

任意配置の $n$ 点をノードとし、全点対をコンダクタンスで結ぶ（`dacore/rom_general.py`）:

$$C_i\frac{dT_i}{dt}=\sum_{j\ne i}K_{ij}(T_j-T_i)+q_i(t)-h\,(T_i-T_\text{air})$$

- $K_{ij}$：ノード間コンダクタンス [W/K]（対称・全点対、 $n$=5 なら 10本）
- $q_i$：ヒータ発熱（**ヒータ最近傍ノード**＝P2 に投入）
- $h$：放熱係数 [W/K]

行列形にすると線形システム $\dot{\mathbf T}=M\mathbf T+\mathbf b$ になり、RK4 で一括前進積分できる。

### 4.2 校正（未知の $C,K,h$ を実データに合わせる）

選んだ5点の**OpenFOAM温度履歴**を目標に、 $C[5],K[10],h$ を最小二乗で同定
（`run/build_calibrate_rom.py`）:

```python
def resid(vec):                       # vec = [C(5), K上三角(10), h(1)] = 16未知
    C = vec[:5]; K = tri_to_matrix(vec[5:15], 5); h = vec[15]
    Ymodel = integrate_rom(C, K, h)   # ROMを前進積分して5点温度履歴を作る
    return (Ymodel - Yobs).ravel()    # OpenFOAMとの差
sol = least_squares(resid, x0, bounds=(lb, ub))
```

### 4.3 結果（検証）

- **当てはめ残差 RMSE = 0.014 K**（最大0.051K）＝OpenFOAMをほぼ完全に再現
- 全熱容量 $\sum C_i=1211$ J/K は鋼2.49kgの $\rho c_p V\approx1195$ J/K とほぼ一致
  （数字合わせでなく物理的にも妥当）

![校正結果](img/rom_calib_fit.png)

*実線=OpenFOAM、破線=POD選定5点ROM。ほぼ重なる。*

---

## 5. まとめ（このフォルダの主張）

1. **代表点は勘でなく POD＋Q-DEIM でデータから選ぶ**（§2, §3）。
   Q-DEIM＝モード行列に列枢軸QRをかけ、場を最もよく覆う点を選ぶアルゴリズム。
2. 選んだ点で一般化ROMを組み、校正すると **0.014K でOpenFOAMを再現**（§4）。
3. この上で拡大状態EnKFにより温度場・ $Q$ ・ $h$ を推定する（→ `02_da_on_pod_rom.md` 予定）。

再現: `run/select_points_qdeim.py` → `run/build_calibrate_rom.py`
