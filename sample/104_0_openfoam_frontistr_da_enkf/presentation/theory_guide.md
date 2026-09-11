# OI・KF・EnKF・PF：理論式の展開と104への接続

共通の推定問題から導出し、具体例・実装へつなげる読み物版。過程ノイズ共分散は、ヒーター発熱量Qと区別してΩと表記する。

## 第0部　予報と観測を統合する理論

**OI → KF → EnKF → PF**

記号 → 仮定 → 式の展開 → 数値例 → 104の実装

理論式と具体例を両方残し、式の各項が何をしているかを追う。

**要点：前提を置き、導出し、数値を代入する。この順で理解する。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## まず、何を推定する問題なのか

観測できない状態xを、時刻kまでの観測から推定する。


$$
x_k=f_k(x_{k-1})+\eta_k,\qquad y_k=h_k(x_k)+\nu_k
$$


x：温度場など、f：時間発展、h：観測に対応する量の計算。

η：モデル・過程の誤差、ν：測定誤差。状態と観測は次元が違ってよい。

**要点：104ではxに固体温度と発熱量Qを含め、hで温度・変位を計算する。**

ここでは加法ノイズの状態空間モデルを導入する。一般の非加法モデルにもベイズフィルタは拡張できる。

参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 記号と単位：発熱量Qと誤差共分散を混同しない

| 記号 | 意味 |
|---|---|
| xᶠ, xᵃ | 観測を使う前の予報／使った後の解析（同化結果） |
| Pᶠ, Pᵃ | 状態の予報／解析誤差共分散 |
| B | OIで与える背景誤差共分散 |
| Ω, R | 過程ノイズ／観測ノイズの共分散 |
| F, H | 線形な時間発展／観測の行列 |
| K | カルマンゲイン（残差から状態補正への係数） |

**要点：過程ノイズ共分散はΩと表記する。104の発熱量Q [W]とは別物。**

文献で過程ノイズ共分散をQと表す場合が多いため、ここでは混同防止のためΩを使用する。

参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 分散と共分散が表すもの




$$
\bar{x}=\mathrm{E}[x],\qquad P=\mathrm{E}[(x-\bar{x})(x-\bar{x})^T]
$$





$$
P_{ij}=\mathrm{E}[(x_i-\bar{x}_i)(x_j-\bar{x}_j)]
$$


対角Piiは成分iの分散。非対角Pijはiとjの変動の関係。

相関係数は Pij / √(Pii Pjj)。共分散には単位があり、相関係数は無次元。

**要点：未観測の温度やQを更新できるのは、観測量との共分散を使うため。**

以下のPは推定誤差についての共分散。EnKFではメンバーのばらつきでこれを近似する。

参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 逐次推定：予報してから、新しい観測で更新する

```text
前時刻の解析 → モデルで予報 → 現時刻の観測を同化 → 現時刻の解析
      ↑                                                    │
      └──────────────── 次の時刻へ ───────────────────────┘
```

KF・EnKF・PFは、状態と不確かさをこのループで伝える。
OIも逐次解析に使えるが、通常は予報誤差共分散を所与のBで近似する。

**要点：「逐次」は観測が来るたびに更新すること。OI単独は時間発展の方程式ではない。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## ベイズ更新：どの状態が観測を説明しやすいか




$$
p(x_k\mid y_{1:k})=\frac{p(y_k\mid x_k)\,p(x_k\mid y_{1:k-1})}{p(y_k\mid y_{1:k-1})}
$$


右辺の分子は **尤度 × 事前分布**、左辺が事後分布。

尤度：その状態なら、この観測が出る確からしさ。
事前分布：観測を使う前の状態の不確かさ。
分母：全体の確率を1にする正規化。

**要点：予報に合うかと、観測に合うかを合わせて状態を評価する。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 予報分布：前の事後分布をモデルで伝える




$$
p(x_k\mid y_{1:k-1})=\int p(x_k\mid x_{k-1})p(x_{k-1}\mid y_{1:k-1})\,dx_{k-1}
$$


全ての前状態について「そこにいる確率×そこから進む確率」を足す。

KF：線形ガウスの平均・共分散を計算する。
EnKF：各メンバーを前進する。PF：各粒子を前進し、重みを引き継ぐ。

**要点：同じ予報・更新の問題を、違う表現と近似で解く。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 第0部A　最適内挿法（OI）

**Optimal Interpolation：最適内挿法／最適補間法**

背景値と観測を、誤差共分散に基づいて統合する。

まず1時刻での「最も妥当な補正」を導出する。

**要点：IDWの距離補間とOIの統計的な解析は、目的も重みの由来も異なる。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OIの前提：背景と観測は誤差を持つ




$$
x^b=x^{true}+e_b,\qquad y=Hx^{true}+e_o
$$





$$
\mathrm{E}[e_b]=\mathrm{E}[e_o]=0,\quad\mathrm{Cov}(e_b)=B,\quad\mathrm{Cov}(e_o)=R
$$


背景誤差と観測誤差は無相関とする。Hは既知の線形観測行列。
BとRは正定値として導出する。

**要点：不偏な誤差と共分散を仮定し、どちらをどれだけ信頼するかを決める。**

ガウス分布を仮定すれば後の解は事後平均かつMAP。ガウス性なしでも線形不偏推定の最小分散として同じゲインを導ける。

参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OI導出①：背景と観測からのずれを評価する




$$
J(x)=\frac{1}{2}(x-x^b)^TB^{-1}(x-x^b)+\frac{1}{2}(y-Hx)^TR^{-1}(y-Hx)
$$


第1項：背景からのずれ。第2項：観測からのずれ。

誤差分散が小さい方向のずれほど、強く罰する。
ガウス仮定なら、これは事後確率の負の対数（定数を除く）。

**要点：背景と観測の両方を、誤差の大きさで重み付けして合わせる。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OI導出②：Jを微分して0と置く




$$
\nabla_xJ=B^{-1}(x-x^b)-H^TR^{-1}(y-Hx)
$$





$$
\nabla_xJ=0\ \Longrightarrow\ (B^{-1}+H^TR^{-1}H)x^a=B^{-1}x^b+H^TR^{-1}y
$$


転置Hᵀは、観測空間の残差を状態空間へ戻す役割を持つ。

**要点：最小化の条件から、解析値xᵃの連立一次方程式が得られる。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OI導出③：背景との差δxで書き直す




$$
\delta x=x^a-x^b,\qquad d=y-Hx^b
$$





$$
(B^{-1}+H^TR^{-1}H)\delta x=H^TR^{-1}d
$$





$$
\delta x=(B^{-1}+H^TR^{-1}H)^{-1}H^TR^{-1}d
$$


dは観測−背景の残差（イノベーション）。

**要点：状態全体を解き直す式を、「背景＋残差に応じた補正」の式へ変える。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OI導出④：観測空間の逆行列へ変形する




$$
A=B^{-1}+H^TR^{-1}H,\quad S=HBH^T+R
$$





$$
A(BH^TS^{-1})=H^TR^{-1}(R+HBH^T)S^{-1}=H^TR^{-1}
$$





$$
A^{-1}H^TR^{-1}=BH^TS^{-1}
$$


候補BHᵀS⁻¹を左辺に代入すると同じ連立方程式を満たす。

**要点：状態数nの逆行列から、観測数mの系へ計算を移せる。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OIの更新式とゲイン




$$
K=BH^T(HBH^T+R)^{-1}
$$





$$
x^a=x^b+K(y-Hx^b)
$$


B：n×n、H：m×n、R：m×m、K：n×m。

Kの第i行は、m個の残差から状態iをどれだけ直すかを表す。

**要点：104のQ更新も「ゲインのQ行×4成分の残差」という同じ構造。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OI導出⑤：更新後の誤差共分散




$$
e_a=(I-KH)e_b+Ke_o
$$





$$
P^a=(I-KH)B(I-KH)^T+KRK^T
$$





$$
P^a=B-BH^T(HBH^T+R)^{-1}HB
$$


無相関なので交差項が消える。第2式はJoseph形で、最適Kなら第3式へ簡約できる。

**要点：観測誤差Rの寄与を残すことが、解析の不確かさを正しく計算する要点。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## 「最適」の意味：線形不偏推定の誤差を最小化する




$$
P^a=B-KHB-BH^TK^T+KSK^T
$$





$$
\frac{\partial\,\mathrm{tr}(P^a)}{\partial K}=-2BH^T+2KS=0
$$





$$
K=BH^TS^{-1}
$$


共分散のトレースは、各状態成分の誤差分散の和。

**要点：ガウス性がなくても、ここで仮定した線形不偏推定の範囲で最適なKを導ける。**

ガウス性がない場合、これが全ての非線形推定を含むMMSE解とは限らない。

参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OIの数値例：温度予報20 ℃、観測22 ℃




$$
x^b=20,\quad B=4,\quad y=22,\quad R=1,\quad H=1
$$





$$
K=\frac{4}{4+1}=0.8,\qquad x^a=20+0.8(22-20)=21.6
$$





$$
P^a=(1-0.8)^2\,4+0.8^2\,1=0.8
$$


予報の標準偏差2 K、観測の標準偏差1 Kなので、観測に近い側へ補正する。

**要点：重み0.8は任意の調整値ではなく、背景と観測の誤差分散から決まる。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## OIは何を省略しているか

通常のOIではBを統計・相関モデルなどから与える。

- 背景の状態はモデルで時間発展させてもよい。
- ただしBそのものを、毎時刻のモデルで厳密に予報するわけではない。
- 場所・季節などでBを変える実装もある。

次のKFは、**誤差共分散の時間発展**も計算する。

**要点：OIとKFの解析式は共通。違いは主に背景・予報共分散の作り方にある。**



参考：[Bouttier & Courtier, Data Assimilation Concepts and Methods, §4–7 (ECMWF)](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)

## 第0部B　カルマンフィルタ（KF）

**状態の平均と誤差共分散を、時刻ごとに予報・更新する。**

線形モデルとガウス誤差なら、ベイズフィルタを平均と共分散だけで厳密に表せる。

**要点：OIの解析式に、状態と共分散の予報式を組み合わせる。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KFの線形状態空間モデル




$$
x_k=F_kx_{k-1}+b_k+\eta_k,\qquad y_k=H_kx_k+\nu_k
$$





$$
\eta_k\sim\mathcal{N}(0,\Omega_k),\quad\nu_k\sim\mathcal{N}(0,R_k)
$$


初期状態もガウス。ノイズは過去の状態推定誤差と独立、時刻間でも独立とする。
bkは既知入力による項。

**要点：未知の発熱量を推定したい場合は、入力を固定せず状態へ含める拡張が必要。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KF導出①：平均の予報




$$
\hat{x}_{k-1}^a=\mathrm{E}[x_{k-1}\mid y_{1:k-1}]
$$





$$
\hat{x}_k^f=\mathrm{E}[F_kx_{k-1}+b_k+\eta_k\mid y_{1:k-1}]
$$





$$
\hat{x}_k^f=F_k\hat{x}_{k-1}^a+b_k
$$


線形変換では期待値とFの積を入れ替えられ、平均0のηは消える。

**要点：観測を受け取る前に、前時刻の解析平均をモデルで進める。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KF導出②：共分散の予報




$$
e_k^f=F_ke_{k-1}^a-\eta_k
$$





$$
P_k^f=\mathrm{E}[(F_ke_{k-1}^a-\eta_k)(F_ke_{k-1}^a-\eta_k)^T]
$$





$$
P_k^f=F_kP_{k-1}^aF_k^T+\Omega_k
$$


独立・平均0のため交差項は0。既存誤差を伝播する項と、新たなモデル誤差の項が残る。

**要点：KFは状態だけでなく、不確かさが時間とともにどう変わるかも予報する。**

誤差を推定値−真値で定義する。真値には過程ノイズηが加わるため予報誤差には−ηが入るが、共分散の寄与は＋Ωとなる。

参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KF導出③：予報を背景としてOIの式を適用




$$
S_k=H_kP_k^fH_k^T+R_k,\qquad K_k=P_k^fH_k^TS_k^{-1}
$$





$$
\hat{x}_k^a=\hat{x}_k^f+K_k(y_k-H_k\hat{x}_k^f)
$$





$$
P_k^a=(I-K_kH_k)P_k^f(I-K_kH_k)^T+K_kR_kK_k^T
$$


OIのBを、その時刻に予報したPᶠへ置き換えた式。

**要点：解析値と解析共分散の両方を、次の時刻へ引き継ぐ。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KFの数値例①：1回目の予報と更新




$$
F=0.9,\ b=2,\ \Omega=0.16,\ H=1,\ R=1
$$





$$
\hat{x}_0^a=20,\ P_0^a=4\ \Longrightarrow\ \hat{x}_1^f=20,\ P_1^f=0.9^2\cdot4+0.16=3.4
$$





$$
y_1=22,\ K_1=\frac{3.4}{4.4}=0.772727
$$





$$
\hat{x}_1^a=21.545455,\qquad P_1^a=0.772727
$$


説明用の1温度モデル。Fとbは20 ℃へ緩和する時間発展を表す。

**要点：同じ22 ℃の観測でも、OI例と予報共分散が違えば重みも変わる。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## KFの数値例②：2回目は前回の結果から再開




$$
\hat{x}_2^f=0.9\cdot21.545455+2=21.390909
$$





$$
P_2^f=0.9^2\cdot0.772727+0.16=0.785909
$$





$$
y_2=21,\qquad K_2=\frac{0.785909}{1.785909}=0.440061
$$





$$
\hat{x}_2^a=21.218885,\qquad P_2^a=0.440061
$$


前回の同化で誤差が減ったため、今回のゲインも変わる。

**要点：同じRでもKは固定とは限らない。予報Pᶠに応じて変わる。**



参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## 非線形・高次元になると何が難しいか

非線形fでは、一般に **E[f(x)] ≠ f(E[x])**。

- 平均と共分散だけで分布の変化を閉じられない。
- 拡張KFはヤコビアンで線形化する近似。
- 104の20,697状態ではPに約4.28億要素。倍精度なら約3.43 GB。

EnKFは、各メンバーの前進と標本共分散を使う。

**要点：高次元の共分散を全要素で保持する代わりに、有限個の状態で近似する。**

メモリは共分散1枚の概算であり、行列積やソルバ全体の使用量ではない。

参考：[Kalman (1960); Särkkä & Svensson (2023)](https://doi.org/10.1115/1.3662552)

## 第0部C　アンサンブルカルマンフィルタ（EnKF）

**多数の状態候補を前進し、そのばらつきからKを計算する。**

非線形モデルは各メンバーで直接実行できる。

ただし、観測更新は共分散を使うカルマン型の近似である。

**要点：非線形モデルを前進できることと、任意の非ガウス事後分布を厳密に求めることは別。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## EnKF導出①：各メンバーを予報する




$$
x_k^{f,(i)}=f_k(x_{k-1}^{a,(i)})+\eta_k^{(i)},\qquad i=1,\ldots,N
$$





$$
\bar{x}_k^f=\frac{1}{N}\sum_i x_k^{f,(i)}
$$


予報分布をN個の状態候補で表す。
104の実ソルバ版は5ケースを個別に前進。明示的な全状態過程ノイズは加えず、初期ばらつき等を使う。

**要点：104の1メンバー＝1つのOpenFOAMケース。EnKF本体はPython側。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## EnKF導出②：KFのPを標本共分散で近似する




$$
a_x^{(i)}=x^{f,(i)}-\bar{x}^f
$$





$$
\widehat{P}^f=\frac{1}{N-1}\sum_i a_x^{(i)}(a_x^{(i)})^T
$$


平均を同じN標本から推定するので、偏差の和は0で自由度が1減る。

独立同分布の標本なら、この分母で共分散推定は不偏。フィルタ中のメンバーは更新により依存し得る。

**要点：N−1は標本共分散の正規化。有限メンバーによる推定誤差は残る。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## EnKF導出③：予報観測の偏差を作る




$$
y_f^{(i)}=h(x^{f,(i)}),\quad\bar{y}_f=\frac{1}{N}\sum_i y_f^{(i)},\quad a_y^{(i)}=y_f^{(i)}-\bar{y}_f
$$





$$
\widehat{C}_{xy}=\frac{1}{N-1}\sum_i a_x^{(i)}(a_y^{(i)})^T
$$





$$
\widehat{C}_{yy}=\frac{1}{N-1}\sum_i a_y^{(i)}(a_y^{(i)})^T
$$


104ではhが温度の抜き出しとFrontISTRによる変位評価を含む。

**要点：状態の全共分散を作らず、観測との交差共分散だけを計算できる。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## EnKF導出④：線形観測ならKFの式に一致する構造




$$
h(x)=Hx+c\ \Longrightarrow\ a_y^{(i)}=Ha_x^{(i)}
$$





$$
\widehat{C}_{xy}=\widehat{P}^fH^T,\qquad\widehat{C}_{yy}=H\widehat{P}^fH^T
$$





$$
\widehat{K}=\widehat{C}_{xy}(\widehat{C}_{yy}+R)^{-1}
$$


定数cは偏差を取ると消える。非線形hでは、同じ式を標本による近似として使う。

**要点：KFのPᶠを標本から作るとEnKFのゲインになる。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## EnKF導出⑤：確率的EnKFのメンバー更新




$$
\varepsilon^{(i)}\sim\mathcal{N}(0,R)
$$





$$
x^{a,(i)}=x^{f,(i)}+\widehat{K}(y+\varepsilon^{(i)}-y_f^{(i)})
$$


実観測yは全メンバー共通。各メンバーに独立の観測摂動εを与える。

これは観測を作る際の模擬測定ノイズとは別。

**要点：各メンバーを更新し、更新後の分布を次の予報へ残す。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## なぜ観測摂動を入れるのか：共分散の展開




$$
a_a^{(i)}=(I-KH)a_f^{(i)}+K(\varepsilon^{(i)}-\bar{\varepsilon})
$$





$$
\mathrm{E}[\widehat{P}^a\mid\mathrm{ensemble}]=(I-KH)\widehat{P}^f(I-KH)^T+KRK^T
$$


線形H・摂動が予報アンサンブルから独立という条件で、摂動について期待値を取る。

εを単に省くとKRKᵀが失われ、ばらつきを小さく見積もる。

**要点：有限Nの1回の更新では理論共分散と厳密一致しない。期待値・近似を区別する。**

決定論的平方根EnKFは、摂動を使わず別の変換で目標共分散を実現する。単にεを削除する方法とは違う。

参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## 理論の列ベクトルを、NumPyの行配列へ直す




$$
Z_f\in\mathbb{R}^{N\times n},\quad Y_f\in\mathbb{R}^{N\times m}
$$





$$
C_{zy}=\frac{dZ^TdY}{N-1},\qquad C_{yy}=\frac{dY^TdY}{N-1}
$$





$$
Z_a=Z_f+(\mathbf{1}y^T+E-Y_f)K^T
$$


Zfの1行が1メンバー。Eの1行が観測摂動。
104：N＝5、n＝20,697、m＝4。

**要点：式のKがコードではK.Tとして右から掛かるのは、行と列の持ち方が違うため。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## 低ランク・インフレーション・局所化




$$
\mathrm{rank}(dZ)\le N-1,\qquad a^{(i)}\leftarrow\lambda a^{(i)}
$$


- 5メンバーでは、表現できる偏差の方向は最大4。
- 偏差をλ倍すると共分散はλ²倍。104はλ＝1.05。
- 局所化は遠隔の見かけの相関を抑える工夫。今回の実ソルバ実装では使用していない。

**要点：高次元の場を更新していても、全ての独立した誤差を表現できるわけではない。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## 未知のQも状態に入れる：拡大状態




$$
z=[T_1,\ldots,T_{N_c},Q]^T
$$





$$
[C_{zy}]_{Q,:}=\frac{1}{N-1}\sum_i(Q^{(i)}-\bar Q)(y_f^{(i)}-\bar y_f)^T
$$





$$
K_Q=[C_{zy}]_{Q,:}(C_{yy}+R)^{-1}
$$





$$
Q_a^{(i)}=Q_f^{(i)}+K_Q(y+\varepsilon^{(i)}-y_f^{(i)})
$$


Qと観測の共分散が0なら、この更新式ではQを補正できない。

**要点：Qの更新方向は、ゲインのQ行と全観測残差の内積で決まる。**



参考：[Evensen (2003); 104/dacore/enkf.py](https://doi.org/10.1007/s10236-003-0036-9)

## 第0部D　粒子フィルタ（PF）

**状態候補に重みを付けて、事後分布を表す。**

EnKF：候補をゲインで移動する。
PF：候補の尤度から重みを変え、必要に応じて再標本化する。

**要点：ガウス分布に限定しない表現ができるが、有限粒子で十分に表せるとは限らない。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PF導出①：確率分布を重み付き粒子で近似




$$
p(x_k\mid y_{1:k})\approx\sum_{i=1}^{N}w_k^{(i)}\delta(x_k-x_k^{(i)})
$$





$$
w_k^{(i)}\ge0,\qquad\sum_iw_k^{(i)}=1
$$





$$
\hat{x}_k=\sum_iw_k^{(i)}x_k^{(i)}
$$


δは粒子位置に集中した確率質量を表す記号。重み付きの分布が主役。

**要点：重みが不均一なら、状態の平均も単純平均ではなく重み付き平均。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PF導出②：提案分布から粒子を生成する




$$
x_k^{(i)}\sim q_k(x_k\mid x_{k-1}^{(i)},y_k)
$$





$$
\widetilde w_k^{(i)}=w_{k-1}^{(i)}\frac{p(y_k\mid x_k^{(i)})p(x_k^{(i)}\mid x_{k-1}^{(i)})}{q_k(x_k^{(i)}\mid x_{k-1}^{(i)},y_k)}
$$


目標の「尤度×遷移密度」と、実際に標本を引いた提案密度の比で補正する。

これは祖先を同じ番号で引き継ぐ逐次重要度サンプリングの式。

**要点：一般のPFは、モデル予報からしか粒子を出せないわけではない。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PF導出③：ブートストラップPFでは遷移密度が消える




$$
q_k(x_k\mid x_{k-1},y_k)=p(x_k\mid x_{k-1})
$$





$$
\widetilde w_k^{(i)}=w_{k-1}^{(i)}p(y_k\mid x_k^{f,(i)})
$$





$$
w_k^{(i)}=\frac{\widetilde w_k^{(i)}}{\sum_j\widetilde w_k^{(j)}}
$$


モデルで前進し、観測の尤度を前の重みに掛ける。
前回に再標本化していれば、開始重みは1/N。

**要点：再標本化しない場合、前の重みを引き継ぐ必要がある。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PF導出④：ガウス観測誤差から尤度を作る




$$
r^{(i)}=y-h(x^{f,(i)})
$$





$$
L^{(i)}=\frac{\exp[-\tfrac12(r^{(i)})^TR^{-1}r^{(i)}]}{(2\pi)^{m/2}|R|^{1/2}}
$$





$$
\log\widetilde w^{(i)}=\log w_{prev}^{(i)}-\tfrac12(r^{(i)})^TR^{-1}r^{(i)}+c
$$


Rが全粒子で同じなら、定数cは正規化で消える。

**要点：ガウスなのはこの例の観測ノイズ。PFの状態事後分布をガウスに限定していない。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PFの数値例①：5つの温度候補に尤度を付ける

候補温度：[18, 19, 20, 21, 22] ℃。観測21.5 ℃、σ＝0.5 K。初期重みは1/5。


$$
\log L_i=c-\frac{(21.5-x_i)^2}{2\cdot0.5^2}
$$


```text
候補       18        19        20       21       22
log L−c  −24.5     −12.5      −4.5     −0.5     −0.5
重み      ≈0      0.000003   0.009075  0.495461  0.495461
```

観測を説明しやすい21 ℃・22 ℃へ、重みが集まる。

**要点：PFでは、この段階で粒子の温度そのものをゲインで変更していない。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PFの数値例②：重み付き平均とESS




$$
\hat{x}=\sum_iw_ix_i=21.486380\ ^\circ\mathrm{C}
$$





$$
N_{eff}=\frac{1}{\sum_iw_i^2}=2.036470
$$


ESSは有効サンプル数の指標。
等重みならN、1粒子に全重みが集まれば1。
閾値をN/2＝2.5とすれば、この例では再標本化する。

**要点：粒子が5個あっても、実質的には約2個に重みが集中している。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PF導出⑤：系統リサンプリングを具体的に追う




$$
c_j=\sum_{i=1}^{j}w_i,\quad u_i=u_0+\frac{i-1}{N},\quad u_0\sim U(0,1/N)
$$


```text
累積重み ≈ [0, 0.000003, 0.009078, 0.504539, 1]
u₀=0.08 とすると位置=[0.08, 0.28, 0.48, 0.68, 0.88]
選ぶ粒子温度       =[21,   21,   21,   22,   22]
新しい重み         =[1/5, 1/5, 1/5, 1/5, 1/5]
```

**要点：累積重みが各位置を初めて超える粒子を選び、重みを等しく戻す。**

この1回の再標本化後の平均は21.4℃。再標本化前の重み付き平均と毎回ぴったり同じになるわけではない。

参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## PFの弱点：重みの退化と粒子の重複

- 多数の観測に同時に合う粒子が少ないと、重みが一部へ集中する。
- 再標本化は重みを均等にするが、同じ粒子が複製される。
- モデル誤差や移動ステップがなければ、多様性が戻らない場合がある。
- ジッタは多様性を増やすが、事後分布も変えるため設定根拠が必要。

**要点：非ガウスを扱える利点と、高次元で必要粒子数が増える難しさを併せて理解する。**



参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 104のPF参考コードと、標準理論の差

`dacore/pf.py` は尤度・ESS・系統再標本化を実装している。

ただし、現状は **前時刻の重みを引数で受け取らない**。
`dacore/twin.py` も返った重みを平均計算に使わず、単純平均する。

再標本化しない時刻では、一般の逐次重要度更新を正しく表すには修正が必要。

**要点：PFの理論は解説するが、現コードをそのまま標準PFの検証済み結果とは扱わない。**

本作業は資料の改訂であり、PFソルバコードや過去の解析結果を変更しない。pf_updateは一様事前重みからの一回更新としては理解できる。

参考：[Särkkä & Svensson (2023); 104/dacore/pf.py（実装との相違は本文参照）](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 4手法の違いを、同じ軸で比較する

| 方法 | 不確かさの表現 | 観測で行う処理 |
|---|---|---|
| OI | 所与の背景共分散B | 共分散からKを作り補正 |
| KF | 予報する平均・共分散 | Kで平均と共分散を更新 |
| EnKF | 等重みメンバーの標本統計 | メンバーをKで移動 |
| PF | 粒子と重み | 尤度で重み更新・再標本化 |

**要点：どれも観測を使うが、不確かさの伝え方と更新の近似が異なる。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 「最適」「非線形」「非ガウス」の意味を整理する

- OI：与えた共分散と線形不偏推定の範囲で最小分散。
- KF：線形ガウスなら事後分布を平均・共分散で厳密に計算。
- EnKF：非線形モデルを使えても、更新は共分散に基づく近似。
- PF：非ガウス・多峰な分布を表現できるが、有限粒子の近似誤差がある。

どの方法でも、モデル・観測誤差の設定が間違えば推定は偏り得る。

**要点：手法の名前だけで精度が保証されるわけではない。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 104の説明を読む順番：式→行列→ファイル

1. EnKFの理論式で、どの共分散と残差を使うかを把握する。
2. 5メンバーの具体行列で、平均→偏差→共分散→K→Q更新を計算する。
3. `dacore/enkf.py` と `daof/of_fem_twin.py` で実装を確認する。
4. 実ソルバ結果とROMの別実験を、各条件の下で読む。

**要点：理論章はここまで。次は同じ式を104の物理量とコードに結び付ける。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

## 理論の参考資料と本資料の位置付け

- [OI：Bouttier & Courtier, ECMWF講義 §4–7](https://www.ecmwf.int/sites/default/files/elibrary/2002/16928-data-assimilation-concepts-and-methods.pdf)
- [KF：Kalman (1960)](https://doi.org/10.1115/1.3662552)
- [EnKF：Evensen (2003)](https://doi.org/10.1007/s10236-003-0036-9)
- [ベイズフィルタ・PF：Särkkä & Svensson (2023), 著者公開版](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)

式は本資料の記号に統一して展開。数値例は説明用で、104の解析結果とは区別する。

**要点：出典の文章や図を転載せず、導出・数値例・実装対応を組み合わせて説明した。**



参考：[Särkkä & Svensson, Bayesian Filtering and Smoothing (2023)](https://users.aalto.fi/~ssarkka/pub/bfs_book_2023_online.pdf)
