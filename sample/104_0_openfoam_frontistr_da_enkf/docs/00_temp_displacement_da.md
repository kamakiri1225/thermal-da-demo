# 温度＋変位観測のデータ同化 ― 仕組みと数式

## 1. 何が新しいか(103 との差分)

103 の観測は温度のみだった:

$$y = [\, T_\mathrm{hot},\ T_\mathrm{cold} \,]$$

104 は **FrontISTR 熱膨張解析の変位**を観測に加える:

$$y = [\, T_\mathrm{hot},\ T_\mathrm{cold},\ u_z^\mathrm{heater},\ u_z^\mathrm{opposite} \,]$$

これは実験(熱電対＋ダイヤルゲージ/変位センサ)と同じ観測構成であり、
「温度計が少なくても、変位計測が温度場の推定を助ける」ことを確かめる設定になっている。

## 2. 非線形観測演算子

温度観測は状態(温度場)の抜き出しなので線形だが、変位は

$$u = g(T) \quad \text{(温度場 → IDW補間 → FrontISTR線形静解析 → 上面 } U_z)$$

という**物理チェーンを通した非線形写像**。EnKF ではヤコビアンを作る必要はなく、
各メンバー $i$ について実際にチェーンを評価した予報観測

$$y_f^{(i)} = h\!\left(x_f^{(i)}\right) = \left[\, T^{(i)}_\mathrm{obs},\ g\!\left(T^{(i)}\right) \right]$$

を並べ、標本共分散で更新すればよい:

$$C_{zy} = \frac{1}{N-1}\sum_i \left(x_f^{(i)}-\bar{x}_f\right)\left(y_f^{(i)}-\bar{y}_f\right)^\top$$

$$x_a^{(i)} = x_f^{(i)} + C_{zy}\left(C_{yy}+R\right)^{-1}\left(y+\varepsilon^{(i)}-y_f^{(i)}\right)$$

実装: `dacore/enkf.py enkf_update(..., Yf=Yf)`(Yf サンプリング対応)、
`daof/of_fem_twin.py member_obs()`(チェーン評価)、
`fem/fem_obs.py displacement_obs()`(102_1 の `run_one_time` を再利用)。

## 3. なぜ変位観測が効くのか

熱膨張変位は温度場の**重み付き積分**(体積的な情報)なので、点measurementの温度と
相補的な情報を持つ:

- でたらめに熱い/冷たいメンバーは、全体が伸び/縮みして上面変位が大きくズレる
  → 変位1点で温度場全体のオフセットを強く拘束
- ヒータ側と反対側の変位差(傾き)は温度の**非対称分布**を反映
  → Q の位置・大きさの情報

例(この設定の実測値): 真値 t=10 s の上面変位は heater 側 +0.36 µm / 反対側 −0.32 µm。
一方 +27 K のでたらめメンバーは一様膨張で ~+32 µm となり、0.1 µm の観測ノイズに対して
圧倒的に識別できる。

## 4. 観測誤差 R の単位混在

$y$ は K と mm が混在するが、$R$ を対角に

$$R = \mathrm{diag}\!\left(\sigma_T^2,\ \sigma_T^2,\ \sigma_u^2,\ \sigma_u^2\right), \quad \sigma_T = 0.3\ \mathrm{K},\ \sigma_u = 10^{-4}\ \mathrm{mm}$$

と与えれば、イノベーションが各 $\sigma$ で正規化されるため単位の違いは問題にならない。

## 5. 実行手順・規模

README のとおり。規模は 103 と同じ最小実行(N=5, 0→20 s, 2サイクル)。
FrontISTR は1評価あたり十数秒で、OpenFOAM(1ウィンドウ数分)に比べ軽い。

コスト内訳(1サイクル):

$$\underbrace{N \times \text{chtMultiRegionFoam}}_{\text{数分} \times 5} + \underbrace{(N+1) \times \text{FrontISTR}}_{\text{約16秒} \times 6}$$
