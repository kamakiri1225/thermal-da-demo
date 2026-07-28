# OpenFOAM温度時刻歴からFrontISTR熱膨張解析を行う手順

## 1. この連成で行うこと

`101_0_openfoam_cht_radiation_box`で計算した固体温度を、同じ保存時刻ごとに
`101_1_frontistr_cht_box_thermal_expansion`へ渡し、熱膨張変形を求めます。

```text
OpenFOAM
  0, 5, 10, ... 600秒の固体セル温度 T
        ↓ セル中心座標と温度を読み取る
Python
  座標原点を合わせ、FrontISTR節点へ温度を補間
        ↓ !TEMPERATUREとして.cntへ書く
FrontISTR
  各時刻を底面固定の線形静解析として計算
        ↓
  時刻別VTK、変位結果、timehistory.csv、timehistory.png
```

これは**OpenFOAMからFrontISTRへの一方向連成**です。温度は変形へ影響しますが、
変形後の形状をOpenFOAMへ戻しません。またFrontISTRでは各時刻を独立した静解析と
して解くため、慣性、振動、塑性履歴、接触履歴は考慮しません。加熱が構造応答より
十分ゆっくり進む熱膨張を対象とした、準静的な時刻歴です。

## 2. 使用するフォルダ

```text
sample/
├── 101_0_openfoam_cht_radiation_box/
│   ├── 0/, 5/, 10/, ... 600/    OpenFOAMの保存時刻
│   ├── constant/                 物性、リージョンメッシュ
│   ├── system/                   計算条件
│   ├── Allrun.pre                メッシュとリージョンの準備
│   └── Allrun                    準備と熱流体計算
└── 101_1_frontistr_cht_box_thermal_expansion/
    ├── config/                   構造材料物性
    ├── python/                   温度変換とFrontISTR連続実行
    ├── case/t_<時刻>/            時刻ごとのFrontISTR入力・結果
    ├── data/                     時刻歴CSVとグラフ
    └── docs/                     この解説
```

## 3. 事前確認

### 3.1 OpenFOAMを有効にする

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc
```

`source`はOpenFOAMの実行ファイル、ライブラリ、環境変数を現在のシェルへ登録します。
インストール場所が違う場合は実環境の`etc/bashrc`を指定します。

```bash
which chtMultiRegionFoam
which postProcess
foamVersion
```

- `which chtMultiRegionFoam`: CHTソルバーの場所を確認します。
- `which postProcess`: セル中心座標を書き出す後処理コマンドを確認します。
- `foamVersion`: 使用中のOpenFOAM版を確認します。

### 3.2 FrontISTRを確認する

```bash
export FISTR1=/home/kamakiri/local/frontistr/bin/fistr1
"$FISTR1" --help
```

`FISTR1`は本ケースのPythonが参照する実行ファイルの環境変数です。`--help`を
受け付けない版でも、ファイルが存在し実行権限があれば解析時に使用できます。

```bash
test -x "$FISTR1" && echo "FrontISTR executable: OK"
```

## 4. OpenFOAMで温度時刻歴を計算する

### 4.1 ケースへ移動する

```bash
cd sample/101_0_openfoam_cht_radiation_box
```

以降のOpenFOAMコマンドが、このケースの`system/`、`constant/`、`0/`を読むように
作業場所を変更します。

### 4.2 計算時間と保存間隔を確認する

`system/controlDict`の主な設定は次のとおりです。

```text
application     chtMultiRegionFoam;
endTime         600;
writeControl    adjustableRunTime;
writeInterval   5;
adjustTimeStep  yes;
```

- `application`: 使用するソルバーです。
- `endTime 600`: 600秒まで計算します。
- `writeInterval 5`: 5秒ごとに結果を保存します。
- `adjustTimeStep yes`: Courant数等を満たすよう内部時間刻みを自動調整します。

FrontISTRの時刻列はここで決まった保存フォルダをそのまま使用します。つまり
`writeInterval`を10秒へ変えれば、FrontISTRも`0, 10, 20, ...`を処理します。

### 4.3 メッシュとリージョンを準備する

```bash
./Allrun.pre
```

`Allrun.pre`は内部で次を順番に実行します。

```bash
blockMesh
```

`system/blockMeshDict`から背景の六面体メッシュを作ります。

```bash
topoSet
```

`system/topoSetDict`に従い、流体、固体、ヒートマットに使うセル集合を作ります。

```bash
splitMeshRegions -cellZones -overwrite
```

セル集合を`fluid`、`solid`、`heaterMat`の独立リージョンへ分割します。
`-overwrite`は同じケース内へリージョンメッシュを書きます。

```bash
changeDictionary -region <リージョン名>
```

各`system/<region>/changeDictionaryDict`の内容を初期境界条件へ反映します。
実際には`Allrun.pre`が全リージョンに対して繰り返します。

### 4.4 CHT計算を実行する

準備を含めて最初から実行する場合は次の1コマンドです。

```bash
./Allrun
```

`Allrun`は`Allrun.pre`を実行後、`system/controlDict`の`application`に書かれた
`chtMultiRegionFoam`を実行します。既に準備済みで、続きを計算する場合は次でも構いません。

```bash
chtMultiRegionFoam | tee log.chtMultiRegionFoam
```

`tee`は画面表示とログ保存を同時に行います。正常終了はログ末尾の`End`で確認します。

```bash
tail -n 30 log.chtMultiRegionFoam
```

### 4.5 保存時刻を確認する

```bash
foamListTimes
```

本ケースの既定値では`0, 5, 10, ... 600`の121時刻です。この一覧がそのまま
FrontISTRの解析時刻になります。

## 5. OpenFOAM温度をFrontISTRへ渡す

### 5.1 なぜセル中心座標Cが必要か

OpenFOAMの`T`はセルごとの値ですが、FrontISTRの`!TEMPERATURE`は節点ごとの値です。
どの温度が空間上のどこにあるかを知るため、OpenFOAMセル中心座標`C`も読みます。

時刻歴スクリプトは不足している`C`を自動生成します。手動で確認したい場合は次を実行します。

```bash
postProcess -func writeCellCentres -region solid -time '0:'
postProcess -func writeCellCentres -region heaterMat -time '0:'
```

- `-func writeCellCentres`: セル中心の`C`、`Cx`、`Cy`、`Cz`を書きます。
- `-region solid`: 固体リージョンを対象にします。
- `-time '0:'`: 0秒以降の全保存時刻を対象にします。シェル展開を防ぐため引用符を付けます。

### 5.2 座標系を合わせる

OpenFOAMの固体は既定で`x=0.4--0.6m, y=0.4--0.6m, z=0--0.4m`です。
FrontISTRモデルは`x=0--0.2m, y=0--0.2m, z=0--0.4m`です。寸法と軸は同じですが
原点が違うため、セル中心点群と節点群の中心を合わせます。既定の平行移動量は
`(-0.4, -0.4, 0)m`です。回転や拡大縮小は行いません。

OpenFOAMセル中心群のバウンディングボックス中心を $\boldsymbol{c}_{\mathrm{OF}}$、
FrontISTR節点群のバウンディングボックス中心を $\boldsymbol{c}_{\mathrm{FEM}}$ とすると、
平行移動量と移動後座標は次式です。

$$
\begin{aligned}
\boldsymbol{c}_{\mathrm{OF}}
&= \frac{\boldsymbol{x}_{\mathrm{OF,min}}+\boldsymbol{x}_{\mathrm{OF,max}}}{2}, \\
\boldsymbol{c}_{\mathrm{FEM}}
&= \frac{\boldsymbol{x}_{\mathrm{FEM,min}}+\boldsymbol{x}_{\mathrm{FEM,max}}}{2}, \\
\Delta\boldsymbol{x}
&= \boldsymbol{c}_{\mathrm{FEM}}-\boldsymbol{c}_{\mathrm{OF}}, \\
\boldsymbol{x}'_i
&= \boldsymbol{x}_i+\Delta\boldsymbol{x}.
\end{aligned}
$$

### 5.3 セル温度を節点温度へ補間する

各FrontISTR節点について、近いOpenFOAMセル中心8点を探し、距離の逆数を重みとして
平均します。

FrontISTR節点 $n$ の座標を $\boldsymbol{x}_n$、近傍OpenFOAMセル $i$ の移動後
セル中心座標を $\boldsymbol{x}'_i$、セル温度を $T_i$ とすると、距離 $d_{ni}$、
重み $w_{ni}$、節点温度 $T_n$ は次式です。

$$
\begin{aligned}
d_{ni} &= \left\lVert \boldsymbol{x}_n-\boldsymbol{x}'_i \right\rVert_2, \\
w_{ni} &= \frac{1}{d_{ni}}, \\
T_n &= \frac{\displaystyle\sum_{i\in\mathcal{N}_8(n)}w_{ni}T_i}
{\displaystyle\sum_{i\in\mathcal{N}_8(n)}w_{ni}}.
\end{aligned}
$$

$\mathcal{N}_8(n)$は節点$n$に近い8セルの集合です。節点とセル中心が一致して
$d_{ni}<10^{-9}\,\mathrm{m}$となる場合は、ゼロ除算を避け、そのセル温度を直接使います。

これにより`solid`と`heaterMat`のセル温度128個からFrontISTR節点225個の温度を作ります。
温度ピークは平均で少し滑らかになるため、局所温度勾配が重要な場合は補間法と
メッシュ解像度の検証が必要です。

## 6. 全時刻のFrontISTR解析を実行する

```bash
cd ../101_1_frontistr_cht_box_thermal_expansion/python
export FISTR1=/home/kamakiri/local/frontistr/bin/fistr1
python3 run_thermal_expansion_timehistory.py \
  --of-case ../../101_0_openfoam_cht_radiation_box
```

コマンドの意味は次のとおりです。

- `python3`: Python 3で連携スクリプトを実行します。
- `run_thermal_expansion_timehistory.py`: OpenFOAM保存時刻を列挙し、時刻ごとに入力作成、`fistr1`実行、結果集計を行います。
- `--of-case`: 温度を読むOpenFOAMケースを指定します。

既定では0秒を含む全保存時刻を計算します。絞り込み例は次のとおりです。

```bash
# 100～300秒だけ
python3 run_thermal_expansion_timehistory.py --start-time 100 --end-time 300

# 保存時刻を2個おきに処理
python3 run_thermal_expansion_timehistory.py --every 2

# 初期状態0秒を除外
python3 run_thermal_expansion_timehistory.py --exclude-t0
```

単一時刻だけ確認する場合は次です。

```bash
python3 run_thermal_expansion.py \
  --of-case ../../101_0_openfoam_cht_radiation_box \
  --time 300
```

## 7. スクリプトが時刻ごとに行う処理

1. `<時刻>/solid/T`と`<時刻>/heaterMat/T`を読む。
2. 同じ場所の`C`を読み、温度と座標をセル順で対応させる。
3. OpenFOAM座標をFrontISTR座標へ平行移動する。
4. FrontISTRの六面体メッシュを作る。
5. セル温度を節点温度へ補間する。
6. `.msh`、`.cnt`、`hecmw_ctrl.dat`を書く。
7. その時刻のフォルダ内で`fistr1`を実行する。
8. FrontISTRのVTKへ節点温度`TEMPERATURE`とOpenFOAM時刻を追加する。
9. 変位を読み、最大変位、上面Z変位、反り指標を集計する。

### 7.1 処理とプログラムの対応

| 処理 | ファイル | 関数 |
|---|---|---|
| OpenFOAMの`C`と`T`を読む | [`openfoam_temperature.py`](../python/openfoam_temperature.py) | `_load_region_cell_temperatures()` |
| `solid`と`heaterMat`を結合 | [`openfoam_temperature.py`](../python/openfoam_temperature.py) | `load_solid_cell_temperatures()` |
| 座標原点を合わせる | [`openfoam_temperature.py`](../python/openfoam_temperature.py) | `align_cell_centers_to_node_mesh()` |
| 逆距離加重補間 | [`openfoam_temperature.py`](../python/openfoam_temperature.py) | `interpolate_to_nodes()` |
| 1時刻分を接続する | [`run_thermal_expansion.py`](../python/run_thermal_expansion.py) | `run_one_time()` |
| `!TEMPERATURE`を書く | [`fistr_case.py`](../python/fistr_case.py) | `write_cnt()` |
| 全保存時刻を繰り返す | [`run_thermal_expansion_timehistory.py`](../python/run_thermal_expansion_timehistory.py) | `main()` |

### 7.2 OpenFOAMの座標と温度を読むコード

`openfoam_temperature.py`では、同じ時刻・同じリージョンの`C`と`T`を読みます。
OpenFOAMは両フィールドを同じセル順で保存するため、配列の同じ行が同じセルです。

```python
def _load_region_cell_temperatures(of_case_dir, time_dir, region):
    t_path = Path(of_case_dir) / time_dir / region / "T"
    c_path = Path(of_case_dir) / time_dir / region / "C"
    centers = _read_of_vector_field(c_path)
    temps = _read_of_scalar_field(t_path)
    if len(temps) == 1 and len(centers) > 1:
        temps = np.full(len(centers), temps[0])
    if len(temps) != len(centers):
        raise ValueError("TとCのセル数が不一致")
    return centers, temps
```

### 7.3 座標原点を合わせるコード

上記5.2節の$\Delta\boldsymbol{x}$を計算している部分です。

```python
source_center = 0.5 * (cell_centers.min(axis=0) + cell_centers.max(axis=0))
target_center = 0.5 * (node_coords.min(axis=0) + node_coords.max(axis=0))
translation = target_center - source_center
aligned_centers = cell_centers + translation
```

### 7.4 逆距離加重補間を行うコード

上記5.3節の$d_{ni}$、$w_{ni}$、$T_n$に対応します。

```python
for idx in range(n_nodes):
    d = np.linalg.norm(cell_centers - node_coords[idx], axis=1)
    nearest = np.argsort(d)[:k_eff]
    dn = d[nearest]
    if dn[0] < 1e-9:
        result[idx] = cell_temperatures[nearest[0]]
        continue
    w = 1.0 / dn
    result[idx] = np.sum(w * cell_temperatures[nearest]) / np.sum(w)
```

### 7.5 FrontISTRの節点温度を書くコード

補間した`node_temps`を、節点IDと組にして`.cnt`の`!TEMPERATURE`へ書きます。

```python
lines.append("!TEMPERATURE, GRPID=1")
for node_id, temperature in zip(node_ids, node_temperatures):
    lines.append(f" {node_id}, {temperature:.10g}")
```

生成結果は例えば次の形式です。

```text
!TEMPERATURE, GRPID=1
1, 293.157548
2, 293.160214
...
```

### 7.6 全時刻を繰り返すコード

OpenFOAMの数値名フォルダを時刻順に並べ、`run_one_time()`を1回ずつ呼びます。

```python
times = list_time_dirs(of_case)
for time_dir in times:
    case_dir = case_root / f"t_{time_dir}"
    summary = run_one_time(
        of_case, time_dir, case_dir, material, nx, ny, nz
    )
```

したがって、FrontISTR側へ別の時間間隔をハードコードしていません。OpenFOAMに
存在する保存時刻が、そのままFrontISTR解析とPVDの時刻になります。

## 8. 出力を確認する

```text
case/t_300/
├── box_thermal_expansion.msh
├── box_thermal_expansion.cnt
├── hecmw_ctrl.dat
├── box_thermal_expansion.res.0.1
├── box_thermal_expansion_vis_psf.0001.pvtu
├── box_thermal_expansion_vis_psf.0001/*.vtu
├── log.fistr1
└── summary.yaml
```

全時刻の一覧は次です。

```bash
column -s, -t < ../data/timehistory.csv | less -S
```

グラフは`data/timehistory.png`です。ParaViewで全時刻を連続表示する場合は、
OpenFOAMの時刻値を登録した`data/thermal_expansion_timehistory.pvd`を開きます。
`DISPLACEMENT`を選択して`Warp By Vector`を適用すると熱変形を確認できます。

```bash
paraview ../data/thermal_expansion_timehistory.pvd
```

### ParaViewで連番温度分布を見る操作

1. `thermal_expansion_timehistory.pvd`を開き、Propertiesの`Apply`を押します。
2. 上部の時刻操作ボタンで再生します。時刻はOpenFOAMと同じ`0, 5, ... 600s`です。
3. Coloringのプルダウンで`TEMPERATURE`を選ぶと節点温度[K]を表示します。
4. ℃表示が必要なら`Calculator`を追加し、式を`TEMPERATURE - 273.15`、結果名を`Temperature_degC`にします。
5. 変形を見る場合は`Warp By Vector`を追加し、Vectorsに`DISPLACEMENT`を選びます。
6. 本結果は変位がマイクロメートル級なので、形状確認用にはScale Factorを大きくします。数値評価には拡大前の値を使います。

同じPVD内で次の配列を選択できます。

| 配列 | 内容 |
|---|---|
| `TEMPERATURE` | OpenFOAMから写像したFrontISTR節点温度[K] |
| `DISPLACEMENT` | FrontISTRで計算した節点変位[m] |
| `NodalSTRESS` | 節点応力成分[Pa] |
| `NodalMISES` | ミーゼス相当応力[Pa] |

FrontISTRは静解析の温度荷重を標準VTKへ出力しないため、連携スクリプトが解析後の
PVTU/VTUへ`TEMPERATURE`を追加しています。温度値は`.cnt`の`!TEMPERATURE`と同じです。

## 9. 条件を変更するときの確認表

| 変更内容 | OpenFOAM | FrontISTR |
|---|---|---|
| 固体寸法 | `system/include/caseSettings` | `run_thermal_expansion.py`の寸法定数 |
| セル寸法 | `cellSize` | `--nx --ny --nz` |
| 初期温度 | 各リージョンの`0/*/T` | 材料JSONの`reference_temperature_K` |
| 材料 | 熱伝導率、密度、比熱 | E、ν、密度、線膨張係数 |
| 固定方法 | 構造条件は持たない | `.cnt`の`!BOUNDARY` |
| 保存時刻 | `controlDict` | 自動的に同じ時刻フォルダを使用 |

寸法、軸方向、単位系が一致しないまま温度を渡すと、計算は終了しても物理的に誤った
結果になります。実材料を使う前に、基準温度、線膨張係数、固定条件を必ず確認します。
