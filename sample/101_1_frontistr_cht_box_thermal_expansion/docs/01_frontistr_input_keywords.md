# FrontISTR熱膨張ケースの入力ファイルとキーワード

## 1. 3つの入力ファイル

| ファイル | 役割 |
|---|---|
| `hecmw_ctrl.dat` | 使用するメッシュ、解析制御、結果ファイル名を結び付ける |
| `box_thermal_expansion.msh` | 節点、要素、節点群、要素群、材料、断面を定義する |
| `box_thermal_expansion.cnt` | 解析種別、拘束、節点温度、材料物性、ソルバー、出力を定義する |

Pythonが各時刻にこれらを生成します。初めは`case/t_300/`など、実際に生成された
ファイルをこの説明と並べて読むと理解しやすくなります。

## 2. `hecmw_ctrl.dat`

```text
!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE
box_thermal_expansion.msh
```

- `!MESH`: 読み込むメッシュを登録します。
- `NAME=fstrMSH`: FrontISTRが既定で参照するメッシュ名です。
- `TYPE=HECMW-ENTIRE`: 1ファイルに全メッシュが入るHEC-MW形式です。

```text
!CONTROL, NAME=fstrCNT
box_thermal_expansion.cnt
```

解析条件ファイルを`fstrCNT`という既定名へ割り当てます。

```text
!RESULT, NAME=fstrRES, IO=OUT
box_thermal_expansion.res
```

節点変位や応力などの解析結果の出力名を指定します。

```text
!RESULT, NAME=vis_out, IO=OUT
box_thermal_expansion_vis
```

ParaView可視化用結果の出力名を指定します。

## 3. `.msh`の主要キーワード

### `!NODE`

```text
!NODE
1, 0.0, 0.0, 0.0
2, 0.05, 0.0, 0.0
```

`節点ID, X, Y, Z`の順です。本ケースの長さ単位はmです。

### `!ELEMENT, TYPE=361`

```text
!ELEMENT, TYPE=361
1, 1, 2, 7, 6, 26, 27, 32, 31
```

要素IDと接続節点を指定します。`TYPE=361`は8節点一次六面体要素です。
節点順序を誤ると要素の体積が負になったり、形がねじれたりします。

### `!NGROUP`

```text
!NGROUP, NGRP=BOTTOM
1, 2, 3, ...
```

節点集合を作ります。`BOTTOM`はz=0の節点で、固定境界に使用します。
`NALL`は全節点、`TOP`は上面節点です。

### `!EGROUP`

```text
!EGROUP, EGRP=EALL
1, 2, 3, ...
```

要素集合です。`EALL`へ全要素を登録し、材料と断面を一括適用します。

### `!MATERIAL`と`!ITEM`

```text
!MATERIAL, NAME=STEEL, ITEM=3
!ITEM=1, SUBITEM=2
2.05e11, 0.3
!ITEM=2, SUBITEM=1
7850
!ITEM=3, SUBITEM=1
1.2e-5
```

- `ITEM=3`: 3種類の材料項目を持つ定義です。
- `ITEM=1`: ヤング率E[Pa]とポアソン比νです。
- `ITEM=2`: 密度[kg/m3]です。
- `ITEM=3`: 線膨張係数[1/K]です。

本ケースではFrontISTRのメッシュ読込時に材料名を解決できるよう、`.msh`にも
材料を記載します。同じ値を`.cnt`にも記載するため、両者を必ず一致させます。

### `!SECTION`

```text
!SECTION, TYPE=SOLID, EGRP=EALL, MATERIAL=STEEL
1.0
```

全要素`EALL`を3次元ソリッドとし、材料`STEEL`を割り当てます。次行の`1.0`は
ソリッド断面の係数です。

## 4. `.cnt`の主要キーワード

### `!SOLUTION, TYPE=STATIC`

線形静解析を指定します。時刻ごとに温度荷重を変えてこの静解析を繰り返します。
OpenFOAMの時間値はケース名と集計値に使われ、FrontISTR内部の動的時間積分には
使われません。

### `!SOLVER`

```text
!SOLVER,METHOD=CG,PRECOND=1,NSET=0,ITERLOG=NO,TIMELOG=NO
5000, 1
1.0e-08, 1.00, 0.0
```

- `METHOD=CG`: 共役勾配法で連立一次方程式を解きます。
- `PRECOND=1`: 前処理を使用します。
- `5000`: 最大反復回数です。
- `1.0e-08`: 収束判定値です。
- `ITERLOG=NO`: 反復履歴の詳細表示を抑えます。

### `!REFTEMP`

```text
!REFTEMP
293.15
```

熱ひずみがゼロとなる基準温度です。本ケースではOpenFOAM初期温度20℃と同じ
293.15 Kです。節点温度`T`に対する自由熱ひずみは概ね次式です。

$$
\boldsymbol{\varepsilon}_{\mathrm{th}}
= \alpha\left(T-T_{\mathrm{ref}}\right)\boldsymbol{I}
$$

ここで$\alpha$は線膨張係数、$T$は節点温度、$T_{\mathrm{ref}}$は
`!REFTEMP`、$\boldsymbol{I}$は単位テンソルです。

### `!INITIAL_CONDITION, TYPE=TEMPERATURE`

```text
!INITIAL_CONDITION, TYPE=TEMPERATURE
NALL, 293.15
```

全節点の初期温度を設定します。その後の`!TEMPERATURE`で各節点値を上書きします。

### `!BOUNDARY`

```text
!BOUNDARY, GRPID=1
BOTTOM,1,1
BOTTOM,2,2
BOTTOM,3,3
```

`節点群, 最初の自由度, 最後の自由度`です。自由度1、2、3はX、Y、Z変位です。
したがって底面をXYZ全方向に固定します。拘束が足りないと剛体移動で解けず、
拘束しすぎると実物より大きな熱応力が出るため、実際の支持方法に合わせます。

### `!TEMPERATURE`

```text
!TEMPERATURE, GRPID=1
1, 293.15
2, 293.18
...
```

`節点ID, 温度[K]`です。PythonがOpenFOAMセル温度から補間して全節点分を書きます。

### `!ELASTIC`

```text
!MATERIAL, NAME=STEEL
!ELASTIC
2.05e11, 0.3
```

ヤング率E[Pa]とポアソン比νを指定します。温度による物性変化はこの設定では
考慮せず、全温度範囲で一定とします。

### `!EXPANSION_COEFF`

```text
!EXPANSION_COEFF
1.2e-5
```

線膨張係数α[1/K]です。熱変形量へ直接効くため、実材料の値へ更新します。

### `!STEP`

```text
!STEP, SUBSTEPS=1, CONVERG=1.000E-07
BOUNDARY,1
LOAD,1
```

- `SUBSTEPS=1`: 温度荷重を1段階で与えます。線形静解析なので1で十分です。
- `BOUNDARY,1`: `GRPID=1`の拘束をこのステップで有効にします。
- `LOAD,1`: `GRPID=1`の温度荷重を有効にします。

### `!VISUAL`

```text
!VISUAL, method=PSR
!surface_num=1
!surface 1
!output_type = VTK
```

表面抽出結果をVTK形式で出力します。生成された`.pvtu`をParaViewで開きます。

本ケースのFrontISTR標準VTKには`DISPLACEMENT`、`NodalSTRESS`、`NodalMISES`が
出力されますが、熱荷重の`TEMPERATURE`は含まれません。そのため
`fistr_case.add_temperature_to_visualization()`が、解析後に同じ節点順序を座標で
確認して`TEMPERATURE`をPVTU/VTUへ追加します。同時にVTK内の`TimeValue`を
OpenFOAMの実時刻へ置き換え、全時刻をPVDで束ねます。

## 5. 単位系

FrontISTRは単位を自動変換しないため、入力全体で整合させます。本ケースはSI単位です。

| 量 | 単位 |
|---|---|
| 座標・変位 | m |
| ヤング率 | Pa |
| 密度 | kg/m3 |
| 温度 | K |
| 線膨張係数 | 1/K |

`timehistory.csv`では読みやすいように変位だけmからmmへ変換しています。

## 6. エラー時に見るファイル

| ファイル | 確認内容 |
|---|---|
| `log.fistr1` | 起動エラー、入力エラー、反復計算の終了状態 |
| `FSTR.msg` | FrontISTRの詳細メッセージ |
| `FSTR.sta` | ステップの進行状況 |
| `summary.yaml` | 入力温度範囲、座標移動量、変位要約 |

まず`log.fistr1`末尾を確認します。

```bash
tail -n 80 case/t_300/log.fistr1
```

結果が不自然な場合は、エラーの有無だけでなく、`summary.yaml`の基準温度、
節点温度範囲、座標移動量、固定面を確認します。
