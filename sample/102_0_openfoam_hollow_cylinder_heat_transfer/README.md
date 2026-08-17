# 102_0: 中空円筒の過渡共役熱伝達解析

FreeCADで作った中空鋼製試験体を`snappyHexMesh`でメッシュ化し、
`chtMultiRegionFoam`で固体熱伝導と周囲空気の自然対流を連成するケースです。
外周の一部へ15 Wを300秒間入力し、その後300秒間冷却します。

## 解析条件

| 項目 | 値 |
|---|---:|
| 外径 / 内径 / 高さ | 75 / 40 / 100.5 mm |
| 実測質量 | 約2.5 kg |
| 固体 | 炭素鋼相当、rho=7850 kg/m3、cp=480 J/(kg K)、k=50 W/(m K) |
| 初期温度 | 293.15 K (20 degC) |
| ヒーター範囲 | 軸方向50 mm、周方向100 mm、+X側中心 |
| 入力 | 0～300秒: 15 W、300～600秒: 0 W |
| 放熱 | 周囲空気とのCHT（重力による自然対流）、放射なし |
| 保存間隔 | 5秒 |
| 背景格子 | 12.5 mm（ヒーター近傍はsnappyHexMeshで追加細分化） |
| snappyHexMesh | 特徴線level 3、試験体表面level (2 3)、ヒーター選択領域level 3 |

ヒーターは材料領域として作っていません。`heaterSelection.stl`と`topoSet`で
固体表面を選び、`createPatch`で作った`heaterPower`へ
`externalWallHeatFluxTemperature`の`mode power`で合計15 Wを与えます。
したがって、選択面の離散面積が変わっても合計入力は15 Wです。

現在のメッシュ検査値は、固体20,696セル、流体46,580セル、固体体積
`0.0003169756 m3`（幾何学値`0.0003177 m3`との差約0.23%）です。
ヒーターパッチは2,694面、面積`0.00507224 m2`で、15 W入力時の平均熱流束は
約`2,957 W/m2`です。全リージョンと固体リージョンの`checkMesh`はいずれも
`Mesh OK`です。

## フォルダ構成

| パス | 役割 |
|---|---|
| `cad/create_model.py` | `config/geometry.yaml`を読み、FreeCAD形状とSTL/STEPを生成 |
| `cad/generated/` | FCStd、STEP、FreeCADが出力したmm単位STL |
| `config/geometry.yaml` | 寸法、ヒーター位置、出力、背景領域、メッシュレベル |
| `constant/triSurface/` | OpenFOAMが使うm単位の`specimen.stl`と`heaterSelection.stl` |
| `system/snappyHexMeshDict` | 試験体表面へのスナップと局所細分化 |
| `system/topoSetDict` | 全セルから固体を除き、空気の`fluid` cellZoneを作成 |
| `system/*/heaterFaceSetDict` | ヒーター選択体と交差する界面をfaceSet化 |
| `system/*/createHeaterPatchDict` | solid側`heaterPower`、fluid側`heaterCover`を生成 |
| `system/solid/changeDictionaryDict` | 鋼材側CHT条件と15 W時刻表 |
| `0/`, `constant/`, `system/` | OpenFOAMの初期場、物性、数値・境界条件 |
| `VTK_heaterPower/` | ParaViewでヒーター選択面だけ確認するVTP |
| `postProcessing/` | 4点の固体温度プローブ（5秒間隔） |

## 実行方法

```bash
cd sample/102_0_openfoam_hollow_cylinder_heat_transfer

# CAD、snappyHexMesh、領域分割、ヒーターパッチ作成まで
./Allrun.pre

# 前処理済みなら再利用するか確認し、chtMultiRegionFoamまで実行
./Allrun
```

`Allrun`は既存の`constant/fluid/polyMesh`、`constant/solid/polyMesh`と
初期場を検出します。端末から実行した場合は前処理を再構築するか質問し、
既定の`N`ではCAD・メッシュを再利用します。無条件に再構築する場合だけ、
次を使います。

```bash
FORCE_PREPROCESS=1 FORCE_CAD_REBUILD=1 ./Allrun
```

`Allrun.pre`を直接実行した場合も、FreeCADデータが存在すればCADを作り直すか
質問します。非対話実行では既存CADを再利用します。

## 前処理の順序

```text
FreeCAD -> STLをmmからmへ変換 -> blockMesh -> snappyHexMesh
        -> solid/fluid cellZone -> splitMeshRegions
        -> heaterSelectionと界面の交差面を抽出
        -> heaterPower/heaterCover patch -> 境界条件適用
```

ParaViewでは`hollowCylinder.foam`で全領域を開きます。ヒーター面だけを
確認する場合は`VTK_heaterPower/solid/102_0_openfoam_hollow_cylinder_heat_transfer_0.vtm`
を開いてください。

## FrontISTRへの受け渡し

計算後、`../102_1_frontistr_hollow_cylinder_thermal_expansion`のPythonが各5秒時刻の
`solid/T`とセル中心`solid/C`を読み、逆距離補間でFrontISTR節点温度へ変換します。
OpenFOAMは温度と流れ、FrontISTRは各時刻の熱変位・応力を担当する一方向連成です。

## 計算結果（2026-08-11）

0～300秒を15 W加熱、300～600秒を0 W冷却として、5秒間隔の全121時刻を
`chtMultiRegionFoam`で計算し、正常終了しました。

- 固体最高温度: 298.9257 K（25.7757 degC、300秒）
- ヒーター側プローブ最高温度: 25.5805 degC（300秒）
- 反ヒーター側プローブ最高温度: 23.5949 degC（530秒）
- 600秒の固体温度範囲: 296.7030～296.7580 K
- `data/temperature_history.csv`: 0～600秒、121行
- `data/temperature_history.png`: 4測定点の温度時刻歴
- `data/summary.yaml`: メッシュ検査値と主要温度結果
- `postProcessing/solidTemperatureProbes/solid/0/T`: OpenFOAMの生プローブ出力

## アニメーション

OpenFOAM温度・流れとFrontISTR熱変形を並べた121枚のPNG連番から、公開用GIFを
生成済みです。

- `ani/ani.gif`: 900 x 460 px、121フレーム、約2.8 MB（GitHub公開用に減色圧縮済み）
- `ani/make_gif.sh`: ImageMagick `convert`による再生成スクリプト

```bash
cd ani
./make_gif.sh
```

PNG連番は約17 MBあるためGit対象外とし、圧縮済み`ani.gif`だけを公開対象にします。

反対側の温度ピークが加熱停止より遅れるのは、鋼材内を熱が拡散する時間遅れによります。

## 仮定と今後の段階

- 内周ねじ山は滑らかな直径40 mm穴へ簡略化しています。
- ヒーター厚さ1.5 mmはCAD表示用の仮定で、基準解析の熱容量には含めません。
- 現在は15 Wがすべて試験体へ入る簡易モデルです。
- 放射は未使用です。次段階で黒皮0.8、切削面0.3～0.5を比較します。
- 接触熱抵抗モデルは、ヒーターパッドの物性・厚さ・接触状態を測定してから追加します。
- 実験同定では、実測電力、室温、ヒーター貼付角度、接触材、温度測定座標が必要です。
