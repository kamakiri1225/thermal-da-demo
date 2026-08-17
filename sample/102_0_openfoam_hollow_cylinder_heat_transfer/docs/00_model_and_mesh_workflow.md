# 102_0 モデル・メッシュ・熱入力の作成手順

## 1. 形状

`cad/create_model.py`が`config/geometry.yaml`を読み、FreeCADで次の形状を作ります。

- `Specimen`: 外径75 mm、内径40 mm、高さ100.5 mmの中空円筒
- `HeaterMat`: 実験配置を確認する50 x 100 mmの曲面ヒーター表示形状
- `HeaterSelection`: OpenFOAMで加熱面を選ぶため、表面へ0.5 mm食い込ませた選択体

出力は`cad/generated/hollow_cylinder_heater.FCStd`、STEP、mm単位STLです。
FreeCADモデルが既にあれば`Allrun.pre`は再作成の要否を質問し、既定では再利用します。

## 2. OpenFOAM用表面

```bash
surfaceTransformPoints -read-scale 0.001 \
  cad/generated/specimen_mm.stl constant/triSurface/specimen.stl
surfaceTransformPoints -read-scale 0.001 \
  cad/generated/heaterSelection_mm.stl constant/triSurface/heaterSelection.stl
```

FreeCADのmmをOpenFOAMのmへ変換します。`constant/triSurface`に置くのはこの2表面だけです。
`heaterSelection.stl`は材料形状ではなく、面選択専用です。

## 3. snappyHexMeshとリージョン分割

1. `blockMesh`: 周囲空気を含む背景直方体を生成
2. `surfaceFeatureExtract`: 試験体の特徴線を`specimen.eMesh`へ出力
3. `snappyHexMesh -overwrite`: 試験体近傍とヒーター範囲を細分化し、試験体を`solid` cellZone化
4. `topoSet`: 全セルからsolidを差し引き、残りを`fluid` cellZone化
5. `splitMeshRegions -cellZonesOnly -overwrite`: `solid`と`fluid`へ分割

特徴線は`level 3`、試験体表面は`level (2 3)`、ヒーター選択領域は`level 3`です。
ヒーター範囲を細かくする理由は、加熱面端部の温度勾配を表現し、面選択を安定させるためです。

## 4. 加熱パッチ

リージョン分割直後、鋼材と空気の界面は全て`solid_to_fluid`です。次の操作で
ヒーター部分だけを別パッチにします。

```bash
topoSet -region solid -dict system/solid/heaterFaceSetDict
createPatch -region solid -dict system/solid/createHeaterPatchDict -overwrite
topoSet -region fluid -dict system/fluid/heaterFaceSetDict
createPatch -region fluid -dict system/fluid/createHeaterPatchDict -overwrite
```

- solid側: `heaterPower`
- fluid側: `heaterCover`
- 残りの界面: `solid_to_fluid` / `fluid_to_solid`

`heaterPower`には次の合計熱量を与えます。

```text
externalWallHeatFluxTemperature
mode power
Q = 15 W (0 <= t < 300 s)
Q = 0 W  (300 <= t <= 600 s)
```

`mode power`の`Q`は熱流束ではなくパッチ全体の合計Wです。公称面積0.005 m2で
平均化した目安は$15/0.005=3000\ \mathrm{W/m^2}$ですが、離散面積に対する分配は
境界条件が行います。`heaterCover`はヒーターパッドで空気から覆われる面として断熱します。

## 5. CHT計算

`chtMultiRegionFoam`は空気の質量・運動量・エネルギーと、鋼材の熱伝導を同じ時間刻みで解きます。
流体・固体界面では温度と熱流束を結合します。基準ケースは放射なしで、自然対流だけを含みます。

保存は5秒間隔です。`postProcessing/solidTemperatureProbes/solid/0/T`には、
ヒーター側、反対側、周方向90度、上面近傍の4点がK単位で保存されます。

## 6. 確認コマンド

```bash
checkMesh -case .
paraFoam -case .
python3 python/export_temperature_history.py
```

ParaViewでヒーター面だけ確認するときは`VTK_heaterPower/solid/*.vtm`を開きます。
