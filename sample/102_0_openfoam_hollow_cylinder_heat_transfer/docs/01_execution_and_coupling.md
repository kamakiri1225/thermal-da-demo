# 実行とFrontISTR連携

## 前提

- ESI/OpenCFD版OpenFOAM v2406以降
- FreeCAD 1.x（既定パス以外なら`FREECADCMD`を設定）
- Python 3、NumPy、Matplotlib、PyYAML
- 102_1を実行する場合はFrontISTRの`fistr1`

## 初回

```bash
cd sample/102_0_openfoam_hollow_cylinder_heat_transfer
export FREECADCMD=/path/to/freecadcmd
./Allrun.pre
chtMultiRegionFoam
python3 python/export_temperature_history.py
```

`Allrun.pre`が行う各コマンドは`README.md`と`docs/00_model_and_mesh_workflow.md`を参照してください。

## CAD・メッシュを再利用する計算

```bash
./Allrun
```

既存前処理を使う場合は質問にEnterまたは`N`で答えます。ヒーター出力だけを変える場合は
`system/solid/changeDictionaryDict`と`0/solid/T`のQテーブルを揃えてから計算します。

## FrontISTR

```bash
cd ../102_1_frontistr_hollow_cylinder_thermal_expansion
python3 python/run_thermal_expansion.py \
  --time 0 --case-dir case/uniform_10K --uniform-delta-t 10
python3 python/run_timehistory.py --interval 5
```

Pythonは各OpenFOAM時刻の`solid/C`と`solid/T`をFrontISTR節点へ補間し、
`!TEMPERATURE`として渡します。詳細は102_1の`README.md`を参照してください。
