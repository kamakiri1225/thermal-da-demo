# 101_0: OpenFOAM chtMultiRegionFoam — CHT + 輻射 + ヒートマット発熱源

更新日: 2026-07-28
OpenFOAM: ESI/OpenCFD版 v2512で動作確認済み(v2406でも構築・検証。詳細は`docs/00_tutorial_derivation.md`)。

## 概要

立方体の流体ドメイン(既定1000x1000x1000mm)の中央、床面に固体の直方体(既定200x200x400mm)が立っている構成で、以下を計算する。

- 固体・流体間の共役熱伝達(CHT, chtMultiRegionFoam)
- 空気の自然対流(浮力駆動、`constant/g`)
- 輻射(fvDOM、固体表面⇔周囲壁のグレー拡散面間の熱交換)
- 固体の片側面(+X面)に貼り付けたヒートマット発熱源(既定20W、時間変化あり: 0-300秒ON→300-600秒OFF、寸法・位置・時間は変更可能)
- 固体・ヒートマット表面の温度をpatchProbesで計測

初期温度は20℃(293.15K)。得られる固体の温度分布は、`../101_1_frontistr_cht_box_thermal_expansion`でFrontISTRへ引き継ぎ、熱膨張(底面固定)を計算する材料として使う(時刻歴でも計算可能)。

## ジオメトリ(自由に変更可能)

すべて`system/include/caseSettings`に集約している。この1ファイルの数値を変えるだけで、`blockMeshDict`(ジオメトリ・メッシュ解像度)と`topoSetDict`(固体/流体/ヒートマット領域の切り分け)の両方に反映される。

| 変数 | 既定値 | 意味 |
|---|---|---|
| `domainLx/Ly/Lz` | 1.0 (=1000mm) | 流体ドメイン(立方体)の寸法 |
| `solidLx/Ly/Lz` | 0.2/0.2/0.4 (=200/200/400mm) | 固体直方体の寸法(床z=0に設置、X,Y中央) |
| `cellSize` | 0.05 (50mm) | メッシュ解像度。domain/solidの寸法の整数分の1にすること |
| `heaterWattage` | 20.0 [W] | ヒートマットON時の出力 |
| `heaterOnTime` / `heaterOffTime` | 300 / 600 [s] | ヒーターON継続時間 / OFFへ切り替わる時刻(既定: 0-300s=ON, 300-600s=OFF) |
| `heaterMatWidth/Height` | 0.10/0.20 [m] | ヒートマットの寸法(cellSizeの偶数倍にすること、`caseSettings`のコメント参照) |

**メッシュを細かくしたい場合は`cellSize`を小さくするだけでよい**(例: 0.025m → 40x40x40の背景メッシュ)。

**時間変化するヒーターの数値を変える場合の注意**: `caseSettings`の`heaterWattage`/`heaterOnTime`/`heaterOffTime`を変えるだけでは反映されない。`system/heaterMat/changeDictionaryDict`のQテーブルはリテラル値で書く必要がある(`#include`によるドット変数展開が`changeDictionary`の書き出し時に効かないため。詳細は同ファイル内のコメントと`docs/00_tutorial_derivation.md`参照)。**両方を手動で揃えること。**

## 全体構成

```mermaid
flowchart LR
    subgraph Specimen["試験片系"]
        ROD["固体ブロック(solid)"]
        HEATMAT["ヒートマット(heaterMat)<br/>独立リージョン、+X面"]
    end
    ROD -- "solid_to_heaterMat<br/>(熱伝導結合)" --> HEATMAT
    ROD -- "solid_to_fluid<br/>(CHT結合パッチ)" --> AIR["流体(空気)"]
    HEATMAT -- "fluid_to_heaterMat<br/>(断熱、非結合)" -.- AIR
    AIR -- "自然対流+輻射" --> WALLS["ドメイン外周壁(20degC固定)"]
    HEATER_BC["externalWallHeatFluxTemperature<br/>mode=power, Q=table(時間変化)"] -->|"heaterMat_to_fluid面"| HEATMAT
```

ヒートマットは**独立リージョン**(solidのサブゾーンではない)。外側面(`heaterMat_to_fluid`)は熱量境界条件(`externalWallHeatFluxTemperature`, mode=power)で、ケーシングにより室内空気とは直接熱交換しない(断熱)想定。投入した熱はすべて内側面(`heaterMat_to_solid`)を通じて固体へ伝導する。

## ドキュメント

- **`docs/01_setup_openfoam_and_frontistr.md`: GitHubからcloneした状態からOpenFOAM/FrontISTRをインストールし、101_0→101_1を実行するまでの手順(コピペで進められる詳細なコマンド付き)**
- `docs/00_tutorial_derivation.md`: どの公式チュートリアルを参考にどう組み立てたか、つまずいた点と対処(splitMeshRegionsのcellZone競合、fvOptionsの警告等)、バージョン間の注意

FrontISTRへ温度時刻歴を渡し、熱膨張を計算する手順は
[`../101_1_frontistr_cht_box_thermal_expansion/README.md`](../101_1_frontistr_cht_box_thermal_expansion/README.md)と、
その[`docs/00_openfoam_frontistr_coupling_workflow.md`](../101_1_frontistr_cht_box_thermal_expansion/docs/00_openfoam_frontistr_coupling_workflow.md)を参照する。

GitHubには再計算に必要な`0.orig/`、`constant/`、`system/`、実行スクリプト、解説を
収録する。0～600秒の時刻フォルダは約1GBあるため収録せず、`./Allrun`で再生成する。

## 使い方

**OpenFOAM未インストールの場合は先に`docs/01_setup_openfoam_and_frontistr.md`を参照。**

```bash
source /usr/lib/openfoam/openfoam2512/etc/bashrc   # 環境に応じて変更
./Allrun            # blockMesh → topoSet → splitMeshRegions → changeDictionary → chtMultiRegionFoam
# 掃除する場合
./Allclean
```

既定の`controlDict`は `endTime 60`(秒)。8000セル程度の背景メッシュであれば数分で完了する。より長時間の加熱で変形を大きくしたい場合(101_1のFrontISTR側でより明瞭な熱膨張を見たい場合)は`endTime`を数百〜数千秒に延ばすとよい。

## 主な設計判断・前提

- **単一ブロックのblockMeshDict + topoSetDict(boxToCell/setToCellZone)+ splitMeshRegions -cellZones** という、公式チュートリアル`multiRegionHeater`と同じ方式で固体/流体/ヒートマットの3領域を分割している(3つのcellZoneは互いに排他的、`system/topoSetDict`参照)。
- 乱流は`laminar`(層流)としている。自然対流のRayleigh数次第では実際は乱流域の可能性があるが、CHT+輻射+発熱源の結合が正しく動くことの実証を優先し、乱流モデルは簡略化した。
- ヒートマットは**独立リージョン**(solidから切り出した厚さ1セルの薄い直方体)とし、その外側面(`heaterMat_to_fluid`)へ`externalWallHeatFluxTemperature`(mode=power)で熱量[W]を直接与える(セルではなく面への熱量入力、`system/heaterMat/changeDictionaryDict`)。ケーシングで室内側とは断熱されている想定のため、fluid側の対応パッチ(`fluid_to_heaterMat`)は単純な断熱壁としている。
- 固体の底面(floor)は熱的には断熱(zeroGradient)とし、構造解析(101_1)側で「固定端」として使う面という位置づけにしている。
- 初期温度は20℃(293.15K)。ドメイン外周壁もすべて20℃固定(輻射のグレー拡散面としても不透明・輻射率0.9)。
- ヒーターは時間変化する熱量として設定(既定: 0-300秒=20W、300-600秒=0W)。`system/heaterMat/changeDictionaryDict`のQフィールドをFunction1の`table`で与えている。

## 実行確認(2026-07-28、heaterMatリージョン化後)

- `Allrun.pre`実行後、`solid`(120セル)・`heaterMat`(8セル)・`fluid`(7872セル)の3リージョンに分割され、`solid_to_heaterMat`(20面)・`heaterMat_to_fluid`(8面)・`heaterMat_to_solid`(20面)・`fluid_to_heaterMat`(8面)が自動生成されることを確認(幾何学的な期待値と一致)。
- `endTime=600`(1サイクル分、ヒーターON 0-300s→OFF 300-600s)で実行し、エラーなく完走することを確認。
- **要注意(実機で発見した制約)**: `changeDictionaryDict`内で`#include`した`caseSettings`の変数(`$heaterWattage`等)をFunction1の`table`の中で使うと、`changeDictionary`が書き出す`0/heaterMat/T`に**変数名の文字列がそのまま(未展開で)書き込まれ**、ソルバー起動時にFATAL ERRORになった。このため`system/heaterMat/changeDictionaryDict`のQテーブルはリテラル値で書いている(上記「ジオメトリ」節の注意も参照)。

## 未確定・要確認事項

- ヒートマットの実際の寸法・貼付位置・時間プロファイル(既定値は仮案)
- 固体の材料物性(鋼材相当の代表値、実材料が決まったら`constant/solid`・`constant/heaterMat`の`thermophysicalProperties`と`101_1/config/material_properties_steel.json`の両方を更新すること)
- 輻射率(壁0.9、固体表面0.85)は代表値
- 乱流モデル(現状laminar、必要ならRAS kEpsilon等へ変更)
- ヒートマット外側面を完全断熱(fluidと非結合)とした近似の妥当性(実際のケーシングの断熱性能次第)
