# 00. このケースの作り方(参考にしたチュートリアルと手順)

更新日: 2026-07-28

## 0.1 使用したOpenFOAMディストリビューション

このケースは **ESI/OpenCFD版OpenFOAM (openfoam.com系)** のチュートリアルを参考にしている。
**OpenFOAM Foundation版(openfoam.org系, `openfoam13`など)ではない。**

判別ポイント:
- `blockMeshDict`等のヘッダコメントに `Website: www.openfoam.com` と入っている → ESI/OpenCFD版
- Foundation版は `Website: www.openfoam.org` と入る
- `WM_PROJECT_VERSION` が `v2512` のように `v`+年+月の形式 → ESI/OpenCFD版
- Foundation版は `13`, `12` のような単純整数バージョン

このマシンには複数のESI/OpenCFD版(`openfoam2406`, `openfoam2506`, `openfoam2512`)と、Foundation版(`/opt/openfoam13`)が入っている。
**最終的な動作確認は `source /usr/lib/openfoam/openfoam2512/etc/bashrc` (ESI v2512)を通した環境で行った**(当初v2406で構築・検証したのち、ユーザー指示によりv2512へ切り替えて再検証している。v2406とv2512の間で本ケースが使う機能(`#include`/`#eval`によるblockMeshDict、topoSet、splitMeshRegions、fvDOM、fvOptions、patchProbes)に構文差異は見つからなかった)。

## 0.2 参考にした公式チュートリアル(すべてESI/OpenCFD版 `$FOAM_TUTORIALS` 配下)

| 参照元 | 使った部分 |
|---|---|
| `heatTransfer/chtMultiRegionFoam/multiRegionHeater` | 単一blockMeshDict + topoSetDict(boxToCell/setToCellZone) + splitMeshRegions -cellZones によるfluid/solid分割の手順全体、Allrun.pre/Allrun/Allclean雛形、changeDictionaryDictの書き方、T/U/p/p_rghの結合パッチ境界条件の型 |
| `heatTransfer/buoyantSimpleFoam/hotRadiationRoomFvDOM` | fvDOM輻射モデルの基本設定(`radiationProperties`のfvDOMCoeffs, constantAbsorptionEmission)、SIMPLE辞書のpRefCell/pRefValue(閉じた箱でのdomain設定) |
| `heatTransfer/chtMultiRegionFoam/solarBeamWithTrees` | **chtMultiRegionFoam + fvDOM + 固体リージョンを実際に組み合わせた唯一の公式チュートリアル**。G/IDefault/qrフィールドの境界条件型(`calculated`, `wideBandDiffusiveRadiation`)、`boundaryRadiationProperties`(`opaqueDiffusive`+`wallAbsorptionEmissionModel`)の書式、T場の結合境界条件(`compressible::turbulentTemperatureRadCoupledMixed`)に`qr`/`qrNbr`を追加する書き方 |

## 0.3 このケース独自に変更・単純化した点

- ジオメトリを単純な「立方体流体ドメイン+床に立つ直方体固体」に単純化し、`system/include/caseSettings` に寸法パラメータを集約(`#include`で`blockMeshDict`と`topoSetDict`の両方から共有)。数値を変えるだけでサイズ変更できるようにした。
- 輻射の吸収放出モデルは`solarBeamWithTrees`の`multiBandAbsorption`(2バンド、太陽光スペクトル用)ではなく、`hotRadiationRoomFvDOM`と同じ単バンドの`constantAbsorptionEmission`/`constantAbsorption`を採用(本ケースでは太陽光スペクトル分割が不要なため)。
- 乱流は`multiRegionHeater`のtopAirと同じく`laminar`とし、k-epsilon等は使っていない(自然対流の可視化・輻射結合の確認が主目的のため簡略化。要件次第でRAS乱流モデルへの変更は可能)。
- 固体の底面(`floor`パッチ)を固定温度350Kの熱源接触面として設定し、ドメイン外周壁は300K固定とすることで、自然対流+輻射が駆動されるシナリオにした(公式チュートリアルにはない、本ケース固有の境界条件設定)。

## 0.4 実際に行った作業手順(このケースをゼロから作った際の手順)

1. `$FOAM_TUTORIALS/heatTransfer/chtMultiRegionFoam/multiRegionHeater` の `system/blockMeshDict`, `system/topoSetDict`, `Allrun.pre`, `Allrun`, `Allclean`, `constant/regionProperties`, `system/<region>/changeDictionaryDict` を読み、cellZoneベースの領域分割の仕組みを把握した。
2. ジオメトリを本ケース用(立方体+直方体)に単純化し、寸法を`#eval`で計算するパラメータファイル(`system/include/caseSettings`)を新規作成した。
3. `blockMesh`を実行し、意図通りのセル数・境界パッチができることを確認した(8000セル、境界6面)。
4. `topoSet`を実行し、`boxToCell`で選択されるセル数が期待値(固体: 4×4×8=128セル)と一致することを確認した。
5. `$FOAM_TUTORIALS/heatTransfer/buoyantSimpleFoam/hotRadiationRoomFvDOM` の `constant/radiationProperties`, `0.orig/G`, `constant/boundaryRadiationProperties` を読み、fvDOMの基本設定値を把握した。
6. `$FOAM_TUTORIALS/heatTransfer/chtMultiRegionFoam/solarBeamWithTrees` の `0.orig/air/{G,T,qr,IDefault}`, `constant/air/boundaryRadiationProperties`, `constant/solid/radiationProperties` を読み、**マルチリージョン+fvDOMという公式チュートリアルの実例**から、結合パッチでの輻射関連フィールドの境界条件の型を確認した。
7. 上記を踏まえて `constant/fluid/*`, `constant/solid/*`, `0.orig/*`, `system/fluid/*`, `system/solid/*` を作成した。
8. `restore0Dir` → `splitMeshRegions -cellZones -overwrite` → 各リージョンの `changeDictionary` を実行し、`fluid_to_solid`/`solid_to_fluid` という結合パッチが自動生成され、面数(144)が幾何学的な期待値と一致することを確認した。
9. `chtMultiRegionFoam` を実際に実行し、以下を確認した:
   - ログに `Selecting radiationModel fvDOM`, `fvDOM : Allocated 60 rays`, `Selecting boundary radiation Model: opaqueDiffusive`(7パッチ分)が出ること
   - NaN・FATALエラーが出ないこと
   - 固体温度が床面(350K固定)から内部へ伝導している妥当な分布になっていること
   - 輻射場`G`の値が、閉じたグレー拡散面の理論値 `4×σ×T^4`(T≈300Kで約1836 W/m^2)と近い値(1826〜1834)になっていること(物理的妥当性の定量チェック)
10. 検証後、`Allclean` でケースを初期状態に戻し、本番用の`controlDict`(endTime=60等)に戻した。

## 0.4b 追加: 発熱源(ヒートマット)とpatchProbes(2026-07-28)

ユーザー要望により、以下を追加した。

1. **固体+X側面のヒートマット発熱源(常時一定20W、寸法変更可)**
   - 参考: `$FOAM_TUTORIALS/heatTransfer/chtMultiRegionSimpleFoam/jouleHeatingSolid/system/solid/fvOptions`(固体リージョンへのfvOptions適用例)
   - `fvOptions`本体の書式(`volumeMode`, `sources`のフラット構成)は、`$WM_PROJECT_DIR/src/fvOptions/sources/general/semiImplicitSource/SemiImplicitSource.H`のUsageコメント、および実例`$FOAM_TUTORIALS/multiphase/twoPhaseEulerFoam/laminar/injection/constant/fvOptions`(`energySource1`, `scalarSemiImplicitSource`, `volumeMode absolute`)を確認して採用した。
   - **つまずいた点**: 発熱源のcellZone(`heaterMatZone`)を、fluid/solid分割用の`system/topoSetDict`に一緒に定義したところ、`splitMeshRegions -cellZones`が「メッシュ上の全cellZoneを領域分割に使おうとする」ため、solidの部分集合であるheaterMatZoneと衝突して`FOAM FATAL ERROR: ... is multiple zones`になった。
     → 対策: heaterMatZoneの作成は`splitMeshRegions`実行**後**に、`system/solid/topoSetDict`を使って`topoSet -region solid`で行うよう`Allrun.pre`を変更した。
   - **つまずいた点2**: `system/solid/fvOptions`の先頭で`#include "../include/caseSettings"`すると、fvOptionsの読み込み時に`Entry 'heaterMatWidth' is not a dictionary`という警告が出た(fvOptions直下の全エントリを「オプション定義=辞書」として走査するため、寸法パラメータのような非辞書エントリがあると警告になる)。
     → 対策: `#include`を`heaterMat{ ... }`サブ辞書の中に移動し、フィールド発熱源の定義そのものだけをfvOptions直下に置くよう修正した。
   - 実行確認: `Selecting finite volume options type scalarSemiImplicitSource` / `Source: heaterMat` / `State: active` / `selecting cells using cellZones (heaterMatZone)` がログに出ることを確認、警告なし。

2. **patchProbesによる固体表面3点の温度計測**
   - 参考: `$FOAM_TUTORIALS/heatTransfer/chtMultiRegionSimpleFoam/cpuCabinet/system/probes`(マルチリージョンでの`region`キー付き関数オブジェクトの書き方)
   - `patchProbes`自体のキー(`probeLocations`, `patches`)は`$WM_PROJECT_DIR/src/sampling/probes/patchProbes.H`のUsageコメントで確認した(`points`ではなく`probeLocations`が正しいキー名)。
   - 計測点は「ヒートマット直近(+X面)」「反対側面(-X面)」「上面中央」の3点とし、`system/solid/probes`を`system/controlDict`の`functions`から`#include`する構成にした。
   - 実行確認: `postProcessing/solidSurfaceProbes/solid/0/T`が生成され、3点とも`solid_to_fluid`パッチに距離0でスナップされ、ヒートマット直近の点だけ時間とともに温度が上昇していく(他の2点は2秒間ではまだ変化がほぼ無い、鋼材の熱拡散率から見て妥当)ことを確認した。

## 0.4c 再設計: ヒートマットを独立リージョン化(2026-07-28)

ユーザーから「セルではなく固体表面に熱量を与えたい」との指摘を受け、0.4b節の
体積発熱源(fvOptions + solidのサブcellZone)方式から、**ヒートマットを
独立リージョンとし、その外側面に熱量境界条件を与える方式**へ設計変更した。

**検討した代替案とその判断:**

1. **createPatchで既存のsolid_to_fluidパッチを分割**: 変更範囲は小さいが、
   fluid側の対応パッチも整合を取って分割する必要があり、`createPatch`の
   マッピング設定が煩雑になりやすい。
2. **blockMeshDictを完全マルチブロック化**: ヒーター足元(Y-Z方向)まで
   ブロック境界を切る必要があり、背景メッシュ全体を作り直す規模の変更に
   なる。リスク・作業量が大きい。
3. **ヒートマットを第3のリージョンにする(採用)**: `system/topoSetDict`に
   cellZoneをもう1つ追加するだけで、`splitMeshRegions`が自動的に
   `heaterMat_to_fluid`/`heaterMat_to_solid`という独立パッチを作ってくれる
   (`multiRegionHeater`チュートリアルで`leftSolid`が`bottomWater`と
   `topAir`の両方に対しそれぞれ別パッチを持っていたのと同じ仕組み)。
   既存の実証済みパイプライン(単一ブロック+cellZone+splitMeshRegions)を
   ほぼそのまま使えるため、この方式を採用した。

**実装のポイント:**

- `system/topoSetDict`: solid/fluid/heaterMatの3つのcellZoneを、互いに
  重ならないように作る(heaterMatをsolidの範囲から`action subtract`で
  差し引く)。
- `constant/regionProperties`: `solid (solid heaterMat)`として、
  heaterMatをsolidファミリーに属させる。
- `constant/heaterMat/`, `system/heaterMat/`: solidと同じ物性・スキームを
  複製。
- `system/heaterMat/changeDictionaryDict`: 外側面`heaterMat_to_fluid`に
  `externalWallHeatFluxTemperature`(mode=power)。内側面
  `heaterMat_to_solid`は通常の熱伝導結合(`compressible::
  turbulentTemperatureRadCoupledMixed`, kappaMethod solidThermo)。
- `system/fluid/changeDictionaryDict`: 新しくできる`fluid_to_heaterMat`は
  結合せず断熱壁として扱う(ヒートマットのケーシングで室内と断熱されている
  という想定)。正規表現`fluid_to_.*`より後に完全一致`fluid_to_heaterMat`を
  書くことで上書きしている(文字列完全一致は正規表現より優先される)。

**実行確認**: `Allrun.pre`後、`solid`(120セル)・`heaterMat`(8セル)・
`fluid`(7872セル)に分割され、境界パッチ面数(`solid_to_fluid`136,
`solid_to_heaterMat`20, `heaterMat_to_solid`20, `heaterMat_to_fluid`8,
`fluid_to_heaterMat`8)がすべて幾何学的な期待値と一致することを確認した。

**つまずいた点(重要、再発防止用に記録)**: `system/heaterMat/
changeDictionaryDict`で`#include "../include/caseSettings"`し、
Q(熱量)の時間テーブルを`$heaterWattage`等の変数で書いたところ、
`changeDictionary`が書き出す`0/heaterMat/T`に**変数名の文字列がそのまま
(未展開で)書き込まれ**、ソルバー起動時に
`primitiveEntry::expandVariable`でFATAL ERRORになった。
`system/solid/fvOptions`(ソルバーが直接読む)では同じ`#include`+`$var`が
問題なく機能していたのとは対照的で、**changeDictionaryDictの内容が
別ファイル(0/<region>/<field>)へ書き出される場合、埋め込まれた
Function1のtable等の中の`$var`は展開されずそのまま文字列として
コピーされる**ことがわかった。対策として、Qテーブルはリテラル値で
直接書くことにした(`system/heaterMat/changeDictionaryDict`内のコメント、
および101_0/README.md参照)。

## 0.5 バージョン間の注意

- `#include`・`#eval #{ ... #}` によるパラメトリックなblockMeshDictは、ESI/OpenCFD版の比較的新しいバージョン(v1906以降)でサポートされている機能。Foundation版(openfoam.org)では書式が異なる場合があるため、Foundation版で使う場合は書き換えが必要になる可能性がある(**要確認**)。
- `compressible::turbulentTemperatureRadCoupledMixed` 等の境界条件名やライブラリ構成は、ESI/OpenCFD版とFoundation版で細部(利用可能なキーワード、必須引数)が異なることがあるため、別バージョンに移植する場合は該当バージョンのソース(`src/thermophysicalModels/radiation`, `src/TurbulenceModels`)またはチュートリアルで再確認すること。
- 2026-07-28時点でこのマシンには `openfoam2406`(現在使用), `openfoam2506`, `openfoam2512` がインストール済み、`openfoam2606` はaptパッケージとしては存在するが未インストール。バージョンを variar する場合は、`source /usr/lib/openfoam/openfoamXXXX/etc/bashrc` のように切り替えてから `Allclean && Allrun` で再検証すること。
