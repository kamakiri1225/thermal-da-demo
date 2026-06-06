# VTK (legacy) → VTU + PVD 変換メモ

## 目的

FrontISTR が出力するレガシー VTK ファイル (`.vtk`) は時刻情報を持たない。
そのため ParaView で時系列表示すると連番 (1, 2, 3, ...) が時刻として扱われ、
OpenFOAM の時刻ラベル (5s, 10s, ...) と同期が取れない。

VTU + PVD 方式に変換することで ParaView の時刻スライダーを秒単位に揃えられる。

---

## ファイル形式の比較

| 形式 | 拡張子 | 時刻管理 | ParaView での時刻表示 |
|------|--------|----------|----------------------|
| レガシー VTK | `.vtk` | 不可 | ファイル連番 (1, 2, ...) |
| VTK XML (UnstructuredGrid) | `.vtu` | PVD で管理 | PVD の timestep 値 (秒) |
| ParaView Data Collection | `.pvd` | 各 .vtu と時刻を対応付ける XML | — |

---

## 変換スクリプト

Python の `vtk` ライブラリ (ParaView に同梱) を使用。

```python
import vtk
import os

VTK_DIR = 'oi/results/with_da/vtk'   # 元の .vtk が置いてあるディレクトリ
VTU_DIR = os.path.join(VTK_DIR, 'vtu')
os.makedirs(VTU_DIR, exist_ok=True)

DT = 5.0          # 1ステップあたりの物理時間 [s]
N_STEPS = 180     # 総ステップ数

reader = vtk.vtkUnstructuredGridReader()
writer = vtk.vtkXMLUnstructuredGridWriter()
pvd_entries = []

for step in range(1, N_STEPS + 1):
    t = step * DT
    src = os.path.join(VTK_DIR, f'with_da_{step:04d}.vtk')
    dst = os.path.join(VTU_DIR, f'with_da_{step:04d}.vtu')

    # 1. レガシー VTK を読む
    reader.SetFileName(src)
    reader.Update()

    # 2. VTU として書き出す
    writer.SetFileName(dst)
    writer.SetInputData(reader.GetOutput())
    writer.Write()

    pvd_entries.append(
        f'    <DataSet timestep="{t:.1f}" file="vtu/with_da_{step:04d}.vtu"/>'
    )

# 3. PVD ファイルを生成
pvd_lines = [
    '<?xml version="1.0"?>',
    '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
    '  <Collection>',
] + pvd_entries + [
    '  </Collection>',
    '</VTKFile>',
]
with open(os.path.join(VTK_DIR, 'with_da_seconds.pvd'), 'w') as f:
    f.write('\n'.join(pvd_lines) + '\n')
```

---

## PVD ファイルの構造

```xml
<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
  <Collection>
    <DataSet timestep="5.0"   file="vtu/with_da_0001.vtu"/>
    <DataSet timestep="10.0"  file="vtu/with_da_0002.vtu"/>
    <DataSet timestep="15.0"  file="vtu/with_da_0003.vtu"/>
    ...
    <DataSet timestep="900.0" file="vtu/with_da_0180.vtu"/>
  </Collection>
</VTKFile>
```

- `timestep` に物理時刻 [秒] を指定する
- `file` は PVD ファイルからの相対パス

---

## 使い方 (ParaView)

1. ParaView で `with_da_seconds.pvd` を開く
2. Pipeline Browser → Apply
3. 時刻スライダーが `5.0, 10.0, ..., 900.0` と秒単位で表示される
4. OpenFOAM の `case.foam` と並べると同じ時刻 (例: `Time: 230`) で比較できる

---

## 出力ファイル構成

```
vtk/
  with_da_0001.vtk  ... with_da_0180.vtk   # 元の FrontISTR 出力 (変更なし)
  with_da_seconds.pvd                       # PVD インデックス (新規)
  vtu/
    with_da_0001.vtu  ... with_da_0180.vtu  # 変換済み VTU (新規)
```

---

## 注意

- `vtk` Python パッケージは ParaView 付属のものか `pip install vtk` で入手可能
- `meshio` ライブラリでも同様の変換が可能 (`meshio.read` / `meshio.write`)
- truth / model_only ケースも同じスクリプトでディレクトリ名を変えて対応できる
