"""データ同化の結果を ParaView で開ける VTU ファイルに書き出す.

出力(paraview/ フォルダ):
  temperature_fields.vtu … 固体実メッシュ(20696セル)に3つの温度をセルデータで格納
        T_garbage_degC / T_da_degC / T_truth_degC
  displacement_fields.vtu … FrontISTR 円筒メッシュに変位ベクトルを点データで格納
        U_garbage_m / U_da_m / U_truth_m (と Uz_*_um)

ParaView での見方:
  - temperature_fields.vtu を開き、Coloring を T_da_degC に。Clip で断面表示。
  - displacement_fields.vtu を開き、Filters > Warp By Vector に U_da_m を指定(倍率大)。
    Coloring を Uz_da_um にすると傾きが見える。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/export_paraview.py
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
F102 = os.path.abspath(os.path.join(
    ROOT, "..", "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
sys.path.insert(0, F102)

K = 273.15
FIELDS = os.path.join(ROOT, "openfoam", "run_fem_enkf", "fields")
CASE = os.path.join(ROOT, "openfoam", "run_fem_enkf", "truth", "case.foam")
OUT = os.path.join(ROOT, "paraview")

CASES = [("garbage", "ensmean_t0.npy"),
         ("da", "ensmean_t20.npy"),
         ("truth", "truth_final.npy")]


def export_temperature():
    import pyvista as pv
    from scipy.spatial import cKDTree
    if not os.path.exists(CASE):
        open(CASE, "w").close()
    r = pv.POpenFOAMReader(CASE)
    r.set_active_time_value(r.time_values[-1])
    solid = r.read()["solid"]
    grid = solid["internalMesh"] if "internalMesh" in solid.keys() else solid.combine()

    centres = np.load(os.path.join(FIELDS, "cell_centres.npy"))
    idx = cKDTree(centres).query(np.asarray(grid.cell_centers().points))[1]
    for tag, fn in CASES:
        grid.cell_data[f"T_{tag}_degC"] = np.load(os.path.join(FIELDS, fn))[idx] - K
    path = os.path.join(OUT, "temperature_fields.vtu")
    grid.save(path)
    print(f"[paraview] wrote {os.path.relpath(path, ROOT)} ({grid.n_cells} cells)")


def export_displacement():
    import pyvista as pv
    import vtk
    from run.plot_da_displacement_field import fem_displacement_field
    import cylinder_mesh  # 102_1

    centres = np.load(os.path.join(FIELDS, "cell_centres.npy"))
    # 円筒メッシュ(節点・六面体要素)を1回作る
    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    node_coords = np.array([xyz for _nid, xyz in mesh["nodes"]], dtype=float)
    id_to_row = {nid: i for i, (nid, _xyz) in enumerate(mesh["nodes"])}
    cells = []
    for _eid, conn in mesh["elements"]:
        cells.append(8)
        cells.extend(id_to_row[n] for n in conn)
    cells = np.array(cells)
    celltypes = np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, dtype=np.uint8)
    ug = pv.UnstructuredGrid(cells, celltypes, node_coords)

    for tag, fn in CASES:
        cell_T = np.load(os.path.join(FIELDS, fn))
        work = os.path.join(ROOT, "openfoam", "fem_disp_" + fn.replace(".npy", ""))
        coords, U, _ = fem_displacement_field(centres, cell_T, work)
        ug.point_data[f"U_{tag}_m"] = U
        ug.point_data[f"Uz_{tag}_um"] = U[:, 2] * 1e6
        print(f"[paraview]  {tag}: Uz {U[:,2].min()*1e6:.2f}..{U[:,2].max()*1e6:.2f} um")
    path = os.path.join(OUT, "displacement_fields.vtu")
    ug.save(path)
    print(f"[paraview] wrote {os.path.relpath(path, ROOT)} ({ug.n_points} nodes)")


def _cylinder_hex_grid():
    """FrontISTR 円筒メッシュ(六面体)を pyvista UnstructuredGrid で返す."""
    import pyvista as pv
    import vtk
    import cylinder_mesh  # 102_1
    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    node_coords = np.array([xyz for _nid, xyz in mesh["nodes"]], dtype=float)
    id_to_row = {nid: i for i, (nid, _xyz) in enumerate(mesh["nodes"])}
    cells = []
    for _eid, conn in mesh["elements"]:
        cells.append(8)
        cells.extend(id_to_row[n] for n in conn)
    celltypes = np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, dtype=np.uint8)
    return pv.UnstructuredGrid(np.array(cells), celltypes, node_coords), node_coords


def export_timehistory():
    """ROM の 0→600s 温度時刻歴を円筒メッシュへ補間し ParaView 時系列(.pvd)で出力."""
    import yaml
    from dacore.calibrate import load_calibrated
    from dacore.displacement import load_operator
    from dacore.node_locations import NODE_XYZ
    from dacore.twin_fem import run_rom_fem_twin

    cfg = yaml.safe_load(open(os.path.join(ROOT, "config", "da_config.yaml")))
    hist = run_rom_fem_twin(cfg, load_calibrated(), load_operator())

    grid, node_coords = _cylinder_hex_grid()
    rom_nodes = np.array(list(NODE_XYZ.values()))
    # IDW 重み(FEM節点 × 5ROMノード)を1回だけ作る
    W = np.zeros((len(node_coords), len(rom_nodes)))
    for i, p in enumerate(node_coords):
        d = np.linalg.norm(rom_nodes - p, axis=1)
        if d.min() < 1e-9:
            W[i, d.argmin()] = 1.0
        else:
            w = 1.0 / d**2
            W[i] = w / w.sum()

    ts_dir = os.path.join(OUT, "timehistory")
    os.makedirs(ts_dir, exist_ok=True)
    datasets = []
    for k, t in enumerate(hist["times"]):
        g = grid.copy()
        g.point_data["T_da_degC"] = W @ (hist["da_T"][k] - K)
        g.point_data["T_free_degC"] = W @ (hist["free_T"][k] - K)
        g.point_data["T_truth_degC"] = W @ (hist["truth_T"][k] - K)
        fn = os.path.join("timehistory", f"da_t{t:g}.vtu")
        g.save(os.path.join(OUT, fn))
        datasets.append(f'    <DataSet timestep="{t:g}" file="{fn}"/>')
    pvd = os.path.join(OUT, "temperature_timehistory.pvd")
    with open(pvd, "w") as f:
        f.write('<?xml version="1.0"?>\n'
                '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n'
                '  <Collection>\n' + "\n".join(datasets) +
                '\n  </Collection>\n</VTKFile>\n')
    print(f"[paraview] wrote {os.path.relpath(pvd, ROOT)} "
          f"({len(hist['times'])} 時刻ステップ, 0→600s)")


def main():
    os.makedirs(OUT, exist_ok=True)
    export_temperature()
    export_displacement()
    export_timehistory()
    print(f"[paraview] ParaView で {os.path.relpath(OUT, ROOT)}/*.vtu / *.pvd を開いてください")


if __name__ == "__main__":
    main()
