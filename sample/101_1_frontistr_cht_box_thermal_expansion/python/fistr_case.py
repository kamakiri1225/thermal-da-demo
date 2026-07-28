"""
101_0(OpenFOAM CHT+輻射)の固体温度分布を取り込み、FrontISTRで
熱膨張(線形静解析、底面固定)を計算するためのケース一式を書き出し、
fistr1を実行し、変位結果を読み込むモジュール。

HEC-MWの.msh/.cnt書式は、本リポジトリの既存ケース
sample/002_2_frontistr_round_bar_da/fistr_interface.py の実装を
そのまま参考にしている(丸棒の熱伝導デモ用に実際に動作確認済みの書式)。
"""

from __future__ import annotations

import os
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from box_mesh import build_box_mesh

FISTR_BIN = Path(os.environ.get("FISTR1", "/home/kamakiri/local/frontistr/bin/fistr1"))

CASE_NAME = "box_thermal_expansion"


def _format_group_ids(ids: list[int], per_line: int = 8) -> list[str]:
    lines = []
    for i in range(0, len(ids), per_line):
        lines.append(" " + ", ".join(str(x) for x in ids[i:i + per_line]))
    return lines


def write_mesh(
    case_dir: Path,
    nx: int,
    ny: int,
    nz: int,
    lx: float,
    ly: float,
    lz: float,
    young_modulus: float,
    poisson_ratio: float,
    density: float,
    thermal_expansion_coeff: float,
) -> dict:
    """FrontISTR用の.mshを書き出す。

    経緯(要確認事項、試行錯誤の記録): !SECTIONのMATERIAL参照は、メッシュ
    読み込み時点(.cntを読む前)に解決される。そのため.mshファイル自身にも
    同名のMATERIALが必要。ただし.msh側のMATERIALパーサは汎用の
    !ITEM=n/!SUBITEM=n形式(HECMWの生の材料定義形式、
    fistr1/src/lib/physics/material.f90 の M_YOUNGS=1,M_POISSON=2,
    M_EXAPNSION=20 という内部配列位置に対応)しか受け付けず、
    !ELASTIC/!EXPANSION_COEFFのような制御ファイル固有の簡易記法では
    "!ITEM required" エラーになる。
    実例: tests/lib/static_LIB_C3D8_Bbar/elastic_beam_thermal.msh は
        !MATERIAL, NAME=MATERIAL1, ITEM=3
        !ITEM=1, SUBITEM=2
         <E>, <nu>
        !ITEM=2, SUBITEM=1
         <density></density>
        !ITEM=3, SUBITEM=1
         <expansion coeff>
    という形式で書かれていた。同じ物性値を .cnt 側にも
    !MATERIAL/!ELASTIC/!EXPANSION_COEFF として重複定義しており
    (write_cnt参照)、実際の解析にはそちらが使われる模様(値を完全に
    一致させておけば、どちらが優先されても結果は変わらない)。
    mesh dict(nodes/elements/bottom_nodes/top_nodes)を返す。"""
    mesh = build_box_mesh(nx, ny, nz, lx, ly, lz)

    lines = ["!HEADER", " FrontISTR box thermal expansion (101_1)", "!NODE"]
    for node_id, (x, y, z) in mesh["nodes"]:
        lines.append(f"{node_id:8d}, {x:.12g}, {y:.12g}, {z:.12g}")

    lines.append("!ELEMENT, TYPE=361")
    for elem_id, conn in mesh["elements"]:
        lines.append(f"{elem_id:8d}, " + ", ".join(f"{n:8d}" for n in conn))

    lines.append("!NGROUP, NGRP=NALL")
    lines.extend(_format_group_ids(mesh["all_nodes"]))
    lines.append("!NGROUP, NGRP=BOTTOM")
    lines.extend(_format_group_ids(mesh["bottom_nodes"]))
    lines.append("!NGROUP, NGRP=TOP")
    lines.extend(_format_group_ids(mesh["top_nodes"]))

    lines.append("!EGROUP, EGRP=EALL")
    lines.extend(_format_group_ids(mesh["all_elements"]))

    lines.append("!MATERIAL, NAME=STEEL, ITEM=3")
    lines.append("!ITEM=1, SUBITEM=2")
    lines.append(f" {young_modulus:.10g}, {poisson_ratio:.10g}")
    lines.append("!ITEM=2, SUBITEM=1")
    lines.append(f" {density:.10g}")
    lines.append("!ITEM=3, SUBITEM=1")
    lines.append(f" {thermal_expansion_coeff:.10g}")

    lines.append("!SECTION, TYPE=SOLID, EGRP=EALL, MATERIAL=STEEL")
    lines.append(" 1.0")

    lines.append("!END")

    (case_dir / f"{CASE_NAME}.msh").write_text("\n".join(lines) + "\n")
    return mesh


def write_hecmw_ctrl(case_dir: Path) -> None:
    (case_dir / "hecmw_ctrl.dat").write_text(
        f"""\
##
## HEC-MW control file for FrontISTR box thermal expansion (101_1)
##
!MESH, NAME=fstrMSH, TYPE=HECMW-ENTIRE
{CASE_NAME}.msh
!CONTROL, NAME=fstrCNT
{CASE_NAME}.cnt
!RESULT, NAME=fstrRES, IO=OUT
{CASE_NAME}.res
!RESULT, NAME=vis_out, IO=OUT
{CASE_NAME}_vis
"""
    )


def write_cnt(
    case_dir: Path,
    node_ids: list[int],
    node_temperatures: np.ndarray,
    reference_temperature: float,
    young_modulus: float,
    poisson_ratio: float,
    thermal_expansion_coeff: float,
) -> None:
    """STATIC解析用.cntを書く。底面(BOTTOM)をXYZ全固定し、全節点へ
    node_temperatures(OpenFOAMから写像した節点温度)を熱荷重として与える。

    熱ひずみは alpha*(TEMPERATURE - REFTEMP) で計算される(FrontISTRの
    !REFTEMP機構、tests/analysis/static/exF/F241.cnt で確認した書式)。
    REFTEMPは101_0の初期温度(293.15K, 20degC)に合わせる。"""
    lines = [
        "!VERSION",
        " 3",
        "!SOLUTION, TYPE=STATIC",
        "!WRITE,RESULT,FREQUENCY=1",
        "!WRITE,VISUAL",
        "!SOLVER,METHOD=CG,PRECOND=1,NSET=0,ITERLOG=NO,TIMELOG=NO",
        " 5000, 1",
        " 1.0e-08, 1.00, 0.0",
        "!REFTEMP",
        f" {reference_temperature:.10g}",
        "!INITIAL_CONDITION, TYPE=TEMPERATURE",
        f" NALL, {reference_temperature:.10g}",
        "!BOUNDARY, GRPID=1",
        " BOTTOM,1,1",
        " BOTTOM,2,2",
        " BOTTOM,3,3",
        "!TEMPERATURE, GRPID=1",
    ]
    for node_id, t in zip(node_ids, node_temperatures):
        lines.append(f" {node_id}, {t:.10g}")

    lines.extend(
        [
            # !SECTION は .msh 側 (write_mesh) で定義済み。
            "!MATERIAL, NAME=STEEL",
            "!ELASTIC",
            f" {young_modulus:.10g}, {poisson_ratio:.10g}",
            "!EXPANSION_COEFF",
            f" {thermal_expansion_coeff:.10g}",
            "!STEP, SUBSTEPS=1, CONVERG=1.000E-07",
            "BOUNDARY,1",
            "LOAD,1",
            "!VISUAL, method=PSR",
            "!surface_num=1",
            "!surface 1",
            "!output_type = VTK",
            "!END",
        ]
    )
    (case_dir / f"{CASE_NAME}.cnt").write_text("\n".join(lines) + "\n")


def run_fistr(case_dir: Path) -> None:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    result = subprocess.run(
        [str(FISTR_BIN)],
        cwd=case_dir,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    (case_dir / "log.fistr1").write_text(result.stdout)
    if result.returncode != 0:
        raise RuntimeError(f"fistr1 failed (exit={result.returncode}). See {case_dir / 'log.fistr1'}\n{result.stdout[-2000:]}")


def add_temperature_to_visualization(
    case_dir: Path,
    node_coords: np.ndarray,
    node_temperatures: np.ndarray,
    time_value: float,
) -> Path:
    """FrontISTRの最終PVTU/VTUへ節点温度と実時刻を追加する。

    FrontISTRのSTATIC解析では、!TEMPERATUREを熱荷重として使用しても標準の
    可視化結果に温度配列が含まれない。この関数はVTK節点座標を入力メッシュの
    節点座標へ照合し、TEMPERATURE[K]をPointDataへ追加する。戻り値は更新した
    PVTUのパス。
    """
    pvtu_candidates = sorted(case_dir.glob(f"{CASE_NAME}_vis_psf.*.pvtu"))
    if not pvtu_candidates:
        raise FileNotFoundError(f"FrontISTR visualization PVTU not found in {case_dir}")
    pvtu_path = pvtu_candidates[-1]

    pvtu_tree = ET.parse(pvtu_path)
    pvtu_root = pvtu_tree.getroot()
    ppoint_data = pvtu_root.find("./PUnstructuredGrid/PPointData")
    if ppoint_data is None:
        raise ValueError(f"PPointDataが見つかりません: {pvtu_path}")
    for old in list(ppoint_data):
        if old.attrib.get("Name") == "TEMPERATURE":
            ppoint_data.remove(old)
    ET.SubElement(
        ppoint_data,
        "PDataArray",
        {"type": "Float32", "Name": "TEMPERATURE", "NumberOfComponents": "1", "format": "ascii"},
    )
    _set_vtk_time_value(pvtu_root, time_value, parallel=True)

    pieces = pvtu_root.findall("./PUnstructuredGrid/Piece")
    if not pieces:
        raise ValueError(f"Pieceが見つかりません: {pvtu_path}")
    for piece in pieces:
        source = piece.attrib.get("Source")
        if not source:
            raise ValueError(f"Piece Sourceがありません: {pvtu_path}")
        vtu_path = (pvtu_path.parent / source).resolve()
        _add_temperature_to_vtu(vtu_path, node_coords, node_temperatures, time_value)

    ET.indent(pvtu_tree, space="  ")
    pvtu_tree.write(pvtu_path, encoding="utf-8", xml_declaration=True)
    return pvtu_path


def _add_temperature_to_vtu(
    vtu_path: Path,
    node_coords: np.ndarray,
    node_temperatures: np.ndarray,
    time_value: float,
) -> None:
    tree = ET.parse(vtu_path)
    root = tree.getroot()
    piece = root.find("./UnstructuredGrid/Piece")
    if piece is None:
        raise ValueError(f"UnstructuredGrid/Pieceが見つかりません: {vtu_path}")
    point_array = piece.find("./Points/DataArray")
    point_data = piece.find("./PointData")
    if point_array is None or point_array.text is None or point_data is None:
        raise ValueError(f"PointsまたはPointDataが見つかりません: {vtu_path}")

    vtk_coords = np.fromstring(point_array.text, sep=" ").reshape(-1, 3)
    source_coords = np.asarray(node_coords, dtype=float)
    source_temps = np.asarray(node_temperatures, dtype=float)
    mapped_temps = np.empty(len(vtk_coords), dtype=float)
    for i, xyz in enumerate(vtk_coords):
        distances = np.linalg.norm(source_coords - xyz, axis=1)
        nearest = int(np.argmin(distances))
        if distances[nearest] > 1.0e-6:
            raise ValueError(
                f"VTK節点を入力メッシュへ対応付けられません: {vtu_path}, "
                f"point={xyz.tolist()}, distance={distances[nearest]}"
            )
        mapped_temps[i] = source_temps[nearest]

    for old in list(point_data):
        if old.attrib.get("Name") == "TEMPERATURE":
            point_data.remove(old)
    temperature_array = ET.SubElement(
        point_data,
        "DataArray",
        {"type": "Float32", "Name": "TEMPERATURE", "NumberOfComponents": "1", "format": "ascii"},
    )
    temperature_array.text = "\n" + "\n".join(f"{t:.9g}" for t in mapped_temps) + "\n"
    _set_vtk_time_value(root, time_value, parallel=False)

    ET.indent(tree, space="  ")
    tree.write(vtu_path, encoding="utf-8", xml_declaration=True)


def _set_vtk_time_value(root: ET.Element, time_value: float, parallel: bool) -> None:
    grid_name = "PUnstructuredGrid" if parallel else "UnstructuredGrid"
    data_array = root.find(f"./{grid_name}/FieldData/DataArray[@Name='TimeValue']")
    if data_array is not None:
        data_array.text = f"{time_value:.12g}"


def _latest_result_file(case_dir: Path) -> Path:
    # "res.0.0"(荷重ステップ前)より "res.0.1"(荷重適用後)を使いたいので、
    # ステップ番号(末尾の整数)でソートする(文字列ソートだと桁数次第で誤る)。
    candidates = list((case_dir).glob(f"{CASE_NAME}.res.0.*"))
    if not candidates:
        raise FileNotFoundError(f"FrontISTR result not found in {case_dir}")
    candidates.sort(key=lambda p: int(p.suffix.lstrip(".")))
    return candidates[-1]


def read_displacement(case_dir: Path) -> dict[int, tuple[float, float, float]]:
    """FrontISTR結果ファイル(*fstrresult 2.0形式)からDISPLACEMENTを読む。

    フォーマット(実際の出力を確認して把握した):
        *data
        <n_nodes> <n_elements>
        <n_node_value_groups> <n_elem_value_groups>
        <各node value groupの成分数> ...          (例: "3 6 1" = DISP(3),STRESS(6),MISES(1))
        <node value groupの名前(1行ずつ)>          (例: DISPLACEMENT / NodalSTRESS / NodalMISES)
        [elem value groupがあれば同様に成分数・名前]
        <node_id>
        <成分数の合計だけの数値(改行で折り返される場合がある)>
        ...(n_nodes回繰り返し)
    DISPLACEMENTは常に最初のnode value groupという前提(本ケースの.cnt構成では
    常にそうなる)。"""
    path = _latest_result_file(case_dir)
    lines = path.read_text().splitlines()

    try:
        data_idx = next(i for i, line in enumerate(lines) if line.strip() == "*data")
    except StopIteration as exc:
        raise ValueError(f"'*data'セクションが見つかりません: {path}") from exc

    n_nodes, n_elements = (int(x) for x in lines[data_idx + 1].split())
    n_node_groups, n_elem_groups = (int(x) for x in lines[data_idx + 2].split())

    widths_line = lines[data_idx + 3].split()
    node_widths = [int(x) for x in widths_line[:n_node_groups]]
    elem_widths = [int(x) for x in widths_line[n_node_groups:n_node_groups + n_elem_groups]]

    cursor = data_idx + 4
    node_group_names = [lines[cursor + i].strip() for i in range(n_node_groups)]
    cursor += n_node_groups
    cursor += n_elem_groups  # elem group名の行数分スキップ(値そのものは後段で読み飛ばす)

    if "DISPLACEMENT" not in node_group_names:
        raise ValueError(f"DISPLACEMENTがnode value groupに見つかりません: {node_group_names}")
    disp_group_idx = node_group_names.index("DISPLACEMENT")
    offset_before = sum(node_widths[:disp_group_idx])
    n_node_values_total = sum(node_widths)

    # 残り全トークンを平坦化して、"nodeID, 値×n_node_values_total" を
    # n_nodes回繰り返す構造として読む(改行位置に依存しない)。
    rest_text = "\n".join(lines[cursor:])
    tokens = rest_text.split()

    values: dict[int, tuple[float, float, float]] = {}
    pos = 0
    for _ in range(n_nodes):
        node_id = int(tokens[pos]); pos += 1
        vals = [float(tokens[pos + k]) for k in range(n_node_values_total)]
        pos += n_node_values_total
        ux, uy, uz = vals[offset_before:offset_before + 3]
        values[node_id] = (ux, uy, uz)

    if len(values) != n_nodes:
        raise ValueError(f"読み込んだ節点数({len(values)})が期待値({n_nodes})と不一致: {path}")
    return values
