#!/usr/bin/env python3
"""
run_thermal_expansion.py
============================================================
101_0(OpenFOAM chtMultiRegionFoam, CHT+輻射+ヒートマット発熱源)の
固体温度分布を取り込み、101_1(FrontISTR)で熱膨張(線形静解析、
底面固定)を計算する一連の処理をまとめて実行する。

前提: OpenFOAM環境がsource済みであること(postProcessを呼ぶため)。
      101_0側で最低1回 Allrun (または Allrun.pre + 本体実行)が
      済んでおり、使いたい時刻のフィールドが書き出されていること。

使い方:
    python3 run_thermal_expansion.py \
        --of-case ../../101_0_openfoam_cht_radiation_box \
        --time latestTime

出力:
    case/box_thermal_expansion.{msh,cnt,res.0.*}
    ../data/summary.yaml (変位の要約)
============================================================
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from box_mesh import build_box_mesh  # noqa: E402
from openfoam_temperature import (  # noqa: E402
    align_cell_centers_to_node_mesh,
    interpolate_to_nodes,
    load_solid_cell_temperatures,
)
import fistr_case  # noqa: E402

DEFAULT_OF_CASE = HERE / ".." / ".." / "101_0_openfoam_cht_radiation_box"
DEFAULT_CASE_DIR = HERE / ".." / "case"
DEFAULT_MATERIAL = HERE / ".." / "config" / "material_properties_steel.yaml"

# 101_0側 system/include/caseSettings のデフォルト値と揃えている(要連動確認)。
SOLID_LX = 0.2
SOLID_LY = 0.2
SOLID_LZ = 0.4
CELL_SIZE = 0.05
DEFAULT_NX = round(SOLID_LX / CELL_SIZE)   # 4
DEFAULT_NY = round(SOLID_LY / CELL_SIZE)   # 4
DEFAULT_NZ = round(SOLID_LZ / CELL_SIZE)   # 8


def resolve_latest_time(of_case: Path) -> str:
    times = []
    for p in of_case.iterdir():
        if p.is_dir() and re.match(r"^\d+(\.\d+)?$", p.name) and p.name != "0.orig":
            times.append(p.name)
    if not times:
        raise FileNotFoundError(f"{of_case} に時刻ディレクトリが見つかりません。先に101_0を実行してください。")
    times.sort(key=lambda s: float(s))
    return times[-1]


def ensure_cell_centres(of_case: Path, time_dir: str) -> None:
    """solid + heaterMat の両リージョンについて、セル中心座標(C)が
    無ければ postProcess -func writeCellCentres を実行して作る。
    heaterMatリージョンが存在しない(旧構成の)ケースはスキップする。"""
    for region in ("solid", "heaterMat"):
        region_dir = of_case / time_dir / region
        if not region_dir.is_dir():
            continue
        c_path = region_dir / "C"
        if c_path.exists():
            continue
        print(f"[INFO] {c_path} が無いため postProcess -func writeCellCentres を実行します。")
        result = subprocess.run(
            ["postProcess", "-func", "writeCellCentres", "-region", region, "-time", time_dir],
            cwd=of_case,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        print(result.stdout[-3000:])
        if result.returncode != 0 or not c_path.exists():
            raise RuntimeError(
                f"postProcess -func writeCellCentres -region {region} に失敗しました。"
                "OpenFOAM環境がsourceされているか確認してください。"
            )


def run_one_time(
    of_case: Path,
    time_dir: str,
    case_dir: Path,
    mat: dict,
    nx: int,
    ny: int,
    nz: int,
) -> dict:
    """101_0の1時刻分についてFrontISTR熱膨張解析を実行し、サマリdictを返す。
    101_1のワークフロー(単一時刻: run_thermal_expansion.py、
    複数時刻: run_thermal_expansion_timehistory.py)の両方から使う共通処理。"""
    ensure_cell_centres(of_case, time_dir)
    centers, temps = load_solid_cell_temperatures(of_case, time_dir)
    print(f"[INFO] t={time_dir}: 固体セル数={len(temps)}  T範囲: {temps.min():.3f} - {temps.max():.3f} K")

    mesh = fistr_case.write_mesh(
        case_dir, nx, ny, nz, SOLID_LX, SOLID_LY, SOLID_LZ,
        young_modulus=mat["young_modulus_Pa"],
        poisson_ratio=mat["poisson_ratio"],
        density=mat["density_kg_m3"],
        thermal_expansion_coeff=mat["thermal_expansion_coeff_per_K"],
    )
    fistr_case.write_hecmw_ctrl(case_dir)

    node_ids = [nid for nid, _xyz in mesh["nodes"]]
    node_coords = np.array([xyz for _nid, xyz in mesh["nodes"]], dtype=float)

    aligned_centers, coordinate_translation = align_cell_centers_to_node_mesh(centers, node_coords)
    print(
        "[INFO] OpenFOAM -> FrontISTR 座標平行移動 [m]: "
        f"({coordinate_translation[0]:.6g}, {coordinate_translation[1]:.6g}, "
        f"{coordinate_translation[2]:.6g})"
    )
    node_temps = interpolate_to_nodes(node_coords, aligned_centers, temps, k=8)

    fistr_case.write_cnt(
        case_dir,
        node_ids,
        node_temps,
        reference_temperature=mat["reference_temperature_K"],
        young_modulus=mat["young_modulus_Pa"],
        poisson_ratio=mat["poisson_ratio"],
        thermal_expansion_coeff=mat["thermal_expansion_coeff_per_K"],
    )

    fistr_case.run_fistr(case_dir)
    fistr_case.add_temperature_to_visualization(
        case_dir,
        node_coords,
        node_temps,
        time_value=float(time_dir),
    )

    disp = fistr_case.read_displacement(case_dir)

    top_nodes = mesh["top_nodes"]
    top_uz = np.array([disp[n][2] for n in top_nodes if n in disp])
    top_ux = np.array([disp[n][0] for n in top_nodes if n in disp])
    all_mag = np.array([np.linalg.norm(v) for v in disp.values()])

    return {
        "of_case": str(of_case),
        "of_time": time_dir,
        "of_time_s": float(time_dir),
        "n_solid_cells": int(len(temps)),
        "solid_T_min_K": float(temps.min()),
        "solid_T_max_K": float(temps.max()),
        "node_T_min_K": float(node_temps.min()),
        "node_T_max_K": float(node_temps.max()),
        "coordinate_translation_x_m": float(coordinate_translation[0]),
        "coordinate_translation_y_m": float(coordinate_translation[1]),
        "coordinate_translation_z_m": float(coordinate_translation[2]),
        "reference_temperature_K": mat["reference_temperature_K"],
        "top_face_mean_Uz_mm": float(top_uz.mean() * 1000.0),
        "top_face_max_Uz_mm": float(top_uz.max() * 1000.0),
        "top_face_mean_Ux_mm": float(top_ux.mean() * 1000.0),
        "top_face_ux_spread_mm": float((top_ux.max() - top_ux.min()) * 1000.0),
        "max_displacement_magnitude_mm": float(all_mag.max() * 1000.0),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="101_0の固体温度分布をFrontISTRへ引き継いで熱膨張を計算する(単一時刻)")
    parser.add_argument("--of-case", default=str(DEFAULT_OF_CASE))
    parser.add_argument("--time", default="latestTime")
    parser.add_argument("--case-dir", default=str(DEFAULT_CASE_DIR))
    parser.add_argument("--material", default=str(DEFAULT_MATERIAL))
    parser.add_argument("--nx", type=int, default=DEFAULT_NX)
    parser.add_argument("--ny", type=int, default=DEFAULT_NY)
    parser.add_argument("--nz", type=int, default=DEFAULT_NZ)
    args = parser.parse_args()

    of_case = Path(args.of_case).resolve()
    case_dir = Path(args.case_dir).resolve()
    case_dir.mkdir(parents=True, exist_ok=True)

    mat = yaml.safe_load(Path(args.material).read_text(encoding="utf-8"))

    time_dir = args.time
    if time_dir == "latestTime":
        time_dir = resolve_latest_time(of_case)
    print(f"[INFO] OpenFOAM時刻: {time_dir}  (case: {of_case})")

    summary = run_one_time(of_case, time_dir, case_dir, mat, args.nx, args.ny, args.nz)

    out_dir = case_dir.parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.yaml"
    summary_text = yaml.safe_dump(summary, allow_unicode=True, sort_keys=False)
    summary_path.write_text(summary_text, encoding="utf-8")

    print(summary_text, end="")
    print(f"[INFO] サマリを保存しました: {summary_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
