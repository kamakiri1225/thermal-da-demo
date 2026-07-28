#!/usr/bin/env python3
"""
run_thermal_expansion_timehistory.py
============================================================
101_0(OpenFOAM)の複数時刻分の固体温度分布を順番にFrontISTRへ引き継ぎ、
熱膨張(変位)の時刻歴を求める。run_thermal_expansion.py(単一時刻)を
101_0の全時刻ディレクトリに対して繰り返し実行するラッパー。

使い方:
    python3 run_thermal_expansion_timehistory.py \
        --of-case ../../101_0_openfoam_cht_radiation_box

    # 時刻を間引きたい場合(全時刻だと計算時間がかかるため)
    python3 run_thermal_expansion_timehistory.py --every 2

出力:
    case/t_<time>/box_thermal_expansion.{msh,cnt,res.0.*}  (時刻ごとのFrontISTRケース)
    ../data/timehistory.csv   (時刻歴サマリ)
    ../data/timehistory.png   (変位・温度の時刻歴プロット)
============================================================
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import yaml  # noqa: E402

from run_thermal_expansion import (  # noqa: E402
    DEFAULT_MATERIAL,
    DEFAULT_NX,
    DEFAULT_NY,
    DEFAULT_NZ,
    DEFAULT_OF_CASE,
    run_one_time,
)

DEFAULT_CASE_ROOT = HERE / ".." / "case"


def list_time_dirs(of_case: Path) -> list[str]:
    times = []
    for p in of_case.iterdir():
        if p.is_dir() and re.match(r"^\d+(\.\d+)?$", p.name) and p.name != "0.orig":
            times.append(p.name)
    times.sort(key=lambda s: float(s))
    return times


def main() -> int:
    parser = argparse.ArgumentParser(description="101_0の時刻歴温度分布からFrontISTR熱膨張の時刻歴を計算する")
    parser.add_argument("--of-case", default=str(DEFAULT_OF_CASE))
    parser.add_argument("--case-root", default=str(DEFAULT_CASE_ROOT),
                         help="時刻ごとのFrontISTRケースを t_<time>/ として作る親フォルダ")
    parser.add_argument("--material", default=str(DEFAULT_MATERIAL))
    parser.add_argument("--nx", type=int, default=DEFAULT_NX)
    parser.add_argument("--ny", type=int, default=DEFAULT_NY)
    parser.add_argument("--nz", type=int, default=DEFAULT_NZ)
    parser.add_argument("--every", type=int, default=1, help="何個おきに時刻を処理するか(間引き、既定1=全時刻)")
    parser.add_argument("--exclude-t0", action="store_true", help="t=0(初期一様温度)を除外する(既定は含める)")
    parser.add_argument("--start-time", type=float, help="処理する最初のOpenFOAM時刻[s]")
    parser.add_argument("--end-time", type=float, help="処理する最後のOpenFOAM時刻[s]")
    args = parser.parse_args()

    of_case = Path(args.of_case).resolve()
    case_root = Path(args.case_root).resolve()
    case_root.mkdir(parents=True, exist_ok=True)
    mat = yaml.safe_load(Path(args.material).read_text(encoding="utf-8"))

    times = list_time_dirs(of_case)
    if args.exclude_t0:
        times = [t for t in times if float(t) > 0]
    if args.start_time is not None:
        times = [t for t in times if float(t) >= args.start_time]
    if args.end_time is not None:
        times = [t for t in times if float(t) <= args.end_time]
    times = times[:: args.every]

    if not times:
        print(f"[エラー] {of_case} に処理対象の時刻ディレクトリがありません。")
        return 1

    print(f"[INFO] 対象時刻 ({len(times)}個): {times}")
    ensure_cell_centres_for_times(of_case, times)

    rows: list[dict] = []
    for time_dir in times:
        case_dir = case_root / f"t_{time_dir}"
        case_dir.mkdir(parents=True, exist_ok=True)
        try:
            summary = run_one_time(of_case, time_dir, case_dir, mat, args.nx, args.ny, args.nz)
        except Exception as e:  # noqa: BLE001
            print(f"[警告] t={time_dir} の処理に失敗しました: {e}")
            continue
        rows.append(summary)
        (case_dir / "summary.yaml").write_text(
            yaml.safe_dump(summary, allow_unicode=True, sort_keys=False), encoding="utf-8"
        )
        print(
            f"  t={time_dir:>8s}s  T=[{summary['solid_T_min_K']:.3f},{summary['solid_T_max_K']:.3f}]K"
            f"  top_Uz_mean={summary['top_face_mean_Uz_mm']:.6f}mm"
            f"  top_Ux_spread={summary['top_face_ux_spread_mm']:.6f}mm"
        )

    if not rows:
        print("[エラー] すべての時刻でFrontISTR実行に失敗しました。")
        return 1

    out_dir = case_root.parent / "data"
    out_dir.mkdir(parents=True, exist_ok=True)

    csv_path = out_dir / "timehistory.csv"
    fieldnames = list(rows[0].keys())
    with open(csv_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    print(f"[INFO] 時刻歴CSVを保存しました: {csv_path}")

    pvd_path = out_dir / "thermal_expansion_timehistory.pvd"
    write_pvd_collection(rows, case_root, pvd_path)
    print(f"[INFO] ParaView時刻歴を保存しました: {pvd_path}")

    plot_path = out_dir / "timehistory.png"
    try:
        _plot_timehistory(rows, plot_path)
        print(f"[INFO] 時刻歴プロットを保存しました: {plot_path}")
    except Exception as e:  # noqa: BLE001
        print(f"[警告] プロット作成に失敗しました(CSVは正常に保存済み): {e}")

    return 0


def ensure_cell_centres_for_times(of_case: Path, times: list[str]) -> None:
    """不足しているセル中心座標Cをリージョンごとに一括生成する。"""
    for region in ("solid", "heaterMat"):
        missing = [
            time_dir
            for time_dir in times
            if (of_case / time_dir / region).is_dir()
            and not (of_case / time_dir / region / "C").exists()
        ]
        if not missing:
            continue

        time_selector = ",".join(missing)
        print(f"[INFO] {region}: {len(missing)}時刻分のセル中心座標Cを一括生成します。")
        result = subprocess.run(
            ["postProcess", "-func", "writeCellCentres", "-region", region, "-time", time_selector],
            cwd=of_case,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"postProcess -region {region} に失敗しました。\n{result.stdout[-3000:]}"
            )
        still_missing = [
            t for t in missing if not (of_case / t / region / "C").exists()
        ]
        if still_missing:
            raise RuntimeError(f"{region} のセル中心Cが生成されませんでした: {still_missing}")


def write_pvd_collection(rows: list[dict], case_root: Path, out_path: Path) -> None:
    """時刻別FrontISTR PVTUをOpenFOAM時刻付きPVDへまとめる。"""
    datasets = []
    for row in rows:
        case_dir = case_root / f"t_{row['of_time']}"
        candidates = sorted(case_dir.glob("box_thermal_expansion_vis_psf.*.pvtu"))
        if not candidates:
            raise FileNotFoundError(f"可視化PVTUがありません: {case_dir}")
        result_path = candidates[-1]
        relative_path = Path(os.path.relpath(result_path, out_path.parent)).as_posix()
        datasets.append(
            f'    <DataSet timestep="{row["of_time_s"]:.12g}" group="" part="0" '
            f'file="{relative_path}"/>'
        )

    content = [
        '<?xml version="1.0"?>',
        '<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">',
        "  <Collection>",
        *datasets,
        "  </Collection>",
        "</VTKFile>",
    ]
    out_path.write_text("\n".join(content) + "\n", encoding="utf-8")


def _plot_timehistory(rows: list[dict], out_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = [r["of_time_s"] for r in rows]
    solid_tmax = [r["solid_T_max_K"] for r in rows]
    uz = [r["top_face_mean_Uz_mm"] for r in rows]
    ux_spread = [r["top_face_ux_spread_mm"] for r in rows]

    fig, (ax_t, ax_d) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

    ax_t.plot(t, solid_tmax, marker="o", ms=3, color="#d62728", label="solid T max [K]")
    ax_t.set_ylabel("Temperature [K]")
    ax_t.legend(fontsize=9)
    ax_t.grid(alpha=0.3)
    ax_t.set_title("FrontISTR thermal expansion time history (driven by 101_0 OpenFOAM T)")

    ax_d.plot(t, uz, marker="o", ms=3, color="#2ca02c", label="top face mean Uz [mm]")
    ax_d.plot(t, ux_spread, marker="s", ms=3, color="#1f77b4", label="top face Ux spread [mm] (warpage)")
    ax_d.set_xlabel("OpenFOAM time [s]")
    ax_d.set_ylabel("Displacement [mm]")
    ax_d.legend(fontsize=9)
    ax_d.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)


if __name__ == "__main__":
    sys.exit(main())
