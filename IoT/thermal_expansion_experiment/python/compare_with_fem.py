#!/usr/bin/env python3
"""
compare_with_fem.py
============================================================
実験CSV(serial_logger.py出力)と、FEM/熱構造解析側の結果CSVを
時間軸で突き合わせ、変位・温度の差分と誤差率を計算する。

FEM側のCSV境界条件・出力形式は本プロジェクトでは未確定(docs/11_fem_comparison.md参照)。
そのため、本スクリプトは以下の「想定フォーマット」を前提とし、
実際のFEM出力に合わせて列名オプション(--fem-time-col等)で調整できるようにしている。

想定するFEM側CSVフォーマット(例、列名は変更可能):
    time_s, disp_mm, T1_C, T2_C, T3_C
    0.0,    0.000,   25.0, 25.0, 25.0
    1.0,    0.001,   26.0, 25.2, 25.1
    ...

使い方:
    python3 compare_with_fem.py \
        --exp ../data/raw/thermal_expansion_XXXX.csv \
        --fem ../data/processed/fem_result.csv

出力:
    ../data/processed/<exp_basename>_vs_fem.png
    ../data/processed/<exp_basename>_vs_fem_summary.json
============================================================
"""
import argparse
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTDIR = os.path.join(HERE, "..", "data", "processed")


def load_experiment_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", na_values=["NaN", "nan", ""])
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "elapsed_ms" not in df.columns:
        raise ValueError(f"実験CSVに elapsed_ms 列がありません。列: {list(df.columns)}")
    df["elapsed_s"] = df["elapsed_ms"] / 1000.0
    return df


def load_fem_csv(path: str, time_col: str) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", na_values=["NaN", "nan", ""])
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if time_col not in df.columns:
        raise ValueError(f"FEM CSVに時間列 '{time_col}' がありません。列: {list(df.columns)}。"
                          f"--fem-time-col で実際の列名を指定してください。")
    return df


def main() -> int:
    parser = argparse.ArgumentParser(description="実験結果とFEM/理論結果の比較")
    parser.add_argument("--exp", required=True, help="実験CSV(serial_logger.py出力)")
    parser.add_argument("--fem", required=True, help="FEM/理論結果CSV(想定フォーマットは本ファイル冒頭コメント参照)")
    parser.add_argument("--fem-time-col", default="time_s")
    parser.add_argument("--fem-disp-col", default="disp_mm")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    for p in (args.exp, args.fem):
        if not os.path.isfile(p):
            print(f"[エラー] ファイルが見つかりません: {p}")
            return 1

    exp = load_experiment_csv(args.exp)
    fem = load_fem_csv(args.fem, args.fem_time_col)

    if args.fem_disp_col not in fem.columns:
        print(f"[エラー] FEM CSVに変位列 '{args.fem_disp_col}' がありません。列: {list(fem.columns)}")
        return 1

    # FEM側の時系列を、実験側の時刻へ線形補間して揃える
    fem_sorted = fem.sort_values(args.fem_time_col)
    fem_disp_interp = np.interp(
        exp["elapsed_s"], fem_sorted[args.fem_time_col], fem_sorted[args.fem_disp_col],
        left=np.nan, right=np.nan,
    )

    exp_disp = exp["disp_mm"] if "disp_mm" in exp.columns else pd.Series(np.nan, index=exp.index)
    diff = exp_disp - fem_disp_interp

    valid = (~exp_disp.isna()) & (~pd.isna(fem_disp_interp))
    idx_peak = exp_disp.abs().idxmax() if exp_disp.notna().any() else None

    summary = {
        "exp_file": os.path.basename(args.exp),
        "fem_file": os.path.basename(args.fem),
        "n_compared_points": int(valid.sum()),
    }
    if idx_peak is not None and valid.loc[idx_peak]:
        meas = float(exp_disp.loc[idx_peak])
        fem_v = float(fem_disp_interp[exp.index.get_loc(idx_peak)])
        summary["peak_measured_mm"] = meas
        summary["peak_fem_mm"] = fem_v
        summary["peak_diff_mm"] = meas - fem_v
        if fem_v != 0:
            summary["peak_error_rate_percent"] = (meas - fem_v) / fem_v * 100.0

    if valid.sum() > 0:
        summary["rmse_mm"] = float(np.sqrt(np.nanmean(diff[valid] ** 2)))
        summary["max_abs_diff_mm"] = float(np.nanmax(np.abs(diff[valid])))

    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.exp))[0]

    summary_path = os.path.join(args.outdir, f"{base}_vs_fem_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    print(f"サマリを保存しました: {summary_path}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    ax1.plot(exp["elapsed_s"], exp_disp, label="measured", color="#2ca02c")
    ax1.plot(exp["elapsed_s"], fem_disp_interp, label="FEM/theoretical", color="#9467bd", linestyle="--")
    ax1.set_ylabel("Displacement [mm]")
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax1.set_title("Measured vs FEM Displacement")

    ax2.plot(exp["elapsed_s"], diff, color="#d62728")
    ax2.axhline(0, color="black", lw=0.8)
    ax2.set_xlabel("elapsed time [s]")
    ax2.set_ylabel("measured - FEM [mm]")
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    out_png = os.path.join(args.outdir, f"{base}_vs_fem.png")
    fig.savefig(out_png, dpi=140)
    print(f"保存しました: {out_png}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
