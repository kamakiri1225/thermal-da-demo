#!/usr/bin/env python3
"""
plot_results.py
============================================================
実験直後の「とりあえず波形を見る」ための軽量プロットスクリプト。
analyze_experiment.py のような理論値比較・サマリ計算は行わず、
生CSVの時系列だけを素早く画像化する。

使い方:
    python3 plot_results.py path/to/thermal_expansion_XXXX.csv
    python3 plot_results.py path/to/thermal_expansion_XXXX.csv --outdir ../data/processed
============================================================
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUTDIR = os.path.join(HERE, "..", "data", "processed")


def main() -> int:
    parser = argparse.ArgumentParser(description="実験CSVの簡易プロット(速報用)")
    parser.add_argument("csv_path")
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    if not os.path.isfile(args.csv_path):
        print(f"[エラー] CSVが見つかりません: {args.csv_path}")
        return 1

    df = pd.read_csv(args.csv_path, comment="#", na_values=["NaN", "nan", ""])
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    if "elapsed_ms" not in df.columns:
        print(f"[エラー] elapsed_ms列がありません。列: {list(df.columns)}")
        return 1
    t = df["elapsed_ms"] / 1000.0

    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.csv_path))[0]

    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    temp_cols = [c for c in ("T1_C", "T2_C", "T3_C") if c in df.columns]
    for c in temp_cols:
        axes[0].plot(t, df[c], label=c)
    axes[0].set_ylabel("Temperature [C]")
    axes[0].legend(fontsize=9)
    axes[0].grid(alpha=0.3)
    axes[0].set_title(f"Quick look: {os.path.basename(args.csv_path)}")

    if "indicator_mm" in df.columns:
        axes[1].plot(t, df["indicator_mm"], color="#555555")
    axes[1].set_ylabel("indicator_mm (raw)")
    axes[1].grid(alpha=0.3)

    if "disp_mm" in df.columns:
        axes[2].plot(t, df["disp_mm"], color="#2ca02c")
    axes[2].set_ylabel("disp_mm (zero-referenced)")
    axes[2].set_xlabel("elapsed time [s]")
    axes[2].grid(alpha=0.3)

    fig.tight_layout()
    out_png = os.path.join(args.outdir, f"{base}_quicklook.png")
    fig.savefig(out_png, dpi=140)
    print(f"保存しました: {out_png}")

    if "fault1" in df.columns:
        n_fault = int(df[["fault1", "fault2", "fault3"]].fillna(0).sum().sum()) if "fault3" in df.columns else None
        if n_fault is not None:
            print(f"[情報] センサ異常フラグの合計: {n_fault}")
    if "spc_valid" in df.columns:
        print(f"[情報] spc_valid=0 の行数: {int((df['spc_valid'] == 0).sum())}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
