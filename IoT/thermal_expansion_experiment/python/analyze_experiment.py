#!/usr/bin/env python3
"""
analyze_experiment.py
============================================================
serial_logger.py が保存したCSVを読み込み、温度・変位の時系列、
理論熱膨張量との比較、誤差率、ヒステリシス等を計算・図示する。

使い方:
    python3 analyze_experiment.py path/to/thermal_expansion_XXXX.csv
    python3 analyze_experiment.py path/to/LOG.csv --config ../config/experiment_config.example.json \
        --material ../config/material_properties_sus304.json

出力:
    ../data/processed/<basename>_summary.json  (数値サマリ)
    ../data/processed/<basename>_timeseries.png
    ../data/processed/<basename>_temp_vs_disp.png

理論式 (docs/11_fem_comparison.md 参照):
    単純式:   dL = alpha * L * dT            (dT = 平均温度の変化量)
    積分近似: dL = alpha * integral_0^L [T(x,t) - T0(x)] dx
              (3点の位置と値から線形補間 + 台形積分。両端は最近傍値で外挿)
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
DEFAULT_CONFIG = os.path.join(HERE, "..", "config", "experiment_config.example.json")
DEFAULT_MATERIAL = os.path.join(HERE, "..", "config", "material_properties_sus304.json")
DEFAULT_OUTDIR = os.path.join(HERE, "..", "data", "processed")

TEMP_VALID_RANGE = (-50.0, 1000.0)


def load_json(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def load_csv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, comment="#", na_values=["NaN", "nan", ""], skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]

    required = ["elapsed_ms", "T1_C", "T2_C", "T3_C", "indicator_mm", "disp_mm", "spc_valid"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"必要な列がありません: {missing}. 実際の列: {list(df.columns)}")

    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    for col in ("T1_C", "T2_C", "T3_C"):
        lo, hi = TEMP_VALID_RANGE
        df.loc[(df[col] < lo) | (df[col] > hi), col] = np.nan

    df["elapsed_s"] = df["elapsed_ms"] / 1000.0
    return df


def average_temp_simple(df: pd.DataFrame) -> pd.Series:
    return df[["T1_C", "T2_C", "T3_C"]].mean(axis=1, skipna=True)


def average_temp_weighted(df: pd.DataFrame, positions_mm, length_mm: float) -> pd.Series:
    """3点の位置(positions_mm=[x1,x2,x3], 0=加熱端)を使い、線形補間+台形積分で
    軸方向の平均温度を近似する。両端(x=0, x=L)は最近傍値で外挿する簡易近似。"""
    x1, x2, x3 = positions_mm
    L = length_mm

    def integrate_row(t1, t2, t3):
        if any(pd.isna(v) for v in (t1, t2, t3)):
            return np.nan
        xs = [0.0, x1, x2, x3, L]
        ts = [t1, t1, t2, t3, t3]  # 両端は最近傍値で外挿(簡易近似、要検証)
        area = 0.0
        for i in range(len(xs) - 1):
            dx = xs[i + 1] - xs[i]
            area += 0.5 * (ts[i] + ts[i + 1]) * dx
        return area / L

    return df.apply(lambda row: integrate_row(row["T1_C"], row["T2_C"], row["T3_C"]), axis=1)


def theoretical_dL(avg_temp: pd.Series, alpha: float, length_m: float) -> pd.Series:
    t0 = avg_temp.iloc[0] if len(avg_temp) else np.nan
    dT = avg_temp - t0
    return alpha * length_m * dT


def detect_phase(df: pd.DataFrame, temp_col: str = "T1_C", window: int = 5) -> pd.Series:
    """T1の移動平均の傾きから、加熱(heating)/保持(hold)/冷却(cooling)を大まかに判定する。"""
    smooth = df[temp_col].rolling(window=window, min_periods=1, center=True).mean()
    slope = smooth.diff().fillna(0.0)
    phase = pd.Series("hold", index=df.index)
    phase[slope > 0.01] = "heating"
    phase[slope < -0.01] = "cooling"
    return phase


def compute_hysteresis(avg_temp: pd.Series, disp_mm: pd.Series, phase: pd.Series):
    """同程度の平均温度における昇温時disp_mmと冷却時disp_mmの差の最大値を概算する。"""
    heating = pd.DataFrame({"T": avg_temp[phase == "heating"], "d": disp_mm[phase == "heating"]}).dropna()
    cooling = pd.DataFrame({"T": avg_temp[phase == "cooling"], "d": disp_mm[phase == "cooling"]}).dropna()
    if heating.empty or cooling.empty:
        return None
    # 冷却側の各温度に最も近い昇温側の変位を線形補間で求め、差を取る
    heating_sorted = heating.sort_values("T")
    interp_d = np.interp(cooling["T"], heating_sorted["T"], heating_sorted["d"])
    diff = cooling["d"].values - interp_d
    return float(np.nanmax(np.abs(diff))) if len(diff) else None


def main() -> int:
    parser = argparse.ArgumentParser(description="熱膨張実験CSVの解析")
    parser.add_argument("csv_path")
    parser.add_argument("--config", default=DEFAULT_CONFIG)
    parser.add_argument("--material", default=DEFAULT_MATERIAL)
    parser.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    if not os.path.isfile(args.csv_path):
        print(f"[エラー] CSVが見つかりません: {args.csv_path}")
        return 1

    cfg = load_json(args.config) if os.path.isfile(args.config) else {}
    mat = load_json(args.material) if os.path.isfile(args.material) else {}

    alpha = mat.get("properties", {}).get("linear_expansion_coefficient", {}).get("value", 17.3e-6)
    length_mm = cfg.get("specimen", {}).get("length_mm", 300.0)
    positions = cfg.get("sensor_positions_mm_from_heated_end", {})
    x1 = positions.get("T1_heated_side", 20)
    x2 = positions.get("T2_center", length_mm / 2)
    x3 = positions.get("T3_free_end_side", length_mm - 20)

    print(f"読み込み: {args.csv_path}")
    df = load_csv(args.csv_path)
    print(f"行数: {len(df)}")

    os.makedirs(args.outdir, exist_ok=True)
    base = os.path.splitext(os.path.basename(args.csv_path))[0]

    avg_simple = average_temp_simple(df)
    avg_weighted = average_temp_weighted(df, [x1, x2, x3], length_mm)

    dL_simple = theoretical_dL(avg_simple, alpha, length_mm / 1000.0) * 1000.0   # mm
    dL_weighted = theoretical_dL(avg_weighted, alpha, length_mm / 1000.0) * 1000.0  # mm

    disp = df["disp_mm"]

    # 実測 vs 理論の誤差(ピーク時)
    idx_peak = disp.abs().idxmax() if disp.notna().any() else None
    summary = {
        "csv_file": os.path.basename(args.csv_path),
        "n_rows": int(len(df)),
        "alpha_used_per_K": alpha,
        "length_mm_used": length_mm,
        "sensor_positions_mm": {"T1": x1, "T2": x2, "T3": x3},
    }

    if idx_peak is not None:
        meas_peak = float(disp.loc[idx_peak])
        th_simple_peak = float(dL_simple.loc[idx_peak]) if not pd.isna(dL_simple.loc[idx_peak]) else None
        th_weighted_peak = float(dL_weighted.loc[idx_peak]) if not pd.isna(dL_weighted.loc[idx_peak]) else None
        summary["peak_measured_disp_mm"] = meas_peak
        summary["peak_theoretical_simple_mm"] = th_simple_peak
        summary["peak_theoretical_weighted_mm"] = th_weighted_peak
        if th_simple_peak not in (None, 0):
            summary["error_rate_vs_simple_percent"] = (meas_peak - th_simple_peak) / th_simple_peak * 100.0
        if th_weighted_peak not in (None, 0):
            summary["error_rate_vs_weighted_percent"] = (meas_peak - th_weighted_peak) / th_weighted_peak * 100.0

    # 欠損・異常カウント
    summary["missing_T1"] = int(df["T1_C"].isna().sum())
    summary["missing_T2"] = int(df["T2_C"].isna().sum())
    summary["missing_T3"] = int(df["T3_C"].isna().sum())
    summary["missing_disp"] = int(disp.isna().sum())
    for col in ("fault1", "fault2", "fault3"):
        if col in df.columns:
            summary[f"{col}_count"] = int(df[col].fillna(0).sum())
    if "spc_valid" in df.columns:
        summary["spc_invalid_count"] = int((df["spc_valid"] == 0).sum())

    # ヒステリシス
    phase = detect_phase(df)
    hyst = compute_hysteresis(avg_simple, disp, phase)
    summary["hysteresis_max_abs_mm"] = hyst

    summary_path = os.path.join(args.outdir, f"{base}_summary.json")
    with open(summary_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    print(f"サマリを保存しました: {summary_path}")
    print(json.dumps(summary, ensure_ascii=False, indent=2))

    # ── 図1: 温度・変位の時系列 ──
    fig, (ax_t, ax_d) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
    ax_t.plot(df["elapsed_s"], df["T1_C"], label="T1 (heated side)", color="#d62728")
    ax_t.plot(df["elapsed_s"], df["T2_C"], label="T2 (center)", color="#ff7f0e")
    ax_t.plot(df["elapsed_s"], df["T3_C"], label="T3 (free end side)", color="#1f77b4")
    ax_t.set_ylabel("Temperature [C]")
    ax_t.legend(fontsize=9)
    ax_t.grid(alpha=0.3)
    ax_t.set_title("Temperature & Displacement Time History")

    ax_d.plot(df["elapsed_s"], disp, label="measured disp_mm", color="#2ca02c")
    ax_d.plot(df["elapsed_s"], dL_simple, label="theoretical (simple avg)", color="#9467bd", linestyle="--")
    ax_d.plot(df["elapsed_s"], dL_weighted, label="theoretical (weighted/integral)", color="#8c564b", linestyle=":")
    ax_d.set_xlabel("elapsed time [s]")
    ax_d.set_ylabel("Displacement [mm]")
    ax_d.legend(fontsize=9)
    ax_d.grid(alpha=0.3)

    fig.tight_layout()
    ts_png = os.path.join(args.outdir, f"{base}_timeseries.png")
    fig.savefig(ts_png, dpi=140)
    print(f"保存しました: {ts_png}")

    # ── 図2: 平均温度 vs 変位 (ヒステリシス確認用) ──
    fig2, ax2 = plt.subplots(figsize=(7, 6))
    colors = {"heating": "#d62728", "cooling": "#1f77b4", "hold": "#7f7f7f"}
    for ph, color in colors.items():
        mask = phase == ph
        ax2.scatter(avg_simple[mask], disp[mask], s=8, color=color, label=ph, alpha=0.7)
    ax2.set_xlabel("Average Temperature (simple mean) [C]")
    ax2.set_ylabel("Measured Displacement [mm]")
    ax2.set_title("Temperature vs Displacement (heating/cooling hysteresis)")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig2.tight_layout()
    hy_png = os.path.join(args.outdir, f"{base}_temp_vs_disp.png")
    fig2.savefig(hy_png, dpi=140)
    print(f"保存しました: {hy_png}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
