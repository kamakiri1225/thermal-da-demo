#!/usr/bin/env python3
"""
plot_temperature_humidity.py
============================================================
ESP32ロガーが保存したCSVを読み込み、時系列グラフをPNG保存する。

使い方:
    python3 plot_temperature_humidity.py                     # 既定: ../data/sample_log.csv
    python3 plot_temperature_humidity.py path/to/LOG.CSV     # ファイル指定

出力:
    ../outputs/temperature_humidity_timeseries.png

特徴:
    - timestamp は「日時文字列」「経過秒（数値）」のどちらでも読める
    - NaN・欠損・明らかな異常値（熱電対の断線スパイク等）があっても落ちない
    - surface_temp_*_C 列は何本あっても自動で全部プロットする
============================================================
"""
import os
import sys

import matplotlib

matplotlib.use("Agg")  # 画面なし環境でも保存できるように
import matplotlib.pyplot as plt
import pandas as pd

# ── 設定 ────────────────────────────────────────────────
HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CSV = os.path.join(HERE, "..", "data", "sample_log.csv")
OUT_DIR = os.path.join(HERE, "..", "outputs")
OUT_PNG = os.path.join(OUT_DIR, "temperature_humidity_timeseries.png")

# 物理的にあり得ない値はNaN扱いにする（断線スパイク対策）
TEMP_VALID_RANGE = (-50.0, 500.0)   # 温度 [°C]
HUMID_VALID_RANGE = (0.0, 100.0)    # 湿度 [%]


def load_csv(path: str) -> pd.DataFrame:
    """CSVを読み込み、timestampと数値列を整える。多少壊れていても落ちない。"""
    df = pd.read_csv(path, na_values=["NaN", "nan", ""], skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]

    if "timestamp" not in df.columns:
        raise ValueError(f"'timestamp' 列がありません。列: {list(df.columns)}")

    # timestamp: 日時文字列 or 経過秒 の両対応
    ts_datetime = pd.to_datetime(df["timestamp"], errors="coerce")
    if ts_datetime.notna().mean() > 0.5:
        df["time"] = ts_datetime
        time_label = "time"
    else:
        # 経過秒（数値）とみなす
        df["elapsed_s"] = pd.to_numeric(df["timestamp"], errors="coerce")
        df["time"] = df["elapsed_s"]
        time_label = "elapsed time [s]"

    # 数値列を強制変換（文字化けなどはNaNへ）
    for col in df.columns:
        if col in ("timestamp", "time"):
            continue
        df[col] = pd.to_numeric(df[col], errors="coerce")

    # 異常値をNaN化
    for col in df.columns:
        if col.endswith("_C"):
            lo, hi = TEMP_VALID_RANGE
            df.loc[(df[col] < lo) | (df[col] > hi), col] = float("nan")
        if col == "humidity_percent":
            lo, hi = HUMID_VALID_RANGE
            df.loc[(df[col] < lo) | (df[col] > hi), col] = float("nan")

    # 時間が読めない行は捨てる
    df = df.dropna(subset=["time"]).reset_index(drop=True)
    df.attrs["time_label"] = time_label
    return df


def main() -> int:
    csv_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_CSV
    if not os.path.isfile(csv_path):
        print(f"[エラー] CSVが見つかりません: {csv_path}")
        return 1

    print(f"読み込み: {csv_path}")
    df = load_csv(csv_path)
    print(f"行数: {len(df)}  列: {list(df.columns)}")

    surface_cols = [c for c in df.columns if c.startswith("surface_temp_")]
    has_ambient = "ambient_temp_C" in df.columns
    has_humidity = "humidity_percent" in df.columns

    n_missing = int(df[surface_cols + (["ambient_temp_C"] if has_ambient else [])]
                    .isna().sum().sum())
    if n_missing:
        print(f"[情報] 欠損/異常値を {n_missing} 個検出しました（グラフでは途切れとして表示）")

    # ── 描画: 上=温度（表面+雰囲気）、下=湿度 ──
    fig, (ax_t, ax_h) = plt.subplots(
        2, 1, figsize=(11, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )

    colors = ["#d62728", "#ff7f0e", "#9467bd", "#8c564b"]
    for i, col in enumerate(surface_cols):
        ax_t.plot(df["time"], df[col], lw=1.6, color=colors[i % len(colors)],
                  label=col.replace("_C", " [°C]"))
    if has_ambient:
        ax_t.plot(df["time"], df["ambient_temp_C"], lw=1.6, color="#1f77b4",
                  label="ambient_temp [°C]")
    ax_t.set_ylabel("Temperature [°C]")
    ax_t.legend(loc="best", fontsize=9)
    ax_t.grid(alpha=0.3)
    ax_t.set_title("Surface / Ambient Temperature and Humidity Time History")

    if has_humidity:
        ax_h.plot(df["time"], df["humidity_percent"], lw=1.6, color="#2ca02c",
                  label="humidity [%]")
        ax_h.set_ylabel("Humidity [%]")
        ax_h.legend(loc="best", fontsize=9)
    ax_h.grid(alpha=0.3)
    ax_h.set_xlabel(df.attrs.get("time_label", "time"))

    fig.autofmt_xdate()
    fig.tight_layout()

    os.makedirs(OUT_DIR, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=140)
    print(f"保存しました: {OUT_PNG}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
