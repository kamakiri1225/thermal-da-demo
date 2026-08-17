#!/usr/bin/env python3
"""Convert OpenFOAM solid probes and field extrema to CSV and PNG."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CASE = Path(__file__).resolve().parents[1]
PROBES = CASE / "postProcessing" / "solidTemperatureProbes" / "solid" / "0" / "T"
OUT_DIR = CASE / "data"


def scalar_values(path: Path) -> list[float]:
    text = path.read_text(encoding="utf-8")
    uniform = re.search(r"internalField\s+uniform\s+([-+0-9.eE]+)\s*;", text)
    if uniform:
        return [float(uniform.group(1))]
    field = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+\d+\s*\((.*?)\)\s*;",
        text,
        re.S,
    )
    if not field:
        raise ValueError(f"Cannot parse internalField: {path}")
    return [float(value) for value in field.group(1).split()]


def main() -> int:
    probe_rows: dict[float, list[float]] = {}
    for line in PROBES.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        probe_rows[values[0]] = values[1:]
    probe_rows.setdefault(0.0, [293.15, 293.15, 293.15, 293.15])

    rows = []
    for time, probes in sorted(probe_rows.items()):
        field = CASE / f"{time:g}" / "solid" / "T"
        if not field.exists():
            continue
        values = scalar_values(field)
        rows.append({
            "time_s": time,
            "T_hot_C": probes[0] - 273.15,
            "T_cold_C": probes[1] - 273.15,
            "T_mid_C": probes[2] - 273.15,
            "T_top_C": probes[3] - 273.15,
            "T_min_C": min(values) - 273.15,
            "T_max_C": max(values) - 273.15,
            "heater_power_W": 15.0 if time < 300.0 else 0.0,
        })

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = OUT_DIR / "temperature_history.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    plt.rcParams.update({
        "font.size": 20,
        "axes.labelsize": 22,
        "xtick.labelsize": 18,
        "ytick.labelsize": 18,
        "legend.fontsize": 18,
    })
    fig, axis = plt.subplots(figsize=(11, 6.5))
    for key, label in (
        ("T_hot_C", "T_hot (heater side)"),
        ("T_cold_C", "T_cold (opposite side)"),
        ("T_mid_C", "T_mid (90 deg)"),
        ("T_top_C", "T_top (near top)"),
    ):
        axis.plot(
            [r["time_s"] for r in rows],
            [r[key] for r in rows],
            linewidth=2.8,
            label=label,
        )
    axis.set_xlabel("Time [s]", labelpad=10)
    axis.set_ylabel("Temperature [degC]", labelpad=12)
    axis.grid(alpha=0.3)
    axis.legend(loc="best", frameon=True)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "temperature_history.png", dpi=160)
    print(f"Wrote {len(rows)} rows to {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
