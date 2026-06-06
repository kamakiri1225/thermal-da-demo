"""
OpenFOAM round-bar data assimilation results:
plot temperature history for every axial node.

Reads the CSV outputs produced by `da_main.py` and writes a compact
small-multiples figure comparing truth / DA-off / DA-on for each node.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


BASE_DIR = Path(__file__).parent
RESULTS_DIR = BASE_DIR / "results"
RESULTS_IMG_DIR = RESULTS_DIR / "img"
TRUTH_CSV = RESULTS_DIR / "truth" / "axial_temperature_history.csv"
MODEL_CSV = RESULTS_DIR / "model_only" / "axial_temperature_history.csv"
DA_CSV = RESULTS_DIR / "with_da" / "axial_temperature_history.csv"
OUT_PNG = RESULTS_IMG_DIR / "all_nodes_temperature_comparison.png"


def load_history(path: Path) -> tuple[np.ndarray, list[str]]:
    with path.open("r", encoding="utf-8") as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data, header


def parse_node_labels(header: list[str]) -> list[tuple[int, float]]:
    labels: list[tuple[int, float]] = []
    for name in header[1:]:
        # N12_x0.094m -> node=12, x=0.094
        left, x_part = name.split("_x", 1)
        node = int(left[1:])
        x_m = float(x_part[:-1])
        labels.append((node, x_m))
    return labels


def main() -> None:
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    truth, header = load_history(TRUTH_CSV)
    model, _ = load_history(MODEL_CSV)
    da, _ = load_history(DA_CSV)

    labels = parse_node_labels(header)
    time_s = truth[:, 0]

    truth_y = truth[:, 1:]
    model_y = model[:, 1:]
    da_y = da[:, 1:]

    all_vals = np.concatenate([truth_y.ravel(), model_y.ravel(), da_y.ravel()])
    ymin = float(all_vals.min())
    ymax = float(all_vals.max())
    pad = max(0.5, 0.08 * (ymax - ymin))

    n_nodes = truth_y.shape[1]
    ncols = 5
    nrows = int(np.ceil(n_nodes / ncols))

    ASSIM_NODES = {4, 16}
    VALIDATION_NODES = {28, 36}

    fig, axes = plt.subplots(nrows, ncols, figsize=(18, 22), sharex=True, sharey=True)
    fig.suptitle(
        "Round Bar Twin Experiment: each axial node\n"
        "truth / DA off / DA on",
        fontsize=13,
        y=0.995,
    )

    for idx in range(nrows * ncols):
        ax = axes.flat[idx]
        if idx >= n_nodes:
            ax.axis("off")
            continue

        node, x_m = labels[idx]
        ax.plot(time_s, truth_y[:, idx], color="black", ls="--", lw=1.4)
        ax.plot(time_s, model_y[:, idx], color="#377eb8", ls=":", lw=1.2)
        ax.plot(time_s, da_y[:, idx], color="crimson", ls="-", lw=1.5)

        if node in ASSIM_NODES:
            ax.set_facecolor("#FFF3CD")
            title = f"N{node} x={x_m:.3f} m [sensor]"
        else:
            title = f"N{node} x={x_m:.3f} m"
        ax.set_title(title, fontsize=8)
        ax.grid(True, alpha=0.25)
        ax.set_ylim(ymin - pad, ymax + pad)

        if idx % ncols == 0:
            ax.set_ylabel("Temperature [degC]", fontsize=8)
        if idx >= (nrows - 1) * ncols:
            ax.set_xlabel("Time [s]", fontsize=8)

    from matplotlib.patches import Patch
    handles = [
        Line2D([0], [0], color="black", ls="--", lw=1.6, label="truth"),
        Line2D([0], [0], color="#377eb8", ls=":", lw=1.6, label="DA off"),
        Line2D([0], [0], color="crimson", ls="-", lw=1.6, label="DA on"),
        Patch(facecolor="#FFF3CD", edgecolor="gray", lw=0.8, label="assimilation sensor"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.978), fontsize=9)

    plt.tight_layout(rect=(0.02, 0.02, 0.98, 0.957))
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {OUT_PNG}")


if __name__ == "__main__":
    main()
