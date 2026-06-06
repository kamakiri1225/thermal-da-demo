"""
Create clear comparison plots from FrontISTR DA CSV results.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE_DIR = Path(__file__).resolve().parent.parent
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_IMG_DIR = RESULTS_DIR / "img"

MEASUREMENT_NODES = [4, 16, 28, 36]
ASSIM_NODES = {4, 16}


def read_csv(path: Path) -> np.ndarray:
    """Read a CSV file with headers into a structured NumPy array."""
    return np.genfromtxt(path, delimiter=",", names=True)


def plot_measurement_points() -> None:
    """Plot truth / model-only / DA histories at the representative sensor nodes."""
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    data = read_csv(RESULTS_DIR / "measurement_node_temperature_comparison.csv")
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    fig.suptitle(
        "FrontISTR Round Bar: Measurement Point Comparison\n"
        "truth / without data assimilation / with data assimilation",
        fontsize=12,
    )
    y_min = min(float(np.min(data[key])) for key in ["truth", "without_da", "with_da"])
    y_max = max(float(np.max(data[key])) for key in ["truth", "without_da", "with_da"])
    y_margin = max((y_max - y_min) * 0.05, 0.05)

    for ax, node in zip(axes.ravel(), MEASUREMENT_NODES):
        rows = data[data["axial_node"].astype(int) == node]
        time_min = rows["time_s"] / 60.0
        role = "assimilation sensor" if node in ASSIM_NODES else "validation point"
        ax.plot(time_min, rows["truth"], color="black", ls="--", lw=2.0, label="truth")
        ax.plot(time_min, rows["without_da"], color="steelblue", ls=":", lw=2.0, label="without DA")
        ax.plot(time_min, rows["with_da"], color="crimson", lw=1.8, label="with DA")
        ax.set_title(f"N{node}  x={rows['x_m'][0]:.3f} m  ({role})", fontsize=10)
        ax.set_ylabel("Temperature [degC]")
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8)

    for ax in axes[-1]:
        ax.set_xlabel("Time [min]")

    plt.tight_layout()
    out = RESULTS_IMG_DIR / "measurement_points_frontistr_da.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")


ASSIM_NODES = {4, 16}
VALIDATION_NODES = {28, 36}


def plot_all_axial_nodes_grid() -> None:
    """Plot every axial node in a dense panel for quick inspection of the time history."""
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    data = read_csv(RESULTS_DIR / "axial_node_temperature_comparison.csv")
    y_min = min(float(np.min(data[key])) for key in ["truth", "without_da", "with_da"])
    y_max = max(float(np.max(data[key])) for key in ["truth", "without_da", "with_da"])
    y_margin = max((y_max - y_min) * 0.05, 0.05)

    fig, axes = plt.subplots(8, 5, figsize=(18, 22), sharex=True, sharey=True)
    fig.suptitle(
        "FrontISTR Round Bar: All Axial Nodes\n"
        "truth / without data assimilation / with data assimilation",
        fontsize=13,
        y=0.995,
    )
    for node, ax in enumerate(axes.ravel()):
        rows = data[data["axial_node"].astype(int) == node]
        time_s = rows["time_s"]
        ax.plot(time_s, rows["truth"], color="black", ls="--", lw=1.4)
        ax.plot(time_s, rows["without_da"], color="steelblue", ls=":", lw=1.4)
        ax.plot(time_s, rows["with_da"], color="crimson", lw=1.2)

        if node in ASSIM_NODES:
            ax.set_facecolor("#FFF3CD")
            label = f"N{node} x={rows['x_m'][0]:.3f}m [sensor]"
        else:
            label = f"N{node} x={rows['x_m'][0]:.3f}m"
        ax.set_title(label, fontsize=8)
        ax.set_ylim(y_min - y_margin, y_max + y_margin)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        ax.grid(True, alpha=0.25)

    handles = [
        plt.Line2D([0], [0], color="black", ls="--", lw=1.8, label="truth"),
        plt.Line2D([0], [0], color="steelblue", ls=":", lw=1.8, label="without DA"),
        plt.Line2D([0], [0], color="crimson", lw=1.8, label="with DA"),
        plt.Rectangle((0, 0), 1, 1, fc="#FFF3CD", ec="gray", lw=0.8, label="assimilation sensor"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 0.980), fontsize=9)
    fig.text(0.5, 0.025, "Time [s]", ha="center", fontsize=11)
    fig.text(0.03, 0.5, "Temperature [degC]", va="center", rotation="vertical", fontsize=11)
    plt.tight_layout(rect=(0.04, 0.04, 0.99, 0.960))
    out = RESULTS_IMG_DIR / "all_axial_nodes_timeseries_frontistr_da.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")


def plot_axial_heatmap() -> None:
    """Plot temperature heatmaps for truth, model-only, and DA histories."""
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    data = read_csv(RESULTS_DIR / "axial_node_temperature_comparison.csv")
    time_values = np.unique(data["time_s"])
    node_values = np.unique(data["axial_node"]).astype(int)
    time_index = {float(t): idx for idx, t in enumerate(time_values)}

    fields = [
        ("truth", "truth"),
        ("without_da", "without DA"),
        ("with_da", "with DA"),
    ]
    arrays = []
    for key, _label in fields:
        arr = np.zeros((len(node_values), len(time_values)))
        for row in data:
            i = int(row["axial_node"])
            k = time_index[float(row["time_s"])]
            arr[i, k] = row[key]
        arrays.append(arr)

    vmin = min(float(arr.min()) for arr in arrays)
    vmax = max(float(arr.max()) for arr in arrays)
    fig, axes = plt.subplots(3, 1, figsize=(13, 9), sharex=True, sharey=True)
    for ax, arr, (_key, label) in zip(axes, arrays, fields):
        im = ax.imshow(
            arr,
            aspect="auto",
            origin="lower",
            extent=[time_values[0] / 60.0, time_values[-1] / 60.0, node_values[0], node_values[-1]],
            cmap="inferno",
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(label, fontsize=10)
        ax.set_ylabel("Axial node")
        ax.set_yticks(np.arange(0, 40, 4))
        ax.grid(False)

    axes[-1].set_xlabel("Time [min]")
    cbar = fig.colorbar(im, ax=axes, shrink=0.92, pad=0.02)
    cbar.set_label("Temperature [degC]")
    fig.suptitle("FrontISTR Round Bar: Axial Node Temperature Map", fontsize=12)
    out = RESULTS_IMG_DIR / "axial_nodes_heatmap_frontistr_da.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")


def plot_axial_error_heatmap() -> None:
    """Plot error heatmaps for model-only and DA relative to truth."""
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    data = read_csv(RESULTS_DIR / "axial_node_temperature_comparison.csv")
    time_values = np.unique(data["time_s"])
    node_values = np.unique(data["axial_node"]).astype(int)
    time_index = {float(t): idx for idx, t in enumerate(time_values)}

    errors = {}
    for key in ["without_da", "with_da"]:
        arr = np.zeros((len(node_values), len(time_values)))
        for row in data:
            i = int(row["axial_node"])
            k = time_index[float(row["time_s"])]
            arr[i, k] = row[key] - row["truth"]
        errors[key] = arr

    vmax = max(float(np.max(np.abs(arr))) for arr in errors.values())
    fig, axes = plt.subplots(2, 1, figsize=(13, 6.5), sharex=True, sharey=True)
    for ax, key, label in [
        (axes[0], "without_da", "without DA - truth"),
        (axes[1], "with_da", "with DA - truth"),
    ]:
        im = ax.imshow(
            errors[key],
            aspect="auto",
            origin="lower",
            extent=[time_values[0] / 60.0, time_values[-1] / 60.0, node_values[0], node_values[-1]],
            cmap="coolwarm",
            vmin=-vmax,
            vmax=vmax,
        )
        ax.set_title(label, fontsize=10)
        ax.set_ylabel("Axial node")
        ax.set_yticks(np.arange(0, 40, 4))

    axes[-1].set_xlabel("Time [min]")
    cbar = fig.colorbar(im, ax=axes, shrink=0.9, pad=0.02)
    cbar.set_label("Error [degC]")
    fig.suptitle("FrontISTR Round Bar: Axial Node Error Map", fontsize=12)
    out = RESULTS_IMG_DIR / "axial_nodes_error_heatmap_frontistr_da.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")


if __name__ == "__main__":
    plot_measurement_points()
    plot_all_axial_nodes_grid()
    plot_axial_heatmap()
    plot_axial_error_heatmap()
