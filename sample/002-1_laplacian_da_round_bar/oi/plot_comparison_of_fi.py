"""
Generate comparison figure: OpenFOAM vs FrontISTR, all axial nodes.
Sensor nodes and validation nodes are clearly distinguished.
"""

from __future__ import annotations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
import numpy as np

# ── paths ────────────────────────────────────────────────────────────────────
BASE_DIR   = Path(__file__).parent
OF_DIR     = BASE_DIR / "results"
FI_DIR     = BASE_DIR.parent.parent / "002-1-FrontISTR_round_bar_da" / "oi" / "results"
OUT_PNG    = BASE_DIR.parent / "assets" / "comparison_of_vs_fistr.png"

OF_TRUTH   = OF_DIR / "truth"      / "axial_temperature_history.csv"
OF_MODEL   = OF_DIR / "model_only" / "axial_temperature_history.csv"
OF_DA      = OF_DIR / "with_da"    / "axial_temperature_history.csv"
FI_TRUTH   = FI_DIR / "truth"      / "axial_temperature_history.csv"
FI_MODEL   = FI_DIR / "model_only" / "axial_temperature_history.csv"
FI_DA      = FI_DIR / "with_da"    / "axial_temperature_history.csv"

# ── DA design ────────────────────────────────────────────────────────────────
ASSIM_NODES      = {4, 16}
VALIDATION_NODES = {28, 36}

# ── colours ──────────────────────────────────────────────────────────────────
C_TRUTH  = "#222222"
C_MOD_OF = "#d97706"   # amber
C_MOD_FI = "#b45309"   # dark amber
C_OI_OF  = "#1d4ed8"   # blue
C_OI_FI  = "#dc2626"   # red

BG_SENSOR     = "#FFF9C4"   # yellow – assimilation sensor
BG_VALIDATION = "#DBEAFE"   # blue   – validation only
BG_NORMAL     = "white"

BORDER_SENSOR     = "#F59E0B"
BORDER_VALIDATION = "#3B82F6"


def load(path: Path):
    with path.open() as f:
        header = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1)
    return data, header


def parse_labels(header):
    out = []
    for h in header[1:]:
        # "N4_x0.034m"  or  "N4" (FrontISTR)
        if "_x" in h:
            left, xp = h.split("_x", 1)
            node = int(left[1:])
            x_m  = float(xp[:-1])
        else:
            node = int(h[1:])
            x_m  = None
        out.append((node, x_m))
    return out


def main():
    of_truth, hdr = load(OF_TRUTH);  of_model, _ = load(OF_MODEL);  of_da, _ = load(OF_DA)
    fi_truth, _   = load(FI_TRUTH);  fi_model, _ = load(FI_MODEL);  fi_da, _ = load(FI_DA)

    labels   = parse_labels(hdr)
    time_s   = of_truth[:, 0]
    n_nodes  = of_truth.shape[1] - 1

    # align FrontISTR time axis (may differ in length)
    fi_time  = fi_truth[:, 0]
    fi_idx   = [np.argmin(np.abs(fi_time - t)) for t in time_s]

    of_truth_y = of_truth[:, 1:];  of_model_y = of_model[:, 1:];  of_da_y = of_da[:, 1:]
    fi_truth_y = fi_truth[fi_idx, 1:];  fi_model_y = fi_model[fi_idx, 1:];  fi_da_y = fi_da[fi_idx, 1:]

    all_vals = np.concatenate([of_truth_y.ravel(), of_da_y.ravel(),
                                fi_truth_y.ravel(), fi_da_y.ravel()])
    ymin = float(all_vals.min());  ymax = float(all_vals.max())
    pad  = max(0.5, 0.08 * (ymax - ymin))

    NCOLS = 4
    NROWS = int(np.ceil(n_nodes / NCOLS))
    fig, axes = plt.subplots(NROWS, NCOLS,
                             figsize=(16, NROWS * 2.2),
                             sharex=True, sharey=True)
    fig.patch.set_facecolor("white")

    for idx in range(NROWS * NCOLS):
        ax = axes.flat[idx]
        ax.set_facecolor("white")
        if idx >= n_nodes:
            ax.axis("off")
            continue

        node, x_m = labels[idx]
        x_label = f"{x_m:.3f} m" if x_m is not None else ""

        # ── background & border colour by role ──────────────────────────────
        if node in ASSIM_NODES:
            bg     = BG_SENSOR
            border = BORDER_SENSOR
            role   = "[S] Sensor"
        elif node in VALIDATION_NODES:
            bg     = BG_VALIDATION
            border = BORDER_VALIDATION
            role   = "[V] Validation"
        else:
            bg     = BG_NORMAL
            border = "#cccccc"
            role   = ""

        ax.set_facecolor(bg)
        for spine in ax.spines.values():
            spine.set_edgecolor(border)
            spine.set_linewidth(2.5 if role else 0.8)

        # ── lines ───────────────────────────────────────────────────────────
        ax.plot(time_s, of_truth_y[:, idx],  color=C_TRUTH,  ls="-",  lw=1.6)
        ax.plot(time_s, fi_truth_y[:, idx],  color=C_TRUTH,  ls="--", lw=1.0, alpha=0.6)
        ax.plot(time_s, of_model_y[:, idx],  color=C_MOD_OF, ls="--", lw=1.1)
        ax.plot(time_s, fi_model_y[:, idx],  color=C_MOD_FI, ls=":",  lw=1.1)
        ax.plot(time_s, of_da_y[:, idx],     color=C_OI_OF,  ls="-",  lw=1.8)
        ax.plot(time_s, fi_da_y[:, idx],     color=C_OI_FI,  ls="-",  lw=1.4, alpha=0.85)

        # ── title ───────────────────────────────────────────────────────────
        title_color = (BORDER_SENSOR if node in ASSIM_NODES
                       else BORDER_VALIDATION if node in VALIDATION_NODES
                       else "#333333")
        title = f"N{node}  {x_label}"
        if role:
            title += f"\n{role}"
        ax.set_title(title, fontsize=13,
                     fontweight="bold" if role else "normal",
                     color=title_color, pad=3)

        ax.grid(True, alpha=0.2, color="#aaaaaa")
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.tick_params(labelsize=10)

        if idx % NCOLS == 0:
            ax.set_ylabel("T [°C]", fontsize=11)
        if idx >= (NROWS - 1) * NCOLS:
            ax.set_xlabel("Time [s]", fontsize=11)

    # ── legend ───────────────────────────────────────────────────────────────
    legend_handles = [
        Line2D([0],[0], color=C_TRUTH,  ls="-",  lw=2.2, label="Truth (OF)"),
        Line2D([0],[0], color=C_TRUTH,  ls="--", lw=1.4, alpha=0.6, label="Truth (FI)"),
        Line2D([0],[0], color=C_MOD_OF, ls="--", lw=1.8, label="Model only (OF)"),
        Line2D([0],[0], color=C_MOD_FI, ls=":",  lw=1.8, label="Model only (FI)"),
        Line2D([0],[0], color=C_OI_OF,  ls="-",  lw=2.5, label="OI DA (OF)"),
        Line2D([0],[0], color=C_OI_FI,  ls="-",  lw=2.0, label="OI DA (FI)"),
        mpatches.Patch(facecolor=BG_SENSOR,     edgecolor=BORDER_SENSOR,
                       linewidth=2, label="[S] Assimilation sensor"),
        mpatches.Patch(facecolor=BG_VALIDATION, edgecolor=BORDER_VALIDATION,
                       linewidth=2, label="[V] Validation node"),
    ]
    fig.legend(handles=legend_handles,
               loc="upper center",
               ncol=4,
               fontsize=13,
               frameon=True,
               framealpha=0.9,
               edgecolor="#cccccc",
               bbox_to_anchor=(0.5, 1.0))

    fig.suptitle("All Axial Nodes: OpenFOAM vs FrontISTR — OI Data Assimilation",
                 fontsize=16, fontweight="bold", y=1.055, color="#111111")

    plt.tight_layout(rect=(0, 0, 1, 0.97))
    OUT_PNG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[saved] {OUT_PNG}")


if __name__ == "__main__":
    main()
