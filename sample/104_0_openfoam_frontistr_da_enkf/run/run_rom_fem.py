"""ROM 版 温度＋変位データ同化(0→600s)を実行し、分かりやすい図を作る.

数秒で完走。図:
  - docs/img/rom_fem_temperature.png : 5ノード温度の時刻歴(真値/同化あり/同化なし)
  - docs/img/rom_fem_displacement.png: 上面変位2点の時刻歴(真値/同化あり/同化なし)
  - docs/img/rom_fem_rmse.png        : 温度RMSEの収束(同化あり vs なし)

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/run_rom_fem.py
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore import plots as _p  # noqa: F401  日本語フォント設定
import matplotlib.pyplot as plt

from dacore.calibrate import calibrate, load_calibrated
from dacore.cht_rom import NODE_NAMES
from dacore.displacement import calibrate_displacement, load_operator
from dacore.twin_fem import run_rom_fem_twin

IMG = os.path.join(ROOT, "docs", "img")
RESULTS = os.path.join(ROOT, "results")
K = 273.15
HEATER_OFF = 300.0


def _heater_span(ax):
    ax.axvspan(0, HEATER_OFF, color="orange", alpha=0.06)
    ax.axvline(HEATER_OFF, color="orange", ls=":", lw=1)


def plot_temperature(hist, out):
    t = hist["times"]
    obs_nodes = set(hist["tnodes"])
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5), sharex=True)
    axes = axes.ravel()
    for i, name in enumerate(NODE_NAMES):
        ax = axes[i]
        _heater_span(ax)
        ax.plot(t, hist["free_T"][:, i]-K, color="0.6", lw=1.6,
                label="同化なし(でたらめのまま)")
        ax.plot(t, hist["da_T"][:, i]-K, color="tab:blue", lw=1.8,
                label="データ同化あり")
        ax.plot(t, hist["truth_T"][:, i]-K, "k--", lw=1.8, label="真値")
        if i in obs_nodes:
            oi = hist["tnodes"].index(i)
            ax.scatter(hist["obs_times"], hist["obs_T"][:, oi]-K, s=10,
                       color="tab:red", alpha=0.5, zorder=5, label="観測")
            tag = "観測点"
        else:
            tag = "未観測点"
        ax.set_title(f"{name} ({tag})")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=14, loc="upper right")
    axes[-1].axis("off")
    for ax in axes[3:]:
        ax.set_xlabel("time [s]")
    for ax in (axes[0], axes[3]):
        ax.set_ylabel("温度 [degC]")
    fig.suptitle("縮約モデルによる固体温度の時刻歴 (0-600 s, ヒータ通電 0-300 s / 遮断 300-600 s)", fontsize=20)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_displacement(hist, out):
    t = hist["times"]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    for j, name in enumerate(hist["disp_names"]):
        ax = axes[j]
        _heater_span(ax)
        ax.plot(t, hist["free_U"][:, j]*1000, color="0.6", lw=1.6,
                label="同化なし")
        ax.plot(t, hist["da_U"][:, j]*1000, color="tab:blue", lw=1.8,
                label="データ同化あり")
        ax.plot(t, hist["truth_U"][:, j]*1000, "k--", lw=1.8, label="真値")
        ax.scatter(hist["obs_times"], hist["obs_U"][:, j]*1000, s=10,
                   color="tab:red", alpha=0.5, zorder=5, label="観測")
        ax.set_title(name)
        ax.set_xlabel("time [s]")
        ax.set_ylabel("上面変位 Uz [um]")
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=14)
    fig.suptitle("縮約モデルによる上面変位の時刻歴 (FrontISTR 校正)",
                 fontsize=20)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    plt.close(fig)


def plot_rmse(hist, out):
    t = hist["times"]
    da = np.sqrt(((hist["da_T"] - hist["truth_T"])**2).mean(axis=1))
    fr = np.sqrt(((hist["free_T"] - hist["truth_T"])**2).mean(axis=1))
    fig, ax = plt.subplots(figsize=(8.5, 4.8))
    _heater_span(ax)
    ax.plot(t, fr, "-o", ms=3, color="0.6", label="同化なし")
    ax.plot(t, da, "-o", ms=3, color="tab:blue", label="データ同化あり")
    ax.set_xlabel("time [s]"); ax.set_ylabel("温度RMSE(全5ノード) [K]")
    ax.set_yscale("log"); ax.grid(alpha=0.3, which="both"); ax.legend()
    ax.set_title("縮約モデルにおける温度 RMSE の時間変化 (データ同化あり / なし)")
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


def main():
    os.makedirs(IMG, exist_ok=True); os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(ROOT, "config", "da_config.yaml")) as f:
        cfg = yaml.safe_load(f)

    if not os.path.exists(os.path.join(ROOT, "config", "rom_calibrated.yaml")):
        calibrate()
    calib = load_calibrated()
    if not os.path.exists(os.path.join(ROOT, "config", "displacement_operator.yaml")):
        calibrate_displacement(calib)
    op = load_operator()

    hist = run_rom_fem_twin(cfg, calib, op)
    plot_temperature(hist, os.path.join(IMG, "rom_fem_temperature.png"))
    plot_displacement(hist, os.path.join(IMG, "rom_fem_displacement.png"))
    plot_rmse(hist, os.path.join(IMG, "rom_fem_rmse.png"))

    da_final = float(np.sqrt(((hist["da_T"][-1]-hist["truth_T"][-1])**2).mean()))
    fr_final = float(np.sqrt(((hist["free_T"][-1]-hist["truth_T"][-1])**2).mean()))
    print(f"[rom-fem] 最終温度RMSE: 同化あり {da_final:.3f} K / 同化なし {fr_final:.3f} K")
    print(f"[rom-fem] 図を {os.path.relpath(IMG)} に出力")


if __name__ == "__main__":
    main()
