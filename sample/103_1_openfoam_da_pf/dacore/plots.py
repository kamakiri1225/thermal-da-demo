"""データ同化結果の可視化(温度時刻歴の補正・RMSE・パラメータ推定)."""

from __future__ import annotations

import os

import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
import numpy as np

K = 273.15


def _setup_japanese_font():
    """日本語ラベル用フォントを登録する(見つからなければ英字フォントのまま)."""
    for name in ("Noto Sans CJK JP", "IPAGothic", "TakaoGothic", "VL Gothic"):
        try:
            path = fm.findfont(name, fallback_to_default=False)
            fm.fontManager.addfont(path)
            plt.rcParams["font.family"] = name
            plt.rcParams["axes.unicode_minus"] = False
            return name
        except Exception:
            pass
    # フォント名で見つからない場合はファイルを直接探して登録
    patterns = [
        os.path.expanduser("~/.fonts/*CJK*jp*.otf"),
        os.path.expanduser("~/.fonts/*CJK*JP*.otf"),
        "/usr/share/fonts/**/NotoSansCJK*.otf",
        "/usr/share/fonts/**/DroidSansFallback*.ttf",
    ]
    for pat in patterns:
        for path in glob.glob(pat, recursive=True):
            try:
                fm.fontManager.addfont(path)
                prop = fm.FontProperties(fname=path)
                plt.rcParams["font.family"] = prop.get_name()
                plt.rcParams["axes.unicode_minus"] = False
                return prop.get_name()
            except Exception:
                continue
    return None


_FONT = _setup_japanese_font()


def plot_state_convergence(hist, title, out_path):
    """5ノードそれぞれについて 真値 vs 解析平均±1σ + 観測 を描く(補正の時刻歴)."""
    names = hist["node_names"]
    t = hist["times"]
    obs_nodes = set(hist["obs_nodes"])
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    axes = axes.ravel()
    for i, name in enumerate(names):
        ax = axes[i]
        truth = hist["truth_T"][:, i] - K
        mean = hist["ana_mean_T"][:, i] - K
        std = hist["ana_std_T"][:, i]
        ax.fill_between(t, mean - std, mean + std, color="tab:blue", alpha=0.2,
                        label="解析 ±1σ")
        ax.plot(t, mean, color="tab:blue", lw=1.8, label="解析平均")
        ax.plot(t, truth, color="k", lw=1.8, ls="--", label="真値")
        if i in obs_nodes:
            oi = hist["obs_nodes"].index(i)
            ax.scatter(hist["obs_times"],
                       np.array(hist["obs_values_C"])[:, oi],
                       s=12, color="tab:red", alpha=0.6, zorder=5, label="観測")
            tag = "観測ノード"
        else:
            tag = "未観測ノード"
        ax.set_title(f"{name}  ({tag})")
        ax.grid(alpha=0.3)
        if i == 0:
            ax.legend(fontsize=8, loc="upper right")
    axes[-1].axis("off")
    for ax in axes[3:]:
        ax.set_xlabel("time [s]")
    for ax in (axes[0], axes[3]):
        ax.set_ylabel("temperature [degC]")
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_rmse(hist, title, out_path):
    t = hist["times"]
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(t, hist["rmse_all"], "-o", ms=3, label="全ノード")
    ax.plot(t, hist["rmse_obs"], "-s", ms=3, label="観測ノード")
    ax.plot(t, hist["rmse_unobs"], "-^", ms=3, label="未観測ノード")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("温度 RMSE [K]")
    ax.set_yscale("log")
    ax.set_title(title)
    ax.grid(alpha=0.3, which="both")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def plot_params(hist, title, out_path):
    t = hist["times"]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, key, tru, lab in [
        (axes[0], "q", hist["q_scale_true"], "q_scale (ヒータ倍率)"),
        (axes[1], "h", hist["h_true"], "h (放熱係数) [W/K]"),
    ]:
        m = hist[f"{key}_mean"]
        s = hist[f"{key}_std"]
        ax.fill_between(t, m - s, m + s, alpha=0.2, color="tab:green")
        ax.plot(t, m, color="tab:green", lw=1.8, label="推定平均 ±1σ")
        ax.axhline(tru, color="k", ls="--", lw=1.5, label="真値")
        ax.set_xlabel("time [s]")
        ax.set_title(lab)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=9)
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def save_history_csv(hist, out_path):
    """サイクルごとの要約(時刻・RMSE・パラメータ推定・各ノード真値/解析)を CSV 保存."""
    names = hist["node_names"]
    cols = ["time_s", "rmse_all_K", "rmse_obs_K", "rmse_unobs_K",
            "q_mean", "q_std", "h_mean", "h_std", "ess"]
    header = cols[:]
    for n in names:
        header += [f"truth_{n}_C", f"ana_{n}_C", f"anastd_{n}_K"]
    rows = []
    for j in range(len(hist["times"])):
        row = [hist["times"][j], hist["rmse_all"][j], hist["rmse_obs"][j],
               hist["rmse_unobs"][j], hist["q_mean"][j], hist["q_std"][j],
               hist["h_mean"][j], hist["h_std"][j], hist["ess"][j]]
        for i in range(len(names)):
            row += [hist["truth_T"][j, i] - K, hist["ana_mean_T"][j, i] - K,
                    hist["ana_std_T"][j, i]]
        rows.append(row)
    arr = np.array(rows)
    np.savetxt(out_path, arr, delimiter=",", header=",".join(header),
               comments="", fmt="%.6g")
