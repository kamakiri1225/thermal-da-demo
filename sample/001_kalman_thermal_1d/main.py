"""
main.py
=======
カルマンフィルタによる1次元熱モデルのデータ同化デモ。

両端2点だけを温度センサで観測し、未観測の内部節点温度と
棒全体の熱変位を推定する。

注意:
    グラフ内のタイトル、凡例、軸ラベルは文字化け回避のため英語表記にする。
    プログラム説明のコメントと docstring は日本語で記述する。
"""

from __future__ import annotations

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from kalman_filter import KalmanFilter
from thermal_model import ThermalBarModel


# シミュレーション設定
N_NODES = 10
DT = 10.0
N_STEPS = 540
SENSOR_NODES = [0, N_NODES - 1]

# 物理パラメータ
ALPHA_DIFF = 1.2e-5
DX = 0.10
H_CONV = 8.0
RHO_CP = 3.5e6
ALPHA_EXP = 12.0e-6
T_INIT = 20.0

# 熱入力条件（真値とモデルで同一 — 熱入力の誤差なし）
Q_ON = 6.5e3        # 発熱量 [W/m]（真値・モデル共通）
CYCLE_ON  = 60
CYCLE_OFF = 30

# シナリオ設定（誤差源: 初期条件誤差のみ）
# ──────────────────────────────────────────────────────────────────
# 真値は T_INIT（20°C）から、No-DA/KF は T_INIT_EST（22°C）から開始する。
# 「システム起動時に初期温度を正確に知ることができない」状況を模擬する。
# 両者の熱入力モデルは同一（誤差源は初期条件のみ）。
# ──────────────────────────────────────────────────────────────────
T_INIT_EST = 22.0   # No-DA・KF の初期推定温度（真値より 2°C 高め）

# KF パラメータ（最適値の根拠コメント付き）
# ──────────────────────────────────────────────────────────────────
# PROC_NOISE_STD : プロセスノイズ σ_Q = 0.03°C — 現実の揺らぎに対応した値
# OBS_NOISE_STD  : 観測ノイズ σ_R = 0.20°C   — センサ精度に合わせた値（変更不要）
# Q_SCALE        : Q の拡大係数 = 1.0
#                  熱入力誤差なし → モデルは正確 → 拡大不要（最適値: 1.0）
# P0_DIAG        : 初期共分散の対角成分 = 4.0°C²
#                  初期温度誤差 2°C → σ_init = 2°C → P0 = σ² = 4.0（最適値）
# L_P0           : 初期共分散の空間相関長 = 4 節点
#                  初期誤差が全節点で均一（+2°C）なので相関長を長めに設定する。
#                  これにより Node 0/9 センサの補正が隣接節点へ波及できる（最適値: 4〜5）。
# ──────────────────────────────────────────────────────────────────
PROC_NOISE_STD = 0.03
OBS_NOISE_STD  = 0.20
Q_SCALE  = 1.0    # Q 拡大係数（モデル誤差なし → 1.0 が最適）
P0_DIAG  = 4.0    # 初期共分散対角成分 [°C²]（= (2°C)²、最適値）
L_P0     = 4.0    # 初期共分散の空間相関長 [節点数]（最適値）


def make_heat_input(n_steps: int) -> np.ndarray:
    """スピンドルの ON/OFF サイクルを模した熱入力系列 [W/m] を作る。"""
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([Q_ON if (k % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def compute_displacement(temps: np.ndarray, t_ref: float = T_INIT) -> np.ndarray:
    """節点温度から棒全体の熱変位 [um] を計算する。"""
    delta_T = temps - t_ref
    return ALPHA_EXP * np.sum(delta_T, axis=1) * DX * 1e6


def _calc_rmse(a: np.ndarray, b: np.ndarray) -> float:
    """RMSE を計算する。"""
    return float(np.sqrt(np.mean((a - b) ** 2)))


def run() -> dict:
    """真値生成、観測生成、カルマンフィルタ推定を一通り実行する。"""
    model = ThermalBarModel(
        n_nodes=N_NODES,
        dt=DT,
        alpha=ALPHA_DIFF,
        dx=DX,
        h_conv=H_CONV,
        rho_cp=RHO_CP,
    )

    # 観測行列 H は、センサが置かれた節点だけを取り出す行列。
    n_obs = len(SENSOR_NODES)
    H = np.zeros((n_obs, N_NODES))
    for i, node in enumerate(SENSOR_NODES):
        H[i, node] = 1.0

    Q_mat = np.eye(N_NODES) * PROC_NOISE_STD**2 * Q_SCALE
    R_mat = np.eye(n_obs) * OBS_NOISE_STD**2

    # 真値の初期温度（正確な値）と推定側の初期温度（2°C 過大評価）を分ける。
    x0_true  = np.ones(N_NODES) * T_INIT      # True: 正しい初期状態
    x0_model = np.ones(N_NODES) * T_INIT_EST  # No-DA・KF: 2°C 高めにスタート

    # 初期共分散 P0: 空間相関付き
    # 全節点が同じ方向（+2°C）にずれているので、空間相関を持たせると
    # センサ節点の補正が隣接節点へ波及しやすくなる。
    P0_mat = np.array([
        [P0_DIAG * np.exp(-abs(i - j) / L_P0) for j in range(N_NODES)]
        for i in range(N_NODES)
    ], dtype=float)

    kf = KalmanFilter(A=model.A_d, B=model.B_d, H=H, Q=Q_mat, R=R_mat,
                      x0=x0_model.copy(), P0=P0_mat)

    # 記録用バッファ
    true_T = np.zeros((N_STEPS, N_NODES))
    model_T = np.zeros((N_STEPS, N_NODES))   # DA なし（モデルのみ）
    kf_T = np.zeros((N_STEPS, N_NODES))
    obs_y = np.zeros((N_STEPS, n_obs))
    kf_sigma = np.zeros((N_STEPS, N_NODES))
    kf_P_trace = np.zeros(N_STEPS)

    u_series = make_heat_input(N_STEPS)
    # No-DA と KF は同一の熱入力・同一の初期条件を使う。
    # 唯一の違いはセンサ観測を使うかどうかのみ。
    x_true  = x0_true.copy()
    x_model = x0_model.copy()

    for k in range(N_STEPS):
        u = np.array([u_series[k]])

        # 真値はプロセスノイズ込みで進める。
        x_true = model.step(x_true, u, noise_std=PROC_NOISE_STD)
        true_T[k] = x_true

        # DA なし: 観測を使わずモデルだけで予測する（同じ熱入力、センサなし）。
        x_model = model.step(x_model, u, noise_std=0.0)
        model_T[k] = x_model

        # センサ値は真値に観測ノイズを重ねて作る。
        y = H @ x_true + np.random.randn(n_obs) * OBS_NOISE_STD
        obs_y[k] = y

        # カルマンフィルタで全節点の温度を推定する（No-DA と全く同じモデルを使用）。
        x_est, P_est = kf.step(u, y)
        kf_T[k] = x_est
        kf_sigma[k] = np.sqrt(np.diag(P_est))
        kf_P_trace[k] = np.trace(P_est)

    return {
        "time": np.arange(N_STEPS) * DT,
        "u_series": u_series,
        "true_T": true_T,
        "model_T": model_T,
        "kf_T": kf_T,
        "obs_y": obs_y,
        "kf_sigma": kf_sigma,
        "kf_P_trace": kf_P_trace,
        "disp_true": compute_displacement(true_T),
        "disp_model": compute_displacement(model_T),
        "disp_kf": compute_displacement(kf_T),
    }


def print_summary(res: dict) -> None:
    """定常域の推定精度を標準出力にまとめる。"""
    start = 2 * N_STEPS // 3
    err = res["true_T"][start:] - res["kf_T"][start:]

    print()
    print("=" * 56)
    print("  カルマンフィルタ データ同化 結果サマリ")
    print("=" * 56)
    print(f"  評価区間: {res['time'][start] / 60:.1f} - {res['time'][-1] / 60:.1f} min")
    print()
    print("  温度推定 RMSE（定常域）")
    for i in range(N_NODES):
        rmse = np.sqrt(np.mean(err[:, i] ** 2))
        tag = "センサ節点" if i in SENSOR_NODES else "未観測節点"
        print(f"    節点 {i} ({tag}) : {rmse:.4f} degC")

    d_err = res["disp_true"][start:] - res["disp_kf"][start:]
    print()
    print(f"  熱変位推定 RMSE: {np.sqrt(np.mean(d_err**2)):.4f} um")
    print("=" * 56)
    print()


def plot_results(res: dict, save_path: str = "results.png") -> None:
    """推定結果を英語表記のグラフとして保存する。"""
    time_min = res["time"] / 60.0
    hidden_nodes = [i for i in range(N_NODES) if i not in SENSOR_NODES]
    cmap = plt.cm.tab10(np.linspace(0, 1, N_NODES))
    t_start = 2 * N_STEPS // 3

    rmse_T = [
        _calc_rmse(res["true_T"][t_start:, i], res["kf_T"][t_start:, i])
        for i in range(N_NODES)
    ]
    rmse_disp = _calc_rmse(res["disp_true"][t_start:], res["disp_kf"][t_start:])
    peak_disp = float(np.max(np.abs(res["disp_true"])))

    fig = plt.figure(figsize=(15, 16))
    fig.suptitle(
        "Kalman Filter Data Assimilation Demo (Case 1: Initial Condition Error)\n"
        f"True starts at {T_INIT:.0f}°C  |  No-DA & KF start at {T_INIT_EST:.0f}°C (2°C offset)  |  "
        f"Same heat model (Q={Q_ON:.0f} W/m)  |  Sensors: Node 0 & {N_NODES-1}",
        fontsize=10,
        y=0.995,
    )
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.35)

    ax1 = fig.add_subplot(gs[0, 0])
    for i in range(N_NODES):
        label = f"Node {i}" + (" (sensor)" if i in SENSOR_NODES else " (hidden)")
        ax1.plot(time_min, res["true_T"][:, i], color=cmap[i], lw=1.5, label=label)
    for s in SENSOR_NODES:
        ax1.plot(time_min, res["obs_y"][:, SENSOR_NODES.index(s)], color=cmap[s], alpha=0.35, lw=0.8, ls=":")
    ax1.set_title("(1) True temperature field (dotted: sensor measurement with noise)", fontsize=10)
    ax1.set_xlabel("Time [min]")
    ax1.set_ylabel("Temperature [degC]")
    ax1.legend(fontsize=7.5, loc="upper left")
    ax1.grid(True, alpha=0.3)

    ax2 = fig.add_subplot(gs[0, 1])
    for i in hidden_nodes:
        sigma = res["kf_sigma"][:, i]
        ax2.plot(time_min, res["true_T"][:, i], color=cmap[i], ls="--", lw=1.5, label=f"Node {i} (true)")
        ax2.plot(time_min, res["kf_T"][:, i], color=cmap[i], ls="-", lw=1.8, label=f"Node {i} (KF est.)")
        ax2.fill_between(time_min, res["kf_T"][:, i] - 2 * sigma, res["kf_T"][:, i] + 2 * sigma, color=cmap[i], alpha=0.15)
    rmse_txt = "\n".join([f"Node {i} RMSE={rmse_T[i]:.3f}degC" for i in hidden_nodes])
    ax2.text(0.97, 0.04, rmse_txt, transform=ax2.transAxes, fontsize=7.5, va="bottom", ha="right",
             bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.8))
    ax2.set_title("(2) Hidden node temperature estimation\n(dashed: true  solid: KF est.  shaded: +/-2sigma)", fontsize=9)
    ax2.set_xlabel("Time [min]")
    ax2.set_ylabel("Temperature [degC]")
    ax2.legend(fontsize=6.5, ncol=2)
    ax2.grid(True, alpha=0.3)

    ax3 = fig.add_subplot(gs[1, 0])
    for i in range(N_NODES):
        err = res["true_T"][:, i] - res["kf_T"][:, i]
        tag = "hidden" if i not in SENSOR_NODES else "sensor"
        ax3.plot(time_min, err, color=cmap[i], lw=2.0 if tag == "hidden" else 1.0,
                 ls="-" if tag == "hidden" else ":", label=f"Node {i} ({tag}) RMSE={rmse_T[i]:.3f}degC")
    ax3.axhline(0, color="k", lw=0.8, ls="--")
    ax3.axvline(time_min[t_start], color="red", lw=1.0, ls="--", alpha=0.6)
    ax3.set_title("(3) Estimation error (true - KF)\nRMSE in steady-state shown in legend", fontsize=9)
    ax3.set_xlabel("Time [min]")
    ax3.set_ylabel("Error [degC]")
    ax3.legend(fontsize=6.5, ncol=2)
    ax3.grid(True, alpha=0.3)

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(time_min, res["disp_true"], "k-", lw=2.0, label="True")
    ax4.plot(time_min, res["disp_kf"], "r--", lw=2.0, label="KF estimate")
    ax4.text(0.97, 0.07, f"Est. RMSE = {rmse_disp:.3f} um\n({rmse_disp / peak_disp * 100:.1f}% of peak {peak_disp:.1f} um)",
             transform=ax4.transAxes, fontsize=8, va="bottom", ha="right",
             bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.8))
    ax4.set_title("(4) Thermal displacement (total bar elongation)\nEstimated from 2 sensors only", fontsize=9)
    ax4.set_xlabel("Time [min]")
    ax4.set_ylabel("Displacement [um]")
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)

    ax5 = fig.add_subplot(gs[2, 0])
    ax5.step(time_min, res["u_series"] / 1e3, "steelblue", where="post", lw=1.5)
    ax5.fill_between(time_min, 0, res["u_series"] / 1e3, step="post", alpha=0.15, color="steelblue")
    ax5.set_title("(5) Heat input pattern\n(machine tool spindle ON/OFF cycle)", fontsize=9)
    ax5.set_xlabel("Time [min]")
    ax5.set_ylabel("Heat input [kW/m]")
    ax5.set_ylim(-1, Q_ON / 1e3 * 1.3)
    ax5.grid(True, alpha=0.3)

    ax6 = fig.add_subplot(gs[2, 1])
    snap_steps = [N_STEPS // 6, N_STEPS // 3, N_STEPS // 2, N_STEPS - 1]
    snap_cols = ["royalblue", "forestgreen", "darkorange", "crimson"]
    x_pos = np.arange(N_NODES) * DX * 100
    for t_idx, color in zip(snap_steps, snap_cols):
        label = f"t={time_min[t_idx]:.0f}min"
        ax6.plot(x_pos, res["true_T"][t_idx], color=color, ls="--", marker="o", ms=5, lw=1.5, label=f"{label} true")
        ax6.plot(x_pos, res["kf_T"][t_idx], color=color, ls="-", marker="s", ms=5, lw=1.5, label=f"{label} KF")
    ax6.set_title("(6) Temperature snapshots (spatial distribution)\ndashed: true  solid: KF", fontsize=9)
    ax6.set_xlabel("Position [cm]")
    ax6.set_ylabel("Temperature [degC]")
    ax6.legend(fontsize=7, ncol=2, loc="upper right")
    ax6.grid(True, alpha=0.3)

    ax7 = fig.add_subplot(gs[3, 0])
    sigma_mean = np.mean(res["kf_sigma"], axis=1)
    ax7.plot(time_min, sigma_mean, "purple", lw=1.8, label="Mean sigma")
    for i in range(N_NODES):
        ax7.plot(time_min, res["kf_sigma"][:, i], color=cmap[i], lw=0.8, alpha=0.5, label=f"Node {i} sigma")
    ax7.axvline(time_min[t_start], color="red", lw=1.0, ls="--", alpha=0.6)
    ax7.set_title("(7) Uncertainty convergence\n(KF learns over time: sigma decreases)", fontsize=9)
    ax7.set_xlabel("Time [min]")
    ax7.set_ylabel("Estimated std dev sigma [degC]")
    ax7.legend(fontsize=6.5, ncol=4)
    ax7.grid(True, alpha=0.3)
    ax7.set_ylim(bottom=0)

    ax8 = fig.add_subplot(gs[3, 1])
    ax8.axis("off")
    summary = (
        "[Data Assimilation Summary]\n\n"
        f"Sensors   : {len(SENSOR_NODES)} / {N_NODES} nodes\n"
        f"Sensor noise: {OBS_NOISE_STD:.2f} degC\n\n"
        "Temperature RMSE (steady-state):\n"
        + "".join([f"  Node {i}: {rmse_T[i]:.4f} degC\n" for i in range(N_NODES)])
        + "\nDisplacement estimation:\n"
        f"  Peak disp. : {peak_disp:.2f} um\n"
        f"  Est. RMSE  : {rmse_disp:.4f} um\n"
    )
    ax8.text(0.03, 0.97, summary, transform=ax8.transAxes, fontsize=8.5, va="top", ha="left",
             family="monospace", bbox=dict(boxstyle="round,pad=0.5", fc="lightyellow", ec="goldenrod", lw=1.5))

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


def plot_all_nodes(res: dict, save_path: str = "results_all_nodes.png") -> None:
    """全 10 節点の真値 vs No-DA vs KF を 2 列 × 5 行で並べて保存する。"""
    time_min = res["time"] / 60.0
    t_start = 2 * N_STEPS // 3

    fig, axes = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
    fig.suptitle(
        "All-Node Comparison: True vs No-DA vs Kalman Filter (DA)\n"
        f"Initial error: True={T_INIT:.0f}°C, No-DA & KF={T_INIT_EST:.0f}°C (+2°C)  |  "
        f"Same heat model (Q={Q_ON:.0f} W/m)  |  Sensors: Node 0 & {N_NODES-1}\n"
        "No-DA and KF start from IDENTICAL conditions — only difference: KF uses sensor observations",
        fontsize=9, y=0.995,
    )
    plt.subplots_adjust(hspace=0.45, wspace=0.3)

    for node in range(N_NODES):
        row, col = divmod(node, 2)
        ax = axes[row, col]
        is_sensor = node in SENSOR_NODES

        ax.plot(time_min, res["true_T"][:, node],
                color="black", ls="-", lw=2.0, label="True")
        ax.plot(time_min, res["model_T"][:, node],
                color="steelblue", ls="--", lw=1.5, label="Model only (no DA)")
        ax.plot(time_min, res["kf_T"][:, node],
                color="crimson", ls="-", lw=1.8, label="Kalman filter (DA)")

        if not is_sensor:
            ax.fill_between(
                time_min,
                res["kf_T"][:, node] - 2 * res["kf_sigma"][:, node],
                res["kf_T"][:, node] + 2 * res["kf_sigma"][:, node],
                color="crimson", alpha=0.12, label="KF ±2σ",
            )
        else:
            obs_idx = SENSOR_NODES.index(node)
            ax.plot(time_min, res["obs_y"][:, obs_idx],
                    color="gray", ls=":", lw=0.9, alpha=0.55, label="Sensor obs.")

        r_m = _calc_rmse(res["true_T"][t_start:, node], res["model_T"][t_start:, node])
        r_k = _calc_rmse(res["true_T"][t_start:, node], res["kf_T"][t_start:, node])
        ax.text(
            0.98, 0.04,
            f"RMSE (last 30min)\nNo DA: {r_m:.3f}°C\nDA:    {r_k:.3f}°C",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.25", fc="lightyellow", ec="gray", alpha=0.85),
        )
        ax.axvline(time_min[t_start], color="red", lw=0.8, ls=":", alpha=0.4)
        tag = "◆ sensor" if is_sensor else "hidden"
        ax.set_title(f"Node {node}  [{tag}]", fontsize=10)
        ax.set_ylabel("Temp [°C]", fontsize=8)
        ax.grid(True, alpha=0.3)
        if row == 4:
            ax.set_xlabel("Time [min]")
        if node == 0:
            ax.legend(fontsize=8, loc="upper left")

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


if __name__ == "__main__":
    np.random.seed(42)
    print("データ同化シミュレーションを開始します...")
    results = run()
    print_summary(results)
    output_path = os.path.join(os.path.dirname(__file__), "results.png")
    plot_results(results, save_path=output_path)
    output_all = os.path.join(os.path.dirname(__file__), "results_all_nodes.png")
    plot_all_nodes(results, save_path=output_all)
