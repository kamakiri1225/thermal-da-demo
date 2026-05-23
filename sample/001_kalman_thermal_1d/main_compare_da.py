"""
main_compare_da.py
==================
データ同化あり/なしの比較。

2つの温度センサを使ったカルマンフィルタの有無で、温度推定がどう変わるかを比較する。
グラフ内のタイトル、凡例、軸ラベルは文字化け回避のため英語表記にする。
"""

from __future__ import annotations

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from kalman_filter import KalmanFilter
from thermal_model import ThermalBarModel


# main.py と同じ条件を使う。
N_NODES = 10
DT = 10.0
N_STEPS = 540
SENSOR_NODES = [0, N_NODES - 1]

ALPHA_DIFF = 1.2e-5
DX = 0.10
H_CONV = 8.0
RHO_CP = 3.5e6
ALPHA_EXP = 12.0e-6
T_INIT = 20.0

Q_ON = 6.5e3
CYCLE_ON = 60
CYCLE_OFF = 30

PROC_NOISE_STD = 0.03
OBS_NOISE_STD = 0.20

# 比較を分かりやすくするため、推定側のモデルには熱入力の見積もり誤差を入れる。
# 実機では発熱量や境界条件を完全には分からないため、データ同化なしではこの誤差が残る。
MODEL_INPUT_SCALE = 0.90


def make_heat_input(n_steps: int) -> np.ndarray:
    """スピンドルの ON/OFF サイクルを模した熱入力系列 [W/m] を作る。"""
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([Q_ON if (k % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def compute_displacement(temps: np.ndarray) -> np.ndarray:
    """節点温度から棒全体の熱変位 [um] を計算する。"""
    return ALPHA_EXP * np.sum(temps - T_INIT, axis=1) * DX * 1e6


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    """RMSE を計算する。"""
    return float(np.sqrt(np.mean((a - b) ** 2)))


def run() -> dict:
    """データ同化なし/ありの両ケースを、同じ真値と観測値で比較する。"""
    model = ThermalBarModel(
        n_nodes=N_NODES,
        dt=DT,
        alpha=ALPHA_DIFF,
        dx=DX,
        h_conv=H_CONV,
        rho_cp=RHO_CP,
    )

    # センサ節点だけを取り出す観測行列。
    H = np.zeros((len(SENSOR_NODES), N_NODES))
    for i, node in enumerate(SENSOR_NODES):
        H[i, node] = 1.0

    # カルマンフィルタの初期化。
    Q_mat = np.eye(N_NODES) * PROC_NOISE_STD**2
    R_mat = np.eye(len(SENSOR_NODES)) * OBS_NOISE_STD**2
    x0 = np.ones(N_NODES) * T_INIT
    P0 = np.eye(N_NODES) * 2.0
    kf = KalmanFilter(A=model.A_d, B=model.B_d, H=H, Q=Q_mat, R=R_mat, x0=x0, P0=P0)

    # 記録用バッファ。
    true_T = np.zeros((N_STEPS, N_NODES))
    model_T = np.zeros((N_STEPS, N_NODES))
    kf_T = np.zeros((N_STEPS, N_NODES))
    obs_y = np.zeros((N_STEPS, len(SENSOR_NODES)))
    kf_sigma = np.zeros((N_STEPS, N_NODES))

    u_series = make_heat_input(N_STEPS)
    x_true = x0.copy()
    x_model = x0.copy()

    for k in range(N_STEPS):
        u_true = np.array([u_series[k]])
        u_model = np.array([u_series[k] * MODEL_INPUT_SCALE])

        # 真値は小さな外乱を含むものとして生成する。
        x_true = model.step(x_true, u_true, noise_std=PROC_NOISE_STD)
        true_T[k] = x_true

        # データ同化なし: センサ観測を使わず、モデルだけで予測し続ける。
        x_model = model.step(x_model, u_model, noise_std=0.0)
        model_T[k] = x_model

        # データ同化あり: 2つのセンサ観測を使ってカルマンフィルタで補正する。
        y = H @ x_true + np.random.randn(len(SENSOR_NODES)) * OBS_NOISE_STD
        obs_y[k] = y
        x_est, P_est = kf.step(u_model, y)
        kf_T[k] = x_est
        kf_sigma[k] = np.sqrt(np.diag(P_est))

    return {
        "time": np.arange(N_STEPS) * DT,
        "u_series": u_series,
        "true_T": true_T,
        "model_T": model_T,
        "kf_T": kf_T,
        "obs_y": obs_y,
        "kf_sigma": kf_sigma,
        "disp_true": compute_displacement(true_T),
        "disp_model": compute_displacement(model_T),
        "disp_kf": compute_displacement(kf_T),
    }


def print_comparison(res: dict) -> None:
    """定常域での RMSE を標準出力に表示する。"""
    t_start = 2 * N_STEPS // 3

    print()
    width = 24 + N_NODES * 10 + 14
    print("=" * width)
    print("  データ同化なし vs あり: 推定精度比較（定常域 RMSE）")
    print("=" * width)
    node_header = "".join([f"  {f'node{i}':>8}" for i in range(N_NODES)])
    print(f"  {'case':<24}{node_header}  {'disp[um]':>10}")
    print("-" * width)

    for label, key in [
        ("Model only (no DA)", "model_T"),
        ("Kalman filter (DA)", "kf_T"),
    ]:
        rmse_nodes = [rmse(res["true_T"][t_start:, i], res[key][t_start:, i]) for i in range(N_NODES)]
        disp_key = "disp_model" if key == "model_T" else "disp_kf"
        rmse_disp = rmse(res["disp_true"][t_start:], res[disp_key][t_start:])
        node_values = "".join([f"  {value:>8.4f}" for value in rmse_nodes])
        print(
            f"  {label:<24}  "
            f"{node_values}"
            f"{rmse_disp:>10.4f}"
        )

    sensor_rmse = float(
        np.mean([
            rmse(res["obs_y"][:, j], res["true_T"][:, SENSOR_NODES[j]])
            for j in range(len(SENSOR_NODES))
        ])
    )
    print("=" * width)
    print(f"  sensor noise RMSE reference: {sensor_rmse:.4f} degC")
    print()


def _plot_temperature_panel(ax, res: dict, node: int, obs_index: int | None, title: str, t_start: int) -> None:
    """1節点分の真値、モデルのみ、データ同化ありを同じ軸に描く。"""
    time_min = res["time"] / 60.0
    style_true = dict(color="black", ls="-", lw=2.2, label="True")
    style_model = dict(color="steelblue", ls="--", lw=1.8, label="Model only (no DA)")
    style_kf = dict(color="crimson", ls="-", lw=1.8, label="Kalman filter (DA)")

    ax.plot(time_min, res["true_T"][:, node], **style_true)
    if obs_index is not None:
        obs_rmse = rmse(res["obs_y"][:, obs_index], res["true_T"][:, node])
        ax.plot(
            time_min,
            res["obs_y"][:, obs_index],
            color="gray",
            ls=":",
            lw=1.0,
            alpha=0.7,
            label=f"Sensor reading (RMSE={obs_rmse:.3f}degC)",
        )
    ax.plot(time_min, res["model_T"][:, node], **style_model)
    ax.plot(time_min, res["kf_T"][:, node], **style_kf)

    # 未観測節点には推定不確かさの幅も重ねる。
    if obs_index is None:
        ax.fill_between(
            time_min,
            res["kf_T"][:, node] - 2 * res["kf_sigma"][:, node],
            res["kf_T"][:, node] + 2 * res["kf_sigma"][:, node],
            color="crimson",
            alpha=0.12,
            label="KF +/-2sigma",
        )

    r_model = rmse(res["true_T"][t_start:, node], res["model_T"][t_start:, node])
    r_kf = rmse(res["true_T"][t_start:, node], res["kf_T"][t_start:, node])
    ax.text(
        0.97,
        0.05,
        f"RMSE (steady-state)\nModel only: {r_model:.3f}degC\nDA:         {r_kf:.3f}degC",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
        ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.85),
    )
    ax.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [degC]")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)


def plot_comparison(res: dict, save_path: str = "results_compare_da.png") -> None:
    """データ同化なし/ありの差を英語表記のグラフとして保存する。"""
    time_min = res["time"] / 60.0
    t_start = 2 * N_STEPS // 3
    hidden = [i for i in range(N_NODES) if i not in SENSOR_NODES]

    fig = plt.figure(figsize=(16, 18))
    fig.suptitle(
        "Data Assimilation: with DA vs without DA\n"
        f"Two sensor nodes (Node 0 and Node {N_NODES - 1}) are used by the Kalman filter",
        fontsize=13,
        y=0.995,
    )
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.35)

    # 未観測節点の代表例。
    _plot_temperature_panel(fig.add_subplot(gs[0, 0]), res, hidden[0], None, "(1) Node 1 hidden temperature", t_start)
    _plot_temperature_panel(fig.add_subplot(gs[0, 1]), res, hidden[1], None, "(2) Node 2 hidden temperature", t_start)

    # 2つのセンサ節点を別々に表示する。
    _plot_temperature_panel(fig.add_subplot(gs[1, 0]), res, SENSOR_NODES[0], 0, "(3) Node 0 sensor temperature", t_start)
    _plot_temperature_panel(
        fig.add_subplot(gs[1, 1]),
        res,
        SENSOR_NODES[1],
        1,
        f"(4) Node {SENSOR_NODES[1]} sensor temperature",
        t_start,
    )

    # 熱変位の比較。
    ax5 = fig.add_subplot(gs[2, 0])
    r_model_d = rmse(res["disp_true"][t_start:], res["disp_model"][t_start:])
    r_kf_d = rmse(res["disp_true"][t_start:], res["disp_kf"][t_start:])
    ax5.plot(time_min, res["disp_true"], color="black", ls="-", lw=2.2, label="True")
    ax5.plot(time_min, res["disp_model"], color="steelblue", ls="--", lw=1.8, label=f"Model only RMSE={r_model_d:.3f}um")
    ax5.plot(time_min, res["disp_kf"], color="crimson", ls="-", lw=1.8, label=f"DA RMSE={r_kf_d:.3f}um")
    ax5.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)
    ax5.set_title("(5) Thermal displacement", fontsize=10)
    ax5.set_xlabel("Time [min]")
    ax5.set_ylabel("Displacement [um]")
    ax5.legend(fontsize=8)
    ax5.grid(True, alpha=0.3)

    # 未観測節点の推定誤差。
    ax6 = fig.add_subplot(gs[2, 1])
    cmap = plt.cm.Set1(np.linspace(0, 0.7, len(hidden)))
    for i, node in enumerate(hidden):
        err_model = res["true_T"][:, node] - res["model_T"][:, node]
        err_kf = res["true_T"][:, node] - res["kf_T"][:, node]
        ax6.plot(time_min, err_model, color=cmap[i], ls="--", lw=1.4, label=f"Node {node} model only")
        ax6.plot(time_min, err_kf, color=cmap[i], ls="-", lw=1.8, label=f"Node {node} DA")
    ax6.axhline(0, color="k", lw=0.8, ls="--")
    ax6.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5, label="Steady-state boundary")
    ax6.set_title("(6) Hidden-node estimation error\n(dashed: no DA, solid: with DA)", fontsize=9)
    ax6.set_xlabel("Time [min]")
    ax6.set_ylabel("Error [degC]")
    ax6.legend(fontsize=7.5, ncol=2)
    ax6.grid(True, alpha=0.3)

    # 全節点の RMSE。センサ節点も含めることで2点センサの効果を見やすくする。
    ax7 = fig.add_subplot(gs[3, :])
    nodes_label = [f"Node {i}" for i in range(N_NODES)]
    rmse_model = [rmse(res["true_T"][t_start:, i], res["model_T"][t_start:, i]) for i in range(N_NODES)]
    rmse_kf = [rmse(res["true_T"][t_start:, i], res["kf_T"][t_start:, i]) for i in range(N_NODES)]
    x = np.arange(N_NODES)
    width = 0.35
    b1 = ax7.bar(x - width / 2, rmse_model, width, label="Model only (no DA)", color="steelblue", alpha=0.8, edgecolor="k")
    b2 = ax7.bar(x + width / 2, rmse_kf, width, label="Kalman filter (DA)", color="crimson", alpha=0.8, edgecolor="k")
    ax7.bar_label(b1, fmt="%.3f", fontsize=8, padding=2)
    ax7.bar_label(b2, fmt="%.3f", fontsize=8, padding=2)
    ax7.axhline(OBS_NOISE_STD, color="gray", ls=":", lw=1.5, label=f"Sensor noise={OBS_NOISE_STD:.2f}degC")
    for node in SENSOR_NODES:
        ax7.text(node, ax7.get_ylim()[1] * 0.92, "sensor", ha="center", va="top", fontsize=8, color="dimgray")
    ax7.set_xticks(x)
    ax7.set_xticklabels(nodes_label, fontsize=10)
    ax7.set_ylabel("RMSE [degC]")
    ax7.set_title("(7) Temperature RMSE comparison for all nodes", fontsize=10)
    ax7.legend(fontsize=8)
    ax7.grid(True, alpha=0.3, axis="y")

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


if __name__ == "__main__":
    np.random.seed(42)
    print("データ同化あり/なしの比較シミュレーションを開始します...")
    results = run()
    print_comparison(results)
    output_path = os.path.join(os.path.dirname(__file__), "results_compare_da.png")
    plot_comparison(results, save_path=output_path)
