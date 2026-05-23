"""
main_da_demo.py
===============
データ同化の効果を示すデモ（現実的な条件）。

シナリオ:
  「実際の工作機械は公称値より 30% 多く発熱している。
   熱モデルは公称値しか知らない（モデル誤差）。
   さらにシステム起動時は温度が完全にわかっていない（初期条件誤差）。」

設定:
  真値  : 50 分暖機後（公称値 × 1.5 の高発熱）から, 60 分本番計算
           ＋ 本番では公称値 × 1.3 で発熱し続ける
  No-DA : 全節点 20°C から開始, 公称発熱（モデル誤差 30%）で計算のみ
  KF    : 全節点 20°C から開始, 同じ誤ったモデル + センサ観測（N0, N9）で補正

効果:
  ・Node 0/9（センサ節点）: センサが直接補正 → 劇的改善
  ・Node 1-4（N0 寄り）  : 空間相関により補正が波及
  ・Node 5-8（N9 寄り）  : 同様に N9 センサから補正波及
  ・全体的に No-DA より大幅に誤差が小さい
"""

from __future__ import annotations

import os

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from kalman_filter import KalmanFilter
from thermal_model import ThermalBarModel


N_NODES = 10
DT = 10.0
N_WARMUP = 300         # 暖機ステップ数（真値のみ先行計算）
N_STEPS  = 540         # 本番シミュレーション（90 分、main.py と同一）

SENSOR_NODES = [0, N_NODES - 1]

ALPHA_DIFF = 1.2e-5
DX = 0.10
H_CONV = 8.0
RHO_CP = 3.5e6
ALPHA_EXP = 12.0e-6
T_AMBIENT = 20.0

Q_NOMINAL = 6.5e3        # W/m ─ モデルが知っている公称発熱量
Q_TRUE_WARMUP = Q_NOMINAL * 1.5   # 暖機時の真値（高め）
Q_TRUE_MAIN   = Q_NOMINAL * 1.3   # 本番時の真値（公称より 30% 多い）
# モデルは Q_NOMINAL しか知らない → 30% 過小評価

CYCLE_ON  = 60
CYCLE_OFF = 30

# main.py / main_compare_da.py と同一の物理パラメータ
PROC_NOISE_STD = 0.03
OBS_NOISE_STD  = 0.20

# 空間相関長 [節点数]: 初期共分散の広がり
L_P0 = 4.0

# Q のスケール倍率: モデル誤差が既知のとき Q を大きくすると
# 中央節点へのセンサ補正の伝播が改善される（観測をより信頼する設定）
Q_SCALE = 5.0


def make_heat_input(n_steps: int, offset: int = 0) -> np.ndarray:
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([1.0 if ((k + offset) % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def compute_displacement(temps: np.ndarray) -> np.ndarray:
    return ALPHA_EXP * np.sum(temps - T_AMBIENT, axis=1) * DX * 1e6


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def run() -> dict:
    model = ThermalBarModel(
        n_nodes=N_NODES, dt=DT, alpha=ALPHA_DIFF,
        dx=DX, h_conv=H_CONV, rho_cp=RHO_CP,
    )
    H = np.zeros((len(SENSOR_NODES), N_NODES))
    for i, node in enumerate(SENSOR_NODES):
        H[i, node] = 1.0

    # ── 暖機フェーズ: 真値だけ高発熱で先行計算 ─────────────────────────────
    sig_warmup = make_heat_input(N_WARMUP, offset=0)
    x_true = np.ones(N_NODES) * T_AMBIENT
    for k in range(N_WARMUP):
        x_true = model.step(x_true, np.array([Q_TRUE_WARMUP * sig_warmup[k]]),
                            noise_std=PROC_NOISE_STD * 0.3)
    x0_true = x_true.copy()
    x0_est  = np.ones(N_NODES) * T_AMBIENT  # 推定側は周囲温度から開始

    delta0 = x0_true - x0_est  # 初期誤差ベクトル
    sigma0 = float(np.mean(np.abs(delta0)))   # 代表初期誤差

    print(f"[暖機後 真値]  Node0={x0_true[0]:.1f}°C  Node5={x0_true[4]:.1f}°C  Node9={x0_true[9]:.1f}°C")
    print(f"[初期誤差]   Node0={delta0[0]:.1f}°C  Node5={delta0[4]:.1f}°C  Node9={delta0[9]:.1f}°C")

    # ── KF 初期化: 空間相関を持つ P0 + 大きめの Q ────────────────────────────
    # 暖機で全節点が昇温している（正相関あり）→ 空間相関 P0 で補正を波及させる
    # Q を大きく設定: モデル誤差(30%)が既知なので観測をより信頼させる
    #   → Kalman ゲインが上がり、センサ補正が中央節点まで届きやすくなる
    Q_mat = np.eye(N_NODES) * PROC_NOISE_STD**2 * Q_SCALE
    R_mat = np.eye(len(SENSOR_NODES)) * OBS_NOISE_STD**2
    P0 = np.array([
        [float(delta0[i] * delta0[j]) * np.exp(-abs(i - j) / L_P0)
         for j in range(N_NODES)]
        for i in range(N_NODES)
    ])
    # 対角成分が小さすぎる節点（Node9 周辺は初期誤差 ≈ 0）は最小値を確保
    for i in range(N_NODES):
        P0[i, i] = max(P0[i, i], OBS_NOISE_STD**2 * 4)

    kf = KalmanFilter(A=model.A_d, B=model.B_d, H=H, Q=Q_mat, R=R_mat,
                      x0=x0_est.copy(), P0=P0)

    true_T   = np.zeros((N_STEPS, N_NODES))
    model_T  = np.zeros((N_STEPS, N_NODES))
    kf_T     = np.zeros((N_STEPS, N_NODES))
    obs_y    = np.zeros((N_STEPS, len(SENSOR_NODES)))
    kf_sigma = np.zeros((N_STEPS, N_NODES))

    sig_main = make_heat_input(N_STEPS, offset=N_WARMUP)
    x_true_k  = x0_true.copy()
    x_model_k = x0_est.copy()

    for k in range(N_STEPS):
        sig = sig_main[k]
        u_true  = np.array([Q_TRUE_MAIN * sig])    # 真値: 30% 多い発熱
        u_model = np.array([Q_NOMINAL   * sig])    # モデル: 公称発熱

        x_true_k  = model.step(x_true_k,  u_true,  noise_std=PROC_NOISE_STD * 0.3)
        x_model_k = model.step(x_model_k, u_model, noise_std=0.0)
        true_T[k]  = x_true_k
        model_T[k] = x_model_k

        y = H @ x_true_k + np.random.randn(len(SENSOR_NODES)) * OBS_NOISE_STD
        obs_y[k] = y
        x_est, P_est = kf.step(u_model, y)
        kf_T[k]     = x_est
        kf_sigma[k] = np.sqrt(np.diag(P_est))

    return {
        "time": np.arange(N_STEPS) * DT,
        "true_T": true_T, "model_T": model_T, "kf_T": kf_T,
        "obs_y": obs_y, "kf_sigma": kf_sigma,
        "x0_true": x0_true, "x0_est": x0_est,
        "disp_true":  compute_displacement(true_T),
        "disp_model": compute_displacement(model_T),
        "disp_kf":    compute_displacement(kf_T),
    }


def _rmse_box(ax, res: dict, node: int, t_start: int, time_min: np.ndarray) -> None:
    r_m = rmse(res["true_T"][t_start:, node], res["model_T"][t_start:, node])
    r_k = rmse(res["true_T"][t_start:, node], res["kf_T"][t_start:, node])
    imp = (r_m - r_k) / r_m * 100 if r_m > 0 else 0.0
    ax.text(
        0.97, 0.05,
        f"RMSE (後半 20min)\nNo DA: {r_m:.1f}°C\nDA:    {r_k:.1f}°C  ({imp:.0f}% better)",
        transform=ax.transAxes, fontsize=8, va="bottom", ha="right",
        bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.9),
    )
    ax.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)


def plot_results(res: dict, save_path: str = "results_da_demo.png") -> None:
    time_min = res["time"] / 60.0
    t_eval = N_STEPS * 2 // 3   # 後半 1/3 で RMSE 評価

    style_t = dict(color="black",     ls="-",  lw=2.2, label="True (30% excess heat)")
    style_m = dict(color="steelblue", ls="--", lw=1.8, label="Model only (nominal heat, no DA)")
    style_k = dict(color="crimson",   ls="-",  lw=1.8, label="Kalman filter (DA)")

    fig = plt.figure(figsize=(16, 20))
    fig.suptitle(
        "Data Assimilation Demo: Model Underestimates Heat by 30%\n"
        f"True Q = {Q_TRUE_MAIN:.0f} W/m,  Model Q = {Q_NOMINAL:.0f} W/m  "
        f"(30% error).  "
        f"True state pre-warmed {N_WARMUP*DT/60:.0f} min; "
        f"estimator starts at {T_AMBIENT:.0f}°C.\n"
        f"Sensors: Node 0 & Node {N_NODES-1}",
        fontsize=11, y=0.998,
    )
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.35)

    for col, (node, obs_idx, title_sfx) in enumerate([
        (0,          0, " [sensor, heat input end]"),
        (N_NODES-1,  1, " [sensor, cool end]"),
    ]):
        ax = fig.add_subplot(gs[0, col])
        ax.plot(time_min, res["true_T"][:, node], **style_t)
        ax.plot(time_min, res["obs_y"][:, obs_idx], color="gray", ls=":", lw=1.0, alpha=0.7, label="Sensor")
        ax.plot(time_min, res["model_T"][:, node], **style_m)
        ax.plot(time_min, res["kf_T"][:, node], **style_k)
        _rmse_box(ax, res, node, t_eval, time_min)
        ax.set_title(f"({'1' if col==0 else '2'}) Node {node}{title_sfx}", fontsize=10)
        ax.set_xlabel("Time [min]"); ax.set_ylabel("Temperature [°C]")
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)

    for col, node in enumerate([3, 6]):
        ax = fig.add_subplot(gs[1, col])
        ax.plot(time_min, res["true_T"][:, node], **style_t)
        ax.plot(time_min, res["model_T"][:, node], **style_m)
        ax.plot(time_min, res["kf_T"][:, node], **style_k)
        ax.fill_between(
            time_min,
            res["kf_T"][:, node] - 2 * res["kf_sigma"][:, node],
            res["kf_T"][:, node] + 2 * res["kf_sigma"][:, node],
            color="crimson", alpha=0.12, label="KF ±2σ",
        )
        _rmse_box(ax, res, node, t_eval, time_min)
        ax.set_title(f"({'3' if col==0 else '4'}) Node {node}  [hidden]", fontsize=10)
        ax.set_xlabel("Time [min]"); ax.set_ylabel("Temperature [°C]")
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)

    # ── 熱変位 ──────────────────────────────────────────────────────────────
    ax5 = fig.add_subplot(gs[2, 0])
    r_m = rmse(res["disp_true"][t_eval:], res["disp_model"][t_eval:])
    r_k = rmse(res["disp_true"][t_eval:], res["disp_kf"][t_eval:])
    ax5.plot(time_min, res["disp_true"],  **style_t)
    ax5.plot(time_min, res["disp_model"], **dict(style_m, label=f"No DA  RMSE={r_m:.1f}μm"))
    ax5.plot(time_min, res["disp_kf"],    **dict(style_k, label=f"DA     RMSE={r_k:.1f}μm"))
    ax5.axvline(time_min[t_eval], color="red", lw=1.0, ls=":", alpha=0.5)
    ax5.set_title("(5) Thermal displacement", fontsize=10)
    ax5.set_xlabel("Time [min]"); ax5.set_ylabel("Displacement [μm]")
    ax5.legend(fontsize=8); ax5.grid(True, alpha=0.3)

    # ── 全節点 RMSE 棒グラフ（後半）──────────────────────────────────────────
    ax6 = fig.add_subplot(gs[2, 1])
    rmse_m = [rmse(res["true_T"][t_eval:, i], res["model_T"][t_eval:, i]) for i in range(N_NODES)]
    rmse_k = [rmse(res["true_T"][t_eval:, i], res["kf_T"][t_eval:, i])    for i in range(N_NODES)]
    x = np.arange(N_NODES); w = 0.35
    b1 = ax6.bar(x - w/2, rmse_m, w, label="No DA", color="steelblue", alpha=0.85, edgecolor="k")
    b2 = ax6.bar(x + w/2, rmse_k, w, label="DA (KF)", color="crimson",  alpha=0.85, edgecolor="k")
    ax6.bar_label(b1, fmt="%.1f", fontsize=8, padding=2)
    ax6.bar_label(b2, fmt="%.1f", fontsize=8, padding=2)
    ax6.axhline(OBS_NOISE_STD, color="gray", ls=":", lw=1.5, label=f"Sensor noise σ={OBS_NOISE_STD:.2f}°C")
    ymax = max(rmse_m)
    for node in SENSOR_NODES:
        ax6.text(node, ymax * 1.1, "sensor", ha="center", va="bottom", fontsize=8, color="dimgray")
    ax6.set_xticks(x)
    ax6.set_xticklabels([f"N{i}" for i in range(N_NODES)], fontsize=9)
    ax6.set_ylabel("RMSE [°C]")
    ax6.set_title("(6) RMSE (last 20 min, quasi-steady state)", fontsize=10)
    ax6.legend(fontsize=8); ax6.grid(True, alpha=0.3, axis="y")

    # ── 全節点 RMS 誤差の時系列 ──────────────────────────────────────────────
    ax7 = fig.add_subplot(gs[3, :])
    all_err_m = np.sqrt(np.mean((res["model_T"] - res["true_T"])**2, axis=1))
    all_err_k = np.sqrt(np.mean((res["kf_T"]    - res["true_T"])**2, axis=1))
    ax7.plot(time_min, all_err_m, color="steelblue", ls="--", lw=2.2, label="Model only (no DA)")
    ax7.plot(time_min, all_err_k, color="crimson",   ls="-",  lw=2.2, label="Kalman filter (DA)")
    ax7.axhline(OBS_NOISE_STD, color="gray", ls=":", lw=1.5, label=f"Sensor noise {OBS_NOISE_STD:.2f}°C")
    ax7.axvline(time_min[t_eval], color="red", lw=1.0, ls=":", alpha=0.5, label="RMSE evaluation start")
    ax7.set_xlabel("Time [min]")
    ax7.set_ylabel("RMS error over all 10 nodes [°C]")
    ax7.set_title(
        "(7) Overall convergence  —  Model error grows without DA; DA suppresses error quickly",
        fontsize=10,
    )
    ax7.legend(fontsize=9); ax7.grid(True, alpha=0.3)

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


def print_summary(res: dict) -> None:
    t_eval = N_STEPS * 2 // 3
    print()
    print("=" * 65)
    print("  データ同化の効果まとめ（後半 20 分 RMSE）")
    print(f"  真値発熱: {Q_TRUE_MAIN} W/m / モデル発熱: {Q_NOMINAL} W/m (30% 過小)")
    print("=" * 65)
    print(f"  {'Node':<10} {'No DA (°C)':>12} {'With DA (°C)':>14} {'改善率':>10}")
    print("-" * 65)
    for i in range(N_NODES):
        r_m = rmse(res["true_T"][t_eval:, i], res["model_T"][t_eval:, i])
        r_k = rmse(res["true_T"][t_eval:, i], res["kf_T"][t_eval:, i])
        imp = (r_m - r_k) / r_m * 100 if r_m > 0 else 0.0
        tag = " [sensor]" if i in SENSOR_NODES else "         "
        print(f"  Node {i}{tag}  {r_m:>10.1f}   {r_k:>12.1f}   {imp:>8.0f}%")
    print("=" * 65)
    r_m_d = rmse(res["disp_true"][t_eval:], res["disp_model"][t_eval:])
    r_k_d = rmse(res["disp_true"][t_eval:], res["disp_kf"][t_eval:])
    imp_d = (r_m_d - r_k_d) / r_m_d * 100
    print(f"  Displacement      {r_m_d:>10.1f}   {r_k_d:>12.1f}   {imp_d:>8.0f}%  [μm]")
    print("=" * 65)


def plot_all_nodes(res: dict, save_path: str = "results_da_demo_all_nodes.png") -> None:
    """全 10 節点の真値 vs No-DA vs KF を 2 列 × 5 行で並べて保存する。"""
    time_min = res["time"] / 60.0
    t_eval = N_STEPS * 2 // 3

    fig, axes = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
    fig.suptitle(
        "All-Node Comparison: True vs No-DA vs Kalman Filter (DA)\n"
        f"True Q={Q_TRUE_MAIN:.0f} W/m,  Model Q={Q_NOMINAL:.0f} W/m (−30%).  "
        f"Sensors: Node 0 & Node {N_NODES-1}.  "
        f"Shaded band: KF ±2σ (hidden nodes only).",
        fontsize=11, y=0.995,
    )
    plt.subplots_adjust(hspace=0.45, wspace=0.3)

    for node in range(N_NODES):
        row, col = divmod(node, 2)
        ax = axes[row, col]
        is_sensor = node in SENSOR_NODES
        tag = "sensor" if is_sensor else "hidden"

        ax.plot(time_min, res["true_T"][:, node],
                color="black", ls="-", lw=2.0, label="True")
        ax.plot(time_min, res["model_T"][:, node],
                color="steelblue", ls="--", lw=1.5, label="No DA")
        ax.plot(time_min, res["kf_T"][:, node],
                color="crimson", ls="-", lw=1.8, label="DA (KF)")

        if not is_sensor:
            ax.fill_between(
                time_min,
                res["kf_T"][:, node] - 2 * res["kf_sigma"][:, node],
                res["kf_T"][:, node] + 2 * res["kf_sigma"][:, node],
                color="crimson", alpha=0.12,
            )

        if is_sensor:
            obs_idx = SENSOR_NODES.index(node)
            ax.plot(time_min, res["obs_y"][:, obs_idx],
                    color="gray", ls=":", lw=0.9, alpha=0.6, label="Sensor obs.")

        r_m = rmse(res["true_T"][t_eval:, node], res["model_T"][t_eval:, node])
        r_k = rmse(res["true_T"][t_eval:, node], res["kf_T"][t_eval:, node])
        imp = (r_m - r_k) / r_m * 100 if r_m > 0 else 0.0
        ax.text(
            0.98, 0.04,
            f"RMSE(last 30min)\nNo DA:{r_m:.1f}°C / DA:{r_k:.1f}°C ({imp:.0f}%↑)",
            transform=ax.transAxes, fontsize=7.5, va="bottom", ha="right",
            bbox=dict(boxstyle="round,pad=0.25", fc="lightyellow", ec="gray", alpha=0.85),
        )
        ax.axvline(time_min[t_eval], color="red", lw=0.8, ls=":", alpha=0.4)
        sensor_label = " ◆sensor" if is_sensor else ""
        ax.set_title(f"Node {node}{sensor_label}  [{tag}]", fontsize=10)
        ax.set_ylabel("Temp [°C]", fontsize=8)
        ax.grid(True, alpha=0.3)
        if row == 4:
            ax.set_xlabel("Time [min]")
        if node == 0:
            ax.legend(fontsize=7.5, loc="upper left")

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 65)
    print("  データ同化デモ: モデル発熱 30% 過小評価 + 暖機後コールドスタート")
    print(f"  真値: 暖機 {N_WARMUP * DT / 60:.0f} min + 本番 {Q_TRUE_MAIN} W/m")
    print(f"  モデル: {Q_NOMINAL} W/m (過小評価)")
    print("=" * 65)
    results = run()
    print_summary(results)
    out = os.path.join(os.path.dirname(__file__), "results_da_demo.png")
    plot_results(results, save_path=out)
    out_all = os.path.join(os.path.dirname(__file__), "results_da_demo_all_nodes.png")
    plot_all_nodes(results, save_path=out_all)
