"""
main_disp_obs.py
================
変位データを観測に追加した場合の比較実験

【比較するケース】
  Case A : 温度センサ 2 点のみ           (H: 2×10)
  Case B : 温度センサ 2 点 + 変位センサ  (H: 3×10)
  Case C : 変位センサのみ                (H: 1×10)

【ポイント】
  変位 δ は全節点温度の線形和：
      δ [μm] = α_exp · Σ_i (θ_i − T_ref) · dx · 1e6

  そのため観測行列 H に変位行を追加できる：
      H_disp = α_exp · dx · 1e6 · [1, 1, 1, 1, 1]

  変位センサを加えると「全体の温度平均」に相当する情報が増え、
  FEM モデルによる空間補間と組み合わさって精度が向上する。
"""

from __future__ import annotations

import os
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from thermal_model import ThermalBarModel
from kalman_filter import KalmanFilter

# ── 日本語フォント設定 ─────────────────────────────────────────
def _setup_jp_font() -> None:
    import matplotlib.font_manager as fm
    candidates = ["IPAexGothic","IPAGothic","Noto Sans CJK JP",
                  "Yu Gothic","Meiryo","MS Gothic"]
    installed = {f.name for f in fm.fontManager.ttflist}
    for name in candidates:
        if name in installed:
            matplotlib.rcParams["font.family"] = name
            return
    for path in ["/mnt/c/Windows/Fonts/YuGothM.ttc",
                 "/mnt/c/Windows/Fonts/meiryo.ttc",
                 "/mnt/c/Windows/Fonts/msgothic.ttc"]:
        if os.path.exists(path):
            fm.fontManager.addfont(path)
            matplotlib.rcParams["font.family"] = fm.FontProperties(fname=path).get_name()
            return

_setup_jp_font()
matplotlib.rcParams["axes.unicode_minus"] = False

# ══════════════════════════════════════════════════════════════
# パラメータ（main.py と同一）
# ══════════════════════════════════════════════════════════════
N_NODES      = 10
DT           = 10.0
N_STEPS      = 540
SENSOR_NODES = [0, N_NODES - 1]

ALPHA_DIFF   = 1.2e-5
DX           = 0.10
H_CONV       = 8.0
RHO_CP       = 3.5e6
ALPHA_EXP    = 12.0e-6
T_INIT       = 20.0

Q_ON         = 6.5e3
CYCLE_ON     = 60
CYCLE_OFF    = 30

PROC_NOISE_STD = 0.03
OBS_NOISE_STD  = 0.20
DISP_NOISE_STD = 0.50   # 変位センサノイズ [μm]（リニアエンコーダ相当）

# 変位行列の係数：δ[μm] = c_disp * Σ(θ_i − T_ref)
C_DISP = ALPHA_EXP * DX * 1e6   # [μm/°C] per node


# ══════════════════════════════════════════════════════════════
# ヘルパー
# ══════════════════════════════════════════════════════════════
def make_heat_input(n_steps: int) -> np.ndarray:
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([Q_ON if (k % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def true_displacement_um(temps: np.ndarray) -> np.ndarray:
    """全節点温度から真の変位 [μm] を計算"""
    return C_DISP * np.sum(temps - T_INIT, axis=1)


def rmse(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


# ══════════════════════════════════════════════════════════════
# 観測行列の構築
# ══════════════════════════════════════════════════════════════
def build_H_and_R(case: str) -> tuple[np.ndarray, np.ndarray]:
    """
    case = 'A' : 温度センサ 2 点
    case = 'B' : 温度センサ 2 点 + 変位センサ 1 点
    case = 'C' : 変位センサ 1 点のみ
    """
    # 温度センサ行
    H_temp = np.zeros((len(SENSOR_NODES), N_NODES))
    for i, node in enumerate(SENSOR_NODES):
        H_temp[i, node] = 1.0
    R_temp = np.eye(len(SENSOR_NODES)) * OBS_NOISE_STD ** 2

    # 変位センサ行：δ[μm] = C_DISP * Σ(θ_i − T_ref)
    #   y_disp_eff = δ_measured + C_DISP*n*T_ref  ← これを観測値として渡す
    #              = C_DISP * Σ θ_i  + v_disp
    H_disp = C_DISP * np.ones((1, N_NODES))
    R_disp = np.array([[DISP_NOISE_STD ** 2]])

    if case == "A":
        return H_temp, R_temp
    elif case == "B":
        H = np.vstack([H_temp, H_disp])
        R = np.block([[R_temp, np.zeros((2, 1))],
                      [np.zeros((1, 2)), R_disp]])
        return H, R
    elif case == "C":
        return H_disp, R_disp
    else:
        raise ValueError(f"Unknown case: {case}")


# ══════════════════════════════════════════════════════════════
# シミュレーション（1ケース分）
# ══════════════════════════════════════════════════════════════
def run_case(case: str, x_true_series: np.ndarray, u_series: np.ndarray) -> dict:
    """
    x_true_series : shape (N_STEPS, N_NODES)  ← あらかじめ生成した真値
    """
    model = ThermalBarModel(
        n_nodes=N_NODES, dt=DT,
        alpha=ALPHA_DIFF, dx=DX,
        h_conv=H_CONV, rho_cp=RHO_CP,
    )
    H, R = build_H_and_R(case)
    Q_mat = np.eye(N_NODES) * PROC_NOISE_STD ** 2
    x0 = np.ones(N_NODES) * T_INIT
    P0 = np.eye(N_NODES) * 2.0

    kf = KalmanFilter(A=model.A_d, B=model.B_d, H=H,
                      Q=Q_mat, R=R, x0=x0, P0=P0)

    kf_T     = np.zeros((N_STEPS, N_NODES))
    kf_sigma = np.zeros((N_STEPS, N_NODES))

    for k in range(N_STEPS):
        u      = np.array([u_series[k]])
        x_true = x_true_series[k]

        # 観測値の生成（ケースに応じて）
        if case == "A":
            y = H @ x_true + np.random.randn(H.shape[0]) * OBS_NOISE_STD
        elif case == "B":
            # 温度部分のノイズ
            y_temp = H[:2] @ x_true + np.random.randn(2) * OBS_NOISE_STD
            # 変位部分（実効観測値 = δ_measured + C_DISP*n*T_ref）
            disp_true_um = C_DISP * np.sum(x_true - T_INIT)
            y_disp_eff   = (disp_true_um + C_DISP * N_NODES * T_INIT
                            + np.random.randn() * DISP_NOISE_STD)
            y = np.append(y_temp, y_disp_eff)
        elif case == "C":
            disp_true_um = C_DISP * np.sum(x_true - T_INIT)
            y_disp_eff   = (disp_true_um + C_DISP * N_NODES * T_INIT
                            + np.random.randn() * DISP_NOISE_STD)
            y = np.array([y_disp_eff])

        x_est, P_est = kf.step(u, y)
        kf_T[k]     = x_est
        kf_sigma[k] = np.sqrt(np.diag(P_est))

    return {
        "kf_T"    : kf_T,
        "kf_sigma": kf_sigma,
        "disp_kf" : true_displacement_um(kf_T),
    }


# ══════════════════════════════════════════════════════════════
# 統計表示
# ══════════════════════════════════════════════════════════════
def print_comparison(cases: dict[str, dict], true_T: np.ndarray, disp_true: np.ndarray) -> None:
    t_start = 2 * N_STEPS // 3
    labels  = {"A": "温度2点のみ     ", "B": "温度2点+変位1点", "C": "変位1点のみ     "}
    hidden  = [i for i in range(N_NODES) if i not in SENSOR_NODES]

    print()
    print("=" * 72)
    print("  観測条件別の推定精度比較（定常域 RMSE）")
    print("=" * 72)
    print(f"  {'ケース':<20}  {'未観測平均':>10}  {'最大RMSE':>10}  {'変位[μm]':>10}")
    print("-" * 72)
    for key, res in cases.items():
        rmse_T    = [rmse(true_T[t_start:, i], res["kf_T"][t_start:, i]) for i in range(N_NODES)]
        rmse_disp = rmse(disp_true[t_start:], res["disp_kf"][t_start:])
        hidden_rmse = [rmse_T[i] for i in hidden]
        print(f"  {labels[key]}  "
              f"  {np.mean(hidden_rmse):>9.4f}  "
              f"  {np.max(hidden_rmse):>9.4f}  "
              f"  {rmse_disp:>9.4f}")
    print("=" * 72)
    print()


# ══════════════════════════════════════════════════════════════
# プロット
# ══════════════════════════════════════════════════════════════
def plot_comparison(cases: dict[str, dict], true_T: np.ndarray,
                    disp_true: np.ndarray, u_series: np.ndarray,
                    save_path: str = "results_disp_obs.png") -> None:
    time_min = np.arange(N_STEPS) * DT / 60.0
    t_start  = 2 * N_STEPS // 3

    case_styles = {
        "A": dict(color="navy",       ls="-",  lw=1.8, label="Case A: Temp x2 only"),
        "B": dict(color="crimson",    ls="-",  lw=2.0, label="Case B: Temp x2 + Disp"),
        "C": dict(color="darkorange", ls="-.", lw=1.5, label="Case C: Disp only"),
    }
    style_true = dict(color="black", ls="-", lw=2.2, label="True")

    fig = plt.figure(figsize=(16, 20))
    fig.suptitle(
        "Effect of Adding Displacement Observation\n"
        "Case A: 2 temp sensors  /  Case B: 2 temp + 1 disp  /  Case C: disp only\n"
        "All cases use identical initial conditions — only sensor configuration differs",
        fontsize=11, y=0.998,
    )
    gs = gridspec.GridSpec(4, 2, figure=fig, hspace=0.55, wspace=0.35)

    def _rmse_box(ax, node):
        lines = []
        for key in ["A", "B", "C"]:
            r = rmse(true_T[t_start:, node], cases[key]["kf_T"][t_start:, node])
            lines.append(f"{key}: {r:.3f} degC")
        ax.text(0.97, 0.05, "RMSE (last 30min)\n" + "\n".join(lines),
                transform=ax.transAxes, fontsize=8, va="bottom", ha="right",
                bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow", ec="gray", alpha=0.9))
        ax.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)

    # Row 0: センサ節点（Node 0 と Node 9）
    for col, (node, sfx) in enumerate([(0, " [sensor, heat input end]"),
                                        (N_NODES - 1, " [sensor, cool end]")]):
        ax = fig.add_subplot(gs[0, col])
        ax.plot(time_min, true_T[:, node], **style_true, zorder=5)
        for key, res in cases.items():
            ax.plot(time_min, res["kf_T"][:, node], **case_styles[key])
        _rmse_box(ax, node)
        ax.set_title(f"({'1' if col == 0 else '2'}) Node {node}{sfx}", fontsize=10)
        ax.set_xlabel("Time [min]"); ax.set_ylabel("Temperature [degC]")
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)

    # Row 1: 未観測節点（Node 1 と Node 5）
    for col, (node, sfx) in enumerate([(1, " [hidden, near heat source]"),
                                        (5, " [hidden, center]")]):
        ax = fig.add_subplot(gs[1, col])
        ax.plot(time_min, true_T[:, node], **style_true, zorder=5)
        for key, res in cases.items():
            ax.plot(time_min, res["kf_T"][:, node], **case_styles[key])
        ax.fill_between(
            time_min,
            cases["B"]["kf_T"][:, node] - 2 * cases["B"]["kf_sigma"][:, node],
            cases["B"]["kf_T"][:, node] + 2 * cases["B"]["kf_sigma"][:, node],
            color="crimson", alpha=0.12, label="Case B ±2sigma",
        )
        _rmse_box(ax, node)
        ax.set_title(f"({'3' if col == 0 else '4'}) Node {node}{sfx}", fontsize=10)
        ax.set_xlabel("Time [min]"); ax.set_ylabel("Temperature [degC]")
        ax.legend(fontsize=7.5); ax.grid(True, alpha=0.3)

    # Row 2: 熱変位 と 全節点RMSE棒グラフ
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.plot(time_min, disp_true, **style_true, zorder=5)
    for key, res in cases.items():
        r = rmse(disp_true[t_start:], res["disp_kf"][t_start:])
        s = case_styles[key].copy()
        s["label"] = f"{s['label']}  RMSE={r:.3f}um"
        ax5.plot(time_min, res["disp_kf"], **s)
    ax5.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)
    ax5.set_title("(5) Thermal displacement estimation", fontsize=10)
    ax5.set_xlabel("Time [min]"); ax5.set_ylabel("Displacement [um]")
    ax5.legend(fontsize=7.5); ax5.grid(True, alpha=0.3)

    ax6 = fig.add_subplot(gs[2, 1])
    x = np.arange(N_NODES)
    w = 0.25
    for idx, key in enumerate(["A", "B", "C"]):
        rmse_list = [rmse(true_T[t_start:, n], cases[key]["kf_T"][t_start:, n])
                     for n in range(N_NODES)]
        bars = ax6.bar(x + (idx - 1) * w, rmse_list, w,
                       color=case_styles[key]["color"], alpha=0.85,
                       edgecolor="k", label=case_styles[key]["label"])
    ax6.set_xticks(x)
    ax6.set_xticklabels([f"N{i}" for i in range(N_NODES)], fontsize=9)
    ax6.set_ylabel("RMSE [degC]")
    ax6.set_title("(6) Node-by-node RMSE (last 30 min)", fontsize=10)
    ax6.legend(fontsize=8); ax6.grid(True, alpha=0.3, axis="y")
    ymax = ax6.get_ylim()[1]
    for node in SENSOR_NODES:
        ax6.text(node, ymax * 0.92, "sensor", ha="center", fontsize=8, color="dimgray")

    # Row 3: 全節点平均RMSE の時系列（main_da_demo.py と同様）
    ax7 = fig.add_subplot(gs[3, :])
    for key, res in cases.items():
        all_err = np.sqrt(np.mean((res["kf_T"] - true_T) ** 2, axis=1))
        ax7.plot(time_min, all_err, **case_styles[key])
    ax7.axvline(time_min[t_start], color="red", lw=1.0, ls=":", alpha=0.5)
    ax7.set_title("(7) All-node mean RMSE over time", fontsize=10)
    ax7.set_xlabel("Time [min]"); ax7.set_ylabel("Mean RMSE [degC]")
    ax7.legend(fontsize=8); ax7.grid(True, alpha=0.3)

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


def plot_all_nodes_comparison(
    cases: dict,
    true_T: np.ndarray,
    save_path: str = "results_disp_obs_all_nodes.png",
) -> None:
    """全 10 節点について Case A / B / C vs 真値を 5×2 パネルで比較する。"""
    time_min = np.arange(N_STEPS) * DT / 60.0
    t_start  = 2 * N_STEPS // 3

    case_styles = {
        "A": dict(color="navy",       ls="-",  lw=1.5, label="Case A: Temp x2"),
        "B": dict(color="crimson",    ls="-",  lw=1.8, label="Case B: Temp x2 + Disp"),
        "C": dict(color="darkorange", ls="-.", lw=1.5, label="Case C: Disp only"),
    }

    fig, axes = plt.subplots(5, 2, figsize=(14, 20), sharex=True)
    fig.suptitle(
        "All-Node Comparison by Observation Case\n"
        "Case A: 2 temp sensors only  |  "
        "Case B: 2 temp + 1 disp  |  "
        "Case C: disp only\n"
        "All cases use identical initial conditions and heat model — only sensor configuration differs",
        fontsize=9, y=0.998,
    )
    plt.subplots_adjust(hspace=0.45, wspace=0.3)

    for node in range(N_NODES):
        row, col = divmod(node, 2)
        ax = axes[row, col]
        is_sensor = node in SENSOR_NODES

        ax.plot(time_min, true_T[:, node],
                color="black", ls="-", lw=2.2, label="True", zorder=5)

        for key, res in cases.items():
            ax.plot(time_min, res["kf_T"][:, node], **case_styles[key])

        # RMSE box
        lines = []
        for key in ["A", "B", "C"]:
            r = rmse(true_T[t_start:, node], cases[key]["kf_T"][t_start:, node])
            lines.append(f"{key}: {r:.2f}°C")
        ax.text(0.98, 0.04, "RMSE (last 30min)\n" + "\n".join(lines),
                transform=ax.transAxes, fontsize=7, va="bottom", ha="right",
                bbox=dict(boxstyle="round,pad=0.25", fc="lightyellow", ec="gray", alpha=0.85))

        ax.axvline(time_min[t_start], color="red", lw=0.8, ls=":", alpha=0.4)
        tag = "sensor" if is_sensor else "hidden"
        ax.set_title(f"Node {node}  [{tag}]", fontsize=10)
        ax.set_ylabel("Temp [degC]", fontsize=8)
        ax.grid(True, alpha=0.3)
        if row == 4:
            ax.set_xlabel("Time [min]")
        if node == 0:
            ax.legend(fontsize=8, loc="upper left")

        # y軸の最小スパンを確保（温度変化が小さい節点で差が誇張されないよう）
        ymin, ymax = ax.get_ylim()
        MIN_SPAN = 6.0  # °C
        if (ymax - ymin) < MIN_SPAN:
            ymid = (ymin + ymax) / 2
            ax.set_ylim(ymid - MIN_SPAN / 2, ymid + MIN_SPAN / 2)

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.show()


# ══════════════════════════════════════════════════════════════
# エントリーポイント
# ══════════════════════════════════════════════════════════════
if __name__ == "__main__":
    np.random.seed(42)

    print("変位データ観測比較シミュレーション 開始 ...")

    # ── 真値の生成（3ケース共通で使う） ──────────────────────
    model = ThermalBarModel(n_nodes=N_NODES, dt=DT, alpha=ALPHA_DIFF,
                            dx=DX, h_conv=H_CONV, rho_cp=RHO_CP)
    u_series  = make_heat_input(N_STEPS)
    x_true_all = np.zeros((N_STEPS, N_NODES))
    x_true     = np.ones(N_NODES) * T_INIT
    for k in range(N_STEPS):
        x_true = model.step(x_true, np.array([u_series[k]]), noise_std=PROC_NOISE_STD)
        x_true_all[k] = x_true

    disp_true = true_displacement_um(x_true_all)

    # ── 各ケースでKFを実行 ────────────────────────────────────
    cases = {}
    for case in ["A", "B", "C"]:
        np.random.seed(42)   # 観測ノイズのシードを揃えてケース間の比較を公平に
        cases[case] = run_case(case, x_true_all, u_series)
        print(f"  Case {case} 完了")

    print_comparison(cases, x_true_all, disp_true)

    save_path = os.path.join(os.path.dirname(__file__), "results_disp_obs.png")
    plot_comparison(cases, x_true_all, disp_true, u_series, save_path=save_path)

    save_all = os.path.join(os.path.dirname(__file__), "results_disp_obs_all_nodes.png")
    plot_all_nodes_comparison(cases, x_true_all, save_path=save_all)
