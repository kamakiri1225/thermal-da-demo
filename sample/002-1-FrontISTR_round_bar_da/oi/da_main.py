"""
FrontISTR の丸棒データ同化デモ。

`sample/002-1_laplacian_da_round_bar` と同じ流れで、
真値、モデル単独、OI 補正後の結果を比較・可視化する。
"""

from __future__ import annotations

import shutil
import sys
import time as time_mod
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).parent.parent))
from fistr_interface import FrontISTRInterface

N_AXIAL = 40
ROD_LENGTH_M = 0.30
ROD_DIAMETER_M = 0.01
ROD_RADIUS_M = ROD_DIAMETER_M / 2.0
# OpenFOAM版 (002-1_laplacian_da_round_bar) と統一: DT=5s, 180ステップ = 900s
DT = 5.0
N_STEPS = 180

CYCLE_ON = 90    # ON継続ステップ数 = 450s
CYCLE_OFF = 30   # OFF継続ステップ数 = 150s
G_TRUE = 120.0
G_MODEL = 80.0

# OpenFOAM版と同じセンサ配置・OIパラメータ
ASSIM_SENSOR_AXIAL_NODES = [4, 16]       # x = 0.034m, 0.124m
VALIDATION_SENSOR_AXIAL_NODES = [28, 36] # x = 0.214m, 0.274m
OBS_NOISE_STD = 0.05
BACKGROUND_ERROR_STD = 3.0
CORRELATION_LENGTH_M = 0.05

T_INIT = 20.0
T_AMB = 20.0

BASE_DIR = Path(__file__).resolve().parent.parent
CASE_DIR = BASE_DIR / "case"
RESULTS_DIR = Path(__file__).resolve().parent / "results"
RESULTS_IMG_DIR = RESULTS_DIR / "img"


def make_gradient_series(n_steps: int, g_on: float) -> np.ndarray:
    """丸棒の加熱 ON/OFF 列を作る。"""
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([g_on if (k % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def make_left_flux_series(gradient_series: np.ndarray) -> np.ndarray:
    """与えた温度勾配列を、左端境界の熱流束列に変換する。"""
    conductivity = 160.0  # W/(mK) aluminum
    return conductivity * gradient_series


def save_history(path: Path, time_s: np.ndarray, values: np.ndarray, prefix: str) -> None:
    """時系列行列を CSV に保存する。各列は軸方向節点を表す。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.column_stack([time_s, values])
    header = "time_s," + ",".join(f"{prefix}{i}" for i in range(values.shape[1]))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.10g")


def load_frontistr_mesh(mesh_path: Path) -> tuple[np.ndarray, list[list[int]]]:
    """生成済みの FrontISTR メッシュから節点座標と要素接続を読む。"""
    points: list[tuple[float, float, float]] = []
    cells: list[list[int]] = []
    mode: str | None = None
    for line in mesh_path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("!NODE"):
            mode = "node"
            continue
        if stripped.startswith("!ELEMENT"):
            mode = "elem"
            continue
        if stripped.startswith("!"):
            mode = None
            continue
        if not stripped:
            continue
        parts = [p.strip() for p in stripped.split(",")]
        if mode == "node":
            points.append((float(parts[1]), float(parts[2]), float(parts[3])))
        elif mode == "elem":
            cells.append([int(p) - 1 for p in parts[1:]])
    return np.array(points, dtype=float), cells


def write_legacy_vtk(path: Path, points: np.ndarray, cells: list[list[int]], temperature: np.ndarray) -> None:
    """ParaView で開ける legacy VTK を 1 枚書き出す。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("FrontISTR round bar twin experiment\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")
        f.write(f"POINTS {len(points)} float\n")
        for x, y, z in points:
            f.write(f"{x:.10g} {y:.10g} {z:.10g}\n")
        f.write(f"CELLS {len(cells)} {len(cells) * 9}\n")
        for cell in cells:
            f.write("8 " + " ".join(str(idx) for idx in cell) + "\n")
        f.write(f"CELL_TYPES {len(cells)}\n")
        for _ in cells:
            f.write("12\n")
        f.write(f"POINT_DATA {len(points)}\n")
        f.write("SCALARS temperature_degC float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for value in temperature:
            f.write(f"{value:.10g}\n")


def export_vtk_series(case_name: str, history: np.ndarray) -> None:
    """時系列全体を VTK 群として出力する。ParaView のアニメーション用。"""
    points, cells = load_frontistr_mesh(CASE_DIR / "round_bar.msh")
    out_dir = RESULTS_DIR / case_name / "vtk"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for step_index, temperature in enumerate(history, start=1):
        write_legacy_vtk(out_dir / f"{case_name}_{step_index:04d}.vtk", points, cells, temperature)




def build_oi_matrices(fistr: FrontISTRInterface) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[int]]:
    """OI の観測行列 H、背景誤差共分散 B、観測誤差共分散 R を作る。"""
    representative = fistr.axial_representative_node_indices()
    n_nodes = fistr.n_nodes
    # 観測行列 H:
    #   y = H x + v
    # ここでは「どの節点をセンサとして読むか」を 1 行ずつ表す。
    H = np.zeros((len(ASSIM_SENSOR_AXIAL_NODES), n_nodes))
    for i, axial_node in enumerate(ASSIM_SENSOR_AXIAL_NODES):
        H[i, representative[axial_node]] = 1.0

    assert fistr.node_x is not None
    distance = np.abs(fistr.node_x[:, None] - fistr.node_x[None, :])
    # 背景誤差共分散 B:
    #   B_ij = σ_b^2 exp(-(d_ij^2)/(2L^2))
    #   σ_b : 背景誤差の大きさ
    #   d_ij: 節点 i と j の距離
    #   L   : 誤差相関長
    # 近い節点ほど一緒に補正され、遠い節点ほど弱く結び付く。
    B = BACKGROUND_ERROR_STD**2 * np.exp(-(distance**2) / (2.0 * CORRELATION_LENGTH_M**2))
    # 観測誤差共分散 R:
    #   R = σ_o^2 I
    # センサノイズは各観測で独立、同じ分散を持つと仮定する。
    R = np.eye(len(ASSIM_SENSOR_AXIAL_NODES)) * OBS_NOISE_STD**2
    return H, B, R, representative


def run_sequence(label: str, fistr: FrontISTRInterface, left_series: np.ndarray) -> np.ndarray:
    """データ同化なしで FrontISTR を回し、状態履歴を記録する。"""
    print(f"--- {label} ---")
    fistr.reset()
    state = np.full(fistr.n_nodes, T_INIT)
    history = np.zeros((N_STEPS, fistr.n_nodes))
    t0 = time_mod.time()
    for k in range(N_STEPS):
        state = fistr.run_step(state, left_flux=left_series[k], right_temp=T_AMB)
        history[k] = state
        if (k + 1) % 30 == 0:
            profile = fistr.axial_profile(state)
            print(f"  step {k + 1}/{N_STEPS}  T_left={profile[0]:.2f}  elapsed={time_mod.time() - t0:.0f}s")
    return history


def run_da(
    fistr: FrontISTRInterface,
    truth_t: np.ndarray,
    left_model: np.ndarray,
    H: np.ndarray,
    B: np.ndarray,
    R: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """FrontISTR の予測ステップの後に OI 補正をかける。"""
    print("--- FrontISTR + OI data assimilation ---")
    fistr.reset()
    # OI の更新式:
    #   x_a = x_f + K (y - H x_f)
    #   K   = B H^T (H B H^T + R)^-1
    # ここで x_f は FrontISTR の予測場、y は観測、x_a は同化後の場。
    gain = B @ H.T @ np.linalg.inv(H @ B @ H.T + R)
    state = np.full(fistr.n_nodes, T_INIT)
    forecast_t = np.zeros((N_STEPS, fistr.n_nodes))
    da_t = np.zeros((N_STEPS, fistr.n_nodes))
    t0 = time_mod.time()
    for k in range(N_STEPS):
        # 真値場から観測を作る。観測ノイズは N(0, sigma^2) を仮定する。
        y = H @ truth_t[k] + np.random.randn(H.shape[0]) * OBS_NOISE_STD
        # 予測ステップ:
        #   x_f^{k+1} = M(x_a^k, u^k)
        # 直前の同化後状態 state を FrontISTR で 1 ステップ進める。
        forecast = fistr.run_step(state, left_flux=left_model[k], right_temp=T_AMB)
        # イノベーション(観測差):
        #   d = y - H x_f
        innovation = y - H @ forecast
        # OI 更新:
        #   x_a = x_f + K d
        # つまり、観測との差分 d を、B と R から作った重み K で
        # 空間全体へ広げる。
        state = forecast + gain @ innovation
        forecast_t[k] = forecast
        da_t[k] = state
        if (k + 1) % 30 == 0:
            p_true = fistr.axial_profile(truth_t[k])
            p_fcst = fistr.axial_profile(forecast)
            p_da = fistr.axial_profile(state)
            print(
                f"  step {k + 1}/{N_STEPS}"
                f"  true={p_true[0]:.2f}  fistr={p_fcst[0]:.2f}  oi={p_da[0]:.2f}"
                f"  elapsed={time_mod.time() - t0:.0f}s"
            )
    return forecast_t, da_t


def axial_history(fistr: FrontISTRInterface, values: np.ndarray) -> np.ndarray:
    """節点温度履歴を軸方向平均へ畳み込む。"""
    return np.vstack([fistr.axial_profile(row) for row in values])


def save_summary(fistr: FrontISTRInterface, truth_t: np.ndarray, model_t: np.ndarray, da_t: np.ndarray) -> None:
    """モデル単独と同化ありの軸方向 RMSE を保存する。"""
    t_start = 2 * N_STEPS // 3
    truth = axial_history(fistr, truth_t)
    model = axial_history(fistr, model_t)
    da = axial_history(fistr, da_t)
    rmse_model = np.sqrt(np.mean((truth[t_start:] - model[t_start:]) ** 2, axis=0))
    rmse_da = np.sqrt(np.mean((truth[t_start:] - da[t_start:]) ** 2, axis=0))
    data = np.column_stack([np.arange(N_AXIAL), fistr.axial_positions(), rmse_model, rmse_da])
    np.savetxt(
        RESULTS_DIR / "summary_rmse.csv",
        data,
        delimiter=",",
        header="axial_node,x_m,rmse_frontistr_only,rmse_oi_da",
        comments="",
        fmt=["%d", "%.10g", "%.10g", "%.10g"],
    )
    print()
    print("=" * 64)
    print("FrontISTR twin experiment RMSE, last third of run")
    print(f"  FrontISTR-only avg RMSE: {float(np.mean(rmse_model)):.4f} degC")
    print(f"  FrontISTR + OI avg RMSE: {float(np.mean(rmse_da)):.4f} degC")
    if float(np.mean(rmse_model)) > 0.0:
        print(f"  Improvement: {(1.0 - float(np.mean(rmse_da)) / float(np.mean(rmse_model))) * 100.0:.1f}%")
    print("=" * 64)


def save_node_temperature_comparison(
    fistr: FrontISTRInterface,
    truth_t: np.ndarray,
    model_t: np.ndarray,
    da_t: np.ndarray,
) -> None:
    """節点ごとの温度比較 CSV を保存する。軸方向表示にも使う。"""
    time_s = (np.arange(N_STEPS) + 1) * DT
    assert fistr.node_x is not None
    coords = fistr._node_coordinates()
    rows = []
    for k, t in enumerate(time_s):
        for i, node_id in enumerate(fistr.node_ids):
            rows.append(
                [
                    t,
                    node_id,
                    coords[i, 0],
                    coords[i, 1],
                    coords[i, 2],
                    truth_t[k, i],
                    model_t[k, i],
                    da_t[k, i],
                ]
            )
    np.savetxt(
        RESULTS_DIR / "frontistr_node_temperature_comparison.csv",
        np.array(rows),
        delimiter=",",
        header="time_s,node_id,x_m,y_m,z_m,truth,without_da,with_da",
        comments="",
        fmt=["%.10g", "%d", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g", "%.10g"],
    )

    truth_axial = axial_history(fistr, truth_t)
    model_axial = axial_history(fistr, model_t)
    da_axial = axial_history(fistr, da_t)
    axial_rows = []
    for k, t in enumerate(time_s):
        for node in range(N_AXIAL):
            axial_rows.append(
                [
                    t,
                    node,
                    fistr.axial_positions()[node],
                    truth_axial[k, node],
                    model_axial[k, node],
                    da_axial[k, node],
                ]
            )
    np.savetxt(
        RESULTS_DIR / "axial_node_temperature_comparison.csv",
        np.array(axial_rows),
        delimiter=",",
        header="time_s,axial_node,x_m,truth,without_da,with_da",
        comments="",
        fmt=["%.10g", "%d", "%.10g", "%.10g", "%.10g", "%.10g"],
    )

    measurement_nodes = ASSIM_SENSOR_AXIAL_NODES + VALIDATION_SENSOR_AXIAL_NODES
    measurement_rows = []
    for k, t in enumerate(time_s):
        for node in measurement_nodes:
            measurement_rows.append(
                [
                    t,
                    node,
                    fistr.axial_positions()[node],
                    truth_axial[k, node],
                    model_axial[k, node],
                    da_axial[k, node],
                ]
            )
    np.savetxt(
        RESULTS_DIR / "measurement_node_temperature_comparison.csv",
        np.array(measurement_rows),
        delimiter=",",
        header="time_s,axial_node,x_m,truth,without_da,with_da",
        comments="",
        fmt=["%.10g", "%d", "%.10g", "%.10g", "%.10g", "%.10g"],
    )


def plot_results(
    fistr: FrontISTRInterface,
    truth_t: np.ndarray,
    model_t: np.ndarray,
    da_t: np.ndarray,
    gradient_true: np.ndarray,
) -> None:
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    time_min = (np.arange(N_STEPS) + 1) * DT / 60.0
    t_start = 2 * N_STEPS // 3
    truth = axial_history(fistr, truth_t)
    model = axial_history(fistr, model_t)
    da = axial_history(fistr, da_t)
    hidden_to_plot = VALIDATION_SENSOR_AXIAL_NODES + [8, 24]
    colors = plt.cm.tab10(np.linspace(0, 1, N_AXIAL))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Round Bar Twin Experiment: FrontISTR + Optimal Interpolation\n"
        f"phi={ROD_DIAMETER_M * 1000:.0f} mm, L={ROD_LENGTH_M * 1000:.0f} mm, "
        f"Truth G={G_TRUE:g} degC/m, Model G={G_MODEL:g} degC/m",
        fontsize=11,
    )

    ax = axes[0, 0]
    for node in ASSIM_SENSOR_AXIAL_NODES:
        c = colors[node]
        ax.plot(time_min, truth[:, node], color=c, ls="--", lw=1.8, label=f"N{node} truth")
        ax.plot(time_min, model[:, node], color=c, ls=":", lw=1.3, label=f"N{node} FrontISTR only")
        ax.plot(time_min, da[:, node], color=c, lw=2.0, label=f"N{node} OI")
    ax.set_title("(1) Assimilation sensor temperatures", fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [degC]")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    for node in hidden_to_plot:
        c = colors[node]
        ax.plot(time_min, truth[:, node], color=c, ls="--", lw=1.8, label=f"N{node} truth")
        ax.plot(time_min, model[:, node], color=c, ls=":", lw=1.3, label=f"N{node} FrontISTR only")
        ax.plot(time_min, da[:, node], color=c, lw=2.0, label=f"N{node} OI")
    ax.set_title("(2) Validation / hidden node temperatures", fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [degC]")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    rmse_model = np.sqrt(np.mean((truth[t_start:] - model[t_start:]) ** 2, axis=0))
    rmse_da = np.sqrt(np.mean((truth[t_start:] - da[t_start:]) ** 2, axis=0))
    x = np.arange(N_AXIAL)
    ax.bar(x - 0.2, rmse_model, 0.4, label="FrontISTR only", color="steelblue", alpha=0.85, edgecolor="k")
    ax.bar(x + 0.2, rmse_da, 0.4, label="FrontISTR + OI", color="crimson", alpha=0.85, edgecolor="k")
    ax.axhline(OBS_NOISE_STD, color="gray", ls=":", lw=1.5, label="sensor noise")
    ax.set_xticks(x[::2])
    ax.set_xticklabels([f"N{i}" for i in x[::2]], fontsize=8)
    ax.set_title("(3) RMSE over axial profile", fontsize=10)
    ax.set_ylabel("RMSE [degC]")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    ax = axes[1, 1]
    for node in hidden_to_plot:
        c = colors[node]
        ax.plot(time_min, truth[:, node] - model[:, node], color=c, ls=":", lw=1.3, label=f"N{node} FrontISTR only")
        ax.plot(time_min, truth[:, node] - da[:, node], color=c, lw=1.8, label=f"N{node} OI")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(time_min[t_start], color="red", lw=1.0, ls="--", alpha=0.5)
    ax.set_title("(4) Hidden-node error, truth - estimate", fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Error [degC]")
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    out = RESULTS_IMG_DIR / "results_frontistr_da.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[saved] {out}")

    save_history(RESULTS_DIR / "truth" / "axial_temperature_history.csv", time_min * 60.0, truth, "N")
    save_history(RESULTS_DIR / "model_only" / "axial_temperature_history.csv", time_min * 60.0, model, "N")
    save_history(RESULTS_DIR / "with_da" / "axial_temperature_history.csv", time_min * 60.0, da, "N")
    save_history(RESULTS_DIR / "input_gradient.csv", time_min * 60.0, gradient_true[:, None], "G_true")


def plot_measurement_points(
    fistr: FrontISTRInterface,
    truth_t: np.ndarray,
    model_t: np.ndarray,
    da_t: np.ndarray,
) -> None:
    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    time_min = (np.arange(N_STEPS) + 1) * DT / 60.0
    truth = axial_history(fistr, truth_t)
    model = axial_history(fistr, model_t)
    da = axial_history(fistr, da_t)
    nodes = ASSIM_SENSOR_AXIAL_NODES + VALIDATION_SENSOR_AXIAL_NODES

    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    fig.suptitle(
        "FrontISTR Round Bar: Measurement Point Comparison\n"
        "truth / without data assimilation / with data assimilation",
        fontsize=12,
    )
    for ax, node in zip(axes.ravel(), nodes):
        x_m = fistr.axial_positions()[node]
        ax.plot(time_min, truth[:, node], color="black", ls="--", lw=2.0, label="truth")
        ax.plot(time_min, model[:, node], color="steelblue", ls=":", lw=2.0, label="without DA")
        ax.plot(time_min, da[:, node], color="crimson", lw=1.8, label="with DA")
        role = "assimilation sensor" if node in ASSIM_SENSOR_AXIAL_NODES else "validation point"
        ax.set_title(f"N{node}  x={x_m:.3f} m  ({role})", fontsize=10)
        ax.set_ylabel("Temperature [degC]")
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


def main() -> None:
    np.random.seed(42)
    print("=" * 64)
    print("FrontISTR round-bar twin experiment + OI")
    print("=" * 64)
    print(f"rod: phi={ROD_DIAMETER_M * 1000:.0f} mm, L={ROD_LENGTH_M * 1000:.0f} mm")
    print(f"N_AXIAL={N_AXIAL}, N_STEPS={N_STEPS}, DT={DT}s")
    print(f"assim sensors={ASSIM_SENSOR_AXIAL_NODES}, validation={VALIDATION_SENSOR_AXIAL_NODES}")

    RESULTS_IMG_DIR.mkdir(parents=True, exist_ok=True)

    fistr = FrontISTRInterface(
        CASE_DIR,
        n_axial=N_AXIAL,
        length_m=ROD_LENGTH_M,
        radius_m=ROD_RADIUS_M,
        t_amb=T_AMB,
        dt=DT,
    )
    fistr.reset()
    print(f"FrontISTR mesh nodes={fistr.n_nodes}")
    H, B, R, _representative = build_oi_matrices(fistr)

    g_true = make_gradient_series(N_STEPS, G_TRUE)
    g_model = make_gradient_series(N_STEPS, G_MODEL)
    left_true = make_left_flux_series(g_true)
    left_model = make_left_flux_series(g_model)

    truth_t = run_sequence("truth run: FrontISTR with G_true", fistr, left_true)
    save_history(RESULTS_DIR / "truth" / "temperature_history.csv", (np.arange(N_STEPS) + 1) * DT, truth_t, "node")

    print()
    model_t = run_sequence("model-only run: FrontISTR with G_model", fistr, left_model)
    save_history(RESULTS_DIR / "model_only" / "temperature_history.csv", (np.arange(N_STEPS) + 1) * DT, model_t, "node")

    print()
    _forecast_t, da_t = run_da(fistr, truth_t, left_model, H, B, R)
    save_history(RESULTS_DIR / "with_da" / "temperature_history.csv", (np.arange(N_STEPS) + 1) * DT, da_t, "node")

    save_summary(fistr, truth_t, model_t, da_t)
    save_node_temperature_comparison(fistr, truth_t, model_t, da_t)
    plot_results(fistr, truth_t, model_t, da_t, g_true)
    plot_measurement_points(fistr, truth_t, model_t, da_t)
    export_vtk_series("truth", truth_t)
    export_vtk_series("model_only", model_t)
    export_vtk_series("with_da", da_t)


if __name__ == "__main__":
    main()
