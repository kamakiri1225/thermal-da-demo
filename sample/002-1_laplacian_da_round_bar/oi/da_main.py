"""
da_main.py  (oi/)
==================
laplacianFoam を使った円柱丸棒モデルの Optimal Interpolation データ同化デモ。

【ツイン実験の構成】
  真値生成:   laplacianFoam で G_true (= 100 degC/m) の熱入力で走らせる
  予測ステップ: laplacianFoam で G_model (= 50 degC/m, 50% 少ない) で走らせる
  更新ステップ: Optimal Interpolation (OI) がセンサ観測で予測を補正する

  真値と予測モデルは同じ OpenFOAM ソルバー・同じメッシュを使う。
  ただし境界勾配 G を意図的にずらして、モデル誤差がある状況を作る。

【実行方法】
    cd sample/002-1_laplacian_da_round_bar/oi
    python da_main.py
  または
    cd sample/002-1_laplacian_da_round_bar/oi
    ./run.sh

【比較内容】
  OF-only (DA なし) : G_model で OpenFOAM を独立に走らせた予測（センサ補正なし）
  OI あり (DA あり) : G_model + OI でセンサ観測を使って補正した推定
  真値              : G_true で OF を走らせた「正解」の温度場
"""

from __future__ import annotations

import sys
import time as time_mod
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# of_interface.py は親ディレクトリに置いてモデル共有
sys.path.insert(0, str(Path(__file__).parent.parent))
from of_interface import OFInterface

# ─────────────────────────────────────────────
# シミュレーション設定
# ─────────────────────────────────────────────
# 円柱丸棒モデル:
# 市販で扱いやすいアルミ丸棒 φ10 mm x 300 mm を想定。
# blockMesh は中央ブロック + 外周4ブロックの円柱断面を軸方向40分割する。
# 断面は中央4セル + 外周8セル = 12セル。
N_AXIAL = 40
N_CROSS_SECTION = 12
N_CELLS = N_AXIAL * N_CROSS_SECTION
ROD_LENGTH_M = 0.30
ROD_DIAMETER_M = 0.01
DX = ROD_LENGTH_M / N_AXIAL

DT = 5.0            # 時間刻み [s]
N_STEPS = 180       # 900s = 15分

# 熱入力サイクル（ON/OFF スイッチング）
CYCLE_ON = 90       # ON 継続ステップ数 (450s)
CYCLE_OFF = 30      # OFF 継続ステップ数 (150s)

# 境界温度勾配 [degC/m]
# 丸棒長さ L=0.3m, 定常状態左端温度: T[0] = T_AMB + G_TRUE * L = 56 degC
G_TRUE = 120.0      # 真値の熱入力境界勾配
G_MODEL = 80.0      # モデルの境界勾配（真値より弱い加熱を模擬）

# センサ・ノイズ設定
ASSIM_SENSOR_AXIAL_NODES = [4, 16]       # OI 同化に使うセンサ: x=0.034, 0.124 m
VALIDATION_SENSOR_AXIAL_NODES = [28, 36] # 検証専用センサ: x=0.214, 0.274 m
OBS_NOISE_STD = 0.05    # 観測ノイズ標準偏差 [degC]
BACKGROUND_ERROR_STD = 3.0      # OpenFOAM予測の背景誤差標準偏差 [degC]
CORRELATION_LENGTH_M = 0.05     # 補正を空間的に広げる相関長 [m]

# 初期・周囲温度
T_INIT = 20.0
T_AMB = 20.0

CASE_DIR = Path(__file__).parent / "case"
RESULTS_DIR = Path(__file__).parent / "results"


# ─────────────────────────────────────────────
# Multi-block cylinder mesh cell mapping
# ─────────────────────────────────────────────
# blockMeshDict: 5 blocks with K (axial) as the SLOWEST index within each block.
#   Block 0 (center): 2×2×40 = 160 cells, global start 0
#   Block 1 (lower) : 2×1×40 =  80 cells, global start 160
#   Block 2 (right) : 2×1×40 =  80 cells, global start 240
#   Block 3 (upper) : 2×1×40 =  80 cells, global start 320
#   Block 4 (left)  : 2×1×40 =  80 cells, global start 400
#
# Cells at axial position k share the same x coordinate:
#   center: 4k, 4k+1, 4k+2, 4k+3
#   lower : 160+2k, 160+2k+1
#   right : 240+2k, 240+2k+1
#   upper : 320+2k, 320+2k+1
#   left  : 400+2k, 400+2k+1

_CELLS_AT_AXIAL: list[np.ndarray] = [
    np.array([
        4*k, 4*k+1, 4*k+2, 4*k+3,
        160+2*k, 160+2*k+1,
        240+2*k, 240+2*k+1,
        320+2*k, 320+2*k+1,
        400+2*k, 400+2*k+1,
    ])
    for k in range(N_AXIAL)
]

# Inverse mapping: cell index → axial position
_CELL_TO_AXIAL = np.empty(N_CELLS, dtype=int)
for _k, _cidx in enumerate(_CELLS_AT_AXIAL):
    _CELL_TO_AXIAL[_cidx] = _k


def axial_cell_index(axial_node: int, section: int = 0) -> int:
    """Return the global cell index for center-block cell at axial_node."""
    return 4 * axial_node + (section % 4)


def axial_positions() -> np.ndarray:
    """軸方向断面中心の x 座標を返す。"""
    return (np.arange(N_AXIAL) + 0.5) * DX


def cell_axial_positions() -> np.ndarray:
    """全セルに対応する軸方向 x 座標を返す。"""
    xpos = axial_positions()
    return xpos[_CELL_TO_AXIAL]


def axial_profile(temperature: np.ndarray) -> np.ndarray:
    """全セル温度を軸方向断面平均温度へ変換する。"""
    if temperature.ndim == 1:
        return np.array([temperature[g].mean() for g in _CELLS_AT_AXIAL])
    return np.stack([temperature[:, g].mean(axis=1) for g in _CELLS_AT_AXIAL], axis=1)


# ─────────────────────────────────────────────
# OI に必要な観測行列・背景誤差共分散
# ─────────────────────────────────────────────
def build_oi_matrices() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optimal Interpolation (OI) に必要な行列 (H, B, R) を作って返す。

    Returns
    -------
    H   : (n_obs, N_CELLS)
    B   : (N_CELLS, N_CELLS)
    R   : (n_obs, n_obs)
    """
    n_obs = len(ASSIM_SENSOR_AXIAL_NODES)
    H = np.zeros((n_obs, N_CELLS))
    for i, node in enumerate(ASSIM_SENSOR_AXIAL_NODES):
        H[i, axial_cell_index(node)] = 1.0

    x = cell_axial_positions()
    distance = np.abs(x[:, None] - x[None, :])
    B = BACKGROUND_ERROR_STD**2 * np.exp(
        -(distance**2) / (2.0 * CORRELATION_LENGTH_M**2)
    )
    R = np.eye(n_obs) * OBS_NOISE_STD**2
    return H, B, R


# ─────────────────────────────────────────────
# 熱入力系列の生成
# ─────────────────────────────────────────────
def make_gradient_series(n_steps: int, g_on: float) -> np.ndarray:
    """ON/OFF サイクルの境界勾配系列 [degC/m] を作る。"""
    cycle = CYCLE_ON + CYCLE_OFF
    return np.array([g_on if (k % cycle) < CYCLE_ON else 0.0 for k in range(n_steps)])


def save_temperature_history(path: Path, temperature: np.ndarray) -> None:
    """温度時系列をCSVで保存する。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    time_s = (np.arange(N_STEPS) + 1) * DT
    data = np.column_stack([time_s, temperature])
    header = "time_s," + ",".join(f"C{i}" for i in range(temperature.shape[1]))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.10g")


def save_axial_temperature_history(path: Path, temperature: np.ndarray) -> None:
    """断面平均温度時系列をCSVで保存する。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    profile = axial_profile(temperature)
    time_s = (np.arange(N_STEPS) + 1) * DT
    data = np.column_stack([time_s, profile])
    header = "time_s," + ",".join(f"N{i}_x{(i + 0.5) * DX:.3f}m" for i in range(N_AXIAL))
    np.savetxt(path, data, delimiter=",", header=header, comments="", fmt="%.10g")


def save_openfoam_time_dirs(case_dir: Path, dest_dir: Path) -> None:
    """
    現在のOpenFOAM時刻ディレクトリを保存する。

    ofi.reset() は次の計算パスの前に時刻ディレクトリを削除するため、
    各パスの直後にコピーして残す。
    constant/ と system/ もコピーし case.foam を作成することで
    ParaView で直接開けるようにする。
    """
    if dest_dir.exists():
        shutil.rmtree(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    for d in case_dir.iterdir():
        if d.is_dir() and _is_numeric_dir(d.name):
            shutil.copytree(d, dest_dir / d.name)

    for subdir in ("constant", "system"):
        src = case_dir / subdir
        if src.exists():
            shutil.copytree(src, dest_dir / subdir)
    (dest_dir / "case.foam").touch()


def save_summary(path: Path, truth_T: np.ndarray, model_only_T: np.ndarray, da_T: np.ndarray) -> None:
    """節点ごとのRMSEと平均RMSEをCSVで保存する。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    t_start = 2 * N_STEPS // 3
    truth_profile = axial_profile(truth_T)
    model_profile = axial_profile(model_only_T)
    da_profile = axial_profile(da_T)
    rmse_of = np.array([
        np.sqrt(np.mean((truth_profile[t_start:, i] - model_profile[t_start:, i]) ** 2))
        for i in range(N_AXIAL)
    ])
    rmse_da = np.array([
        np.sqrt(np.mean((truth_profile[t_start:, i] - da_profile[t_start:, i]) ** 2))
        for i in range(N_AXIAL)
    ])
    data = np.column_stack([np.arange(N_AXIAL), axial_positions(), rmse_of, rmse_da])
    np.savetxt(
        path,
        data,
        delimiter=",",
        header="axial_node,x_m,rmse_of_only,rmse_oi_da",
        comments="",
        fmt=["%d", "%.10g", "%.10g", "%.10g"],
    )


def _is_numeric_dir(name: str) -> bool:
    try:
        float(name)
        return True
    except ValueError:
        return False


# ─────────────────────────────────────────────
# 真値生成: OF を G_true で走らせる
# ─────────────────────────────────────────────
def generate_truth(ofi: OFInterface, g_true_series: np.ndarray) -> np.ndarray:
    """
    laplacianFoam を G_true で N_STEPS 走らせて真の温度場を生成する。

    Returns
    -------
    truth_T : (N_STEPS, N_CELLS)
    """
    print("─── 真値生成 (G_true で OF を走らせる) ───")
    ofi.reset(initial_T=T_INIT)
    x_current = np.ones(N_CELLS) * T_INIT
    truth_T = np.zeros((N_STEPS, N_CELLS))

    t0 = time_mod.time()
    for k in range(N_STEPS):
        ofi.write_field(k * DT, x_current, gradient=g_true_series[k])
        ofi.run_step(k * DT, (k + 1) * DT)
        x_current = ofi.read_field((k + 1) * DT)
        truth_T[k] = x_current
        if (k + 1) % 30 == 0:
            elapsed = time_mod.time() - t0
            print(f"  step {k+1}/{N_STEPS}  T[N0]={axial_profile(x_current)[0]:.2f}  elapsed {elapsed:.0f}s")

    print(f"真値生成完了  所要時間: {time_mod.time()-t0:.1f}s")
    return truth_T


def generate_model_only(ofi: OFInterface, g_model_series: np.ndarray) -> np.ndarray:
    """
    データ同化を使わず、G_model だけで OpenFOAM を走らせた温度場を生成する。

    これは比較用の「DA なし」ケースであり、OI 補正後の温度場は一切使わない。
    """
    print("─── DAなし予測 (G_model で OF を独立に走らせる) ───")
    ofi.reset(initial_T=T_INIT)
    x_current = np.ones(N_CELLS) * T_INIT
    model_T = np.zeros((N_STEPS, N_CELLS))

    t0 = time_mod.time()
    for k in range(N_STEPS):
        ofi.write_field(k * DT, x_current, gradient=g_model_series[k])
        ofi.run_step(k * DT, (k + 1) * DT)
        x_current = ofi.read_field((k + 1) * DT)
        model_T[k] = x_current
        if (k + 1) % 30 == 0:
            elapsed = time_mod.time() - t0
            print(f"  step {k+1}/{N_STEPS}  T[N0]={axial_profile(x_current)[0]:.2f}  elapsed {elapsed:.0f}s")

    print(f"DAなし予測完了  所要時間: {time_mod.time()-t0:.1f}s")
    return model_T


# ─────────────────────────────────────────────
# データ同化ループ
# ─────────────────────────────────────────────
def run_da(
    ofi: OFInterface,
    truth_T: np.ndarray,
    g_model_series: np.ndarray,
    H: np.ndarray,
    B: np.ndarray,
    R: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    laplacianFoam (G_model) + Optimal Interpolation でデータ同化を実行する。

    Returns
    -------
    of_T  : (N_STEPS, N_CELLS)  OI更新直前の OF 予測
    da_T  : (N_STEPS, N_CELLS)  OI 補正後の推定値
    """
    print("─── データ同化ループ (G_model + OI) ───")
    ofi.reset(initial_T=T_INIT)

    n_obs = len(ASSIM_SENSOR_AXIAL_NODES)
    gain = B @ H.T @ np.linalg.inv(H @ B @ H.T + R)

    x_current = np.ones(N_CELLS) * T_INIT
    of_T = np.zeros((N_STEPS, N_CELLS))
    da_T = np.zeros((N_STEPS, N_CELLS))

    t0 = time_mod.time()
    for k in range(N_STEPS):
        # センサ観測（真値にノイズを加える）
        y = H @ truth_T[k] + np.random.randn(n_obs) * OBS_NOISE_STD

        # ① OF 予測ステップ: x_current を IC にして G_model で1ステップ
        ofi.write_field(k * DT, x_current, gradient=g_model_series[k])
        ofi.run_step(k * DT, (k + 1) * DT)
        x_pred = ofi.read_field((k + 1) * DT)

        # ② OI 更新ステップ: y で x_pred を補正
        innovation = y - H @ x_pred
        x_da = x_pred + gain @ innovation

        # 次ステップのICを OI 補正値に設定
        x_current = x_da
        of_T[k] = x_pred
        da_T[k] = x_da

        if (k + 1) % 30 == 0:
            elapsed = time_mod.time() - t0
            print(
                f"  step {k+1}/{N_STEPS}"
                f"  true={axial_profile(truth_T[k])[0]:.2f}"
                f"  of={axial_profile(x_pred)[0]:.2f}"
                f"  oi={axial_profile(x_da)[0]:.2f}"
                f"  elapsed {elapsed:.0f}s"
            )

    print(f"DAループ完了  所要時間: {time_mod.time()-t0:.1f}s")
    return of_T, da_T


# ─────────────────────────────────────────────
# 結果の表示
# ─────────────────────────────────────────────
def print_summary(truth_T: np.ndarray, of_T: np.ndarray, da_T: np.ndarray) -> None:
    """定常域（後半 1/3）の RMSE をターミナルに表示する。"""
    t_start = 2 * N_STEPS // 3
    truth_profile = axial_profile(truth_T)
    of_profile = axial_profile(of_T)
    da_profile = axial_profile(da_T)
    print()
    print("=" * 68)
    print("  ツイン実験 RMSE（定常域）  [対 truth: G_true で OF を走らせた真値]")
    print("=" * 68)
    header = f"  {'case':<28}" + "".join(f" N{i:>3}" for i in range(N_AXIAL)) + "  avg"
    print(header)
    print("-" * 68)
    for label, arr in [
        ("OF-only (DA なし)    ", of_profile),
        ("OI あり (DA あり)    ", da_profile),
    ]:
        rmse = [
            float(np.sqrt(np.mean((truth_profile[t_start:, i] - arr[t_start:, i]) ** 2)))
            for i in range(N_AXIAL)
        ]
        vals = "".join(f" {v:>4.2f}" for v in rmse)
        print(f"  {label}{vals}  {np.mean(rmse):.3f}")
    avg_of = float(np.sqrt(np.mean((truth_profile[t_start:] - of_profile[t_start:]) ** 2)))
    avg_da = float(np.sqrt(np.mean((truth_profile[t_start:] - da_profile[t_start:]) ** 2)))
    if avg_of > 0:
        print(f"  平均RMSE改善率: {(1.0 - avg_da / avg_of) * 100:.1f}%")
    print("=" * 68)
    for n in ASSIM_SENSOR_AXIAL_NODES:
        print(f"  ※ Node {n} は同化センサ節点 x={(n + 0.5) * DX:.3f}m")
    for n in VALIDATION_SENSOR_AXIAL_NODES:
        print(f"  ※ Node {n} は検証専用センサ節点 x={(n + 0.5) * DX:.3f}m")
    print()


def plot_results(
    truth_T: np.ndarray,
    of_T: np.ndarray,
    da_T: np.ndarray,
    g_true_series: np.ndarray,
    save_path: str = "results_of_da.png",
) -> None:
    """結果グラフを描画して保存する。"""
    time_min = np.arange(N_STEPS) * DT / 60.0
    t_start = 2 * N_STEPS // 3
    truth_profile = axial_profile(truth_T)
    of_profile = axial_profile(of_T)
    da_profile = axial_profile(da_T)
    hidden = [i for i in range(N_AXIAL) if i not in ASSIM_SENSOR_AXIAL_NODES]
    hidden_to_plot = VALIDATION_SENSOR_AXIAL_NODES + [
        i for i in hidden if i not in VALIDATION_SENSOR_AXIAL_NODES
    ][:2]
    cmap = plt.cm.tab10(np.linspace(0, 1, N_AXIAL))

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(
        "Round Bar Twin Experiment: laplacianFoam + Optimal Interpolation\n"
        f"Truth: G={G_TRUE}degC/m,  Model: G={G_MODEL}degC/m ({G_MODEL/G_TRUE*100:.0f}%),  "
        f"Assim sensors: nodes {ASSIM_SENSOR_AXIAL_NODES}, validation: nodes {VALIDATION_SENSOR_AXIAL_NODES}",
        fontsize=11,
    )

    # (1) センサ節点の温度推定
    ax = axes[0, 0]
    for node in ASSIM_SENSOR_AXIAL_NODES:
        c = cmap[node]
        ax.plot(time_min, truth_profile[:, node], color=c, ls="--", lw=1.8, label=f"N{node} truth")
        ax.plot(time_min, of_profile[:, node], color=c, ls=":", lw=1.3, alpha=0.8, label=f"N{node} OF only")
        ax.plot(time_min, da_profile[:, node], color=c, ls="-", lw=2.0, label=f"N{node} OI (DA)")
    ax.set_title("(1) Assimilation sensors: truth vs OF-only vs OI", fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [degC]")
    ax.legend(fontsize=7.5)
    ax.grid(True, alpha=0.3)

    # (2) 隠れ節点の温度推定（代表4節点）
    ax = axes[0, 1]
    for node in hidden_to_plot:
        c = cmap[node]
        ax.plot(time_min, truth_profile[:, node], color=c, ls="--", lw=1.8, label=f"N{node} truth")
        ax.plot(time_min, of_profile[:, node], color=c, ls=":", lw=1.3, alpha=0.8, label=f"N{node} OF only")
        ax.plot(time_min, da_profile[:, node], color=c, ls="-", lw=2.0, label=f"N{node} OI (DA)")
    node_label = ", ".join(f"N{node}" for node in hidden_to_plot)
    ax.set_title(f"(2) Validation/hidden nodes ({node_label}): truth vs OF-only vs OI", fontsize=10)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Temperature [degC]")
    ax.legend(fontsize=7)
    ax.grid(True, alpha=0.3)

    # (3) 全節点の RMSE 棒グラフ（定常域）
    ax = axes[1, 0]
    rmse_of = [
        float(np.sqrt(np.mean((truth_profile[t_start:, i] - of_profile[t_start:, i]) ** 2)))
        for i in range(N_AXIAL)
    ]
    rmse_da = [
        float(np.sqrt(np.mean((truth_profile[t_start:, i] - da_profile[t_start:, i]) ** 2)))
        for i in range(N_AXIAL)
    ]
    x = np.arange(N_AXIAL)
    b1 = ax.bar(x - 0.2, rmse_of, 0.4, label=f"OF only (G={G_MODEL})", color="steelblue", alpha=0.8, edgecolor="k")
    b2 = ax.bar(x + 0.2, rmse_da, 0.4, label="OF + OI (DA)", color="crimson", alpha=0.8, edgecolor="k")
    ax.bar_label(b1, fmt="%.2f", fontsize=7, padding=2)
    ax.bar_label(b2, fmt="%.2f", fontsize=7, padding=2)
    ax.axhline(OBS_NOISE_STD, color="gray", ls=":", lw=1.5,
               label=f"Sensor noise = {OBS_NOISE_STD} degC")
    y_label = max(max(rmse_of), max(rmse_da), OBS_NOISE_STD) * 1.1
    for node in ASSIM_SENSOR_AXIAL_NODES:
        ax.text(node, y_label, "assim", ha="center", fontsize=8, color="dimgray")
    for node in VALIDATION_SENSOR_AXIAL_NODES:
        ax.text(node, y_label, "valid", ha="center", fontsize=8, color="dimgray")
    ax.set_xticks(x)
    ax.set_xticklabels([f"N{i}" for i in range(N_AXIAL)], fontsize=9)
    ax.set_ylabel("RMSE [degC]")
    ax.set_title("(3) Temperature RMSE (steady-state region)", fontsize=10)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3, axis="y")

    # (4) 隠れ節点の誤差時系列
    ax = axes[1, 1]
    for node in hidden_to_plot:
        c = cmap[node]
        ax.plot(time_min, truth_profile[:, node] - of_profile[:, node],
                color=c, ls=":", lw=1.3, label=f"N{node} OF only")
        ax.plot(time_min, truth_profile[:, node] - da_profile[:, node],
                color=c, ls="-", lw=1.8, label=f"N{node} OI (DA)")
    ax.axhline(0, color="k", lw=0.8, ls="--")
    ax.axvline(time_min[t_start], color="red", lw=1.0, ls="--", alpha=0.5,
               label="Steady-state boundary")
    ax.set_title("(4) Error for hidden nodes (truth - estimate)\ndotted: OF only, solid: OI", fontsize=9)
    ax.set_xlabel("Time [min]")
    ax.set_ylabel("Error [degC]")
    ax.legend(fontsize=7.5, ncol=2)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    print(f"[保存] {save_path}")
    plt.close(fig)


# ─────────────────────────────────────────────
# メイン
# ─────────────────────────────────────────────
if __name__ == "__main__":
    np.random.seed(42)
    print("=" * 60)
    print("  丸棒 laplacianFoam ツイン実験 + OI データ同化デモ")
    print("=" * 60)
    print(f"  rod L={ROD_LENGTH_M}m, D={ROD_DIAMETER_M}m, axial nodes={N_AXIAL}, cells={N_CELLS}, dx={DX:.4f}m")
    print(f"  N_STEPS={N_STEPS}, DT={DT}s, 総時間={N_STEPS*DT}s={N_STEPS*DT/60:.1f}min")
    print(f"  G_true={G_TRUE}, G_model={G_MODEL} ({G_MODEL/G_TRUE*100:.0f}% of truth)")
    print(f"  assim sensors={ASSIM_SENSOR_AXIAL_NODES}, validation sensors={VALIDATION_SENSOR_AXIAL_NODES}")
    print()

    # 熱入力勾配系列を生成
    g_true_series = make_gradient_series(N_STEPS, G_TRUE)
    g_model_series = make_gradient_series(N_STEPS, G_MODEL)

    # OI 行列を準備
    H, B_mat, R_mat = build_oi_matrices()

    # OpenFOAM インターフェース
    ofi = OFInterface(CASE_DIR, n_cells=N_CELLS, dt=DT, t_amb=T_AMB)

    # PASS 1: 真値生成（G_true で OF を走らせる）
    truth_T = generate_truth(ofi, g_true_series)
    save_temperature_history(RESULTS_DIR / "truth" / "temperature_history.csv", truth_T)
    save_axial_temperature_history(RESULTS_DIR / "truth" / "axial_temperature_history.csv", truth_T)
    save_openfoam_time_dirs(CASE_DIR, RESULTS_DIR / "truth" / "openfoam_time_dirs")

    # PASS 2: 比較用の DA なし予測（G_model だけで OF を独立に走らせる）
    print()
    model_only_T = generate_model_only(ofi, g_model_series)
    save_temperature_history(RESULTS_DIR / "model_only" / "temperature_history.csv", model_only_T)
    save_axial_temperature_history(RESULTS_DIR / "model_only" / "axial_temperature_history.csv", model_only_T)
    save_openfoam_time_dirs(CASE_DIR, RESULTS_DIR / "model_only" / "openfoam_time_dirs")

    # PASS 3: データ同化（G_model + OI）
    print()
    da_forecast_T, da_T = run_da(ofi, truth_T, g_model_series, H, B_mat, R_mat)
    save_temperature_history(RESULTS_DIR / "with_da" / "temperature_history.csv", da_T)
    save_temperature_history(RESULTS_DIR / "with_da" / "openfoam_forecast_before_oi.csv", da_forecast_T)
    save_axial_temperature_history(RESULTS_DIR / "with_da" / "axial_temperature_history.csv", da_T)
    save_axial_temperature_history(RESULTS_DIR / "with_da" / "axial_openfoam_forecast_before_oi.csv", da_forecast_T)
    save_openfoam_time_dirs(CASE_DIR, RESULTS_DIR / "with_da" / "openfoam_time_dirs")

    # 結果の表示と保存
    print_summary(truth_T, model_only_T, da_T)
    save_summary(RESULTS_DIR / "summary_rmse.csv", truth_T, model_only_T, da_T)
    out_path = str(RESULTS_DIR / "results_of_da.png")
    plot_results(truth_T, model_only_T, da_T, g_true_series, save_path=out_path)
    legacy_out_path = str(Path(__file__).parent / "results_of_da.png")
    if legacy_out_path != out_path:
        shutil.copyfile(out_path, legacy_out_path)
