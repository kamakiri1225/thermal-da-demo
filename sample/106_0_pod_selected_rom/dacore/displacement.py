"""ROM 節点温度 → 上面変位(2点)の線形観測オペレータ.

FrontISTR の熱膨張は線形静解析なので、上面変位は温度場の線形汎関数
(温度上昇の重み付き和)。ROM では固体を5ノードに縮約しているので、変位を

    uz = D (T_nodes - T_ref)        D: (2, 5)

という線形写像で近似する。係数 D は 102_1(FrontISTR)の変位時刻歴と、
102_0 に校正済み ROM 真値ランのノード温度を突き合わせて最小二乗で同定する。
こうすると ROM の変位が本物の FrontISTR の大きさと一致する。
"""

from __future__ import annotations

import os

import numpy as np
import yaml

from .cht_rom import N_NODES, T_AIR_K, integrate_single
from .observations import true_params

HERE = os.path.dirname(os.path.abspath(__file__))
CONFIG = os.path.join(HERE, "..", "config")
OUT = os.path.join(CONFIG, "displacement_operator.yaml")
# 102_1 FrontISTR の変位時刻歴
FEM_CSV = os.path.join(
    HERE, "..", "..", "102_1_frontistr_hollow_cylinder_thermal_expansion",
    "data", "timehistory.csv")
DISP_NAMES = ["uz_heater_mm", "uz_opposite_mm"]


def calibrate_displacement(calib, verbose=True):
    """D (2,5) を同定して displacement_operator.yaml に保存し、辞書を返す."""
    # 一部の列に "[a, b, c]" 形式(カンマ入り)があり素朴な CSV パースが壊れるため
    # pandas で必要列だけ取り出す
    import pandas as pd
    df = pd.read_csv(FEM_CSV)
    t_fem = df["of_time_s"].to_numpy(dtype=float)
    uz = np.column_stack([df["heater_side_top_Uz_mm"].to_numpy(dtype=float),
                          df["opposite_side_top_Uz_mm"].to_numpy(dtype=float)])

    # ROM 真値ランのノード温度を FEM サンプル時刻で取得
    p = true_params(calib)
    ts, traj = integrate_single(np.full(N_NODES, T_AIR_K), p, 0.0,
                                float(t_fem.max()), dt=1.0)
    idx = np.clip(np.searchsorted(ts, t_fem), 0, len(ts) - 1)
    dT = traj[idx] - T_AIR_K                    # (nt, 5) 温度上昇

    # 最小二乗 uz ≈ dT @ D^T  (切片なし: 一様 T_ref で変位ゼロ)
    D, *_ = np.linalg.lstsq(dT, uz, rcond=None)  # (5, 2)
    D = D.T                                       # (2, 5)
    pred = dT @ D.T
    rmse_um = float(np.sqrt(np.mean((pred - uz) ** 2)) * 1000.0)

    result = {
        "D": [[float(x) for x in row] for row in D],
        "T_ref_K": float(T_AIR_K),
        "disp_names": DISP_NAMES,
        "fit_rmse_um": rmse_um,
        "note": ("ROM節点温度→上面変位[mm]の線形写像 uz=D(T-T_ref)。"
                 "102_1 FrontISTR の変位時刻歴に最小二乗で同定。"),
    }
    with open(OUT, "w") as f:
        yaml.safe_dump(result, f, sort_keys=False, allow_unicode=True)
    if verbose:
        print(f"[disp] D fit RMSE = {rmse_um:.4f} um")
        print(f"[disp] D =\n{np.array2string(D, precision=4)}")
        print(f"[disp] wrote {os.path.relpath(OUT)}")
    return result


def load_operator():
    with open(OUT) as f:
        return yaml.safe_load(f)


def displacement(T_nodes, op):
    """ノード温度 (…,5)[K] から変位 (…,2)[mm] を返す."""
    D = np.asarray(op["D"])
    return (np.asarray(T_nodes) - op["T_ref_K"]) @ D.T


if __name__ == "__main__":
    from .calibrate import load_calibrated
    calibrate_displacement(load_calibrated())
