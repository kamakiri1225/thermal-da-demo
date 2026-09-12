"""FrontISTR 変位観測演算子: OpenFOAM の固体温度場 → 熱膨張変位(観測点).

104 の観測演算子 h(x) は「メンバーの固体温度場を FrontISTR の線形静解析に通して
上面変位を取り出す」という物理チェーンそのもの。102_1 の連成コード
(cylinder_mesh / fistr_case / openfoam_temperature / run_thermal_expansion)を
そのまま再利用する。1回の評価は FrontISTR 込みで数秒。

観測量(実験のダイヤルゲージ/変位センサに対応):
    uz_heater_mm   : ヒータ側上面ノード列の平均 Uz [mm]
    uz_opposite_mm : 反ヒータ側上面ノード列の平均 Uz [mm]
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
# 102_1 の連成モジュールを利用(隣に clone されている前提。docs 参照)
F102 = os.path.abspath(os.path.join(
    ROOT, "..", "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
if F102 not in sys.path:
    sys.path.insert(0, F102)

import yaml  # noqa: E402

import run_thermal_expansion as rte  # noqa: E402  (102_1)

OBS_NAMES = ["uz_heater_mm", "uz_opposite_mm"]

# 材料物性(102_1 と同じ鋼材)を dict として読む
MATERIAL = yaml.safe_load(Path(str(rte.DEFAULT_MATERIAL)).read_text(encoding="utf-8"))


def displacement_obs(of_case_dir, time_name, work_dir):
    """OpenFOAM ケースの時刻 time_name の固体温度場から変位観測 (2,) を返す [mm].

    of_case_dir : メンバー(または真値)の OpenFOAM ケースディレクトリ
    time_name   : 時刻ディレクトリ名(例 "10")
    work_dir    : FrontISTR の作業ディレクトリ(呼び出しごとに上書き)
    """
    row = rte.run_one_time(
        of_case=Path(of_case_dir),
        time_dir=str(time_name),
        case_dir=Path(work_dir),
        mat=MATERIAL,
    )
    return np.array([float(row["heater_side_top_Uz_mm"]),
                     float(row["opposite_side_top_Uz_mm"])])


if __name__ == "__main__":
    # スモークテスト: 103_0 の真値ラン t=10 で変位を評価
    of_case = os.path.join(ROOT, "..", "103_0_openfoam_da_enkf",
                           "openfoam", "run_enkf", "truth")
    work = os.path.join(ROOT, "openfoam", "fem_smoke")
    u = displacement_obs(of_case, "10", work)
    print(f"[fem-smoke] uz_heater={u[0]*1000:.3f} um  uz_opposite={u[1]*1000:.3f} um")
