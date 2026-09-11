"""同化なしフリーラン: でたらめな初期状態のまま OpenFOAM を回す(比較用).

EnKF 実行(run_openfoam_enkf.py)と同じ seed で同じ「でたらめ初期アンサンブル」
(温度・Q とも同一)を作り、**解析更新なし**で 0→t_end を回す。
データ同化あり/なしの温度時刻歴比較(plot_timehistory.py)の材料になる。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 \
      nohup python3 run/run_free_run.py > openfoam/run_free.log 2>&1 &
"""

from __future__ import annotations

import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from daof import of_case

WORKDIR = os.path.join(ROOT, "openfoam", "run_free")


def main():
    with open(os.path.join(ROOT, "openfoam", "da_openfoam_config.yaml")) as f:
        cfg = yaml.safe_load(f)

    seed = cfg["experiment"]["seed"]
    # run_openfoam_twin と同一の乱数列: rng_filt = seed+1, T0 -> Q の順で draw
    rng = np.random.default_rng(seed + 1)
    ec = cfg["ensemble"]
    N = ec["n_members"]
    T0 = rng.uniform(ec["init_T_min_K"], ec["init_T_max_K"], N)
    Q = rng.uniform(ec["init_Q_min_W"], ec["init_Q_max_W"], N)
    t_end = cfg["experiment"]["t_end_s"]

    print(f"[free] no-DA free run: T0={np.round(T0-273.15,1)}C Q={np.round(Q,1)}W "
          f"(EnKF実行と同一のでたらめ初期状態)", flush=True)

    os.makedirs(WORKDIR, exist_ok=True)
    for i in range(N):
        m = os.path.join(WORKDIR, f"member_{i:02d}")
        of_case.prepare_member(m)
        of_case.set_solid_state(m, 0.0, float(T0[i]), float(Q[i]))
        ok, log = of_case.run_window(m, 0.0, t_end, logname="log.free")
        print(f"[free]  member {i} 0->{t_end:g}s ok={ok}", flush=True)
        if not ok:
            raise RuntimeError(f"free run member {i} failed, see {log}")
    print("[free] DONE", flush=True)


if __name__ == "__main__":
    main()
