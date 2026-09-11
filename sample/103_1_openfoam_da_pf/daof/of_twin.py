"""OpenFOAM-in-the-loop 双子実験ドライバ(field-space EnKF / PF).

各メンバーは実際に chtMultiRegionFoam(既存メッシュ・輻射あり)を短ウィンドウで回し、
固体温度場(全セル)＋ヒータ発熱 Q を状態として EnKF/PF で更新する。

流れ:
    1. 真値ラン(Q=Q_true, 初期一様 293.15 K)を各観測時刻まで回し、真の固体温度場を得る
    2. 真値場の観測セル値にノイズを乗せて合成観測を作る
    3. でたらめな初期状態(各メンバー一様ランダム温度＋ランダム Q)から出発し、
       サイクルごとに [全メンバーを前進積分 → 観測で解析更新 → 場を書き戻す]
    4. アンサンブル平均場と真値場の RMSE 等を記録

解析(EnKF/PF)のコアは ROM 版と共通(dacore.enkf / dacore.pf)。モデル非依存。
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore.enkf import enkf_update           # noqa: E402
from dacore.pf import pf_update               # noqa: E402
from daof import of_case                      # noqa: E402

K = 273.15


def _cycle_times(cfg):
    t_end = cfg["experiment"]["t_end_s"]
    dt = cfg["experiment"]["obs_interval_s"]
    return list(np.round(np.arange(dt, t_end + 1e-9, dt), 6))


def _log(msg):
    print(msg, flush=True)


def run_openfoam_twin(cfg, filter_kind, workdir):
    os.makedirs(workdir, exist_ok=True)
    fields_dir = os.path.join(workdir, "fields")
    os.makedirs(fields_dir, exist_ok=True)

    centres = of_case.solid_cell_centres()
    Nc = len(centres)
    np.save(os.path.join(fields_dir, "cell_centres.npy"), centres)

    probes = cfg["observation"]["probes"]
    pnames = list(probes)
    obs_cells = of_case.nearest_cells([probes[n] for n in pnames], centres)
    n_obs = len(obs_cells)
    n_aug = Nc + 1
    iQ = Nc

    H = np.zeros((n_obs, n_aug))
    H[np.arange(n_obs), obs_cells] = 1.0
    noise = cfg["observation"]["noise_C"]
    R = np.eye(n_obs) * noise ** 2

    seed = cfg["experiment"]["seed"]
    rng_obs = np.random.default_rng(seed)
    rng = np.random.default_rng(seed + 1)

    fb = cfg["filter"]
    Tclip = tuple(fb["T_clip_K"])
    Qb = tuple(fb["Q_bounds_W"])
    times = _cycle_times(cfg)

    # ---------------- 真値ラン ----------------
    _log(f"[twin] truth run over cycles {times} (Q_true={cfg['truth']['Q_true_W']} W)")
    truth_dir = os.path.join(workdir, "truth")
    of_case.prepare_member(truth_dir)
    of_case.set_solid_state(truth_dir, 0.0, cfg["truth"]["T0_K"],
                            cfg["truth"]["Q_true_W"])
    truth_fields = {0.0: np.full(Nc, cfg["truth"]["T0_K"])}
    obs_values = {}
    t0 = 0.0
    for t1 in times:
        ok, log = of_case.run_window(truth_dir, t0, t1, logname=f"log.truth_{t0:g}_{t1:g}")
        if not ok:
            raise RuntimeError(f"truth run failed [{t0},{t1}] see {log}")
        Tf = of_case.read_solid_T(truth_dir, t1)
        truth_fields[t1] = Tf
        obs_values[t1] = Tf[obs_cells] + rng_obs.normal(0.0, noise, n_obs)
        _log(f"[twin]  truth t={t1:g}s  Tmax={Tf.max()-K:.3f}C  "
             f"obs={np.round(obs_values[t1]-K,3)}")
        t0 = t1
    np.save(os.path.join(fields_dir, "truth_final.npy"), truth_fields[times[-1]])

    # ---------------- でたらめな初期アンサンブル ----------------
    N = cfg["ensemble"]["n_members"]
    ec = cfg["ensemble"]
    T0 = rng.uniform(ec["init_T_min_K"], ec["init_T_max_K"], N)
    Q = rng.uniform(ec["init_Q_min_W"], ec["init_Q_max_W"], N)
    members = []
    for i in range(N):
        m = os.path.join(workdir, f"member_{i:02d}")
        of_case.prepare_member(m)
        of_case.set_solid_state(m, 0.0, float(T0[i]), float(Q[i]))
        members.append(m)
    _log(f"[twin] {N} members init: T0={np.round(T0-K,1)}C  Q={np.round(Q,1)}W")

    def rmse(a, b):
        return float(np.sqrt(np.mean((a - b) ** 2)))

    # t=0 の記録(でたらめ状態 vs 真値 293.15)
    mean0 = np.full(Nc, T0.mean())
    hist = {
        "filter": filter_kind, "Nc": Nc, "obs_cells": obs_cells.tolist(),
        "Q_true": cfg["truth"]["Q_true_W"],
        "times": [0.0],
        "rmse_field": [rmse(mean0, truth_fields[0.0])],
        "rmse_obs": [rmse(mean0[obs_cells], truth_fields[0.0][obs_cells])],
        "q_mean": [float(Q.mean())], "q_std": [float(Q.std())],
        "ess": [float(N)],
    }
    np.save(os.path.join(fields_dir, "ensmean_t0.npy"), mean0)

    # ---------------- 同化サイクル ----------------
    t0 = 0.0
    for ci, t1 in enumerate(times):
        _log(f"[twin] === cycle {ci+1}/{len(times)}: forecast [{t0:g},{t1:g}] "
             f"for {N} members (serial) ===")
        for i, m in enumerate(members):
            ok, log = of_case.run_window(m, t0, t1)
            if not ok:
                raise RuntimeError(f"member {i} failed [{t0},{t1}] see {log}")
            _log(f"[twin]   member {i} done [{t0:g},{t1:g}]")

        # 予報アンサンブルを集める
        Z = np.empty((N, n_aug))
        for i, m in enumerate(members):
            Z[i, :Nc] = of_case.read_solid_T(m, t1)
            Z[i, iQ] = Q[i]
        # 解析前の予報(補正前)平均場を保存(イノベーション可視化用)
        np.save(os.path.join(fields_dir, f"foremean_t{t1:g}.npy"),
                Z[:, :Nc].mean(axis=0))

        y = obs_values[t1]

        # 解析更新
        if filter_kind == "enkf":
            Za = enkf_update(Z, y, H, R, rng, inflation=fb["inflation"])
            ess = np.nan
        elif filter_kind == "pf":
            Za, _w, ess, resampled = pf_update(
                Z, y, H, R, rng, ess_frac=fb["pf_ess_frac"],
                jitter_T=fb["pf_jitter_T_K"], n_nodes=Nc)
            if resampled:
                Za[:, iQ] += rng.normal(0.0, fb["pf_jitter_Q_W"], N)
        else:
            raise ValueError(filter_kind)

        Za[:, :Nc] = np.clip(Za[:, :Nc], *Tclip)
        Za[:, iQ] = np.clip(Za[:, iQ], *Qb)

        # 解析場を t1 に書き戻す(次ウィンドウはここから再開)
        for i, m in enumerate(members):
            of_case.set_solid_state(m, t1, Za[i, :Nc], float(Za[i, iQ]))
        Q = Za[:, iQ].copy()

        ens_mean = Za[:, :Nc].mean(axis=0)
        np.save(os.path.join(fields_dir, f"ensmean_t{t1:g}.npy"), ens_mean)
        hist["times"].append(float(t1))
        hist["rmse_field"].append(rmse(ens_mean, truth_fields[t1]))
        hist["rmse_obs"].append(rmse(ens_mean[obs_cells], truth_fields[t1][obs_cells]))
        hist["q_mean"].append(float(Q.mean()))
        hist["q_std"].append(float(Q.std()))
        hist["ess"].append(float(ess))
        _log(f"[twin]  analysis t={t1:g}s  RMSE(field)={hist['rmse_field'][-1]:.3f}K  "
             f"Q={Q.mean():.2f}±{Q.std():.2f}W (true {hist['Q_true']})")
        t0 = t1

    _log(f"[twin] DONE. field RMSE: {hist['rmse_field'][0]:.2f} -> "
         f"{hist['rmse_field'][-1]:.3f} K")
    return hist
