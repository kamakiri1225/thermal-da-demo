"""OpenFOAM×FrontISTR 双子実験ドライバ: 温度＋変位観測の field-space EnKF.

103 の of_twin.py を拡張し、観測ベクトルを

    y = [ T_hot, T_cold,  uz_heater_mm, uz_opposite_mm ]

の4成分にする。変位は「固体温度場 → FrontISTR 線形静解析 → 上面 Uz」という
非線形観測演算子なので、H 行列ではなく **各メンバーごとに実際に FrontISTR を
回して予報観測 Yf をサンプリング**し、enkf_update(Yf=...) に渡す。

状態は 103 と同じ z = [固体温度場(20696セル), Q]。書き戻しも同じ。
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore.enkf import enkf_update       # noqa: E402
from daof import of_case                  # noqa: E402
from fem.fem_obs import displacement_obs  # noqa: E402

K = 273.15


def _log(msg):
    print(msg, flush=True)


def _tname(t):
    return str(int(t)) if float(t) == int(t) else ("%g" % float(t))


def run_fem_twin(cfg, workdir):
    os.makedirs(workdir, exist_ok=True)
    fields_dir = os.path.join(workdir, "fields")
    os.makedirs(fields_dir, exist_ok=True)

    centres = of_case.solid_cell_centres()
    Nc = len(centres)
    np.save(os.path.join(fields_dir, "cell_centres.npy"), centres)

    oc = cfg["observation"]
    tnames = list(oc["temp_probes"])
    obs_cells = of_case.nearest_cells([oc["temp_probes"][n] for n in tnames], centres)
    n_t = len(obs_cells)                     # 温度観測数
    n_d = len(oc["disp_names"])              # 変位観測数
    n_obs = n_t + n_d
    n_aug = Nc + 1
    iQ = Nc

    R = np.diag([oc["temp_noise_C"] ** 2] * n_t +
                [oc["disp_noise_mm"] ** 2] * n_d)

    seed = cfg["experiment"]["seed"]
    rng_obs = np.random.default_rng(seed)
    rng = np.random.default_rng(seed + 1)
    fb = cfg["filter"]

    dt_obs = cfg["experiment"]["obs_interval_s"]
    times = list(np.round(np.arange(dt_obs, cfg["experiment"]["t_end_s"] + 1e-9,
                                    dt_obs), 6))

    def member_obs(case_dir, t, fem_work):
        """メンバー/真値の観測ベクトル [T_hot,T_cold,uz_heater,uz_opposite] を評価."""
        T = of_case.read_solid_T(case_dir, t)
        u = displacement_obs(case_dir, _tname(t), fem_work)
        return np.concatenate([T[obs_cells], u]), T

    # ---------------- 真値ラン ----------------
    _log(f"[twin] truth run over cycles {times} (Q_true={cfg['truth']['Q_true_W']} W)")
    truth_dir = os.path.join(workdir, "truth")
    of_case.prepare_member(truth_dir)
    of_case.set_solid_state(truth_dir, 0.0, cfg["truth"]["T0_K"],
                            cfg["truth"]["Q_true_W"])
    truth_fields = {0.0: np.full(Nc, cfg["truth"]["T0_K"])}
    obs_values, truth_obs_clean = {}, {}
    t0 = 0.0
    for t1 in times:
        ok, log = of_case.run_window(truth_dir, t0, t1,
                                     logname=f"log.truth_{t0:g}_{t1:g}")
        if not ok:
            raise RuntimeError(f"truth run failed [{t0},{t1}] see {log}")
        clean, Tf = member_obs(truth_dir, t1,
                               os.path.join(truth_dir, f"fem_t{t1:g}"))
        truth_fields[t1] = Tf
        noise_vec = np.sqrt(np.diag(R))
        obs_values[t1] = clean + rng_obs.normal(0.0, noise_vec)
        truth_obs_clean[t1] = clean
        _log(f"[twin]  truth t={t1:g}s Tmax={Tf.max()-K:.3f}C "
             f"obs(T)={np.round(obs_values[t1][:n_t]-K,3)} "
             f"obs(uz_um)={np.round(obs_values[t1][n_t:]*1000,3)}")
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

    mean0 = np.full(Nc, T0.mean())
    hist = {
        "filter": "enkf+fem", "Nc": Nc, "obs_cells": obs_cells.tolist(),
        "Q_true": cfg["truth"]["Q_true_W"],
        "times": [0.0],
        "rmse_field": [rmse(mean0, truth_fields[0.0])],
        "q_mean": [float(Q.mean())], "q_std": [float(Q.std())],
        "uz_truth": [np.zeros(n_d)],
        "uz_mean": [np.zeros(n_d)],
        "uz_std": [np.zeros(n_d)],
    }
    np.save(os.path.join(fields_dir, "ensmean_t0.npy"), mean0)

    # ---------------- 同化サイクル ----------------
    t0 = 0.0
    for ci, t1 in enumerate(times):
        _log(f"[twin] === cycle {ci+1}/{len(times)}: forecast [{t0:g},{t1:g}] ===")
        Z = np.empty((N, n_aug))
        Yf = np.empty((N, n_obs))
        for i, m in enumerate(members):
            ok, log = of_case.run_window(m, t0, t1)
            if not ok:
                raise RuntimeError(f"member {i} failed [{t0},{t1}] see {log}")
            yv, Tm = member_obs(m, t1, os.path.join(m, f"fem_t{t1:g}"))
            Z[i, :Nc] = Tm
            Z[i, iQ] = Q[i]
            Yf[i] = yv
            _log(f"[twin]   member {i} done: T@hot={Tm[obs_cells[0]]-K:.2f}C "
                 f"uz_heater={yv[n_t]*1000:.2f}um")

        y = obs_values[t1]
        Za = enkf_update(Z, y, None, R, rng, inflation=fb["inflation"], Yf=Yf)
        Za[:, :Nc] = np.clip(Za[:, :Nc], *fb["T_clip_K"])
        Za[:, iQ] = np.clip(Za[:, iQ], *fb["Q_bounds_W"])

        for i, m in enumerate(members):
            of_case.set_solid_state(m, t1, Za[i, :Nc], float(Za[i, iQ]))
        Q = Za[:, iQ].copy()

        ens_mean = Za[:, :Nc].mean(axis=0)
        np.save(os.path.join(fields_dir, f"ensmean_t{t1:g}.npy"), ens_mean)
        hist["times"].append(float(t1))
        hist["rmse_field"].append(rmse(ens_mean, truth_fields[t1]))
        hist["q_mean"].append(float(Q.mean()))
        hist["q_std"].append(float(Q.std()))
        hist["uz_truth"].append(truth_obs_clean[t1][n_t:].copy())
        hist["uz_mean"].append(Yf[:, n_t:].mean(axis=0))
        hist["uz_std"].append(Yf[:, n_t:].std(axis=0))
        _log(f"[twin]  analysis t={t1:g}s RMSE(field)={hist['rmse_field'][-1]:.3f}K "
             f"Q={Q.mean():.2f}±{Q.std():.2f}W (true {hist['Q_true']})")
        t0 = t1

    _log(f"[twin] DONE. field RMSE: {hist['rmse_field'][0]:.2f} -> "
         f"{hist['rmse_field'][-1]:.3f} K")
    for k in ["uz_truth", "uz_mean", "uz_std"]:
        hist[k] = np.array(hist[k])
    return hist
