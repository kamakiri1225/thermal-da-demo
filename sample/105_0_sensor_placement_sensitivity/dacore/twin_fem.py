"""ROM 版 双子実験: 0→600s の全ストーリー(assim_disp=True で変位もデータ同化).

既定(assim_disp=False)は温度2点のみ同化し、変位は補正後温度からの予測。
assim_disp=True にすると、IDW→FrontISTR で構築した演算子 M(µm/K,
config/M_frontistr_operator.npy) によるアフィン観測 uz=M(T-T_ref) を
温度と一緒に同化する(温度のみ0.076K→温度+変位0.004K、docs/11参照)。
注意: 校正線形写像Dでの同化や、絶対温度 M@T での評価は発散する(過去のバグ)。

さらに「同化なし(free run)」も同じ初期アンサンブルで走らせ、
温度・変位の時刻歴を [真値 / 同化あり / 同化なし] で比較する。
ROM なので 0→600s(ヒータON 300s→OFF)が数秒で完走する。
"""

from __future__ import annotations

import numpy as np

from .cht_rom import N_NODES, T_AIR_K
from .displacement import displacement, load_operator
from .enkf import enkf_update
from .ensemble import I_H, I_QSCALE, N_AUG, clip_params, forecast, init_ensemble
from .observations import (
    generate_truth,
    obs_node_indices,
    true_params,
)


def build_H_R(cfg):
    """温度観測 H (n_t,7) と観測誤差 R を作る(温度のみ同化)."""
    tnodes = obs_node_indices(cfg)                 # 例 [hot, cold]
    n_t = len(tnodes)
    H = np.zeros((n_t, N_AUG))
    for r, c in enumerate(tnodes):
        H[r, c] = 1.0
    sigT = cfg["observation"]["noise_C"]
    R = np.diag([sigT**2]*n_t)
    return H, R, n_t


def _disp_traj(traj_T, op):
    """ノード温度時刻歴 (nt,5) から変位時刻歴 (nt,2)[mm]."""
    return displacement(traj_T, op)


def run_rom_fem_twin(cfg, calib, op, assim_disp=False, M_um=None, sigU_um=0.3):
    """戻り値 hist: 温度(5)・変位(2)の [真値/同化あり/同化なし] 時刻歴."""
    seed = cfg["experiment"]["seed"]
    rng_obs = np.random.default_rng(seed)
    rng_da = np.random.default_rng(seed + 1)
    rng_free = np.random.default_rng(seed + 1)     # 同一の初期アンサンブル

    H, R, n_t = build_H_R(cfg)
    tnodes = obs_node_indices(cfg)
    if assim_disp and M_um is None:
        import os as _os
        M_um = np.load(_os.path.join(_os.path.dirname(__file__), "..", "config",
                                     "M_frontistr_operator.npy")) * 1e6  # µm/K
    n_d = np.asarray(op["D"]).shape[0]
    sigU = cfg["observation"].get("disp_noise_mm", 1e-4)
    fb = cfg["filter"]

    times, truth = generate_truth(calib, cfg)      # (nt,), (nt,5)
    dt_obs = cfg["experiment"]["obs_interval_s"]
    cyc_t = np.arange(dt_obs, cfg["experiment"]["t_end_s"] + 1e-9, dt_obs)
    cyc_idx = np.clip(np.searchsorted(times, cyc_t), 0, len(times) - 1)

    # 合成観測: 同化に使うのは温度2点。変位は表示用に真値＋ノイズを別途作る。
    truth_disp = _disp_traj(truth, op)
    obs = {}
    obs_U = {}
    for k, ci in enumerate(cyc_idx):
        obs[k] = truth[ci][tnodes] + rng_obs.normal(0.0, np.sqrt(np.diag(R)))
        obs_U[k] = truth_disp[ci] + rng_obs.normal(0.0, sigU, n_d)

    Zda = init_ensemble(cfg, rng_da)
    Zfree = init_ensemble(cfg, rng_free)           # 同じでたらめ初期状態

    rec_t = [0.0]
    da_T = [Zda[:, :N_NODES].mean(0)]
    free_T = [Zfree[:, :N_NODES].mean(0)]
    tr_T = [truth[0]]

    t_prev = 0.0
    for k, (t1, ci) in enumerate(zip(cyc_t, cyc_idx)):
        Zda = forecast(Zda, t_prev, t1, calib, cfg, rng_da)
        Zfree = forecast(Zfree, t_prev, t1, calib, cfg, rng_free)  # 同化しない
        if assim_disp:
            # 変位も同化: 観測=温度2点+変位2点(アフィン M(T-Tref) [µm])
            ci_t = cyc_idx[k]
            du_true = M_um @ (truth[ci_t] - T_AIR_K)
            y4 = np.concatenate([obs[k], du_true + rng_obs.normal(0, sigU_um, len(du_true))])
            R4 = np.diag(list(np.diag(R)) + [sigU_um**2]*len(du_true))
            Yf = np.zeros((len(Zda), len(y4)))
            Yf[:, :n_t] = Zda[:, tnodes]
            Yf[:, n_t:] = (Zda[:, :N_NODES] - T_AIR_K) @ M_um.T
            Zda = enkf_update(Zda, y4, None, R4, rng_da, inflation=fb["inflation"], Yf=Yf)
        else:
            Zda = enkf_update(Zda, obs[k], H, R, rng_da, inflation=fb["inflation"])
        clip_params(Zda, cfg)
        rec_t.append(float(t1))
        da_T.append(Zda[:, :N_NODES].mean(0))
        free_T.append(Zfree[:, :N_NODES].mean(0))
        tr_T.append(truth[ci])
        t_prev = t1

    rec_t = np.array(rec_t)
    da_T = np.array(da_T); free_T = np.array(free_T); tr_T = np.array(tr_T)
    hist = {
        "times": rec_t,
        "truth_T": tr_T, "da_T": da_T, "free_T": free_T,
        "truth_U": _disp_traj(tr_T, op),
        "da_U": _disp_traj(da_T, op),
        "free_U": _disp_traj(free_T, op),
        "obs_times": cyc_t,
        "obs_T": np.array([obs[k] for k in range(len(cyc_t))]),
        "obs_U": np.array([obs_U[k] for k in range(len(cyc_t))]),
        "tnodes": tnodes.tolist(),
        "disp_names": op["disp_names"],
    }
    return hist
