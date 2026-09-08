"""ROM 版 双子実験: 0→600s の分かりやすい全ストーリー(温度同化＋変位の予測).

温度2点を同化して固体温度場(5ノード)を真値へ補正し、補正後の温度から
FrontISTR 校正済みの線形写像 uz=D(T-T_ref) で上面変位を予測する。
「温度を直せば変位も自動的に真値に一致する」ことを見せる。

なぜ変位を「同化」しないか(重要な教訓):
    縮約した5ノードROMでは、上面変位は一様加熱モードにほぼ比例し、既に観測している
    温度2点とほとんど一次従属(collinear)になる。そのため変位を観測に足しても
    新しい情報がほぼ増えず、逆に条件数が悪化して推定を乱す。変位が真に効くのは、
    全セル温度場を積分する OpenFOAM 版(run_openfoam_fem_enkf.py)の方。
    → ROM: 温度を同化 / 変位は予測、 OpenFOAM: 変位も同化。

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


def run_rom_fem_twin(cfg, calib, op):
    """戻り値 hist: 温度(5)・変位(2)の [真値/同化あり/同化なし] 時刻歴."""
    seed = cfg["experiment"]["seed"]
    rng_obs = np.random.default_rng(seed)
    rng_da = np.random.default_rng(seed + 1)
    rng_free = np.random.default_rng(seed + 1)     # 同一の初期アンサンブル

    H, R, n_t = build_H_R(cfg)
    tnodes = obs_node_indices(cfg)
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
