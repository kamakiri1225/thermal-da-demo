"""温度時刻歴の比較図: でたらめのまま(同化なし) vs データ同化あり vs 真値.

各ケースの postProcessing/solidTemperatureProbes に残る4点プローブ
(hot/cold=観測点, mid/top=未観測点)の時刻歴を読み、2x2 パネルで

  - 真値ラン(黒破線)
  - 同化なしフリーラン(灰色) … でたらめな初期状態のまま(openfoam/run_free)
  - EnKF メンバー(青)        … 解析更新で真値へジャンプ(openfoam/run_enkf)
  - 観測値(赤点, hot/cold のみ)

を重ねる。解析時刻では書き戻された解析場(<t>/solid/T)から各プローブ最近傍セルの
値を「解析直後の点」として挿入し、直前に NaN を入れてジャンプを線で結ばない。

使い方:
    OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/plot_timehistory.py
出力: docs/img/openfoam_enkf_timehistory.png
"""

from __future__ import annotations

import glob
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from dacore import plots as _p  # noqa: F401  日本語フォント設定
import matplotlib.pyplot as plt
import yaml

from daof import of_case

K = 273.15
# base_case の system/solid/probes と同じ順序・座標
PROBES = {
    "hot":  (0.032, 0.0, 0.05025),
    "cold": (-0.032, 0.0, 0.05025),
    "mid":  (0.0, 0.024, 0.05025),
    "top":  (0.025, 0.0, 0.095),
}
PROBE_NAMES = list(PROBES)
OBSERVED = {"hot", "cold"}


def read_probe_history(case_dir, analysis_times=(), probe_cells=None):
    """postProcessing の全ウィンドウを連結して (times, T[nt,4]) を返す [K].

    analysis_times の各時刻では、書き戻された解析場から4プローブ全ての
    「解析直後の値」を挿入し、直前に NaN を入れて線を切る。
    """
    files = sorted(
        glob.glob(os.path.join(case_dir, "postProcessing",
                               "solidTemperatureProbes", "solid", "*", "T")),
        key=lambda p: float(os.path.basename(os.path.dirname(p))))
    times, rows = [], []
    for f in files:
        win_start = float(os.path.basename(os.path.dirname(f)))
        if probe_cells is not None and any(abs(win_start - a) < 1e-6
                                           for a in analysis_times):
            try:
                Ta = of_case.read_solid_T(case_dir, win_start)[probe_cells]
                times.append(win_start)
                rows.append([np.nan] * 4)          # 解析ジャンプで線を切る
                times.append(win_start + 1e-6)
                rows.append(list(Ta))              # 解析直後の値
            except Exception:
                pass
        for line in open(f):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            v = s.split()
            t = float(v[0])
            if times and t <= times[-1] and not np.isnan(rows[-1][0]):
                continue  # ウィンドウ境界の重複を除去
            times.append(t)
            rows.append([float(x) for x in v[1:5]])
    return np.array(times), np.array(rows)


def main():
    with open(os.path.join(ROOT, "openfoam", "da_openfoam_config.yaml")) as f:
        cfg = yaml.safe_load(f)
    noise = cfg["observation"]["noise_C"]
    seed = cfg["experiment"]["seed"]

    truth_dir = os.path.join(ROOT, "openfoam", "run_enkf", "truth")
    da_members = sorted(glob.glob(os.path.join(ROOT, "openfoam", "run_enkf", "member_*")))
    free_members = sorted(glob.glob(os.path.join(ROOT, "openfoam", "run_free", "member_*")))

    centres = of_case.solid_cell_centres()
    probe_cells = of_case.nearest_cells(list(PROBES.values()), centres)
    analysis_times = [10.0, 20.0]

    t_tr, T_tr = read_probe_history(truth_dir)

    # 観測の再現(run_openfoam_twin と同じ rng_obs=seed, hot/cold の順で draw)
    rng_obs = np.random.default_rng(seed)
    obs_times = analysis_times
    obs_vals = []
    for t in obs_times:
        i = int(np.argmin(np.abs(t_tr - t)))
        clean = np.array([T_tr[i, 0], T_tr[i, 1]])  # hot, cold
        obs_vals.append(clean + rng_obs.normal(0.0, noise, 2))
    obs_vals = np.array(obs_vals)

    # 各ケースの時刻歴を先に読む
    da_hists = [read_probe_history(m, analysis_times, probe_cells)
                for m in da_members]
    free_hists = []
    for m in free_members:
        try:
            tm, Tm = read_probe_history(m)
            if Tm.ndim == 2 and len(tm) >= 2:
                free_hists.append((tm, Tm))
        except Exception:
            pass

    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), sharex=True)
    axes = axes.ravel()
    for pi, name in enumerate(PROBE_NAMES):
        ax = axes[pi]
        lbl = pi == 0
        for j, (tm, Tm) in enumerate(free_hists):
            ax.plot(tm, Tm[:, pi] - K, color="0.65", lw=1.2,
                    label="同化なし(でたらめのまま)" if lbl and j == 0 else None)
        for j, (tm, Tm) in enumerate(da_hists):
            ax.plot(tm, Tm[:, pi] - K, color="tab:blue", lw=1.2, alpha=0.85,
                    label="EnKFメンバー(同化あり)" if lbl and j == 0 else None)
        ax.plot(t_tr, T_tr[:, pi] - K, "k--", lw=2.2,
                label="真値" if lbl else None)
        if name in OBSERVED:
            oi = ["hot", "cold"].index(name)
            ax.scatter(obs_times, obs_vals[:, oi] - K, s=55, color="tab:red",
                       zorder=6, label="観測(ノイズ入り)" if lbl else None)
            tag = "観測点"
        else:
            tag = "未観測点"
        for t in analysis_times:
            ax.axvline(t, color="tab:red", ls=":", lw=0.8, alpha=0.4)
        ax.set_title(f"{name}  ({tag})")
        ax.grid(alpha=0.3)
        if lbl:
            ax.legend(loc="upper right", fontsize=9)
    for ax in axes[2:]:
        ax.set_xlabel("time [s]")
    for ax in (axes[0], axes[2]):
        ax.set_ylabel("温度 [degC]")
    fig.suptitle("温度時刻歴: でたらめな初期状態のまま(灰) vs EnKFデータ同化(青) vs 真値(黒破線)",
                 fontsize=13)
    fig.tight_layout()
    out = os.path.join(ROOT, "docs", "img", "openfoam_enkf_timehistory.png")
    fig.savefig(out, dpi=130)
    print(f"wrote {os.path.relpath(out, ROOT)}  "
          f"(free members plotted: {len(free_hists)}, DA members: {len(da_hists)})")


if __name__ == "__main__":
    main()
