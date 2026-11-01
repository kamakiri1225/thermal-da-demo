"""ROM 単体（データ同化なし）から gappy-POD で全体の温度分布を戻し、代表点以外の場所で OpenFOAM と比べる。

  1. ROM（代表 5 点の熱回路）を、発熱量とヒータの ON/OFF を OpenFOAM と同じに与えて 0→600 s 解く
     （run/make_rom_applicability_fig.py と同じ計算）
  2. 5 点の温度から全 20,696 セルの温度を戻す：T = 平均場 + Φ Φ_P⁻¹ (T_P − 平均場_P)
  3. OpenFOAM の温度（results/limit_truth_*.npz、30 秒ごと）と、代表 5 点を除いた 20,691 セルで比べる
  比べる場所：代表点以外の 3 セル（上面のヒータ側 T2・上面の反対側 T3 は §11 と同じ、上面の横 T4 は §11 と別）と、全セルの RMSE・最大の差
  参考として「OpenFOAM の温度を 5 つの型で一番よく表したときの差」（型が 5 つしかないことによる限界）も描く。

出力: docs/img/rom_field_vs_openfoam.png, results/rom_field_vs_openfoam.json
再現: OMP_NUM_THREADS=4 python3 run/rom_field_vs_openfoam.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401  （日本語フォント）
from dacore import rom_general as rg
import matplotlib.pyplot as plt
from improve_heating_schedule import integrate, SCHED, QTRUE

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
CASES = [("learned", 1.0, "15 W・0〜300 s\n（ROM を作った条件）"), ("q25", 25 / 15, "25 W・0〜300 s\n（新しい条件）"),
         ("intermittent", 1.0, "15 W 間欠加熱\n（新しい条件）")]
PLACES = [("T2 上面・ヒータ側", (0.028, 0.0, 0.0975), "#C0392B"),
          ("T3 上面・反対側", (-0.028, 0.0, 0.0975), "#1f77b4"),
          ("T4 上面・横（+y 側）", (0.0, 0.028, 0.0975), "#2ca02c")]


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], 5); h = float(d["h"]); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPinv = np.linalg.inv(U[pod])
    others = np.setdiff1d(np.arange(len(mean)), pod)                     # 代表 5 点以外のセル
    cells = [int(np.linalg.norm(Cc - np.array(x), axis=1).argmin()) for _, x, _ in PLACES]
    out = {"note": "ROM 単体（同化なし）→ gappy-POD。代表 5 点以外の 20,691 セルで OpenFOAM と比較。単位 K", "cases": {}}
    res = {}
    for case, q, _ in CASES:
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]; Tf = z["Tfield"].astype(float)
        T = np.full((1, 5), rg.T_AIR_K); tr = [T[0].copy()]
        for k in range(1, len(t)):
            T = integrate(T, C, Km, np.array([h]), np.array([q]), heat, SCHED[case], t[k - 1], t[k]); tr.append(T[0].copy())
        TP = np.array(tr)
        field = mean[None] + ((TP - mean[pod]) @ UPinv.T) @ U.T          # ROM → 全セル
        proj = mean[None] + ((Tf - mean[None]) @ U) @ U.T                # OpenFOAM を 5 つの型で一番よく表した温度
        e = (field - Tf)[:, others]; ep = (proj - Tf)[:, others]
        rmse = np.sqrt((e ** 2).mean(1)); emax = np.abs(e).max(1); rmse_p = np.sqrt((ep ** 2).mean(1))
        res[case] = dict(t=t, truth=Tf[:, cells] - 273.15, rom=field[:, cells] - 273.15, rmse=rmse, emax=emax, rmse_p=rmse_p,
                         rise=float(Tf.max() - rg.T_AIR_K))
        out["cases"][case] = dict(rmse_max=float(rmse.max()), rmse_mean=float(rmse[1:].mean()), abs_err_max=float(emax.max()),
                                  pod5_rmse_max=float(rmse_p.max()), max_rise=res[case]["rise"],
                                  places_abs_err_max={n: float(np.abs(res[case]["rom"][:, i] - res[case]["truth"][:, i]).max())
                                                      for i, (n, _, _) in enumerate(PLACES)})
        print(case, {k: round(v, 3) if isinstance(v, float) else v for k, v in out["cases"][case].items()}, flush=True)

    plt.rcParams.update({"font.size": 16, "axes.titlesize": 19, "axes.labelsize": 18, "xtick.labelsize": 16,
                         "ytick.labelsize": 16, "legend.fontsize": 16})
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.6), sharex=True, gridspec_kw=dict(height_ratios=[1.4, 1], hspace=0.12, wspace=0.42))
    for j, (case, _, lab) in enumerate(CASES):
        r = res[case]; t = r["t"]
        ax = axs[0, j]
        # ヒータの発熱（黄色、右の目盛り）を温度の後ろに重ねる（図 8b と同じ描き方）
        ax2 = ax.twinx(); ax.set_zorder(ax2.get_zorder() + 1); ax.patch.set_visible(False)
        tq = np.linspace(0, 600, 1201); Pw = np.array([QTRUE[case](x) for x in tq], float)
        ax2.fill_between(tq, Pw, step="post", color="#FFF4B8", lw=0); ax2.step(tq, Pw, where="post", color="#D4A800", lw=2.5)
        ax2.set_ylim(0, 40); ax2.set_yticks([0, 15, 25]); ax2.tick_params(axis="y", colors="#9A7B00")
        if j == 2: ax2.set_ylabel("ヒータの発熱 [W]", color="#9A7B00")
        for i, (n, _, c) in enumerate(PLACES):
            ax.plot(t, r["truth"][:, i], "-", color=c, lw=9, alpha=.28, solid_capstyle="round")
            ax.plot(t, r["rom"][:, i], "--", color=c, lw=2.2, label=n)
        ax.set_title(lab, pad=8); ax.grid(alpha=.3); ax.set_ylim(19.5, np.ceil(r["truth"].max() + 0.5))
        ax = axs[1, j]
        ax.plot(t, r["emax"], "-", color="#555", lw=2.2)
        ax.plot(t, r["rmse"], "-", color="#8E44AD", lw=2.8)
        ax.plot(t, r["rmse_p"], ":", color="#8E44AD", lw=2.2)
        ax.set_ylim(0, 0.21); ax.set_yticks([0, 0.05, 0.1, 0.15, 0.2])   # 差が小さいので縦軸を広げて見せる; ax.grid(alpha=.3)
        ax.set_xlabel("時刻 [s]"); ax.set_xticks([0, 200, 400, 600])
        ax.text(0.97, 0.93, f"RMSE 最大 {r['rmse'].max():.2f} K\n最大の差 {r['emax'].max():.2f} K\n（ノイズ 0.3 K より小）", transform=ax.transAxes,
                ha="right", va="top", fontsize=16, bbox=dict(fc="white", ec="#bbb"))
    axs[0, 0].set_ylabel("温度 [°C]"); axs[1, 0].set_ylabel("代表点以外の全セルでの\nROM と OpenFOAM の差 [K]")
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    hs, ls = axs[0, 0].get_legend_handles_labels()
    style = [Line2D([0], [0], color="#555", lw=9, alpha=.28, label="OpenFOAM"), Line2D([0], [0], color="#555", lw=2.2, ls="--", label="ROM から戻した温度"),
             Patch(fc="#FFF4B8", ec="#D4A800", label="ヒータの発熱")]
    err = [Line2D([0], [0], color="#8E44AD", lw=2.8, label="RMSE（代表点以外の全セル）"), Line2D([0], [0], color="#555", lw=2.2, label="最大の差"),
           Line2D([0], [0], color="#8E44AD", lw=2.2, ls=":", label="参考：5 つの型で表せる限界（RMSE）")]
    fig.legend(hs + style, ls + [s.get_label() for s in style], loc="lower center", ncol=6, columnspacing=0.8, handlelength=1.6, bbox_to_anchor=(0.5, 0.05), frameon=False)
    fig.legend(err, [s.get_label() for s in err], loc="lower center", ncol=3, bbox_to_anchor=(0.5, 0.0), frameon=False)
    fig.subplots_adjust(left=.07, right=.93, top=.90, bottom=.17)
    fig.savefig(os.path.join(IMG, "rom_field_vs_openfoam.png"), dpi=130, bbox_inches="tight", facecolor="white"); plt.close(fig)
    json.dump(out, open(os.path.join(RES, "rom_field_vs_openfoam.json"), "w"), ensure_ascii=False, indent=1)
    print("wrote docs/img/rom_field_vs_openfoam.png")


if __name__ == "__main__":
    main()
