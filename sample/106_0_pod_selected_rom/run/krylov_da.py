"""クリロフ ROM にデータ同化（EnKF）を組み合わせ、熱伝達率のずれを観測で直せるかを確かめる。

  状態：z（クリロフ ROM の座標 r=10 ＋ 一様な温度の成分 1 つ）、発熱の倍率 q、熱伝達率の倍率 s
  ROM：Cr dz/dt + (Kr + s·Hr) z = fr · 15 W · q · on(t)        （Kr = VᵀKV、Hr = VᵀHV、基底 V は s=1 で作り固定）
  観測：温度計 P2・P0（OpenFOAM の値＋ノイズ 0.3 K、30 秒ごと）。今の ROM の「温度 2 点のみ」と同じ条件
  初期のばらつき：温度 20 °C の −3〜+12 K（一様）、q 0.3〜1.8（一様）、s は対数で 0.1〜1.5（一様）
  60 メンバー、膨張係数 1.02、乱数 5 通りの平均。比べるもの：全節点の温度（OpenFOAM を節点へ移した値）

出力：results/krylov_da.json、docs/img/krylov_da.png
再現：OMP_NUM_THREADS=4 python3 run/krylov_da.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import krylov_rom as kr
from dacore import plots as _p  # noqa: F401
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from dacore.temperature_mapping import build_idw_mapping, map_temperature
from improve_heating_schedule import SCHED
import matplotlib.pyplot as plt

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
R_K = 10; N_ENS = 60; SIG = 0.30; INFL = 1.02; DT = 2.0; OBS = 30.0
SEEDS = [20260913, 20260914, 20260915, 20260916, 20260917]
QTRUE = {"learned": 1.0, "q25": 25 / 15, "intermittent": 1.0}


def basis_with_const(fem, r):
    """クリロフ基底に「一様な温度」の成分を足し、C の内積で正規直交にする（初期温度のずれを表せるように）"""
    C = fem["C"]; V = kr.krylov_basis(fem, r)
    v = np.ones(V.shape[0])
    for _ in range(2): v = v - V @ (V.T @ (C @ v))
    return np.column_stack([v / np.sqrt(v @ (C @ v)), V])


def main():
    fem1 = kr.build_fem(1.0); fem0 = kr.build_fem(1e-12)
    Hm = fem1["G"] - fem0["G"]; Km = fem0["G"]; C = fem1["C"]; f = fem1["f"]; X = fem1["X"]
    V = basis_with_const(fem1, R_K); m = V.shape[1]
    Cr = V.T @ (C @ V); Kr = V.T @ (Km @ V); Hr = V.T @ (Hm @ V); fr = V.T @ f
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    nodeP = [int(np.linalg.norm(X - Cc[c], axis=1).argmin()) for c in pod]
    Hobs = V[[nodeP[2], nodeP[0]]]                                   # P2・P0 の温度 = Hobs @ z
    Mmap = build_idw_mapping(Cc, X); one = V.T @ (C @ np.ones(len(X)))  # 一様な温度 1 K の座標

    def forecast(Z, t0, t1, on):
        Z = Z.copy()
        for e in range(len(Z)):
            z = Z[e, :m]; q = Z[e, m]; s = Z[e, m + 1]; G = Kr + s * Hr
            Al = np.linalg.inv(Cr / DT + G / 2); Ar = Cr / DT - G / 2; t = t0
            while t < t1 - 1e-9:
                z = Al @ (Ar @ z + fr * 15.0 * q * (float(on(t)) + float(on(t + DT))) / 2); t += DT
            Z[e, :m] = z
        return Z

    out = {"note": "クリロフ ROM（r=10＋一様成分）＋EnKF（温度計 P2・P0）。全節点の温度の RMSE [K]、乱数 5 通りの平均", "cases": {}}
    hist = {}
    for case in ("learned", "q25", "intermittent"):
        z_ = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z_["times"]; Tf = z_["Tfield"]
        Tt = np.array([map_temperature(Mmap, Tf[k]) for k in range(len(t))]) - rg.T_AIR_K
        yobs = Tf[:, pod][:, [2, 0]] - rg.T_AIR_K
        on = SCHED[case]; rms, qs, ss = [], [], []
        for seed in SEEDS:
            rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
            Z = np.zeros((N_ENS, m + 2))
            Z[:, :m] = rng.uniform(-3, 12, N_ENS)[:, None] * one[None]     # 一様な初期温度のずれ
            Z[:, m] = rng.uniform(0.3, 1.8, N_ENS); Z[:, m + 1] = np.exp(rng.uniform(np.log(0.1), np.log(1.5), N_ENS))
            R = np.diag([SIG ** 2] * 2); rm = [np.nan]; q_ = [Z[:, m].mean()]; s_ = [Z[:, m + 1].mean()]; tp = 0.0
            for k, tb in enumerate(t[1:], 1):
                Z = forecast(Z, tp, tb, on); tp = tb
                y = yobs[k] + ro.normal(0, SIG, 2)
                Z = enkf_update(Z, y, None, R, rng, inflation=INFL, Yf=Z[:, :m] @ Hobs.T)
                Z[:, m] = np.clip(Z[:, m], 0, 3); Z[:, m + 1] = np.clip(Z[:, m + 1], 0.01, 3)
                th = V @ Z[:, :m].mean(0); rm.append(float(np.sqrt(((th - Tt[k]) ** 2).mean())))
                q_.append(Z[:, m].mean()); s_.append(Z[:, m + 1].mean())
            rms.append(rm); qs.append(q_); ss.append(s_)
        rms = np.nanmean(np.array(rms, float), 0) if False else np.mean(np.array(rms, float), 0)
        qs = np.mean(qs, 0); ss = np.mean(ss, 0); heat = (t > 0) & (t <= 300)
        hist[case] = (t, rms, qs, ss)
        out["cases"][case] = dict(rmse_heat_mean_K=float(np.nanmean(rms[heat])), rmse_after100_max_K=float(np.nanmax(rms[t >= 100])),
                                  rmse_end_K=float(rms[-1]), q_end=float(qs[-1]), q_true=QTRUE[case], h_scale_end=float(ss[-1]))
        print(case, {k: round(v, 3) for k, v in out["cases"][case].items()}, flush=True)
    json.dump(out, open(os.path.join(RES, "krylov_da.json"), "w"), ensure_ascii=False, indent=1)
    np.savez(os.path.join(RES, "krylov_da.npz"), **{f"{c}_{k}": v for c, h in hist.items() for k, v in zip(("t", "rmse", "q", "s"), h)})

    ho = np.load(os.path.join(RES, "holdout_unobserved_cells.npz"))
    plt.rcParams.update({"font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
    fig, axs = plt.subplots(2, 3, figsize=(16, 8.6), sharex=True, gridspec_kw=dict(hspace=0.14, wspace=0.27))
    titles = {"learned": "15 W・0〜300 s", "q25": "25 W・0〜300 s", "intermittent": "15 W 間欠加熱"}
    for j, case in enumerate(("learned", "q25", "intermittent")):
        t, rms, qs, ss = hist[case]
        ax = axs[0, j]
        ax.plot(ho[f"{case}_t"], ho[f"{case}_temp_rmse"], "-", color="#999999", lw=2.5, label="今の ROM＋同化（温度 2 点）")
        ax.plot(ho[f"{case}_t"], ho[f"{case}_sel_rmse"], "-", color="#C0392B", lw=2.5, label="今の ROM＋同化（温度 2 点＋変位 2 点）")
        ax.plot(t, rms, "-", color="#1f77b4", lw=2.8, marker="o", ms=5, label="クリロフ ROM＋同化（温度 2 点）")
        ax.set_yscale("log"); ax.set_ylim(0.02, 5); ax.grid(alpha=.3, which="both"); ax.set_title(titles[case])
        ax = axs[1, j]
        ax.plot(t, qs, "-", color="#d62728", lw=2.4, label="発熱の倍率 q（推定）"); ax.axhline(QTRUE[case], color="#d62728", ls=":", lw=1.6)
        ax.plot(t, ss, "-", color="#2ca02c", lw=2.4, label="熱伝達率の倍率 s（推定）")
        ax.axhline(0.21, color="#2ca02c", ls=":", lw=1.6)
        ax.set_ylim(0, 2.4); ax.grid(alpha=.3); ax.set_xlabel("時刻 [s]"); ax.set_xticks([0, 200, 400, 600])
    axs[0, 0].set_ylabel("全体の温度の RMSE [K]"); axs[1, 0].set_ylabel("推定した倍率")
    
    h1, l1 = axs[0, 0].get_legend_handles_labels(); h2, l2 = axs[1, 0].get_legend_handles_labels()
    from matplotlib.lines import Line2D
    h2 += [Line2D([], [], color="#d62728", ls=":", lw=1.6), Line2D([], [], color="#2ca02c", ls=":", lw=1.6)]
    l2 += ["q の正解", "s=0.21（15 W に合わせた値）"]
    fig.legend(h1 + h2, l1 + l2, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 0.02), frameon=False)
    fig.savefig(os.path.join(IMG, "krylov_da.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote docs/img/krylov_da.png")


if __name__ == "__main__":
    main()
