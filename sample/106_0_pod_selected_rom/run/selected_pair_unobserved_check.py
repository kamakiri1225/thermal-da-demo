"""選んだセンサ配置で、観測していない温度・変位がどう合うか（blog_008 §5-2・§6）.

真値は OpenFOAM（温度場）＋ FrontISTR（全節点の変位）。条件は 15 W・25 W・15 W 間欠加熱。
同化は ROM＋EnKF（60メンバー、30 秒ごと、5 seed）。
評価（すべて観測に使っていない）
  温度：代表点 P1・P3・P4、代表点以外のセル4点、全20,696セルのRMSE
  変位：反り A−O、Uz(A)、Uz(O)、ほかの未観測点の Uz（D・E・F）

構成のセット
  set=s  : 温度 P2・P0 に、変位2点（なし／底面近く／中段 B・C／選定 S1・S2）を足す
  set=n  : 代表点以外から σ で選んだ温度2点＋変位2点（run/select_nonrep_temp_disp.py）と、今の配置の比較

出力: results/unobserved_check_<set>_<case>.json, docs/img/unobserved_<set>_<case>.png
再現: OMP_NUM_THREADS=4 python3 run/selected_pair_unobserved_check.py [s|n] [case ...]
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from improve_heating_schedule import integrate as integ_on

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NPT = 5; IQ = 5; IH = 6; NAUG = 7; N_ENS = 60; SIG_T = 0.30; SIG_U = 0.30; INFL = 1.02
SEEDS = [20260913, 20260914, 20260915, 20260916, 20260917]
H_ = 0.1005
A_XYZ = (0.028, 0.0, H_); O_XYZ = (-0.028, 0.0, H_)
P2P0 = [("rep", 2), ("rep", 0)]
SEL = [((0.0097, -0.0362, H_), 0), ((-0.0346, 0.0144, H_), 2)]          # §5 で選んだ S1(Ux)・S2(Uz)
LOW = [((0.0375, 0.0, 0.005), 2), ((-0.0375, 0.0, 0.005), 2)]           # 底面近く（悪い例）
BC = [((0.028, 0.0, 0.075), 2), ((-0.028, 0.0, 0.075), 2)]              # 中段 B・C の Uz（人が決めた例）
FREE_T = [("ヒータ側・上部", (0.0351, -0.0008, 0.0946)), ("反対側・上部", (-0.0375, 0.0, 0.090)),
          ("横（+y）・中段", (0.0, 0.0355, 0.050)), ("反対側・下部", (-0.0329, -0.0016, 0.0297))]
FREE_U = [("D（上面・+y 側）", (0.0, 0.0362, H_)), ("E（ヒータ側・中段 z=50）", (0.0375, 0.0, 0.050)),
          ("F（反対側・中段 z=50）", (-0.0375, 0.0, 0.050))]
SCHED = {"learned": lambda t: (t < 300.0), "q25": lambda t: (t < 300.0),
         "intermittent": lambda t: (t < 150.0) | ((t >= 300.0) & (t < 450.0))}
CASES = {"learned": "15 W", "q25": "25 W", "intermittent": "15 W 間欠加熱"}
COMPI = {"Ux": 0, "Uy": 1, "Uz": 2}


def configs(setname):
    if setname == "t2":
        js = json.load(open(os.path.join(RES, "two_nonrep_selection.json"), encoding="utf-8"))
        tsel = [("xyz", tuple(np.array(p["xyz_mm"]) / 1000)) for p in js["modes"]["T2"]]
        usel = [(tuple(np.array(p["xyz_mm"]) / 1000), COMPI[p["comp"]]) for p in js["modes"]["U2"]]
        return ("センサ2本だけで反り A−O を当てる：代表点以外から σ で選んだ2本と、代表点の温度計2本",
                [("データ同化なし（ROM だけで予測）", None, [], "#E1A100", "-"),
                 ("温度 P2・P0（代表点・片側）", P2P0, [], "#9AA5B1", "--"),
                 ("温度 P2・P4（代表点・両側）", [("rep", 2), ("rep", 4)], [], "#6B8E23", "--"),
                 ("変位 S1・S2 のみ（§5 で選んだ2点）", [], SEL, "#8E44AD", "--"),
                 ("温度2本（代表点以外、σ で選んだ）", tsel, [], "#2E86C1", "-"),
                 ("変位2本（代表点以外、σ で選んだ）", [], usel, "#C0392B", "-")])
    if setname == "v":
        return ("検証3：A・O を一度も観測せずに、3条件で確かめる",
                [("データ同化なし（ROM だけで予測）", None, [], "#E1A100", "-"),
                 ("温度 P2・P0 だけでデータ同化", P2P0, [], "#9AA5B1", "-"),
                 ("A・O 自身を測る（参考）", P2P0, [(A_XYZ, 2), (O_XYZ, 2)], "#2E86C1", "-"),
                 ("選んだ別の2点 S1・S2（本研究）", P2P0, SEL, "#C0392B", "-")])
    if setname == "s":
        return ("温度計 P2・P0 に変位計2点を足す：置き場所による違い",
                [("データ同化なし（ROM だけで予測）", None, [], "#E1A100", "-"),
                 ("温度 P2・P0 のみ", P2P0, [], "#9AA5B1", "--"),
                 ("＋底面近くの変位2点（悪い例）", P2P0, LOW, "#8E44AD", "--"),
                 ("＋中段 B・C の変位2点", P2P0, BC, "#E67E22", "--"),
                 ("＋選定 S1・S2（σ で選んだ）", P2P0, SEL, "#C0392B", "-")])
    js = json.load(open(os.path.join(RES, "nonrep_selection.json"), encoding="utf-8"))
    tsel = [("xyz", tuple(np.array(p["xyz_mm"]) / 1000)) for p in js["picked"] if p["kind"] == "温度"]
    usel = [(tuple(np.array(p["xyz_mm"]) / 1000), COMPI[p["comp"]]) for p in js["picked"] if p["kind"] == "変位"]
    return ("温度計も代表点以外から選ぶ：σ で選んだ温度2点＋変位2点と、今の配置",
            [("データ同化なし（ROM だけで予測）", None, [], "#E1A100", "-"),
             ("温度 P2・P0 のみ", P2P0, [], "#9AA5B1", "--"),
             ("温度 P2・P0 ＋ 変位 S1・S2（今の配置）", P2P0, SEL, "#E67E22", "--"),
             ("代表点以外の温度2点のみ（σ で選んだ）", tsel, [], "#2E86C1", "--"),
             ("代表点以外の温度2点＋変位2点（σ で選んだ）", tsel, usel, "#C0392B", "-")])


def main(setname, case):
    on = SCHED[case]
    title, CFG = configs(setname)
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); mp = mean[pod]; G = U @ UPp            # 全セル温度 = mean + G (T5 - mp)
    op = np.load(os.path.join(RES, "dispop_allnodes.npz"))
    um = op["u_mean"]; Wall = np.einsum("nck,kj->ncj", op["D"], UPp); co = op["coords"]
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    nc = lambda x: int(np.linalg.norm(Cc - np.array(x), axis=1).argmin())
    cell_of = lambda o: int(pod[o[1]]) if o[0] == "rep" else nc(o[1])
    iA, iO = nn(A_XYZ), nn(O_XYZ)
    z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]; Tf = z["Tfield"]
    ut = np.load(os.path.join(RES, f"truth_disp_all_{case}.npz"))["u"]   # (時刻, 節点, 3) µm

    tcells = [("P1（代表点）", int(pod[1])), ("P3（代表点）", int(pod[3])), ("P4（代表点）", int(pod[4]))] + \
             [(f"{lab}（代表点以外）", nc(x)) for lab, x in FREE_T]
    unodes = [("反り A−O", None), ("Uz(A)", iA), ("Uz(O)", iO)] + [(lab, nn(x)) for lab, x in FREE_U]
    truth_u = lambda k, i: ut[k, iA, 2] - ut[k, iO, 2] if i is None else ut[k, i, 2]

    def est_u(dm, i):
        f = lambda j: um[j, 2] + Wall[j, 2] @ dm
        return f(iA) - f(iO) if i is None else f(i)

    def run(tobs, dobs, seed):
        noda = tobs is None
        tc = [] if noda else [cell_of(o) for o in tobs]
        drow = [(nn(x), c) for x, c in dobs]
        rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
        Z = np.zeros((N_ENS, NAUG)); Z[:, :NPT] = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N_ENS, NPT))
        Z[:, IQ] = rng.uniform(0.3, 1.8, N_ENS); Z[:, IH] = np.clip(rng.normal(0.02, 0.01, N_ENS), 1e-3, 0.1)
        Rd = np.diag([SIG_T**2] * len(tc) + [SIG_U**2] * len(drow))
        Tc, Uc, rms = [], [], []; tp = 0.0
        for ci, tb in enumerate(t[1:], 1):
            Z = Z.copy(); Z[:, :NPT] = integ_on(Z[:, :NPT], C, Km, Z[:, IH], Z[:, IQ], heat, on, tp, tb); tp = tb
            dZ = Z[:, :NPT] - mp
            if noda:
                dm = Z[:, :NPT].mean(0) - mp
                Tc.append([mean[j] + G[j] @ dm for _, j in tcells]); Uc.append([est_u(dm, i) for _, i in unodes])
                rms.append(np.sqrt(np.mean((mean + G @ dm - Tf[ci]) ** 2))); continue
            yv = [Tf[ci, j] for j in tc] + [ut[ci, i, c] for i, c in drow]
            Yf = np.column_stack([mean[j] + dZ @ G[j] for j in tc] + [um[i, c] + dZ @ Wall[i, c] for i, c in drow])
            y = np.array(yv) + ro.normal(0, np.sqrt(np.diag(Rd)))
            Z = enkf_update(Z, y, None, Rd, rng, inflation=INFL, Yf=Yf)
            Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
            dm = Z[:, :NPT].mean(0) - mp
            Tc.append([mean[j] + G[j] @ dm for _, j in tcells]); Uc.append([est_u(dm, i) for _, i in unodes])
            rms.append(np.sqrt(np.mean((mean + G @ dm - Tf[ci]) ** 2)))
        return np.array(Tc), np.array(Uc), np.array(rms)

    tt = t[1:]
    Ttrue = np.array([[Tf[k, j] for _, j in tcells] for k in range(1, len(t))])
    Utrue = np.array([[truth_u(k, i) for _, i in unodes] for k in range(1, len(t))])
    ev = (tt >= 30) & (tt <= 300)
    out = {"note": f"真値＝OpenFOAM＋FrontISTR（{CASES[case]}）。誤差は30〜300秒の平均絶対誤差（温度場は全セルRMSEの時間平均）、5 seed平均。",
           "configs": {}}
    traj = {}
    for nm, tobs, dobs, _, _ in CFG:
        r = [run(tobs, dobs, s) for s in SEEDS]
        traj[nm] = tuple(np.mean([x[i] for x in r], 0) for i in range(3))
        eT = np.mean([np.abs(x[0] - Ttrue)[ev].mean(0) for x in r], 0)
        eU = np.mean([np.abs(x[1] - Utrue)[ev].mean(0) for x in r], 0)
        out["configs"][nm] = {"field_rmse_K": round(float(np.mean([x[2][ev].mean() for x in r])), 3),
                              "temperature_K": {lab: round(float(v), 3) for (lab, _), v in zip(tcells, eT)},
                              "displacement_um": {lab: round(float(v), 3) for (lab, _), v in zip(unodes, eU)}}
        c_ = out["configs"][nm]
        print(f"[{setname} {case}] {nm}: 温度場 {c_['field_rmse_K']} K  反り {c_['displacement_um']['反り A−O']} µm  "
              f"Uz(A) {c_['displacement_um']['Uz(A)']}  Uz(O) {c_['displacement_um']['Uz(O)']}", flush=True)
    json.dump(out, open(os.path.join(RES, f"unobserved_check_{setname}_{case}.json"), "w"), ensure_ascii=False, indent=1)
    np.savez_compressed(os.path.join(RES, f"unobserved_traj_{setname}_{case}.npz"), t=tt, Ttrue=Ttrue, Utrue=Utrue,
                        tlabels=np.array([l for l, _ in tcells]), ulabels=np.array([l for l, _ in unodes]),
                        names=np.array([c[0] for c in CFG]), colors=np.array([c[3] for c in CFG]),
                        T=np.array([traj[c[0]][0] for c in CFG]), U=np.array([traj[c[0]][1] for c in CFG]),
                        R=np.array([traj[c[0]][2] for c in CFG]))
    if setname == "v":
        return

    # ---- 図 ----------------------------------------------------------------
    fig = plt.figure(figsize=(21, 26))
    gs = fig.add_gridspec(5, 4, height_ratios=[1.25, 1, 1, 1, 1], hspace=0.45, wspace=0.27, top=(0.885 if len(CFG) > 5 else 0.90), bottom=0.03)
    th = np.linspace(0, 2 * np.pi, 200)
    axm = fig.add_subplot(gs[0, 0:2]); axs = fig.add_subplot(gs[0, 2:4])
    for r_ in (20, 37.5): axm.plot(r_ * np.cos(th), r_ * np.sin(th), color="#999", lw=1)
    axs.add_patch(plt.Rectangle((-37.5, 0), 75, 100.5, fc="#f2f2f2", ec="#999"))
    def mark(x, mk, col, lab, s=170):
        X = np.array(x) * 1000
        axm.scatter(X[0], X[1], marker=mk, s=s, c=col, edgecolor="k", zorder=5, label=lab)
        axs.scatter(X[0], X[2], marker=mk, s=s, c=col, edgecolor="k", zorder=5, label=lab)
    seen = set()
    for nm, tobs, dobs, col, _ in CFG:
        for o in (tobs or []):
            key = ("T", str(o))
            if key in seen: continue
            seen.add(key); mark(Cc[cell_of(o)], "o", col, f"温度計：{nm}")
        for x, c in dobs:
            key = ("U", x, c)
            if key in seen: continue
            seen.add(key); mark(x, "s", col, f"変位計：{nm}")
    for k, (_, j) in enumerate(tcells): mark(Cc[j], "^", "#7FB3D5", "温度を評価（観測しない）" if k == 0 else None, s=130)
    for k, x in enumerate([A_XYZ, O_XYZ] + [x for _, x in FREE_U]):
        mark(x, "*", "#111", "変位を評価（観測しない）" if k == 0 else None, s=260)
    for ax in (axm, axs): ax.set_aspect("equal"); ax.grid(alpha=.3)
    axm.set_title("上から見た位置（x–y）", fontsize=16); axm.set_xlabel("x [mm]"); axm.set_ylabel("y [mm]")
    axm.text(44, 0, "ヒータ側\n（+x）", fontsize=13, color="#a0522d", va="center")
    axs.set_title("横から見た位置（x–z）", fontsize=16); axs.set_xlabel("x [mm]"); axs.set_ylabel("z [mm]")
    h_, l_ = axs.get_legend_handles_labels(); uniq = dict(zip(l_, h_))
    axs.legend(uniq.values(), uniq.keys(), fontsize=11.5, loc="center left", bbox_to_anchor=(1.02, 0.5))

    def shade(ax):
        for a0, a1 in ([(0, 150), (300, 450)] if case == "intermittent" else [(0, 300)]):
            ax.axvspan(a0, a1, color="#FDEBD0", alpha=.45)

    def errbox(ax, vals, unit):
        txt = "\n".join(f"{v:.2f} {unit}" for v in vals)
        ys = 0.97
        for (nm, _, _, col, _), v in zip(CFG, vals):
            ax.text(0.98, ys, f"{v:.2f} {unit}", transform=ax.transAxes, ha="right", va="top", fontsize=11.5,
                    color=col, weight="bold", bbox=dict(fc="white", ec="none", alpha=.75, pad=1)); ys -= 0.085

    def series(ax, k, kind, title):
        tv = Ttrue[:, k] - 273.15 if kind == "T" else Utrue[:, k]
        shade(ax); ax.plot(tt, tv, "-", color="#111", lw=6, alpha=.30)
        for nm, _, _, col, ls in CFG:
            v = traj[nm][0][:, k] - 273.15 if kind == "T" else traj[nm][1][:, k]
            ax.plot(tt, v, ls, color=col, lw=2.3)
        key = "temperature_K" if kind == "T" else "displacement_um"
        errbox(ax, [out["configs"][nm][key][title] for nm, *_ in CFG], "K" if kind == "T" else "µm")
        ax.set_title(title, fontsize=15); ax.grid(alpha=.3); ax.set_xlabel("時刻 [s]")
        ax.set_ylabel("温度 [°C]" if kind == "T" else "変位 [µm]")

    tslots = [(1, 0), (1, 1), (1, 2), (2, 0), (2, 1), (2, 2), (2, 3)]
    for k, ((lab, _), (r_, c_)) in enumerate(zip(tcells, tslots)): series(fig.add_subplot(gs[r_, c_]), k, "T", lab)
    ax = fig.add_subplot(gs[1, 3]); shade(ax)
    for nm, _, _, col, ls in CFG: ax.plot(tt, traj[nm][2], ls, color=col, lw=2.3)
    errbox(ax, [out["configs"][nm]["field_rmse_K"] for nm, *_ in CFG], "K")
    ax.set_title("温度場全体（20,696セル）の RMSE", fontsize=15); ax.set_ylabel("RMSE [K]"); ax.set_xlabel("時刻 [s]")
    ax.set_yscale("log"); ax.grid(alpha=.3)
    uslots = [(3, 0), (3, 1), (3, 2), (4, 0), (4, 1), (4, 2)]
    for k, ((lab, _), (r_, c_)) in enumerate(zip(unodes, uslots)): series(fig.add_subplot(gs[r_, c_]), k, "U", lab)
    axt = fig.add_subplot(gs[3:5, 3]); axt.axis("off")
    yy = 1.0
    axt.text(0, yy, "観測していない量の誤差\n（30〜300 秒の平均、5 seed）", va="top", fontsize=15, weight="bold"); yy -= 0.10
    for nm, _, _, col, _ in CFG:
        c_ = out["configs"][nm]
        axt.text(0, yy, f"■ {nm}", va="top", fontsize=13, color=col, weight="bold"); yy -= 0.045
        axt.text(0.04, yy, f"温度場 {c_['field_rmse_K']:.2f} K、反り A−O {c_['displacement_um']['反り A−O']:.2f} µm\n"
                 f"Uz(A) {c_['displacement_um']['Uz(A)']:.2f} µm、Uz(O) {c_['displacement_um']['Uz(O)']:.2f} µm",
                 va="top", fontsize=12.5, linespacing=1.5); yy -= (0.135 if len(CFG) > 5 else 0.17)
    handles = [Line2D([], [], color="#111", lw=6, alpha=.3, label="真値（OpenFOAM＋FrontISTR）")] + \
              [Line2D([], [], color=col, lw=2.5, ls=ls, label=nm) for nm, _, _, col, ls in CFG]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.945), ncol=3, fontsize=13.5, frameon=False)
    on_txt = "0〜150 秒と 300〜450 秒" if case == "intermittent" else "0〜300 秒"
    fig.text(0.5, 0.985, f"{title}（真値＝OpenFOAM＋FrontISTR、{CASES[case]}）", ha="center", fontsize=20, weight="bold")
    fig.text(0.5, 0.965, f"グラフの点はすべて観測に使っていない。数字＝30〜300 秒の平均誤差（線の色と対応）。橙の帯＝ヒータ ON（{on_txt}）。1回目の補正（30 秒）以降を表示",
             ha="center", fontsize=13.5, color="#444")
    fig.savefig(os.path.join(IMG, f"unobserved_{setname}_{case}.png"), dpi=90, facecolor="white", bbox_inches="tight")
    plt.close(fig); print("wrote", setname, case, flush=True)


if __name__ == "__main__":
    sets = [a for a in sys.argv[1:] if a in ("s", "n", "v", "t2")] or ["s", "n", "v"]
    cases = [a for a in sys.argv[1:] if a in CASES] or list(CASES)
    for s_ in sets:
        for c in cases:
            main(s_, c)
