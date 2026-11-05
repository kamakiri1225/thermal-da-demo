"""クリロフ部分空間法で作る低次元モデル（ROM）を試し、今の ROM（POD・Q-DEIM＋熱回路）と OpenFOAM で比べる。

  1. 固体の熱伝導の有限要素モデル（FEM）を作る（FrontISTR と同じ六面体メッシュ 4×48×20、5,040 節点）
       C dθ/dt + (K + H) θ = f P(t)      θ = T − T_air
     C：熱容量、K：熱伝導、H：表面の熱伝達（面ごとの h は blog_001 の CHT の実効値）、f：ヒータ 1 W あたりの熱流
  2. クリロフ部分空間法（s = 0 のまわりのモーメント合わせ）で r 次元に縮める
       V = span{ G⁻¹f, (G⁻¹C) G⁻¹f, (G⁻¹C)² G⁻¹f, … }（G = K + H、C の内積で正規直交化）
       縮めた式：Vᵀ C V dz/dt + Vᵀ G V z = Vᵀ f P(t)、θ ≈ V z（全節点の温度がそのまま出る）
  3. OpenFOAM（CHT）の温度を節点へ移した正解と、3 条件（15 W・25 W・間欠加熱）で比べる
     比べるもの：全 FEM（縮める前）、クリロフ ROM（r = 5・10・20）、今の ROM（5 点の熱回路＋gappy-POD）
  4. 参考：表面の h を一律に何倍かする 1 つの係数だけを、15 W の OpenFOAM に合わせた版も作る

出力：results/krylov_rom.json、docs/img/krylov_rom_compare.png
再現：OMP_NUM_THREADS=4 python3 run/krylov_rom.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy.optimize import minimize_scalar
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
import cylinder_mesh
from dacore import plots as _p  # noqa: F401
from dacore import rom_general as rg
from dacore.temperature_mapping import build_idw_mapping, map_temperature
from improve_heating_schedule import integrate, SCHED, QTRUE
import matplotlib.pyplot as plt

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NR, NTH, NZ = 4, 48, 20; RI, RO, HH = 0.020, 0.0375, 0.1005
RHO, CP, KAP = 7850.0, 480.0, 50.0                     # 102_0 の固体の物性（鋼）
HTC = dict(outer_px=15.8, outer_mx=10.5, inner=10.8, top=11.9, bottom=19.5)   # blog_001 の実効 h [W/(m²K)]
HEAT_Z = (0.02525, 0.07525); HEAT_HALF = 0.100 / RO / 2  # ヒータ：高さ 25〜75 mm、周長 100 mm（+x 中心）
DT = 2.0
CASES = [("learned", "15 W・0〜300 s\n（今の ROM を作った条件）"), ("q25", "25 W・0〜300 s"), ("intermittent", "15 W 間欠加熱")]
G8 = np.array([[a, b, c] for c in (-1, 1) for b in (-1, 1) for a in (-1, 1)]) / np.sqrt(3)
XI = np.array([[-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1], [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]], float)


def hex_matrices(X):
    """六面体 8 節点の熱伝導・熱容量の要素行列（ne,8,8）"""
    ne = len(X); Ke = np.zeros((ne, 8, 8)); Ce = np.zeros((ne, 8, 8))
    for g in G8:
        N = np.prod(1 + XI * g, axis=1) / 8
        dN = np.stack([XI[:, d] * np.prod(np.delete(1 + XI * g, d, axis=1), axis=1) / 8 for d in range(3)], axis=1)  # (8,3)
        J = np.einsum("ai,eaj->eij", dN, X); dJ = np.abs(np.linalg.det(J))
        dNdx = np.einsum("eij,aj->eai", np.linalg.inv(J), dN)
        Ke += KAP * np.einsum("eai,ebi->eab", dNdx, dNdx) * dJ[:, None, None]
        Ce += RHO * CP * np.einsum("a,b->ab", N, N)[None] * dJ[:, None, None]
    return Ke, Ce


def quad_matrices(Y):
    """四角形 4 節点の面の ∫N Nᵀ dA（nq,4,4）と ∫N dA（nq,4）"""
    s = np.array([-1, 1, 1, -1], float); t = np.array([-1, -1, 1, 1], float)
    M = np.zeros((len(Y), 4, 4)); L = np.zeros((len(Y), 4))
    for gs in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
        for gt in (-1 / np.sqrt(3), 1 / np.sqrt(3)):
            N = (1 + s * gs) * (1 + t * gt) / 4; Ns = s * (1 + t * gt) / 4; Nt = t * (1 + s * gs) / 4
            dA = np.linalg.norm(np.cross(np.einsum("a,qai->qi", Ns, Y), np.einsum("a,qai->qi", Nt, Y)), axis=1)
            M += np.einsum("a,b->ab", N, N)[None] * dA[:, None, None]; L += N[None] * dA[:, None]
    return M, L


def build_fem(h_scale=1.0):
    m = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, RI, RO, HH)
    X = np.array([x for _, x in m["nodes"]]); pos = {n: i for i, (n, _) in enumerate(m["nodes"])}
    E = np.array([[pos[n] for n in c] for _, c in m["elements"]]); n = len(X)
    Ke, Ce = hex_matrices(X[E])
    rows = np.repeat(E, 8, axis=1).ravel(); cols = np.tile(E, (1, 8)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(n, n)); C = sp.csr_matrix((Ce.ravel(), (rows, cols)), shape=(n, n))
    ir = np.arange(len(E)) % NR; iz = np.arange(len(E)) // (NR * NTH)
    faces = {"outer": (ir == NR - 1, [1, 2, 6, 5]), "inner": (ir == 0, [0, 3, 7, 4]),
             "bottom": (iz == 0, [0, 1, 2, 3]), "top": (iz == NZ - 1, [4, 5, 6, 7])}
    Hr, Hc, Hv = [], [], []; f = np.zeros(n); A_heat = 0.0; qs = []
    for name, (mask, loc) in faces.items():
        Q = E[mask][:, loc]; Mq, Lq = quad_matrices(X[Q]); cen = X[Q].mean(1)
        if name == "outer":
            h = np.where(cen[:, 0] > 0, HTC["outer_px"], HTC["outer_mx"])
            th = np.arctan2(cen[:, 1], cen[:, 0])
            heat = (np.abs(th) <= HEAT_HALF) & (cen[:, 2] >= HEAT_Z[0]) & (cen[:, 2] <= HEAT_Z[1])
            qs.append((Q[heat], Lq[heat])); A_heat += Lq[heat].sum()
        else:
            h = np.full(len(Q), HTC[name])
        h = h * h_scale
        Hr.append(np.repeat(Q, 4, axis=1).ravel()); Hc.append(np.tile(Q, (1, 4)).ravel()); Hv.append((Mq * h[:, None, None]).ravel())
    H = sp.csr_matrix((np.concatenate(Hv), (np.concatenate(Hr), np.concatenate(Hc))), shape=(n, n))
    for Q, Lq in qs: np.add.at(f, Q.ravel(), Lq.ravel() / A_heat)   # 1 W を面積に比例して配る
    return dict(X=X, C=C.tocsc(), G=(K + H).tocsc(), f=f)


def run_full(fem, Pfun, t_out):
    """全 FEM をクランク・ニコルソン法で 0→600 s。t_out の時刻の θ を返す"""
    C, G, f = fem["C"], fem["G"], fem["f"]
    Al = splu((C / DT + G / 2).tocsc()); Ar = (C / DT - G / 2).tocsr()
    th = np.zeros(len(f)); out = [th.copy()]; t = 0.0
    while t < t_out[-1] - 1e-9:
        th = Al.solve(Ar @ th + f * (Pfun(t) + Pfun(t + DT)) / 2); t += DT
        if np.any(np.isclose(t, t_out[1:])): out.append(th.copy())
    return np.array(out)


def krylov_basis(fem, r):
    """G⁻¹C のクリロフ部分空間（G⁻¹f から始める）を C の内積で正規直交化"""
    C, G, f = fem["C"], fem["G"], fem["f"]; lu = splu(G)
    V = np.zeros((len(f), r)); v = lu.solve(f)
    for j in range(r):
        for _ in range(2):
            v = v - V[:, :j] @ (V[:, :j].T @ (C @ v))
        V[:, j] = v / np.sqrt(v @ (C @ v)); v = lu.solve(C @ V[:, j])
    return V


def run_rom(fem, V, Pfun, t_out):
    Cr = V.T @ (fem["C"] @ V); Gr = V.T @ (fem["G"] @ V); fr = V.T @ fem["f"]
    Al = np.linalg.inv(Cr / DT + Gr / 2); Ar = Cr / DT - Gr / 2
    z = np.zeros(V.shape[1]); out = [z.copy()]; t = 0.0
    while t < t_out[-1] - 1e-9:
        z = Al @ (Ar @ z + fr * (Pfun(t) + Pfun(t + DT)) / 2); t += DT
        if np.any(np.isclose(t, t_out[1:])): out.append(z.copy())
    return np.array(out) @ V.T


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz")); Cn = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], 5); hh = float(d["h"]); heat = int(d["heat_node"])
    fem = build_fem(); X = fem["X"]; M = build_idw_mapping(Cc, X)
    nodeP = [int(np.linalg.norm(X - Cc[c], axis=1).argmin()) for c in pod]
    truth = {}
    for case, _ in CASES:
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]
        truth[case] = (t, np.array([map_temperature(M, z["Tfield"][k]) for k in range(len(t))]) - rg.T_AIR_K)
    # 参考：h を一律に s 倍する係数を、15 W の OpenFOAM に合わせる
    t15, T15 = truth["learned"]
    def loss(s):
        return float(np.sqrt(((run_full(build_fem(s), QTRUE["learned"], t15) - T15) ** 2).mean()))
    opt = minimize_scalar(loss, bounds=(0.01, 3.0), method="bounded", options=dict(xatol=1e-3))
    s_fit = float(opt.x); fem_fit = build_fem(s_fit); print("h の倍率（15 W に合わせた）", round(s_fit, 3), "RMSE", round(opt.fun, 4))
    bases = {r: krylov_basis(fem, r) for r in (5, 10, 20)}; Vfit = krylov_basis(fem_fit, 10)
    out = {"note": "クリロフ ROM と今の ROM を、OpenFOAM（節点へ移した温度）と比較。単位 K。時刻 0〜600 s（30 s ごと）",
           "h_table_W_m2K": HTC, "h_scale_fit": s_fit, "n_nodes": int(len(X)), "cases": {}}
    res = {}
    for case, _ in CASES:
        t, Tt = truth[case]; P = QTRUE[case]; r = {}
        r["全 FEM（h は CHT の実効値）"] = run_full(fem, P, t)
        for k, V in bases.items():
            t0 = time.perf_counter(); r[f"クリロフ ROM r={k}"] = run_rom(fem, V, P, t); dt_ = time.perf_counter() - t0
            if k == 10: out.setdefault("time_rom_r10_s", dt_)
        r["クリロフ ROM r=10（h を合わせた版）"] = run_rom(fem_fit, Vfit, P, t)
        q = {"learned": 1.0, "q25": 25 / 15, "intermittent": 1.0}[case]
        Tp = np.full((1, 5), rg.T_AIR_K); tr = [Tp[0].copy()]
        for k in range(1, len(t)):
            Tp = integrate(Tp, Cn, Km, np.array([hh]), np.array([q]), heat, SCHED[case], t[k - 1], t[k]); tr.append(Tp[0].copy())
        field = mean[None] + ((np.array(tr) - mean[pod]) @ np.linalg.inv(U[pod]).T) @ U.T
        r["今の ROM（5 点の熱回路＋gappy-POD）"] = np.array([map_temperature(M, fk) for fk in field]) - rg.T_AIR_K
        res[case] = (t, Tt, r)
        out["cases"][case] = {nm: dict(rmse_max_K=float(np.sqrt(((v - Tt) ** 2).mean(1)).max()),
                                       abs_max_K=float(np.abs(v - Tt).max()),
                                       P_abs_max_K=float(np.abs(v[:, nodeP] - Tt[:, nodeP]).max())) for nm, v in r.items()}
        for nm, v in out["cases"][case].items(): print(case, nm, {k: round(x, 3) for k, x in v.items()})
    json.dump(out, open(os.path.join(RES, "krylov_rom.json"), "w"), ensure_ascii=False, indent=1)

    plt.rcParams.update({"font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
    # 線の描き方：全 FEM は太い灰、r=5 は丸の印、r=10 は四角の印（3 本はほぼ重なるので、印で見分ける）
    show = [("全 FEM（h は CHT の実効値）", dict(color="#8c564b", lw=5, alpha=.45)),
            ("クリロフ ROM r=5", dict(color="#1f77b4", lw=1.6, ls="--", marker="o", ms=7, mfc="white", markevery=(0, 3))),
            ("クリロフ ROM r=10", dict(color="#2ca02c", lw=1.6, ls=":", marker="s", ms=6, markevery=(1, 3))),
            ("クリロフ ROM r=10（h を合わせた版）", dict(color="#9467bd", lw=2.2)),
            ("今の ROM（5 点の熱回路＋gappy-POD）", dict(color="#C0392B", lw=2.2))]
    fig, axs = plt.subplots(3, 3, figsize=(16, 12.5), sharex=True, gridspec_kw=dict(hspace=0.16, wspace=0.27))
    for j, (case, lab) in enumerate(CASES):
        t, Tt, r = res[case]
        ax = axs[0, j]; ax.plot(t, Tt[:, nodeP[2]] + 20, "-", color="k", lw=9, alpha=.18, label="OpenFOAM（正解）")
        for nm, kw in show: ax.plot(t, r[nm][:, nodeP[2]] + 20, label=nm, **kw)
        ax.set_title(lab); ax.grid(alpha=.3)
        ax = axs[1, j]
        for nm, kw in show: ax.plot(t, np.sqrt(((r[nm] - Tt) ** 2).mean(1)), **kw)
        ax.set_yscale("log"); ax.set_ylim(3e-3, 5); ax.grid(alpha=.3, which="both")
        ax.text(0.03, 0.95, "全 FEM・r=5・r=10 は重なる", transform=ax.transAxes, va="top", fontsize=14, color="#555")
        ax = axs[2, j]                     # 縮める誤差：全 FEM との差（全節点の最大）
        full = r["全 FEM（h は CHT の実効値）"]
        for nm, kw in show[1:3]: ax.plot(t, np.abs(r[nm] - full).max(1) + 1e-6, **kw)
        V20 = r["クリロフ ROM r=20"]; ax.plot(t, np.abs(V20 - full).max(1) + 1e-6, color="#ff7f0e", lw=2, label="クリロフ ROM r=20")
        ax.set_yscale("log"); ax.set_ylim(1e-5, 1); ax.grid(alpha=.3, which="both"); ax.set_xlabel("時刻 [s]"); ax.set_xticks([0, 200, 400, 600])
    axs[0, 0].set_ylabel("P2 の温度 [°C]"); axs[1, 0].set_ylabel("OpenFOAM との差\n全節点の RMSE [K]"); axs[2, 0].set_ylabel("縮める誤差\n|ROM − 全 FEM| の最大 [K]")
    h_, l_ = axs[0, 0].get_legend_handles_labels(); h2, l2 = axs[2, 0].get_legend_handles_labels()
    fig.legend(h_ + h2[-1:], l_ + l2[-1:], loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.06), frameon=False)
    fig.savefig(os.path.join(IMG, "krylov_rom_compare.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote docs/img/krylov_rom_compare.png")


if __name__ == "__main__":
    main()
