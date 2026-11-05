"""クリロフ ROM の温度から熱変形を出す行列を作り、変位の観測も使ってデータ同化する。

  1. 熱弾性の行列を作る（FrontISTR と同じ六面体 8 節点、鋼 E=205 GPa・ν=0.3・α=1.2e-5、3-2-1 の拘束）
       Ku u = L θ        Ku：剛性行列（3n×3n）、L：温度から熱荷重を作る行列（3n×n）
     → u = Ku⁻¹ L θ。クリロフ ROM では θ = V z なので、u = (Ku⁻¹ L V) z = Wz z（3n×m を前もって計算）
  2. 確かめ：OpenFOAM の温度（節点へ移したもの）を入れた変位を、FrontISTR の変位（正解）と比べる
  3. データ同化：温度計 P2・P0 ＋ 変位計 S1（Ux）・S2（Uz）（目的指向で選んだ 2 本）で z・q・s を直す
     比べるもの：全体の温度の RMSE、測らない反り A−O の誤差（今の ROM の結果と同じ指標）

出力：results/krylov_thermoelastic.json、docs/img/krylov_thermoelastic.png
再現：OMP_NUM_THREADS=4 python3 run/krylov_thermoelastic.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import splu
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import krylov_rom as kr
from krylov_da import basis_with_const
from dacore import plots as _p  # noqa: F401
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from dacore.temperature_mapping import build_idw_mapping, map_temperature
from improve_heating_schedule import SCHED
import cylinder_mesh
import matplotlib.pyplot as plt

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
E_, NU, ALPHA = 2.05e11, 0.3, 1.2e-5
R_K = 10; N_ENS = 60; SIG_T = 0.30; SIG_U = 0.30; INFL = 1.02; DT = 2.0
SEEDS = [20260913, 20260914, 20260915, 20260916, 20260917]
QTRUE = {"learned": 1.0, "q25": 25 / 15, "intermittent": 1.0}
S_OBS = [((0.0097, -0.0362, 0.1005), 0), ((-0.0346, 0.0144, 0.1005), 2)]   # 目的指向の S1（Ux）・S2（Uz）
A_XYZ = (0.028, 0.0, 0.1005); O_XYZ = (-0.028, 0.0, 0.1005)


def thermoelastic(X, E):
    """剛性 Ku（3n×3n）と熱荷重の行列 L（3n×n）"""
    lam = E_ * NU / ((1 + NU) * (1 - 2 * NU)); mu = E_ / (2 * (1 + NU))
    D = np.zeros((6, 6)); D[:3, :3] = lam; D[np.arange(3), np.arange(3)] += 2 * mu; D[3:, 3:] = np.eye(3) * mu
    m = np.array([1, 1, 1, 0, 0, 0.0]); Dm = D @ m * ALPHA
    ne = len(E); Xe = X[E]; Ke = np.zeros((ne, 24, 24)); Le = np.zeros((ne, 24, 8))
    for g in kr.G8:
        N = np.prod(1 + kr.XI * g, axis=1) / 8
        dN = np.stack([kr.XI[:, d] * np.prod(np.delete(1 + kr.XI * g, d, axis=1), axis=1) / 8 for d in range(3)], axis=1)
        J = np.einsum("ai,eaj->eij", dN, Xe); dJ = np.abs(np.linalg.det(J))
        dx = np.einsum("eij,aj->eai", np.linalg.inv(J), dN)                  # (ne,8,3)
        B = np.zeros((ne, 6, 24))
        for a in range(8):
            B[:, 0, 3 * a] = dx[:, a, 0]; B[:, 1, 3 * a + 1] = dx[:, a, 1]; B[:, 2, 3 * a + 2] = dx[:, a, 2]
            B[:, 3, 3 * a] = dx[:, a, 1]; B[:, 3, 3 * a + 1] = dx[:, a, 0]
            B[:, 4, 3 * a + 1] = dx[:, a, 2]; B[:, 4, 3 * a + 2] = dx[:, a, 1]
            B[:, 5, 3 * a] = dx[:, a, 2]; B[:, 5, 3 * a + 2] = dx[:, a, 0]
        Ke += np.einsum("eki,kl,elj->eij", B, D, B) * dJ[:, None, None]
        Le += np.einsum("eki,k,b->eib", B, Dm, N) * dJ[:, None, None]
    n = len(X); dof = (3 * E[:, :, None] + np.arange(3)).reshape(ne, 24)
    Ku = sp.csr_matrix((Ke.ravel(), (np.repeat(dof, 24, 1).ravel(), np.tile(dof, (1, 24)).ravel())), shape=(3 * n, 3 * n))
    L = sp.csr_matrix((Le.ravel(), (np.repeat(dof, 8, 1).ravel(), np.tile(E, (1, 24)).ravel())), shape=(3 * n, n))
    return Ku, L


def main():
    mesh = cylinder_mesh.build_cylinder_mesh(kr.NR, kr.NTH, kr.NZ, kr.RI, kr.RO, kr.HH)
    pos = {nid: i for i, (nid, _) in enumerate(mesh["nodes"])}
    X = np.array([x for _, x in mesh["nodes"]]); E = np.array([[pos[n] for n in c] for _, c in mesh["elements"]]); n = len(X)
    Ku, L = thermoelastic(X, E)
    fixed = [3 * pos[i] + c for i in mesh["fix_xyz"] for c in range(3)] + [3 * pos[i] + c for i in mesh["fix_yz"] for c in (1, 2)] + [3 * pos[i] + 2 for i in mesh["fix_z"]]
    free = np.setdiff1d(np.arange(3 * n), fixed); lu = splu(Ku[free][:, free].tocsc())

    def disp(theta):                                                # θ（n,…）→ 変位 µm（3n,…）
        Fv = L @ theta; u = np.zeros((3 * n,) + Fv.shape[1:]); u[free] = lu.solve(Fv[free]); return u * 1e6

    kv = np.load(os.path.join(RES, "qdeim_points.npz")); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    M = build_idw_mapping(Cc, X)
    near = lambda x: int(np.linalg.norm(X - np.array(x), axis=1).argmin()); iA, iO = near(A_XYZ), near(O_XYZ)
    nodeP = [near(Cc[c]) for c in pod]
    srow = [3 * near(x) + c for x, c in S_OBS]; rA, rO = 3 * iA + 2, 3 * iO + 2
    # ---- 2. 確かめ：OpenFOAM の温度 → この行列の変位 と FrontISTR の変位 ----
    check = {}; truth = {}
    for case in ("learned", "q25", "intermittent"):
        z = np.load(os.path.join(RES, f"limit_truth_{case}.npz")); t = z["times"]
        Tt = np.array([map_temperature(M, z["Tfield"][k]) for k in range(len(t))]) - rg.T_AIR_K
        uf = np.load(os.path.join(RES, f"truth_disp_all_{case}.npz"))["u"].reshape(len(t), -1)   # FrontISTR [µm]
        um = disp(Tt.T).T
        check[case] = dict(max_abs_um=float(np.abs(um - uf).max()), umax_um=float(np.abs(uf).max()),
                           warp_max_diff_um=float(np.abs((um[:, rA] - um[:, rO]) - (uf[:, rA] - uf[:, rO])).max()))
        truth[case] = (t, z["Tfield"], Tt, uf)
        print("確かめ", case, {k: round(v, 4) for k, v in check[case].items()}, flush=True)

    # ---- 3. クリロフ ROM ＋ 同化（温度 2 点／温度 2 点＋変位 2 点） ----
    fem1 = kr.build_fem(1.0); fem0 = kr.build_fem(1e-12)
    Hm = fem1["G"] - fem0["G"]; Km = fem0["G"]; C = fem1["C"]; f = fem1["f"]
    V = basis_with_const(fem1, R_K); m = V.shape[1]
    Cr = V.T @ (C @ V); Kr = V.T @ (Km @ V); Hr = V.T @ (Hm @ V); fr = V.T @ f
    Wz = disp(V)                                                    # 変位 = Wz @ z（3n×m）
    one = V.T @ (C @ np.ones(n))
    Ht = V[[nodeP[2], nodeP[0]]]; Hu = Wz[srow]

    def forecast(Z, t0, t1, on):
        Z = Z.copy()
        for e in range(len(Z)):
            zz = Z[e, :m]; q = Z[e, m]; s = Z[e, m + 1]; G = Kr + s * Hr
            Al = np.linalg.inv(Cr / DT + G / 2); Ar = Cr / DT - G / 2; tt = t0
            while tt < t1 - 1e-9:
                zz = Al @ (Ar @ zz + fr * 15.0 * q * (float(on(tt)) + float(on(tt + DT))) / 2); tt += DT
            Z[e, :m] = zz
        return Z

    out = {"note": "クリロフ ROM（r=10＋一様成分）の温度から熱弾性の行列で変位を出し、EnKF で同化。温度は全節点の RMSE [K]、反りは 30〜300 s の絶対誤差の平均 [µm]（5 seed）",
           "check_vs_frontistr": check, "cases": {}}
    hist = {}
    for case in ("learned", "q25", "intermittent"):
        t, Tf, Tt, uf = truth[case]; on = SCHED[case]; heat = (t > 0) & (t <= 300)
        yT = Tf[:, pod][:, [2, 0]] - rg.T_AIR_K; yU = uf[:, srow]; warp_true = uf[:, rA] - uf[:, rO]
        res = {}
        for key, use_u in (("temp", False), ("temp_disp", True)):
            rms, werr = [], []
            for seed in SEEDS:
                rng = np.random.default_rng(seed); ro = np.random.default_rng(seed + 7)
                Z = np.zeros((N_ENS, m + 2)); Z[:, :m] = rng.uniform(-3, 12, N_ENS)[:, None] * one[None]
                Z[:, m] = rng.uniform(0.3, 1.8, N_ENS); Z[:, m + 1] = np.exp(rng.uniform(np.log(0.1), np.log(1.5), N_ENS))
                R = np.diag([SIG_T ** 2] * 2 + ([SIG_U ** 2] * 2 if use_u else [])); tp = 0.0; rm = [np.nan]; we = [np.nan]
                for k, tb in enumerate(t[1:], 1):
                    Z = forecast(Z, tp, tb, on); tp = tb
                    y = list(yT[k] + ro.normal(0, SIG_T, 2)); Yf = Z[:, :m] @ Ht.T
                    if use_u: y += list(yU[k] + ro.normal(0, SIG_U, 2)); Yf = np.column_stack([Yf, Z[:, :m] @ Hu.T])
                    Z = enkf_update(Z, np.array(y), None, R, rng, inflation=INFL, Yf=Yf)
                    Z[:, m] = np.clip(Z[:, m], 0, 3); Z[:, m + 1] = np.clip(Z[:, m + 1], 0.01, 3)
                    zm = Z[:, :m].mean(0); rm.append(float(np.sqrt(((V @ zm - Tt[k]) ** 2).mean())))
                    u = Wz @ zm; we.append(float(abs((u[rA] - u[rO]) - warp_true[k])))
                rms.append(rm); werr.append(we)
            rms = np.mean(rms, 0); werr = np.mean(werr, 0)
            res[key] = dict(rmse_heat_mean_K=float(np.nanmean(rms[heat])), warp_err_heat_mean_um=float(np.nanmean(werr[heat])))
            hist[(case, key)] = (t, rms, werr)
        out["cases"][case] = res; print(case, res, flush=True)
    json.dump(out, open(os.path.join(RES, "krylov_thermoelastic.json"), "w"), ensure_ascii=False, indent=1)

    ho = np.load(os.path.join(RES, "holdout_unobserved_cells.npz"))
    plt.rcParams.update({"font.size": 16, "axes.titlesize": 18, "axes.labelsize": 17, "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15})
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.6), sharey=True, gridspec_kw=dict(wspace=0.1))
    titles = {"learned": "15 W・0〜300 s", "q25": "25 W・0〜300 s", "intermittent": "15 W 間欠加熱"}
    for ax, case in zip(axs, ("learned", "q25", "intermittent")):
        ax.plot(ho[f"{case}_t"], ho[f"{case}_temp_rmse"], "-", color="#999999", lw=2.4, label="今の ROM＋同化（温度 2 点）")
        ax.plot(ho[f"{case}_t"], ho[f"{case}_sel_rmse"], "-", color="#C0392B", lw=2.4, label="今の ROM＋同化（温度 2 点＋変位 2 点）")
        for key, c, mk, lab in (("temp", "#1f77b4", "o", "クリロフ ROM＋同化（温度 2 点）"), ("temp_disp", "#2ca02c", "s", "クリロフ ROM＋同化（温度 2 点＋変位 2 点）")):
            t, rms, _ = hist[(case, key)]; ax.plot(t, rms, "-", color=c, lw=2.6, marker=mk, ms=5, label=lab)
        ax.set_yscale("log"); ax.set_ylim(0.02, 5); ax.grid(alpha=.3, which="both"); ax.set_title(titles[case])
        ax.set_xlabel("時刻 [s]"); ax.set_xticks([0, 200, 400, 600])
    axs[0].set_ylabel("全体の温度の RMSE [K]")
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.0), frameon=False)
    fig.savefig(os.path.join(IMG, "krylov_thermoelastic.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote docs/img/krylov_thermoelastic.png")


if __name__ == "__main__":
    main()
