"""blog_012（Fisher 情報行列でセンサの置き場所を選ぶ）の図と数値。すべて本番のデータから計算する。

  分かりやすくするため、未知数を「POD の型 1・型 2 の重み a1, a2」の 2 つに絞る（型 1・2 で温度の変化の 99.89 %）。
  型は「一番大きく変わるセルで 1」になるように大きさをそろえ、重み a の単位を K にする。
    温度計 c の感度   h_c = Φ'[c, :2]           （その場所の温度が a1, a2 でどう変わるか）
    変位計 n,成分 の感度 h = D'[n, 成分, :2]     （FrontISTR で前もって解いた、型ごとの変位 [µm/K]）
    反り A−O の係数    w = D'[A, z, :2] − D'[O, z, :2]
  測る前の不確かさ B：本番（run/placement_method_comparison.py）と同じ 4,000 本の事前アンサンブル（初期温度・発熱・放熱を
  ばらつかせて ROM で進めたもの）の t = 150 秒の散らばりを、型 1・2 の重みに直した共分散。
  ノイズ：温度計 0.3 K、変位計 0.3 µm。温度計は P2・P0（本番と同じ）。

  出力：docs/img/blog012_fisher_1d.png、blog012_fisher_ellipse.png、blog012_fisher_score_map.png、
        results/blog012_fisher.json
  再現：OMP_NUM_THREADS=4 python3 run/make_blog012_fisher_figs.py
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p  # noqa: F401
from dacore import rom_general as rg
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
K = 273.15; SIG_T = 0.30; SIG_U = 0.30; TN = [2, 0]; COMP = ["Ux", "Uy", "Uz"]
A_XYZ = (0.028, 0.0, 0.1005); O_XYZ = (-0.028, 0.0, 0.1005); EXCL_R = 0.012; T_EVAL = 150.0


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], 5); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); D = op["D"]; coords = op["coords"]
    s = np.abs(U).max(0)                                   # 型ごとの「一番大きく変わるセル」の値
    Us = U / s; Ds = D / s                                # 大きさをそろえた型（a の単位が K になる）
    near = lambda x: int(np.linalg.norm(coords - np.array(x), axis=1).argmin())
    iA, iO = near(A_XYZ), near(O_XYZ)
    w = Ds[iA, 2, :2] - Ds[iO, 2, :2]                     # 反り A−O = w・a  [µm]

    # ---- 測る前の不確かさ B（本番と同じ事前アンサンブル、t = 150 秒）----
    rng = np.random.default_rng(20260930)
    Tp = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (4000, 5))
    qp = rng.uniform(0.3, 1.8, 4000); hp = np.clip(rng.normal(0.02, 0.01, 4000), 1e-3, 0.1)
    Tp = rg.integrate_ensemble(Tp, C, Km, hp, qp, heat, 0.0, T_EVAL, 2.0)
    a_ens = (Tp - mean[pod]) @ np.linalg.inv(U[pod]).T * s  # 5 つの重み（K）
    a2 = a_ens[:, :2]; B = np.cov(a2.T); a0 = a2.mean(0)

    # 正解（OpenFOAM、t = 150 秒）の重み
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); kt = int(np.argmin(abs(z["times"] - T_EVAL)))
    a_true = ((z["Tfield"][kt] - mean) @ U * s)[:2]

    def post(Bm, hs, rs):
        F = sum(np.outer(h, h) / r for h, r in zip(hs, rs)); return np.linalg.inv(np.linalg.inv(Bm) + F), F

    hT = [Us[pod[j], :2] for j in TN]
    PT, FT = post(B, hT, [SIG_T ** 2] * 2)
    sig = lambda P: float(np.sqrt(w @ P @ w))

    # ---- 候補（表面の全節点 × 3 成分、A・O の近くを除く）----
    far = (np.linalg.norm(coords - np.array(A_XYZ), axis=1) > EXCL_R) & (np.linalg.norm(coords - np.array(O_XYZ), axis=1) > EXCL_R)
    H = Ds[:, :, :2].reshape(-1, 2); cand = np.where(np.repeat(far, 3))[0]
    r = SIG_U ** 2
    q = np.einsum("ij,jk,ik->i", H, PT, H)                 # h^T P h
    goal = (H @ PT @ w) ** 2 / (q + r)                     # 目的指向の得（反りの分散の減り）
    Aopt = np.einsum("ij,jk,ik->i", H, PT @ PT, H) / (q + r)
    Dopt = np.log(1 + q / r)
    big = np.linalg.norm(H, axis=1)
    pick = lambda v: int(cand[np.argmax(v[cand])])
    top = cand[big[cand] >= np.percentile(big[cand], 90)]          # よく動く上位 10 % の候補
    bad = int(top[np.argmin(goal[top])])                              # その中で反りに一番効かない点
    sel = {"目的指向": pick(goal), "A 最適": pick(Aopt), "D 最適": pick(Dopt), "よく動く点": pick(big), "よく動くが反りに効かない点": bad}
    out = {"t_eval_s": T_EVAL, "scale_s": s[:2].tolist(), "B": B.tolist(), "a_prior_mean": a0.tolist(), "a_true": a_true.tolist(),
           "w": w.tolist(), "hT": [h.tolist() for h in hT], "FT": FT.tolist(), "PT": PT.tolist(),
           "sigma_AO_prior": sig(B), "sigma_AO_temp": sig(PT), "picks": {}}
    for nm, k in sel.items():
        h = H[k]; P1, _ = post(PT, [h], [r])
        out["picks"][nm] = dict(xyz_mm=(coords[k // 3] * 1000).round(1).tolist(), comp=COMP[k % 3], h=h.tolist(),
                                goal=float(goal[k]), sigma_AO=sig(P1), trace=float(np.trace(P1)), P=P1.tolist())
    # 1 次元（型 1 だけ）の例：P2 の温度計 1 本
    b1 = B[0, 0]; h1 = Us[pod[2], 0]; f1 = h1 ** 2 / SIG_T ** 2; p1 = 1 / (1 / b1 + f1)
    yP2 = z["Tfield"][kt][pod[2]] - K
    a1_hat = p1 * (a0[0] / b1 + h1 * (yP2 - (mean[pod[2]] - K) - Us[pod[2], 1] * a0[1]) / SIG_T ** 2)
    out["one_d"] = dict(b=b1, h=h1, F=f1, P=p1, y=yP2, a1_prior=a0[0], a1_post=a1_hat, a1_true=a_true[0])
    json.dump(out, open(os.path.join(RES, "blog012_fisher.json"), "w"), ensure_ascii=False, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k in ("B", "a_prior_mean", "a_true", "w", "hT", "FT", "PT", "sigma_AO_prior", "sigma_AO_temp", "one_d")}, ensure_ascii=False))
    for nm, v in out["picks"].items(): print(nm, v["xyz_mm"], v["comp"], "h", np.round(v["h"], 3), "得", round(v["goal"], 4), "σAO", round(v["sigma_AO"], 3))

    plt.rcParams.update({"font.size": 16, "axes.labelsize": 18, "xtick.labelsize": 16, "ytick.labelsize": 16, "legend.fontsize": 15})
    # ---- 図 1：1 次元 ----
    x = np.linspace(-3.2 * np.sqrt(b1), 3.2 * np.sqrt(b1), 600)
    g = lambda v: np.exp(-x ** 2 / (2 * v)) / np.sqrt(2 * np.pi * v)
    fig, ax = plt.subplots(figsize=(16, 6.2))
    ax.plot(x, g(b1), "--", color="#555", lw=2.5, label=f"測る前（分散 B = {b1:.2f}）：標準偏差 {np.sqrt(b1):.2f} K")
    ax.plot(x, g(1 / f1), "-", color="#2E8B57", lw=2.5, label=f"P2 の温度計だけから（分散 1/F = {1 / f1:.2f}）：標準偏差 {np.sqrt(1 / f1):.2f} K")
    ax.plot(x, g(p1), "-", color="#C0392B", lw=3.2, label=f"合わせた後（分散 P = {p1:.3f}）：標準偏差 {np.sqrt(p1):.2f} K")
    ax.set_xlabel("型 1 の重み $a_1$ のずれ [K]（幅を比べるため、中心を 0 にそろえた）"); ax.set_ylabel("確率密度"); ax.set_yticks([])
    ax.set_ylim(0, g(p1).max() * 1.55); ax.legend(loc="upper left", frameon=False); ax.grid(alpha=.3)
    fig.savefig(os.path.join(IMG, "blog012_fisher_1d.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)

    # ---- 図 2：2 次元の不確かさの楕円 ----
    def ell(ax, P, c, ls, lab, lw=2.5):
        vals, vecs = np.linalg.eigh(P); ang = np.degrees(np.arctan2(vecs[1, 1], vecs[0, 1]))
        ax.add_patch(Ellipse(a0, 2 * np.sqrt(vals[1]), 2 * np.sqrt(vals[0]), angle=ang, fill=False, ec=c, ls=ls, lw=lw, label=lab))
    fig, axs = plt.subplots(1, 2, figsize=(16, 7.4), gridspec_kw=dict(wspace=0.28))
    for ax, zoom in zip(axs, [False, True]):
        ell(ax, B, "#555", "--", f"測る前 B（反りの幅 ±{sig(B):.2f} µm）")
        ell(ax, PT, "#2E8B57", "-", f"温度計 P2・P0 で測った後（±{sig(PT):.2f} µm）")
        cs = {"目的指向": "#C0392B", "よく動くが反りに効かない点": "#1f77b4"}
        for nm in cs:
            ell(ax, np.array(out["picks"][nm]["P"]), cs[nm], "-", f"＋変位計 1 本：{nm}（±{out['picks'][nm]['sigma_AO']:.2f} µm）", lw=3)
        # 反りが一定の線（w・a = 一定）
        L = 3.2 * np.sqrt(np.linalg.eigvalsh(B).max()) if not zoom else 1.2 * np.sqrt(PT.diagonal().max())
        n = np.array([-w[1], w[0]]) / np.linalg.norm(w)
        for kk in np.linspace(-1, 1, 5):
            c0 = a0 + kk * L * 0.5 * w / np.linalg.norm(w)
            ax.plot(*np.array([c0 - 3 * L * n, c0 + 3 * L * n]).T, color="#999", lw=0.8, ls=":")
        ax.annotate("", xy=a0 + 0.45 * L * w / np.linalg.norm(w), xytext=a0, arrowprops=dict(arrowstyle="-|>", color="k", lw=2))
        ax.text(*(a0 + 0.47 * L * w / np.linalg.norm(w) + np.array([0.06 * L, 0])), "反りが増える向き", fontsize=15, ha="left", va="center",
                bbox=dict(fc="white", ec="none", pad=1))
        ax.set_xlim(a0[0] - L, a0[0] + L); ax.set_ylim(a0[1] - L, a0[1] + L); ax.set_aspect("equal")
        ax.set_xlabel("型 1 の重み $a_1$ [K]"); ax.set_ylabel("型 2 の重み $a_2$ [K]"); ax.grid(alpha=.3)
        ax.set_title("全体" if not zoom else "拡大（温度計で測った後）", fontsize=19)
    h_, l_ = axs[0].get_legend_handles_labels()
    fig.legend(h_, l_, loc="upper center", bbox_to_anchor=(0.5, 1.07), ncol=2, frameon=False)
    fig.savefig(os.path.join(IMG, "blog012_fisher_ellipse.png"), dpi=110, bbox_inches="tight", facecolor="white"); plt.close(fig)

    # ---- 図 3：目的指向の得の分布（表面、3 成分のうち一番大きいもの）----
    import pyvista as pv
    sys.path.insert(0, HERE); from make_unobserved_points_fig import grid
    import vtk
    pv.OFF_SCREEN = True
    g3 = grid(); surf = g3.extract_surface(algorithm="dataset_surface"); oid = surf["vtkOriginalPointIds"]
    gm = goal.reshape(-1, 3).max(1); gm[~far] = 0
    surf.point_data["g"] = gm[oid] / gm.max()
    pl = pv.Plotter(off_screen=True, window_size=(1000, 950)); pl.set_background("white")
    pl.add_mesh(surf, scalars="g", cmap="viridis", clim=[0, 1], show_scalar_bar=False, lighting=False)
    pl.add_mesh(surf.extract_feature_edges(), color="#333333", line_width=1.5)
    marks = {"A": (np.array(A_XYZ) * 1000, "black"), "O": (np.array(O_XYZ) * 1000, "black"),
             "P2": (kv["xyz"][2] * 1000, "#C0392B"), "P0": (kv["xyz"][0] * 1000, "#C0392B"),
             "best": (coords[sel["目的指向"] // 3] * 1000, "#E8452C"), "move": (coords[sel["よく動くが反りに効かない点"] // 3] * 1000, "#1f77b4")}
    for nm, (xx, c) in marks.items(): pl.add_mesh(pv.Sphere(radius=2.6, center=xx), color=c)
    pl.camera_position = [(230, -260, 230), (0, 0, 50), (0, 0, 1)]; pl.camera.zoom(1.25)
    img = pl.screenshot(return_img=True)
    xy = {}
    for nm, (xx, _) in marks.items():
        cc = vtk.vtkCoordinate(); cc.SetCoordinateSystemToWorld(); cc.SetValue(*map(float, xx))
        u, v = cc.GetComputedDisplayValue(pl.renderer); xy[nm] = (u, 950 - v)
    pl.close()
    fig = plt.figure(figsize=(12, 9.2)); ax = fig.add_axes([0, 0.1, 1, 0.9]); ax.imshow(img); ax.axis("off")
    lab = {"A": "A（反りの評価点）", "O": "O（反りの評価点）", "P2": "P2 温度計", "P0": "P0 温度計（内面の近く）",
           "best": f"目的指向で選んだ点（{out['picks']['目的指向']['comp']}）", "move": f"よく動くが反りに効かない点（{out['picks']['よく動くが反りに効かない点']['comp']}）"}
    off = dict(A=(40, -150), O=(-150, -40), P2=(170, 30), P0=(-200, 60), best=(250, -60), move=(220, 30))
    for nm, (u, v) in xy.items():
        ax.annotate(lab[nm], (u, v), xytext=(u + off[nm][0], v + off[nm][1]), fontsize=16, ha="center", va="center",
                    arrowprops=dict(arrowstyle="-", color="black", lw=1.2), bbox=dict(fc="white", ec="none", pad=1))
    cax = fig.add_axes([0.2, 0.06, 0.6, 0.025])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap="viridis"), cax=cax, orientation="horizontal")
    cb.set_label("変位計を 1 本置いたときの、反りの分散の減り（一番大きい所を 1）", fontsize=16)
    fig.savefig(os.path.join(IMG, "blog012_fisher_score_map.png"), dpi=100, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote figs")


if __name__ == "__main__":
    main()
