"""熱感度 W の A・O の行（全節点の温度に対する感度）を VTK に出し、温度分布の時刻歴と掛け算したアニメーションを作る.

W の行は、FrontISTR がダンプした剛性行列 K と熱荷重行列 H（105 の dumpw_case）から随伴法で求める
（W[i,:] = (K^{-1} e_i)^T H。run/W_rows_by_adjoint.py と同じ）。
温度は OpenFOAM の真値（15 W、21 時刻）を、各節点に一番近いセルから渡す。

  U_z(A)(t) = Σ_j W_A,j (T_j(t) − T_ref)       （T_ref = 20 °C）

出力
  results/vtk/W_AO.vtu                    W_A, W_O, W_A−O（µm/K）
  results/vtk/contrib_XXX.vtu, contrib.pvd 時刻ごとの ΔT、W_A×ΔT、W_O×ΔT、(W_A−W_O)×ΔT（ParaView で時系列として開ける）
  docs/img/W_AO_times_temp_anim.gif       アニメーション
  results/W_AO_times_temp.json            掛け算の合計と FrontISTR の比較
再現: OMP_NUM_THREADS=4 python3 run/make_W_AO_vtk_anim.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import scipy.io as sio
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import pyvista as pv
import cylinder_mesh
pv.OFF_SCREEN = True

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img"); VTK = os.path.join(RES, "vtk")
CASE = os.path.join(SAMPLE, "105_0_sensor_placement_sensitivity", "openfoam", "dumpw_case")
NR, NTH, NZ = 4, 48, 20; R_IN, R_OUT, H_ = 0.020, 0.0375, 0.1005
TREF = 293.15


def main():
    os.makedirs(VTK, exist_ok=True)
    mesh = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, R_IN, R_OUT, H_)
    ids = [n for n, _ in mesh["nodes"]]; coords = np.array([x for _, x in mesh["nodes"]])
    pos = {n: i for i, n in enumerate(ids)}
    cells = np.array([[8] + [pos[n] for n in conn] for _, conn in mesh["elements"]]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(mesh["elements"]), pv.CellType.HEXAHEDRON), coords * 1000)
    iA = int(np.linalg.norm(coords - np.array((0.028, 0, H_)), axis=1).argmin())
    iO = int(np.linalg.norm(coords - np.array((-0.028, 0, H_)), axis=1).argmin())
    # ---- 随伴法で W の A・O の行 ----
    K = sio.mmread(os.path.join(CASE, "dump_matrix_1_0.mm")).tocsc().astype(np.float64)
    Hm = sio.mmread(os.path.join(CASE, "H_matrix.mtx")).tocsc().astype(np.float64)
    lu = splu(K); E = np.zeros((K.shape[0], 2)); E[3 * iA + 2, 0] = 1; E[3 * iO + 2, 1] = 1
    Wr = (lu.solve(E).T @ Hm) * 1e6                       # (2, 5040) µm/K
    WA, WO = Wr[0], Wr[1]; WAO = WA - WO
    print(f"sum W_A = {WA.sum():.3f} µm/K, sum W_O = {WO.sum():.3f}, sum W_A-O = {WAO.sum():.4f}", flush=True)
    g = grid.copy(); g.point_data["W_A_um_per_K"] = WA; g.point_data["W_O_um_per_K"] = WO; g.point_data["W_AminusO_um_per_K"] = WAO
    g.save(os.path.join(VTK, "W_AO.vtu"))
    # ---- 温度の時刻歴（真値）と掛け算 ----
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); Cc = kv["cell_centres"]
    _, near = cKDTree(Cc).query(coords)
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; Tn = z["Tfield"][:, near]   # (21, 5040)
    dT = Tn - TREF
    cA, cO = dT * WA, dT * WO
    sA, sO = cA.sum(1), cO.sum(1)
    ut = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]
    trA, trO = ut[:, iA, 2], ut[:, iO, 2]
    pvd = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1"><Collection>']
    for k, tk in enumerate(t):
        gk = grid.copy()
        gk.point_data["T_K"] = Tn[k]; gk.point_data["dT_K"] = dT[k]
        gk.point_data["WA_x_dT_um"] = cA[k]; gk.point_data["WO_x_dT_um"] = cO[k]; gk.point_data["WAO_x_dT_um"] = cA[k] - cO[k]
        fn = f"contrib_{k:03d}.vtu"; gk.save(os.path.join(VTK, fn))
        pvd.append(f'<DataSet timestep="{tk:g}" file="{fn}"/>')
    pvd.append("</Collection></VTKFile>"); open(os.path.join(VTK, "contrib.pvd"), "w").write("\n".join(pvd))
    out = {"note": "U_z = Σ W·(T−20°C)。W は 105 のダンプ（四面体分割メッシュ）から随伴法、真値は FrontISTR（六面体）に OpenFOAM の全温度場。",
           "sum_W_A": float(WA.sum()), "sum_W_O": float(WO.sum()),
           "t": t.tolist(), "UzA_W": sA.round(4).tolist(), "UzA_FrontISTR": trA.round(4).tolist(),
           "UzO_W": sO.round(4).tolist(), "UzO_FrontISTR": trO.round(4).tolist(),
           "max_abs_diff_UzA": float(np.abs(sA - trA).max()), "max_abs_diff_AO": float(np.abs((sA - sO) - (trA - trO)).max())}
    json.dump(out, open(os.path.join(RES, "W_AO_times_temp.json"), "w"), ensure_ascii=False, indent=1)
    print("max diff UzA", out["max_abs_diff_UzA"], "AO", out["max_abs_diff_AO"], flush=True)

    # ---- アニメーション：pyvista で3つの円筒、下に時刻歴 ----
    surf = grid.extract_surface(algorithm="dataset_surface") if hasattr(grid, "extract_surface") else grid
    lim_w = np.percentile(np.abs(WA), 99.5); lim_t = float(dT.max()); lim_c = np.percentile(np.abs(cA), 99.8)
    cam = [(30, -255, 120), (0, 0, 50), (0, 0, 1)]

    def render(scal, name, clim, cmap, title):
        p = pv.Plotter(off_screen=True, window_size=(620, 560))
        gg = grid.copy(); gg.point_data[name] = scal
        p.add_mesh(gg, scalars=name, cmap=cmap, clim=clim, show_scalar_bar=True,
                   scalar_bar_args=dict(title=title, vertical=False, title_font_size=18, label_font_size=15,
                                        position_x=0.15, position_y=0.04, width=0.7, height=0.08, fmt="%.3g"))
        p.add_point_labels([coords[iA] * 1000 + [0, 0, 6], coords[iO] * 1000 + [0, 0, 6]], ["A", "O"], font_size=26,
                           text_color="black", shape=None, show_points=True, point_size=14, point_color="black", always_visible=True)
        p.add_point_labels([(R_OUT * 1000 + 14, 0, 50)], ["HEATER"], font_size=18, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.05)
        img = p.screenshot(return_img=True); p.close(); return img

    # 色の範囲は、支持点（底面3節点）と A・O のすぐ近くを除いて決める（点拘束・点評価による局所的な値のため）
    fixn = [pos[n] for n in mesh["fix_xyz"] + mesh["fix_yz"] + mesh["fix_z"]]
    dist = lambda idx: np.linalg.norm(coords - coords[idx], axis=1)
    mask = np.ones(len(coords), bool)
    for f in fixn: mask &= dist(f) > 0.015
    for f in (iA, iO): mask &= dist(f) > 0.008
    rob = lambda v: float(np.percentile(np.abs(v[..., mask]), 98))
    rows = [("A の上下変位 U_z(A)", WA, cA, sA, trA, "W_A"),
            ("O の上下変位 U_z(O)", WO, cO, sO, trO, "W_O"),
            ("反り A−O", WAO, cA - cO, sA - sO, trA - trO, "W_A - W_O")]
    statics, anims = [], []
    for name, Wv, cv, _, _, lab in rows:
        lw = rob(Wv); lc = rob(cv)
        statics.append(render(Wv, "w", [-lw, lw], "RdBu_r", f"{lab} [um/K]"))
        anims.append([render(cv[k], "c", [-lc, lc], "RdBu_r", f"{lab} x (T-20C) [um]") for k in range(len(t))])
    fig = plt.figure(figsize=(19, 17.5))
    gs = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1.15], hspace=0.18, wspace=0.05, top=0.90, bottom=0.06, left=0.02, right=0.98)
    fig.suptitle("左：熱感度 W（その節点が 1 K 上がると何 µm 動くか）　右：W × 温度の上がり分（その節点の温度が実際にどれだけ効いているか）",
                 fontsize=17, weight="bold")
    txt = fig.text(0.5, 0.925, "", ha="center", fontsize=15)
    fig.text(0.5, 0.005, "色の範囲は、底面の支持点（3節点）と A・O の真上の局所的な値を除いて決めた（その付近は色が飽和）。右手前がヒータ側（+x）、左奥が反対側（−x）。",
             ha="center", fontsize=12, color="#555")
    axL = [fig.add_subplot(gs[r, 0]) for r in range(3)]; axR = [fig.add_subplot(gs[r, 1]) for r in range(3)]
    axT = [fig.add_subplot(gs[r, 2]) for r in range(3)]; curs = []
    for r, (name, Wv, cv, sv, tv, lab) in enumerate(rows):
        axL[r].imshow(statics[r]); axL[r].axis("off")
        axL[r].set_title(f"{name}：熱感度 {lab}（時間で変わらない）", fontsize=13)
        a = axT[r]; a.axvspan(0, 300, color="#FDEBD0", alpha=.5)
        a.plot(t, tv, "o", color="#C0392B", mfc="none", ms=7, label="FrontISTR（全温度場）")
        a.plot(t, sv, "-", color="#14213D", lw=2.3, label="右の図を全節点で足したもの")
        curs.append(a.axvline(0, color="#555", lw=2)); a.grid(alpha=.3); a.set_ylabel("変位 [µm]")
        a.set_title(f"{name} の時刻歴", fontsize=13)
        if r == 0: a.legend(fontsize=11, loc="lower right")
        if r == 2: a.set_xlabel("時刻 [s]")

    def draw(k):
        for r, (name, Wv, cv, sv, tv, lab) in enumerate(rows):
            axR[r].clear(); axR[r].imshow(anims[r][k]); axR[r].axis("off")
            axR[r].set_title(f"{name}：{lab} × 温度の上がり分（t = {t[k]:.0f} 秒）", fontsize=13)
            curs[r].set_xdata([t[k], t[k]])
        txt.set_text(f"t = {t[k]:.0f} 秒　　U_z(A) = {sA[k]:.2f} µm（FrontISTR {trA[k]:.2f}）　"
                     f"U_z(O) = {sO[k]:.2f} µm（{trO[k]:.2f}）　反り A−O = {sA[k]-sO[k]:.2f} µm（{trA[k]-trO[k]:.2f}）")
        return []
    anim = FuncAnimation(fig, draw, frames=len(t), interval=700)
    anim.save(os.path.join(IMG, "W_AO_times_temp_anim.gif"), writer=PillowWriter(fps=1.5), dpi=44)
    draw(10); fig.savefig(os.path.join(IMG, "W_AO_times_temp_t300.png"), dpi=70, facecolor="white"); plt.close(fig)
    print("wrote gif", os.path.getsize(os.path.join(IMG, "W_AO_times_temp_anim.gif")) // 1024, "KB")


if __name__ == "__main__":
    main()
