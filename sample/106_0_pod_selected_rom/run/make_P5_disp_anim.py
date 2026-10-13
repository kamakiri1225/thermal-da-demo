"""代表点5点の温度から変位場ができる様子のアニメーション（GIF）.

変位（全節点・3方向）は、代表点5点の温度だけで次のように書ける（blog_008・011 の H_u）：
  u(t) = Σ_j H_u[:, :, j] (T_Pj(t) − 20 °C) + c0
  H_u = D (Φ_P)^{-1}   （D：POD の型を FrontISTR に入れた変位、Φ_P：代表5セルの POD の型）
  c0 = u_mean − H_u (mean_P − 20 °C)   （定数。大きさ 0.001 µm で、ほぼ 0）
各項 H_u[:, :, j] (T_Pj − 20 °C) が「代表点 Pj の温度から来る変位」。5つを足すと変位場になる。

温度 T_Pj(t) は OpenFOAM の真値（15 W、代表5セル）。比較に FrontISTR の真値（全温度場から計算）を並べる。
表示：変形を拡大した形（全パネル同じ倍率）を、U_z の色で塗る。

出力: docs/img/P5_disp_anim.gif, docs/img/P5_disp_t300.png, results/P5_disp_anim.json
再現: OMP_NUM_THREADS=4 python3 run/make_P5_disp_anim.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pyvista as pv
import cylinder_mesh
from make_W_AO_vtk_anim import FixedPaletteWriter
pv.OFF_SCREEN = True

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NR, NTH, NZ = 4, 48, 20; R_IN, R_OUT, H_ = 0.020, 0.0375, 0.1005
TREF = 293.15; SCALE = 2000.0          # 変形の拡大倍率（µm → mm の表示で ×2000）
CMAP = "RdBu_r"
NAMES = ["P0（ヒータ側・内面・中段）", "P1（+y側・中段）", "P2（ヒータ側・中段）", "P3（ヒータ側・下端）", "P4（反対側・下端）"]


def main():
    mesh = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, R_IN, R_OUT, H_)
    ids = [n for n, _ in mesh["nodes"]]; coords = np.array([x for _, x in mesh["nodes"]])
    pos = {n: i for i, n in enumerate(ids)}
    cells = np.array([[8] + [pos[n] for n in conn] for _, conn in mesh["elements"]]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(mesh["elements"]), pv.CellType.HEXAHEDRON), coords * 1000)
    surf = grid.extract_surface(algorithm="dataset_surface")
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); assert np.abs(op["coords"] - coords).max() < 1e-12
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; mean = kv["mean"]
    Pxyz = kv["xyz"] * 1000
    Hu = np.einsum("nck,kj->ncj", op["D"], np.linalg.inv(U[pod]))              # (5040,3,5) µm/K
    c0 = op["u_mean"] - np.einsum("ncj,j->nc", Hu, mean[pod] - TREF)
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; TP = z["Tfield"][:, pod]   # (21,5)
    dT = TP - TREF
    contrib = np.einsum("ncj,tj->jtnc", Hu, dT)                                 # (5,21,5040,3)
    total = contrib.sum(0) + c0[None]
    truth = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]
    iA = int(np.linalg.norm(coords - np.array((0.028, 0, H_)), axis=1).argmin())
    iO = int(np.linalg.norm(coords - np.array((-0.028, 0, H_)), axis=1).argmin())
    out = {"note": "u = Σ_j H_u[:,:,j](T_Pj−20°C) + c0。T_Pj は OpenFOAM 真値（15 W）。µm。",
           "c0_max_um": float(np.abs(c0).max()), "max_abs_total_minus_FrontISTR_um": float(np.abs(total - truth).max()),
           "t": t.tolist(), "T_P_degC": (TP - 273.15).round(3).tolist(),
           "UzA_contrib_um": contrib[:, :, iA, 2].round(4).tolist(), "UzO_contrib_um": contrib[:, :, iO, 2].round(4).tolist(),
           "UzA_total": total[:, iA, 2].round(4).tolist(), "UzA_FrontISTR": truth[:, iA, 2].round(4).tolist(),
           "UzO_total": total[:, iO, 2].round(4).tolist(), "UzO_FrontISTR": truth[:, iO, 2].round(4).tolist()}
    json.dump(out, open(os.path.join(RES, "P5_disp_anim.json"), "w"), ensure_ascii=False, indent=1)
    print("c0 max", out["c0_max_um"], "total-truth max", out["max_abs_total_minus_FrontISTR_um"], flush=True)

    lim = float(np.percentile(np.abs(np.concatenate([contrib[..., 2].ravel(), truth[..., 2].ravel()])), 99.5))
    cam = [(190, -230, 210), (0, 0, 50), (0, 0, 1)]

    def render(u, hi=None):
        p = pv.Plotter(off_screen=True, window_size=(640, 640))
        s = surf.copy(); s.points = surf.points + u[surf["vtkOriginalPointIds"]] * 1e-3 * SCALE
        s.point_data["uz"] = u[surf["vtkOriginalPointIds"], 2]
        p.add_mesh(surf, style="wireframe", color="#BBBBBB", opacity=0.25, line_width=1)
        p.add_mesh(s, scalars="uz", cmap=CMAP, clim=[-lim, lim], show_scalar_bar=False)
        for j in range(5):
            col = "#E31A1C" if j == hi else "#333333"
            p.add_point_labels([Pxyz[j]], [f"P{j}"], font_size=30 if j == hi else 20, text_color=col, shape=None,
                               show_points=True, point_size=18 if j == hi else 9, point_color=col, always_visible=True)
        p.add_point_labels([(R_OUT * 1000 + 16, 0, 50)], ["HEATER"], font_size=22, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.2)
        img = p.screenshot(return_img=True); p.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2))
        r0, r1 = max(nz[0].min() - 8, 0), min(nz[0].max() + 8, img.shape[0])
        c_0, c_1 = max(nz[1].min() - 8, 0), min(nz[1].max() + 8, img.shape[1])
        return np.ascontiguousarray(img[r0:r1, c_0:c_1, :3])

    frames = []
    for k in range(len(t)):
        frames.append([render(contrib[j, k], hi=j) for j in range(5)] + [render(total[k]), render(truth[k])])
        print("rendered", k, flush=True)

    plt.rcParams.update({"font.size": 17})
    fig = plt.figure(figsize=(24, 13.5))
    gs = fig.add_gridspec(2, 4, hspace=0.28, wspace=0.04, top=0.86, bottom=0.12, left=0.01, right=0.99)
    fig.suptitle("代表点5点の温度から、変位がどうできるか（15 W、OpenFOAM の真値の温度を使用）", fontsize=28, weight="bold")
    fig.text(0.5, 0.905, "変位 ＝ Σ（Pj の重み H_u）×（Pj の温度 − 20 ℃）。各パネルは1つの項。5つを足すと「合計」になる。"
             f"変形は {SCALE:.0f} 倍に拡大（灰色の線＝変形前）。色＝上下変位 Uz。", ha="center", fontsize=18, color="#444")
    axs = [fig.add_subplot(gs[i // 4, i % 4]) for i in range(7)]
    for a in axs: a.axis("off")
    ims = [axs[i].imshow(frames[0][i]) for i in range(7)]
    axH = fig.add_subplot(gs[1, 3]); axH.axvspan(0, 300, color="#FDEBD0", alpha=.6)
    cols = ["#1f77b4", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]
    for j in range(5): axH.plot(t, contrib[j, :, iA, 2], color=cols[j], lw=2, label=f"P{j} から")
    axH.plot(t, total[:, iA, 2], color="k", lw=3.5, label="5点の合計")
    axH.plot(t, truth[:, iA, 2], "o", color="#E1A100", mfc="none", ms=9, mew=2, label="FrontISTR 真値")
    axH.axhline(0, color="#999", lw=1); cur = axH.axvline(0, color="#555", lw=2)
    axH.set_title("点 A（ヒータ側上面）の Uz の内訳", fontsize=19); axH.set_xlabel("時刻 [s]"); axH.set_ylabel("Uz(A) [µm]")
    axH.legend(fontsize=12, ncol=2, loc="upper left"); axH.grid(alpha=.3)
    cax = fig.add_axes([0.30, 0.045, 0.40, 0.022])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(-lim, lim), cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_label("上下変位 Uz [µm]（全パネル共通）", fontsize=18)
    tk = fig.text(0.01, 0.93, "", fontsize=22, weight="bold")

    def draw(k):
        for i in range(7): ims[i].set_data(frames[k][i])
        for j in range(5):
            axs[j].set_title(f"{NAMES[j]} から\nT = {TP[k, j] - 273.15:.1f} ℃（+{dT[k, j]:.1f} K）　Uz(A) {contrib[j, k, iA, 2]:+.2f} µm", fontsize=16)
        axs[5].set_title(f"合計（5点の温度から）\nUz(A) {total[k, iA, 2]:.2f} µm", fontsize=17, weight="bold")
        axs[6].set_title(f"FrontISTR 真値（全温度場から）\nUz(A) {truth[k, iA, 2]:.2f} µm", fontsize=17)
        cur.set_xdata([t[k], t[k]]); tk.set_text(f"t = {t[k]:.0f} 秒" + ("（加熱中）" if t[k] <= 300 else "（加熱停止後）"))
        return []
    anim = FuncAnimation(fig, draw, frames=len(t), interval=700)
    anim.save(os.path.join(IMG, "P5_disp_anim.gif"), writer=FixedPaletteWriter(fps=1.5), dpi=60)
    draw(10); fig.savefig(os.path.join(IMG, "P5_disp_t300.png"), dpi=75, facecolor="white"); plt.close(fig)
    print("wrote gif", os.path.getsize(os.path.join(IMG, "P5_disp_anim.gif")) // 1024, "KB")


if __name__ == "__main__":
    main()
