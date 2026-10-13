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
from dacore.temperature_mapping import build_idw_mapping, map_temperature, mapping_cache_valid, MAPPING_METHOD
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



class FixedPaletteWriter(PillowWriter):
    """固定したWの色がGIFのフレーム別量子化で変化しないようにする。"""
    def finish(self):
        from PIL import Image
        rgb = [frame.convert("RGB") for frame in self._frames]
        # 256 色の表は全コマを並べた画像から作る（最初のコマだけだと、後のコマの色が化ける）。
        # 表は全コマ共通なので、固定した W の色はコマ間で変わらない。
        w, h = rgb[0].size
        mosaic = Image.new("RGB", (w, h * len(rgb)))
        for i, frame in enumerate(rgb): mosaic.paste(frame, (0, h * i))
        palette = mosaic.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
        frames = [frame.quantize(palette=palette, dither=Image.Dither.FLOYDSTEINBERG) for frame in rgb]
        frames[0].save(self.outfile, save_all=True, append_images=frames[1:],
                       duration=int(1000 / self.fps), loop=0, optimize=False)


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
    # DUMPW=YES が指定した A/O 節点から直接出力した感度差を照合する。
    point_file = os.path.join(CASE, "sensitivity_points.dat")
    configured_ids = [int(x) for line in open(point_file) if not line.lstrip().startswith("#") for x in line.split()]
    assert configured_ids == [ids[iA], ids[iO]], (configured_ids, ids[iA], ids[iO])
    native = os.path.join(CASE, "sensitivity_Wdiff.vtk")
    if not os.path.exists(native):
        native = os.path.join(SAMPLE, "105_0_sensor_placement_sensitivity", "paraview", "sensitivity_Wdiff_cylinder.vtk")
    ng = pv.read(native)
    distance, native_index = cKDTree(ng.points).query(coords)
    assert distance.max() < 1e-8, distance.max()
    native_diff = np.asarray(ng.point_data["Sensitivity"])[native_index, 2] * 1e6
    native_rel_error = float(np.linalg.norm(WAO - native_diff) / np.linalg.norm(native_diff))
    print(f"FrontISTR native Wdiff_z vs Python: relative error={native_rel_error:.3g}; A/O node IDs={configured_ids}", flush=True)
    assert native_rel_error < 1e-5, native_rel_error
    ng.save(os.path.join(VTK, "FrontISTR_sensitivity_Wdiff.vtu"))
    # A/O個別の行はPython随伴法、差の行はFrontISTR直接出力を使う。
    WAO = native_diff
    print(f"sum W_A = {WA.sum():.3f} µm/K, sum W_O = {WO.sum():.3f}, sum W_A-O = {WAO.sum():.4f}", flush=True)
    g = grid.copy(); g.point_data["W_A_um_per_K"] = WA; g.point_data["W_O_um_per_K"] = WO; g.point_data["W_AminusO_um_per_K"] = WAO
    g.save(os.path.join(VTK, "W_AO.vtu"))
    # ---- 温度の時刻歴（真値）と掛け算 ----
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); Cc = kv["cell_centres"]
    temperature_map = build_idw_mapping(Cc, coords)  # 近傍8セルのIDW
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; Tn = map_temperature(temperature_map, z["Tfield"], axis=1)   # (21, 5040)
    dT = Tn - TREF
    cA, cO = dT * WA, dT * WO
    cAO = dT * WAO; sAO = cAO.sum(1)
    sA, sO = cA.sum(1), cO.sum(1)
    ut = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]
    trA, trO = ut[:, iA, 2], ut[:, iO, 2]
    pvd = ['<?xml version="1.0"?>', '<VTKFile type="Collection" version="0.1"><Collection>']
    for k, tk in enumerate(t):
        gk = grid.copy()
        gk.point_data["T_K"] = Tn[k]; gk.point_data["dT_K"] = dT[k]
        gk.point_data["WA_x_dT_um"] = cA[k]; gk.point_data["WO_x_dT_um"] = cO[k]; gk.point_data["WAO_x_dT_um"] = cAO[k]
        fn = f"contrib_{k:03d}.vtu"; gk.save(os.path.join(VTK, fn))
        pvd.append(f'<DataSet timestep="{tk:g}" file="{fn}"/>')
    pvd.append("</Collection></VTKFile>"); open(os.path.join(VTK, "contrib.pvd"), "w").write("\n".join(pvd))
    out = {"note": "U_z = Σ W·(T−20°C)。W は 105 のダンプ（四面体分割メッシュ）から随伴法、真値は FrontISTR（六面体）に OpenFOAM の全温度場。",
           "A_node_id": ids[iA], "O_node_id": ids[iO], "native_Wdiff_source": native,
           "native_Wdiff_relative_error": native_rel_error, "sum_W_A": float(WA.sum()), "sum_W_O": float(WO.sum()),
           "t": t.tolist(), "UzA_W": sA.round(4).tolist(), "UzA_FrontISTR": trA.round(4).tolist(),
           "UzO_W": sO.round(4).tolist(), "UzO_FrontISTR": trO.round(4).tolist(),
           "max_abs_diff_UzA": float(np.abs(sA - trA).max()), "max_abs_diff_AO": float(np.abs((sA - sO) - (trA - trO)).max())}
    json.dump(out, open(os.path.join(RES, "W_AO_times_temp.json"), "w"), ensure_ascii=False, indent=1)
    print("max diff UzA", out["max_abs_diff_UzA"], "AO", out["max_abs_diff_AO"], flush=True)

    # ---- アニメーション：pyvista で3つの円筒、下に時刻歴 ----
    surf = grid.extract_surface(algorithm="dataset_surface") if hasattr(grid, "extract_surface") else grid
    lim_w = np.percentile(np.abs(WA), 99.5); lim_t = float(dT.max()); lim_c = np.percentile(np.abs(cA), 99.8)
    cam = [(180, -220, 200), (0, 0, 50), (0, 0, 1)]

    def render(scal, name, clim, cmap, title, fmt="%.1f"):
        import hashlib
        from PIL import Image
        cache = "/tmp/thermal_W_render_cache"
        os.makedirs(cache, exist_ok=True)
        digest = hashlib.sha256(np.asarray(scal).tobytes() + repr((name, clim, cmap, title, cam, fmt, "separate-legend-v7-nolight")).encode()).hexdigest()
        cache_file = os.path.join(cache, digest + ".png")
        if os.path.exists(cache_file):
            return np.asarray(Image.open(cache_file).convert("RGB"))
        p = pv.Plotter(off_screen=True, window_size=(700, 600))
        gg = grid.copy(); gg.point_data[name] = scal
        p.add_mesh(gg, scalars=name, cmap=cmap, clim=clim, show_scalar_bar=False, lighting=False)   # 陰影なし＝カラーバーと同じ色
        p.add_mesh(gg.extract_surface(algorithm="dataset_surface").extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False, feature_angle=30), color="#333333", line_width=2)
        p.add_point_labels([coords[iA] * 1000 + [0, 0, 4], coords[iO] * 1000 + [0, 0, 4]], ["A", "O"], font_size=42,
                           text_color="black", shape=None, show_points=True, point_size=14, point_color="black", always_visible=True)
        p.add_point_labels([(R_OUT * 1000 + 14, 0, 50)], ["HEATER"], font_size=34, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.15)
        model_img = p.screenshot(return_img=True); p.close()
        # 白い余白を切り取って、円筒を大きく見せる
        nonwhite = np.where((model_img[:, :, :3] < 245).any(axis=2))
        r0, r1 = max(nonwhite[0].min() - 12, 0), min(nonwhite[0].max() + 12, model_img.shape[0])
        c0, c1 = max(nonwhite[1].min() - 12, 0), min(nonwhite[1].max() + 12, model_img.shape[1])
        model_img = np.ascontiguousarray(model_img[r0:r1, c0:c1, :3])
        w = model_img.shape[1]
        if w < 700:                                    # 凡例（幅700）と同じ幅に白で広げる
            pad = np.full((model_img.shape[0], 700 - w, 3), 255, np.uint8)
            model_img = np.concatenate([pad[:, :(700 - w) // 2], model_img, pad[:, (700 - w) // 2:]], axis=1)
        elif w > 700:
            from PIL import Image as _I
            model_img = np.asarray(_I.fromarray(model_img).resize((700, int(model_img.shape[0] * 700 / w))))
        # 凡例専用の画像を下に連結し、モデルとは重ならない領域に描く。
        from matplotlib.colors import Normalize
        from matplotlib.cm import ScalarMappable
        legend_fig = plt.figure(figsize=(7, 2.1), dpi=100, facecolor="white")
        legend_ax = legend_fig.add_axes([0.10, 0.40, 0.80, 0.15])
        cb = legend_fig.colorbar(ScalarMappable(norm=Normalize(*clim), cmap=cmap), cax=legend_ax, orientation="horizontal")
        cb.set_ticks(np.linspace(clim[0], clim[1], 3))
        cb.ax.tick_params(labelsize=34, pad=6)
        from matplotlib.ticker import FormatStrFormatter
        cb.ax.xaxis.set_major_formatter(FormatStrFormatter(fmt))
        cb.ax.set_title(title, fontsize=34, pad=12)
        legend_fig.canvas.draw()
        legend_img = np.asarray(legend_fig.canvas.buffer_rgba())[:, :, :3].copy()
        plt.close(legend_fig)
        img = np.concatenate([model_img[:, :, :3], legend_img], axis=0)
        Image.fromarray(img).save(cache_file)
        return img

    # 色の範囲は、支持点（底面3節点）と A・O のすぐ近くを除いて決める（点拘束・点評価による局所的な値のため）
    fixn = [pos[n] for n in mesh["fix_xyz"] + mesh["fix_yz"] + mesh["fix_z"]]
    dist = lambda idx: np.linalg.norm(coords - coords[idx], axis=1)
    mask = np.ones(len(coords), bool)
    for f in fixn: mask &= dist(f) > 0.015
    for f in (iA, iO): mask &= dist(f) > 0.008
    rob = lambda v: float(np.percentile(np.abs(v[..., mask]), 98))
    rows = [("A の上下変位に対する感度・寄与", WA, cA, sA, trA, "W_A"),
            ("O の上下変位に対する感度・寄与", WO, cO, sO, trO, "W_O"),
            ("反り A−O に対する感度・寄与", WAO, cAO, sAO, trA - trO, "W_A - W_O")]
    statics, anims = [], []
    for name, Wv, cv, _, _, lab in rows:
        lw = rob(Wv); lc = rob(cv)
        statics.append(render(Wv * 1000, "w", [-lw * 1000, lw * 1000], "RdBu_r", f"{lab} [nm/K]", "%.1f"))
        anims.append([render(cv[k] * 1000, "c", [-lc * 1000, lc * 1000], "RdBu_r", f"{lab} x dT [nm]", "%.1f") for k in range(len(t))])
    temp_images = [render(Tn[k] - 273.15, "temperature", [20.0, float(Tn.max() - 273.15)], "rainbow", "Temperature [degC]", "%.0f") for k in range(len(t))]
    plt.rcParams.update({"font.size": 20, "axes.labelsize": 22, "xtick.labelsize": 20, "ytick.labelsize": 20})
    fig = plt.figure(figsize=(23, 16))
    gs = fig.add_gridspec(3, 4, width_ratios=[1, 1, 1, 1.15], hspace=0.32, wspace=0.05, top=0.84, bottom=0.10, left=0.02, right=0.98)
    fig.suptitle("① 温度場 → ② W〈固定〉 → ③ WΔT〈変化〉 → ④ 変位時刻歴",
                 fontsize=30, weight="bold")
    txt = fig.text(0.02, 0.905, "", ha="left", fontsize=23)
    fig.text(0.5, 0.008, "色の範囲は、底面の支持点（3節点）と A・O の真上の局所的な値を除いて決めた（その付近は色が飽和）。\n"
             "斜め上から表示。A＝節点4803（+x側）、O＝節点4923（−x側）。色は温度感度・寄与であり、各位置の変位ではない。",
             ha="center", fontsize=17, color="#555")
    axTemp = fig.add_subplot(gs[:, 0]); axTemp.axis("off")
    temp_artist = axTemp.imshow(temp_images[0])
    axL = [fig.add_subplot(gs[r, 1]) for r in range(3)]; axR = [fig.add_subplot(gs[r, 2]) for r in range(3)]
    axT = [fig.add_subplot(gs[r, 3]) for r in range(3)]; curs = []
    for r, (name, Wv, cv, sv, tv, lab) in enumerate(rows):
        axL[r].imshow(statics[r]); axL[r].axis("off")
        axL[r].set_title(f"② {lab}：温度感度〈固定〉", fontsize=24)
        a = axT[r]; a.axvspan(0, 300, color="#FDEBD0", alpha=.5)
        a.plot(t, tv, "o", color="#C0392B", mfc="none", ms=9, mew=2, label="FrontISTR（全温度場）")
        a.plot(t, sv, "-", color="#14213D", lw=3, label="③を全節点で足したもの")
        curs.append(a.axvline(0, color="#555", lw=2)); a.grid(alpha=.3); a.set_ylabel("変位 [µm]")
        a.set_title(["④ Aの上下変位", "④ Oの上下変位", "④ 反り A−O"][r], fontsize=24)
        if r == 0: fig.legend(*a.get_legend_handles_labels(), loc="upper right", bbox_to_anchor=(0.995, 0.935), fontsize=19, frameon=False, ncol=1)
        if r == 2: a.set_xlabel("時刻 [s]")

    def draw(k):
        temp_artist.set_data(temp_images[k])
        axTemp.set_title(f"① 温度：時間で変化\nt = {t[k]:.0f} 秒（3行に共通）", fontsize=25)
        for r, (name, Wv, cv, sv, tv, lab) in enumerate(rows):
            axR[r].clear(); axR[r].imshow(anims[r][k]); axR[r].axis("off")
            axR[r].set_title(f"③ {lab} × ΔT：寄与〈変化〉", fontsize=24)
            curs[r].set_xdata([t[k], t[k]])
        txt.set_text(f"t = {t[k]:.0f} 秒　　U_z(A) = {sA[k]:.2f} µm（FrontISTR {trA[k]:.2f}）　"
                     f"U_z(O) = {sO[k]:.2f} µm（{trO[k]:.2f}）　反り A−O = {sA[k]-sO[k]:.2f} µm（{trA[k]-trO[k]:.2f}）")
        return []
    anim = FuncAnimation(fig, draw, frames=len(t), interval=700)
    anim.save(os.path.join(IMG, "W_AO_times_temp_anim.gif"), writer=FixedPaletteWriter(fps=1.5), dpi=70)
    draw(10); fig.savefig(os.path.join(IMG, "W_AO_times_temp_t300.png"), dpi=85, facecolor="white"); plt.close(fig)
    print("wrote gif", os.path.getsize(os.path.join(IMG, "W_AO_times_temp_anim.gif")) // 1024, "KB")


if __name__ == "__main__":
    main()
