"""温度（色）と熱変形（形）のアニメーション（GIF）：OpenFOAM＋FrontISTR の計算結果そのもの（15 W、0〜600 秒）。

  温度：OpenFOAM（CHT）の固体の温度（results/limit_truth_learned.npz の Tfield、30 秒ごと）を、
        FrontISTR の節点へ IDW（近傍 8 セル）で渡したもの
  変形：その温度を FrontISTR に入れて計算した全節点の変位（results/truth_disp_all_learned.npz、holdout_AO_realsolver.py が作成）
  形は変位を 2000 倍にして描く。ROM・データ同化は使っていない。
  引数 rom を付けると、正解の代わりに「ROM を解いて gappy-POD で戻した温度と変位」を描く：
    ROM（代表 5 点の熱回路、results/rom_calibrated_pod.npz）を 15 W・0〜300 秒加熱で 0→600 秒解く（データ同化なし、
    初期 20 °C）→ 代表 5 点の温度から gappy-POD で全 20,696 セルの温度（T = 平均 + Φ Φ_P⁻¹ (T_P − 平均_P)）
    → 変位は、POD の型ごとに FrontISTR で前もって解いた変位の足し合わせ（u = u_平均 + D Φ_P⁻¹ (T_P − 平均_P)）。
  出力：docs/img/temp_deform_anim.gif、temp_deform_t300.png（正解）
        docs/img/temp_deform_anim_rom.gif、temp_deform_t300_rom.png（ROM ＋ gappy-POD）
  再現：python3 run/make_temp_deform_anim.py [rom]
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p  # noqa: F401
from dacore.temperature_mapping import build_idw_mapping, map_temperature
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import pyvista as pv
import cylinder_mesh
from make_P5_deform_anim import AllFramePaletteWriter
from dacore import rom_general as rg

pv.OFF_SCREEN = True
RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NR, NTH, NZ = 4, 48, 20; R_IN, R_OUT, H_ = 0.020, 0.0375, 0.1005
SCALE = 2000.0; CMAP = "turbo"


def main():
    mesh = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, R_IN, R_OUT, H_)
    coords = np.array([x for _, x in mesh["nodes"]]); pos = {n: i for i, (n, _) in enumerate(mesh["nodes"])}
    cells = np.array([[8] + [pos[n] for n in conn] for _, conn in mesh["elements"]]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(mesh["elements"]), pv.CellType.HEXAHEDRON), coords * 1000)
    surf = grid.extract_surface(algorithm="dataset_surface"); oid = surf["vtkOriginalPointIds"]
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; Tf = z["Tfield"]
    M = build_idw_mapping(kv["cell_centres"], coords)
    Tn_true = np.array([map_temperature(M, Tf[k]) for k in range(len(t))]) - 273.15  # 正解の節点温度 [°C]
    u_true = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]           # 正解の変位 (時刻, 節点, 3) µm
    assert u_true.shape[0] == len(t) and u_true.shape[1] == len(coords)
    ROM = "rom" in sys.argv[1:]
    if ROM:
        d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
        C = d["C"]; K = rg.tri_to_matrix(d["K_upper"], 5); h = float(d["h"]); heat = int(d["heat_node"])
        TP = [np.full(5, rg.T_AIR_K)]
        for a, b in zip(t[:-1], t[1:]):                                              # ROM を 30 秒ずつ進める（中は 2 秒刻み）
            TP.append(rg.integrate_ensemble(TP[-1][None], C, K, np.array([h]), np.array([1.0]), heat, a, b, 2.0)[0])
        TP = np.array(TP)                                                            # (時刻, 5) K
        U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]
        a_mode = (TP - mean[pod]) @ np.linalg.inv(U[pod]).T                          # 型の重み（5 点から）
        Tcell = mean[None] + a_mode @ U.T                                            # gappy-POD：全セルの温度
        Tn = np.array([map_temperature(M, Tcell[k]) for k in range(len(t))]) - 273.15
        op = np.load(os.path.join(RES, "dispop_allnodes.npz"))
        u = op["u_mean"][None] + np.einsum("nck,tk->tnc", op["D"], a_mode)            # 型ごとの変位の足し合わせ
        eT = np.sqrt(((Tn - Tn_true) ** 2).mean(1)); eU = np.abs(np.linalg.norm(u, axis=2).max(1) - np.linalg.norm(u_true, axis=2).max(1))
        print("ROM vs 正解：節点温度 RMSE 最大", eT.max().round(3), "K、変位の最大の差 最大", eU.max().round(3), "µm")
    else:
        Tn, u = Tn_true, u_true
    tmin, tmax = 20.0, float(np.ceil(Tn_true.max()))                              # 正解と同じ色の範囲
    cam = [(190, -230, 210), (0, 0, 50), (0, 0, 1)]

    def render(k):
        p = pv.Plotter(off_screen=True, window_size=(900, 900))
        s = surf.copy(); s.points = surf.points + u[k][oid] * 1e-3 * SCALE; s.point_data["T"] = Tn[k][oid]
        p.add_mesh(surf, style="wireframe", color="#AAAAAA", opacity=0.25, line_width=1)
        p.add_mesh(s, scalars="T", cmap=CMAP, clim=[tmin, tmax], show_scalar_bar=False, lighting=False)
        p.add_mesh(s.extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False, feature_angle=30),
                   color="#333333", line_width=2)
        p.add_point_labels([(R_OUT * 1000 + 18, 0, 50)], ["HEATER"], font_size=26, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.15)
        img = p.screenshot(return_img=True); p.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2))
        return np.ascontiguousarray(img[max(nz[0].min() - 10, 0):nz[0].max() + 10, max(nz[1].min() - 10, 0):nz[1].max() + 10, :3])

    frames = [render(k) for k in range(len(t))]
    umax = np.linalg.norm(u, axis=2).max(axis=1)
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(9, 10.5))
    ax = fig.add_axes([0.02, 0.15, 0.96, 0.69]); ax.axis("off"); im = ax.imshow(frames[0])
    if ROM:
        fig.text(0.5, 0.955, "低次元モデル（ROM）から戻した温度と熱変形", ha="center", fontsize=22, weight="bold")
        fig.text(0.5, 0.915, f"代表 5 点の ROM を解く → gappy-POD で全体の温度・変位。15 W。変形は {SCALE:.0f} 倍",
                 ha="center", fontsize=14, color="#444")
    else:
        fig.text(0.5, 0.955, "温度（色）と熱変形（形）", ha="center", fontsize=22, weight="bold")
        fig.text(0.5, 0.915, f"OpenFOAM（温度）→ FrontISTR（熱変形）。15 W。変形は {SCALE:.0f} 倍（灰色の線＝変形前）",
                 ha="center", fontsize=14, color="#444")
    tx = fig.text(0.5, 0.875, "", ha="center", fontsize=17)
    cax = fig.add_axes([0.15, 0.07, 0.70, 0.025])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(tmin, tmax), cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_label("温度 [°C]", fontsize=15)

    def draw(k):
        im.set_data(frames[k])
        tx.set_text(f"t = {t[k]:.0f} 秒" + ("（加熱中）" if t[k] <= 300 else "（加熱停止後）") +
                    f"　最高 {Tn[k].max():.1f} °C　変位の最大 {umax[k]:.1f} µm")
        return []
    anim = FuncAnimation(fig, draw, frames=len(t), interval=700)
    suf = "_rom" if ROM else ""
    anim.save(os.path.join(IMG, f"temp_deform_anim{suf}.gif"), writer=AllFramePaletteWriter(fps=1.5), dpi=85)
    draw(int(np.argmin(abs(t - 300)))); fig.savefig(os.path.join(IMG, f"temp_deform_t300{suf}.png"), dpi=75, facecolor="white")
    plt.close(fig)
    print(f"wrote temp_deform_anim{suf}.gif", "Tmax", Tn.max(), "umax", umax.max())


if __name__ == "__main__":
    main()
