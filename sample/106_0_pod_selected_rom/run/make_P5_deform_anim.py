"""代表点5点の温度から求めた変形のアニメーション（GIF、1画面）.

  u(t) = Σ_j H_u[:, :, j] (T_Pj(t) − 20 °C) + c0     （H_u = D (Φ_P)^{-1}、c0 は 0.001 µm 以下の定数）
T_Pj(t) は OpenFOAM の真値（15 W）を代表5セルで取り出したもの。変形を拡大して描き、色は変位の大きさ |u|。

引数 b009 / b010 を付けると、その記事の観測点（温度 T1・T2、変位 S1・S2 と測る向き）と評価点 A・O も描く。
出力: docs/img/P5_deform_anim.gif（引数なし）、P5_deform_anim_b009.gif、P5_deform_anim_b010.gif
再現: OMP_NUM_THREADS=4 python3 run/make_P5_deform_anim.py [b009|b010]
"""
from __future__ import annotations
import os, sys
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
from matplotlib.animation import PillowWriter


class AllFramePaletteWriter(PillowWriter):
    """GIF の 256 色の表を、全コマを並べた画像から作る（最初のコマだけから作ると後のコマの色が化ける）。"""
    def finish(self):
        from PIL import Image
        rgb = [f.convert("RGB") for f in self._frames]
        w, h = rgb[0].size
        mosaic = Image.new("RGB", (w, h * len(rgb)))
        for i, f in enumerate(rgb): mosaic.paste(f, (0, h * i))
        pal = mosaic.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
        frames = [f.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG) for f in rgb]
        frames[0].save(self.outfile, save_all=True, append_images=frames[1:], duration=int(1000 / self.fps), loop=0)
pv.OFF_SCREEN = True

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NR, NTH, NZ = 4, 48, 20; R_IN, R_OUT, H_ = 0.020, 0.0375, 0.1005
TREF = 293.15; SCALE = 2000.0; CMAP = "turbo"
# 記事ごとの観測点 [mm]。変位計は測る向き（0=x, 1=y, 2=z）も持つ。
SETS = {
    "b009": {"title": "blog009 の観測点", "A": (28.75, 0, 100.5), "O": (-28.75, 0, 100.5),
             "T": {"T1": (36.7, 0.8, 49.2), "T2": (-35.9, -1.5, 51.6)},
             "S": {"S1 (Ux)": ((4.9, -37.2, 100.5), 0), "S2 (Ux)": ((4.9, 37.2, 100.5), 0)}},
    "b010": {"title": "blog010 の観測点", "A": (0, 28.75, 100.5), "O": (0, -28.75, 100.5),
             "T": {"T1": (-36.6, 0.8, 99.3), "T2": (21.1, 0.8, 0.8)},
             "S": {"S1 (Uy)": ((-17.3, -10.0, 100.5), 1), "S2 (Uy)": ((21.1, 12.2, 100.5), 1)}},
}


def main():
    key = next((a for a in sys.argv[1:] if a in SETS), None); cfg = SETS.get(key)
    mesh = cylinder_mesh.build_cylinder_mesh(NR, NTH, NZ, R_IN, R_OUT, H_)
    ids = [n for n, _ in mesh["nodes"]]; coords = np.array([x for _, x in mesh["nodes"]])
    pos = {n: i for i, n in enumerate(ids)}
    cells = np.array([[8] + [pos[n] for n in conn] for _, conn in mesh["elements"]]).ravel()
    grid = pv.UnstructuredGrid(cells, np.full(len(mesh["elements"]), pv.CellType.HEXAHEDRON), coords * 1000)
    surf = grid.extract_surface(algorithm="dataset_surface"); oid = surf["vtkOriginalPointIds"]
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); assert np.abs(op["coords"] - coords).max() < 1e-12
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; mean = kv["mean"]
    Hu = np.einsum("nck,kj->ncj", op["D"], np.linalg.inv(U[pod]))
    c0 = op["u_mean"] - np.einsum("ncj,j->nc", Hu, mean[pod] - TREF)
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; TP = z["Tfield"][:, pod]
    u = np.einsum("ncj,tj->tnc", Hu, TP - TREF) + c0[None]                     # (21,5040,3) µm
    mag = np.linalg.norm(u, axis=2); vmax = float(np.ceil(mag.max()))
    Pxyz = kv["xyz"] * 1000
    cam = [(190, -230, 210), (0, 0, 50), (0, 0, 1)]

    def render(k):
        p = pv.Plotter(off_screen=True, window_size=(900, 900))
        s = surf.copy(); s.points = surf.points + u[k][oid] * 1e-3 * SCALE; s.point_data["m"] = mag[k][oid]
        p.add_mesh(surf, style="wireframe", color="#AAAAAA", opacity=0.25, line_width=1)
        p.add_mesh(s, scalars="m", cmap=CMAP, clim=[0, vmax], show_scalar_bar=False, lighting=False)   # 陰影なし＝カラーバーと同じ色
        p.add_mesh(s.extract_feature_edges(boundary_edges=True, feature_edges=True, manifold_edges=False, feature_angle=30), color="#333333", line_width=2)
        if cfg is None:
            p.add_point_labels(Pxyz, [f"P{j}" for j in range(5)], font_size=26, text_color="black", shape=None,
                               show_points=True, point_size=14, point_color="black", always_visible=True)
        else:
            def moved(x):   # 一番近い節点の変位だけ動かした位置（変形後の形の上に描く）
                i = int(np.linalg.norm(coords * 1000 - np.array(x), axis=1).argmin())
                return np.array(x) + u[k][i] * 1e-3 * SCALE
            # 印は球（形状として描くので、ラベルの重なりで消えない）。文字は別に付ける。
            def mark(x, col, r, name, fs, dz=5.0):
                p.add_mesh(pv.Sphere(radius=r, center=x), color=col, lighting=False)
                p.add_point_labels([np.asarray(x) + [0, 0, dz]], [name], font_size=fs, text_color=col, shape="rounded_rect",
                                   shape_color="white", shape_opacity=0.75, show_points=False, always_visible=True)
            for j in range(5): mark(moved(Pxyz[j]), "#555555", 1.2, f"P{j}", 16, 3.5)
            mark(moved(cfg["A"]), "black", 2.4, "A", 30)
            mark(moved(cfg["O"]), "black", 2.4, "O", 30)
            for name, x in cfg["T"].items(): mark(moved(x), "#C2185B", 2.6, name, 28)
            for name, (x, comp) in cfg["S"].items():
                xm = moved(x); d = np.zeros(3); d[comp] = 1.0
                p.add_mesh(pv.Arrow(start=xm - 7 * d, direction=d, scale=14.0, tip_length=0.35, tip_radius=0.12, shaft_radius=0.05),
                           color="#1565C0", lighting=False)
                mark(xm, "#1565C0", 2.2, name, 26, 7.0)
        p.add_point_labels([(R_OUT * 1000 + 18, 0, 50)], ["HEATER"], font_size=26, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.15)
        img = p.screenshot(return_img=True); p.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2))
        return np.ascontiguousarray(img[max(nz[0].min() - 10, 0):nz[0].max() + 10, max(nz[1].min() - 10, 0):nz[1].max() + 10, :3])

    frames = [render(k) for k in range(len(t))]
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(9, 10.5))
    ax = fig.add_axes([0.02, 0.15, 0.96, 0.69]); ax.axis("off"); im = ax.imshow(frames[0])
    fig.text(0.5, 0.955, "代表点5点の温度から求めた変形" + (f"（{cfg['title']}）" if cfg else ""), ha="center", fontsize=22, weight="bold")
    if cfg:
        fig.text(0.5, 0.135, "赤＝温度センサ T1・T2　青＝変位センサ S1・S2（矢印＝測る向き）\n黒＝評価点 A・O（測らない）　灰＝代表点 P0〜P4　（裏側・内側の点は文字だけ見える）",
                 ha="center", va="top", fontsize=13, color="#333")
    fig.text(0.5, 0.915, f"15 W。変形は {SCALE:.0f} 倍に拡大（灰色の線＝変形前）", ha="center", fontsize=14, color="#444")
    tx = fig.text(0.5, 0.875, "", ha="center", fontsize=17)
    cax = fig.add_axes([0.15, 0.07, 0.70, 0.025])
    cb = fig.colorbar(ScalarMappable(norm=Normalize(0, vmax), cmap=CMAP), cax=cax, orientation="horizontal")
    cb.set_label("変位の大きさ [µm]", fontsize=15)

    def draw(k):
        im.set_data(frames[k])
        tx.set_text(f"t = {t[k]:.0f} 秒" + ("（加熱中）" if t[k] <= 300 else "（加熱停止後）") +
                    "　　" + "  ".join(f"P{j} {TP[k, j] - 273.15:.1f}℃" for j in range(5)))
        return []
    anim = FuncAnimation(fig, draw, frames=len(t), interval=700)
    suf = f"_{key}" if cfg else ""
    anim.save(os.path.join(IMG, f"P5_deform_anim{suf}.gif"), writer=AllFramePaletteWriter(fps=1.5), dpi=85)
    draw(10); fig.savefig(os.path.join(IMG, f"P5_deform_t300{suf}.png"), dpi=75, facecolor="white"); plt.close(fig)
    print("wrote gif", suf, os.path.getsize(os.path.join(IMG, f"P5_deform_anim{suf}.gif")) // 1024, "KB", "max |u|", mag.max())


if __name__ == "__main__":
    main()
