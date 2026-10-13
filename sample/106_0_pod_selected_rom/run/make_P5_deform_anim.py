"""代表点5点の温度から求めた変形のアニメーション（GIF、1画面）.

  u(t) = Σ_j H_u[:, :, j] (T_Pj(t) − 20 °C) + c0     （H_u = D (Φ_P)^{-1}、c0 は 0.001 µm 以下の定数）
T_Pj(t) は OpenFOAM の真値（15 W）を代表5セルで取り出したもの。変形を拡大して描き、色は変位の大きさ |u|。

出力: docs/img/P5_deform_anim.gif
再現: OMP_NUM_THREADS=4 python3 run/make_P5_deform_anim.py
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


def main():
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
        p.add_point_labels(Pxyz, [f"P{j}" for j in range(5)], font_size=26, text_color="black", shape=None,
                           show_points=True, point_size=14, point_color="black", always_visible=True)
        p.add_point_labels([(R_OUT * 1000 + 18, 0, 50)], ["HEATER"], font_size=26, text_color="#c0392b", shape=None,
                           show_points=False, always_visible=True)
        p.camera_position = cam; p.camera.zoom(1.15)
        img = p.screenshot(return_img=True); p.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2))
        return np.ascontiguousarray(img[max(nz[0].min() - 10, 0):nz[0].max() + 10, max(nz[1].min() - 10, 0):nz[1].max() + 10, :3])

    frames = [render(k) for k in range(len(t))]
    plt.rcParams.update({"font.size": 16})
    fig = plt.figure(figsize=(9, 10.5))
    ax = fig.add_axes([0.02, 0.13, 0.96, 0.74]); ax.axis("off"); im = ax.imshow(frames[0])
    fig.text(0.5, 0.955, "代表点5点の温度から求めた変形", ha="center", fontsize=22, weight="bold")
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
    anim.save(os.path.join(IMG, "P5_deform_anim.gif"), writer=AllFramePaletteWriter(fps=1.5), dpi=85)
    draw(10); fig.savefig(os.path.join(IMG, "P5_deform_t300.png"), dpi=75, facecolor="white"); plt.close(fig)
    print("wrote gif", os.path.getsize(os.path.join(IMG, "P5_deform_anim.gif")) // 1024, "KB", "max |u|", mag.max())


if __name__ == "__main__":
    main()
