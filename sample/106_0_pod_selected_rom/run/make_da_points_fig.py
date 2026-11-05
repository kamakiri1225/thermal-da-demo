"""データ同化の時刻歴の図（da_timeseries.png、krylov_da_timeseries.png）で使った点の場所の図。
  温度計 P2・P0（観測）、変位計 S1（Ux）・S2（Uz）（観測、矢印＝測る向き）、
  測らない点 T3（上面・反対側、温度の時刻歴を描いた点）、評価点 A・O（反り A−O）、ヒータ。
  python3 run/make_da_points_fig.py → docs/img/da_points.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
import pyvista as pv
import vtk
from make_unobserved_points_fig import grid
from make_rom_field_points_fig import heater_surface

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); P = kv["xyz"] * 1000; Cc = kv["cell_centres"] * 1000
    T3 = Cc[int(np.linalg.norm(Cc - np.array([-28.0, 0.0, 97.5]), axis=1).argmin())]
    pts = {"P2": (P[2], "#C0392B", None), "P0": (P[0], "#C0392B", None),
           "S1": ((9.7, -36.2, 100.5), "#1565C0", (1, 0, 0)), "S2": ((-34.6, 14.4, 100.5), "#1565C0", (0, 0, 1)),
           "T3": (T3, "#7B2CBF", None), "A": ((28.7, 0, 100.5), "black", None), "O": ((-28.7, 0, 100.5), "black", None)}
    labels = {"P2": "P2 温度計", "P0": "P0 温度計", "S1": "S1 変位計（Ux）", "S2": "S2 変位計（Uz）",
              "T3": "T3 測らない点（温度の時刻歴）", "A": "A（反り A−O）", "O": "O（反り A−O）"}
    pv.OFF_SCREEN = True; g = grid()
    pl = pv.Plotter(off_screen=True, window_size=(900, 900)); pl.set_background("white")
    pl.add_mesh(g.extract_surface(), color="#c8d0da", opacity=0.28)
    pl.add_mesh(g.extract_feature_edges(), color="#555555", line_width=1.5)
    pl.add_mesh(heater_surface(), color="#F26B5B", opacity=0.45, lighting=False)
    for n, (x, c, d) in pts.items():
        pl.add_mesh(pv.Sphere(radius=3.0, center=x), color=c)
        if d: pl.add_mesh(pv.Arrow(start=x, direction=d, scale=18, tip_length=0.3, tip_radius=0.12, shaft_radius=0.05), color=c, lighting=False)
    pl.camera_position = [(230, -260, 230), (0, 0, 50), (0, 0, 1)]; pl.camera.zoom(1.25)
    img = pl.screenshot(return_img=True); xy = {}
    for n, (x, _, _) in pts.items():
        cc = vtk.vtkCoordinate(); cc.SetCoordinateSystemToWorld(); cc.SetValue(*map(float, x))
        u, v = cc.GetComputedDisplayValue(pl.renderer); xy[n] = (u, 900 - v)
    pl.close()
    nz = np.where((img[:, :, :3] < 245).any(axis=2)); r0, c0 = max(nz[0].min() - 40, 0), nz[1].min()
    img = img[r0:nz[0].max() + 30, c0:nz[1].max() + 1]; M = 300
    img = np.pad(img, ((0, 0), (M, M), (0, 0)), constant_values=255); W = img.shape[1]
    xy = {n: (u - c0 + M, v - r0) for n, (u, v) in xy.items()}
    plt.rcParams.update({"font.size": 18})
    fig = plt.figure(figsize=(11, 8.6)); ax = fig.add_axes([0, 0.08, 1, 0.92]); ax.imshow(img); ax.axis("off")
    for side in (-1, 1):
        items = sorted([(v, n, u) for n, (u, v) in xy.items() if (u < W / 2) == (side < 0)]); last = -1e9
        for v, n, u in items:
            tv = max(v, last + 55); last = tv; tu = M * 0.5 if side < 0 else W - M * 0.5
            ax.annotate(labels[n], (u, v), xytext=(tu, tv), fontsize=17, color=pts[n][1], ha="center", va="center",
                        weight="bold" if n == "T3" else "normal", arrowprops=dict(arrowstyle="-", color=pts[n][1], lw=1.3),
                        bbox=dict(fc="white", ec="none", pad=1))
    fig.text(0.5, 0.03, "赤＝温度計（観測）　青＝変位計（観測、矢印＝測る向き）　紫＝測らない点　黒＝反りの評価点（測らない）　薄い赤の面＝ヒータ",
             ha="center", fontsize=14)
    fig.savefig(os.path.join(IMG, "da_points.png"), dpi=100, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote da_points.png")


if __name__ == "__main__":
    main()
