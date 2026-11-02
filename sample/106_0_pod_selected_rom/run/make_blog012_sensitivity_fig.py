"""blog_012 §2 の図：型 1・型 2 の形と、温度計 P2・P0 の場所での値（＝温度計の感度 h）。

  型は「一番大きく変わるセルで 1」にそろえた（blog_012 と同じ）。色は型の値、球は温度計 P2・P0。
  python3 run/make_blog012_sensitivity_fig.py → docs/img/blog012_sensitivity.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p  # noqa: F401
from dacore.temperature_mapping import build_idw_mapping, map_temperature
import matplotlib.pyplot as plt
import pyvista as pv
import vtk
from make_unobserved_points_fig import grid

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    Us = U / np.abs(U).max(0)
    g = grid(); M = build_idw_mapping(Cc, g.points / 1000)
    pv.OFF_SCREEN = True
    P = {"P2": kv["xyz"][2] * 1000, "P0": kv["xyz"][0] * 1000}
    shots = []
    for k in range(2):
        gg = g.copy(); gg.point_data["v"] = map_temperature(M, Us[:, k])
        pl = pv.Plotter(off_screen=True, window_size=(800, 820)); pl.set_background("white")
        pl.add_mesh(gg, scalars="v", cmap="coolwarm", clim=[-1, 1], show_scalar_bar=False, opacity=0.55, lighting=False)
        pl.add_mesh(gg.extract_feature_edges(), color="#333333", line_width=1.2)
        for x in P.values(): pl.add_mesh(pv.Sphere(radius=3.0, center=x), color="black")
        pl.camera_position = [(260, -230, 200), (0, 0, 50), (0, 0, 1)]; pl.camera.zoom(1.3)
        img = pl.screenshot(return_img=True)
        xy = {}
        for n, x in P.items():
            c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*map(float, x))
            u, v = c.GetComputedDisplayValue(pl.renderer); xy[n] = (u, 820 - v)
        pl.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2)); r0, c0 = nz[0].min() - 10, nz[1].min() - 150
        img = img[max(r0, 0):nz[0].max() + 10, max(c0, 0):nz[1].max() + 150]
        xy = {n: (u - max(c0, 0), v - max(r0, 0)) for n, (u, v) in xy.items()}
        shots.append((img, xy))
    plt.rcParams.update({"font.size": 16})
    fig, axs = plt.subplots(1, 2, figsize=(14, 7.6), gridspec_kw=dict(wspace=0.05))
    titles = ["型 1（全体が温まる）", "型 2（ヒータ側と反対側の差）"]
    off = {"P2": (150, 40), "P0": (-160, 70)}
    for k, (ax, (img, xy)) in enumerate(zip(axs, shots)):
        ax.imshow(img); ax.axis("off"); ax.set_title(titles[k], fontsize=19)
        for n, (u, v) in xy.items():
            j = 2 if n == "P2" else 0
            ax.annotate(f"{n}：{Us[pod[j], k]:+.2f}", (u, v), xytext=(u + off[n][0], v + off[n][1]), fontsize=18, weight="bold",
                        ha="center", va="center", arrowprops=dict(arrowstyle="-", color="black", lw=1.4),
                        bbox=dict(fc="white", ec="black", lw=1, pad=2))
        ax.text(0.98, 0.02, "ヒータは右手前（+x 側）", transform=ax.transAxes, ha="right", fontsize=15)
    cax = fig.add_axes([0.25, 0.04, 0.5, 0.025])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap="coolwarm"), cax=cax, orientation="horizontal")
    cb.set_label("型の値（一番大きく変わるセルで ±1）", fontsize=16)
    fig.savefig(os.path.join(IMG, "blog012_sensitivity.png"), dpi=100, bbox_inches="tight", facecolor="white")
    print("wrote blog012_sensitivity.png", {n: Us[pod[2 if n == 'P2' else 0], :2].round(3).tolist() for n in P})


if __name__ == "__main__":
    main()
