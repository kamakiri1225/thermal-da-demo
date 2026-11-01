"""未観測点（温度を比べた 4 か所）と、観測点・評価点の場所の図（斜めからと真上から）。

  run/holdout_unobserved_cells.py の 4 か所（results/holdout_unobserved_cells.json の cells）を描く。
  python3 run/make_unobserved_points_fig.py   → docs/img/unobserved_points.png
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
import cylinder_mesh
import pyvista as pv
import vtk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
TEMP = {"P2": (36.5, 6.9, 50.8), "P0": (21.3, -1.6, 51.6)}                      # 温度計
DISP = {"S1": (9.7, -36.2, 100.5), "S2": (-34.6, 14.4, 100.5)}                  # 変位計
EVAL = {"A": (28.0, 0.0, 100.5), "O": (-28.0, 0.0, 100.5)}                      # 反りの評価点
COL = dict(temp="#d62728", disp="#2ca02c", eval="#111111", unobs="#7b2cbf")


def grid():
    m = cylinder_mesh.build_cylinder_mesh(4, 48, 20, .020, .0375, .1005)
    xyz = np.array([p for _, p in m["nodes"]]) * 1000
    idx = {n: i for i, (n, _) in enumerate(m["nodes"])}
    cells = []
    for _, conn in m["elements"]:
        cells += [8] + [idx[n] for n in conn[:8]]
    return pv.UnstructuredGrid(np.array(cells), np.full(len(m["elements"]), 12, dtype=np.uint8), xyz)


def shot(g, unobs, cam, zoom, size=(900, 860)):
    p = pv.Plotter(off_screen=True, window_size=size); p.set_background("white")
    p.add_mesh(g.extract_surface(), color="#c8d0da", opacity=0.30, show_edges=False)
    p.add_mesh(g.extract_feature_edges(), color="#555555", line_width=1.5)
    for grp, pts, col in [("temp", TEMP, COL["temp"]), ("disp", DISP, COL["disp"]), ("eval", EVAL, COL["eval"]),
                          ("unobs", unobs, COL["unobs"])]:
        for name, x in pts.items():
            p.add_mesh(pv.Sphere(radius=3.2, center=x), color=col)
    p.camera_position = cam; p.camera.zoom(zoom)
    img = p.screenshot(return_img=True)
    xy = {}
    for pts in (TEMP, DISP, EVAL, unobs):
        for name, x in pts.items():
            c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*map(float, x))
            u, v = c.GetComputedDisplayValue(p.renderer); xy[name] = (u, size[1] - v)
    p.close(); return img, xy


def main():
    cells = json.load(open(os.path.join(RES, "holdout_unobserved_cells.json")))["cells"]
    keys = list(cells)
    unobs = {f"T{i + 1}": tuple(np.array(cells[k]) * 1000) for i, k in enumerate(keys)}
    pv.OFF_SCREEN = True
    g = grid()
    obl, xo = shot(g, unobs, [(230, -260, 230), (0, 0, 50), (0, 0, 1)], 1.25)
    top, xt = shot(g, unobs, [(0, 0, 400), (0, 0, 50), (0, 1, 0)], 1.55)
    colof = {**{k: COL["temp"] for k in TEMP}, **{k: COL["disp"] for k in DISP}, **{k: COL["eval"] for k in EVAL},
             **{k: COL["unobs"] for k in unobs}}
    # 名前の置き場所（点からのずれ [px]）。重ならないように手で決めた
    OFF_O = dict(T1=(-120, 10), T2=(90, 70), T3=(-150, 60), T4=(-130, -10), P0=(-30, -70), P2=(70, -40),
                 S1=(-130, 20), S2=(-80, -60), A=(90, -40), O=(-110, -50))
    OFF_T = dict(T1=(-60, -90), T2=(40, 90), T3=(-20, 120), T4=(-150, 40), P0=(-90, 100), P2=(60, -60),
                 S1=(70, 50), S2=(-80, -60), A=(80, 10), O=(-120, -10))

    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16})
    fig = plt.figure(figsize=(16, 7.4))
    a1 = fig.add_axes([0.0, 0.03, 0.33, 0.86]); a1.imshow(obl); a1.axis("off"); a1.set_title("斜めから見た図", fontsize=19)
    a2 = fig.add_axes([0.32, 0.03, 0.33, 0.86]); a2.imshow(top); a2.axis("off"); a2.set_title("真上から見た図（ヒータは右の +x 側）", fontsize=19)
    for ax, xy, off in [(a1, xo, OFF_O), (a2, xt, OFF_T)]:
        for n, (u, v) in xy.items():
            du, dv = off[n]
            ax.annotate(n, (u, v), xytext=(u + du, v + dv), color=colof[n], fontsize=17, weight="bold",
                        ha="center", va="center", arrowprops=dict(arrowstyle="-", color=colof[n], lw=1.5),
                        bbox=dict(fc="white", ec=colof[n], lw=1, pad=1.5))
    a3 = fig.add_axes([0.66, 0.0, 0.34, 1.0]); a3.axis("off")
    rows = [("■ 温度を比べた未観測点（観測しない）", COL["unobs"], True)]
    for i, k in enumerate(keys):
        x = np.array(cells[k]) * 1000
        rows.append((f"T{i + 1}：{k}\n　 ({x[0]:.1f}, {x[1]:.1f}, {x[2]:.1f}) mm", COL["unobs"], False))
    rows += [("■ 温度計：P2・P0（ヒータ側の中段）", COL["temp"], True),
             ("■ 変位計：S1（Ux）・S2（Uz）（上面）", COL["disp"], True),
             ("■ 反りの評価点：A・O（上面、観測しない）", COL["eval"], True)]
    y = 0.86
    for txt, col, bold in rows:
        a3.text(0.0, y, txt, color=col, fontsize=16, weight="bold" if bold else "normal", va="top", transform=a3.transAxes)
        y -= 0.11 if "\n" in txt else 0.075
    a3.text(0.0, 0.03, "座標は (x, y, z)、z が高さ。円筒の高さ 100.5 mm", fontsize=15, color="0.3", transform=a3.transAxes)
    fig.suptitle("温度を比べた未観測点 T1〜T4 と、観測点・評価点の場所", fontsize=20, y=0.99)
    fig.savefig(os.path.join(IMG, "unobserved_points.png"), dpi=100, bbox_inches="tight", facecolor="white")
    print("wrote unobserved_points.png")


if __name__ == "__main__":
    main()
