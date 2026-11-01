"""図 8c（run/rom_field_vs_openfoam.py）で代表点以外の温度を比べた 3 か所（＋代表点をつなぐ熱回路の K） T2・T3・T4 と、代表点 P0〜P4、ヒータの場所の図。

  色は図 8c の線と同じ（T2 赤・T3 青・T4 緑）。代表点は灰色。ヒータ（外周 +x 側、高さ 25〜75 mm、周長 100 mm）は黄色。
  python3 run/make_rom_field_points_fig.py   → docs/img/rom_field_points.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
import pyvista as pv
import vtk
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from make_unobserved_points_fig import grid
from rom_field_vs_openfoam import PLACES

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
GRAY = "#555555"; HEAT = "#E6B800"


def heater_surface():
    """外周 r = 37.5 mm の少し外側に、ヒータの面（周長 100 mm、高さ 25〜75 mm）を作る。"""
    r = 38.6; half = 100.0 / 37.5 / 2
    th, z = np.meshgrid(np.linspace(-half, half, 40), np.linspace(25.25, 75.25, 12))
    return pv.StructuredGrid(r * np.cos(th), r * np.sin(th), z)


def shot(g, pts, cam, zoom, size=(900, 860), links=()):
    p = pv.Plotter(off_screen=True, window_size=size); p.set_background("white")
    p.add_mesh(g.extract_surface(), color="#c8d0da", opacity=0.30, show_edges=False)
    p.add_mesh(g.extract_feature_edges(), color="#555555", line_width=1.5)
    p.add_mesh(heater_surface(), color="#F26B5B", opacity=0.55, lighting=False)   # ヒータは赤っぽく（陰影なし）
    for a, b, w in links:   # 代表点どうしをつなぐ熱コンダクタンス K（熱回路の「ばね」）。太さは K に合わせる
        p.add_mesh(pv.Tube(pointa=a, pointb=b, radius=w, n_sides=16), color="#444444", lighting=False)
    for name, (x, col, r) in pts.items():
        p.add_mesh(pv.Sphere(radius=r, center=x), color=col)
    p.camera_position = cam; p.camera.zoom(zoom)
    img = p.screenshot(return_img=True)
    xy = {}
    for name, (x, _, _) in pts.items():
        c = vtk.vtkCoordinate(); c.SetCoordinateSystemToWorld(); c.SetValue(*map(float, x))
        u, v = c.GetComputedDisplayValue(p.renderer); xy[name] = (u, size[1] - v)
    p.close(); return img, xy


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); Cc = kv["cell_centres"] * 1000; P = kv["xyz"] * 1000
    pts = {}
    for lab, x, col in PLACES:
        name = lab.split()[0]
        c = Cc[int(np.linalg.norm(Cc - np.array(x) * 1000, axis=1).argmin())]   # 図 8c と同じセル
        pts[name] = (c, col, 3.4)
    for j in range(5):
        pts[f"P{j}"] = (P[j], GRAY, 2.4)
    pv.OFF_SCREEN = True
    g = grid()
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz")); from dacore import rom_general as rg
    K = rg.tri_to_matrix(d["K_upper"], 5)
    links = [(P[i], P[j], 0.35 + 0.55 * np.sqrt(K[i, j] / K.max())) for i in range(5) for j in range(i + 1, 5) if K[i, j] > 1e-3]
    obl, xo = shot(g, pts, [(230, -260, 230), (0, 0, 50), (0, 0, 1)], 1.25, links=links)
    # 白い余白を切る（名前を置く分だけ残す）。点の画面上の位置も同じだけずらす
    nz = np.where((obl[:, :, :3] < 245).any(axis=2)); m = 62
    r0, c0 = max(nz[0].min() - 8, 0), max(nz[1].min() - m, 0)
    obl = obl[r0:nz[0].max() + 28, c0:nz[1].max() + m]
    xo = {k: (u - c0, v - r0) for k, (u, v) in xo.items()}
    # 名前の置き場所（点からのずれ [px]）。重ならないように手で決めた
    OFF_O = dict(T2=(70, -30), T3=(-75, -25), T4=(-80, 22), P0=(-60, -70), P1=(-110, 40), P2=(55, 30), P3=(55, 25), P4=(-60, 35))
    col = {k: v[1] for k, v in pts.items()}

    # 文字はすべて黒の太字（色付きの字・枠は使わない）。色は点と、説明の印だけに付ける
    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 16, "text.color": "black", "font.weight": "bold"})   # 字は濃い黒（太字）
    # スライドの狭い場所に置けるよう、凡例は図の上にまとめ、図を縦長にする（座標はブログの図の説明に書く）
    ih, iw = obl.shape[:2]; W = 7.6; imgH = W * ih / iw; legH = 1.1; H = imgH + legH
    fig = plt.figure(figsize=(W, H))
    a1 = fig.add_axes([0.0, 0.0, 1.0, imgH / H]); a1.imshow(obl); a1.axis("off")
    for ax, xy, off in [(a1, xo, OFF_O)]:
        for n, (u, v) in xy.items():
            du, dv = off[n]
            if n.startswith("T"):   # 代表点以外：大きい太字、黒枠で囲む
                ax.annotate(n, (u, v), xytext=(u + du, v + dv), color="black", fontsize=19, weight="bold",
                            ha="center", va="center", arrowprops=dict(arrowstyle="-", color="black", lw=1.6),
                            bbox=dict(fc="white", ec="black", lw=1.2, pad=2.0))
            else:                   # 代表点：小さめの普通の太さ、枠なし
                ax.annotate(n, (u, v), xytext=(u + du, v + dv), color="#333333", fontsize=15, weight="normal",
                            ha="center", va="center", arrowprops=dict(arrowstyle="-", color="#555555", lw=1.0),
                            bbox=dict(fc="white", ec="none", pad=0.5))
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    hs = [Line2D([], [], ls="", marker="o", ms=13, color=c, label=lab) for lab, _, c in PLACES]
    hs += [Line2D([], [], ls="", marker="o", ms=10, color=GRAY, label="代表点 P0〜P4"),
           Line2D([], [], color="#444444", lw=4, label="熱のばね K（太いほど大）"),
           Patch(fc="#F26B5B", ec="black", lw=0.6, label="ヒータ（+x 側）")]
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=16,
               columnspacing=1.2, handletextpad=0.5, labelspacing=0.45, prop={"weight": "bold", "size": 16})
    fig.savefig(os.path.join(IMG, "rom_field_points.png"), dpi=100, bbox_inches="tight", pad_inches=0.04, facecolor="white")
    print("wrote rom_field_points.png")


if __name__ == "__main__":
    main()
