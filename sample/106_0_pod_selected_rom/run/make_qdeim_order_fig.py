"""Q-DEIMが代表点を1点ずつ選んでいく様子をParaView風3D図で示す.

5コマ：1点目 → 2点まで → … → 5点まで。新しく追加された点を強調する。
出力: docs/img/qdeim_order.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_qdeim_order_fig.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
IMG = os.path.join(ROOT, "docs", "img"); RES = os.path.join(ROOT, "results")
# 貪欲選択の順序と残差（run/qdeim_explained.py の出力）
RESID = [0.0369, 0.0266, 0.0192, 0.0173, 0.0142]


def trim(img, pad=6):
    m = np.any(img[..., :3] < 246, axis=-1)
    ys, xs = np.where(m)
    return img[max(0, ys.min()-pad):ys.max()+pad, max(0, xs.min()-pad):xs.max()+pad]


def main():
    d = np.load(os.path.join(RES, "qdeim_points.npz"))
    xyz = d["xyz"]           # 選ばれた順に並んでいる
    import pyvista as pv, vtk
    pv.OFF_SCREEN = True
    try: pv.start_xvfb()
    except Exception: pass
    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    co = np.array([c for _n, c in mesh["nodes"]])
    idr = {nid: i for i, (nid, _x) in enumerate(mesh["nodes"])}
    cl = []
    for _e, conn in mesh["elements"]: cl.append(8); cl.extend(idr[n] for n in conn)
    ug = pv.UnstructuredGrid(np.array(cl),
                             np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, np.uint8), co)

    shots = []
    for k in range(5):
        pl = pv.Plotter(off_screen=True, window_size=(620, 740))
        pl.add_mesh(ug, color="lightsteelblue", opacity=0.30, show_edges=False)
        for j in range(k + 1):
            newest = (j == k)
            pl.add_mesh(pv.Sphere(radius=0.0040 if newest else 0.0026, center=xyz[j]),
                        color="#c0392b" if newest else "#7f8c9a")
            pl.add_point_labels([xyz[j] + np.array([0.005, 0, 0.007])],
                                [f"{j+1}"], font_size=26 if newest else 18,
                                text_color="#c0392b" if newest else "#5a6672",
                                shape=None, always_visible=True)
        pl.camera_position = [(0.26, -0.24, 0.20), (0, 0, 0.05), (0, 0, 1)]
        pl.set_background("white"); pl.camera.zoom(1.25)
        f = os.path.join(tempfile.mkdtemp(), f"s{k}.png")
        pl.screenshot(f); pl.close()
        shots.append(trim(plt.imread(f)))
        print(f"  {k+1}点目: cell座標 {np.round(xyz[k]*1000,1)} mm")

    fig, axes = plt.subplots(1, 5, figsize=(17.0, 5.2))
    for k, (ax, im) in enumerate(zip(axes, shots)):
        ax.imshow(im); ax.axis("off")
        ax.set_title(f"{k+1}点目を追加\n"
                     f"残差 {RESID[k]:.4f}" + ("  → 0.0000" if k == 4 else ""),
                     fontsize=12.5, weight="bold",
                     color="#c0392b" if k == 4 else "#1b2430")
    fig.suptitle("Q-DEIMが代表点を選んでいく順序：赤が新しく追加された点、灰色がすでに選んだ点\n"
                 "1点目はヒータ最近傍ではなく中央高さ。2点目は反対側、最後に底面の両側へ ― 場を広く覆うように広がっていく",
                 fontsize=13.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.86])
    out = os.path.join(IMG, "qdeim_order.png")
    fig.savefig(out, dpi=130, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
