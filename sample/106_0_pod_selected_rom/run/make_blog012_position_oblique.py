"""blog_012 §8-2 の「選ばれた位置」を、斜め上から見た 3D の図にする（6 方式、1 方式 1 枚）。

  位置と向きは results/placement_method_comparison.json（本番の選定結果）から読む。
  半透明の円筒、ヒータ（薄い赤）、温度計 P2・P0（赤の球）、評価点 A・O（黒の球）、
  選んだ変位計 S1・S2（青の球＋測る向きの矢印）を描く。
  出力：docs/img/blog012_positions_{goal,magnitude,qdeim,aopt,dopt,direct}.png（上から見た旧版を置き換える）
  再現：python3 run/make_blog012_position_oblique.py
"""
import os, sys, json
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
METHODS = [("目的指向（本研究）", "goal", "目的指向"), ("よく動く点に置く", "magnitude", "感度の大きさ"),
           ("Q-DEIM（QR列ピボット）", "qdeim", "Q-DEIM"), ("A最適（温度の trace 最小）", "aopt", "A最適"),
           ("D最適（温度の det 最小）", "dopt", "D最適"), ("A・O 自身（循環・参考）", "direct", "A・O直接（参考）")]
DIRV = {"Ux": (1, 0, 0), "Uy": (0, 1, 0), "Uz": (0, 0, 1)}
SC = "#1565C0"; CAM = [(230, -260, 230), (0, 0, 50), (0, 0, 1)]


def main():
    j = json.load(open(os.path.join(RES, "placement_method_comparison.json")))["methods"]
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); Pt = kv["xyz"] * 1000
    g = grid(); pv.OFF_SCREEN = True
    fixed = {"A": ((28.7, 0, 100.5), "black"), "O": ((-28.7, 0, 100.5), "black"),
             "P2": (tuple(Pt[2]), "#C0392B"), "P0": (tuple(Pt[0]), "#C0392B")}
    for key, fn, label in METHODS:
        m = j[key]; S = {f"S{i + 1}": (tuple(x), c) for i, (x, c) in enumerate(zip(m["xyz_mm"], m["comp"]))}
        pl = pv.Plotter(off_screen=True, window_size=(900, 900)); pl.set_background("white")
        pl.add_mesh(g.extract_surface(), color="#c8d0da", opacity=0.28)
        pl.add_mesh(g.extract_feature_edges(), color="#555555", line_width=1.5)
        pl.add_mesh(heater_surface(), color="#F26B5B", opacity=0.45, lighting=False)
        for _, (x, c) in fixed.items(): pl.add_mesh(pv.Sphere(radius=2.6, center=x), color=c)
        for _, (x, comp) in S.items():
            pl.add_mesh(pv.Sphere(radius=3.0, center=x), color=SC)
            pl.add_mesh(pv.Arrow(start=x, direction=DIRV[comp], scale=18, tip_length=0.3, tip_radius=0.12, shaft_radius=0.05),
                        color=SC, lighting=False)
        pl.camera_position = CAM; pl.camera.zoom(1.25)
        img = pl.screenshot(return_img=True); xy = {}
        for n, (x, _) in {**fixed, **S}.items():
            cc = vtk.vtkCoordinate(); cc.SetCoordinateSystemToWorld(); cc.SetValue(*map(float, x))
            u, v = cc.GetComputedDisplayValue(pl.renderer); xy[n] = (u, 900 - v)
        pl.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2))
        img = img[max(nz[0].min() - 50, 0):nz[0].max() + 30, nz[1].min():nz[1].max() + 1]
        r0, c0 = max(nz[0].min() - 50, 0), nz[1].min()
        M = 260                                                    # 名前を書く左右の余白 [px]
        img = np.pad(img, ((0, 0), (M, M), (0, 0)), constant_values=255)
        xy = {n: (u - c0 + M, v - r0) for n, (u, v) in xy.items()}
        W = img.shape[1]; cx = W / 2
        plt.rcParams.update({"font.size": 18})
        fig = plt.figure(figsize=(9, 8.6)); ax = fig.add_axes([0, 0.1, 1, 0.84]); ax.imshow(img); ax.axis("off")
        fig.text(0.5, 0.975, label, ha="center", va="top", fontsize=24, weight="bold")
        lab = {**{n: n for n in fixed}, **{n: f"{n}（{c}）" for n, (_, c) in S.items()}}
        col = {"A": "black", "O": "black", "P2": "#C0392B", "P0": "#C0392B", "S1": SC, "S2": SC}
        # 名前は、点が中心より左なら左の余白、右なら右の余白に置く。同じ側で近いときは縦にずらす
        for side in (-1, 1):
            items = sorted([(v, n, u) for n, (u, v) in xy.items() if (u < cx) == (side < 0)])
            last = -1e9
            for v, n, u in items:
                tv = max(v, last + 55); last = tv
                tu = M * 0.55 if side < 0 else W - M * 0.55
                ax.annotate(lab[n], (u, v), xytext=(tu, tv), fontsize=19 if n.startswith("S") else 18,
                            weight="bold" if n.startswith("S") else "normal", color=col[n], ha="center", va="center",
                            arrowprops=dict(arrowstyle="-", color=col[n], lw=1.3), bbox=dict(fc="white", ec="none", pad=1))
        fig.text(0.5, 0.06, "青＝選んだ変位計（矢印＝測る向き）　赤＝温度計 P2・P0\n黒＝評価点 A・O（測らない）　薄い赤の面＝ヒータ",
                 ha="center", va="center", fontsize=16)
        fig.savefig(os.path.join(IMG, f"blog012_positions_{fn}.png"), dpi=100, bbox_inches="tight", facecolor="white")
        plt.close(fig); print("wrote", fn)


if __name__ == "__main__":
    main()
