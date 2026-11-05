"""クリロフ基底の最初の 4 本を円筒の上に色で描く（blog_013 §3「縮めるとは」用）。
  python3 run/make_krylov_basis_fig.py → docs/img/krylov_basis.png
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import krylov_rom as kr
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
import pyvista as pv
from make_unobserved_points_fig import grid

IMG = os.path.join(ROOT, "docs", "img")


def main():
    fem = kr.build_fem(1.0); V = kr.krylov_basis(fem, 4)
    g = grid(); pv.OFF_SCREEN = True; shots = []
    for j in range(4):
        v = V[:, j] / np.abs(V[:, j]).max()
        if v[np.argmax(np.abs(v))] < 0: v = -v
        gg = g.copy(); gg.point_data["v"] = v
        pl = pv.Plotter(off_screen=True, window_size=(600, 700)); pl.set_background("white")
        pl.add_mesh(gg, scalars="v", cmap="coolwarm", clim=[-1, 1], show_scalar_bar=False, lighting=False)
        pl.add_mesh(gg.extract_feature_edges(), color="#333333", line_width=1.2)
        pl.camera_position = [(230, -260, 230), (0, 0, 50), (0, 0, 1)]; pl.camera.zoom(1.3)
        img = pl.screenshot(return_img=True); pl.close()
        nz = np.where((img[:, :, :3] < 245).any(axis=2)); shots.append(img[nz[0].min() - 5:nz[0].max() + 5, nz[1].min() - 5:nz[1].max() + 5])
    titles = ["1 本目：$G^{-1}f$\n（入れ続けた後の形）", "2 本目：$G^{-1}C$×1 本目\n（そこへ向かう遅れ）", "3 本目", "4 本目"]
    plt.rcParams.update({"font.size": 16})
    fig, axs = plt.subplots(1, 4, figsize=(16, 5.6), gridspec_kw=dict(wspace=0.25))
    for ax, im, tl in zip(axs, shots, titles):
        ax.imshow(im); ax.axis("off"); ax.set_title(tl, fontsize=17)
    cax = fig.add_axes([0.3, 0.04, 0.4, 0.035])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-1, 1), cmap="coolwarm"), cax=cax, orientation="horizontal")
    cb.set_ticks([-1, 0, 1]); cb.set_label("基底の値（一番大きい所を 1）", fontsize=16)
    fig.savefig(os.path.join(IMG, "krylov_basis.png"), dpi=100, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("wrote krylov_basis.png")


if __name__ == "__main__":
    main()
