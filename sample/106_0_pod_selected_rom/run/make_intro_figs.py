"""オープンCAE発表「はじめに」スライド用の2枚（左:温度 / 右:変位アニメ）.

左 intro_temp_history.png : 熱電対相当の測温点における温度の時間変化
右 intro_deform_anim.gif  : 変形の様子（誇張表示）と変位グラフを並べたアニメーション

出力: docs/img/intro_temp_history.png, docs/img/intro_deform_anim.gif
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_intro_figs.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG = os.path.join(ROOT, "docs", "img"); RES = os.path.join(ROOT, "results")
VTU = os.path.join(SAMPLE, "104_0_openfoam_frontistr_da_enkf", "paraview", "displacement_fields.vtu")
KC = 273.15; EXAG = 6000.0


def temp_fig():
    ev = np.load(os.path.join(RES, "fullsolver_evidence.npz"))
    t = ev["time"]; T = ev["temperature_truth_K"] - KC
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    ax.axvspan(0, 300, color="orange", alpha=.10)
    ax.text(150, T.max(), "加熱中", ha="center", va="top", fontsize=11, color="#b36a12", weight="bold")
    ax.text(450, T.min()+0.25, "冷却", ha="center", va="bottom", fontsize=11, color="#3b74b8", weight="bold")
    ax.plot(t, T[:, 0], "-o", color="tab:red", lw=2.6, ms=6, label="測温点A（熱源側）")
    ax.plot(t, T[:, 1], "-s", color="tab:blue", lw=2.6, ms=6, label="測温点B（反対側）")
    ax.fill_between(t, T[:, 1], T[:, 0], color="gray", alpha=.15)
    i = np.argmax(T[:, 0] - T[:, 1])
    ax.annotate(f"温度差 {T[i,0]-T[i,1]:.1f} K",
                xy=(t[i], (T[i, 0] + T[i, 1]) / 2), xytext=(t[i] + 60, T[:, 0].min() + 0.6),
                fontsize=11, weight="bold",
                arrowprops=dict(arrowstyle="-|>", color="k", lw=1.4))
    ax.set_xlabel("時間 [s]"); ax.set_ylabel("温度 [℃]")
    ax.set_title("温度の時間変化（測温点）", fontsize=13.5, weight="bold")
    ax.legend(fontsize=10.5, loc="lower right"); ax.grid(alpha=.3)
    fig.tight_layout()
    out = os.path.join(IMG, "intro_temp_history.png")
    fig.savefig(out, dpi=150); plt.close(fig); print("wrote", out)


def deform_anim():
    import pyvista as pv
    from PIL import Image
    pv.OFF_SCREEN = True
    try: pv.start_xvfb()
    except Exception: pass
    m = pv.read(VTU)
    m.set_active_vectors("U_truth_m"); m.set_active_scalars("Uz_truth_um")
    uz = np.asarray(m.point_data["Uz_truth_um"]); clim = [float(uz.min()), float(uz.max())]
    ev = np.load(os.path.join(RES, "fullsolver_evidence.npz"))
    tt = ev["time"].astype(float); uA = ev["u_truth_um"][:, 0]; uO = ev["u_truth_um"][:, 1]
    diff = uA - uO; amax = float(diff.max())
    tf = np.linspace(0, 600, 36)
    aA = np.interp(tf, tt, uA); aO = np.interp(tf, tt, uO); aD = np.interp(tf, tt, diff)

    frames = []
    for k, (tv, a) in enumerate(zip(tf, aD / amax)):
        pl = pv.Plotter(off_screen=True, window_size=(620, 700))
        pl.add_mesh(m.warp_by_vector("U_truth_m", factor=EXAG * a), scalars="Uz_truth_um",
                    cmap="coolwarm", clim=clim, show_edges=False,
                    scalar_bar_args={"title": "Uz [um]", "fmt": "%.1f"})
        pl.add_mesh(m, color="lightgray", opacity=0.10, show_edges=False)
        pl.add_text(f"t = {tv:4.0f} s\ndeformation x{EXAG:.0f}", position="upper_left",
                    font_size=12, color="black")
        pl.camera_position = [(0.34, -0.30, 0.22), (0, 0, 0.05), (0, 0, 1)]
        pl.set_background("white"); pl.camera.zoom(1.55)
        img3d = pl.screenshot(return_img=True); pl.close()
        msk = np.any(img3d[:, :, :3] < 245, axis=-1)      # 白余白をトリム
        ys, xs = np.where(msk); pad = 8
        img3d = img3d[max(0, ys.min()-pad):ys.max()+pad, max(0, xs.min()-pad):xs.max()+pad]

        fig, (a0, a1) = plt.subplots(1, 2, figsize=(11.4, 4.8),
                                     gridspec_kw={"width_ratios": [0.95, 1.25]})
        a0.imshow(img3d); a0.axis("off")
        a0.set_title("熱変形の様子（誇張表示）", fontsize=13, weight="bold")
        a1.axvspan(0, 300, color="orange", alpha=.10)
        a1.plot(tf, aA, color="tab:red", lw=2.4, label="Uz(A) 熱源側")
        a1.plot(tf, aO, color="tab:blue", lw=2.4, label="Uz(O) 反対側")
        a1.plot(tf, aD, "--", color="tab:green", lw=2.4, label="差 A−O（反り）")
        a1.plot(tf[k], aA[k], "o", color="tab:red", ms=10, zorder=5)
        a1.plot(tf[k], aO[k], "s", color="tab:blue", ms=9, zorder=5)
        a1.plot(tf[k], aD[k], "D", color="tab:green", ms=9, zorder=5)
        a1.axvline(tf[k], color="0.5", lw=1.0, ls=":")
        a1.set_xlim(0, 600); a1.set_ylim(min(aO.min(), aD.min()) - .4, aA.max() + .5)
        a1.set_xlabel("時間 [s]"); a1.set_ylabel("上向き変位 [µm]")
        a1.set_title("変位の時間変化", fontsize=13, weight="bold")
        a1.legend(fontsize=10, loc="upper left"); a1.grid(alpha=.3)
        fig.tight_layout()
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
        frames.append(Image.fromarray(buf.copy()))
        plt.close(fig)

    out = os.path.join(IMG, "intro_deform_anim.gif")
    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=140, loop=0, disposal=2)
    print("wrote", out, f"({len(frames)}フレーム)")


if __name__ == "__main__":
    temp_fig(); deform_anim()
