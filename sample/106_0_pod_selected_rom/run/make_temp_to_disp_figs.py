"""温度→FrontISTR変位の可視化（blog_004用）: 変形アニメGIF＋変位の時系列グラフ.

- 変形形状: 104のFrontISTR節点変位 U_truth（displacement_fields.vtu）を使用。
- 振幅の時間変化: fullsolver_evidence の A/O 変位時系列（0-600s）で駆動（加熱で反り→冷却で戻る）。
- 変位は数µmなので誇張表示（EXAG倍）。底面(z=0)は固定なので変位ゼロ。

出力: docs/img/blog_deform_anim.gif, docs/img/blog_disp_timeseries.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_temp_to_disp_figs.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
VTU=os.path.join(SAMPLE,"104_0_openfoam_frontistr_da_enkf","paraview","displacement_fields.vtu")
EXAG=6000.0   # 変位の誇張倍率（µmオーダーを可視化）


def disp_timeseries():
    ev=np.load(os.path.join(RES,"fullsolver_evidence.npz"))
    t=ev["time"]; uA=ev["u_truth_um"][:,0]; uO=ev["u_truth_um"][:,1]
    fig,ax=plt.subplots(figsize=(8.6,5.0))
    ax.axvspan(0,300,color="orange",alpha=.07,label="加熱期(ヒータON)")
    ax.plot(t,uA,"-o",color="tab:red",lw=2.4,label="A（ヒータ側）上向き変位 Uz")
    ax.plot(t,uO,"-o",color="tab:blue",lw=2.4,label="O（反対側）上向き変位 Uz")
    ax.plot(t,uA-uO,"--s",color="tab:green",lw=2.2,label="A−O 差（反りの大きさ）")
    ax.axhline(0,color="k",lw=1.0)
    ax.text(150,5.4,"底面(z=0)は固定＝変位ゼロ。動くのは上ほど大きい",
            fontsize=10,color="dimgray",
            bbox=dict(boxstyle="round,pad=0.3",fc="white",ec="lightgray"))
    ax.set_xlabel("time [s]"); ax.set_ylabel("変位 [µm]")
    ax.set_title("FrontISTRで計算した熱変位の時系列（片側加熱→冷却）",fontsize=13,weight="bold")
    ax.grid(alpha=.3); ax.legend(fontsize=10)
    fig.tight_layout(); out=os.path.join(IMG,"blog_disp_timeseries.png")
    fig.savefig(out,dpi=140); plt.close(fig); print("wrote",out)


def deform_anim():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    m=pv.read(VTU)
    m.set_active_vectors("U_truth_m"); m.set_active_scalars("Uz_truth_um")
    clim=[float(np.asarray(m.point_data["Uz_truth_um"]).min()),
          float(np.asarray(m.point_data["Uz_truth_um"]).max())]
    ev=np.load(os.path.join(RES,"fullsolver_evidence.npz"))
    tt=ev["time"].astype(float); diff=ev["u_truth_um"][:,0]-ev["u_truth_um"][:,1]
    amax=diff.max()
    # 滑らかに: 0-600sを40フレームに内挿
    tf=np.linspace(0,600,40)
    ampf=np.interp(tf,tt,diff)/amax    # 0..1（反りの相対量）

    frames=[]
    for tv,a in zip(tf,ampf):
        pl=pv.Plotter(off_screen=True,window_size=(760,860))
        warped=m.warp_by_vector("U_truth_m", factor=EXAG*a)
        pl.add_mesh(warped, scalars="Uz_truth_um", cmap="coolwarm", clim=clim,
                    show_edges=False, scalar_bar_args={"title":"Uz [um]","fmt":"%.2f"})
        pl.add_mesh(m, color="lightgray", opacity=0.12, show_edges=False)  # 元形状を薄く重ねる
        phase="heating (heater ON)" if tv<=300 else "cooling (heater OFF)"
        pl.add_text(f"t = {tv:4.0f} s   ({phase})\ndeformation x{EXAG:.0f}",
                    position="upper_left", font_size=11, color="black")
        pl.camera_position=[(0.34,-0.30,0.22),(0,0,0.05),(0,0,1)]
        pl.set_background("white"); pl.camera.zoom(1.25)
        img=pl.screenshot(return_img=True); pl.close()
        frames.append(img)
    # PILでGIF化（imageio不要）
    from PIL import Image
    ims=[Image.fromarray(f) for f in frames]
    out=os.path.join(IMG,"blog_deform_anim.gif")
    ims[0].save(out, save_all=True, append_images=ims[1:], duration=125, loop=0, disposal=2)
    print("wrote",out)


def main():
    disp_timeseries()
    deform_anim()


if __name__=="__main__": main()
