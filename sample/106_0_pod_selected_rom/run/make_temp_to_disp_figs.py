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


def disp_da_timeseries():
    """EnKF/OI比較用：FrontISTR真値を太線で重ねた変位差の時系列。"""
    d=np.load(os.path.join(RES,"da_compare_displacement.npz"), allow_pickle=True)
    t=d["time"]; truth=d["truth_u_um"][:,0]-d["truth_u_um"][:,1]
    est=d["estimate_u_um"]
    names=[str(x) for x in d["names"]]
    # 凡例をわかりやすく（どの点か・何の感度かを明記）
    rename={"同化なし(free run)":"同化なし",
            "温度1点(低感度)":"温度1点:P4底(dT/dQ小)",
            "温度1点(高感度)":"温度1点:P2ヒータ側(dT/dQ大)",
            "温度2点(高感度)":"温度2点:P2+P0",
            "温度2点+変位2点":"温度2点+変位2点(高W上端)"}
    names=[rename.get(n,n) for n in names]
    # 各設定について5 seedの平均を表示（ばらつきは薄い帯）
    colors=["0.45","tab:orange","tab:blue","tab:green","tab:red"]
    fig,ax=plt.subplots(figsize=(11.6,6.0))
    ax.axvspan(0,300,color="orange",alpha=.07,label="加熱期（ヒータON）")
    ax.plot(t,truth,color="black",lw=4.4,zorder=10,label="FrontISTR真値（A−O）")
    ax.annotate("真値（黒太線）",xy=(370,float(np.interp(370,t,truth))),xytext=(430,1.7),
                fontsize=12,weight="bold",
                arrowprops=dict(arrowstyle="-|>",color="black",lw=1.6))
    for i,(name,c) in enumerate(zip(names,colors)):
        y=est[i]  # (seed,time,2)  ※軸は(config,seed,time,2)
        diff=y[:,:,0]-y[:,:,1]
        mu=diff.mean(axis=0); sd=diff.std(axis=0)
        ax.plot(t,mu,lw=1.8,color=c,label=name)
        ax.fill_between(t,mu-sd,mu+sd,color=c,alpha=.06,linewidth=0)
    ax.axhline(0,color="k",lw=.8)
    ax.set_xlabel("時間 [s]"); ax.set_ylabel("変位差 Uz(A)−Uz(O) [µm]")
    ax.grid(alpha=.3)
    fig.suptitle("FrontISTR真値とデータ同化後の変位差（5 seed平均）\n"
                 "変位差RMSE: 温度1点 低感度1.19／高感度0.99µm（温度センサの場所差はほぼ出ない）"
                 "→ 変位2点を足すと0.17µm（約6倍改善）",
                 fontsize=12.5,weight="bold",y=0.99)
    hd,lb=ax.get_legend_handles_labels()
    fig.legend(hd,lb,loc="upper center",bbox_to_anchor=(0.5,0.945),ncol=4,fontsize=9,frameon=False)
    fig.tight_layout(rect=[0,0,1,0.865])
    out=os.path.join(IMG,"blog_disp_timeseries_truth_vs_da.png")
    fig.savefig(out,dpi=160,bbox_inches="tight"); plt.close(fig); print("wrote",out)


def disp_da_timeseries_points():
    """個別の Uz(A)・Uz(O) も並べる版（差だけだと誤差の相殺で構成差が見えないため）。"""
    d=np.load(os.path.join(RES,"da_compare_displacement.npz"), allow_pickle=True)
    t=d["time"]; tru=d["truth_u_um"]              # (time,2)
    est=d["estimate_u_um"]                         # (seed,cfg,time,2)
    names=[str(x) for x in d["names"]]
    rename={"同化なし(free run)":"同化なし",
            "温度1点(低感度)":"温度1点:P4底(dT/dQ小)",
            "温度1点(高感度)":"温度1点:P2ヒータ側(dT/dQ大)",
            "温度2点(高感度)":"温度2点:P2+P0",
            "温度2点+変位2点":"温度2点+変位2点(高W上端)"}
    names=[rename.get(n,n) for n in names]
    colors=["0.45","tab:orange","tab:blue","tab:green","tab:red"]
    ht=(t>0)&(t<=300)
    fig,axes=plt.subplots(1,3,figsize=(16.5,5.6))
    panels=[("Uz(A) ヒータ側・上面",0),("Uz(O) 反対側・上面",1),("差 Uz(A)−Uz(O)",None)]
    for ax,(ti,k) in zip(axes,panels):
        ax.axvspan(0,300,color="orange",alpha=.07)
        trv=tru[:,0]-tru[:,1] if k is None else tru[:,k]
        ax.plot(t,trv,color="black",lw=4.0,zorder=10)
        for i,c in enumerate(colors):
            y=est[i]
            v=(y[:,:,0]-y[:,:,1]) if k is None else y[:,:,k]
            mu=v.mean(axis=0)
            rm=np.sqrt(((v-trv[None,:])**2)[:,ht].mean())
            ax.plot(t,mu,lw=1.8,color=c,label=f"{names[i]}  RMSE {rm:.2f}µm")
        ax.set_title(ti,fontsize=12.5,weight="bold"); ax.grid(alpha=.3)
        ax.set_xlabel("時間 [s]"); ax.set_ylabel("変位 [µm]")
        ax.legend(fontsize=7.2,loc="lower center")
    fig.suptitle("個別の変位で見ると構成差は大きい：差(A−O)だけでは相殺で見えにくい\n"
                 "温度1点は“センサに近い側”しか合わない（P2ヒータ側→A良/O悪、P4底→O良/A悪）。変位2点(赤)は両点とも合う",
                 fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.88])
    out=os.path.join(IMG,"blog_disp_timeseries_truth_vs_da_points.png")
    fig.savefig(out,dpi=150,bbox_inches="tight"); plt.close(fig); print("wrote",out)


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
    disp_da_timeseries()
    disp_da_timeseries_points()
    deform_anim()


if __name__=="__main__": main()
