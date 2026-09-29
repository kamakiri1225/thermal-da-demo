"""熱感度分布が「何を基準に、どちらが動いた量か」を説明する図.

行感度 |W_i| = Σ_j |Wz[i,j]| は
  「底面（固定・変位ゼロ）を基準に、全節点の温度が1 K上がったとき
    その節点が上下方向にどれだけ動くか [µm/K]」。
つまり基準は固定面であって、隣の節点ではない。

① 基準の説明：底面が拘束されているので変位は下から積み上がる（高さ方向プロファイル）
② 空間分布：どこがよく動くか（底面固定を明示）
③ 実際に測るのは2点の差：A−O の相対変位（片方だけだと全体の伸びも拾う）

出力: docs/img/sensitivity_reference.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_sensitivity_reference_fig.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
IMG=os.path.join(ROOT,"docs","img")
KINVH=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","kinvh_sensitivity.npz")
H=0.1005; A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    kv=np.load(KINVH, allow_pickle=True)
    rs=kv["row_sens"]; valid=kv["valid"]
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    co=np.array([p for _n,p in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,cc in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in cc)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),co)
    field=np.where(valid, rs, np.nan)
    vmax=float(np.nanpercentile(field,99.5))
    ug.point_data["W"]=np.nan_to_num(field, nan=0.0)

    # ── ① 高さ方向プロファイル（底面基準で積み上がる）──
    z=co[:,2]; zs_all=np.unique(np.round(z,6))
    prof=[]; zs=[]
    for zz in zs_all:                       # 固定部近傍は除外域なので NaN になる層がある
        v=field[np.isclose(z,zz)]
        if np.isfinite(v).any():
            prof.append(float(np.nanmean(v))); zs.append(float(zz))
    prof=np.array(prof); zs=np.array(zs)
    pmax=float(prof.max())

    # ── ② 3D分布（底面固定を明示）──
    pl=pv.Plotter(off_screen=True, window_size=(720,820))
    pl.add_mesh(ug, scalars="W", cmap="turbo", clim=[0,vmax], n_colors=14,
                scalar_bar_args={"title":"|W| [um/K]","title_font_size":18,
                                 "label_font_size":15,"fmt":"%.1f"})
    pl.add_mesh(pv.Disc(center=(0,0,-0.0015), inner=0.020, outer=0.0375, normal=(0,0,1),
                        r_res=2, c_res=48), color="#333333")
    pl.add_point_labels([np.array([0.0,0.0,-0.012])],["FIXED (Uz = 0)"],
                        font_size=30, text_color="#111111", shape=None,
                        always_visible=True, bold=True)
    for xyz,lab,col in [(A_XYZ,"A","#C0392B"),(O_XYZ,"O","#2e7d32")]:
        pl.add_mesh(pv.Sphere(radius=0.0036, center=xyz), color=col)
        pl.add_point_labels([xyz+np.array([0,0,0.008])],[lab],font_size=34,
                            text_color=col, shape=None, always_visible=True, bold=True)
    pl.camera_position=[(0.28,-0.25,0.17),(0,0,0.05),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.18)
    img=pl.screenshot(return_img=True); pl.close()
    m=np.any(img[:,:,:3]<246,axis=-1); ys,xs=np.where(m); p=6
    img=img[max(0,ys.min()-p):ys.max()+p, max(0,xs.min()-p):xs.max()+p]

    # ── 図 ──
    fig=plt.figure(figsize=(18.0,6.6))
    gs=fig.add_gridspec(1,3,width_ratios=[1.02,0.86,1.12],wspace=.30)

    a0=fig.add_subplot(gs[0])
    a0.plot(prof, zs*1000, "-o", color="#1F4E9C", lw=2.8, ms=5)
    a0.axhline(0, color="#333", lw=3)
    a0.text(pmax*0.98, 2.5, "底面＝固定（$U_z$=0）＝ここが基準",
            ha="right", fontsize=12, weight="bold", color="#333")
    a0.fill_betweenx([0,5],0,pmax*1.05, color="#888", alpha=.18)
    a0.set_xlim(0, pmax*1.05); a0.set_ylim(-6, 108)
    a0.set_xlabel("その高さの平均 $|W|$ [µm/K]", fontsize=12.5)
    a0.set_ylabel("高さ z [mm]", fontsize=12.5)
    a0.set_title("① 基準は「固定した底面」\n"
                 "膨張は下から積み上がるので、上へ行くほど大きく動く",
                 fontsize=13, weight="bold")
    a0.grid(alpha=.3)

    a1=fig.add_subplot(gs[1]); a1.imshow(img); a1.axis("off")
    a1.set_title("② 全節点が1 K上がったときの\n各点の上下動（固定面基準）",
                 fontsize=13, weight="bold")

    a2=fig.add_subplot(gs[2]); a2.axis("off")
    a2.set_title("③ 実際に測るのは「2点の相対変位」", fontsize=13, weight="bold")
    a2.text(0.02,0.86,"①②は「固定面を基準にした絶対変位」。\n"
                      "ただし機械で測れるのはふつう2点の差です。",
            fontsize=12.5, transform=a2.transAxes)
    a2.text(0.02,0.655,r"$\Delta u = U_z(A)-U_z(O)$",
            fontsize=18, transform=a2.transAxes, color="#C0392B")
    a2.text(0.02,0.50,"A＝上面のヒータ側(+X)、O＝上面の反対側(−X)。\n"
                      "どちらも同じ高さなので、"
                      "全体が一様に伸びた分は引き算で消える。",
            fontsize=12.5, transform=a2.transAxes)
    a2.text(0.02,0.335,"→ 残るのは「ヒータ側だけが余計に伸びた分」＝上面の傾き（反り）。\n"
                       "　 加工精度に直結するのはこちら。",
            fontsize=12.5, transform=a2.transAxes, color="#C0392B", weight="bold")
    a2.text(0.02,0.155,"片方だけ測ると、室温が上がっただけの一様な伸びまで\n"
                       "拾ってしまい、「反り」と区別できない。",
            fontsize=12, transform=a2.transAxes, color="#556")
    a2.text(0.02,0.02,"※ 本研究の注目量 $U_z(A)-U_z(O)$ はこの相対変位。\n"
                      "　 センサ選定に使う $|W|$ は②の絶対変位の大きさ（動ける点かどうか）。",
            fontsize=11.5, transform=a2.transAxes, color="#556")

    fig.suptitle("熱感度 $W$ の分布は「何を基準に、どちらが動いた量」か",
                 fontsize=15.5, weight="bold")
    fig.tight_layout(rect=[0,0,1,0.905])
    out=os.path.join(IMG,"sensitivity_reference.png")
    fig.savefig(out, dpi=130); plt.close(fig); print("wrote", out)
    print(f"  底面付近 {prof[0]:.3f} → 上端 {prof[-1]:.3f} µm/K")


if __name__=="__main__":
    main()
