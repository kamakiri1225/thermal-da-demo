"""校正した5点ROMは、校正に使っていない「無関係な点」でも温度が合うのか検証する.

校正図(rom_calib_fit.png)は5点しか見ていない。ここでは
  ① 5点ROMを0→600s前進（同化なし。発熱は真値）
  ② その5点温度を gappy-POD で全20,696セルへ復元
  ③ OpenFOAM(CHT)の答えと、5点から遠い「無関係な点」で突き合わせる
出力: docs/img/rom_unrelated_points.png, results/rom_unrelated_points.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/rom_unrelated_points_check.py
"""
from __future__ import annotations
import os, sys, json, tempfile, importlib.util
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
spec=importlib.util.spec_from_file_location('q', os.path.join(HERE,'select_points_qdeim.py'))
q=importlib.util.module_from_spec(spec); spec.loader.exec_module(q)
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KC=273.15; NPT=5; DT=2.0; NPICK=6


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"])
    heat_node=int(d["heat_node"]); xyz5=d["xyz"]
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float)
    Cc=kv["cell_centres"]; pod=kv["cell_idx"]
    ts, _Cc2, X = q.load_snapshots()               # OpenFOAMの答え (20696, 121)

    # ① 5点ROMを前進（同化なし・発熱は真値の倍率1.0）
    T=np.full(NPT, rg.T_AIR_K); Yrom=[T.copy()]
    for a,b in zip(ts[:-1], ts[1:]):
        _,tr=rg.integrate_single(T,C,Kmat,h,1.0,heat_node,a,b,DT); T=tr[-1]
        Yrom.append(T.copy())
    Yrom=np.array(Yrom).T                          # (5, nt)

    # ② gappy-PODで全セルへ復元
    UPp=np.linalg.pinv(U[pod,:])
    Xhat=mean[:,None] + U @ (UPp @ (Yrom - mean[pod,None]))

    # ③ 全セル / 5点以外 / 5点から遠い点 で突き合わせ
    err=Xhat-X
    others=np.setdiff1d(np.arange(len(X)), pod)
    rmse_all=float(np.sqrt((err**2).mean()))
    rmse_5=float(np.sqrt((err[pod]**2).mean()))
    rmse_oth=float(np.sqrt((err[others]**2).mean()))
    dmin=np.linalg.norm(Cc[:,None,:]-xyz5[None,:,:],axis=2).min(axis=1)   # 5点までの最短距離
    far=others[np.argsort(dmin[others])[::-1]]
    print(f"全20,696セル RMSE {rmse_all:.4f} K / 校正に使った5点 {rmse_5:.4f} K / "
          f"それ以外20,691セル {rmse_oth:.4f} K")
    print(f"5点から最も遠いセルの距離 {dmin[far[0]]*1000:.1f} mm")

    # 「無関係な点」を、遠い順かつ互いに離れた6点で選ぶ
    picks=[]
    for c in far:
        if all(np.linalg.norm(Cc[c]-Cc[p])>0.018 for p in picks):
            picks.append(int(c))
        if len(picks)==NPICK: break
    cell_rmse=np.sqrt((err**2).mean(axis=1))

    # ── ParaView風：どこの点かを示す ──
    import pyvista as pv, vtk, cylinder_mesh
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([p for _n,p in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,cc in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in cc)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    _,near=cKDTree(Cc).query(coords)
    ug.point_data["e"]=cell_rmse[near]*1000        # mK
    TMP=tempfile.mkdtemp(); shots=[]
    for cam,lab in [([(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)],"ヒータ側(+X)から"),
                    ([(-0.26,0.24,0.20),(0,0,0.05),(0,0,1)],"反対側(−X)から")]:
        pl=pv.Plotter(off_screen=True,window_size=(600,640))
        pl.add_mesh(ug,scalars="e",cmap="turbo",n_colors=14,opacity=0.94,
            scalar_bar_args={"title":"ROM+復元の誤差 [mK]","title_font_size":15,"label_font_size":12})
        for i in range(NPT):
            pl.add_mesh(pv.Sphere(radius=0.0022,center=xyz5[i]),color="white")
        for k,c in enumerate(picks):
            pl.add_mesh(pv.Sphere(radius=0.0022,center=Cc[c]),color="magenta")
            pl.add_point_labels([Cc[c]],[f"Q{k+1}"],font_size=17,text_color="magenta",
                                shape=None,always_visible=True)
        pl.camera_position=cam; pl.set_background("white"); pl.camera.zoom(1.35)
        p=os.path.join(TMP,lab[:3]+".png"); pl.screenshot(p); pl.close(); shots.append((p,lab))

    # ── 図 ── 上2段＝時刻歴6点、下段＝どこの点かの3D
    fig=plt.figure(figsize=(16.0,13.2))
    gs=fig.add_gridspec(3,6,height_ratios=[1,1,1.5],hspace=.52,wspace=.55)
    for k,c in enumerate(picks):
        ax=fig.add_subplot(gs[k//3, (k%3)*2:(k%3)*2+2])
        ax.axvspan(0,300,color="orange",alpha=.07)
        ax.plot(ts,X[c]-KC,color="black",lw=5.0,alpha=.40,label="OpenFOAM（真値）")
        ax.plot(ts,Xhat[c]-KC,"--",color="tab:red",lw=2.2,dashes=(4,3),label="5点ROM＋復元")
        ax.set_title(f"Q{k+1}  ({Cc[c][0]*1000:.0f}, {Cc[c][1]*1000:.0f}, {Cc[c][2]*1000:.0f}) mm\n"
                     f"5点から {dmin[c]*1000:.0f} mm 離れた点　誤差 {cell_rmse[c]*1000:.1f} mK",
                     fontsize=12.5, weight="bold")
        ax.set_xlabel("時刻 [s]",fontsize=11); ax.set_ylabel("温度 [℃]",fontsize=11)
        ax.grid(alpha=.3); ax.tick_params(labelsize=10)
        if k==0: ax.legend(fontsize=11,loc="lower right")
    for j,(p2,lab) in enumerate(shots):
        ax=fig.add_subplot(gs[2, j*3:j*3+3])
        ax.imshow(plt.imread(p2)); ax.axis("off")
        ax.set_title(f"{lab}　白球＝校正に使った5点／マゼンタ＝無関係な検証点 Q1〜Q6",
                     fontsize=12.5, weight="bold")
    fig.suptitle("校正した5点ROMは、校正に使っていない点でも温度が合うか\n"
                 f"全20,696セル {rmse_all:.3f} K　｜　校正に使った5点 {rmse_5:.3f} K　｜　"
                 f"それ以外20,691セル {rmse_oth:.3f} K",fontsize=16,weight="bold")
    fig.text(0.5,0.008,"※ データ同化なし。校正済みROMを0→600sそのまま前進させ、"
             "5点の温度をgappy-PODで全場へ復元してOpenFOAMと比較した。",
             ha="center",fontsize=12,color="#444")
    fig.tight_layout(rect=[0,0.022,1,0.935])
    out=os.path.join(IMG,"rom_unrelated_points.png")
    fig.savefig(out,dpi=125); plt.close(fig)
    with open(os.path.join(RES,"rom_unrelated_points.json"),"w") as f:
        json.dump(dict(rmse_all=rmse_all,rmse_calib5=rmse_5,rmse_others=rmse_oth,
                       max_err_K=float(np.abs(err[others]).max()),
                       picks=[dict(cell=int(c),xyz_mm=[round(v*1000,1) for v in Cc[c]],
                                   dist_to_5pt_mm=round(float(dmin[c])*1000,1),
                                   rmse_mK=round(float(cell_rmse[c])*1000,2)) for c in picks]),
                  f,ensure_ascii=False,indent=2)
    print("wrote",out)


if __name__=="__main__":
    main()
