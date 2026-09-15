"""5点→全温度場の復元アニメ（blog_004用）.

各時刻で、校正済みROMが出す「5点の温度」だけから、PODモードを使って
全20,696セルの温度分布を復元する（gappy-POD）。遠いセルも"近くの点からの内挿"でなく
"モード経由"で決まることを、動く絵で見せる。

  a(t) = U5P^{-1} (T5(t) - mean5)          # 5点で5個のモード係数を解く
  T_hat(t) = mean + U @ a(t)               # 全セルをmean+モードの重ねで復元

出力: docs/img/blog_temp_recon_anim.gif
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_temp_recon_anim.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
sys.path.insert(0, ROOT)
from dacore import plots as _p
from dacore import rom_general as rg
import cylinder_mesh, vtk
from scipy.spatial import cKDTree
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results"); K=273.15


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    idx=kv["cell_idx"]
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])

    # ROMを0-600s回して5点温度履歴（ヒータは300sで自動OFF）
    T0=np.full(5, rg.T_AIR_K)
    ts,traj=rg.integrate_single(T0,C,Kmat,h,1.0,heat,0.0,600.0,2.0)   # (nt,5) [K]

    # 5点→係数→全場（gappy-POD 復元）
    U5P=U[idx,:]                       # (5,5) 選点でのモード値
    U5P_inv=np.linalg.inv(U5P)
    mean5=mean[idx]
    A=(traj-mean5)@U5P_inv.T           # (nt,5) 各時刻のモード係数 a(t)
    field=mean[None,:]+A@U.T           # (nt,20696) 復元温度場 [K]

    # 描画メッシュ（FEM円筒）＋最寄りセル対応
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,c in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in c)
    ug=pv.UnstructuredGrid(np.array(conn),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    _,near=cKDTree(Cc).query(coords)
    sensors=Cc[idx]                    # 5点の座標

    clim=[float(field.min()-K), float(field.max()-K)]
    frames=[]; sub=range(0,len(ts),8)  # 間引いてコマ数を減らす(約38コマ)
    for fi in sub:
        tv=ts[fi]; fld=field[fi]-K      # ℃
        g=ug.copy(); g.point_data["T"]=fld[near]
        pl=pv.Plotter(off_screen=True,window_size=(720,760))
        # 5点が透けて見えるよう温度コンタは半透明に
        pl.add_mesh(g,scalars="T",cmap="turbo",clim=clim,n_colors=16,opacity=0.55,
                    scalar_bar_args={"title":"T [degC]","fmt":"%.1f"})
        # 復元に使った5点（マゼンタ球＝turboに無い色で目立たせる、少し大きめ）
        for s in sensors:
            pl.add_mesh(pv.Sphere(radius=0.0038,center=s),color="magenta")
        phase="heating (heater ON)" if tv<=300 else "cooling (heater OFF)"
        pl.add_text(f"t = {tv:4.0f} s   ({phase})\n5 sensors -> full field (POD)",
                    position="upper_left",font_size=11,color="black")
        pl.camera_position=[(0.26,-0.24,0.20),(0,0,0.05),(0,0,1)]
        pl.set_background("white"); pl.camera.zoom(0.85)
        frames.append(pl.screenshot(return_img=True)); pl.close()

    from PIL import Image
    ims=[Image.fromarray(f) for f in frames]
    out=os.path.join(IMG,"blog_temp_recon_anim.gif")
    ims[0].save(out,save_all=True,append_images=ims[1:],duration=140,loop=0,disposal=2)
    print("wrote",out,"frames",len(ims),"clim",np.round(clim,2))


if __name__=="__main__": main()
