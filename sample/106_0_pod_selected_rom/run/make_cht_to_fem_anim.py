"""OpenFOAM(CHT温度) → FrontISTR(熱変形) の一気通貫を1枚のGIFで見せる.

左  : OpenFOAMのCHTで解いた固体温度場（節点へIDW補間して表示）
右  : その温度をFrontISTRに渡して得た熱変形（誇張表示）
同じ時刻で左右を進めるので「温度を引き継いで変形まで出る」ことが一目で伝わる。

出力: docs/img/cht_to_fem_anim.gif
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_cht_to_fem_anim.py
"""
from __future__ import annotations
import os, sys, tempfile, importlib.util
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
from scipy.spatial import cKDTree
from PIL import Image
spec=importlib.util.spec_from_file_location('q', os.path.join(HERE,'select_points_qdeim.py'))
q=importlib.util.module_from_spec(spec); spec.loader.exec_module(q)
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
VTU=os.path.join(SAMPLE,"104_0_openfoam_frontistr_da_enkf","paraview","displacement_fields.vtu")
KC=273.15; EXAG=6000.0; NFR=30


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass

    ts, Cc, X = q.load_snapshots()                 # OpenFOAM CHT スナップショット
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,c in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in c)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    _,near=cKDTree(Cc).query(coords)               # CFDセル中心 → FEM節点（IDWの最近傍版）
    Tn=(X-KC)[near,:]                              # (nnode, ntime) [degC]
    clim_T=[float(Tn.min()), float(Tn.max())]

    m=pv.read(VTU); m.set_active_vectors("U_truth_m")
    uz=np.asarray(m.point_data["Uz_truth_um"]); clim_U=[float(uz.min()),float(uz.max())]
    ev=np.load(os.path.join(RES,"fullsolver_evidence.npz"))
    tt=ev["time"].astype(float); diff=ev["u_truth_um"][:,0]-ev["u_truth_um"][:,1]
    amax=float(diff.max())

    tf=np.linspace(0,600,NFR)
    CAM=[(0.30,-0.27,0.20),(0,0,0.05),(0,0,1)]
    frames=[]
    for tv in tf:
        j=int(np.abs(ts-tv).argmin())
        g=ug.copy(); g.point_data["T"]=Tn[:,j]
        pl=pv.Plotter(off_screen=True,window_size=(560,640))
        pl.add_mesh(g,scalars="T",cmap="turbo",clim=clim_T,n_colors=16,
                    scalar_bar_args={"title":"T [degC]","fmt":"%.1f",
                                     "title_font_size":16,"label_font_size":13})
        pl.camera_position=CAM; pl.set_background("white"); pl.camera.zoom(1.5)
        left=pl.screenshot(return_img=True); pl.close()

        a=float(np.interp(tv,tt,diff))/amax
        w=m.warp_by_vector("U_truth_m",factor=EXAG*a)
        w.point_data["Uz"]=uz*a          # 変形量と同じ倍率で色も動かす
        pl=pv.Plotter(off_screen=True,window_size=(560,640))
        pl.add_mesh(w,scalars="Uz",cmap="coolwarm",clim=clim_U,
                    scalar_bar_args={"title":"Uz [um]","fmt":"%.1f",
                                     "title_font_size":16,"label_font_size":13})
        pl.add_mesh(m,color="lightgray",opacity=0.12)
        pl.camera_position=CAM; pl.set_background("white"); pl.camera.zoom(1.5)
        right=pl.screenshot(return_img=True); pl.close()

        fig,(a0,a1)=plt.subplots(1,2,figsize=(11.6,5.6))
        a0.imshow(left);  a0.axis("off")
        a0.set_title("① OpenFOAM（CHT）で解いた温度",fontsize=14,weight="bold")
        a1.imshow(right); a1.axis("off")
        a1.set_title(f"② FrontISTRの熱変形（×{EXAG:.0f}誇張）",fontsize=14,weight="bold")
        fig.suptitle(f"OpenFOAM → FrontISTR 一気通貫　t = {tv:4.0f} s"
                     + ("（加熱中）" if tv<=300 else "（冷却中）"),
                     fontsize=15,weight="bold")
        fig.text(0.5,0.028,"温度をIDWで節点へ写像し、そのままFrontISTRの熱弾性に渡している",
                 ha="center",fontsize=11.5,color="#333")
        fig.tight_layout(rect=[0,0.05,1,0.94]); fig.canvas.draw()
        frames.append(Image.fromarray(
            np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
        plt.close(fig)
        print(f"  frame t={tv:.0f}s",flush=True)

    out=os.path.join(IMG,"cht_to_fem_anim.gif")
    frames[0].save(out,save_all=True,append_images=frames[1:],duration=160,loop=0,disposal=2)
    print("wrote",out)


if __name__=="__main__":
    main()
