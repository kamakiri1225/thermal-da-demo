"""「モードを足すと真値に近づく」を、全部同じ実温度スケールで見せる（直感的なPOD図）.

平均場 → +モード1 → +モード1,2 → +モード1,2,3 → 真値 を、共通のカラーバーで並べる。
各段の再構成誤差(RMSE)も表示。「2モードでほぼ真値」が一目で分かる。

出力: docs/img/pod_reconstruction.png
再現: OMP_NUM_THREADS=4 python3 run/make_pod_reconstruction_fig.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
from scipy.spatial import cKDTree
from select_points_qdeim import load_snapshots
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img"); K=273.15


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    ts,_,X=load_snapshots()
    # 代表時刻＝加熱ピーク付近(t=300)
    kt=int(np.argmin(np.abs(ts-300))); u=X[:,kt]
    a=U.T@(u-mean)                       # モード振幅 a_k = φ_k^T (u - mean)
    # 段階再構成（真値を左端の基準に置き、右へ「平均場→+モード」で近づける）
    recons=[("真値(OpenFOAM)",u),
            ("平均場のみ",mean.copy()),
            ("＋モード1",mean+U[:,:1]@a[:1]),
            ("＋モード1,2",mean+U[:,:2]@a[:2]),
            ("＋モード1,2,3",mean+U[:,:3]@a[:3])]
    errs={lab:float(np.sqrt(((f-u)**2).mean())) for lab,f in recons}

    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    _,near=cKDTree(Cc).query(coords)
    clim=[float(u.min()-K),float(u.max()-K)]         # 全部共通スケール(真値の範囲)
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; TMP=tempfile.mkdtemp(); shots=[]
    for lab,f in recons:
        g=ug.copy(); g.point_data["T"]=f[near]-K
        pl=pv.Plotter(off_screen=True,window_size=(460,540))
        pl.add_mesh(g,scalars="T",cmap="turbo",clim=clim,n_colors=14,show_scalar_bar=False)
        pl.camera_position=CAM; pl.set_background("white")
        p=os.path.join(TMP,lab[:3]+".png"); pl.screenshot(p); pl.close(); shots.append((p,lab))
    fig,axes=plt.subplots(1,5,figsize=(18,5))
    for ax,(p,lab) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off")
        e=errs[lab]
        sub=f"\n(誤差RMSE {e:.3f}K)" if "真値" not in lab else "\n(左端＝基準)"
        ax.set_title(lab+sub,fontsize=13)
    # 共通カラーバー
    sm=plt.cm.ScalarMappable(cmap="turbo",norm=plt.Normalize(*clim))
    cb=fig.colorbar(sm,ax=axes,fraction=0.02,pad=0.01); cb.set_label("温度 [degC]（全パネル共通）",fontsize=11)
    fig.suptitle(f"モードを足すと真値に近づく（t={ts[kt]:g}s, 全パネル同じ温度スケール）"
                 "— 平均場＋2モードでほぼ真値",fontsize=15,weight="bold")
    fig.savefig(os.path.join(IMG,"pod_reconstruction.png"),dpi=130,bbox_inches="tight")
    plt.close(fig)
    print("[recon] 段階再構成の誤差:",{k:round(v,3) for k,v in errs.items()})
    print("[recon] wrote docs/img/pod_reconstruction.png")


if __name__=="__main__": main()
