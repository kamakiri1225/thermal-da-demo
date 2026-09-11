"""QoIの定義点(Point A / Point O)の位置図と、点名を明記したQoI時刻歴を作る.

QoI = Uz(Point A) − Uz(Point O)
  Point A: ヒータ側 上面 (+X, x=+28mm, z=100.5mm)
  Point O: 反ヒータ側 上面 (−X, x=−28mm, z=100.5mm)
(KinvH検証の Point_A − Point_O と同型の定義)

出力:
  docs/img/qoi_points.png              … A/Oの位置図(FEMメッシュ面)
  docs/img/selection_qoi_timehistory.png … 点名入りタイトルで再生成
  docs/img/selection_qoi_error.png       … 同
"""
from __future__ import annotations
import os, sys
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, vtk
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore.calibrate import load_calibrated
from dacore.cht_rom import N_NODES, T_AIR_K
from run.displacement_point_selection import build_M_all, NR,NTH,NZ,R_IN,R_OUT,H
from run.make_selection_gif import run_da_traj

K=273.15; IMG=os.path.join(ROOT,"docs","img")
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])
TITLE_QOI="QoI = Uz(Point A: ヒータ側上面+X) − Uz(Point O: 反対側上面−X)"

def main():
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated()
    coords,M=build_M_all()
    qa=int(np.linalg.norm(coords-A_XYZ,axis=1).argmin())
    qo=int(np.linalg.norm(coords-O_XYZ,axis=1).argmin())
    # 選点は FrontISTR(KinvH W=K^-1 H) 行感度による確定値 (run/kinvh_sensitivity.py)
    kv=np.load(os.path.join(ROOT,"results","kinvh_sensitivity.npz"))
    hi=[int(i) for i in kv["hi"]]; lo=[int(i) for i in kv["lo"]]
    ts,da_hi,tr=run_da_traj(cfg,calib,M[hi]); _,da_lo,_=run_da_traj(cfg,calib,M[lo])
    wq=M[qa]-M[qo]                      # QoI行 [µm/K]
    qoi=lambda A:(A-T_AIR_K)@wq

    # --- A/O 位置図 ---
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    pl=pv.Plotter(off_screen=True,window_size=(950,820))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.55,show_edges=False)
    pl.add_mesh(pv.Sphere(radius=0.005,center=coords[qa]),color="red")
    pl.add_mesh(pv.Sphere(radius=0.005,center=coords[qo]),color="blue")
    pl.add_point_labels([coords[qa]],["Point A (heater side, +X, top)"],font_size=22,
                        text_color="red",shape=None,always_visible=True)
    pl.add_point_labels([coords[qo]],["Point O (opposite side, -X, top)"],font_size=22,
                        text_color="blue",shape=None,always_visible=True)
    pl.add_text("QoI = Uz(Point A) - Uz(Point O)   [top-face tilt]",
                font_size=17,color="black")
    pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    pl.screenshot(os.path.join(IMG,"qoi_points.png")); pl.close()
    print("[qoi] wrote qoi_points.png")

    # --- QoI時刻歴(点名入り) ---
    fig,ax=plt.subplots(figsize=(11,5.4)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,qoi(tr),"k--",lw=2.5,label="真値")
    ax.plot(ts,qoi(da_hi),"-o",ms=4,color="tab:red",label="高感度2点でDA")
    ax.plot(ts,qoi(da_lo),"-s",ms=4,color="tab:blue",label="低感度2点でDA")
    ax.set_xlabel("time [s]"); ax.set_ylabel("Uz(A) − Uz(O) [µm]")
    ax.set_title(TITLE_QOI)
    ax.grid(alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"selection_qoi_timehistory.png"),dpi=140); plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,5)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,np.abs(qoi(da_hi)-qoi(tr)),"-o",ms=4,color="tab:red",label="高感度2点")
    ax.plot(ts,np.abs(qoi(da_lo)-qoi(tr)),"-s",ms=4,color="tab:blue",label="低感度2点")
    ax.set_yscale("log"); ax.set_xlabel("time [s]")
    ax.set_ylabel("|Uz(A)−Uz(O) の誤差| [µm]")
    ax.set_title("QoI誤差の時刻歴 — "+TITLE_QOI)
    ax.grid(alpha=0.3,which="both"); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"selection_qoi_error.png"),dpi=140); plt.close(fig)
    print("[qoi] regenerated QoI figures with point names")

if __name__=="__main__": main()
