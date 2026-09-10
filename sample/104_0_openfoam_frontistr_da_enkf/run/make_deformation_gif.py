"""真値の熱膨張変形アニメ(0→600s)。ROM真値温度→FrontISTR→上面変位、単一パネル。
出力: docs/img/deformation.gif（presentation/assets にもコピー）
"""
from __future__ import annotations
import glob, os, subprocess, sys
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0,ROOT)
F102=os.path.abspath(os.path.join(ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
sys.path.insert(0,F102)
import cylinder_mesh, fistr_case, vtk
from dacore.calibrate import load_calibrated
from dacore.node_locations import NODE_XYZ
from dacore.cht_rom import N_NODES, T_AIR_K, integrate_single
from dacore.observations import true_params
from fem.fem_obs import MATERIAL
IMG=os.path.join(ROOT,"docs","img"); TMP=os.path.join(ROOT,"openfoam","def_frames")
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005; K=273.15

def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    os.makedirs(TMP,exist_ok=True); os.makedirs(IMG,exist_ok=True)
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]],dtype=float)
    node_ids=[nid for nid,_ in mesh["nodes"]]
    id_row={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(id_row[n] for n in conn)
    grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    romn=np.array(list(NODE_XYZ.values()))
    W=np.zeros((len(coords),len(romn)))
    for i,p in enumerate(coords):
        dd=np.linalg.norm(romn-p,axis=1)
        if dd.min()<1e-9:
            W[i,dd.argmin()]=1.0
        else:
            w=1/dd**2; W[i]=w/w.sum()
    # ROM真値温度(5点)を0-600s
    calib=load_calibrated(); p=true_params(calib)
    ts,traj=integrate_single(np.full(N_NODES,T_AIR_K),p,0.0,600.0,dt=1.0)
    times=list(range(0,601,20))
    work=os.path.join(ROOT,"openfoam","def_fistr"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    # 全フレームの変位を先に計算(スケール決定)
    Us=[]
    for t in times:
        Tnode=traj[int(round(t))]           # 5点温度[K]
        Tmesh=W@Tnode
        fistr_case.write_cnt(Path(work),node_ids,Tmesh,
            reference_temperature=MATERIAL["reference_temperature_K"],
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work))
        disp=fistr_case.read_displacement(Path(work))
        Us.append(np.array([disp[n] for n in node_ids]))
        print(f"[def] FrontISTR {len(Us)}/{len(times)} t={t}s",flush=True)
    Us=np.array(Us); uzmax=np.abs(Us[...,2]).max()*1e6
    clim=[float((Us[...,2].min())*1e6),float((Us[...,2].max())*1e6)]
    scale=0.02/max(uzmax,1e-9)*1e6
    for k,t in enumerate(times):
        pl=pv.Plotter(off_screen=True,window_size=(720,760))
        g=grid.copy(); g.points=coords+Us[k]*scale; g["Uz [um]"]=Us[k][:,2]*1e6
        pl.add_mesh(g,scalars="Uz [um]",cmap="coolwarm",clim=clim,
                    scalar_bar_args={"title":"Uz [um]","title_font_size":20,"label_font_size":16})
        phase="加熱中 (0-300s)" if t<300 else "遮断後 (300-600s)"
        pl.add_text(f"t = {t} s   {phase}  (warp x{scale:.0f})",position="upper_edge",font_size=15,color="black")
        pl.set_background("white"); pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.screenshot(os.path.join(TMP,f"d_{k:03d}.png")); pl.close()
    frames=sorted(glob.glob(os.path.join(TMP,"d_*.png")))
    out=os.path.join(IMG,"deformation.gif")
    subprocess.run(["convert","-delay","18","-loop","0",*frames,"-resize","560x",
                    "-layers","Optimize","-colors","96",out],check=True)
    for f in frames: os.remove(f)
    subprocess.run(["cp",out,os.path.join(ROOT,"presentation","assets","deformation.gif")],check=False)
    print(f"[def] wrote {os.path.relpath(out,ROOT)}")
if __name__=="__main__": main()
