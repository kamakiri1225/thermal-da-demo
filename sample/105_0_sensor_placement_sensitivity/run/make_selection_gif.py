"""変位観測点の選び方(高感度 vs 低感度)でDAの収束がどう違うかをGIFで見せる.

3パネル(左=高感度2点でDA / 中=低感度2点でDA / 右=真値)の温度分布を0→600sで動かす。
選点の違いが時間とともに効いてくる様子を1本で伝える。
出力: docs/img/selection_da_compare.gif
"""
from __future__ import annotations
import glob, os, subprocess, sys
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, vtk
from dacore.calibrate import load_calibrated
from dacore.enkf import enkf_update
from dacore.ensemble import init_ensemble, clip_params, forecast, N_AUG
from dacore.observations import generate_truth
from dacore.cht_rom import N_NODES, T_AIR_K, NODE_NAMES
from dacore.node_locations import NODE_XYZ
from run.displacement_point_selection import build_M_all, SIG_U, HOT

K=273.15
IMG=os.path.join(ROOT,"docs","img"); TMP=os.path.join(ROOT,"openfoam","selgif")


def run_da_traj(cfg,calib,Mrows_um,seed=20260908):
    """同化の平均5点温度の時刻歴を返す (nt,5)."""
    sigT=cfg["observation"]["noise_C"]
    rng_o=np.random.default_rng(seed); rng=np.random.default_rng(seed+1)
    times,truth=generate_truth(calib,cfg); dt=cfg["experiment"]["obs_interval_s"]
    cyc=np.arange(dt,cfg["experiment"]["t_end_s"]+1e-9,dt)
    idx=np.clip(np.searchsorted(times,cyc),0,len(times)-1)
    R=np.diag([sigT**2,SIG_U**2,SIG_U**2])
    Z=init_ensemble(cfg,rng); tp=0
    rec=[Z[:,:N_NODES].mean(0)]; tr=[truth[0]]; ts=[0.0]
    for t1,ci in zip(cyc,idx):
        Z=forecast(Z,tp,t1,calib,cfg,rng); Tt=truth[ci]
        Yf=np.zeros((len(Z),3)); Yf[:,0]=Z[:,HOT]
        Yf[:,1:]=(Z[:,:N_NODES]-T_AIR_K)@Mrows_um.T
        y=np.concatenate([[Tt[HOT]],Mrows_um@(Tt-T_AIR_K)])+rng_o.normal(0,np.sqrt(np.diag(R)))
        Z=enkf_update(Z,y,None,R,rng,inflation=cfg["filter"]["inflation"],Yf=Yf)
        clip_params(Z,cfg); tp=t1
        rec.append(Z[:,:N_NODES].mean(0)); tr.append(truth[ci]); ts.append(float(t1))
    return np.array(ts),np.array(rec),np.array(tr)


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    os.makedirs(TMP,exist_ok=True); os.makedirs(IMG,exist_ok=True)
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated()
    coords,M=build_M_all()
    sens=np.abs(M).sum(axis=1); free=coords[:,2]>0.005
    order=np.argsort(sens)[::-1]; of=[i for i in order if free[i]]
    hi=[of[0]]; v0=M[of[0]]/np.linalg.norm(M[of[0]])
    for i in of[1:]:
        v=M[i]/max(np.linalg.norm(M[i]),1e-12)
        if abs(v@v0)<0.9: hi.append(i); break
    lo=[i for i in order[::-1] if free[i] and sens[i]>1e-4][:2]
    ts,da_hi,tr=run_da_traj(cfg,calib,M[hi])
    _,da_lo,_=run_da_traj(cfg,calib,M[lo])

    # 円筒メッシュ + IDW (make_da_gifs と同じ)
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    mcoords=np.array([xyz for _n,xyz in mesh["nodes"]],dtype=float)
    id_row={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(id_row[n] for n in conn)
    grid=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mcoords)
    romn=np.array(list(NODE_XYZ.values())); W=np.zeros((len(mcoords),5))
    for i,p in enumerate(mcoords):
        d=np.linalg.norm(romn-p,axis=1)
        if d.min()<1e-9: W[i,d.argmin()]=1
        else: w=1/d**2; W[i]=w/w.sum()

    clim=[float(tr.min()-K),float(tr.max()-K)]
    panels=[("high-sensitivity 2pts",da_hi,hi,"red"),
            ("low-sensitivity 2pts",da_lo,lo,"blue"),("truth",tr,None,None)]
    for k,t in enumerate(ts):
        pl=pv.Plotter(off_screen=True,shape=(1,3),window_size=(1500,620))
        for j,(lab,traj,pts,col) in enumerate(panels):
            pl.subplot(0,j)
            g=grid.copy(); g.point_data["T"]=W@(traj[k]-K)
            clip=g.clip(normal="y",origin=(0,0,0.05025))
            pl.add_mesh(clip,scalars="T",cmap="turbo",clim=clim,
                        scalar_bar_args={"title":"T [degC]","title_font_size":18,"label_font_size":14})
            if pts is not None:
                for i in pts: pl.add_mesh(pv.Sphere(radius=0.0035,center=coords[i]),color=col)
            pl.add_text(lab,font_size=17,color="black")
            pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.add_text(f"t = {t:g} s",position="lower_edge",font_size=16,color="black")
        pl.set_background("white")
        pl.screenshot(os.path.join(TMP,f"s_{k:03d}.png")); pl.close()
        if k%8==0: print(f"[gif] frame {k+1}/{len(ts)}",flush=True)
    frames=sorted(glob.glob(os.path.join(TMP,"s_*.png")))
    out=os.path.join(IMG,"selection_da_compare.gif")
    subprocess.run(["convert","-delay","18","-loop","0",*frames,"-resize","1000x",
                    "-layers","Optimize","-colors","80",out],check=True)
    for f in frames: os.remove(f)
    print(f"[gif] wrote {os.path.relpath(out,ROOT)} ({len(frames)} frames)")

    # ---------- 変位で見せる: QoI時刻歴・誤差時刻歴・変形アニメ ----------
    import fistr_case
    from pathlib import Path as _P
    from fem.fem_obs import MATERIAL
    node_ids=[nid for nid,_x in mesh["nodes"]]
    work=os.path.join(ROOT,"openfoam","selgif_fistr"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(_P(work),4,48,20,0.020,0.0375,0.1005,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(_P(work))
    def fistr_U(nodeT5):
        fistr_case.write_cnt(_P(work),node_ids,W@nodeT5,
            reference_temperature=MATERIAL["reference_temperature_K"],
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(_P(work))
        disp=fistr_case.read_displacement(_P(work))
        return np.array([disp[n] for n in node_ids],dtype=float)
    # QoI: 上面の傾き = uz(+X上面) - uz(-X上面)  (並行セッションのA-O差と同型)
    qa=int(np.linalg.norm(mcoords-np.array([0.028,0,0.1005]),axis=1).argmin())
    qo=int(np.linalg.norm(mcoords-np.array([-0.028,0,0.1005]),axis=1).argmin())
    U_hi=[];U_lo=[];U_tr=[]
    for k in range(len(ts)):
        U_hi.append(fistr_U(da_hi[k])); U_lo.append(fistr_U(da_lo[k])); U_tr.append(fistr_U(tr[k]))
        if k%8==0: print(f"[disp] FrontISTR {k+1}/{len(ts)} x3",flush=True)
    U_hi=np.array(U_hi);U_lo=np.array(U_lo);U_tr=np.array(U_tr)
    qoi=lambda U: (U[:,qa,2]-U[:,qo,2])*1e6   # 傾きQoI [µm]
    import matplotlib.pyplot as plt
    # (1) QoI時刻歴
    fig,ax=plt.subplots(figsize=(10,5.2)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,qoi(U_tr),"k--",lw=2.5,label="真値")
    ax.plot(ts,qoi(U_hi),"-o",ms=4,color="tab:red",label="高感度2点でDA")
    ax.plot(ts,qoi(U_lo),"-s",ms=4,color="tab:blue",label="低感度2点でDA")
    ax.set_xlabel("time [s]"); ax.set_ylabel("上面の傾き Uz(+X)−Uz(−X) [µm]")
    ax.set_title("変位QoI(傾き)の時刻歴: 観測点の選び方で真値への追従が変わる")
    ax.grid(alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"selection_qoi_timehistory.png"),dpi=140); plt.close(fig)
    # (2) QoI誤差時刻歴
    fig,ax=plt.subplots(figsize=(10,4.8)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,np.abs(qoi(U_hi)-qoi(U_tr)),"-o",ms=4,color="tab:red",label="高感度2点")
    ax.plot(ts,np.abs(qoi(U_lo)-qoi(U_tr)),"-s",ms=4,color="tab:blue",label="低感度2点")
    ax.set_yscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("|QoI誤差| [µm]")
    ax.set_title("変位QoI誤差の時刻歴（対数）"); ax.grid(alpha=0.3,which="both"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"selection_qoi_error.png"),dpi=140); plt.close(fig)
    # (3) 変形アニメ 3パネル
    uzall=np.concatenate([U_tr[...,2].ravel()])*1e6
    clim2=[float(uzall.min()),float(uzall.max())]
    scale=0.02/max(abs(np.array(clim2)).max(),1e-9)*1e6
    for k,t in enumerate(ts):
        pl=pv.Plotter(off_screen=True,shape=(1,3),window_size=(1500,640))
        for j,(lab,U,pts,col) in enumerate([("high-sens 2pts DA",U_hi[k],hi,"red"),
                                            ("low-sens 2pts DA",U_lo[k],lo,"blue"),
                                            ("truth",U_tr[k],None,None)]):
            pl.subplot(0,j)
            g=grid.copy(); g.points=mcoords+U*scale; g["Uz [um]"]=U[:,2]*1e6
            pl.add_mesh(g,scalars="Uz [um]",cmap="coolwarm",clim=clim2,
                        scalar_bar_args={"title":"Uz [um]","title_font_size":18,"label_font_size":14})
            if pts is not None:
                for i in pts: pl.add_mesh(pv.Sphere(radius=0.0035,center=coords[i]),color=col)
            pl.add_text(lab,font_size=17,color="black")
            pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.add_text(f"t = {t:g} s  (warp x{scale:.0f})",position="lower_edge",font_size=16,color="black")
        pl.set_background("white")
        pl.screenshot(os.path.join(TMP,f"d_{k:03d}.png")); pl.close()
    frames=sorted(glob.glob(os.path.join(TMP,"d_*.png")))
    out2=os.path.join(IMG,"selection_deform_compare.gif")
    subprocess.run(["convert","-delay","18","-loop","0",*frames,"-resize","1000x",
                    "-layers","Optimize","-colors","80",out2],check=True)
    for f in frames: os.remove(f)
    print(f"[gif] wrote {os.path.relpath(out2,ROOT)}")

if __name__=="__main__": main()
