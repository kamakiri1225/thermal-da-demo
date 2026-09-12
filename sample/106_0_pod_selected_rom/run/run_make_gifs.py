"""熱感度違い(高感度/低感度センサ)のデータ同化を、温度分布と変形のGIFで見せる.

3パネル（左=高感度センサ同化 / 中=低感度センサ同化 / 右=真値）を0→600sで動かす。
- 温度GIF: 各時刻の5点温度を gappy-POD で全場復元して描く
- 変形GIF: 各時刻の全場をFrontISTRに通した本物の熱変形（誇張表示）

出力(docs/img/): da_temp_compare.gif, da_deform_compare.gif
再現: OMP_NUM_THREADS=4 python3 run/run_make_gifs.py
"""
from __future__ import annotations
import os, sys, glob, subprocess, tempfile
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh, vtk, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img"); K=273.15
NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; INFL=1.02; SEED=20260913
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]


def da_traj(nodes, C, Kmat, heat_node, Ttr, cyc):
    rng=np.random.default_rng(SEED); rng_o=np.random.default_rng(SEED+7)
    Z=np.zeros((N_ENS,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    R=np.diag([SIG_T**2]*len(nodes)); rec=[Z[:,:NPT].mean(0).copy()]; tp=0.0
    for ci,tb in enumerate(cyc,1):
        Tn=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
        Z=Z.copy(); Z[:,:NPT]=Tn; tp=tb
        y=np.array(Ttr[ci][nodes])+rng_o.normal(0,SIG_T,len(nodes))
        Z=enkf_update(Z,y,None,R,rng,inflation=INFL,Yf=Z[:,nodes])
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        rec.append(Z[:,:NPT].mean(0).copy())
    return np.array(rec)


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat_node=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]; pod_cells=kv["cell_idx"]
    UP_pinv=np.linalg.pinv(U[pod_cells,:])
    def gappy(T5): a=(T5-mean[pod_cells])@UP_pinv.T; return mean+U@a

    # dT/dQ で高/低感度ノード
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens=(pert-base)/0.1; hi=int(np.argmax(sens)); lo=int(np.argmin(sens))

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tgrid=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    da_hi=da_traj([hi],C,Kmat,heat_node,Ttr,cyc)
    da_lo=da_traj([lo],C,Kmat,heat_node,Ttr,cyc)

    # メッシュ + gappy場をFEM節点へ
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    mco=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mco)
    _,near=cKDTree(Cc).query(mco)
    TMP=tempfile.mkdtemp()

    # ===== 温度GIF =====
    fields={"hi":[gappy(da_hi[k])[near] for k in range(len(tgrid))],
            "lo":[gappy(da_lo[k])[near] for k in range(len(tgrid))],
            "tr":[gappy(Ttr[k])[near] for k in range(len(tgrid))]}
    clim=[float(min(f.min() for f in fields["tr"])-K),float(max(f.max() for f in fields["tr"])-K)]
    for k,t in enumerate(tgrid):
        pl=pv.Plotter(off_screen=True,shape=(1,3),window_size=(1500,600))
        for j,(key,lab) in enumerate([("hi","高感度センサ同化"),("lo","低感度センサ同化"),("tr","真値")]):
            pl.subplot(0,j); g=grid.copy(); g.point_data["T"]=fields[key][k]-K
            clip=g.clip(normal="y",origin=(0,0,0.05025))
            pl.add_mesh(clip,scalars="T",cmap="turbo",clim=clim,n_colors=14,
                        scalar_bar_args={"title":"T [degC]","title_font_size":16,"label_font_size":12})
            pl.add_text(lab,font_size=15,color="black")
            pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.add_text(f"t = {t:g} s",position="lower_edge",font_size=15,color="black")
        pl.set_background("white"); pl.screenshot(os.path.join(TMP,f"t_{k:03d}.png")); pl.close()
    frames=sorted(glob.glob(os.path.join(TMP,"t_*.png")))
    subprocess.run(["convert","-delay","22","-loop","0",*frames,"-resize","1000x","-layers","Optimize",
                    "-colors","80",os.path.join(IMG,"da_temp_compare.gif")],check=True)
    for f in frames: os.remove(f)
    print("[gif] da_temp_compare.gif")

    # ===== 変形GIF (各フレームFrontISTR) =====
    work=os.path.join(ROOT,"openfoam","gif_fistr"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    def fistr_U(field):
        fistr_case.write_cnt(Path(work),node_ids,field,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        return np.array([disp[n] for n in node_ids],dtype=float)
    # フレーム間引き（変形は重いので5サイクルおき）
    ks=list(range(0,len(tgrid),2))
    Us={key:[] for key in ["hi","lo","tr"]}
    for k in ks:
        for key,fld in [("hi",fields["hi"][k]),("lo",fields["lo"][k]),("tr",fields["tr"][k])]:
            Us[key].append(fistr_U(fld))
        if k%6==0: print(f"[gif] deform FrontISTR {k}/{len(tgrid)}",flush=True)
    uzmax=max(abs(np.array(Us["tr"])[...,2]).max(),1e-12)
    scale=0.015/uzmax
    clim2=[float(np.array(Us["tr"])[...,2].min()*1e6),float(np.array(Us["tr"])[...,2].max()*1e6)]
    for fi,k in enumerate(ks):
        pl=pv.Plotter(off_screen=True,shape=(1,3),window_size=(1500,620))
        for j,(key,lab) in enumerate([("hi","高感度センサ同化"),("lo","低感度センサ同化"),("tr","真値")]):
            pl.subplot(0,j); Uk=Us[key][fi]; g=grid.copy(); g.points=mco+Uk*scale; g["Uz[um]"]=Uk[:,2]*1e6
            pl.add_mesh(g,scalars="Uz[um]",cmap="coolwarm",clim=clim2,n_colors=12,
                        scalar_bar_args={"title":"Uz [um]","title_font_size":16,"label_font_size":12})
            pl.add_text(lab,font_size=15,color="black")
            pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.add_text(f"t = {tgrid[k]:g} s  (warp x{scale:.0f})",position="lower_edge",font_size=15,color="black")
        pl.set_background("white"); pl.screenshot(os.path.join(TMP,f"d_{fi:03d}.png")); pl.close()
    frames=sorted(glob.glob(os.path.join(TMP,"d_*.png")))
    subprocess.run(["convert","-delay","28","-loop","0",*frames,"-resize","1000x","-layers","Optimize",
                    "-colors","80",os.path.join(IMG,"da_deform_compare.gif")],check=True)
    for f in frames: os.remove(f)
    print("[gif] da_deform_compare.gif")


if __name__=="__main__": main()
