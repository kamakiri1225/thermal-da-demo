"""証拠計算: 温度センサを「変位熱感度の高い場所」に置くと真値の分布を再現でき、
「低い場所」では再現できないことを、分布そのもので示す。

設定: 変位2点(M演算子)は常に同化。温度センサ1点の位置だけを
  (a) 感度マップ最大の場所(底面近くヒータ側)   (b) 感度マップ最小の場所
で置き分ける。任意位置の温度センサは IDW行 w(x)·T_node で観測演算子化する。

出力:
  docs/img/evidence_sensor_points.png   … 感度マップ上に2つのセンサ位置
  docs/img/evidence_field_compare.gif   … 温度分布の時刻歴(真値/高感度センサ/低感度センサ)
  docs/img/evidence_rmse.png            … RMSE時刻歴
  docs/img/evidence_qoi.png             … 変位QoI(傾き)の時刻歴
  results/evidence.csv
"""
from __future__ import annotations
import glob, os, subprocess, sys
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, vtk
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore.calibrate import load_calibrated
from dacore.enkf import enkf_update
from dacore.ensemble import init_ensemble, clip_params, forecast, N_AUG
from dacore.observations import generate_truth
from dacore.cht_rom import N_NODES, T_AIR_K
from dacore.node_locations import NODE_XYZ, HEIGHT, OUTER_RADIUS, INNER_RADIUS

K=273.15
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
TMP=os.path.join(ROOT,"openfoam","evgif")
Mum=np.load(os.path.join(ROOT,"config","M_frontistr_operator.npy"))*1e6
SIG_U=0.3
ROMN=np.array(list(NODE_XYZ.values()))


def idw_row(p):
    d=np.linalg.norm(ROMN-np.asarray(p),axis=1)
    if d.min()<1e-9:
        w=np.zeros(5); w[d.argmin()]=1; return w
    w=1/d**2; return w/w.sum()


def pick_sensor_points():
    d=np.load(os.path.join(RES,"sensitivity_uz.npz"))
    S=d["S_um"]; NTH_P=int(d["NTH_P"]); NZ_P=int(d["NZ_P"])
    mag=np.abs(S).mean(axis=0)               # (NTH_P,NZ_P)
    ah,bh=np.unravel_index(mag.argmax(),mag.shape)
    al,bl=np.unravel_index(mag.argmin(),mag.shape)
    r=0.5*(OUTER_RADIUS+INNER_RADIUS)
    def xyz(a,b):
        th=(a+0.5)/NTH_P*2*np.pi; z=(b+0.5)/NZ_P*HEIGHT
        return np.array([r*np.cos(th),r*np.sin(th),z]), th, z
    ph,thh,zh=xyz(ah,bh); plo,thl,zl=xyz(al,bl)
    return (ph,mag[ah,bh],thh,zh),(plo,mag[al,bl],thl,zl),(S,NTH_P,NZ_P,mag)


def run_traj(cfg,calib,wrow,seed=20260908):
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
        Yf=np.zeros((len(Z),3))
        Yf[:,0]=Z[:,:N_NODES]@wrow
        Yf[:,1:]=(Z[:,:N_NODES]-T_AIR_K)@Mum.T
        y=np.concatenate([[wrow@Tt],Mum@(Tt-T_AIR_K)])+rng_o.normal(0,np.sqrt(np.diag(R)))
        Z=enkf_update(Z,y,None,R,rng,inflation=cfg["filter"]["inflation"],Yf=Yf)
        clip_params(Z,cfg); tp=t1
        rec.append(Z[:,:N_NODES].mean(0)); tr.append(truth[ci]); ts.append(float(t1))
    return np.array(ts),np.array(rec),np.array(tr)


def main():
    os.makedirs(IMG,exist_ok=True); os.makedirs(TMP,exist_ok=True)
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated()
    (ph,sh,thh,zh),(plo,sl_,thl,zl),(S,NTH_P,NZ_P,mag)=pick_sensor_points()
    print(f"[ev] 高感度センサ位置: θ={np.degrees(thh):.0f}° z={zh*1000:.0f}mm 感度{sh:.3f}µm/K")
    print(f"[ev] 低感度センサ位置: θ={np.degrees(thl):.0f}° z={zl*1000:.0f}mm 感度{sl_:.4f}µm/K")
    w_hi=idw_row(ph); w_lo=idw_row(plo)
    ts,da_hi,tr=run_traj(cfg,calib,w_hi)
    _,da_lo,_=run_traj(cfg,calib,w_lo)

    rmse=lambda A: np.sqrt(((A-tr)**2).mean(axis=1))
    print(f"[ev] 最終RMSE: 高感度 {rmse(da_hi)[-1]:.4f} K / 低感度 {rmse(da_lo)[-1]:.4f} K")
    with open(os.path.join(RES,"evidence.csv"),"w") as f:
        f.write("case,theta_deg,z_mm,sens_um_per_K,final_rmse_K\n")
        f.write(f"high,{np.degrees(thh):.0f},{zh*1000:.0f},{sh:.4f},{rmse(da_hi)[-1]:.5f}\n")
        f.write(f"low,{np.degrees(thl):.0f},{zl*1000:.0f},{sl_:.4f},{rmse(da_lo)[-1]:.5f}\n")

    # --- センサ位置図(感度面+マーカー) ---
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    nth,nz=181,81
    th=np.linspace(0,2*np.pi,nth); zzl=np.linspace(0,HEIGHT,nz)
    TH,ZZ=np.meshgrid(th,zzl,indexing="ij")
    X=OUTER_RADIUS*np.cos(TH); Y=OUTER_RADIUS*np.sin(TH)
    a=np.minimum((TH.ravel(order='F')%(2*np.pi)/(2*np.pi)*NTH_P).astype(int),NTH_P-1)
    b=np.minimum((ZZ.ravel(order='F')/HEIGHT*NZ_P).astype(int),NZ_P-1)
    pl=pv.Plotter(off_screen=True,window_size=(950,800))
    surf=pv.StructuredGrid(X,Y,ZZ); surf["sens"]=mag[a,b]
    pl.add_mesh(surf,scalars="sens",cmap="viridis",
                scalar_bar_args={"title":"mean|dUz/dT| [um/K]","title_font_size":18,"label_font_size":14})
    pl.add_mesh(pv.Sphere(radius=0.0045,center=ph),color="red")
    pl.add_mesh(pv.Sphere(radius=0.0045,center=plo),color="deepskyblue")
    pl.add_text("temperature sensor candidates on sensitivity map\nred=HIGH sens  blue=LOW sens",
                font_size=16,color="black")
    pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    pl.screenshot(os.path.join(IMG,"evidence_sensor_points.png")); pl.close()

    # --- 温度分布GIF(真値/高感度/低感度) ---
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    mcoords=np.array([xyz for _n,xyz in mesh["nodes"]],dtype=float)
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    grid=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mcoords)
    Wm=np.zeros((len(mcoords),5))
    for i,p in enumerate(mcoords):
        d=np.linalg.norm(ROMN-p,axis=1)
        if d.min()<1e-9: Wm[i,d.argmin()]=1
        else: w=1/d**2; Wm[i]=w/w.sum()
    clim=[float(tr.min()-K),float(tr.max()-K)]
    marks=[("truth",tr,None,None),
           ("HIGH-sens temp sensor + 2 disp",da_hi,ph,"red"),
           ("LOW-sens temp sensor + 2 disp",da_lo,plo,"deepskyblue")]
    for k,t in enumerate(ts):
        pl=pv.Plotter(off_screen=True,shape=(1,3),window_size=(1500,620))
        for j,(lab,traj,pt,col) in enumerate(marks):
            pl.subplot(0,j)
            g=grid.copy(); g.point_data["T"]=Wm@(traj[k]-K)
            clip=g.clip(normal="y",origin=(0,0,0.05025))
            pl.add_mesh(clip,scalars="T",cmap="turbo",clim=clim,
                        scalar_bar_args={"title":"T [degC]","title_font_size":18,"label_font_size":14})
            if pt is not None: pl.add_mesh(pv.Sphere(radius=0.004,center=pt),color=col)
            pl.add_text(lab,font_size=15,color="black")
            pl.camera_position=[(0.24,-0.20,0.20),(0,-0.01,0.05),(0,0,1)]
        pl.add_text(f"t = {t:g} s",position="lower_edge",font_size=16,color="black")
        pl.set_background("white")
        pl.screenshot(os.path.join(TMP,f"e_{k:03d}.png")); pl.close()
    frames=sorted(glob.glob(os.path.join(TMP,"e_*.png")))
    out=os.path.join(IMG,"evidence_field_compare.gif")
    subprocess.run(["convert","-delay","18","-loop","0",*frames,"-resize","1000x",
                    "-layers","Optimize","-colors","80",out],check=True)
    for f in frames: os.remove(f)
    print(f"[ev] wrote {os.path.relpath(out,ROOT)}")

    # --- RMSE時刻歴 / 変位QoI時刻歴 ---
    fig,ax=plt.subplots(figsize=(10,5))
    ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,rmse(da_hi),"-o",ms=4,color="red",label=f"高感度位置の温度センサ (最終 {rmse(da_hi)[-1]:.3f} K)")
    ax.plot(ts,rmse(da_lo),"-s",ms=4,color="deepskyblue",label=f"低感度位置の温度センサ (最終 {rmse(da_lo)[-1]:.3f} K)")
    ax.set_yscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("温度分布RMSE [K]")
    ax.set_title("温度センサ位置(感度マップ最大/最小)による分布再現性 — 変位2点は常に同化")
    ax.grid(alpha=0.3,which="both"); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"evidence_rmse.png"),dpi=140); plt.close(fig)

    qoi=lambda A: (Mum@(A-T_AIR_K).T)[0]-(Mum@(A-T_AIR_K).T)[1]   # 傾き[µm]
    fig,ax=plt.subplots(figsize=(10,5))
    ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,qoi(tr),"k--",lw=2.5,label="真値")
    ax.plot(ts,qoi(da_hi),"-o",ms=4,color="red",label="高感度位置センサ")
    ax.plot(ts,qoi(da_lo),"-s",ms=4,color="deepskyblue",label="低感度位置センサ")
    ax.set_xlabel("time [s]"); ax.set_ylabel("上面の傾き uz_heater−uz_opp [µm]")
    ax.set_title("変位QoI(傾き)の時刻歴 — センサ位置による追従の差")
    ax.grid(alpha=0.3); ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"evidence_qoi.png"),dpi=140); plt.close(fig)
    print("[ev] wrote evidence_rmse.png / evidence_qoi.png")

if __name__=="__main__": main()
