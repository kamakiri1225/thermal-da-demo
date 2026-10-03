"""blog_008 §6 用：A・O を測らずに当てる様子を、温度分布と変形の3Dアニメーションで示す.

3本並べ
  左  = 真値（ROM真値から復元した温度場と、FrontISTR由来の変形）
  中  = 温度2点のみで同化（変位を測らない）
  右  = 温度2点 ＋ 選定した変位2点（A・O は一度も測らない）
色＝温度、形＝熱変形（誇張表示）。観測点と評価点 A・O も球で描く。

出力: docs/img/blog008_3d_anim.gif
再現: OMP_NUM_THREADS=4 python3 run/make_blog008_3d_anim.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
KC=273.15; NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02; SEED=20260913
TN=[2,0]
SEL=[((0.0097,-0.0362,0.1005),0),((-0.0346,0.0144,0.1005),2)]
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005)
EXAG=None          # 変形の誇張倍率（最大変位が5 mmになるよう自動設定）


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; D=op["D"]; coords=op["coords"]
    Wall=np.einsum("nck,kj->ncj",D,UPp)
    near_node=lambda x: int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    iA=near_node(A_XYZ); iO=near_node(O_XYZ); sel=[(near_node(x),c) for x,c in SEL]
    field=lambda T5: mean+U@(UPp@(T5-mean[pod]))
    # 全節点の変位（µm）→ m
    dispall=lambda T5: (um+np.einsum("nck,k->nc",Wall,(T5-mean[pod])))*1e-6
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    def run(obs):
        ws=[Wall[i,c] for i,c in obs]; u0=[um[i,c] for i,c in obs]
        rng=np.random.default_rng(SEED); ro=np.random.default_rng(SEED+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*len(TN)+[SIG_U**2]*len(ws)); rec=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        mp=mean[pod]
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][TN]); Yf=Z[:,TN]
            for w,b0 in zip(ws,u0):
                yv.append(b0+w@(Ttr[ci]-mp)); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
            y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            rec.append(Z[:,:NPT].mean(0).copy())
        return np.array(rec)
    tOnly=run([]); tDisp=run(sel)
    # 誇張倍率：全フレームの最大変位が 5 mm に見えるよう自動設定
    global EXAG
    mx=max(np.abs(dispall(x)).max() for arr in (Ttr,tOnly,tDisp) for x in arr)
    EXAG=0.005/mx
    print(f"[anim] DA done. 最大変位 {mx*1e6:.2f} µm → 誇張倍率 {EXAG:.0f}倍",flush=True)

    import pyvista as pv, vtk
    from PIL import Image
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    co=np.array([c for _n,c in mesh["nodes"]]); idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cl=[]
    for _e,conn in mesh["elements"]: cl.append(8); cl.extend(idr[n] for n in conn)
    cells=np.array(cl); ctypes=np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8)
    _,near=cKDTree(Cc).query(co)
    clim=[19.5,26.5]
    def shot(T5,title,col,marks):
        f=field(T5)-KC; u=dispall(T5)
        g=pv.UnstructuredGrid(cells,ctypes,co+u*EXAG); g.point_data["T"]=f[near]
        pl=pv.Plotter(off_screen=True,window_size=(620,760))
        pl.add_mesh(g,scalars="T",cmap="turbo",clim=clim,n_colors=18,show_scalar_bar=False)
        for p0,c0,r0,kind in marks:
            q0=p0+(u[near_node(p0)]*EXAG)
            if kind=="temp":      # 温度計：赤い球
                pl.add_mesh(pv.Sphere(radius=r0,center=q0),color=c0)
            elif kind=="disp":    # 変位計：緑の立方体＋測る方向の矢印
                pl.add_mesh(pv.Cube(center=q0,x_length=2*r0,y_length=2*r0,z_length=2*r0),color=c0)
            else:                 # 評価点：黒い球（測らない）
                pl.add_mesh(pv.Sphere(radius=r0,center=q0),color=c0)
        pl.camera_position=[(0.26,-0.24,0.21),(0,0,0.05),(0,0,1)]
        pl.set_background("white"); pl.camera.zoom(1.5)
        im=pl.screenshot(return_img=True); pl.close()
        m=np.any(im[...,:3]<246,axis=-1); ys,xs=np.where(m)
        if ys.size==0: return im          # 全部白なら切り取らずそのまま返す
        return im[max(0,ys.min()-4):ys.max()+4,max(0,xs.min()-4):xs.max()+4]
    marks_obs=[(np.array(A_XYZ),"#111111",0.0040,"eval"),(np.array(O_XYZ),"#111111",0.0040,"eval"),
               (coords[sel[0][0]],"#17c04f",0.0040,"disp"),(coords[sel[1][0]],"#17c04f",0.0040,"disp"),
               (kv["xyz"][2],"red",0.0042,"temp"),(kv["xyz"][0],"red",0.0042,"temp")]
    marks_t=[(np.array(A_XYZ),"#111111",0.0040,"eval"),(np.array(O_XYZ),"#111111",0.0040,"eval"),
             (kv["xyz"][2],"red",0.0042,"temp"),(kv["xyz"][0],"red",0.0042,"temp")]
    frames=[]
    for k in range(len(tg)):
        ut=dispall(Ttr[k]); u1=dispall(tOnly[k]); u2=dispall(tDisp[k])
        aoT=(ut[iA,2]-ut[iO,2])*1e6; ao1=(u1[iA,2]-u1[iO,2])*1e6; ao2=(u2[iA,2]-u2[iO,2])*1e6
        ims=[shot(Ttr[k],"TRUTH","black",marks_t),
             shot(tOnly[k],"TEMP 2 ONLY","#8a1c1c",marks_t),
             shot(tDisp[k],"TEMP 2 + DISP 2 (elsewhere)","#14459c",marks_obs)]
        fig,axes=plt.subplots(1,3,figsize=(13.8,6.4))
        ttls=[f"真値（正解）\n反り A−O = {aoT:5.2f} µm",
              f"観測＝温度2点（●赤）だけ\n反り {ao1:5.2f} µm（誤差 {abs(ao1-aoT):4.2f} µm）",
              f"観測＝温度2点（●赤）＋変位2点（■緑）\n反り {ao2:5.2f} µm（誤差 {abs(ao2-aoT):4.2f} µm）"]
        cols=["#1b2430","#c0392b","#14459c"]
        for ax,im,t,c in zip(axes,ims,ttls,cols):
            ax.imshow(im); ax.axis("off"); ax.set_title(t,fontsize=12,weight="bold",color=c)
        # 温度のカラーバー
        sm=plt.cm.ScalarMappable(cmap="turbo",norm=plt.Normalize(*clim))
        cax=fig.add_axes([0.25,0.055,0.50,0.022])
        cb=fig.colorbar(sm,cax=cax,orientation="horizontal"); cb.set_label("温度 [℃]",fontsize=10.5)
        fig.suptitle(f"A・O を一度も測らずに熱変形を当てる（t = {tg[k]:.0f} s）\n"
                     f"色＝温度　形＝熱変形（{EXAG:.0f}倍に誇張）　　●赤＝温度センサ（観測）　■緑＝変位計（観測）　●黒＝評価点 A・O（一度も測らない）",fontsize=12)
        fig.tight_layout(rect=[0,0.09,1,0.88])
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:,:,:3].copy()))
        plt.close(fig)
        if k%5==0: print(f"  t={tg[k]:.0f}s  温度のみ誤差 {abs(ao1-aoT):.2f} / 変位あり {abs(ao2-aoT):.2f} µm",flush=True)
    out=os.path.join(IMG,"blog008_3d_anim.gif")
    frames[0].save(out,save_all=True,append_images=frames[1:],duration=260,loop=0,disposal=2)
    print("wrote",out,f"({len(frames)}フレーム)")


if __name__=="__main__": main()
