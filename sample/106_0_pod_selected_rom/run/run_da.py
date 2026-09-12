"""POD選定5点ROMで拡大状態EnKFデータ同化（双子実験）＋gappy-POD全場復元.

拡大状態 Z = [T0..T4, q_scale, h]（7次元）。
少数の温度観測から、5点温度・発熱量Q・放熱hを同時推定し、
gappy-PODで全温度場を復元する。温度センサ位置は dT/dQ 最大のノードに置く。

出力(docs/img/):
  da_rmse.png     … 5点温度RMSEの収束（でたらめ初期→真値）
  da_qh.png       … Q と h の推定収束
  da_field.png    … gappy-POD復元の全場（初期/同化後/真値）
出力: results/da_summary.yaml
再現: OMP_NUM_THREADS=4 python3 run/run_da.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh, vtk
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
K=273.15; NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; INFL=1.02
SEED=20260913


def load_rom():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    return (d["cells"], d["xyz"], int(d["heat_node"]), d["C"],
            rg.tri_to_matrix(d["K_upper"],NPT), float(d["h"]))


def integrate(Z, C, Kmat, heat_node, t0, t1):
    """アンサンブルZ(n,7)を t0->t1 前進（T,q,hを使う）."""
    Tn=rg.integrate_ensemble(Z[:,:NPT], C, Kmat, Z[:,IH], Z[:,IQ], heat_node, t0, t1, DT)
    Z2=Z.copy(); Z2[:,:NPT]=Tn; return Z2


def dTdQ(C, Kmat, heat_node):
    """各ノード温度の発熱Qへの感度 dT_i/dq_scale @300s（有限差分）."""
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,0.0154,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,0.0154,1.1,heat_node,0,300,DT)[1][-1]
    return (pert-base)/0.1


def main():
    import yaml
    cells,xyz,heat_node,C,Kmat,h_true=load_rom()
    rng=np.random.default_rng(SEED); rng_o=np.random.default_rng(SEED+7)

    # --- 温度センサ選定: dT/dQ 最大のノード ---
    sens=dTdQ(C,Kmat,heat_node)
    obs_node=int(np.argmax(sens))
    print(f"[da] dT/dQ = {np.round(sens,3)} -> 温度センサは P{obs_node}(最大感度)")

    # --- 真値(twin): 校正ROMを q=1,h=h_true で走らせる ---
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    ttr=[0.0]; Ttr=[np.full(NPT,rg.T_AIR_K)]
    T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]
        ttr.append(b); Ttr.append(T.copy())
    Ttr=np.array(Ttr)   # (ncyc+1,5)

    # --- アンサンブル初期化(でたらめ) ---
    Z=np.zeros((N_ENS,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3, rg.T_AIR_K+12, (N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS)          # q_scale 事前
    Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)  # h 事前
    R=np.array([[SIG_T**2]])

    rec_rmse=[np.sqrt(((Z[:,:NPT].mean(0)-Ttr[0])**2).mean())]
    rec_q=[(Z[:,IQ].mean(),Z[:,IQ].std())]; rec_h=[(Z[:,IH].mean(),Z[:,IH].std())]
    rec_T=[Z[:,:NPT].mean(0).copy()]; rec_Ts=[Z[:,:NPT].std(0).copy()]   # 各ノード推定の平均・ばらつき
    tp=0.0
    for ci,tb in enumerate(cyc,1):
        Z=integrate(Z,C,Kmat,heat_node,tp,tb); tp=tb
        y=np.array([Ttr[ci][obs_node]+rng_o.normal(0,SIG_T)])
        Yf=Z[:,obs_node:obs_node+1]
        Z=enkf_update(Z,y,None,R,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0.0,3.0); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        rec_rmse.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean()))
        rec_q.append((Z[:,IQ].mean(),Z[:,IQ].std())); rec_h.append((Z[:,IH].mean(),Z[:,IH].std()))
        rec_T.append(Z[:,:NPT].mean(0).copy()); rec_Ts.append(Z[:,:NPT].std(0).copy())
    rec_T=np.array(rec_T); rec_Ts=np.array(rec_Ts)
    rec_q=np.array(rec_q); rec_h=np.array(rec_h)
    tgrid=np.r_[0,cyc]
    print(f"[da] 5点温度RMSE: {rec_rmse[0]:.3f}K -> {rec_rmse[-1]:.4f}K")
    print(f"[da] Q推定: 最終 q_scale={rec_q[-1,0]:.3f} (真値1.0) = {rec_q[-1,0]*15:.2f}W")
    print(f"[da] h推定: 最終 h={rec_h[-1,0]:.4f} (真値{h_true:.4f})")

    # --- 図-1: 変位QoIの追従（基本設定: 温度1点だけの同化でも変位が真値に追従） ---
    dop=os.path.join(RES,"disp_operator.npz")
    if os.path.exists(dop):
        do=np.load(dop); uz_mean=do["uz_mean"]; Dmode=do["Dmode"]
        kvq=np.load(os.path.join(RES,"qdeim_points.npz"))
        Uq=kvq["pod_modes"].astype(float); meanq=kvq["mean"].astype(float); pc=kvq["cell_idx"]
        UPp=np.linalg.pinv(Uq[pc,:])
        def qoi(T5):
            a=(T5-meanq[pc])@UPp.T; u=uz_mean+a@Dmode.T; return u[...,0]-u[...,1]
        qoi_tr=np.array([qoi(Ttr[k]) for k in range(len(tgrid))])
        qoi_da=np.array([qoi(rec_T[k]) for k in range(len(tgrid))])
        fig,ax=plt.subplots(figsize=(10,5)); ax.axvspan(0,300,color="orange",alpha=0.06)
        ax.plot(tgrid,qoi_tr,"-",color="k",lw=4,alpha=0.35,label="真値")
        ax.plot(tgrid,qoi_da,"--",color="tab:purple",lw=2.2,dashes=(4,3),label="同化推定(温度1点のみ)")
        ax.set_xlabel("time [s]"); ax.set_ylabel("変位QoI Uz(A)−Uz(O) [µm]")
        ax.set_title("変位の推定は真値を捉えているか（温度1点の同化→gappy-POD場→FrontISTR）\n"
                     "温度の推定誤差が変位QoIにどう現れるかを確認")
        ax.grid(alpha=0.3); ax.legend()
        fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_disp_basic.png"),dpi=140); plt.close(fig)
        print("[da] wrote da_disp_basic.png")

    # 各点を分け、同化後平均と真値の差を明示する。
    np.savez(os.path.join(RES,"da_history.npz"), time=tgrid, truth_T=Ttr,
             analysis_T=rec_T, analysis_std_T=rec_Ts, q=rec_q, h=rec_h,
             obs_node=obs_node, rmse=np.asarray(rec_rmse))
    fig,axes=plt.subplots(3,2,figsize=(13,11),sharex=True)
    for i,ax in enumerate(axes.flat):
        ax.axvspan(0,300,color="orange",alpha=0.08)
        if i<NPT:
            ax.plot(tgrid,Ttr[:,i]-K,color="black",lw=2.5,label="真値（正解）")
            ax.plot(tgrid,rec_T[:,i]-K,"--o",ms=3,color="tab:blue",label="同化推定（60メンバー平均）")
            ax.fill_between(tgrid,rec_T[:,i]-rec_Ts[:,i]-K,
                            rec_T[:,i]+rec_Ts[:,i]-K,color="tab:blue",alpha=0.16,label="メンバーのばらつき±1σ")
            ax.set_title(f"P{i}："+("← 唯一の観測点" if i==obs_node else "未観測（相関で推定）"),
                         fontsize=12,weight=("bold" if i==obs_node else "normal"))
            ax.set_ylabel("温度 [℃]")
            # t=0 のズレ＝でたらめ初期。これが縮むのが見どころ、と1枚だけ注記
            if i==0:
                ax.annotate("わざと外した初期推定\n(t=0で誤差 数K)",xy=(0,rec_T[0,i]-K),
                            xytext=(75,rec_T[0,i]-K+1.3),fontsize=10,color="crimson",
                            arrowprops=dict(arrowstyle="->",color="crimson"))
                ax.annotate("観測を取り込むと\n真値の線に重なる",xy=(300,Ttr[np.argmin(np.abs(tgrid-300)),i]-K),
                            xytext=(340,Ttr[np.argmin(np.abs(tgrid-300)),i]-K-2.6),fontsize=10,color="green",
                            arrowprops=dict(arrowstyle="->",color="green"))
        else:
            for j in range(NPT):
                ax.plot(tgrid,rec_T[:,j]-Ttr[:,j],label=f"P{j}")
            ax.axhline(0,color="black",lw=1)
            ax.set_title("推定の残差（同化推定 − 真値）＝この図の要点",fontsize=12)
            ax.set_ylabel("温度誤差 [K]")
            ax.annotate("初期は数K外れている",xy=(0,4.5),xytext=(110,4.2),fontsize=10,color="dimgray",
                        arrowprops=dict(arrowstyle="->",color="dimgray"))
            ax.annotate("加熱期の終わりには\nほぼ0（小さな残差は残る）",xy=(300,0),xytext=(315,2.4),fontsize=10,color="green",
                        arrowprops=dict(arrowstyle="->",color="green"))
        ax.grid(alpha=0.3); ax.legend(fontsize=9); ax.set_xlabel("時刻 [s]")
    fig.suptitle(f"温度センサ1点（P{obs_node}）だけの観測で、わざと外した初期推定がどれだけ真値に近づくか\n"
                 "黒=真値／青破線=同化推定（60メンバー平均）／水色帯=メンバーのばらつき。"
                 "見どころは「t=0 の大きなズレ→観測で縮む」（要点は右下の残差パネル）",fontsize=13)
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(os.path.join(IMG,"da_track.png"),dpi=140); plt.close(fig)

    # --- 図1: RMSE収束 ---
    fig,ax=plt.subplots(figsize=(9,5)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.semilogy(tgrid,rec_rmse,"-o",ms=4,color="tab:blue")
    ax.set_xlabel("time [s]"); ax.set_ylabel("5点温度の推定RMSE [K]")
    ax.set_title("POD選定5点ROMのデータ同化: でたらめ初期→真値に収束"); ax.grid(alpha=0.3,which="both")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_rmse.png"),dpi=140); plt.close(fig)

    # --- 図2: Q,h 推定 ---
    fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5))
    a1.axvspan(0,300,color="orange",alpha=0.06)
    a1.plot(tgrid,rec_q[:,0]*15,"-o",ms=4,color="tab:red")
    a1.fill_between(tgrid,(rec_q[:,0]-rec_q[:,1])*15,(rec_q[:,0]+rec_q[:,1])*15,color="tab:red",alpha=0.2)
    a1.axhline(15,ls="--",color="k",label="真値 15W"); a1.set_ylabel("発熱量 Q [W]")
    a1.set_xlabel("time [s]"); a1.set_title("発熱量Qの推定"); a1.legend(); a1.grid(alpha=0.3)
    a2.axvspan(0,300,color="orange",alpha=0.06)
    a2.plot(tgrid,rec_h[:,0],"-o",ms=4,color="tab:green")
    a2.fill_between(tgrid,rec_h[:,0]-rec_h[:,1],rec_h[:,0]+rec_h[:,1],color="tab:green",alpha=0.2)
    a2.axhline(h_true,ls="--",color="k",label=f"真値 {h_true:.4f}")
    a2.set_ylabel("放熱係数 h [W/K]"); a2.set_xlabel("time [s]")
    a2.set_title("放熱hの推定（冷却期に決まる）"); a2.legend(); a2.grid(alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_qh.png"),dpi=140); plt.close(fig)

    # --- 図3: gappy-POD 全場復元(初期/同化後/真値) @ 加熱ピーク近傍(t=300) ---
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    pod_cells=kv["cell_idx"]
    ci300=int(np.argmin(np.abs(tgrid-300)))
    def gappy(T5):   # 5点温度 -> 全場(gappy-POD)
        a,_,_,_=np.linalg.lstsq(U[pod_cells,:], T5-mean[pod_cells], rcond=None)
        return mean+U@a
    # 保存した実際の解析平均を再構成する。真値＋乱数による代用は禁止。
    field_truth=gappy(Ttr[ci300])
    field_init=gappy(rec_T[0])
    field_da=gappy(rec_T[ci300])
    np.savez(os.path.join(RES,"da_field_values.npz"), time=tgrid[ci300],
             initial=field_init, analysis=field_da, truth=field_truth)
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz2 for _n,xyz2 in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cellsL=[]
    for _e,conn in mesh["elements"]: cellsL.append(8); cellsL.extend(idr[n] for n in conn)
    from scipy.spatial import cKDTree
    _,near=cKDTree(Cc).query(coords)
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    ug=pv.UnstructuredGrid(np.array(cellsL),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    clim=[float(min(f.min() for f in (field_init,field_da,field_truth))-K),
          float(max(f.max() for f in (field_init,field_da,field_truth))-K)]
    import tempfile; TMP=tempfile.mkdtemp(); shots=[]
    for lab,fld in [("初期推定 t=0s",field_init),("同化後 t=300s",field_da),("ROM真値の復元 t=300s",field_truth)]:
        g=ug.copy(); g.point_data["T"]=fld[near]-K
        pl=pv.Plotter(off_screen=True,window_size=(520,560))
        pl.add_mesh(g,scalars="T",cmap="turbo",clim=clim,n_colors=14,show_scalar_bar=False)
        pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; pl.set_background("white")
        p=os.path.join(TMP,lab+".png"); pl.screenshot(p); pl.close(); shots.append((p,lab))
    fig,axes=plt.subplots(1,3,figsize=(14,5.2))
    for ax,(p,lab) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(lab,fontsize=14)
    sm=plt.cm.ScalarMappable(cmap="turbo",norm=plt.Normalize(*clim))
    fig.colorbar(sm,ax=list(axes),fraction=0.025,pad=0.02,label="温度 [℃]（全パネル共通）")
    fig.suptitle("gappy-POD全場復元：実際の同化後平均を使用（参照もROMから復元）",fontsize=14)
    fig.savefig(os.path.join(IMG,"da_field.png"),dpi=130); plt.close(fig)

    yaml.safe_dump({"n_ens":N_ENS,"obs_node":obs_node,"dTdQ":[float(x) for x in sens],
        "rmse_init_K":float(rec_rmse[0]),"rmse_final_K":float(rec_rmse[-1]),
        "q_final":float(rec_q[-1,0]),"h_final":float(rec_h[-1,0]),"h_true":float(h_true)},
        open(os.path.join(RES,"da_summary.yaml"),"w"),allow_unicode=True)
    print("[da] wrote da_rmse.png, da_qh.png, da_field.png, da_summary.yaml")


if __name__=="__main__": main()
