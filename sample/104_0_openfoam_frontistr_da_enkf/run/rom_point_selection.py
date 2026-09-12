"""ROMの代表点を「流用/POD最適/W感度」で選び分け、全セル温度場の復元精度で比較する.

同じ土俵（102_0の121スナップショットを、5点の温度だけから復元したときのRMSE）で、
代表点の選び方3種を突き合わせる。「流用した測定点」が最適だったのか偶然だったのかを検証。

選び方:
  probes … 現行(102_0の測定プローブ4点+肉厚中心core) = 流用ROM
  qdeim  … PODモードにQ-DEIM(枢軸QR)を適用した最適補間点(データ駆動の理論最適)
  wsens  … 105のW=K⁻¹H行感度が高い上位点(構造の変位応答が大きい点)
  random … 参照(下限)
復元法:
  gappy  … Gappy-POD(5モードの最小二乗当てはめ) = 点の情報量そのものを測る
  idw    … 距離逆二乗補間 = 実際の5点ROMがやっている復元

出力: docs/img/rom_point_selection.png, results/rom_point_selection.yaml
再現: OMP_NUM_THREADS=4 python3 run/rom_point_selection.py
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy.linalg import qr
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore.node_locations import NODE_XYZ
from dacore.cht_rom import NODE_NAMES
import pod_analysis as pa
import cylinder_mesh

IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
K=273.15; NPT=5; RMODE=5


def idw_matrix(C, pts_xyz):
    W=np.zeros((len(C),len(pts_xyz)))
    for i,p in enumerate(C):
        d=np.linalg.norm(pts_xyz-p,axis=1)
        if d.min()<1e-9: W[i,d.argmin()]=1
        else: w=1/d**2; W[i]=w/w.sum()
    return W


def gappy_rmse(U, mean, Xc, X, P):
    """点集合Pの温度から5モード最小二乗で全場復元したRMSE(時刻ごと)."""
    M=U[P,:RMODE]
    errs=[]
    for j in range(X.shape[1]):
        a,_,_,_=np.linalg.lstsq(M,Xc[P,j],rcond=None)
        rec=mean+U[:,:RMODE]@a
        errs.append(np.sqrt(((rec-X[:,j])**2).mean()))
    return np.array(errs)


def idw_rmse(C, X, P):
    W=idw_matrix(C, C[P])
    errs=[np.sqrt(((W@X[P,j]-X[:,j])**2).mean()) for j in range(X.shape[1])]
    return np.array(errs)


def main():
    import yaml
    os.makedirs(IMG,exist_ok=True); os.makedirs(RES,exist_ok=True)
    ts, C, X = pa.load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,Vt=np.linalg.svd(Xc,full_matrices=False)

    # --- 選び方ごとに5つのセル添字を決める ---
    sel={}
    # probes: 現行5点(最寄りセル)
    romn=np.array([NODE_XYZ[n] for n in NODE_NAMES])
    sel["probes"]=[int(np.linalg.norm(C-p,axis=1).argmin()) for p in romn]
    # qdeim: PODモードU[:,:5]に枢軸QR → 最初の5枢軸
    _,_,piv=qr(U[:,:RMODE].T, pivoting=True)
    sel["qdeim"]=list(piv[:NPT])
    # wsens: 105のW行感度 上位(有効節点)を空間分散させて5点→最寄りセル
    kv=np.load(os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity",
                            "results","kinvh_sensitivity.npz"))
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    fem=np.array([xyz for _n,xyz in mesh["nodes"]])
    order=[i for i in np.argsort(kv["row_sens"])[::-1] if kv["valid"][i]]
    picks=[order[0]]
    for i in order[1:]:
        if all(np.linalg.norm(fem[i]-fem[j])>0.02 for j in picks): picks.append(i)
        if len(picks)==NPT: break
    sel["wsens"]=[int(np.linalg.norm(C-fem[i],axis=1).argmin()) for i in picks]
    # random: 3seed平均で評価(点自体はseed0を図示)
    rng=np.random.default_rng(0)
    rand_sets=[list(rng.choice(len(C),NPT,replace=False)) for _ in range(3)]
    sel["random"]=rand_sets[0]

    # --- 評価 ---
    res={}
    span=float(X.max()-X.min())
    for name,P in sel.items():
        g=gappy_rmse(U,mean,Xc,X,P); d=idw_rmse(C,X,P)
        res[name]=dict(gappy_max=float(g.max()),gappy_mean=float(g.mean()),
                       idw_max=float(d.max()),idw_mean=float(d.mean()))
    # randomは3組平均
    gm=[gappy_rmse(U,mean,Xc,X,Pp).max() for Pp in rand_sets]
    dm=[idw_rmse(C,X,Pp).max() for Pp in rand_sets]
    res["random"]["gappy_max"]=float(np.mean(gm)); res["random"]["idw_max"]=float(np.mean(dm))
    for name in sel:
        print(f"[sel] {name:7s} gappy_max={res[name]['gappy_max']:.4f}K "
              f"idw_max={res[name]['idw_max']:.4f}K (span {span:.1f}K)")

    # --- 図1: 棒グラフ(復元RMSE) ---
    labels={"qdeim":"POD最適\n(Q-DEIM)","probes":"流用\n(測定プローブ)",
            "wsens":"W行感度\n上位","random":"ランダム\n(参照)"}
    order_plot=["qdeim","probes","wsens","random"]
    fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5.4))
    for ax,key,ttl in [(a1,"gappy_max","Gappy-POD復元(点の情報量)"),
                       (a2,"idw_max","IDW復元(実ROMの復元法)")]:
        vals=[res[n][key] for n in order_plot]
        bars=ax.bar([labels[n] for n in order_plot],vals,
                    color=["tab:green","tab:orange","tab:red","tab:gray"])
        ax.set_ylabel("全セル温度場の最大RMSE [K]"); ax.set_title(ttl)
        ax.grid(alpha=0.3,axis="y")
        for b,v in zip(bars,vals): ax.text(b.get_x()+b.get_width()/2,v,f"{v:.3f}",
                                           ha="center",va="bottom",fontsize=11)
    fig.suptitle(f"ROM代表点5個の選び方 vs 温度分布の復元精度(102_0, 温度スパン{span:.1f}K)",
                 fontsize=15)
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(IMG,"rom_point_selection.png"),dpi=140); plt.close(fig)

    # --- 図2: 各手法の点配置(3D) ---
    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
    import tempfile; TMPD=tempfile.mkdtemp(); shots=[]
    for name,col in [("qdeim","tab:green"),("probes","orange"),("wsens","red")]:
        cc={"tab:green":"green","orange":"orange","red":"red"}[col]
        pl=pv.Plotter(off_screen=True,window_size=(560,640))
        pl.add_mesh(ug,color="lightsteelblue",opacity=0.4)
        for P in sel[name]:
            pl.add_mesh(pv.Sphere(radius=0.005,center=C[P]),color=cc)
        pl.camera_position=CAM; pl.set_background("white")
        p=os.path.join(TMPD,f"{name}.png"); pl.screenshot(p); pl.close(); shots.append((p,name))
    fig,axes=plt.subplots(1,3,figsize=(14,5.4))
    tt={"qdeim":"POD最適(Q-DEIM):場を覆うよう分散",
        "probes":"流用(測定プローブ):広く配置",
        "wsens":"W行感度上位:上端に集中(変位用の点)"}
    for ax,(p,name) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(tt[name],fontsize=13)
    fig.suptitle("5点をどこに置いたか（球=代表点）",fontsize=15)
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(IMG,"rom_point_placement.png"),dpi=130); plt.close(fig)

    yaml.safe_dump({"span_K":span,"n_points":NPT,"n_modes":RMODE,"results":res},
                   open(os.path.join(RES,"rom_point_selection.yaml"),"w"),allow_unicode=True)
    print("[sel] wrote 2 figs + results/rom_point_selection.yaml")


if __name__=="__main__": main()
