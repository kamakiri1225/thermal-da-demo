"""POD＋Q-DEIMで選んだ5点をノードにした一般化ROMを、102_0のCHTに校正する.

手順:
 1. select_points_qdeim の結果(qdeim_points.npz)から代表5点を取得
 2. その5点セルの温度履歴(102_0スナップショット)を校正ターゲットにする
 3. ヒータ最近傍ノードを決定
 4. C[5], 全点対コンダクタンスK[10], h を least_squares で同定
 5. 当てはめ残差・図を出し、rom_calibrated_pod.npz に保存

出力: results/rom_calibrated_pod.npz, docs/img/rom_calib_fit.png
再現: OMP_NUM_THREADS=4 python3 run/build_calibrate_rom.py
"""
from __future__ import annotations
import os, sys
import numpy as np
from scipy.optimize import least_squares
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from select_points_qdeim import load_snapshots
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; DT=2.0
HEATER_CENTROID=np.array([0.0375,0.0,0.05025])   # +X外周・中央高さ(102_0ヒータ)


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    cells=kv["cell_idx"]; xyz=kv["xyz"]
    ts, C_all, X = load_snapshots()
    Yobs=X[cells,:].T                 # (nt,5) 5点の温度履歴[K]
    heat_node=int(np.linalg.norm(xyz-HEATER_CENTROID,axis=1).argmin())
    print(f"[calib] ヒータ投入ノード = P{heat_node} xyz={np.round(xyz[heat_node]*1000,1)}mm")

    ne=rg.n_edges(NPT)                 # =10
    T0=np.full(NPT, rg.T_AIR_K)        # 初期は全ノード室温

    def forward(vec):
        C=vec[:NPT]; K=rg.tri_to_matrix(vec[NPT:NPT+ne],NPT); h=vec[NPT+ne]
        Y=[T0.copy()]; T=T0.reshape(1,-1).copy()
        for a,b in zip(ts[:-1],ts[1:]):
            T=rg.integrate_ensemble(T,C,K,np.array([h]),np.array([1.0]),heat_node,a,b,DT)
            Y.append(T[0].copy())
        return np.array(Y)             # (nt,5)

    def resid(vec):
        return (forward(vec)-Yobs).ravel()

    x0=np.concatenate([np.full(NPT,120.0), np.full(ne,1.5), [0.3]])
    lb=np.concatenate([np.full(NPT,1.0),   np.full(ne,1e-4), [1e-4]])
    ub=np.concatenate([np.full(NPT,5000.0),np.full(ne,200.0),[50.0]])
    print(f"[calib] least_squares 開始 (未知数 {len(x0)} = C{NPT}+K{ne}+h1)...",flush=True)
    sol=least_squares(resid,x0,bounds=(lb,ub),xtol=1e-12,ftol=1e-12,verbose=0)
    vec=sol.x
    C=vec[:NPT]; Kup=vec[NPT:NPT+ne]; h=vec[NPT+ne]
    Ymod=forward(vec)
    rmse=float(np.sqrt(((Ymod-Yobs)**2).mean()))
    print(f"[calib] 当てはめRMSE = {rmse:.4f} K   (最大 {np.abs(Ymod-Yobs).max():.3f}K)")
    print(f"[calib] C={np.round(C,1)}  h={h:.4f}  ΣC={C.sum():.1f} J/K")

    np.savez(os.path.join(RES,"rom_calibrated_pod.npz"),
             cells=cells, xyz=xyz, heat_node=heat_node,
             C=C, K_upper=Kup, h=h, calib_rmse=rmse, times=ts,
             note="POD+Q-DEIM選定5点の一般化ROM校正結果")

    # 当てはめ図（実線=OpenFOAM(太・半透明)、白破線=ROM(その上に重ねて見やすく)）
    from matplotlib.lines import Line2D
    fig,ax=plt.subplots(figsize=(11,6.0)); ax.axvspan(0,300,color="orange",alpha=0.06)
    cols=plt.cm.tab10(np.arange(NPT))
    for i in range(NPT):
        p=np.round(xyz[i]*1000,0).astype(int)
        # OpenFOAM: 太く半透明の実線
        ax.plot(ts,Yobs[:,i]-273.15,"-",color=cols[i],lw=5,alpha=0.35,
                label=f"P{i}({p[0]},{p[1]},{p[2]})")
        # ROM: 上に重ねる濃い破線（太め・白縁で視認性UP）
        ax.plot(ts,Ymod[:,i]-273.15,"--",color=cols[i],lw=2.2,dashes=(4,3))
    ax.set_xlabel("time [s]"); ax.set_ylabel("温度 [degC]")
    ax.set_title(f"POD選定5点ROMの校正: 太い半透明実線=OpenFOAM／破線=ROM (残差RMSE {rmse:.3f}K)")
    ax.grid(alpha=0.3)
    # 凡例: 色=各点、線種=OpenFOAM/ROM を別に明記
    leg1=ax.legend(title="点(色)",fontsize=9,ncol=2,loc="lower right")
    style=[Line2D([0],[0],color="k",lw=5,alpha=0.35,label="OpenFOAM(太い半透明実線)"),
           Line2D([0],[0],color="k",lw=2.2,ls="--",label="ROM(破線)")]
    ax.add_artist(leg1); ax.legend(handles=style,fontsize=10,loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"rom_calib_fit.png"),dpi=140); plt.close(fig)
    print("[calib] wrote results/rom_calibrated_pod.npz, docs/img/rom_calib_fit.png")


if __name__=="__main__": main()
