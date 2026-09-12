"""docs/14用の追加解説図: SVD分解の可視化と、条件変化に強い実用計算フロー.

出力(docs/img/):
  svd_decomposition.png  … X ≈ σ1·(型1)·(時間1) + σ2·(型2)·(時間2) を絵で示す
  practical_workflow.png … 一度だけ計算(条件によらない) vs 条件ごと(安い) の実用フロー
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import pod_analysis as pa
import cylinder_mesh
IMG=os.path.join(ROOT,"docs","img"); K=273.15


def fig_svd():
    ts,C,X=pa.load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
    energy=S**2/(S**2).sum()
    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    from scipy.spatial import cKDTree
    _,near=cKDTree(C).query(coords)
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
    import tempfile; TMPD=tempfile.mkdtemp()
    def cyl(field,cmap,fn,sym=True):
        g=ug.copy(); g.point_data["v"]=field[near]
        pl=pv.Plotter(off_screen=True,window_size=(430,470))
        if sym: v=float(np.abs(field).max()); clim=[-v,v]
        else: clim=[float(field.min()),float(field.max())]
        pl.add_mesh(g,scalars="v",cmap=cmap,clim=clim,n_colors=12,show_scalar_bar=False)
        pl.camera_position=CAM; pl.set_background("white")
        p=os.path.join(TMPD,fn); pl.screenshot(p); pl.close(); return p
    pmean=cyl(mean-K,"turbo","mean.png",sym=False)
    pm1=cyl(U[:,0],"coolwarm","m1.png"); pm2=cyl(U[:,1],"coolwarm","m2.png")

    fig=plt.figure(figsize=(15,8.4))
    fig.suptitle("SVD/POD分解:温度スナップショット = 平均場 ＋ Σ(型 × 強さ × 時間変化)",
                 fontsize=17,weight="bold")
    gs=fig.add_gridspec(2,4,width_ratios=[1.1,0.25,1.4,0.05],height_ratios=[1,1],
                        hspace=0.35,wspace=0.15)
    # 平均場
    axm=fig.add_subplot(gs[:,0]); axm.imshow(plt.imread(pmean)); axm.axis("off")
    axm.set_title("平均場\n(いつも共通の下地)",fontsize=13)
    axp=fig.add_subplot(gs[:,1]); axp.axis("off"); axp.text(0.5,0.5,"＋",fontsize=40,ha="center",va="center")
    # モード1行
    for r,(pm,k) in enumerate([(pm1,0),(pm2,1)]):
        sub=gs[r,2].subgridspec(1,3,width_ratios=[1,0.3,1.6],wspace=0.1)
        a1=fig.add_subplot(sub[0]); a1.imshow(plt.imread(pm)); a1.axis("off")
        a1.set_title(f"型{k+1} ({energy[k]*100:.1f}%)\n"+
                     ("全体の昇温" if k==0 else "ヒータ側↔反対側"),fontsize=12)
        ax_x=fig.add_subplot(sub[1]); ax_x.axis("off")
        ax_x.text(0.5,0.5,"×",fontsize=26,ha="center",va="center")
        a2=fig.add_subplot(sub[2])
        a2.plot(ts,S[k]*Vt[k,:],color="tab:red" if k==0 else "tab:blue",lw=2)
        a2.axvspan(0,300,color="orange",alpha=0.06)
        a2.set_title(f"強さ×時間変化 (σ{k+1}·v{k+1})",fontsize=12)
        a2.set_xlabel("time [s]",fontsize=10); a2.grid(alpha=0.3)
        a2.tick_params(labelsize=9)
    fig.text(0.5,0.02,"型(U列)=空間パターン、σ(特異値)=強さ、時間変化(V列)。"
             "上位2型だけで場の99%を説明できる = だからROMは少自由度でよい",
             ha="center",fontsize=13,
             bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
    fig.savefig(os.path.join(IMG,"svd_decomposition.png"),dpi=125,bbox_inches="tight")
    plt.close(fig); print("[svd] svd_decomposition.png")


def fig_workflow():
    fig,ax=plt.subplots(figsize=(14,8.6)); ax.set_xlim(0,14); ax.set_ylim(0,8.6); ax.axis("off")
    ax.text(7,8.3,"条件が変わっても作り直さない実用フロー",ha="center",fontsize=17,weight="bold")
    ax.text(7,7.8,"― 何が「構造だけで決まる(一度きり)」で、何が「条件ごと(安い)」かを分ける ―",
            ha="center",fontsize=12,color="dimgray")

    # 左: 一度だけ(オフライン, 条件によらない)
    ax.add_patch(FancyBboxPatch((0.4,2.4),6.3,4.7,boxstyle="round,pad=0.15",
                                fc="#eef5ec",ec="tab:green",lw=2))
    ax.text(3.55,6.75,"① 一度だけ計算(オフライン)",fontsize=14,ha="center",
            color="tab:green",weight="bold")
    ax.text(3.55,6.35,"形状・材料・拘束だけで決まる → 条件が変わっても不変",
            fontsize=10.5,ha="center",color="dimgray")
    items=[
        ("W = K^-1 H (感度行列)","構造のみ。ヒータ位置・出力に依存しない",5.7),
        ("ROMの回路係数 C, k, h","熱容量・熱抵抗・放熱=材料/形状の物性",4.75),
        ("PODモード(場の型)","運転範囲を張る数ケースのスナップショットから",3.8),
        ("代表点/観測点の選定","POD最適点・W行感度で決める(105,14章)",2.9),
    ]
    for txt,sub,y in items:
        ax.add_patch(FancyBboxPatch((0.7,y-0.32),5.7,0.72,boxstyle="round,pad=0.05",
                                    fc="white",ec="tab:green"))
        ax.text(1.0,y+0.06,txt,fontsize=11.5,weight="bold",va="center")
        ax.text(1.0,y-0.19,sub,fontsize=9,color="dimgray",va="center")

    ax.annotate("",xy=(7.5,4.6),xytext=(6.75,4.6),
                arrowprops=dict(arrowstyle="-|>",color="k",lw=2.5))

    # 右: 条件ごと(オンライン, 安い)
    ax.add_patch(FancyBboxPatch((7.5,2.4),6.1,4.7,boxstyle="round,pad=0.15",
                                fc="#eef2fb",ec="tab:blue",lw=2))
    ax.text(10.55,6.75,"② 条件ごと(オンライン, 安い)",fontsize=14,ha="center",
            color="tab:blue",weight="bold")
    ax.text(10.55,6.35,"ヒータ位置・出力・運転が変わるたびに、ここだけ回す",
            fontsize=10.5,ha="center",color="dimgray")
    items2=[
        ("観測を取得(温度数点+変位)","実機センサ or その条件の1回計算",5.7),
        ("ROM前進予測(数秒)","校正済み回路に、その条件の入力Qを与える",4.75),
        ("データ同化(EnKF)で Q・場を推定","①で決めた観測演算子をそのまま使用",3.8),
        ("必要なら分布/変位を復元","PODモード or FrontISTRで場を再構成",2.9),
    ]
    for txt,sub,y in items2:
        ax.add_patch(FancyBboxPatch((7.8,y-0.32),5.5,0.72,boxstyle="round,pad=0.05",
                                    fc="white",ec="tab:blue"))
        ax.text(8.1,y+0.06,txt,fontsize=11.5,weight="bold",va="center")
        ax.text(8.1,y-0.19,sub,fontsize=9,color="dimgray",va="center")

    ax.text(7,1.55,"条件が「①で想定した運転範囲の中」なら②だけで済む(作り直しゼロ)。",
            ha="center",fontsize=12.5,color="black",
            bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
    ax.text(7,0.75,"範囲を外れた新条件(例:ヒータを全く別の位置へ)だけ、①のPODに"
            "その1ケースを足して再SVD(数秒)。W・回路係数はそれでも不変。",
            ha="center",fontsize=11,color="dimgray")
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"practical_workflow.png"),dpi=135); plt.close(fig)
    print("[svd] practical_workflow.png")


def fig_how_modes():
    """モード(型)はどう出てくるか: 実スナップショットの重み付き和として現れることを示す."""
    ts,C,X=pa.load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
    # 型1の重み: 各時刻スナップショットにかかる係数に比例する量 = 相関行列の固有ベクトル v1
    # (U[:,0] = Xc @ (Vt[0,:]/S[0]) なので、重みは v1[t] に比例)
    w1=Vt[0,:]   # 相関行列 X^T X の第1固有ベクトル(O(0.1)、読みやすい)
    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    from scipy.spatial import cKDTree
    _,near=cKDTree(C).query(coords)
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
    import tempfile; TMPD=tempfile.mkdtemp()
    def cyl(field,cmap,fn,sym=True):
        g=ug.copy(); g.point_data["v"]=field[near]
        pl=pv.Plotter(off_screen=True,window_size=(360,400))
        if sym: v=float(np.abs(field).max()); clim=[-v,v]
        else: clim=[float(field.min()),float(field.max())]
        pl.add_mesh(g,scalars="v",cmap=cmap,clim=clim,n_colors=12,show_scalar_bar=False)
        pl.camera_position=CAM; pl.set_background("white")
        p=os.path.join(TMPD,fn); pl.screenshot(p); pl.close(); return p
    # 代表4時刻のスナップショット(変動成分)と、型1
    tsel=[10,30,60,90]   # index (t=50,150,300,450 s付近)
    snaps=[(cyl(Xc[:,i],"coolwarm",f"s{i}.png"), ts[i], w1[i]) for i in tsel]
    pmode=cyl(U[:,0],"coolwarm","mode1.png")

    fig=plt.figure(figsize=(15,6.6))
    fig.suptitle("「型」はどう出てくるか ― 実スナップショットの重み付き和として現れる",
                 fontsize=16,weight="bold")
    n=len(snaps)
    gs=fig.add_gridspec(1,2*n+2,width_ratios=[1,0.35]*n+[0.6,1.4],wspace=0.05)
    for j,(p,t,w) in enumerate(snaps):
        ax=fig.add_subplot(gs[2*j]); ax.imshow(plt.imread(p)); ax.axis("off")
        ax.set_title(f"t={t:.0f}s の場\n(平均を引いた変動)",fontsize=10.5)
        axo=fig.add_subplot(gs[2*j+1]); axo.axis("off")
        sign="＋" if (j==n-1 or True) else ""
        axo.text(0.5,0.62,f"×{w:+.3f}",fontsize=12,ha="center",va="center",color="crimson")
        axo.text(0.5,0.44,"(に比例)",fontsize=8,ha="center",va="center",color="crimson")
        axo.text(0.5,0.20,"＋" if j<n-1 else "＋…",fontsize=20,ha="center",va="center")
    axe=fig.add_subplot(gs[2*n]); axe.axis("off"); axe.text(0.5,0.5,"=",fontsize=34,ha="center",va="center")
    axm=fig.add_subplot(gs[2*n+1]); axm.imshow(plt.imread(pmode)); axm.axis("off")
    axm.set_title("型1(モード1)\n= 全121枚の重み付き和",fontsize=12,color="tab:red")
    fig.text(0.5,0.03,"重み w_t は「相関行列 X^T X の第1固有ベクトル ÷ 特異値」に比例。"
             "型は勝手に発明されるのでなく、実データの blend として出てくる",
             ha="center",fontsize=12.5,bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
    fig.savefig(os.path.join(IMG,"svd_how_modes.png"),dpi=125,bbox_inches="tight")
    plt.close(fig); print("[svd] svd_how_modes.png")


if __name__=="__main__":
    fig_svd(); fig_how_modes(); fig_workflow()
