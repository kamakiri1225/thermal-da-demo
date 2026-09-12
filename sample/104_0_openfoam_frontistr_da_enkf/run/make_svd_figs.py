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
    fig,ax=plt.subplots(figsize=(13.5,10)); ax.set_xlim(0,13.5); ax.set_ylim(0,10); ax.axis("off")
    ax.text(6.75,9.6,"実用の計算フロー ― 「1回だけの準備」と「毎回の本番」に分ける",
            ha="center",fontsize=17,weight="bold")

    # ===== フェーズ1: 準備(装置を動かす前に、計算機で1回だけ) =====
    ax.add_patch(FancyBboxPatch((0.4,6.0),12.7,3.0,boxstyle="round,pad=0.1",
                                fc="#eef5ec",ec="tab:green",lw=2.5))
    ax.text(0.75,8.55,"STEP 1｜準備（装置を動かす“前”に、計算機で1回だけ）",fontsize=14.5,
            color="tab:green",weight="bold")
    ax.text(0.75,8.15,"下の3つは別々のデータから独立に作る(順番ではない)。それを使って観測点を決める。"
            "形状・材料由来なので条件が変わっても作り直さない",
            fontsize=10.0,color="dimgray")
    # 独立な3つの準備(横並び・矢印なし)。各箱に「どのデータ源から」を明記
    prep=[("W = K^-1 H","温度→変形の感度行列","← 構造FEM (剛性K,H)"),
          ("ROM係数 C,k,h","熱容量・熱抵抗・放熱","← 熱解析に校正"),
          ("POD (任意)","何自由度要るか等の診断","← 省略可。使うなら熱ｽﾅｯﾌﾟｼｮｯﾄ")]
    centers=[]
    x=0.9
    for i,(t,s,src) in enumerate(prep):
        ax.add_patch(FancyBboxPatch((x,6.9),3.0,1.05,boxstyle="round,pad=0.06",fc="white",ec="tab:green"))
        ax.text(x+1.5,7.68,t,fontsize=12,weight="bold",ha="center",va="center")
        ax.text(x+1.5,7.34,s,fontsize=9.2,color="dimgray",ha="center",va="center")
        ax.text(x+1.5,7.05,src,fontsize=8.8,color="tab:green",ha="center",va="center",style="italic")
        centers.append(x+1.5)
        x+=3.25
    ax.text(11.9,7.4,"※互いに独立\n(Wから\nPODは出せない)\nPODは診断用で\n無くてもROMは\n作れる",fontsize=8.4,
            color="crimson",ha="center",va="center")
    # 3つ → 観測点を決める箱へ集約 (各箱の真下から短く降ろし、水平線で束ねる)
    ybus=6.6
    for xc in centers:
        ax.annotate("",xy=(xc,ybus),xytext=(xc,6.88),
                    arrowprops=dict(arrowstyle="-",color="tab:green",lw=1.4))
    ax.plot([centers[0],centers[-1]],[ybus,ybus],color="tab:green",lw=1.4)
    ax.add_patch(FancyBboxPatch((3.5,5.95),6.0,0.5,boxstyle="round,pad=0.05",fc="#f6fbf4",ec="tab:green"))
    ax.annotate("",xy=(6.5,6.47),xytext=(6.5,ybus),
                arrowprops=dict(arrowstyle="-|>",color="tab:green",lw=1.6))
    ax.text(6.5,6.2,"→ 観測点を決める（温度＝POD最適点／変位＝W行感度）",fontsize=10.5,
            weight="bold",ha="center",va="center",color="tab:green")

    ax.annotate("",xy=(6.75,5.95),xytext=(6.75,5.35),
                arrowprops=dict(arrowstyle="-|>",color="k",lw=3))
    ax.text(6.95,5.62,"準備完了。あとは本番で使い回すだけ",fontsize=11,va="center",style="italic")

    # ===== フェーズ2: 本番(運転しながら、何度も繰り返す) =====
    ax.add_patch(FancyBboxPatch((0.4,1.25),12.7,4.05,boxstyle="round,pad=0.1",
                                fc="#eef2fb",ec="tab:blue",lw=2.5))
    ax.text(0.75,4.85,"STEP 2｜本番（装置を動かし“ながら”、一定間隔で繰り返す）",fontsize=14.5,
            color="tab:blue",weight="bold")
    ax.text(0.75,4.45,"温度を数点測るだけで、全体の温度分布と熱変位・発熱量を毎サイクル推定する（1周＝数秒）",
            fontsize=10.3,color="dimgray")
    loop=[("① 測る","温度2点＋変位\n(実機センサ)",1.5),
          ("② ROMで予測","校正済みROMを\n数秒回す",4.2),
          ("③ データ同化","測定で予測を補正\nQ・全体温度を推定",6.9),
          ("④ 使う","分布/熱変位を出力\n→補正に反映",9.6)]
    for i,(t,s,x) in enumerate(loop):
        ax.add_patch(FancyBboxPatch((x,2.2),2.5,1.5,boxstyle="round,pad=0.06",fc="white",ec="tab:blue"))
        ax.text(x+1.25,3.35,t,fontsize=12.5,weight="bold",ha="center",va="center")
        ax.text(x+1.25,2.65,s,fontsize=9.3,color="dimgray",ha="center",va="center")
        if i<3: ax.annotate("",xy=(x+2.8,2.95),xytext=(x+2.52,2.95),
                            arrowprops=dict(arrowstyle="-|>",color="tab:blue",lw=2.2))
    # ループの戻り矢印(④→①)
    ax.annotate("",xy=(2.6,2.15),xytext=(11.0,2.15),
                arrowprops=dict(arrowstyle="-|>",color="tab:blue",lw=2,
                                connectionstyle="arc3,rad=-0.06",ls="--"))
    ax.text(6.75,1.62,"④まで済んだら、次のタイミングでまた①へ（何度でも繰り返す＝リアルタイム推定）",
            fontsize=10.5,ha="center",color="tab:blue",style="italic")

    ax.text(6.75,0.62,"条件（ヒータ位置・出力・運転）が変わっても、変えるのはSTEP2の入力だけ。"
            "STEP1は作り直さない。",ha="center",fontsize=12,
            bbox=dict(boxstyle="round",fc="#fffbe6",ec="orange"))
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
