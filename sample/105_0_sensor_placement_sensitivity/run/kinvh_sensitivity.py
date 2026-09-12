"""KinvH流の熱感度(W=K⁻¹H)を唯一のソースとして全感度成果物を再生成する.

W は DUMPW実行がダンプした K(dump_matrix_1_0.mm) と H(H_matrix.mtx) から構築済み
(results/Wz_full.npy, /tmp/build_W.py)。DUMPWのWdiff_zと最大相対差1.9e-7で一致検証済み。

生成物:
  results/kinvh_sensitivity.npz  … row_sens(変位点選定用)/colA/colO(温度影響)/マスク
  docs/img/sensitivity_frontistr.png … W行ノルム + W列(∂uz(A)/∂T) の面表示
  docs/img/disp_point_sensitivity_3d.png … 実験2の選定点(各点ラベル付き)
実験2の高/低感度点もここで W 行感度から確定し、選定結果を npz に保存する。
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, vtk
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
H=0.1005

def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    Wz=np.load(os.path.join(RES,"Wz_full.npy")).astype(np.float64)  # (5040,5040) m/K
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    qa=int(np.linalg.norm(coords-np.array([0.028,0,H]),axis=1).argmin())
    qo=int(np.linalg.norm(coords-np.array([-0.028,0,H]),axis=1).argmin())

    row=np.abs(Wz).sum(axis=1)*1e6                 # 行感度 [µm/K]
    med=np.median(row)
    artifact=(row>100*med)                          # 拘束DOFの見かけ感度(KinvHの教訓)
    valid=(~artifact)&(coords[:,2]>0.005)           # 底面近傍も除外
    print(f"[W] row-sens 有効{valid.sum()}節点 / アーティファクト除外{artifact.sum()}節点")
    print(f"[W] 有効範囲: {row[valid].min():.3f}..{row[valid].max():.3f} µm/K")
    colA=Wz[qa]*1e6; colO=Wz[qo]*1e6                # 列感度 ∂uz(A or O)/∂T(x) [µm/K]

    # 実験2の選定(方向分散つき)
    order=[i for i in np.argsort(row)[::-1] if valid[i]]
    hi=[order[0]]; v0=Wz[order[0]]/np.linalg.norm(Wz[order[0]])
    for i in order[1:]:
        v=Wz[i]/max(np.linalg.norm(Wz[i]),1e-30)
        if abs(v@v0)<0.9: hi.append(i); break
    lo=[i for i in order[::-1][:len(order)] if valid[i]][:0]
    lo=[i for i in np.argsort(row) if valid[i]][:2]
    cur=[qa,qo]
    np.savez(os.path.join(RES,"kinvh_sensitivity.npz"),
             row_sens=row,colA=colA,colO=colO,valid=valid,
             hi=np.array(hi),lo=np.array(lo),cur=np.array(cur),qa=qa,qo=qo,
             note="全てFrontISTR(DUMPW)ダンプのK,Hから構築したW=K^-1Hに基づく感度")
    for tag,pts in [("高感度2点",hi),("低感度2点",lo),("現行(A/O)",cur)]:
        for i in pts:
            print(f"[W] {tag}: node{i} xyz={np.round(coords[i]*1000,1)}mm row_sens={row[i]:.3f}µm/K")

    # メッシュ(面表示)
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
    import tempfile
    TMPD=tempfile.mkdtemp()

    def mark_AO(pl):
        """Point A(黄)/Point O(緑)の球とラベルを描く."""
        pl.add_mesh(pv.Sphere(radius=0.0045,center=coords[qa]),color="yellow")
        pl.add_point_labels([coords[qa]+np.array([0,0,0.007])],["Point A (+X, heater side)"],
                            font_size=19,text_color="black",shape=None,always_visible=True)
        pl.add_mesh(pv.Sphere(radius=0.0045,center=coords[qo]),color="limegreen")
        pl.add_point_labels([coords[qo]+np.array([0,0,0.032])],["Point O (-X, opposite)"],
                            font_size=19,text_color="darkgreen",shape=None,always_visible=True)

    # --- 図1(表・裏の2アングル): W行ノルム + W_diff の2パネル合成 ---
    from dacore import plots as _plots  # 日本語フォント登録
    import matplotlib.pyplot as plt
    g1=ug.copy(); g1.point_data["row"]=np.where(valid,row,np.nan)
    wdiff=colA-colO
    g2=ug.copy(); v=float(np.percentile(np.abs(wdiff[valid]),90)); g2.point_data["col"]=wdiff
    CAM_BACK=[(-0.24,0.22,0.20),(0,0,0.05),(0,0,1)]   # 反対側(−X)から
    def shot(g,scalars,cmap,clim,nc,bar,cam,fname):
        pl=pv.Plotter(off_screen=True,window_size=(780,780))
        pl.add_mesh(g,scalars=scalars,cmap=cmap,clim=clim,n_colors=nc,nan_color="gray",
                    scalar_bar_args={"title":bar,"title_font_size":20,"label_font_size":16})
        mark_AO(pl)
        pl.add_axes(line_width=4,labels_off=False,xlabel="X",ylabel="Y",zlabel="Z")
        pl.camera_position=cam; pl.set_background("white")
        p=os.path.join(TMPD,fname); pl.screenshot(p); pl.close(); return p
    clim_row=[0,float(row[valid].max())]
    p11=shot(g1,"row","turbo",clim_row,8,"row-sens [um/K]",CAM,"p11.png")
    p12=shot(g2,"col","coolwarm",[-v,v],11,"dQoI/dT [um/K]",CAM,"p12.png")
    p21=shot(g1,"row","turbo",clim_row,8,"row-sens [um/K]",CAM_BACK,"p21.png")
    p22=shot(g2,"col","coolwarm",[-v,v],11,"dQoI/dT [um/K]",CAM_BACK,"p22.png")
    fig,axes=plt.subplots(2,2,figsize=(15,15.6))
    col_titles=["Wの「行」ノルム: どこの変位を測ると情報が多いか\n（灰色=候補外。固定部のすぐ近くは計算の都合で非物理な巨大値が出るため）",
                "Wの「列」の差: QoI=Uz(A)−Uz(O) はどこの温度に敏感か\n（=DUMPWが出力するW_diff。赤=温めるとQoI+、青=QoI−）"]
    row_titles=["【表】Point A側(+X)から見る","【裏】Point O側(−X)から見る"]
    for r,(pa,pb) in enumerate([(p11,p12),(p21,p22)]):
        for c,img in enumerate([pa,pb]):
            ax=axes[r][c]; ax.imshow(plt.imread(img)); ax.axis("off")
            if r==0: ax.set_title(col_titles[c],fontsize=15)
        axes[r][0].text(-0.04,0.5,row_titles[r],transform=axes[r][0].transAxes,
                        fontsize=16,weight="bold",rotation=90,va="center",ha="center")
    fig.suptitle("FrontISTR(KinvH)の熱感度行列 $W=K^{-1}H$ — 上段=+X側から/下段=−X側から（各パネル左下に座標軸）。"
                 "黄球=Point A、緑球=Point O",fontsize=14.5)
    fig.tight_layout(rect=[0.02,0,1,0.96])
    fig.savefig(os.path.join(IMG,"sensitivity_frontistr.png"),dpi=120); plt.close(fig)

    # --- 図2: 実験2の選定点(全点ラベル、離散バンド) ---
    pl=pv.Plotter(off_screen=True,window_size=(1050,880))
    pl.add_mesh(g1,scalars="row",cmap="turbo",clim=[0,float(row[valid].max())],
                n_colors=8,nan_color="gray",
                scalar_bar_args={"title":"W row-sens [um/K]","title_font_size":18,"label_font_size":14})
    labels={"hi":("red","HIGH"),"lo":("deepskyblue","LOW"),"cur":("orange","A/O")}
    # ラベル重なり回避の縦オフセット(点ごとに手動指定)
    DZ={("hi",0):0.005,("hi",1):0.018,("cur",0):-0.006,("cur",1):0.024,
        ("lo",0):-0.008,("lo",1):-0.016}
    for key,pts in [("hi",hi),("lo",lo),("cur",cur)]:
        col,tag=labels[key]
        for n,i in enumerate(pts):
            pl.add_mesh(pv.Sphere(radius=0.0045,center=coords[i]),color=col)
            xyz=np.round(coords[i]*1000,0).astype(int)
            lp=coords[i]+np.array([0,0,DZ[(key,n)]])
            pl.add_point_labels([lp],[f"{tag}{n+1} ({xyz[0]},{xyz[1]},{xyz[2]})mm"],
                                font_size=19,text_color=col,shape=None,always_visible=True)
    pl.camera_position=CAM; pl.set_background("white")
    p3=os.path.join(TMPD,"p3.png"); pl.screenshot(p3); pl.close()
    fig,ax=plt.subplots(figsize=(10.2,9.4))
    ax.imshow(plt.imread(p3)); ax.axis("off")
    ax.set_title("変位観測点の候補（色=FrontISTR $W=K^{-1}H$ の行感度）\n"
                 "赤=高感度2点(HIGH1,2)  青=低感度2点(LOW1,2)  橙=104の現行観測点(A/O1,2)",
                 fontsize=15)
    fig.tight_layout()
    fig.savefig(os.path.join(IMG,"disp_point_sensitivity_3d.png"),dpi=130); plt.close(fig)
    print("[W] wrote sensitivity_frontistr.png / disp_point_sensitivity_3d.png")

if __name__=="__main__": main()
