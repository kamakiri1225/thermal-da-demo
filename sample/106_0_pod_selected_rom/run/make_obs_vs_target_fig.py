"""「観測に使う点」と「本当に知りたい点」を分けた実験の構図を図にする.

循環を避けるための設計:
  観測に使う  : 温度2点（P2, P0）＋ 変位2点（上面 A/O, z=100.5mm）
  本当に知りたい: 変位2点（中高さ B/C, z=75.0mm）← 同化に一切使わない
A/O を同化に使って、B/C が正しくなるかを見る。

出力: docs/img/obs_vs_target.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_obs_vs_target_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
H=0.1005


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    Cc=kv["cell_centres"]; pod=kv["cell_idx"]
    d=json.load(open(os.path.join(RES,"xyz_disp_da.json")))
    rows={(r["point"],r["comp"]):r for r in d["rows"]}

    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    co=np.array([p for _n,p in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,cc in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in cc)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),co)

    OBS=[("A",( 0.0287,0,0.1005)),("O",(-0.0287,0,0.1005))]
    TGT=[("B",( 0.0287,0,0.0754)),("C",(-0.0287,0,0.0754))]
    TEMP=[2,0]                       # 観測に使う温度点（発熱感度 上位2）

    pl=pv.Plotter(off_screen=True, window_size=(760,860))
    pl.add_mesh(ug,color="#c3cedb",opacity=0.22)
    for i in TEMP:
        x=Cc[pod[i]]
        pl.add_mesh(pv.Sphere(radius=0.0034,center=x),color="#1F4E9C")
        pl.add_point_labels([x+np.array([0.004,0,0.004])],[f"P{i} temp OBS"],
            font_size=25,text_color="#1F4E9C",shape=None,always_visible=True,bold=True)
    for lab,x in OBS:
        x=np.array(x)
        pl.add_mesh(pv.Sphere(radius=0.0042,center=x),color="#C0392B")
        pl.add_mesh(pv.Arrow(start=x,direction=(0,0,1),scale=0.030),color="#C0392B")
        pl.add_point_labels([x+np.array([0,0,0.015])],[f"{lab} disp OBS"],
            font_size=30,text_color="#C0392B",shape=None,always_visible=True,bold=True)
    for lab,x in TGT:
        x=np.array(x)
        pl.add_mesh(pv.Sphere(radius=0.0042,center=x),color="#2e7d32")
        pl.add_point_labels([x+np.array([0,0,-0.011])],[f"{lab} EVAL"],
            font_size=31,text_color="#2e7d32",shape=None,always_visible=True,bold=True)
    pl.camera_position=[(0.30,-0.27,0.19),(0,0,0.055),(0,0,1)]
    pl.set_background("white"); pl.camera.zoom(1.16)
    img=pl.screenshot(return_img=True); pl.close()
    m=np.any(img[:,:,:3]<246,axis=-1); ys,xs=np.where(m); p=6
    img=img[max(0,ys.min()-p):ys.max()+p, max(0,xs.min()-p):xs.max()+p]

    fig=plt.figure(figsize=(16.6,8.4))
    gs=fig.add_gridspec(1,2,width_ratios=[0.92,1.18],wspace=.16)
    ax=fig.add_subplot(gs[0]); ax.imshow(img); ax.axis("off")
    ax.set_title("観測点と評価点（未観測点）を分ける\n"
                 "青・赤＝観測点（OBS）／緑＝評価点＝未観測点（EVAL）\n温度P2・P0と上面A/Oを観測し、中高さB/Cは同化に一切使わない",
                 fontsize=13,weight="bold")

    ax=fig.add_subplot(gs[1]); ax.axis("off")
    ax.set_title("手続きと結果", fontsize=14, weight="bold")
    steps=["① 観測点を決める\n　 温度 P2・P0 ＋ 変位 A・O（上面 z=100.5）",
           "② 評価点（未観測点）を決める\n　 変位 B・C（中高さ z=75.4）← 同化に一切使わない",
           "③ でたらめな初期温度・誤ったQ,h からアンサンブル60本",
           "④ 30秒ごとに前進 → ①の4つだけ観測して状態を更新（20回）",
           "⑤ 同化後の温度場から 評価点 B・C の変位を計算する",
           "⑥ 評価点の真値と比べる ← ここが本当の成績"]
    for i,t in enumerate(steps):
        y=0.93-i*0.108
        c="#C0392B" if i in (1,5) else "#1F4E9C"
        ax.add_patch(plt.Rectangle((0.01,y-0.072),0.97,0.095,transform=ax.transAxes,
                     fc=("#fdecea" if i in (1,5) else "#eef3f9"),ec=c,lw=1.7,clip_on=False))
        ax.text(0.035,y-0.024,t,transform=ax.transAxes,fontsize=11.6,va="center")

    tb=[("B  $U_z$",rows[("B","Uz")]),("C  $U_z$",rows[("C","Uz")]),
        ("B  $U_x$",rows[("B","Ux")]),("C  $U_x$",rows[("C","Ux")])]
    ax.text(0.02,0.265,"⑥の結果（加熱期60–300 s のMAE、5 seed平均）",
            transform=ax.transAxes,fontsize=12.5,weight="bold",color="#C0392B")
    ax.text(0.03,0.198,"評価点の量　　　　同化なし　　　同化後　　　改善",
            transform=ax.transAxes,fontsize=11.5)
    for i,(lab,r) in enumerate(tb):
        ax.text(0.03,0.155-i*0.040,
                f"{lab:<12}{r['mae_noDA_after1cyc_um']:>9.3f} µm{r['mae_DA_after1cyc_um']:>9.3f} µm"
                f"{r[chr(39)+chr(39)] if False else r['improve_after1cyc']:>6.0f} x",
                transform=ax.transAxes,fontsize=11.5,family="monospace",weight="bold")
    ax.text(0.02,-0.055,"※ $U_x$（半径方向）は x 方向の観測を一切していないのに直る。\n"
            "　 温度場が正しくなれば $u=W(T-T_{\\mathrm{ref}})$ で全方向が従うため。",
            transform=ax.transAxes,fontsize=11,color="#556")

    fig.suptitle("観測点 A/O から、評価点（未観測点）B/C を当てられるか ― 循環のない検証",
                 fontsize=16,weight="bold",y=0.985)
    fig.tight_layout(rect=[0,0.055,1,0.885])
    o=os.path.join(IMG,"obs_vs_target.png"); fig.savefig(o,dpi=130); plt.close(fig)
    print("wrote",o)


if __name__=="__main__":
    main()
