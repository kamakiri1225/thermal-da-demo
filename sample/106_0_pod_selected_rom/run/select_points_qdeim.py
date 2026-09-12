"""POD＋Q-DEIM で、ROMの代表点をOpenFOAM計算から系統的に選ぶ（106の出発点）.

「代表点を勘で置く」のをやめ、元の102_0のCHT温度場スナップショットにPODをかけ、
Q-DEIM（枢軸QR）で「場を最もよく覆う補間点」をデータから決める。

手順:
 1. 102_0 の固体温度場(20,696セル)を0→600s・5秒間隔で121スナップショット読込
 2. 平均を引いてSVD→PODモード U
 3. 先頭 r モードに Q-DEIM(枢軸QR)を適用 → r個のセル(=代表点)を選ぶ
 4. 選点の座標を保存(results/qdeim_points.npz) + 3D図

出力: results/qdeim_points.npz, docs/img/qdeim_points.png
再現: OMP_NUM_THREADS=4 python3 run/select_points_qdeim.py
"""
from __future__ import annotations
import os, re, sys
import numpy as np
from scipy.linalg import qr
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
OF=os.path.join(SAMPLE,"102_0_openfoam_hollow_cylinder_heat_transfer")
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
K=273.15; NPT=5; RMODE=5


def read_foam_field(path, ncell=None):
    txt=open(path).read(); i=txt.find("internalField"); seg=txt[i:]
    if "uniform" in seg[:40] and "nonuniform" not in seg[:40]:
        v=[float(x) for x in re.findall(r"uniform\s+([-\d.eE ]+);",seg)[0].split()]
        return np.full(ncell,v[0]) if (ncell and len(v)==1) else np.array(v)
    n=int(re.search(r"List<\w+>\s*\n?\s*(\d+)",seg).group(1)); body=seg[seg.find("(")+1:]
    if "vector" in seg[:60] or body.lstrip().startswith("("):
        vals=re.findall(r"\(([-\d.eE ]+)\)",body)[:n]
        return np.array([[float(t) for t in v.split()] for v in vals])
    return np.array([float(x) for x in re.findall(r"[-+]?\d[\d.eE+-]*",body)[:n]])


def load_snapshots():
    times=sorted(int(d) for d in os.listdir(OF)
                 if d.isdigit() and os.path.isdir(os.path.join(OF,d,"solid")))
    C=read_foam_field(os.path.join(OF,str(times[1]),"solid","C")); ncell=len(C)
    cols=[]; kept=[]
    for t in times:
        f=os.path.join(OF,str(t),"solid","T")
        if os.path.exists(f): cols.append(read_foam_field(f,ncell)); kept.append(t)
    X=np.column_stack(cols)
    print(f"[qdeim] snapshots X={X.shape}, t={kept[0]}..{kept[-1]}s")
    return np.array(kept,float), C, X


def main():
    os.makedirs(IMG,exist_ok=True); os.makedirs(RES,exist_ok=True)
    ts, C, X = load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
    energy=S**2/(S**2).sum(); cum=np.cumsum(energy)
    print(f"[qdeim] 累積エネルギー: 1モード{cum[0]*100:.1f}% 2モード{cum[1]*100:.1f}% "
          f"3モード{cum[2]*100:.1f}% 5モード{cum[4]*100:.1f}%")
    # Q-DEIM: 先頭RMODEモード U[:, :r] の転置に枢軸QR → 最初のr枢軸が補間点
    _,_,piv=qr(U[:,:RMODE].T, pivoting=True)
    pts=list(piv[:NPT])
    coords=C[pts]
    print(f"[qdeim] 選ばれた代表点(Q-DEIM, r={RMODE}):")
    for j,(p,xyz) in enumerate(zip(pts,coords)):
        print(f"   点{j}: cell{p} xyz={np.round(xyz*1000,1)}mm")
    np.savez(os.path.join(RES,"qdeim_points.npz"),
             cell_idx=np.array(pts), xyz=coords, cell_centres=C,
             pod_modes=U[:,:RMODE].astype(np.float32), mean=mean.astype(np.float32),
             energy=energy[:10], times=ts,
             note="POD+Q-DEIMで102_0のCHT温度場から選んだROM代表点")

    # 3D図: メッシュ面 + 選点
    import pyvista as pv, vtk
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    mco=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),mco)
    pl=pv.Plotter(off_screen=True,window_size=(950,860))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.5,show_edges=False)
    for j,xyz in enumerate(coords):
        pl.add_mesh(pv.Sphere(radius=0.005,center=xyz),color="crimson")
        p=np.round(xyz*1000,0).astype(int)
        pl.add_point_labels([xyz+np.array([0,0,0.006*(1 if j%2 else -1)*(j//2+1)])],
            [f"P{j} ({p[0]},{p[1]},{p[2]})mm"],font_size=17,text_color="crimson",
            shape=None,always_visible=True)
    pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    pl.camera.zoom(1.25)
    p4=os.path.join(IMG,"_tmp.png"); pl.screenshot(p4); pl.close()
    # 円柱まわりの白い余白を自動トリミング（内容=円柱・球・ラベルの外接矩形＋わずかな余白）
    img=plt.imread(p4)
    mask=np.any(img[...,:3]<0.96,axis=-1)
    ys,xs=np.where(mask)
    pad=14
    y0,y1=max(0,ys.min()-pad),min(img.shape[0],ys.max()+pad)
    x0,x1=max(0,xs.min()-pad),min(img.shape[1],xs.max()+pad)
    img=img[y0:y1,x0:x1]
    h,w=img.shape[:2]
    fig,ax=plt.subplots(figsize=(9.5,9.5*h/w)); ax.imshow(img); ax.axis("off")
    ax.set_title("POD＋Q-DEIM で選んだROM代表点（勘でなくデータから系統的に選定）\n"
                 f"温度場は2モードで{cum[1]*100:.0f}%説明→少数点で表せる",fontsize=14)
    fig.savefig(os.path.join(IMG,"qdeim_points.png"),dpi=130,bbox_inches="tight",pad_inches=0.05)
    plt.close(fig)
    os.remove(p4)
    print("[qdeim] wrote results/qdeim_points.npz, docs/img/qdeim_points.png")


if __name__=="__main__": main()
