"""PODで「ROMに何自由度必要か」「5点ROMで温度分布は失われないか」を実データで検証する.

元データ: 102_0 の chtMultiRegionFoam 固体温度場(20,696セル)を0→600s・5秒間隔で
121スナップショット読み込み、スナップショット行列 X(20696×121) を作る。

生成物(docs/img/):
  pod_spectrum.png       … 特異値スペクトルと累積エネルギー(何モードで何%説明できるか)
  pod_modes.png          … 先頭POD モードの空間分布(何を表すモードか)
  pod_reconstruction.png … モード数ごとの場の再構成誤差
  rom_5points.png        … 5点ROMの代表点が円筒のどこか(3D)
  rom_idw_error.png      … 5点→IDW復元した温度分布の、全セル真値に対する誤差の時刻歴
出力: results/pod_summary.yaml
"""
from __future__ import annotations
import os, re, sys, glob
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p   # 日本語フォント登録
import matplotlib.pyplot as plt
from dacore.node_locations import NODE_XYZ
from dacore.cht_rom import NODE_NAMES

OF=os.path.join(SAMPLE,"102_0_openfoam_hollow_cylinder_heat_transfer")
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
K=273.15


def read_foam_field(path,ncell=None):
    """OpenFOAMのvolScalarField/vectorFieldのinternalFieldを配列で返す.
    一様場は ncell 指定時に全セルへ展開する."""
    txt=open(path).read()
    i=txt.find("internalField")
    seg=txt[i:]
    if "uniform" in seg[:40] and "nonuniform" not in seg[:40]:
        v=re.findall(r"uniform\s+([-\d.eE ]+);",seg)[0].split()
        val=np.array([float(x) for x in v])
        if ncell is not None and val.size==1: return np.full(ncell,val[0])
        return val
    n=int(re.search(r"List<\w+>\s*\n?\s*(\d+)",seg).group(1))
    body=seg[seg.find("(")+1:]
    if "(" in body[:body.find("\n2")] if False else False:
        pass
    # scalar: 1値/行, vector: (x y z)/行
    if "vector" in seg[:60] or body.lstrip().startswith("("):
        vals=re.findall(r"\(([-\d.eE ]+)\)",body)[:n]
        return np.array([[float(t) for t in v.split()] for v in vals])
    nums=re.findall(r"[-+]?\d[\d.eE+-]*",body)
    return np.array([float(x) for x in nums[:n]])


def load_snapshots():
    times=sorted([int(d) for d in os.listdir(OF)
                  if d.isdigit() and os.path.isdir(os.path.join(OF,d,"solid"))])
    C=read_foam_field(os.path.join(OF,str(times[1]),"solid","C"))  # セル中心(m)
    ncell=len(C)
    cols=[]; kept=[]
    for t in times:
        f=os.path.join(OF,str(t),"solid","T")
        if os.path.exists(f):
            cols.append(read_foam_field(f,ncell=ncell)); kept.append(t)
    times=kept
    X=np.column_stack(cols)   # (Ncell, Ntime) [K]
    print(f"[pod] snapshots X={X.shape}, t={times[0]}..{times[-1]}s")
    return np.array(times,float), C, X


def main():
    os.makedirs(IMG,exist_ok=True); os.makedirs(RES,exist_ok=True)
    import yaml
    ts, C, X = load_snapshots()
    Xc=X-X.mean(axis=1,keepdims=True)          # 平均場を引く(変動成分)
    U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
    energy=S**2/ (S**2).sum()
    cum=np.cumsum(energy)
    n90=int(np.searchsorted(cum,0.90))+1
    n99=int(np.searchsorted(cum,0.99))+1
    n999=int(np.searchsorted(cum,0.999))+1
    print(f"[pod] modes for 90%={n90}, 99%={n99}, 99.9%={n999}")

    # --- 図1: 特異値スペクトルと累積エネルギー ---
    fig,(a1,a2)=plt.subplots(1,2,figsize=(13,5))
    a1.semilogy(np.arange(1,len(S)+1),energy,"o-",ms=5,color="tab:blue")
    a1.set_xlabel("モード番号"); a1.set_ylabel("各モードのエネルギー比(対数)")
    a1.set_title("特異値スペクトル:少数モードに集中"); a1.grid(alpha=0.3,which="both")
    a1.set_xlim(0,20); a1.set_ylim(1e-14,1)
    a2.plot(np.arange(1,len(cum)+1),cum*100,"o-",ms=5,color="tab:red")
    for n,c,xy,dxy in [(n90,90,(45,14),None),(n99,99,(60,-30),None),(n999,99.9,(80,-52),None)]:
        a2.axhline(c,ls="--",color="gray",lw=1)
        a2.annotate(f"{c}% ← {n}モードで達成",(n,c),textcoords="offset points",
                    xytext=xy,fontsize=11,
                    arrowprops=dict(arrowstyle="->",color="gray",lw=1))
    a2.set_xlabel("使うモード数"); a2.set_ylabel("説明できる場の変動 [%]")
    a2.set_title("累積エネルギー:数モードで場をほぼ説明",pad=12); a2.grid(alpha=0.3)
    a2.set_xlim(0,20); a2.set_ylim(88,100.2)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"pod_spectrum.png"),dpi=140); plt.close(fig)

    # --- 図2: 先頭POD モードの空間分布(面表示) ---
    import pyvista as pv, vtk, cylinder_mesh
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    # 102_1の六面体メッシュ節点へ最近傍でモードを写して面表示
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
    import tempfile; TMPD=tempfile.mkdtemp(); shots=[]
    meanfield=X.mean(axis=1)
    fields=[("平均場",meanfield-K,"turbo",None)]+[
        (f"モード{k+1} ({energy[k]*100:.1f}%)",U[:,k],"coolwarm",k) for k in range(3)]
    for mi,(lab,fld,cmap,k) in enumerate(fields):
        g=ug.copy(); g.point_data["v"]=fld[near]
        pl=pv.Plotter(off_screen=True,window_size=(600,620))
        if cmap=="coolwarm":
            v=float(np.abs(fld).max()); clim=[-v,v]
        else: clim=[float((meanfield-K).min()),float((meanfield-K).max())]
        pl.add_mesh(g,scalars="v",cmap=cmap,clim=clim,n_colors=12,show_scalar_bar=False)
        pl.camera_position=CAM; pl.set_background("white")
        p=os.path.join(TMPD,f"mode{mi}.png"); pl.screenshot(p); pl.close(); shots.append((p,lab))
    fig,axes=plt.subplots(1,4,figsize=(16,5))
    for ax,(p,lab) in zip(axes,shots):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(lab,fontsize=13)
    fig.suptitle("PODモード:場の変動を作っている「基本パターン」(平均場＋モード1,2,3で大半を説明)",fontsize=15)
    fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(os.path.join(IMG,"pod_modes.png"),dpi=130); plt.close(fig)

    # --- 図3: モード数ごとの場の再構成誤差 ---
    errs=[]
    for r in range(1,16):
        Xr=X.mean(axis=1,keepdims=True)+U[:,:r]@np.diag(S[:r])@Vt[:r,:]
        errs.append(np.sqrt(((Xr-X)**2).mean()))
    fig,ax=plt.subplots(figsize=(9,5))
    ax.semilogy(range(1,16),errs,"o-",color="tab:purple",ms=6)
    ax.axvline(5,ls="--",color="tab:green",lw=2)
    ax.annotate("5点ROMの自由度=5",(5,errs[4]),textcoords="offset points",
                xytext=(12,20),fontsize=12,color="tab:green")
    ax.set_xlabel("使うモード数"); ax.set_ylabel("全セル温度場のRMSE [K]")
    ax.set_title("何モードで温度分布をどこまで復元できるか"); ax.grid(alpha=0.3,which="both")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"pod_reconstruction.png"),dpi=140); plt.close(fig)

    # --- 図4: 5点ROMの位置(3D, ラベル付き) ---
    romn=np.array([NODE_XYZ[n] for n in NODE_NAMES])
    pl=pv.Plotter(off_screen=True,window_size=(950,860))
    pl.add_mesh(ug,color="lightsteelblue",opacity=0.45,show_edges=False)
    cols=["red","tab:green","blue","orange","purple"]
    colmap={"red":"red","tab:green":"green","blue":"blue","orange":"orange","purple":"purple"}
    dz={"hot":0.006,"mid":0.02,"cold":0.006,"top":0.014,"core":-0.01}
    for nm,p,c in zip(NODE_NAMES,romn,cols):
        cc=colmap[c]
        pl.add_mesh(pv.Sphere(radius=0.005,center=p),color=cc)
        xyz=np.round(np.array(p)*1000,0).astype(int)
        pl.add_point_labels([np.array(p)+np.array([0,0,dz[nm]])],
            [f"{nm} ({xyz[0]},{xyz[1]},{xyz[2]})mm"],font_size=18,text_color=cc,
            shape=None,always_visible=True)
    pl.camera_position=CAM; pl.set_background("white")
    p4=os.path.join(TMPD,"rom5.png"); pl.screenshot(p4); pl.close()
    fig,ax=plt.subplots(figsize=(9.5,9.0)); ax.imshow(plt.imread(p4)); ax.axis("off")
    ax.set_title("5点ROMの代表点:これだけで円筒全体の温度分布を表す\n"
                 "hot=ヒータ側(+X)  cold=反対側(−X)  mid=+Y  core=−Y  top=上部(+X)",
                 fontsize=14)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"rom_5points.png"),dpi=130); plt.close(fig)

    # --- 図5: 5点→IDW復元 vs 全セル真値の誤差時刻歴(「分布は失われるか」への答え) ---
    # 5点の位置に最も近いセルの温度を「ROMの5温度」とし、IDWで全セルへ復元して真値と比較
    romcell=[int(np.linalg.norm(C-np.array(p),axis=1).argmin()) for p in romn]
    Wd=np.zeros((len(C),5))
    for i,p in enumerate(C):
        d=np.linalg.norm(romn-p,axis=1)
        if d.min()<1e-9: Wd[i,d.argmin()]=1
        else: w=1/d**2; Wd[i]=w/w.sum()
    idw_err=[]; nmodes_err5=[]
    for j in range(X.shape[1]):
        t5=X[romcell,j]
        rec=Wd@t5
        idw_err.append(np.sqrt(((rec-X[:,j])**2).mean()))
    # 比較用: 5モードPOD の各時刻誤差
    Xr5=X.mean(axis=1,keepdims=True)+U[:,:5]@np.diag(S[:5])@Vt[:5,:]
    pod5_err=np.sqrt(((Xr5-X)**2).mean(axis=0))
    fig,ax=plt.subplots(figsize=(11,5.2)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(ts,idw_err,"-o",ms=3,color="tab:blue",label="5点ROM→IDW復元 の誤差")
    ax.plot(ts,pod5_err,"-s",ms=3,color="tab:purple",label="5モードPOD復元 の誤差(理論上の下限)")
    ax.set_xlabel("time [s]"); ax.set_ylabel("全セル温度分布のRMSE [K]")
    ax.set_title("5点ROMで温度分布は失われるか:全20,696セルとの誤差(橙帯=加熱期)")
    ax.grid(alpha=0.3); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"rom_idw_error.png"),dpi=140); plt.close(fig)

    dT=X.max()-X.min()
    yaml.safe_dump({
        "snapshots":int(X.shape[1]),"cells":int(X.shape[0]),
        "modes_for_90pct":n90,"modes_for_99pct":n99,"modes_for_99_9pct":n999,
        "field_temp_span_K":float(dT),
        "idw_5pt_rmse_max_K":float(np.max(idw_err)),
        "idw_5pt_rmse_mean_K":float(np.mean(idw_err)),
        "pod5_rmse_max_K":float(np.max(pod5_err)),
    },open(os.path.join(RES,"pod_summary.yaml"),"w"),allow_unicode=True)
    print(f"[pod] IDW 5pt rmse max={np.max(idw_err):.4f}K mean={np.mean(idw_err):.4f}K "
          f"(span {dT:.1f}K); POD5 max={np.max(pod5_err):.4f}K")
    print("[pod] wrote 5 figs + results/pod_summary.yaml")


if __name__=="__main__": main()
