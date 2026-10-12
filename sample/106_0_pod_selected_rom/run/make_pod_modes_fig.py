"""PODで出た「型（モード）」を絵にする（平均場＋モード1〜5）."""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity"))
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
sys.path.insert(0, ROOT)
from dacore.temperature_mapping import build_idw_mapping, map_temperature, mapping_cache_valid, MAPPING_METHOD
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img"); K=273.15


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float)
    Cc=kv["cell_centres"]; energy=kv["energy"]
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,0.1005)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    temperature_map = build_idw_mapping(Cc, coords)  # 近傍8セルのIDW
    CAM=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; TMP=tempfile.mkdtemp()
    panels=[("平均場",mean-K,"turbo",False)]+[
        (f"モード{k+1} ({energy[k]*100:.3g}%)",U[:,k],"coolwarm",True) for k in range(5)]
    shots=[]
    for i,(lab,fld,cmap,sym) in enumerate(panels):
        fld=fld.copy()
        if sym:
            # PODモードの符号は任意。最大振幅成分が正になるよう向きを揃える
            # （モード1＝一様昇温は正＝暖色で「全体が温まる」を直感的に）。
            if fld[np.argmax(np.abs(fld))]<0: fld=-fld
            # 一様に近い型でも「うすい単色」に見えないよう、分布の主部が飽和するclim
            # （外れ値のp85で頭打ち。型1は一様なので全面が濃い暖色、型2・3は±構造が残る）。
            a=float(np.abs(fld).max()) or float(np.abs(fld).max())
            clim=[-a,a]
        else:
            clim=[float((mean-K).min()),float((mean-K).max())]
        g=ug.copy(); g.point_data["v"]=map_temperature(temperature_map, fld)
        pl=pv.Plotter(off_screen=True,window_size=(560,640))
        pl.add_mesh(g,scalars="v",cmap=cmap,clim=clim,n_colors=12,show_scalar_bar=False)
        pl.camera_position=CAM; pl.set_background("white")
        pl.camera.zoom(1.18)      # 少し引いて上下の見切れ（円柱が枠に接する）を防ぐ
        p=os.path.join(TMP,f"m{i}.png"); pl.screenshot(p); pl.close(); shots.append((p,lab,clim,cmap,sym))
    fig,axes=plt.subplots(2,3,figsize=(13,10))
    axes=axes.ravel()
    tt=["実際の温度[℃]","全体が同じだけ昇温(一様→全面が同じ暖色)","ヒータ側(赤)↔反対側(青)","高次の細部", "高次の細部", "高次の細部"]
    for ax,(p,lab,clim,cmap,sym),sub in zip(axes,shots,tt):
        ax.imshow(plt.imread(p)); ax.axis("off"); ax.set_title(lab,fontsize=14)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*clim),cmap=cmap),ax=ax,orientation="horizontal",fraction=.045,pad=.01)
        cb.set_label("正規化モード成分 [無次元]" if sym else "平均温度 [℃]",fontsize=10)
        cb.ax.tick_params(labelsize=9)
    fig.suptitle("PODで出た温度場の「型」: 平均場とモード1〜5",fontsize=16,weight="bold")
    fig.text(0.5,0.055,"※ 各パネルは独立スケール（各モードの最大絶対値を色範囲に使用）。モード(型)は単位長さに正規化した『形』で絶対温度ではない。"
             "重要度は色の濃さでなく %（型1=91.8%が最も支配的）。",ha="center",fontsize=10.5,color="dimgray")
    fig.text(0.5,0.02,"型1は『全体が同じだけ温まる』一様な型なので全面がほぼ同じ暖色（凹凸のある型2・3と対照的）。"
             "モードの符号は任意で、ここでは昇温を暖色に揃えた。",ha="center",fontsize=10.5,color="dimgray")
    fig.tight_layout(rect=[0,0.06,1,0.93]); fig.savefig(os.path.join(IMG,"pod_modes.png"),dpi=130); plt.close(fig)
    print("[modes] wrote docs/img/pod_modes.png")


if __name__=="__main__": main()
