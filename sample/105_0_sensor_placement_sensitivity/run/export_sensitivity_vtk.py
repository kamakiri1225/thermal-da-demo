"""感度をVTKで出力する(ParaView用, KinvHのsensitivity_Wdiff.vtkと同型の使い方).

FrontISTRのFEMメッシュ(5040節点)に、2種類の感度を節点フィールドとして格納:
  Sens_RowNorm      … Wの行感度 |∂uz(node)/∂mode| の合計 [µm/K]
                       (変位観測点の選定用。FrontISTR 5回の実結果 M_all)
  dUzHeater_dT      … Wの列感度 ∂uz_heater/∂T(x) [µm/K] (パッチ計算をノードへ展開)
  dUzOpp_dT         … 同 ∂uz_opp/∂T(x)
出力:
  paraview/sensitivity_nodal.vtu … ParaViewでそのまま開ける
  docs/img/sensitivity_frontistr.png … 上記から描いた図(FEMメッシュ面表示)
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, vtk
from dacore.node_locations import HEIGHT
from run.displacement_point_selection import build_M_all, NR,NTH,NZ,R_IN,R_OUT,H
IMG=os.path.join(ROOT,"docs","img"); PV=os.path.join(ROOT,"paraview"); RES=os.path.join(ROOT,"results")

def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    os.makedirs(PV,exist_ok=True); os.makedirs(IMG,exist_ok=True)

    coords,M=build_M_all()                     # FrontISTR実計算(5回)
    sens_row=np.abs(M).sum(axis=1)             # 行感度 [µm/K]
    # 列感度(パッチ)をノードへ展開
    d=np.load(os.path.join(RES,"sensitivity_uz.npz"))
    S=d["S_um"]; NTH_P=int(d["NTH_P"]); NZ_P=int(d["NZ_P"])
    th=np.arctan2(coords[:,1],coords[:,0])%(2*np.pi)
    a=np.minimum((th/(2*np.pi)*NTH_P).astype(int),NTH_P-1)
    b=np.minimum((coords[:,2]/HEIGHT*NZ_P).astype(int),NZ_P-1)

    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    cells=[]
    for _e,conn in mesh["elements"]: cells.append(8); cells.extend(idr[n] for n in conn)
    ug=pv.UnstructuredGrid(np.array(cells),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    ug.point_data["Sens_RowNorm_umK"]=sens_row
    ug.point_data["dUzHeater_dT_umK"]=S[0,a,b]
    ug.point_data["dUzOpp_dT_umK"]=S[1,a,b]
    out=os.path.join(PV,"sensitivity_nodal.vtu"); ug.save(out)
    print(f"[vtk] wrote {os.path.relpath(out,ROOT)} (fields: Sens_RowNorm, dUzHeater_dT, dUzOpp_dT)")

    # 図: FrontISTR節点感度(行感度)をFEMメッシュ面で
    pl=pv.Plotter(off_screen=True,shape=(1,2),window_size=(1500,760))
    for j,(name,fld,cmap) in enumerate([
        ("W row-norm: pick DISPLACEMENT obs points","Sens_RowNorm_umK","viridis"),
        ("W column: dUz_heater/dT(x) (temp influence)","dUzHeater_dT_umK","coolwarm")]):
        pl.subplot(0,j)
        clim=None
        if "coolwarm"==cmap:
            v=float(np.abs(ug[fld]).max()); clim=[-v,v]
        pl.add_mesh(ug.copy(),scalars=fld,cmap=cmap,clim=clim,show_edges=False,
                    scalar_bar_args={"title":fld,"title_font_size":17,"label_font_size":13})
        pl.add_text(name,font_size=15,color="black")
        pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]
    pl.set_background("white")
    fig=os.path.join(IMG,"sensitivity_frontistr.png"); pl.screenshot(fig); pl.close()
    print(f"[vtk] wrote {os.path.relpath(fig,ROOT)}")

if __name__=="__main__": main()
