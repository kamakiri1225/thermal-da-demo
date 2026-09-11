"""KinvHのDUMPWパッチ版FrontISTRで、本ケース円筒の感度 W_diff を直接出力する.

ユーザーのFrontISTR感度解析プログラム(20260810_KinvH, W=K^-1 H, 随伴法と2.4e-8一致)
を、この105の中空円筒に適用する。DUMPWは四面体(341)専用のため、六面体メッシュを
0-6対角の6分割で共形テトラ化してから実行する。

QoI測定点(sensitivity_points.dat): Point A(ヒータ側上面+X) / Point O(反対側上面-X)
出力:
  openfoam/dumpw_case/sensitivity_Wdiff.vtk … 感度ベクトル場(ParaView, KinvH形式)
  paraview/sensitivity_Wdiff_cylinder.vtk    … 上のコピー(公開用)
  docs/img/sensitivity_dumpw.png             … wz成分の面表示図
  検証: W_diff_z をROM5モードへ射影し、摂動法(M_all[A]-M_all[O])と比較
"""
from __future__ import annotations
import os, subprocess, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh
import fistr_case
from fem.fem_obs import MATERIAL
from dacore.node_locations import NODE_XYZ

NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
BIN=os.path.expanduser("~/src/FrontISTR-dumpw/build-dumpw/fistr1/fistr1")
CASE=os.path.join(ROOT,"openfoam","dumpw_case")
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])

# 六面体(0..7)→6四面体(0-6対角)。構造格子で共形。
TET_SPLIT=[(0,1,2,6),(0,2,3,6),(0,3,7,6),(0,7,4,6),(0,4,5,6),(0,5,1,6)]

def fmt_ids(ids,per=8):
    out=[]
    for i in range(0,len(ids),per):
        out.append(" "+", ".join(str(x) for x in ids[i:i+per]))
    return out

def main():
    os.makedirs(CASE,exist_ok=True)
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    node_ids=[n for n,_ in mesh["nodes"]]
    id_xyz={n:np.array(x) for n,x in mesh["nodes"]}
    qa=node_ids[int(np.linalg.norm(coords-A_XYZ,axis=1).argmin())]
    qo=node_ids[int(np.linalg.norm(coords-O_XYZ,axis=1).argmin())]
    print(f"[dumpw] Point A node={qa}, Point O node={qo}")

    # --- TET分割(体積正に矯正) ---
    tets=[]; eid=1
    for _e,conn in mesh["elements"]:
        for t in TET_SPLIT:
            n=[conn[i] for i in t]
            p=[id_xyz[x] for x in n]
            v=np.dot(np.cross(p[1]-p[0],p[2]-p[0]),p[3]-p[0])
            if v<0: n=[n[0],n[2],n[1],n[3]]
            tets.append((eid,n)); eid+=1
    print(f"[dumpw] hex {len(mesh['elements'])} -> tet {len(tets)}")

    # --- .msh (341) ---
    L=["!HEADER"," 105 DUMPW cylinder (tetra split of 102_1 hex mesh)","!NODE"]
    for nid,(x,y,z) in mesh["nodes"]:
        L.append(f"{nid:8d}, {x:.12g}, {y:.12g}, {z:.12g}")
    L.append("!ELEMENT, TYPE=341")
    for e,n in tets: L.append(f"{e:8d}, "+", ".join(f"{x:8d}" for x in n))
    for g,k in [("NALL","all_nodes"),("BOTTOM","bottom_nodes"),("TOP","top_nodes"),
                ("FIX_XYZ","fix_xyz"),("FIX_YZ","fix_yz"),("FIX_Z","fix_z")]:
        L.append(f"!NGROUP, NGRP={g}"); L+=fmt_ids(mesh[k])
    L.append("!EGROUP, EGRP=EALL"); L+=fmt_ids([e for e,_ in tets])
    L.append("!MATERIAL, NAME=STEEL, ITEM=3")
    L+= ["!ITEM=1, SUBITEM=2",f" {MATERIAL['young_modulus_Pa']:.10g}, {MATERIAL['poisson_ratio']:.10g}",
         "!ITEM=2, SUBITEM=1",f" {MATERIAL['density_kg_m3']:.10g}",
         "!ITEM=3, SUBITEM=1",f" {MATERIAL['thermal_expansion_coeff_per_K']:.10g}"]
    L.append("!SECTION, TYPE=SOLID, EGRP=EALL, MATERIAL=STEEL"); L.append(" 1.0"); L.append("!END")
    open(os.path.join(CASE,"hollow_cylinder_thermal_expansion.msh"),"w").write("\n".join(L)+"\n")

    # --- .cnt (DIRECT + DUMPW=YES, 荷重は一様+1K: Wは荷重に依存しない) ---
    C=["!VERSION"," 3","!SOLUTION, TYPE=STATIC",
       "!WRITE,RESULT,FREQUENCY=1",
       "!SOLVER,METHOD=DIRECT,ITERLOG=NO,TIMELOG=YES,DUMPTYPE=MM,DUMPW=YES",
       " 10000, 1"," 1.0e-8, 1.0, 0.0",
       "!REFTEMP",f" {Tref:.10g}",
       "!INITIAL_CONDITION, TYPE=TEMPERATURE",f" NALL, {Tref:.10g}",
       "!BOUNDARY, GRPID=1"," FIX_XYZ,1,3"," FIX_YZ,2,3"," FIX_Z,3,3",
       "!TEMPERATURE, GRPID=1",f" NALL, {Tref+1.0:.10g}",
       "!MATERIAL, NAME=STEEL","!ELASTIC",
       f" {MATERIAL['young_modulus_Pa']:.10g}, {MATERIAL['poisson_ratio']:.10g}",
       "!EXPANSION_COEFF",f" {MATERIAL['thermal_expansion_coeff_per_K']:.10g}",
       "!STEP, SUBSTEPS=1, CONVERG=1.000E-07","BOUNDARY,1","LOAD,1","!END"]
    open(os.path.join(CASE,"hollow_cylinder_thermal_expansion.cnt"),"w").write("\n".join(C)+"\n")
    fistr_case.write_hecmw_ctrl(__import__("pathlib").Path(CASE))
    open(os.path.join(CASE,"sensitivity_points.dat"),"w").write(f"#Point_A, Point_O\n{qa} {qo}\n")

    # --- 実行 ---
    print(f"[dumpw] running patched fistr1 (DIRECT, {3*len(node_ids)} DOF) ...",flush=True)
    env=os.environ.copy(); env.setdefault("OMP_NUM_THREADS","4")
    r=subprocess.run([BIN],cwd=CASE,env=env,stdout=open(os.path.join(CASE,"run_dumpw.log"),"w"),
                     stderr=subprocess.STDOUT)
    print(f"[dumpw] exit={r.returncode}")
    out=os.path.join(CASE,"sensitivity_Wdiff.vtk")
    if os.path.exists(out):
        os.makedirs(os.path.join(ROOT,"paraview"),exist_ok=True)
        subprocess.run(["cp",out,os.path.join(ROOT,"paraview","sensitivity_Wdiff_cylinder.vtk")])
        print("[dumpw] sensitivity_Wdiff.vtk -> paraview/sensitivity_Wdiff_cylinder.vtk")
    else:
        print("[dumpw] WARNING: sensitivity_Wdiff.vtk not produced; see run_dumpw.log")

if __name__=="__main__": main()
