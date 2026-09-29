"""結果③の「変位差 Uz(A)-Uz(O)」がどこの2点なのかを図で示す.

左 : 中空円筒の上面にある変位観測点 A(+X, ヒータ側) と O(-X, 反対側)、
     および温度の代表点 P0〜P4 の位置（ParaView風）
右 : 熱変形（誇張）に A・O を重ね、Uz(A)-Uz(O) が「反り」であることを示す

出力: docs/img/disp_obs_points.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_disp_obs_points_fig.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh, vtk
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
VTU=os.path.join(SAMPLE,"104_0_openfoam_frontistr_da_enkf","paraview","displacement_fields.vtu")
H=0.1005; A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H]); EXAG=6000.0


def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    Cc=kv["cell_centres"]; pod=kv["cell_idx"]
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,0.020,0.0375,H)
    coords=np.array([p for _n,p in mesh["nodes"]])
    idr={nid:i for i,(nid,_x) in enumerate(mesh["nodes"])}
    conn=[]
    for _e,cc in mesh["elements"]: conn.append(8); conn.extend(idr[n] for n in cc)
    ug=pv.UnstructuredGrid(np.array(conn),
        np.full(len(mesh["elements"]),vtk.VTK_HEXAHEDRON,np.uint8),coords)
    CAM=[(0.27,-0.25,0.20),(0,0,0.05),(0,0,1)]

    def trim(img):
        m=np.any(img[:,:,:3]<246,axis=-1); ys,xs=np.where(m); p=6
        return img[max(0,ys.min()-p):ys.max()+p, max(0,xs.min()-p):xs.max()+p]

    # 左：観測点の位置
    pl=pv.Plotter(off_screen=True,window_size=(720,800))
    pl.add_mesh(ug,color="#c3cedb",opacity=0.26)
    for i,c in enumerate(pod):
        pl.add_mesh(pv.Sphere(radius=0.0028,center=Cc[c]),color="#1F4E9C")
        pl.add_point_labels([Cc[c]],[f"P{i}"],font_size=30,text_color="#1F4E9C",
                            shape=None,always_visible=True,bold=True)
    for xyz,lab,col in [(A_XYZ,"A  (+X ヒータ側)","#C0392B"),(O_XYZ,"O  (−X 反対側)","#2e7d32")]:
        pl.add_mesh(pv.Sphere(radius=0.0040,center=xyz),color=col)
        pl.add_point_labels([xyz+np.array([0,0,0.007])],[lab],font_size=34,text_color=col,
                            shape=None,always_visible=True,bold=True)
        pl.add_mesh(pv.Arrow(start=xyz,direction=(0,0,1),scale=0.030),color=col)
    pl.camera_position=CAM; pl.set_background("white"); pl.camera.zoom(1.30)
    left=trim(pl.screenshot(return_img=True)); pl.close()

    # 右：熱変形にA・Oを重ねる
    m=pv.read(VTU); m.set_active_vectors("U_truth_m")
    uz=np.asarray(m.point_data["Uz_truth_um"])
    w=m.warp_by_vector("U_truth_m",factor=EXAG)
    pl=pv.Plotter(off_screen=True,window_size=(720,800))
    pl.add_mesh(w,scalars="Uz_truth_um",cmap="coolwarm",
                scalar_bar_args={"title":"Uz [um]","fmt":"%.1f",
                                 "title_font_size":18,"label_font_size":15})
    pl.add_mesh(m,color="lightgray",opacity=0.12)
    for xyz,lab,col in [(A_XYZ,"A","#C0392B"),(O_XYZ,"O","#2e7d32")]:
        j=int(np.linalg.norm(m.points-xyz,axis=1).argmin())
        pt=w.points[j]
        pl.add_mesh(pv.Sphere(radius=0.0040,center=pt),color=col)
        pl.add_point_labels([pt+np.array([0,0,0.008])],[lab],font_size=38,text_color=col,
                            shape=None,always_visible=True,bold=True)
    pl.camera_position=CAM; pl.set_background("white"); pl.camera.zoom(1.30)
    right=trim(pl.screenshot(return_img=True)); pl.close()

    fig,(a0,a1)=plt.subplots(1,2,figsize=(11.2,6.6))
    a0.imshow(left); a0.axis("off")
    a0.set_title("変位観測点 A・O の位置（上面 z=100.5 mm）\n青 P0〜P4＝温度の代表点",
                 fontsize=13, weight="bold")
    a1.imshow(right); a1.axis("off")
    a1.set_title(f"熱変形（×{EXAG:.0f}誇張）に重ねると\nA側が持ち上がり O側が下がる＝「反り」",
                 fontsize=13, weight="bold")
    fig.suptitle("「変位差 $U_z(A)-U_z(O)$」とはどこの2点か",
                 fontsize=15, weight="bold")
    fig.text(0.5,0.018,"A＝上面のヒータ側(+X)、O＝上面の反対側(−X)。"
             "2点の上下方向変位の差が、上面の傾き（反り）を表す。",
             ha="center",fontsize=11.5,color="#444")
    fig.tight_layout(rect=[0,0.035,1,0.92])
    out=os.path.join(IMG,"disp_obs_points.png")
    fig.savefig(out,dpi=135); plt.close(fig); print("wrote",out)


if __name__=="__main__":
    main()
