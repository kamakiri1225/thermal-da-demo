"""102_0 のCHTから、固体温度場の0→600sアニメーション(断面)を作る。
出力: docs/img/solid_temperature.gif（presentation/assets にもコピー）
"""
from __future__ import annotations
import glob, os, subprocess
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
FOAM=os.path.abspath(os.path.join(ROOT,"..","102_0_openfoam_hollow_cylinder_heat_transfer","hollowCylinder.foam"))
IMG=os.path.join(ROOT,"docs","img"); TMP=os.path.join(ROOT,"openfoam","temp_frames")
K=273.15
def main():
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    os.makedirs(TMP,exist_ok=True); os.makedirs(IMG,exist_ok=True)
    r=pv.POpenFOAMReader(FOAM)
    times=[t for t in r.time_values if abs(t%20)<1e-6]
    clim=[20.0,26.0]  # degC
    for k,t in enumerate(times):
        r.set_active_time_value(t)
        sol=r.read()["solid"]
        g=sol["internalMesh"] if "internalMesh" in sol.keys() else sol.combine()
        g=g.cell_data_to_point_data()
        g["T_degC"]=g["T"]-K
        clip=g.clip(normal="y",origin=(0,0,0.05025))
        pl=pv.Plotter(off_screen=True,window_size=(720,820))
        pl.add_mesh(clip,scalars="T_degC",cmap="turbo",clim=clim,
                    scalar_bar_args={"title":"T [degC]","title_font_size":20,"label_font_size":16})
        phase="加熱中 (0-300s)" if t<300 else "遮断後 (300-600s)"
        pl.add_text(f"t = {t:g} s   {phase}",position="upper_edge",font_size=16,color="black")
        pl.set_background("white"); pl.camera_position="xz"; pl.camera.zoom(1.35)
        pl.screenshot(os.path.join(TMP,f"t_{k:03d}.png")); pl.close()
        print(f"[temp] frame {k+1}/{len(times)} t={t:g}s",flush=True)
    frames=sorted(glob.glob(os.path.join(TMP,"t_*.png")))
    out=os.path.join(IMG,"solid_temperature.gif")
    subprocess.run(["convert","-delay","18","-loop","0",*frames,"-resize","600x",
                    "-layers","Optimize","-colors","96",out],check=True)
    for f in frames: os.remove(f)
    subprocess.run(["cp",out,os.path.join(ROOT,"presentation","assets","solid_temperature.gif")],check=False)
    print(f"[temp] wrote {os.path.relpath(out,ROOT)}")
if __name__=="__main__": main()
