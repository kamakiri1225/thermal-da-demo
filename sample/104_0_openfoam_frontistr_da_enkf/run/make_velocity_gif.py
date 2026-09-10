"""102_0 のCHTから、流体の自然対流の流速分布アニメーション(0→600s)を作る。

chtMultiRegionFoam は流体側でナビエ・ストークス+エネルギー+浮力を解いている。
その速度場 U を、鉛直断面(y=0, ヒータ+X〜反対-X)で |U| として時系列表示し、
矢印で流れの向きも重ねてGIFにする。0-300s加熱→300-600s遮断。

出力: docs/img/fluid_velocity.gif（と presentation/assets へコピー）
使い方: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_velocity_gif.py
"""
from __future__ import annotations
import glob, os, subprocess
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FOAM = os.path.abspath(os.path.join(
    ROOT, "..", "102_0_openfoam_hollow_cylinder_heat_transfer", "hollowCylinder.foam"))
IMG = os.path.join(ROOT, "docs", "img")
TMP = os.path.join(ROOT, "openfoam", "vel_frames")


def main():
    import pyvista as pv
    pv.OFF_SCREEN = True
    try:
        pv.start_xvfb()
    except Exception:
        pass
    os.makedirs(TMP, exist_ok=True); os.makedirs(IMG, exist_ok=True)

    r = pv.POpenFOAMReader(FOAM)
    times = [t for t in r.time_values if abs(t % 20) < 1e-6]  # 20sごと=31枚
    clim = [0.0, 0.08]  # m/s

    for k, t in enumerate(times):
        r.set_active_time_value(t)
        blk = r.read()
        fl = blk["fluid"]
        g = fl["internalMesh"] if "internalMesh" in fl.keys() else fl.combine()
        g = g.cell_data_to_point_data()
        g["Umag"] = np.linalg.norm(g["U"], axis=1)
        sl = g.slice(normal="y", origin=(0, 0, 0.05025))
        # 固体を半透明で文脈表示(温度3D図と同じアングル)
        sol = blk["solid"]
        sg = sol["internalMesh"] if "internalMesh" in sol.keys() else sol.combine()
        sclip = sg.clip(normal="y", origin=(0, 0, 0.05025))

        pl = pv.Plotter(off_screen=True, window_size=(780, 840))
        pl.add_mesh(sclip, color="lightgray", opacity=0.30)
        pl.add_mesh(sl, scalars="Umag", cmap="turbo", clim=clim, opacity=0.92,
                    scalar_bar_args={"title": "|U| [m/s]", "title_font_size": 20,
                                     "label_font_size": 16})
        arrows = sl.glyph(orient="U", scale=False, factor=0.010, tolerance=0.025)
        pl.add_mesh(arrows, color="black")
        phase = "加熱中 (0-300s)" if t < 300 else "遮断後 (300-600s)"
        pl.add_text(f"t = {t:g} s   {phase}", position="upper_edge",
                    font_size=16, color="black")
        pl.set_background("white")
        pl.camera_position = [(0.26, -0.24, 0.22), (0.0, -0.01, 0.05), (0, 0, 1)]
        pl.camera.zoom(1.25)
        pl.screenshot(os.path.join(TMP, f"v_{k:03d}.png"))
        pl.close()
        print(f"[vel] frame {k+1}/{len(times)}  t={t:g}s", flush=True)

    frames = sorted(glob.glob(os.path.join(TMP, "v_*.png")))
    out = os.path.join(IMG, "fluid_velocity.gif")
    subprocess.run(["convert", "-delay", "18", "-loop", "0", *frames,
                    "-resize", "620x", "-layers", "Optimize", "-colors", "96", out], check=True)
    for f in frames:
        os.remove(f)
    # プレゼンにもコピー
    passets = os.path.join(ROOT, "presentation", "assets", "fluid_velocity.gif")
    subprocess.run(["cp", out, passets], check=False)
    print(f"[vel] wrote {os.path.relpath(out, ROOT)} ({len(frames)} frames)")


if __name__ == "__main__":
    main()
