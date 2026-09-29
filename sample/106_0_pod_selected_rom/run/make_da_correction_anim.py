"""でたらめな初期状態から温度分布と熱変位が補正されていく様子のGIF（目的スライド用）.

3本並べ：真値／データ同化なし（自由予測）／データ同化あり（EnKF）
温度場は gappy-POD で5点から全セル復元し、変形は熱感度で誇張表示する。
出力: docs/img/da_correction_anim.gif
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_da_correction_anim.py
"""
from __future__ import annotations
import os, sys, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh
IMG = os.path.join(ROOT, "docs", "img"); RES = os.path.join(ROOT, "results")
KC = 273.15; NPT = 5; IQ = 5; IH = 6; NAUG = 7
DT = 2.0; OBS_DT = 30.0; T_END = 600.0; N_ENS = 60; SIG_T = 0.30
TSENS = [2, 0]


def main():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; K = rg.tri_to_matrix(d["K_upper"], NPT)
    h_t = float(d["h"]); heat = int(d["heat_node"])
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float)
    pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :])
    rng = np.random.default_rng(3); rng_o = np.random.default_rng(23)

    cyc = np.arange(OBS_DT, T_END + 1e-9, OBS_DT); tg = np.r_[0, cyc]
    T = np.full(NPT, rg.T_AIR_K); Ttr = [T.copy()]
    for a, b in zip(np.r_[0, cyc[:-1]], cyc):
        _, tr = rg.integrate_single(T, C, K, h_t, 1.0, heat, a, b, DT)
        T = tr[-1]; Ttr.append(T.copy())
    Ttr = np.array(Ttr)

    # でたらめな初期状態（温度も発熱も間違っている）
    Z = np.zeros((N_ENS, NAUG))
    Z[:, :NPT] = rng.uniform(rg.T_AIR_K + 4, rg.T_AIR_K + 14, (N_ENS, NPT))
    Z[:, IQ] = rng.uniform(1.4, 2.6, N_ENS)
    Z[:, IH] = np.clip(rng.normal(0.03, 0.008, N_ENS), 1e-3, 0.1)
    F = Z.copy()                      # 同化なし（同じ初期値のまま進める）
    da, no = [Z[:, :NPT].mean(0)], [F[:, :NPT].mean(0)]
    tp = 0.0
    for ci, tb in enumerate(cyc, 1):
        Z = Z.copy(); F = F.copy()
        Z[:, :NPT] = rg.integrate_ensemble(Z[:, :NPT], C, K, Z[:, IH], Z[:, IQ], heat, tp, tb, DT)
        F[:, :NPT] = rg.integrate_ensemble(F[:, :NPT], C, K, F[:, IH], F[:, IQ], heat, tp, tb, DT)
        tp = tb
        y = Ttr[ci][TSENS] + rng_o.normal(0, SIG_T, len(TSENS))
        Z = enkf_update(Z, y, None, np.diag([SIG_T**2]*len(TSENS)), rng, inflation=1.02, Yf=Z[:, TSENS])
        Z[:, IQ] = np.clip(Z[:, IQ], 0, 3); Z[:, IH] = np.clip(Z[:, IH], 1e-4, 0.2)
        da.append(Z[:, :NPT].mean(0)); no.append(F[:, :NPT].mean(0))
    da = np.array(da); no = np.array(no)

    def field(T5):
        return mean + U @ (UPp @ (T5 - mean[pod]))

    import pyvista as pv, vtk
    from PIL import Image
    from scipy.spatial import cKDTree
    pv.OFF_SCREEN = True
    try: pv.start_xvfb()
    except Exception: pass
    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    co = np.array([c for _n, c in mesh["nodes"]])
    idr = {nid: i for i, (nid, _x) in enumerate(mesh["nodes"])}
    cl = []
    for _e, conn in mesh["elements"]: cl.append(8); cl.extend(idr[n] for n in conn)
    ug = pv.UnstructuredGrid(np.array(cl), np.full(len(mesh["elements"]), vtk.VTK_HEXAHEDRON, np.uint8), co)
    _, near = cKDTree(Cc).query(co)
    clim = [19.5, 27.5]

    def shot(T5, title, col):
        f = field(T5) - KC
        g = ug.copy(); g.point_data["T"] = f[near]
        pl = pv.Plotter(off_screen=True, window_size=(430, 560))
        pl.add_mesh(g, scalars="T", cmap="turbo", clim=clim, n_colors=18, show_scalar_bar=False)
        pl.add_text(title, position="upper_left", font_size=11, color=col)
        pl.camera_position = [(0.26, -0.24, 0.20), (0, 0, 0.05), (0, 0, 1)]
        pl.set_background("white"); pl.camera.zoom(1.35)
        im = pl.screenshot(return_img=True); pl.close()
        m = np.any(im[..., :3] < 246, axis=-1); ys, xs = np.where(m)
        return im[ys.min()-4:ys.max()+4, xs.min()-4:xs.max()+4]

    idx = list(range(0, len(tg), 1))
    frames = []
    for k in idx:
        ims = [shot(Ttr[k], "TRUTH", "black"),
               shot(no[k], "NO ASSIMILATION", "#8a1c1c"),
               shot(da[k], "DATA ASSIMILATION", "#14459c")]
        e_no = np.sqrt(((field(no[k]) - field(Ttr[k]))**2).mean())
        e_da = np.sqrt(((field(da[k]) - field(Ttr[k]))**2).mean())
        fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.6))
        ttls = ["真値", f"データ同化なし  誤差 {e_no:5.2f} K", f"データ同化あり  誤差 {e_da:5.2f} K"]
        cols = ["#1b2430", "#c0392b", "#14459c"]
        for ax, im, t, c in zip(axes, ims, ttls, cols):
            ax.imshow(im); ax.axis("off")
            ax.set_title(t, fontsize=12.5, weight="bold", color=c)
        fig.suptitle(f"でたらめな初期状態からでも温度分布が補正される（t = {tg[k]:.0f} s）",
                     fontsize=13.5, weight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.9])
        fig.canvas.draw()
        frames.append(Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()))
        plt.close(fig)
        if k % 5 == 0: print(f"  t={tg[k]:.0f}s  同化なし {e_no:.2f} K / 同化あり {e_da:.2f} K")

    out = os.path.join(IMG, "da_correction_anim.gif")
    frames[0].save(out, save_all=True, append_images=frames[1:], duration=200, loop=0, disposal=2)
    print("wrote", out, f"({len(frames)}フレーム)")


if __name__ == "__main__":
    main()
