"""Q-DEIM の選び方の解説アニメーション（実際の円筒、POD モード 5 本 → 代表 5 点）。

  列の入れ替え付き QR（Φ_rᵀ Π = Q R）は、次の貪欲法と同じ点を同じ順に選ぶ（ここで確かめる）。
    各セル i に、5 つの型の値を並べたベクトル c_i = (φ1(i), …, φ5(i)) を持たせる。
    1. 長さ |c_i| が最大のセルを選ぶ
    2. 選んだセルの向き q = c_sel/|c_sel| の成分を、全セルの c_i から引く（似ている分を消す）
    3. 残りの長さが最大のセルを次に選ぶ。これを 5 回くりかえす
  色＝その時点の残りの長さ（最大で 1）。

  python3 run/make_qdeim_animation.py → docs/img/qdeim_animation.gif、.mp4、qdeim_animation_last.png
"""
import os, sys, subprocess
import numpy as np
from scipy.linalg import qr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
import cylinder_mesh
from dacore.temperature_mapping import build_idw_mapping
import pyvista as pv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img"); WORK = os.path.join(ROOT, "results", "qdeim_anim_work")
NAMES = ["P0", "P1", "P2", "P3", "P4"]


def grid():
    m = cylinder_mesh.build_cylinder_mesh(4, 48, 20, .020, .0375, .1005)
    xyz = np.array([p for _, p in m["nodes"]])
    idx = {n: i for i, (n, _) in enumerate(m["nodes"])}
    cells = []
    for _, conn in m["elements"]:
        cells += [8] + [idx[n] for n in conn[:8]]
    return pv.UnstructuredGrid(np.array(cells), np.full(len(m["elements"]), 12, dtype=np.uint8), xyz * 1000), xyz


def render(g, val, picks, Cc, new=None):
    p = pv.Plotter(off_screen=True, window_size=(900, 900)); p.set_background("white")
    p.add_mesh(g.extract_feature_edges(), color="#666666", line_width=1.2)
    cloud = pv.PolyData(Cc * 1000); cloud["v"] = val
    p.add_mesh(cloud, scalars="v", cmap="viridis", clim=(0, 1), point_size=6, render_points_as_spheres=True,
               show_scalar_bar=False, opacity=0.55)
    for k, c in enumerate(picks):
        col = "#ff3030" if c == new else "#ffffff"
        p.add_mesh(pv.Sphere(radius=3.6 if c == new else 3.0, center=Cc[c] * 1000), color=col)
    p.camera_position = [(170, -210, 190), (0, 0, 48), (0, 0, 1)]; p.camera.zoom(1.25)
    img = p.screenshot(return_img=True); p.close(); return img


def main():
    os.makedirs(WORK, exist_ok=True)
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); Cc = kv["cell_centres"]
    _, _, piv = qr(U.T, pivoting=True)
    Rm = U.copy(); picks = []; steps = []
    for k in range(5):
        s = np.linalg.norm(Rm, axis=1); i = int(s.argmax())
        steps.append((s / s.max(), i, s.max()))
        picks.append(i); q = Rm[i] / np.linalg.norm(Rm[i]); Rm = Rm - np.outer(Rm @ q, q)
    assert picks == list(piv[:5]), (picks, piv[:5])
    resid = []; Rm = U.copy()
    for i in picks:
        resid.append(Rm.copy()); q = Rm[i] / np.linalg.norm(Rm[i]); Rm = Rm - np.outer(Rm @ q, q)
    rng = np.random.default_rng(0)
    others = rng.choice(np.setdiff1d(np.arange(len(U)), picks), 35, replace=False)
    sample = list(others)
    for pos, c in zip([3, 10, 17, 24, 31], picks):   # 選ぶ列を見やすく間をあけて並べる（並びに意味はない）
        sample.insert(pos, c)
    g, xyz = grid()
    M = build_idw_mapping(Cc, xyz)

    plt.rcParams.update({"font.family": ["Noto Sans CJK JP", "DejaVu Sans"], "font.size": 18})
    text1 = ["各セルに、5 つの型の値を並べた\nベクトル c = (φ1, …, φ5) を持たせる。",
             "色 = ベクトルの長さ |c|。\n型の値が大きいセルほど明るい。"]
    textk = ["すでに選んだ点と「同じ向き」の成分を、\n全セルのベクトルから引く。",
             "色 = 引いた後に残った長さ。\n選んだ点と似ているセルは暗くなる。"]
    frames = []
    for k, (val, i, smax) in enumerate(steps):
        for stage in (0, 1):
            img = render(g, val, picks[:k] + ([i] if stage else []), Cc, new=i if stage else None)
            fig = plt.figure(figsize=(16, 8.4))
            a = fig.add_axes([0.0, 0.02, 0.47, 0.9]); a.imshow(img); a.axis("off")
            cb = fig.add_axes([0.46, 0.18, 0.012, 0.55])
            fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap="viridis"), cax=cb)
            cb.set_ylabel("色 = 各セルの列の長さ（最大で 1）", fontsize=16); cb.tick_params(labelsize=15)
            # 行列 Φᵀ（5 × セル数）の一部の列：残りの成分と、列の長さ
            Rk = resid[k]; cols = list(sample)
            hm = fig.add_axes([0.58, 0.52, 0.40, 0.30])
            vmax = np.abs(U[cols]).max()
            hm.imshow(Rk[cols].T, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
            hm.set_yticks(range(5)); hm.set_yticklabels([f"φ{j + 1}" for j in range(5)], fontsize=15)
            hm.set_xticks([]); hm.set_title("行列 Φᵀ（5 行 × 20,696 列）の一部の列。1 列 = 1 セル", fontsize=16, pad=6)
            for kk, c in enumerate(picks[:k + stage]):
                x = cols.index(c); hm.add_patch(plt.Rectangle((x - .5, -.5), 1, 5, fill=False, ec="#cc0000", lw=2.5))
                pass
            bar = fig.add_axes([0.58, 0.31, 0.40, 0.15])
            nrm = np.linalg.norm(Rk[cols], axis=1) / smax
            bar.bar(range(len(cols)), nrm, color=["#cc0000" if (stage and c == i) else "#4a7bb7" for c in cols], width=0.8)
            bar.set_xlim(-0.5, len(cols) - 0.5); bar.set_ylim(0, 1.05); bar.set_xticks([])
            bar.set_ylabel("列の長さ", fontsize=15); bar.tick_params(labelsize=14)
            for kk, c in enumerate(picks[:k + stage]):
                bar.text(cols.index(c), -0.08, NAMES[kk], ha="center", va="top", color="#cc0000", fontsize=15, weight="bold",
                         transform=bar.get_xaxis_transform())
            t = fig.add_axes([0.58, 0.0, 0.41, 0.26]); t.axis("off")
            t.text(0, 0.86, f"ステップ {k + 1} / 5：" + ("各列（セル）の長さを比べる" if k == 0 else "選んだ列の向きを引いた後の長さ"),
                   fontsize=18, weight="bold", va="top")
            msg = (f"→ 一番長い列を選んで先頭へ（Π）：{NAMES[k]}" if stage else
                   ("型の値が大きいセルほど、列が長い" if k == 0 else "選んだ点と似ている列は短くなる"))
            t.text(0, 0.58, msg, fontsize=17, va="top", color="#cc0000" if stage else "black")
            t.text(0, 0.32, "「一番長い列を選ぶ → その向きを引く」の 5 回が、\n列の入れ替え付き QR 分解 Φᵀ Π = Q R",
                   fontsize=16, va="top", color="0.25", linespacing=1.4)
            fig.suptitle("Q-DEIM：型を見分けるのに一番よい点を 1 つずつ選ぶ", fontsize=24, y=0.995)
            fn = os.path.join(WORK, f"f{2 * k + stage:02d}.png"); fig.savefig(fn, dpi=100, facecolor="white"); plt.close(fig)
            frames.append(Image.open(fn).convert("RGB"))
    # 最後：5 点そろった図
    img = render(g, np.zeros(len(Cc)) + 0.15, picks, Cc)
    fig = plt.figure(figsize=(16, 8.4))
    a = fig.add_axes([0.0, 0.02, 0.5, 0.9]); a.imshow(img); a.axis("off")
    t = fig.add_axes([0.55, 0.02, 0.44, 0.9]); t.axis("off")
    t.text(0, 0.97, "選んだ代表 5 点", fontsize=26, weight="bold", va="top")
    t.text(0, 0.84, "P0 → P1 → P2 → P3 → P4 の順に選ばれた。\n\n型の値が大きく、互いに似ていない点なので、\n5 点の温度から 5 つの型の重みを\nしっかり決められる。",
           fontsize=19, va="top", linespacing=1.5)
    fig.suptitle("Q-DEIM：型を見分けるのに一番よい点を 1 つずつ選ぶ", fontsize=24, y=0.995)
    fn = os.path.join(WORK, "f99.png"); fig.savefig(fn, dpi=100, facecolor="white"); plt.close(fig)
    frames.append(Image.open(fn).convert("RGB"))
    frames[-1].save(os.path.join(IMG, "qdeim_animation_last.png"))
    durs = [2600] * (len(frames) - 1) + [4000]
    pal = frames[3].quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    gif = [f.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG) for f in frames]
    gif[0].save(os.path.join(IMG, "qdeim_animation.gif"), save_all=True, append_images=gif[1:], duration=durs, loop=0)
    import imageio_ffmpeg
    lst = os.path.join(WORK, "list.txt")
    with open(lst, "w") as f:
        for i, d in enumerate(durs):
            f.write(f"file '{os.path.join(WORK, f'f{i:02d}.png' if i < len(durs) - 1 else 'f99.png')}'\nduration {d / 1000}\n")
        f.write(f"file '{os.path.join(WORK, 'f99.png')}'\n")
    subprocess.run([imageio_ffmpeg.get_ffmpeg_exe(), "-y", "-loglevel", "error", "-f", "concat", "-safe", "0", "-i", lst,
                    "-vf", "fps=10,scale=1600:-2", "-pix_fmt", "yuv420p", os.path.join(IMG, "qdeim_animation.mp4")], check=True)
    print("picks", picks, "frames", len(frames))


if __name__ == "__main__":
    main()
