"""発表「はじめに」の真ん中に置く GIF：立形マシニングセンタの熱変形（ParaView）と、工具とワークのずれ（相対変位）のグラフ。

  計算は 107（sample/107_machining_center_rom）の S1：4 時間運転 → 4 時間停止（OpenFOAM の熱伝導 → FrontISTR の熱変形）。
  工作機械を想定した計算例で、実機ではない。
  python3 presentation/make_intro_machine_gif.py → docs/img/intro_machine_deform.gif、intro_machine_deform_4h.png
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent
M107 = ROOT.parent / '107_machining_center_rom'
sys.path.insert(0, str(M107 / 'common'))
import pvlabel                          # noqa: E402  （107 の ParaView 用ラベル）
from vtu import write_vtu, write_points, XYZ, CONN   # noqa: E402
from scenarios import SOURCES                      # noqa: E402
import numpy as _np


OUT = ROOT / 'docs' / 'img'; WORK = ROOT / 'presentation' / 'intro_machine_work'; WORK.mkdir(parents=True, exist_ok=True)
INFO = json.loads((M107 / 'common' / 'mesh' / 'mesh_info.json').read_text())
TIP, TAB = INFO['tip_node'] - 1, INFO['table_node'] - 1
SCALE = 2000
FONT = '/home/kamakiri/.fonts/NotoSansCJKjp-Bold.otf'
SCEN = 'S1_run4h_stop4h'


JP = {'bearing': '主軸の軸受 150 W', 'motor': '主軸モータ 250 W', 'zfeed': 'Z 軸の送り 40 W', 'yfeed': 'Y 軸の送り 40 W'}


def source_flag():
    """107 の 01_openfoam_thermal/make_figs.py と同じ：発熱源の箱の中の鋳鉄セルに番号を付ける"""
    cen = XYZ[CONN].mean(axis=1)
    part = np.load(M107 / 'common' / 'mesh' / 'mesh.npz')['part']
    cast = np.isin(part, ['bed', 'table', 'column', 'head', 'spindle'])
    flag = np.zeros(len(cen))
    for k, src in enumerate(SOURCES.values(), 1):
        lo, hi = np.array(src['box'][0]), np.array(src['box'][1])
        flag[np.all((cen >= lo) & (cen <= hi), axis=1) & cast] = k
    return flag, cen


def main():
    d = np.load(M107 / '02_frontistr_thermal_deform' / 'data' / f'{SCEN}.npz')
    t, U = d['times'], d['U'].astype(float); h = t / 3600
    ks = list(range(0, len(t), 3))
    flag, cen = source_flag()
    # 発熱源の場所（箱の中のセルの中心）。名前は変形後の位置に付ける
    src_cells = {n: np.where(flag == i)[0] for i, n in enumerate(SOURCES, 1)}
    umax = float(np.linalg.norm(U, axis=2).max() * 1e6); rng = [0, float(np.ceil(umax / 10) * 10)]
    specs = []
    for k in ks:
        write_vtu(WORK / f'U_{k:03d}.vtu', cell={'source': flag}, point={'U': U[k], 'Umag_um': np.linalg.norm(U[k], axis=1) * 1e6})
        write_points(WORK / f'p_{k:03d}.vtp', XYZ[[TIP, TAB]] + SCALE * U[k, [TIP, TAB]])
        specs.append(dict(vtu=str(WORK / f'U_{k:03d}.vtu'), out=str(WORK / f'f_{k:03d}.png'), color='Umag_um',
                          assoc='POINTS', range=rng, cmap='Jet', cbar='変位の大きさ [µm]', fmt='%.0f',
                          ticks=list(np.linspace(0, rng[1], 5)), edges=False, warp=dict(vector='U', scale=SCALE),
                          points=str(WORK / f'p_{k:03d}.vtp'), point_size=0.06, point_color=[0.95, 0.0, 0.6],
                          size=[1500, 1300], threshold=dict(array='source', range=[0.5, 9], opacity=0.5, on_warped=True, color=[0.2, 0.2, 0.2]),
                          hide_bar=True))
    p = WORK / 'specs.json'; p.write_text(json.dumps(specs, ensure_ascii=False))
    if not all((WORK / f'f_{k:03d}.png').exists() for k in ks) or '--render' in sys.argv:
      subprocess.run([sys.executable, str(M107 / 'common' / 'pvpython_launch.py'), str(M107 / 'common' / 'render_pv.py'),
                    str(p)], check=True, stdout=subprocess.DEVNULL)

    rel = (U[:, TIP] - U[:, TAB]) * 1e6
    plt.rcParams.update({'font.family': ['Noto Sans CJK JP', 'DejaVu Sans'], 'axes.unicode_minus': False,
                         'font.size': 24, 'axes.titlesize': 26, 'axes.labelsize': 26, 'xtick.labelsize': 24,
                         'ytick.labelsize': 24, 'legend.fontsize': 23})
    font = ImageFont.truetype(FONT, 52)
    frames = []
    for k in ks:
        fig, ax = plt.subplots(figsize=(10, 11.5))
        ax.axvspan(0, 4, color='#fff2a8', lw=0, label='運転中（主軸・送りが発熱）')
        for i, (n, c) in [(1, ('Y（前後）', '#d95f02')), (2, ('Z（上下）', '#7570b3'))]:
            ax.plot(h, rel[:, i], color=c, lw=4, label=n)
            ax.plot(h[k], rel[k, i], 'o', color=c, ms=16)
        ax.axvline(h[k], color='k', lw=2); ax.axhline(0, color='k', lw=1)
        ax.set_xlim(0, 8); ax.set_xticks(range(0, 9, 2)); ax.set_ylim(-75, 15); ax.grid(alpha=.3)
        ax.set_xlabel('時間 [h]'); ax.set_ylabel('工具とワークのずれ [µm]')
        ax.set_title('工具とワークのずれ（熱変位）')
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=21)
        fig.tight_layout(); fig.savefig(WORK / 'plot.png', dpi=100); plt.close(fig)
        pvlabel.label_png(WORK / f'f_{k:03d}.png', XYZ[[TIP, TAB]] + SCALE * U[k, [TIP, TAB]], ['工具', 'ワーク'],
                          offset=[(110, -150), (110, 110)], size=48, color=(180, 0, 110))
        cs = [cen[c].mean(0) + SCALE * U[k][CONN[c]].reshape(-1, 3).mean(0) for c in src_cells.values()]
        pvlabel.label_png(WORK / f'f_{k:03d}.png', cs, [JP[n] for n in src_cells],
                          offset=[(330, 40), (330, -60), (300, -80), (-470, 60)], size=40, color=(40, 40, 40))
        left = Image.open(WORK / f'f_{k:03d}.png').convert('RGB')
        right = Image.open(WORK / 'plot.png').convert('RGB')
        im = Image.new('RGB', (2550, 1420), 'white')
        im.paste(left, (0, 140)); im.paste(right.resize((1000, 1150)), (1520, 160))
        state = '運転中' if h[k] <= 4 else '停止中'
        dr = ImageDraw.Draw(im)
        dr.text((40, 15), f'立形マシニングセンタの熱変形（変形 {SCALE} 倍）　{h[k]:.1f} 時間（{state}）', font=font, fill='black')
        dr.text((40, 85), '温度：OpenFOAM、熱変形：FrontISTR で計算。色＝変位の大きさ（青 0 → 赤 120 µm）、黒い箱＝発熱源', font=ImageFont.truetype(FONT, 40), fill=(60, 60, 60))
        frames.append(im.resize((1275, 710), Image.LANCZOS))
    k4 = min(ks, key=lambda k: abs(t[k] - 4 * 3600))
    frames[ks.index(k4)].save(OUT / 'intro_machine_deform_4h.png')
    big = Image.new('RGB', (1275, 710 * 8))
    for i, f in enumerate(frames[::max(1, len(frames) // 8)][:8]):
        big.paste(f, (0, 710 * i))
    pal = big.quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    gif = [f.quantize(palette=pal, dither=Image.Dither.FLOYDSTEINBERG) for f in frames]
    gif[0].save(OUT / 'intro_machine_deform.gif', save_all=True, append_images=gif[1:], duration=250, loop=0)
    print('frames', len(frames), 'range', rng)


if __name__ == '__main__':
    main()
