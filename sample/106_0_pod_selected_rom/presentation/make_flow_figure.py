"""発表スライド「全体の流れ」の図：5 つの段階を矢印でつなぎ、各段階に既存の結果の図を小さく入れる。

  python3 presentation/make_flow_figure.py   → docs/img/overall_flow.png
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from PIL import Image

ROOT = Path(__file__).resolve().parent.parent
IMG = ROOT / 'docs' / 'img'
E104 = ROOT.parent / '104_0_openfoam_frontistr_da_enkf' / 'docs' / 'img'

plt.rcParams.update({'font.family': ['Noto Sans CJK JP', 'DejaVu Sans'], 'axes.unicode_minus': False})


def thumb(path, crop=None, frame=None):
    im = Image.open(path)
    if frame is not None:
        im.seek(frame)
    im = im.convert('RGB')
    return im.crop(crop) if crop else im


STEPS = [
    dict(no='①', title='正解を作る', body='OpenFOAM（熱連成）\n→ FrontISTR（変形）',
         img=thumb(IMG / 'cht_to_fem_anim.gif', frame=14, crop=(60, 95, 1100, 500)), group=0),
    dict(no='②', title='実ソルバで同化', body='OI・EnKF\n当たるが約 13 時間',
         img=thumb(E104 / 'openfoam_fem_enkf_Q.png'), group=0),
    dict(no='③', title='軽いモデル（ROM）', body='代表 5 点の ROM\n600 秒が 31 ms',
         img=thumb(IMG / 'qdeim_points.png', crop=(0, 60, 771, 1048)), group=0),
    dict(no='④', title='変位も観測に', body='温度 2＋変位 2 で同化\n反りの誤差 0.076 µm',
         img=thumb(IMG / 'blog008_obs_points.png', crop=(0, 130, 950, 1340)), group=1),
    dict(no='⑤', title='置き場所を選ぶ', body='反りの不確かさが\n最小の点を測る前に',
         img=thumb(IMG / 'sensor_score_idw_oblique.png', crop=(0, 230, 853, 960)), group=1),
]
COL = [('#e8eef7', '#3b6db3'), ('#fdeedd', '#d9661f')]


def main():
    fig = plt.figure(figsize=(16, 7.6))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 16); ax.set_ylim(0, 7.6); ax.axis('off')
    w, gap, x0 = 2.9, 0.25, 0.25
    for g, (label, i0, i1) in enumerate([('道具をそろえる', 0, 2), ('本研究の中心', 3, 4)]):
        xa, xb = x0 + i0 * (w + gap), x0 + i1 * (w + gap) + w
        ax.add_patch(FancyBboxPatch((xa - 0.1, 0.95), xb - xa + 0.2, 6.15, boxstyle='round,pad=0.02,rounding_size=0.15',
                                    fc='none', ec=COL[g][1], lw=2.5, ls='--'))
        ax.text((xa + xb) / 2, 7.1, label, ha='center', va='center', fontsize=22, color=COL[g][1], weight='bold',
                bbox=dict(fc='white', ec='none', pad=4))
    for i, s in enumerate(STEPS):
        x = x0 + i * (w + gap); fc, ec = COL[s['group']]
        ax.add_patch(FancyBboxPatch((x, 1.1), w, 5.85, boxstyle='round,pad=0.02,rounding_size=0.12', fc=fc, ec=ec, lw=2))
        ax.text(x + w / 2, 6.55, f"{s['no']} {s['title']}", ha='center', va='center', fontsize=20, weight='bold')
        im = s['img']; r = im.height / im.width
        bw = w - 0.2; bh = min(3.3, bw * r); bw = bh / r
        iax = fig.add_axes([(x + (w - bw) / 2) / 16, (2.55 + (3.3 - bh) / 2) / 7.6, bw / 16, bh / 7.6])
        iax.imshow(im); iax.axis('off')
        ax.text(x + w / 2, 1.85, s['body'], ha='center', va='center', fontsize=16, linespacing=1.5)
        if i < len(STEPS) - 1:
            ax.add_patch(FancyArrowPatch((x + w + 0.03, 4.0), (x + w + gap - 0.03, 4.0), arrowstyle='-|>',
                                         mutation_scale=28, color='0.25', lw=2.5))
    ax.text(8, 0.45, 'すべてオープンソース：OpenFOAM ＋ FrontISTR ＋ Python（真値は計算。実機は未検証）',
            ha='center', va='center', fontsize=18)
    fig.savefig(IMG / 'overall_flow.png', dpi=100, facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    main()
