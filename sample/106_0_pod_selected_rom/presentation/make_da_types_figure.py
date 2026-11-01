"""発表スライド「データ同化とは」の模式図：逐次型（予報 → 観測で修正 → 予報…）と変分型（期間全体の軌道を観測に合わせる）。

  計算結果ではなく説明用の模式図（曲線は式で描いた例）。
  python3 presentation/make_da_types_figure.py   → docs/img/da_sequential_vs_variational.png
"""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family': ['Noto Sans CJK JP', 'DejaVu Sans'], 'axes.unicode_minus': False,
                     'font.size': 16, 'axes.titlesize': 20, 'axes.labelsize': 18, 'xtick.labelsize': 16,
                     'ytick.labelsize': 16, 'legend.fontsize': 15})


def truth(t):
    return 1 - np.exp(-t / 3) + 0.15 * np.sin(1.3 * t)


def main():
    rng = np.random.default_rng(1)
    tobs = np.arange(1, 10, 1.5); yobs = truth(tobs) + rng.normal(0, 0.05, len(tobs))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(16, 6.2), sharey=True)
    t = np.linspace(0, 10, 400)
    for a in (ax, bx):
        a.plot(t, truth(t), color='0.6', lw=5, alpha=0.6, label='本当の状態')
        a.plot(tobs, yobs, 'o', color='tab:red', ms=11, label='観測')
        a.set_xlabel('時間'); a.set_xticks([]); a.set_yticks([]); a.set_xlim(0, 10); a.set_ylim(-0.2, 1.5)
    # 逐次型：観測のたびに修正して次を予報（モデルは少しずれている＝伸びが遅い）
    x0, t0 = 0.0, 0.0
    for k, (to, yo) in enumerate(zip(list(tobs) + [10], list(yobs) + [None])):
        tt = np.linspace(t0, to, 50)
        xf = x0 + (truth(tt) - truth(t0)) * 0.6
        ax.plot(tt, xf, color='tab:blue', lw=2.5, label='予報（モデル）' if k == 0 else None)
        if yo is not None:
            xa = xf[-1] + 0.7 * (yo - xf[-1])
            ax.annotate('', xy=(to, xa), xytext=(to, xf[-1]), arrowprops=dict(arrowstyle='-|>', color='tab:green', lw=2.5))
            x0, t0 = xa, to
    ax.plot([], [], color='tab:green', lw=2.5, label='観測で修正')
    ax.set_title('逐次型：予報 → 観測で修正 → 予報 …')
    ax.legend(loc='lower right')
    # 変分型：期間全体の初期値（とパラメータ）を変えて、軌道全体を観測に合わせる
    bx.plot(t, (truth(t) - truth(0)) * 0.6, color='tab:blue', lw=2, ls='--', label='合わせる前の軌道')
    bx.plot(t, truth(t) + 0.03, color='tab:purple', lw=2.8, label='期間全体を合わせた軌道')
    bx.axvspan(0.5, 9.5, color='#fff2a8', alpha=0.6, lw=0, label='同化する期間（窓）')
    bx.set_title('変分型：期間全体を観測に合わせる')
    bx.legend(loc='lower right')
    fig.text(0.5, -0.02, '模式図（計算結果ではない）', ha='center', fontsize=15, color='0.3')
    fig.tight_layout()
    fig.savefig(ROOT / 'docs' / 'img' / 'da_sequential_vs_variational.png', dpi=100, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
