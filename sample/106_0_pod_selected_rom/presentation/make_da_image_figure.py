"""発表スライド「データ同化とは」のイメージ図（逐次型・変分型によらない）：CAE と実験を、ばらつきに応じて混ぜる。

  左：時間とともに。CAE（条件が分からずずれ、ばらつきが広がる）と、ときどきの実験（誤差あり）を混ぜると、
      本当の温度に近づき、ばらつきが狭くなる。
  右：ある時刻での 3 つの分布。CAE（広い）× 実験（誤差あり）→ データ同化（狭く、両者の間）。
    正規分布どうしなら  μ_a = (σ_o² μ_f + σ_f² μ_o)/(σ_f² + σ_o²),  σ_a² = σ_f² σ_o²/(σ_f² + σ_o²)
  計算結果ではなく説明用の模式図（数値は例）。
  python3 presentation/make_da_spread_figure.py → docs/img/da_image.png
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
C_F, C_O, C_A, C_T = '#1f77b4', '#d62728', '#2ca02c', '0.55'


def truth(t):
    return 20 + 5 * (1 - np.exp(-t / 250))


def gauss(x, m, s):
    return np.exp(-0.5 * ((x - m) / s) ** 2) / (s * np.sqrt(2 * np.pi))


def main():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(16, 6.6), gridspec_kw=dict(width_ratios=[1.35, 1]))
    # 左：時間とともに。CAE（条件が分からずずれ、ばらつきが広がる）、実験（ときどき・誤差あり）、データ同化（両方を使う）
    t = np.linspace(0, 600, 601)
    true = 20 + 4.5 * (1 - np.exp(-t / 200)) + 1.1 * np.sin(2 * np.pi * t / 180)   # 運転の入り切りで波打つ
    cae = 20 + 3.2 * (1 - np.exp(-t / 260)) + 0.8 * np.sin(2 * np.pi * t / 180 - 0.7)   # 条件が分からず、大きさも時刻もずれた CAE
    sf = 0.25 + 1.0 * (t / 600) ** 0.8                              # CAE のばらつきは時間とともに広がる
    to = np.arange(40, 601, 40); so = 0.5
    yo = np.interp(to, t, true) + np.random.default_rng(5).normal(0, 0.25, len(to))
    # データ同化：CAE と実験を混ぜた推定（模式。実験の近くでは実験側へ、ばらつきは狭い）
    wgt = sf ** 2 / (sf ** 2 + so ** 2)
    ana = cae + wgt * (true - cae) + 0.06 * np.sin(t / 25)
    sa = np.sqrt(sf ** 2 * so ** 2 / (sf ** 2 + so ** 2))
    ax.plot(t, true, color=C_T, lw=5, alpha=.6, label='本当の温度')
    ax.fill_between(t, cae - 2 * sf, cae + 2 * sf, color=C_F, alpha=.15, lw=0, label='CAE のばらつき（±2σ）')
    ax.plot(t, cae, color=C_F, lw=2.2, label='CAE（解析だけ）')
    ax.fill_between(t, ana - 2 * sa, ana + 2 * sa, color=C_A, alpha=.22, lw=0, label='データ同化のばらつき（±2σ）')
    ax.plot(t, ana, color=C_A, lw=2.8, label='データ同化（CAE＋実験）')
    ax.errorbar(to, yo, yerr=2 * so, fmt='o', color=C_O, ms=8, capsize=4, lw=1.8, label='実験（測定値）と誤差（±2σ）')
    T0 = 440; ax.axvline(T0, color='0.3', ls='--', lw=1.5); ax.text(T0 + 6, 16.4, '右の図の時刻', fontsize=15, color='0.3')
    ax.set_xlim(0, 620); ax.set_ylim(16, 29); ax.set_xlabel('時間 [s]'); ax.set_ylabel('温度 [°C]')
    ax.set_title('CAE だけではずれていく。実験と混ぜると本当の温度に近づく'); ax.grid(alpha=.3)
    ax.legend(loc='upper left', fontsize=14, ncol=2)
    i0 = int(np.argmin(abs(t - T0)))
    MF, SF, MO = float(cae[i0]), float(sf[i0]), float(yo[list(to).index(T0)])
    # 右：1 回の修正での分布
    mf, sf, mo = MF, SF, MO
    ma = (so ** 2 * mf + sf ** 2 * mo) / (sf ** 2 + so ** 2); sa = np.sqrt(sf ** 2 * so ** 2 / (sf ** 2 + so ** 2))
    x = np.linspace(mf - 4 * sf, mo + 3, 600)
    for mu, sg, col, lab in [(mf, sf, C_F, f'CAE（σ = {sf:.2f}）'), (mo, so, C_O, f'実験（σ = {so}）'),
                             (ma, sa, C_A, f'データ同化（σ = {sa:.2f}）')]:
        bx.plot(x, gauss(x, mu, sg), color=col, lw=3, label=lab); bx.fill_between(x, gauss(x, mu, sg), color=col, alpha=.15, lw=0)
        bx.axvline(mu, color=col, ls=':', lw=1.5)
    bx.set_xlabel('温度 [°C]'); bx.set_ylabel('確からしさ'); bx.set_yticks([]); bx.set_title(f'時刻 {T0} s での、確からしさの分布')
    bx.legend(loc='upper left', fontsize=14); bx.grid(alpha=.3); bx.set_ylim(0, gauss(ma, ma, sa) * 1.7)
    bx.text(0.98, 0.97, 'データ同化の平均は、\nばらつきの小さい\n実験の側に寄り、\n幅は両方より狭い', transform=bx.transAxes,
            ha='right', va='top', fontsize=15, color=C_A)
    fig.text(0.5, -0.02, '模式図（計算結果ではない。数値は例）', ha='center', fontsize=15, color='0.3')
    fig.tight_layout()
    fig.savefig(ROOT / 'docs' / 'img' / 'da_image.png', dpi=100, bbox_inches='tight'); plt.close(fig)


if __name__ == '__main__':
    main()
