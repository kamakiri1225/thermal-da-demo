"""σ の求め方を、5つの世界の数値例で4段階の絵にする（blog_008 §5-2 (b′)）.

① 5つの世界を点で描く ② 直線を引き、外れ（縦の棒）を測る ③ 傾きを変えると外れの二乗平均が変わる
④ 一番小さくなる傾きのときの外れが σ。良い場所 X と、悪い場所 Y を並べる。
出力: docs/img/sigma_steps.png
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt

IMG = os.path.join(ROOT, "docs", "img")
plt.rcParams.update({"font.size": 15})
d = np.array([1.0, 2.0, 3.0, 2.0, 1.0])
S = {"X": (np.array([2.0, 4.0, 6.0, 4.1, 1.9]), "#C0392B", "良い場所 X"),
     "Y": (np.array([5.1, 1.3, 4.8, 3.0, 2.2]), "#2E6FD8", "悪い場所 Y")}


def main():
    dd = d - d.mean()
    fig, axs = plt.subplots(2, 4, figsize=(24, 11))
    for r, (key, (y, col, name)) in enumerate(S.items()):
        dy = y - y.mean(); vd = (dd ** 2).mean(); c = (dd * dy).mean(); vy = (dy ** 2).mean(); a0 = c / vy
        sig = np.sqrt(vd - c * c / vy)
        # ① 点
        ax = axs[r, 0]
        ax.scatter(dy, dd, s=160, color=col, zorder=3)
        offs = [(-10, 14), (-50, 14), (8, 6), (10, -22), (-20, -26)] if key == "X" else [(8, 6)] * 5
        for i in range(5): ax.annotate(f"世界{i+1}", (dy[i], dd[i]), xytext=offs[i], textcoords="offset points", fontsize=13)
        ax.set_title(f"① {name}：5つの世界を点で描く", fontsize=16)
        # ② 直線と外れ（試しに a=0.3）
        ax2 = axs[r, 1]; atry = 0.3 if key == "X" else 0.3
        xs = np.linspace(-2.6, 2.6, 10)
        ax2.scatter(dy, dd, s=160, color=col, zorder=3); ax2.plot(xs, atry * xs, color="#555", lw=2, ls="--")
        for i in range(5): ax2.plot([dy[i], dy[i]], [dd[i], atry * dy[i]], color="#E07A1F", lw=3)
        m_try = ((dd - atry * dy) ** 2).mean()
        ax2.set_title(f"② 試しに傾き {atry} の直線を引く\n外れ（橙の棒）の二乗平均 ＝ {m_try:.3f}", fontsize=16)
        # ③ 傾きを変えたときの外れの二乗平均
        ax3 = axs[r, 2]; aa = np.linspace(-0.4, 1.0, 200); mm = [((dd - a * dy) ** 2).mean() for a in aa]
        ax3.plot(aa, mm, color=col, lw=3); ax3.axvline(atry, color="#555", ls="--", lw=1.5)
        ax3.plot(a0, sig ** 2, "o", ms=14, color="#14213D")
        ax3.annotate(f"一番小さい：傾き {a0:.2f}\n外れの二乗平均 {sig**2:.3f}", (a0, sig ** 2), xytext=(15, 40),
                     textcoords="offset points", fontsize=14, arrowprops=dict(arrowstyle="->"))
        ax3.set_xlabel("直線の傾き a"); ax3.set_ylabel("外れの二乗平均")
        ax3.set_title("③ 傾きを変えて、外れが一番小さい直線を探す", fontsize=16)
        # ④ 最良の直線と σ
        ax4 = axs[r, 3]
        ax4.scatter(dy, dd, s=160, color=col, zorder=3); ax4.plot(xs, a0 * xs, color="#14213D", lw=3)
        ax4.fill_between(xs, a0 * xs - sig, a0 * xs + sig, color="#14213D", alpha=.12)
        for i in range(5): ax4.plot([dy[i], dy[i]], [dd[i], a0 * dy[i]], color="#E07A1F", lw=3)
        ax4.set_title(f"④ 一番良い直線で残る外れの幅\nσ ＝ √{sig**2:.4f} ＝ {sig:.2f}（測る前は {np.sqrt(vd):.2f}）",
                      fontsize=16, color=col)
        for ax in (axs[r, 0], ax2, ax4):
            ax.set_xlabel(f"場所 {key} のセンサの値（平均からのずれ）"); ax.set_ylabel("反り（平均からのずれ）")
            ax.set_xlim(-2.8, 2.8); ax.set_ylim(-1.6, 1.8); ax.grid(alpha=.3); ax.axhline(0, color="#bbb", lw=1); ax.axvline(0, color="#bbb", lw=1)
        ax3.grid(alpha=.3)
    fig.suptitle("σ の求め方：センサの値から反りを直線で予想し、一番うまく予想できたときに残る外れの幅が σ（5つの世界の作り話の数値例）",
                 fontsize=19, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(IMG, "sigma_steps.png"), dpi=85, facecolor="white"); plt.close(fig)
    print("wrote sigma_steps.png")


if __name__ == "__main__":
    main()
