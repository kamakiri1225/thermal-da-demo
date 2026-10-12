"""FrontISTR を6回解くだけで、全節点・全方向の重みがそろう理由（重ね合わせ）の図（blog_008 §5-2）.

出力: docs/img/fem6_superposition.png
"""
from __future__ import annotations
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

IMG = os.path.join(ROOT, "docs", "img")
INK, RED, BLUE, GRAY, GOLD = "#14213D", "#C0392B", "#2E6FD8", "#6B7280", "#B8860B"


def rbox(ax, x, y, w, h, fc, ec, txt, fs=17, col=INK, bold=False):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012", fc=fc, ec=ec, lw=2))
    ax.text(x + w / 2, y + h / 2, txt, ha="center", va="center", fontsize=fs, color=col, weight="bold" if bold else "normal",
            linespacing=1.4)


def arrow(ax, x0, y0, x1, y1, col=INK):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", lw=2.5, color=col, mutation_scale=22))


def main():
    fig = plt.figure(figsize=(20, 11.5))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.965, "FrontISTR を6回解くだけで、全節点・全方向の「重み」がそろう理由 ＝ 重ね合わせ",
            ha="center", fontsize=24, weight="bold", color=INK)

    # 段1：6つの温度の型
    ax.text(0.03, 0.885, "① 温度分布を「6つの型」の足し算で表す（POD）", fontsize=19, weight="bold", color=RED)
    labs = ["平均の\n温度分布\n$\\bar T(x)$", "型1\n$\\varphi_1(x)$\n全体が温まる", "型2\n$\\varphi_2(x)$\nヒータ側⇔反対側",
            "型3\n$\\varphi_3(x)$", "型4\n$\\varphi_4(x)$", "型5\n$\\varphi_5(x)$"]
    xs = [0.04 + i * 0.155 for i in range(6)]
    for x, l in zip(xs, labs):
        rbox(ax, x, 0.70, 0.135, 0.16, "#FDF0EE", RED, l, fs=16)
    ax.text(0.5, 0.655, "実際の温度分布 $T(x)\\;\\approx\\;\\bar T(x)+a_1\\varphi_1(x)+a_2\\varphi_2(x)+\\cdots+a_5\\varphi_5(x)$"
            "　（ $a_k$ ＝各型の強さ）", ha="center", fontsize=19, color=INK)

    # 段2：FrontISTR 6回
    for x in xs:
        arrow(ax, x + 0.0675, 0.635, x + 0.0675, 0.598, GRAY)
    ax.text(0.03, 0.568, "② 型ごとに FrontISTR を1回ずつ解く（計6回）", fontsize=19, weight="bold", color=BLUE)
    ulabs = ["$u[\\bar T]$", "$u[\\varphi_1]$", "$u[\\varphi_2]$", "$u[\\varphi_3]$",
             "$u[\\varphi_4]$", "$u[\\varphi_5]$"]
    for x, l in zip(xs, ulabs):
        rbox(ax, x, 0.415, 0.135, 0.13, "#EDF3FC", BLUE, l + "\n全 5,040 節点 × 3方向\nの変位", fs=15)
    ax.text(0.5, 0.375, "熱弾性は線形なので、変位も同じ足し算になる：　"
            "$u\\;\\approx\\;u[\\bar T]+a_1u[\\varphi_1]+a_2u[\\varphi_2]+\\cdots+a_5u[\\varphi_5]$",
            ha="center", fontsize=19, color=INK)

    # 段3：a を代表5点の温度から
    ax.text(0.03, 0.30, "③ 型の強さ $a_k$ は、代表5点の温度から決まる（gappy-POD）", fontsize=19, weight="bold", color=GOLD)
    ax.text(0.5, 0.24, "代表5点の温度 $T_P=\\bar{T}_P+\\Phi_P\\,a$　→　$a=\\Phi_P^{-1}(T_P-\\bar{T}_P)$"
            "　（ $\\Phi_P$ ＝5つの型の、代表5点での値。5×5）", ha="center", fontsize=19, color=INK)

    # 段4：まとめ
    rbox(ax, 0.06, 0.04, 0.88, 0.14, "#14213D", "#14213D",
         "②と③を合わせると、どの節点 $c$ のどの方向の変位も、代表5点の温度の一次式になる\n"
         "$u_c=u_c[\\bar T]+h_c^{T}(T_P-\\bar{T}_P)$,　　重み $h_c^{T}=(u_c[\\varphi_1],\\dots,u_c[\\varphi_5])\\,\\Phi_P^{-1}$"
         "　→　6回の結果から、15,120 通りすべての重みが掛け算で出る", fs=18, col="#F4F1EA")
    fig.savefig(os.path.join(IMG, "fem6_superposition.png"), dpi=95, facecolor="white"); plt.close(fig)
    print("wrote fem6_superposition.png")


if __name__ == "__main__":
    main()
