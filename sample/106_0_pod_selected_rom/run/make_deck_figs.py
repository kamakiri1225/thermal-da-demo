"""スライド（blog_008 の内容）用の図を、大きな文字で作る.

出力（presentation/deck_img/）
  deck_map.png        センサと評価点の位置（真上から）
  deck_noda.png       データ同化なし／温度だけ／＋変位2点 の時刻歴（15 W、実ソルバ真値）
  deck_placement.png  置き場所の決め方の比較（反り A−O の誤差、15 W）
入力: results/unobserved_traj_v_learned.npz, blog_008 §7 の表の値
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Circle, FancyArrow

RES = os.path.join(ROOT, "results"); OUT = os.path.join(ROOT, "presentation", "deck_img")
INK, MUTED, BLUE, ORANGE, RED, GOLD = "#14213D", "#5B6475", "#2E6FD8", "#E07A1F", "#C0392B", "#D9A400"
plt.rcParams.update({"font.size": 22, "axes.titlesize": 26, "axes.labelsize": 22, "xtick.labelsize": 19, "ytick.labelsize": 19})


def fig_map():
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.add_patch(Wedge((0, 0), 37.5, 0, 360, width=17.5, fc="#E6EBF2", ec="#7A8699", lw=2))
    ax.add_patch(Wedge((0, 0), 41.5, -60, 60, width=3, fc="#F05A28", ec="none"))
    ax.text(45, 0, "ヒータ", color="#F05A28", fontsize=26, va="center", weight="bold")
    for x, y, lab, dx in [(36.5, 6.9, "温度 P2", 2), (21.3, -1.6, "温度 P0", -2)]:
        ax.add_patch(Circle((x, y), 2.2, fc=RED, ec="white", lw=2, zorder=5))
    ax.text(17, 9, "温度計 P2・P0", color=RED, fontsize=24, weight="bold", ha="right")
    for x, y in [(28.7, 0), (-28.7, 0)]:
        ax.plot(x, y, marker="*", ms=34, color=INK, mec="white", zorder=6)
    ax.text(28.7, -7.5, "A", fontsize=28, weight="bold", ha="center", color=INK)
    ax.text(-28.7, -7.5, "O", fontsize=28, weight="bold", ha="center", color=INK)
    ax.text(0, -4, "★ A・O は測らない\n（知りたい反り A−O）", fontsize=21, ha="center", va="center", color=INK)
    ax.add_patch(plt.Rectangle((9.7 - 2.5, -36.2 - 2.5), 5, 5, fc=BLUE, ec="white", lw=2, zorder=5))
    ax.add_patch(FancyArrow(9.7, -36.2, 8, 0, width=0.9, head_width=3, head_length=3, color=BLUE, zorder=5))
    ax.text(11, -46, "変位計 S1（横方向 Ux）", color=BLUE, fontsize=23, weight="bold", ha="center")
    ax.add_patch(plt.Rectangle((-34.6 - 2.5, 14.4 - 2.5), 5, 5, fc=BLUE, ec="white", lw=2, zorder=5))
    ax.text(-34.6, 22, "変位計 S2\n（上下方向 Uz）", color=BLUE, fontsize=23, weight="bold", ha="center")
    ax.set_xlim(-50, 58); ax.set_ylim(-52, 52); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("真上から見た図（変位計・A・O は上面、温度計は中段 z ≈ 51 mm）", fontsize=21, color=MUTED)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "deck_map.png"), dpi=110, facecolor="#FBFBF8"); plt.close(fig)


def fig_noda():
    z = np.load(os.path.join(RES, "unobserved_traj_v_learned.npz"))
    t = z["t"]; names = list(z["names"]); tl = list(z["tlabels"]); ul = list(z["ulabels"])
    pick = {"データ同化なし（ROM だけ）": (names.index([n for n in names if n.startswith("データ同化なし")][0]), GOLD, "-"),
            "温度計 P2・P0 だけ": (names.index([n for n in names if n.startswith("温度 P2")][0]), BLUE, "--"),
            "＋変位計 S1・S2": (names.index([n for n in names if n.startswith("選んだ")][0]), RED, "-")}
    panels = [("T", tl.index("P4（代表点）"), "観測していない温度（P4）", "温度 [°C]"),
              ("U", ul.index("Uz(O)"), "観測していない変位 Uz(O)", "変位 [µm]"),
              ("U", ul.index("Uz(A)"), "観測していない変位 Uz(A)", "変位 [µm]"),
              ("U", ul.index("反り A−O"), "観測していない反り A−O", "変位 [µm]")]
    fig, axs = plt.subplots(1, 4, figsize=(24, 6.6))
    for ax, (kind, i, title, ylab) in zip(axs, panels):
        ax.axvspan(0, 300, color="#FDEBD0", alpha=.6)
        tv = z["Ttrue"][:, i] - 273.15 if kind == "T" else z["Utrue"][:, i]
        ax.plot(t, tv, color="#111111", lw=8, alpha=.28, label="真値（OpenFOAM＋FrontISTR）")
        for lab, (k, col, ls) in pick.items():
            v = z["T"][k][:, i] - 273.15 if kind == "T" else z["U"][k][:, i]
            ax.plot(t, v, ls, color=col, lw=3.6, label=lab)
        ax.set_title(title); ax.set_xlabel("時刻 [s]"); ax.set_ylabel(ylab); ax.grid(alpha=.3)
        ax.set_xticks([0, 200, 400, 600])
    h, l = axs[0].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=4, fontsize=22, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.tight_layout(rect=(0, 0, 1, 0.88)); fig.savefig(os.path.join(OUT, "deck_noda.png"), dpi=90, facecolor="#FBFBF8")
    plt.close(fig)


def fig_placement():
    rows = [("よく動く点の上位2つ（間隔 5 mm）", 0.564, "#9AA5B1"), ("ランダム（20 通りの中央値）", 0.315, "#9AA5B1"),
            ("A・O 自身を測る（参考）", 0.137, "#9AA5B1"), ("Q-DEIM", 0.137, "#9AA5B1"),
            ("D 最適", 0.108, "#9AA5B1"), ("A 最適", 0.090, "#9AA5B1"), ("本研究（σ で選ぶ）", 0.075, RED)]
    fig, ax = plt.subplots(figsize=(15, 7.2))
    y = np.arange(len(rows))[::-1]
    ax.barh(y, [r[1] for r in rows], color=[r[2] for r in rows], height=0.62)
    for yi, (lab, v, c) in zip(y, rows):
        ax.text(v + 0.008, yi, f"{v:.3f} µm", va="center", fontsize=23, color=INK if c != RED else RED, weight="bold")
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=23)
    ax.set_xlabel("観測していない反り A−O の誤差 [µm]（小さいほど良い）"); ax.set_xlim(0, 0.68)
    ax.spines[["top", "right"]].set_visible(False); ax.grid(axis="x", alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "deck_placement.png"), dpi=100, facecolor="#FBFBF8"); plt.close(fig)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    fig_map(); fig_noda(); fig_placement(); print("wrote deck figures")
