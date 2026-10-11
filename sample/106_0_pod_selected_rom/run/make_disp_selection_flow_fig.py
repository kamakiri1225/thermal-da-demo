"""変位計の置き場所を選ぶ手順（blog_008 §5-2）を、文字式でフロー図にする.

出力: docs/img/disp_selection_flow.png
数値は使わず、どの量をどこから作り、どう採点して選ぶかだけを示す。
計算そのものは run/disp_elsewhere_predict_AO.py。
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from dacore import plots as _p  # noqa: F401  (日本語フォント設定)
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

IMG = os.path.join(os.path.dirname(__file__), "..", "docs", "img")
NAVY, ORANGE, GREEN, RED, GRAY = "#1f3b5c", "#d9822b", "#2e7d32", "#c0392b", "#6b7280"


def box(ax, x, y, w, h, title, body, color, fs=15, tfs=17):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=0.012",
                                fc="white", ec=color, lw=2.5))
    ax.add_patch(FancyBboxPatch((x, y + h - 0.034), w, 0.034, boxstyle="round,pad=0,rounding_size=0.012",
                                fc=color, ec=color, lw=2.5))
    ax.text(x + w / 2, y + h - 0.017, title, ha="center", va="center", fontsize=tfs, color="white", weight="bold")
    ax.text(x + 0.015, y + h - 0.048, body, ha="left", va="top", fontsize=fs, color="#111", linespacing=1.55)


def arrow(ax, x0, y0, x1, y1, color=NAVY):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                arrowprops=dict(arrowstyle="-|>", lw=3, color=color, mutation_scale=28))


def main():
    fig = plt.figure(figsize=(17, 21))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    ax.text(0.5, 0.982, "変位計の置き場所の選び方（データ同化を回す前に、計算だけで選ぶ）",
            ha="center", va="center", fontsize=23, weight="bold", color=NAVY)

    # ---- 準備（1回だけ） -------------------------------------------------
    ax.text(0.03, 0.947, "準備（1回だけ計算する）", fontsize=19, weight="bold", color=GRAY)
    box(ax, 0.03, 0.725, 0.30, 0.205, "① FrontISTR：重み",
        "候補 $c$（節点×向き）の変位を\n代表5点の温度の一次式で書く：\n"
        "$y_c = \\mathbf{h}_c^{\\mathsf{T}}\\,\\mathbf{T}$\n"
        "$\\mathbf{h}_c=(h_{c,0},\\dots,h_{c,4})$\n"
        "$h_{c,j}$：$P_j$ が 1 K 上がると\n　　$c$ が何 µm 動くか", ORANGE)
    box(ax, 0.35, 0.725, 0.30, 0.205, "② ROM 4000本：ばらつき",
        "初期温度・発熱量・放熱を振って\nROM を観測なしで回す\n"
        "時刻 $t_k=30k$ 秒（$k=1..10$）で\n"
        "$B(t_k)=\\mathrm{Cov}(\\mathbf{T})$（5×5）\n"
        "＝5点の温度が「どれだけ・\n　どう一緒に」ばらつくか", GREEN)
    box(ax, 0.67, 0.725, 0.30, 0.205, "③ 目的とノイズ",
        "知りたい反り：\n$d=U_z(A)-U_z(O)=\\mathbf{w}^{\\mathsf{T}}\\mathbf{T}$\n"
        "（$\\mathbf{w}$ も①から作る）\n"
        "センサのノイズ：\n$R=0.3^2\\,I$\n（温度 0.3 K、変位 0.3 µm）", NAVY)
    for x in (0.18, 0.50, 0.82):
        arrow(ax, x, 0.723, x, 0.697)

    # ---- 採点ループ ---------------------------------------------------------
    ax.add_patch(FancyBboxPatch((0.03, 0.305), 0.94, 0.390, boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc="#f4f7fb", ec=NAVY, lw=2.5, ls="--"))
    ax.text(0.05, 0.678, "採点：候補 $c$ を1つずつ取り出して点数 $J_c$ を付ける（14,754 候補すべて）",
            fontsize=19, weight="bold", color=NAVY, va="center")

    steps = [
        ("1", "観測行列を作る",
         "$H_c=\\left[\\ \\mathbf{e}_{P2}^{\\mathsf{T}}\\ ;\\ \\mathbf{e}_{P0}^{\\mathsf{T}}\\ ;\\ \\mathbf{h}_c^{\\mathsf{T}}\\ \\right]$　（3行×5列）",
         "1・2行目＝温度計 P2・P0（固定）、3行目＝候補 $c$ の重み"),
        ("2", "測った後に残る温度の共分散",
         "$P_c(t_k)=B-BH_c^{\\mathsf{T}}\\left(H_cBH_c^{\\mathsf{T}}+R\\right)^{-1}H_cB$",
         "カルマンフィルタの更新式。測定値 $y$ は入らない（$B$ は時刻 $t_k$ のもの）"),
        ("3", "測った後に残る反りの分からなさ",
         "$\\sigma_c(t_k)=\\sqrt{\\mathbf{w}^{\\mathsf{T}}P_c(t_k)\\,\\mathbf{w}}$",
         "これを $k=1,\\dots,10$（30〜300 秒）の各時刻で別々に計算"),
        ("4", "時刻平均して点数にする",
         "$J_c=\\frac{1}{10}\\sum_{k=1}^{10}\\sigma_c(t_k)$",
         "小さいほど良い候補"),
    ]
    y = 0.625
    for i, (no, title, eq, note) in enumerate(steps):
        ax.add_patch(plt.Circle((0.075, y - 0.012), 0.018, color=NAVY))
        ax.text(0.075, y - 0.012, no, ha="center", va="center", fontsize=18, color="white", weight="bold")
        ax.text(0.11, y, title, fontsize=18, weight="bold", color="#111", va="center")
        ax.text(0.11, y - 0.033, eq, fontsize=19, color="#111", va="center")
        ax.text(0.11, y - 0.062, note, fontsize=14.5, color=GRAY, va="center")
        if i < len(steps) - 1:
            arrow(ax, 0.075, y - 0.033, 0.075, y - 0.063, GRAY)
        y -= 0.083
    ax.text(0.95, 0.322, "次の候補へ（14,754 回くり返す）", ha="right", fontsize=16, color=NAVY, weight="bold")

    # ---- 選ぶ ---------------------------------------------------------------
    arrow(ax, 0.50, 0.303, 0.50, 0.278)
    box(ax, 0.03, 0.170, 0.455, 0.105, "1点目を選ぶ",
        "$S_1=\\arg\\min_c\\ J_c$\n点数が最小の候補 → 上面の $U_x$", RED, fs=17)
    box(ax, 0.515, 0.170, 0.455, 0.105, "2点目を選ぶ",
        "$H_c$ に $S_1$ の行を足して4行にし、\n残りの候補で採点をやり直す → $S_2=\\arg\\min_c\\ J_c$", RED, fs=17)
    arrow(ax, 0.487, 0.215, 0.513, 0.215, RED)

    # ---- 検証 ---------------------------------------------------------------
    arrow(ax, 0.74, 0.168, 0.74, 0.140)
    box(ax, 0.03, 0.030, 0.94, 0.105, "ここで初めてデータ同化（選んだ配置の検証）",
        "温度計 P2・P0 ＋ 変位計 $S_1$・$S_2$ を観測にして EnKF を回す（60メンバー、30 秒ごとに 20 回、5 seed）。\n"
        "A・O は一度も観測せず、反り $d$ の誤差を真値と比べる。上の採点では測定値もデータ同化も使っていない。",
        GRAY, fs=16)

    os.makedirs(IMG, exist_ok=True)
    out = os.path.join(IMG, "disp_selection_flow.png")
    fig.savefig(out, dpi=110, facecolor="white"); plt.close(fig)
    print("wrote", out)


if __name__ == "__main__":
    main()
