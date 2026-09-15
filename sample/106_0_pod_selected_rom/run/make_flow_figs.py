"""blog_002/004 の mermaid フロー図を PNG 化（WordPress貼り付け用）.

GitHubは mermaid をそのまま描画するので .md 側は mermaid のまま残す。
WordPress用HTMLではコードブロックになってしまうので、その2図だけ画像に差し替える。

出力: docs/img/blog_flow_oi.png, docs/img/blog_flow_pod.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_flow_figs.py
"""
from __future__ import annotations
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p        # 日本語フォント設定
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
IMG=os.path.join(ROOT,"docs","img")


def draw_flow(steps, colors, out, title, cyclic=False):
    n=len(steps)
    fig,ax=plt.subplots(figsize=(2.9*n, 3.2))
    ax.set_xlim(0, n); ax.set_ylim(0, 1); ax.axis("off")
    bw=0.86; bh=0.42; y=0.5
    centers=[]
    for i,(txt,c) in enumerate(zip(steps,colors)):
        x=i+0.5; centers.append(x)
        box=FancyBboxPatch((x-bw/2, y-bh/2), bw, bh,
                           boxstyle="round,pad=0.02,rounding_size=0.04",
                           linewidth=1.8, edgecolor="#33475b", facecolor=c, alpha=0.95)
        ax.add_patch(box)
        ax.text(x, y, txt, ha="center", va="center", fontsize=11.5, weight="bold", color="#12212e")
    for i in range(n-1):
        a=FancyArrowPatch((centers[i]+bw/2, y),(centers[i+1]-bw/2, y),
                          arrowstyle="-|>", mutation_scale=22, linewidth=2.2, color="#33475b")
        ax.add_patch(a)
    if cyclic:
        # 最後→最初へ戻る円環矢印（下側を回す）
        a=FancyArrowPatch((centers[-1], y-bh/2),(centers[0], y-bh/2),
                          connectionstyle="arc3,rad=0.32", arrowstyle="-|>",
                          mutation_scale=22, linewidth=2.2, color="#c0392b", linestyle="--")
        ax.add_patch(a)
        ax.text((centers[0]+centers[-1])/2, y-bh/2-0.24, "次サイクルへ繰り返し",
                ha="center", va="center", fontsize=11, color="#c0392b", weight="bold")
    ax.set_title(title, fontsize=13.5, weight="bold", pad=12)
    fig.tight_layout()
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("wrote", out)


def main():
    # blog_002 §6：実ソルバOIの連成ループ（循環）
    draw_flow(
        ["OpenFOAM CHT\n固体20,696セルの温度",
         "セル温度を\nFEM節点へ補間",
         "FrontISTR熱弾性\n温度上昇→熱変位",
         "A/O変位・A−O差\n未観測C/D・C−D差",
         "OIで温度とQを補正\n次の60秒計算へ"],
        ["#dbeafe","#e0e7ff","#fde7d6","#dcfce7","#fee2e2"],
        os.path.join(IMG,"blog_flow_oi.png"),
        "実ソルバOIの連成ループ（OpenFOAM＋FrontISTR＋OI）",
        cyclic=True)
    # blog_004：POD→Q-DEIM→ROM→EnKF/OI（一方向）
    draw_flow(
        ["OpenFOAM CHT\n20,696セルの温度履歴",
         "① POD(SVD)\n場の型(モード)抽出",
         "② Q-DEIM\n代表5点を選定",
         "③ 5点でROM校正\n残差0.014K",
         "④ 軽いROMで\nEnKF/OI (約31ms)"],
        ["#dbeafe","#e0e7ff","#ede9fe","#dcfce7","#fef9c3"],
        os.path.join(IMG,"blog_flow_pod.png"),
        "POD→Q-DEIM→ROM→データ同化 の流れ",
        cyclic=False)


if __name__=="__main__": main()
