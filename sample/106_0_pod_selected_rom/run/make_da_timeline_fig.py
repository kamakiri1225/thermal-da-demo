"""データ同化の1サイクルが「いつ・何を」しているかを、時間軸つきで図解する.

ROM（代表5点）で予報 → 観測点での予測値を作る（gappy-POD／熱感度 W の1行との内積）
→ EnKF で 5点温度・Q・h を補正 → 補正後の値から次の予報、を 30 秒ごとに繰り返す。

出力: docs/img/da_cycle_timeline.png
再現: python3 run/make_da_timeline_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
IMG=os.path.join(ROOT,"docs","img")
NAVY="#0B2545"; BLUE="#2E6FD8"; RED="#C0392B"; GREEN="#1F9D62"; ORANGE="#E67E22"; GRAY="#8696a7"


def box(ax,x,y,w,h,title,body,fc,ec):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.012,rounding_size=0.018",
                                fc=fc,ec=ec,lw=2.2,transform=ax.transAxes))
    ax.text(x+w/2,y+h-0.035,title,transform=ax.transAxes,ha="center",va="top",fontsize=14.5,
            fontweight="bold",color=ec)
    ax.text(x+w/2,y+h-0.105,body,transform=ax.transAxes,ha="center",va="top",fontsize=11.6,
            color="#1b2430",linespacing=1.55)


def arrow(ax,x0,y0,x1,y1,c=NAVY,rad=0.0):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),transform=ax.transAxes,arrowstyle="-|>",
                                 mutation_scale=22,lw=2.4,color=c,connectionstyle=f"arc3,rad={rad}"))


def main():
    fig=plt.figure(figsize=(17,11.2))
    # ---------- 上段：600 秒の時間軸 ----------
    ax=fig.add_axes([0.05,0.80,0.90,0.13])
    ax.axvspan(0,300,color="#FDEBD0",alpha=.7); ax.text(150,1.32,"ヒータ ON（0〜300 秒）",ha="center",fontsize=12.5,color="#a0522d")
    ax.plot([0,600],[0.5,0.5],color=GRAY,lw=2)
    for t in np.arange(30,601,30):
        ax.plot([t,t],[0.18,0.82],color=RED,lw=2.6)
    ax.plot([0],[0.5],"o",color=NAVY,ms=10)
    ax.text(0,-0.15,"0 秒\n初期値はばらばら\n（60 メンバー）",ha="center",va="top",fontsize=10.5)
    for t in (30,60,90):
        ax.text(t,1.0,f"{t}",ha="center",fontsize=11,color=RED,fontweight="bold")
    ax.text(600,1.0,"600",ha="center",fontsize=11,color=RED,fontweight="bold")
    ax.text(345,1.0,"…　30 秒ごとに観測・補正（全 20 回）　…",ha="center",fontsize=11.5,color=RED)
    ax.annotate("",xy=(60,0.5),xytext=(30,0.5),arrowprops=dict(arrowstyle="<->",color=BLUE,lw=2))
    ax.text(45,-0.15,"この 30 秒を\n下で拡大",ha="center",va="top",fontsize=10.5,color=BLUE,fontweight="bold")
    ax.set_xlim(-25,625); ax.set_ylim(-0.9,1.55); ax.axis("off")
    ax.text(-25,1.55,"① 全体の時間の流れ（600 秒）",fontsize=15,fontweight="bold",color=NAVY,va="bottom")

    # ---------- 中段：1 サイクル（t=30 → 60 秒）----------
    ax2=fig.add_axes([0.03,0.29,0.94,0.45]); ax2.axis("off")
    ax2.text(0.0,1.03,"② 1 サイクルの中身（例：t = 30 秒 → 60 秒）",transform=ax2.transAxes,
             fontsize=15,fontweight="bold",color=NAVY)
    W=0.215; H=0.80; Y=0.14
    box(ax2,0.015,Y,0.205,H,"(a) 予報　30 → 60 秒",
        "ROM で 60 メンバーを\n30 秒進める\n（2 秒刻み × 15 ステップ）\n\n進めるのは\n代表 5 点の温度だけ\n\n各メンバーの Q・h を使う\n観測はまだ使わない",
        "#eef3fb",BLUE)
    box(ax2,0.255,Y,W,H,"(b) 観測点での予測値",
        "各メンバーの 5 点温度から\n「観測点ではいくつか」を作る\n\n温度計（代表点）→ そのまま\n温度計（代表点以外）\n　→ gappy-POD の 1 行との内積\n変位計 → 熱感度 W の 1 行との内積\n\n全セルの復元はしない\n（必要な行だけ。FrontISTR も呼ばない）",
        "#f2f8f4",GREEN)
    box(ax2,0.505,Y,W,H,"(c) 補正（t = 60 秒）",
        "センサの読みと (b) を比べる\n\nEnKF で状態を引き寄せる\n直すのは\n　代表 5 点の温度\n　発熱量 Q\n　放熱係数 h\nの 7 つだけ\n\n60 メンバーの散らばりが\n「どれをどれだけ直すか」を決める",
        "#fff4e6",ORANGE)
    box(ax2,0.755,Y,0.225,H,"(d) 次のサイクルへ",
        "補正後の 5 点温度・Q・h を\n初期値にして (a) に戻る\n（60 → 90 秒の予報）\n\n━━━━━━━━━━━━\n評価用（ループの外）\n補正後の 5 点温度から\n　全 20,696 セル → gappy-POD\n　全 5,040 節点の変位 → W\nで復元して真値と比べる",
        "#fdeeec",RED)
    for x0,x1 in [(0.222,0.253),(0.472,0.503),(0.722,0.753)]:
        arrow(ax2,x0,Y+H/2,x1,Y+H/2)
    # 戻りの矢印：(d) の下 → 左へ → (a) の下（コの字）
    yb=Y-0.07
    ax2.plot([0.872,0.872],[Y-0.005,yb],color=RED,lw=2.6,transform=ax2.transAxes)
    ax2.plot([0.872,0.112],[yb,yb],color=RED,lw=2.6,transform=ax2.transAxes)
    arrow(ax2,0.112,yb,0.112,Y-0.005,c=RED)
    ax2.text(0.49,yb-0.075,"これを 30 秒ごとに 20 回くりかえす（予報 → 予測値 → 補正 → 予報 …）",
             transform=ax2.transAxes,ha="center",fontsize=13,color=RED,fontweight="bold")

    # ---------- 下段：状態と従属量 ----------
    ax3=fig.add_axes([0.03,0.02,0.94,0.20]); ax3.axis("off")
    ax3.text(0.0,1.0,"③ 何が「状態」で、何が「そこから計算するもの」か",transform=ax3.transAxes,
             fontsize=15,fontweight="bold",color=NAVY,va="top")
    ax3.text(0.0,0.70,"状態（EnKF が直すもの）",transform=ax3.transAxes,fontsize=13,fontweight="bold",color=ORANGE)
    ax3.text(0.0,0.44,"$z=(T_0,T_1,T_2,T_3,T_4,\\;Q,\\;h)$　… 代表 5 点の温度と、発熱量・放熱係数の 7 つ",
             transform=ax3.transAxes,fontsize=13)
    ax3.text(0.55,0.70,"そこから計算するもの（直接は直さない）",transform=ax3.transAxes,fontsize=13,fontweight="bold",color=GREEN)
    ax3.text(0.55,0.50,"・代表点以外の温度 ＝ 平均 ＋ (gappy-POD の行)・(5点温度 − 平均)\n・変位 ＝ 平均 ＋ (W の行)・(5点温度 − 平均)",
             transform=ax3.transAxes,fontsize=12.5,va="top",linespacing=1.6)
    ax3.text(0.0,0.02,"→ 代表点以外の温度計や変位計の読みも、5 点温度の一次式で予測できるので、そのまま観測に使える。"
             "補正が入るのは 5 点温度・Q・h で、ほかはそれに従って変わる。",
             transform=ax3.transAxes,fontsize=12.5,color=NAVY,fontweight="bold")
    fig.savefig(os.path.join(IMG,"da_cycle_timeline.png"),dpi=150,facecolor="white"); plt.close(fig)
    print("wrote docs/img/da_cycle_timeline.png")


if __name__=="__main__": main()
