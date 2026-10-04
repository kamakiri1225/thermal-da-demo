"""データ同化の1サイクルが「いつ・何を」しているかを、時間軸つきで図解する.

ROM（代表5点）で予報 → 観測点での予測値を作る（gappy-POD／熱感度 W の1行との内積）
→ EnKF で 5点温度・Q・h を補正 → 補正後の値から次の予報、を 30 秒ごとに繰り返す。
画面で読めるよう、文字は大きく・文は短くしている。

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
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.01,rounding_size=0.02",
                                fc=fc,ec=ec,lw=3,transform=ax.transAxes))
    ax.text(x+w/2,y+h-0.06,title,transform=ax.transAxes,ha="center",va="top",fontsize=21,
            fontweight="bold",color=ec)
    ax.text(x+w/2,y+h/2-0.07,body,transform=ax.transAxes,ha="center",va="center",fontsize=18.5,
            color="#1b2430",linespacing=1.5)


def arrow(ax,x0,y0,x1,y1,c=NAVY):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),transform=ax.transAxes,arrowstyle="-|>",
                                 mutation_scale=34,lw=3.5,color=c))


def main():
    fig=plt.figure(figsize=(14,10.5))

    # ---------- ① 600 秒の時間軸 ----------
    ax=fig.add_axes([0.04,0.78,0.92,0.17])
    ax.text(-30,1.75,"① 時間の流れ（0〜600 秒）",fontsize=22,fontweight="bold",color=NAVY,va="bottom")
    ax.axvspan(0,300,ymin=.12,ymax=.62,color="#FDEBD0",alpha=.9)
    ax.text(150,1.15,"ヒータ ON",ha="center",fontsize=17,color="#a0522d",fontweight="bold")
    ax.plot([0,600],[0.5,0.5],color=GRAY,lw=3)
    for t in np.arange(30,601,30):
        ax.plot([t,t],[0.15,0.85],color=RED,lw=3.2)
    ax.text(0,-0.25,"0",ha="center",va="top",fontsize=17)
    for t in (300,600): ax.text(t,-0.25,f"{t} 秒",ha="center",va="top",fontsize=17)
    ax.text(450,1.15,"赤線＝30 秒ごとに補正（20 回）",ha="center",fontsize=17,color=RED,fontweight="bold")
    ax.set_xlim(-30,630); ax.set_ylim(-0.9,1.75); ax.axis("off")

    # ---------- ② 1 サイクル ----------
    ax2=fig.add_axes([0.02,0.30,0.96,0.43]); ax2.axis("off")
    ax2.text(0.01,1.02,"② 30 秒ごとの1サイクル（例：30 → 60 秒）",transform=ax2.transAxes,
             fontsize=22,fontweight="bold",color=NAVY)
    Y=0.20; H=0.72; W=0.225
    box(ax2,0.010,Y,W,H,"(a) 予報","ROM で\n30 秒進める\n\n（5 点の温度\nだけ計算）","#eef3fb",BLUE)
    box(ax2,0.258,Y,W,H,"(b) 読みを予想","このメンバーが\n正しいなら、\nセンサは何を\n示すはずか\n\n5 点の温度に\n重みを掛けて足す","#f2f8f4",GREEN)
    box(ax2,0.506,Y,W,H,"(c) 補正","予想と実測の\nずれを、\n5 点温度・Q・h\nの直す量に\n換算して直す\n\n（換算表＝K）","#fff4e6",ORANGE)
    box(ax2,0.754,Y,0.228,H,"(d) 次へ","直した 5 点温度\nと Q・h が、\n次の ROM の\n初期値になる","#fdeeec",RED)
    for x0,x1 in [(0.236,0.256),(0.484,0.504),(0.732,0.752)]:
        arrow(ax2,x0,Y+H/2,x1,Y+H/2)
    yb=0.07
    ax2.plot([0.87,0.87],[Y-0.005,yb],color=RED,lw=3.5,transform=ax2.transAxes)
    ax2.plot([0.87,0.12],[yb,yb],color=RED,lw=3.5,transform=ax2.transAxes)
    arrow(ax2,0.12,yb,0.12,Y-0.005,c=RED)
    ax2.text(0.495,yb-0.10,"20 回くりかえす",transform=ax2.transAxes,ha="center",fontsize=19,
             color=RED,fontweight="bold")

    # ---------- ③ 状態と従属量 ----------
    ax3=fig.add_axes([0.02,0.02,0.96,0.22]); ax3.axis("off")
    ax3.text(0.01,0.92,"③ 直すもの と 計算で出すもの",transform=ax3.transAxes,
             fontsize=22,fontweight="bold",color=NAVY,va="top")
    ax3.add_patch(FancyBboxPatch((0.01,0.05),0.47,0.58,boxstyle="round,pad=0.01",fc="#fff4e6",
                                 ec=ORANGE,lw=2.5,transform=ax3.transAxes))
    ax3.text(0.245,0.50,"EnKF が直すもの",transform=ax3.transAxes,ha="center",fontsize=19,
             fontweight="bold",color=ORANGE)
    ax3.text(0.245,0.22,"代表 5 点の温度・Q・h\n（この 7 つだけ）",transform=ax3.transAxes,ha="center",
             va="center",fontsize=18)
    ax3.add_patch(FancyBboxPatch((0.52,0.05),0.465,0.58,boxstyle="round,pad=0.01",fc="#f2f8f4",
                                 ec=GREEN,lw=2.5,transform=ax3.transAxes))
    ax3.text(0.755,0.50,"5 点温度から計算で出すもの",transform=ax3.transAxes,ha="center",fontsize=19,
             fontweight="bold",color=GREEN)
    ax3.text(0.755,0.22,"代表点以外の温度（gappy-POD）\n変位（熱感度 W）",transform=ax3.transAxes,
             ha="center",va="center",fontsize=18)
    fig.savefig(os.path.join(IMG,"da_cycle_timeline.png"),dpi=130,facecolor="white"); plt.close(fig)
    print("wrote docs/img/da_cycle_timeline.png")


if __name__=="__main__": main()
