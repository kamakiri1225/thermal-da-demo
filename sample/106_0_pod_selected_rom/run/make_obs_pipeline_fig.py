"""観測データをどう作り、どう同化に使うかの手続きを1枚の流れ図にする.

左  : 観測データの作り方（双子実験）。真値にノイズを載せて y を作る
右  : 予報観測 Yf の作り方（メンバーごとに温度→変位を計算）
下  : EnKF の更新式。y と Yf の差が状態を動かす

出力: docs/img/obs_pipeline.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_obs_pipeline_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
IMG=os.path.join(ROOT,"docs","img")
BLUE="#1F4E9C"; RED="#C0392B"; GREEN="#2e7d32"; GRAY="#556"


def box(ax,x,y,w,h,title,body,fc,ec,fs=17):
    ax.add_patch(plt.Rectangle((x,y),w,h,fc=fc,ec=ec,lw=2.6,zorder=2))
    ax.text(x+w/2,y+h-0.012,title,ha="center",va="top",fontsize=fs+2.5,
            weight="bold",color=ec,zorder=3)
    ax.text(x+w/2,y+h-0.050,body,ha="center",va="top",fontsize=fs,
            linespacing=1.55,zorder=3)


def arrow(ax,x0,y0,x1,y1,col="#333",lw=3.0):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",
                 mutation_scale=30,color=col,lw=lw,zorder=1))


def main():
    # 文字を大きくするため、図の大きさに対して文字を大きめにとる（画面・スライドで読める大きさ）
    fig,ax=plt.subplots(figsize=(16,14))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")
    BW=0.44; BH=0.115; LX=0.035; RX=0.525
    ys=[0.795,0.655,0.515,0.355]

    # ── 左：観測データの作り方（双子実験）──
    ax.text(LX+BW/2,0.935,"① 観測データの作り方（双子実験）",ha="center",fontsize=23,weight="bold",color=BLUE)
    box(ax,LX,ys[0],BW,BH,"OpenFOAM（CHT）",r"流体＋固体を解いて真の温度場"+"\n"+r"20,696 セル × 121 時刻","#eef3f9",BLUE)
    box(ax,LX,ys[1],BW,BH,"FEM 節点へ写す（IDW）",r"近い 8 セルの距離の逆数で平均"+"\n"+r"20,696 セル → 5,040 節点","#eef3f9",BLUE)
    box(ax,LX,ys[2],BW,BH,"FrontISTR（熱弾性）",r"$K_s u=f_\mathrm{th}(T)$ を解いて"+"\n"+r"真の変位","#eef3f9",BLUE)
    box(ax,LX,ys[3],BW,BH+0.03,"ノイズを足して「観測」に",r"$y=$ 真値 $+\ \varepsilon$"+"\n"+r"温度 0.30 K、変位 0.30 µm","#fdecea",RED)
    for y0,y1 in [(ys[0],ys[1]+BH),(ys[1],ys[2]+BH),(ys[2],ys[3]+BH+0.03)]:
        arrow(ax,LX+BW/2,y0,LX+BW/2,y1+0.004)

    # ── 右：予報観測 Yf の作り方 ──
    ax.text(RX+BW/2,0.935,"② メンバーごとの「予想」 $Y_f$",ha="center",fontsize=23,weight="bold",color=GREEN)
    box(ax,RX,ys[0],BW,BH,"各メンバーの状態",r"$z=(T_0,\dots,T_4,\ q,\ h)$"+"\n"+r"60 メンバー（変位は入れない）","#eaf6ee",GREEN)
    box(ax,RX,ys[1],BW,BH,"gappy-POD でモード係数へ",r"$a=U_{rP}^{+}\,(T_5-\bar T_P)$"+"\n"+r"5×5 を解くだけ","#eaf6ee",GREEN)
    box(ax,RX,ys[2],BW,BH,"モード係数から変位へ",r"$u=\bar u+D\,a$"+"\n"+r"$D$ は FrontISTR を事前に 6 回","#eaf6ee",GREEN)
    box(ax,RX,ys[3],BW,BH+0.03,"温度と変位を並べる",r"$Y_f=(T_{P2},\ T_{P0},\ u_1,\ u_2)$"+"\n"+r"60 メンバー × 4 個","#eaf6ee",GREEN)
    for y0,y1 in [(ys[0],ys[1]+BH),(ys[1],ys[2]+BH),(ys[2],ys[3]+BH+0.03)]:
        arrow(ax,RX+BW/2,y0,RX+BW/2,y1+0.004)

    # ── 下：EnKF 更新 ──
    arrow(ax,LX+BW/2,ys[3]-0.005,0.40,0.265,col=RED)
    arrow(ax,RX+BW/2,ys[3]-0.005,0.60,0.265,col=GREEN)
    ax.add_patch(plt.Rectangle((0.035,0.075),0.93,0.185,fc="#fff8f0",ec="#e58f2a",lw=3,zorder=2))
    ax.text(0.50,0.245,"③ EnKF で直す（30 秒ごと）",ha="center",va="top",fontsize=23,weight="bold",color="#b3720f",zorder=3)
    ax.text(0.50,0.175,r"$z^{a}=z^{f}+K\,(y-Y_f),\qquad K=C_{zy}\,(C_{yy}+R)^{-1}$",
            ha="center",va="center",fontsize=26,zorder=3)
    ax.text(0.50,0.105,"ずれ $(y-Y_f)$ を、5 点温度・q・h の直す量に換算する",
            ha="center",va="center",fontsize=18,zorder=3,color="#333")
    ax.text(0.50,0.030,"ループの中で FrontISTR は一度も呼ばない（$D$ は事前計算）",
            ha="center",fontsize=17,color=GRAY,weight="bold")

    fig.suptitle("観測データの作り方と、同化での使い方",fontsize=26,weight="bold",y=0.995)
    fig.tight_layout(rect=[0,0,1,0.965])
    o=os.path.join(IMG,"obs_pipeline.png"); fig.savefig(o,dpi=110); plt.close(fig)
    print("wrote",o)


if __name__=="__main__":
    main()
