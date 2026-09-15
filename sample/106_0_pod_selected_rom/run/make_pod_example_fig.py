"""blog_004 §3 の「3セル×4時刻」最小例が何をイメージしているかの図解.

左：中空円筒の断面イメージ。ヒータ側からA(強)・B(中)・C(弱)。
右：3セルの温度の時間変化。全セル共通の増え方 g_k を、場所ごとの強さ a_i でスケール
    → u_i(t_k)=ū_i + a_i g_k（右のグラフがまさにこれ）。

出力: docs/img/pod_example_geometry.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_pod_example_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrowPatch
IMG=os.path.join(ROOT,"docs","img")

# §3の例そのままの数値
T={"A":[20,22,24,26],"B":[20,21,22,23],"C":[20,20.5,21,21.5]}
COL={"A":"tab:red","B":"tab:orange","C":"tab:blue"}
AMP={"A":"a_A (大)","B":"a_B (中)","C":"a_C (小)"}


def left_panel(ax):
    R,r=1.0,0.42; Rm=(R+r)/2
    ax.add_patch(Circle((0,0),R,facecolor="#cfd8e3",edgecolor="#33475b",lw=2,zorder=1))
    ax.add_patch(Circle((0,0),r,facecolor="white",edgecolor="#33475b",lw=2,zorder=2))
    # ヒータ（右側外周）
    ax.add_patch(Rectangle((R-0.02,-0.28),0.14,0.56,facecolor="#e8420c",edgecolor="k",lw=1,zorder=3))
    ax.annotate("ヒータ\n(片側加熱)",(R+0.12,0),xytext=(R+0.30,0.0),fontsize=10,color="#e8420c",
                weight="bold",va="center",
                arrowprops=dict(arrowstyle="-|>",color="#e8420c"))
    # A(右=ヒータ側), B(上=中間), C(左=反対側)
    pos={"A":(Rm,0),"B":(0,Rm),"C":(-Rm,0)}
    for k,(x,y) in pos.items():
        ax.scatter(x,y,s=320,color=COL[k],edgecolor="k",lw=1.2,zorder=5)
        ax.annotate(k,(x,y),fontsize=15,color="white",weight="bold",ha="center",va="center",zorder=6)
        ax.annotate(AMP[k],(x,y),textcoords="offset points",
                    xytext=(0,-30 if k=="B" else (24 if y>=0 else -24)),
                    fontsize=11,color=COL[k],weight="bold",ha="center")
    ax.text(0,-R-0.28,"ヒータに近い順に強く振れる：a_A > a_B > a_C",
            ha="center",fontsize=10.5,color="dimgray")
    ax.set_xlim(-1.6,1.7); ax.set_ylim(-1.5,1.5); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title("① 場所（中空円筒の断面イメージ）\n場所ごとの「変動の強さ」a_i",fontsize=12.5,weight="bold")


def right_panel(ax):
    tt=[1,2,3,4]
    for k in ["A","B","C"]:
        ax.plot(tt,T[k],"-o",color=COL[k],lw=2.6,ms=8,label=f"セル{k}")
        ax.annotate(k,(tt[-1],T[k][-1]),textcoords="offset points",xytext=(8,0),
                    fontsize=13,color=COL[k],weight="bold",va="center")
    ax.axhline(20,color="gray",ls=":",lw=1)
    ax.text(1.0,20.15,"全セル 20℃ から出発",fontsize=9.5,color="gray")
    ax.set_xticks(tt); ax.set_xticklabels([f"$t_{i}$" for i in tt])
    ax.set_xlim(0.8,4.6); ax.set_ylim(19.5,26.8)
    ax.set_xlabel("時刻"); ax.set_ylabel("温度 [℃]")
    ax.grid(alpha=.3); ax.legend(loc="upper left",fontsize=10)
    ax.set_title("② 時間（全セル共通の増え方 g_k）\n傾きの違い＝a_i の違い",fontsize=12.5,weight="bold")


def main():
    fig,(axL,axR)=plt.subplots(1,2,figsize=(13.2,5.8),gridspec_kw={"width_ratios":[1,1.15]})
    left_panel(axL); right_panel(axR)
    fig.suptitle("最小例のイメージ：温度 = 下地 ū_i ＋ 場所の強さ a_i × 共通の時間変化 g_k",
                 fontsize=13.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.94])
    out=os.path.join(IMG,"pod_example_geometry.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("wrote",out)


if __name__=="__main__": main()
