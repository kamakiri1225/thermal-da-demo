"""Q-DEIMの直感を説明する概念図（blog_004 §4用）.

各セルは「モード指紋」(φ1_i, φ2_i)を持つ。r点から場を復元するには、
選点でのモード行列 Φ_P が良条件（行列式が大きい＝可逆）である必要がある。
Q-DEIMは「信号が大きく・互いに独立な点」を貪欲に選ぶ（＝Φ_Pの体積を大きくする）。

左：Q-DEIMの良い選び方（直交気味の2点、行列式大）
右：似た2点（平行、行列式小＝復元が壊れる）

出力: docs/img/qdeim_concept.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_qdeim_concept_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
IMG=os.path.join(ROOT,"docs","img")

# 候補セルの「モード指紋」(mode1, mode2)
CELLS={"a":(0.9,0.15),"b":(0.75,-0.45),"c":(0.55,0.55),"d":(0.35,0.10),
       "e":(-0.15,0.75),"f":(0.15,-0.35),"g":(-0.55,-0.25),"h":(0.20,0.85)}


def panel(ax, pick, color, title, note):
    # 全候補（灰）
    for k,(x,y) in CELLS.items():
        ax.scatter(x,y,s=60,color="lightgray",zorder=2,edgecolor="gray",linewidth=0.6)
        ax.annotate(k,(x,y),textcoords="offset points",xytext=(6,4),fontsize=9,color="gray")
    # 選んだ2点＋原点からのベクトル＋平行四辺形
    p1=np.array(CELLS[pick[0]]); p2=np.array(CELLS[pick[1]])
    for p,k in zip((p1,p2),pick):
        ax.annotate("",xy=p,xytext=(0,0),arrowprops=dict(arrowstyle="-|>",color=color,lw=2.2))
        ax.scatter(*p,s=180,marker="*",color=color,zorder=5,edgecolor="k",linewidth=0.7)
        ax.annotate(k,p,textcoords="offset points",xytext=(7,5),fontsize=12,color=color,weight="bold")
    poly=Polygon([(0,0),p1,p1+p2,p2],closed=True,facecolor=color,alpha=0.16,edgecolor=color,lw=1.2)
    ax.add_patch(poly)
    det=abs(p1[0]*p2[1]-p1[1]*p2[0])
    ax.text(0.5,0.04,f"選点のモード行列式 |det Φ_P| = {det:.3f}",transform=ax.transAxes,
            fontsize=11.5,ha="center",color=color,weight="bold",
            bbox=dict(boxstyle="round,pad=0.25",fc="white",ec=color,alpha=0.9))
    ax.axhline(0,color="k",lw=0.8); ax.axvline(0,color="k",lw=0.8)
    ax.set_xlim(-0.9,1.05); ax.set_ylim(-0.75,1.05)
    ax.set_xlabel("モード1 の値 φ1"); ax.set_ylabel("モード2 の値 φ2")
    ax.set_title(title,fontsize=12.5,weight="bold")
    ax.text(0.5,0.98,note,transform=ax.transAxes,ha="center",va="top",fontsize=10,
            bbox=dict(boxstyle="round,pad=0.3",fc="#fffbe6",ec="#e0c86a"))
    ax.set_aspect("equal")


def main():
    fig,(axL,axR)=plt.subplots(1,2,figsize=(13.6,6.2))
    panel(axL,("a","h"),"tab:green",
          "良い選び方（Q-DEIM）",
          "信号が大きく互いに独立な2点\n→平行四辺形が広い＝復元が安定")
    panel(axR,("a","d"),"tab:red",
          "悪い選び方（似た2点）",
          "モードから見てそっくりな2点\n→平行四辺形がぺちゃんこ＝ノイズ爆発")
    fig.suptitle("Q-DEIMの気持ち：各セルの『モード指紋』(φ1,φ2)から、指紋が最も張り合う点を選ぶ",
                 fontsize=13.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.95])
    out=os.path.join(IMG,"qdeim_concept.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("wrote",out)


if __name__=="__main__": main()
