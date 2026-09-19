"""blog_002用の図解：変位観測が温度を直す仕組み（射影と逆流）＋高W/低Wの信号とノイズ.

左：前向き＝温度誤差が重みw経由で変位計に射影 / 逆向き＝ハズレがゲインK=Bw/(..)で各温度へ逆流
右：信号 g·σa とノイズ √r の比較（高Wは信号がノイズを超える、低Wは埋もれる）
出力: docs/img/disp_to_temp_mechanism.png
"""
from __future__ import annotations
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
IMG=os.path.join(ROOT,"docs","img")

def box(ax,x,y,w,h,txt,fc,fs=12):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.02,rounding_size=0.03",
                 fc=fc,ec="#33475b",lw=1.8))
    ax.text(x+w/2,y+h/2,txt,ha="center",va="center",fontsize=fs,weight="bold")

def arrow(ax,p,q,col,lw=2.4,style="-|>",ls="-",rad=0.0):
    ax.add_patch(FancyArrowPatch(p,q,arrowstyle=style,mutation_scale=20,
                 color=col,lw=lw,linestyle=ls,connectionstyle=f"arc3,rad={rad}"))

def left(ax):
    ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
    box(ax,0.4,6.6,2.6,1.5,"温度誤差\n$\\delta T_1$（ヒータ側）","#fde3d6")
    box(ax,0.4,2.2,2.6,1.5,"温度誤差\n$\\delta T_2$","#dceafb")
    box(ax,6.6,4.4,3.0,1.7,"変位計\n$u-u^b=w_1\\delta T_1+w_2\\delta T_2$","#e8f4e0",11)
    arrow(ax,(3.0,7.3),(6.6,5.8),"#33475b",2.8)
    ax.text(4.6,7.15,"重み $w_1$（大）",fontsize=12,color="#33475b",weight="bold")
    arrow(ax,(3.0,3.0),(6.6,4.7),"#33475b",1.6)
    ax.text(4.6,3.15,"重み $w_2$（小）",fontsize=12,color="#33475b",weight="bold")
    # 逆流（同化）
    arrow(ax,(6.7,4.4),(3.0,6.6),"#c0392b",2.8,ls="--",rad=0.25)
    arrow(ax,(6.7,4.3),(3.0,2.4),"#c0392b",1.8,ls="--",rad=-0.25)
    ax.text(4.7,5.35,"補正の逆流\n$K_i \\propto Cov(T_i,u)$",fontsize=12,
            color="#c0392b",weight="bold",ha="center")
    ax.text(5.0,0.9,"前向き（実線）：誤差が変位へ射影 ／ 逆向き（破線）：ハズレが共分散に比例して温度へ配られる",
            fontsize=10.5,ha="center",color="dimgray")
    ax.set_title("① 仕組み：変位は温度誤差の「射影」を測り、共分散で補正が逆流する",fontsize=12.5,weight="bold")

def right(ax):
    ax.set_xlim(0,10); ax.set_ylim(0,10); ax.axis("off")
    # ノイズ帯
    ax.axhspan(0,2.2,xmin=0.06,xmax=0.94,color="#f3c7c7",alpha=.8)
    ax.text(9.2,1.1,"ノイズ $\\sqrt{r}$",fontsize=12,color="#a33",va="center",ha="right",weight="bold")
    # 高W: 大きな信号バー
    ax.add_patch(plt.Rectangle((1.6,0),1.6,7.6,fc="#3b74b8",ec="k"))
    ax.text(2.4,7.9,"高W点\n信号 $g\\,\\sigma_a$（大）",ha="center",fontsize=12,weight="bold",color="#3b74b8")
    # 低W: 小さな信号バー
    ax.add_patch(plt.Rectangle((6.4,0),1.6,1.5,fc="#e58f2a",ec="k"))
    ax.text(7.2,2.6,"低W点\n信号（ノイズに埋没）",ha="center",fontsize=12,weight="bold",color="#b36a12")
    ax.text(5.0,9.4,"$g=w\\cdot\\varphi$（誤差の型と重みの整合）",ha="center",fontsize=12)
    ax.set_title("② 高W/低Wの差＝信号がノイズを超えるか",fontsize=12.5,weight="bold")

def main():
    fig,(a1,a2)=plt.subplots(1,2,figsize=(13.6,5.6),gridspec_kw={"width_ratios":[1.25,1]})
    left(a1); right(a2)
    fig.tight_layout()
    out=os.path.join(IMG,"disp_to_temp_mechanism.png")
    fig.savefig(out,dpi=150,bbox_inches="tight"); plt.close(fig); print("wrote",out)

if __name__=="__main__": main()
