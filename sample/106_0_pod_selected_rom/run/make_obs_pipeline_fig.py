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


def box(ax,x,y,w,h,title,body,fc,ec,fs=10.5):
    ax.add_patch(plt.Rectangle((x,y),w,h,fc=fc,ec=ec,lw=2.0,zorder=2))
    ax.text(x+w/2,y+h-0.014,title,ha="center",va="top",fontsize=fs+1.4,
            weight="bold",color=ec,zorder=3)
    ax.text(x+w/2,y+h-0.048,body,ha="center",va="top",fontsize=fs,
            linespacing=1.70,zorder=3)


def arrow(ax,x0,y0,x1,y1,col="#333",lw=2.2):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle="-|>",
                 mutation_scale=20,color=col,lw=lw,zorder=1))


def main():
    fig,ax=plt.subplots(figsize=(18.4,9.6))
    ax.set_xlim(0,1); ax.set_ylim(0,1); ax.axis("off")

    # ── 左：観測データの作り方（双子実験）──
    ax.text(0.235,0.965,"① 観測データをどう作ったか（双子実験）",
            ha="center",fontsize=15.5,weight="bold",color=BLUE)
    box(ax,0.06,0.828,0.35,0.112,"OpenFOAM（CHT）",
        "流体＋固体を解いて真の温度場\n20,696セル × 121時刻",  "#eef3f9",BLUE)
    arrow(ax,0.235,0.828,0.235,0.800)
    box(ax,0.06,0.686,0.35,0.112,"IDW でFEM節点へ写像",
        "$T_p=\\sum w_c T_c/\\sum w_c,\\ \\ w_c=1/\\|x_p-x_c\\|$\n20,696セル → 5,040節点（$k=8$）",
        "#eef3f9",BLUE)
    arrow(ax,0.235,0.686,0.235,0.658)
    box(ax,0.06,0.544,0.35,0.112,"FrontISTR（熱弾性）",
        "$K_s u=f_\\mathrm{th}(T)$ を解いて真の変位\n上面A/O・中高さB/C の $U_z$",
        "#eef3f9",BLUE)
    arrow(ax,0.235,0.544,0.235,0.510)
    box(ax,0.06,0.368,0.35,0.140,"ノイズを載せて「観測」にする",
        "$y=h(x^{\\mathrm{true}})+\\varepsilon,\\qquad \\varepsilon\\sim N(0,R)$\n"
        "$R=\\mathrm{diag}(\\sigma_T^2,\\sigma_T^2,\\sigma_u^2,\\sigma_u^2)$\n"
        "$\\sigma_T=0.30$ K,  $\\sigma_u=0.30$ µm",
        "#fdecea",RED)
    ax.text(0.235,0.344,"→ 4つの数字（温度2＋変位2）が30秒ごとに届く",
            ha="center",fontsize=11.5,color=RED,weight="bold")

    # ── 右：予報観測 Yf の作り方 ──
    ax.text(0.745,0.965,"② メンバーごとに「予報観測」$Y_f$ を作る",
            ha="center",fontsize=15.5,weight="bold",color=GREEN)
    box(ax,0.57,0.828,0.35,0.112,"各メンバーの状態",
        "$z^{(m)}=(T_1,\\dots,T_5,\\ Q,\\ h)$\n$m=1,\\dots,60$（変位は入っていない）",
        "#eaf6ee",GREEN)
    arrow(ax,0.745,0.828,0.745,0.800)
    box(ax,0.57,0.686,0.35,0.112,"gappy-POD でモード係数へ",
        "$a^{(m)}=U_{r,P}^{+}(T_5^{(m)}-\\bar T_P)$\n5×5 の連立を解くだけ",
        "#eaf6ee",GREEN)
    arrow(ax,0.745,0.686,0.745,0.658)
    box(ax,0.57,0.544,0.35,0.112,"モード係数から変位へ",
        "$u^{(m)}=u_\\mathrm{mean}+D\\,a^{(m)}$\n$D$ は事前にFrontISTRを6回呼んで作成",
        "#eaf6ee",GREEN)
    arrow(ax,0.745,0.544,0.745,0.510)
    box(ax,0.57,0.368,0.35,0.140,"温度と横に並べる",
        "$Y_f^{(m)}=(T_{i_1}^{(m)},T_{i_2}^{(m)},\\ u_1^{(m)},u_2^{(m)})$\n"
        "（60メンバー × 4成分）",
        "#eaf6ee",GREEN)
    ax.text(0.745,0.344,"→ 同化ループ中にFrontISTRは一度も呼ばない",
            ha="center",fontsize=11.5,color=GREEN,weight="bold")

    # ── 下：EnKF 更新 ──
    arrow(ax,0.235,0.332,0.40,0.288,col=RED)
    arrow(ax,0.745,0.332,0.60,0.288,col=GREEN)
    ax.add_patch(plt.Rectangle((0.09,0.105),0.82,0.175,fc="#fff8f0",ec="#e58f2a",lw=2.6,zorder=2))
    ax.text(0.50,0.258,"③ EnKF の更新 ― 温度だけを観測したときと同じ式",
            ha="center",fontsize=15,weight="bold",color="#b3720f",zorder=3)
    ax.text(0.50,0.196,
        r"$z^{a}=z^{b}+K\,(y-Y_f),\qquad "
        r"K=C_{zy}(C_{yy}+R)^{-1},\qquad "
        r"C_{zy}=\frac{dZ^{T}dY}{N-1}$",
        ha="center",fontsize=17,zorder=3)
    ax.text(0.50,0.135,
        "$C_{zy}$ に $\\mathrm{Cov}(T_i,u_j)$ の列が入るので、"
        "変位の残差が温度と発熱量 $Q$ を動かす",
        ha="center",fontsize=12.5,zorder=3,color="#333")

    ax.text(0.50,0.055,
        "直しているのは変位ではなく温度場。温度が合えば $u=W(T-T_\\mathrm{ref})$ で"
        "どの点のどの方向も自動的に従う（結果③・③-2）",
        ha="center",fontsize=12,color=GRAY)
    ax.text(0.50,0.015,
        "実装: dacore/enkf.py ― $Y_f$ を直接受け取る設計なので、非線形な観測演算子にもそのまま使える",
        ha="center",fontsize=11,color=GRAY)

    fig.suptitle("観測データはどう作られ、どう同化に使われるのか ― 数式の手続き",
                 fontsize=17,weight="bold",y=0.995)
    fig.tight_layout(rect=[0,0,1,0.97])
    o=os.path.join(IMG,"obs_pipeline.png"); fig.savefig(o,dpi=125); plt.close(fig)
    print("wrote",o)


if __name__=="__main__":
    main()
