"""Q-DEIMを説明する16:9スライド画像。模式図のみで、解析結果は変更しない。"""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dacore import plots as _plots
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

IMG=ROOT/"docs/img"
BLUE="#1455bd"; INK="#132c48"; RED="#b8382d"; GREEN="#16704e"
plt.rcParams.update({"font.size":32,"mathtext.fontset":"stix"})

def base(title):
    fig=plt.figure(figsize=(16,9),dpi=100,facecolor="white")
    ax=fig.add_axes([0,0,1,1]); ax.set(xlim=(0,1),ylim=(0,1)); ax.axis("off")
    ax.add_patch(plt.Rectangle((0,.865),1,.135,color=BLUE))
    ax.text(.035,.935,title,color="#fff3cc",fontsize=36,weight="bold",va="center")
    return fig,ax

def box(ax,x,y,w,h,color="#eff4fc",edge="#b6cbe4"):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.009,rounding_size=0.012",
                              facecolor=color,edgecolor=edge,lw=2))

def text(ax,x,y,s,size=32,**kw):
    ax.text(x,y,s,fontsize=size,ha="center",va="center",color=kw.pop("color",INK),**kw)

def save(fig,name):
    fig.savefig(IMG/name,facecolor="white")
    plt.close(fig)
    print(IMG/name)


def matrix(ax, xs, ys, rows, left, right, label, label_x, size=28, column_colors=None, row_colors=None):
    """行列を成分・括弧・寸法ラベルで描く。数式フォントはmathtext。"""
    for row_index, (yy, row) in enumerate(zip(ys, rows)):
        for col_index, (xx, value) in enumerate(zip(xs, row)):
            color = (row_colors[row_index] if row_colors else
                     column_colors[col_index] if column_colors else INK)
            text(ax,xx,yy,"$"+value+"$",size,color=color)
    top=max(ys)+.022; bottom=min(ys)-.022
    for x,sign in ((left,1),(right,-1)):
        ax.plot([x+sign*.007,x,x,x+sign*.007],
                [top,top,bottom,bottom],color=INK,lw=2)
    if label:
        text(ax,label_x,.574,"$"+label+"$",25)


def main():
    fig,ax=base("① 5点の温度から、5つの係数を求める")
    text(ax,.5,.815,"各セルの温度 ＝ そのセルの平均温度 ＋ 5モードの寄与",30)
    ys=[.75,.71,.67,.63]
    matrix(ax,[.085],ys,[[r"T[1]"],[r"T[2]"],[r"\vdots"],[r"T[N]"]],
           .03,.14,r"\mathbf{T}:N\times1",.085)
    text(ax,.18,.69,r"$\approx$",34)
    matrix(ax,[.275],ys,[[r"\overline{T}[1]"],[r"\overline{T}[2]"],
           [r"\vdots"],[r"\overline{T}[N]"]],.215,.335,
           r"\overline{\mathbf{T}}:N\times1",.275)
    text(ax,.37,.69,r"$+$",34)

    for x, color in ((.46,BLUE),(.56,"#c46212"),(.77,GREEN)):
        ax.add_patch(plt.Rectangle((x-.043,.602),.086,.172,
                                  facecolor=color,alpha=.09,edgecolor="none"))
    matrix(ax,[.46,.56,.665,.77],ys,[
        [r"\varphi_1[1]",r"\varphi_2[1]",r"\cdots",r"\varphi_5[1]"],
        [r"\varphi_1[2]",r"\varphi_2[2]",r"\cdots",r"\varphi_5[2]"],
        [r"\vdots",r"\vdots",r"\ddots",r"\vdots"],
        [r"\varphi_1[N]",r"\varphi_2[N]",r"\cdots",r"\varphi_5[N]"]],
        .405,.825,r"\Phi:N\times5",.615,column_colors=[BLUE,"#c46212",INK,GREEN])
    matrix(ax,[.93],[.77,.73,.69,.65,.61],
           [[rf"a_{j}"] for j in range(1,6)],.90,.96,
           r"\mathbf{a}:5\times1",.93,row_colors=[BLUE,"#c46212",INK,INK,GREEN])
    cards=[(.035,"代表5点の温度",r"$T[P0],\ldots,T[P4]$","ROMで計算する"),
           (.375,"5本の式を解く",r"$a_1,\ldots,a_5$","モードの係数が分かる"),
           (.715,"全温度場を復元",r"$\overline{\mathbf{T}}+\Phi\mathbf{a}$","全Nセルへ戻す")]
    for x,title,eq,sub in cards:
        box(ax,x,.315,.25,.205)
        text(ax,x+.125,.475,title,30,weight="bold")
        text(ax,x+.125,.407,eq,32)
        text(ax,x+.125,.347,sub,26)
    for x in (.30,.64):
        ax.add_patch(FancyArrowPatch((x,.42),(x+.055,.42),arrowstyle="-|>",mutation_scale=35,lw=3,color=BLUE))
    text(ax,.5,.262,"Φの1列が1モードのベクトル。同じ色の係数を掛けて足す。",29)
    text(ax,.5,.193,r"$T[P0]"
         r"=\overline{T}[P0]+\varphi_1[P0]a_1+\cdots+\varphi_5[P0]a_5$",34)
    box(ax,.035,.035,.93,.105,color="#fff4df",edge="#e3b55d")
    text(ax,.5,.088,"ただし、似た式ばかりになる5点では、係数を見分けにくい。",30,weight="bold")
    save(fig,"presen_qdeim_why_five.png")


    # Q-DEIMが比較する「点ごとのモード値のベクトル」を示す。
    fig,ax=base("選ぶのは、モードへの反応が異なる点")
    text(ax,.5,.807,"矢印1本 ＝ ある点でのモード値の組。P0とP1を比べる。",30)
    text(ax,.25,.741,"似た向きの2点",31,color=RED,weight="bold")
    text(ax,.75,.741,"異なる向きの2点",31,color=GREEN,weight="bold")
    for rect,good in [([.085,.435,.35,.245],False),([.585,.435,.35,.245],True)]:
        ap=fig.add_axes(rect)
        ap.set(xlim=(-.06,1.05),ylim=(-.06,1.05),aspect="equal")
        ap.axis("off")
        for endpoint in [(1.02,0),(0,1.02)]:
            ap.annotate("",endpoint,(0,0),arrowprops=dict(arrowstyle="->",color="#9aa8b8",lw=1.8))
        ap.text(1.02,-.11,r"$\varphi_1[P]$",ha="right",fontsize=26,color=INK)
        ap.text(-.1,1.01,r"$\varphi_2[P]$",ha="right",fontsize=26,color=INK)
        if good:
            values=[((.85,.22),BLUE,r"$\mathbf{c}_{P0}$",(.86,.10)),
                    ((.30,.85),"#c46212",r"$\mathbf{c}_{P1}$",(.43,.88))]
        else:
            # 長い矢印を先に描き、比例する短い矢印を上へ重ねる。
            values=[((.90,.50),"#c46212",r"$\mathbf{c}_{P1}$",(.92,.64)),
                    ((.54,.30),BLUE,r"$\mathbf{c}_{P0}$",(.59,.16))]
        for endpoint,color,label,label_xy in values:
            ap.annotate("",endpoint,(0,0),
                        arrowprops=dict(arrowstyle="-|>",color=color,lw=4.5,mutation_scale=24))
            ap.text(*label_xy,label,color=color,fontsize=28,ha="center",va="center")
    text(ax,.25,.383,"同じ係数の組合せに反応する",28,color=RED)
    text(ax,.75,.383,"別の係数の組合せに反応する",28,color=GREEN)
    text(ax,.5,.305,
         r"$\mathbf{c}_{Pj}=[\,\varphi_1[Pj],\ \varphi_2[Pj]\,]^{\mathsf{T}}$",32)
    text(ax,.5,.242,r"$T[P0]=\overline{T}[P0]+\varphi_1[P0]a_1+\varphi_2[P0]a_2$",32)
    text(ax,.5,.18,r"$T[P1]=\overline{T}[P1]+\varphi_1[P1]a_1+\varphi_2[P1]a_2$",32)
    box(ax,.035,.045,.93,.085,color="#fff4df",edge="#e3b55d")
    text(ax,.5,.087,"$a_1,a_2$ は両点に共通。係数に掛かるモード値の比が異なる点を選ぶ。",28,weight="bold")
    text(ax,.5,.015,"2モードでの模式例。実際は各点の5モードの値を比較する。",24)
    save(fig,"presen_qdeim_vector_view.png")

    fig,ax=base("② Q-DEIMは、異なる情報を持つ点を選ぶ")
    text(ax,.5,.80,"まず2モード・2点に減らして、位置を選ぶ理由を見る。",30)
    for x,color,edge in ((.035,"#fff0ee","#d89088"),(.52,"#edf8f1","#84b69b")):
        box(ax,x,.36,.445,.355,color,edge)
    text(ax,.257,.655,"似た情報の2点",32,color=RED,weight="bold")
    text(ax,.742,.655,"違う情報の2点",32,color=GREEN,weight="bold")
    text(ax,.257,.567,r"$T[P0]=\overline{T}[P0]+b_1a_1+b_2a_2$",32)
    text(ax,.257,.492,r"$T[P1]=\overline{T}[P1]+\lambda(b_1a_1+b_2a_2)$",32)
    text(ax,.257,.407,"同じ式が2本 → 区別できない",26,color=RED)
    text(ax,.742,.567,r"$T[P0]=\overline{T}[P0]+b_1a_1+b_2a_2$",32)
    text(ax,.742,.492,r"$T[P1]=\overline{T}[P1]+c_1a_1+c_2a_2$",32)
    text(ax,.742,.407,"異なる式 → 係数を見分けられる",26,color=GREEN)
    text(ax,.5,.307,"モード値の比が異なり、ほぼ比例する組合せも避けたい。",27)
    for x,label in ((.035,"1  モード値が大きい点"),(.36,"2  重複する成分を除く"),(.685,"3  残りが大きい点")):
        box(ax,x,.16,.28,.09)
        text(ax,x+.14,.205,label,27,weight="bold")
    text(ax,.5,.102,"この選択を5点まで繰り返す処理が、列ピボット付きQR分解。",29)
    text(ax,.5,.045,"選ぶのはROMの代表点。観測センサの配置は、後半で別に考える。",26,color="#586a80")
    save(fig,"presen_qdeim_why_positions.png")

if __name__=="__main__": main()
