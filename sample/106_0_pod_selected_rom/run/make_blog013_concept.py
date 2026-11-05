"""クリロフROMの手順を示す模式図。数値解析の分布図ではない。"""
from make_presen_qdeim_explainer import base, text, box, save, BLUE
from matplotlib.patches import FancyArrowPatch

def arrow(ax,start,end):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle="-|>",mutation_scale=32,lw=3,color=BLUE))

def main():
    fig,ax=base("クリロフROM：熱方程式から、分布の基底を作る")
    text(ax,.5,.805,"全N点の温度を、r個の係数で表して計算する",32)
    for x,title,formula,detail in [
        (.035,"① 最初の分布",r"$G\mathbf{m}_0=\mathbf{f}$","定常加熱の応答"),
        (.375,"② 次の分布",r"$G\mathbf{m}_k=C\mathbf{m}_{k-1}$","行列を使って繰り返す"),
        (.715,"③ 基底を作る",r"$V=[\mathbf{v}_1,\ldots,\mathbf{v}_r]$","重複する成分を除く")]:
        box(ax,x,.515,.25,.215)
        text(ax,x+.125,.679,title,30,weight="bold")
        text(ax,x+.125,.610,formula,29)
        text(ax,x+.125,.551,detail,26)
    arrow(ax,(.30,.62),(.355,.62)); arrow(ax,(.64,.62),(.695,.62))
    text(ax,.5,.450,"Vの1列が全N点の分布。係数zを変え、温度の時間変化を表す。",29)
    box(ax,.055,.18,.39,.205)
    box(ax,.555,.18,.39,.205)
    text(ax,.25,.333,"小さな熱方程式を解く",30,weight="bold")
    text(ax,.25,.26,r"$C_r\dot{\mathbf{z}}+G_r\mathbf{z}=\mathbf{f}_rP(t)$",31)
    text(ax,.25,.208,"時間積分する未知量はr個",26)
    arrow(ax,(.465,.28),(.535,.28))
    text(ax,.75,.333,"全節点の温度へ戻す",30,weight="bold")
    text(ax,.75,.26,r"$\mathbf{T}\approx T_{\rm air}\mathbf{1}+V\mathbf{z}$",32)
    text(ax,.75,.208,"元のN点で温度を評価できる",26)
    text(ax,.5,.09,"部分空間 ＝ r本の基底を足し合わせて表せる、温度分布の集合",30,weight="bold")
    text(ax,.5,.033,"基底の作成は事前に行い、加熱履歴P(t)に対する係数z(t)を計算する。",25)
    save(fig,"blog013_krylov_flow.png")

if __name__=="__main__": main()
