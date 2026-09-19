"""計算力学用の中心図：推定対象ごとに最適な『センサ種類』『配置・領域』が変わる（まとめ表）.

すべて本リポジトリの実計算（ROM＋EnKF/OI, 5seed平均）に基づく。
  全温度場  : 温度センサを空間に散らす（感度最大=近接冗長は2.6倍悪化）
  変位差    : 変位センサを高感度W点へ追加（0.151→0.077µm）
  熱量Q    : 場所はほぼ不問（空間差1.6倍）、時間＝加熱期が担う（冷却のみで誤差7倍）
  放熱h    : 温度では応答がノイズ以下で困難。冷却期で相対改善（7.6→5.3mW/K）

出力: docs/img/placement_summary_table.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_placement_summary.py
"""
from __future__ import annotations
import os, sys
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img")

ROWS=[
 ("全温度場\n[K]",
  "温度",
  "空間：散らす\n（冗長を避ける＝Q-DEIMの原理）",
  "散らす 0.35 K vs 感度最大(近接2点) 0.91 K\n→ 冗長配置は2.6倍悪化",
  "#dff3e3"),
 ("未観測の変位差\n[µm]",
  "温度＋変位\n（種類を足す）",
  "空間：熱感度 W（FrontISTR）の高い点",
  "温度2点 0.151 → ＋高W変位2点 0.077 µm（2.0倍改善）\n低W点だと 0.104 µm（効きが半減）",
  "#dceafb"),
 ("発熱量 Q\n[W]",
  "温度\n（場所はほぼ不問）",
  "時間：加熱期の観測が担う",
  "加熱期のみ 0.57 W ≈ 全期間 0.58 W\n冷却期のみ 4.23 W（7倍悪化）。空間差は1.6倍のみ",
  "#fdeeda"),
 ("放熱係数 h\n[mW/K]",
  "温度では困難\n（応答0.05K＜ノイズ0.30K）",
  "時間：冷却期で相対改善\n（∂T/∂h は冷却期に最大）",
  "加熱のみ 7.6 → 全期間 5.3 mW/K\n※依然難しい＝観測の工夫（長時間・別種類）が今後の課題",
  "#fbe3e3"),
]


def main():
    fig,ax=plt.subplots(figsize=(13.6,6.4)); ax.axis("off")
    cols=["推定したい量","最適なセンサ種類","最適な配置・領域","実データの証拠（ROM＋EnKF, 5seed平均）"]
    x=[0.0,0.155,0.34,0.60,1.0]
    ax.text(0.5,1.06,"最適なセンサは『位置』だけでなく『種類』と『領域(空間/時間)』で変わる",
            ha="center",fontsize=15,weight="bold",transform=ax.transAxes)
    ax.text(0.5,1.005,"推定対象ごとのまとめ（すべて本研究の実計算に基づく）",
            ha="center",fontsize=11,color="dimgray",transform=ax.transAxes)
    yh=0.92/ (len(ROWS)+0.6)
    # header
    for j in range(4):
        ax.add_patch(plt.Rectangle((x[j],0.92-yh*0.6),x[j+1]-x[j],yh*0.6,
                     transform=ax.transAxes,facecolor="#33475b",edgecolor="white",lw=1.5))
        ax.text((x[j]+x[j+1])/2,0.92-yh*0.3,cols[j],ha="center",va="center",
                fontsize=11.5,color="white",weight="bold",transform=ax.transAxes)
    y=0.92-yh*0.6
    for row in ROWS:
        target,kind,dom,ev,color=row
        y-=yh
        cells=[target,kind,dom,ev]
        for j,c in enumerate(cells):
            fc = color if j>0 else "#f2f5f8"
            ax.add_patch(plt.Rectangle((x[j],y),x[j+1]-x[j],yh,
                         transform=ax.transAxes,facecolor=fc,edgecolor="white",lw=1.5))
            ax.text((x[j]+x[j+1])/2,y+yh/2,c,ha="center",va="center",
                    fontsize=10.6,transform=ax.transAxes,
                    weight=("bold" if j in (1,2) else "normal"))
    out=os.path.join(IMG,"placement_summary_table.png")
    fig.savefig(out,dpi=150,bbox_inches="tight"); plt.close(fig)
    print("wrote",out)


if __name__=="__main__": main()
