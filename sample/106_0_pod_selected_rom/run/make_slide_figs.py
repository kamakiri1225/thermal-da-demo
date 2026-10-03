"""発表スライド用に、文字を大きくした「1枚＝1メッセージ」の図を作る.

ブログ用の図は3パネル並べたものが多く、会場のスクリーンでは軸ラベルが読めない。
ここでは計算済みの results/*.json から、スライドに載せる図だけを作り直す。
計算はしない（数値はブログの表と同じもの）。

出力: docs/img/slide_*.png
再現: OMP_NUM_THREADS=4 python3 run/make_slide_figs.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
plt.rcParams.update({"font.size":18,"axes.titlesize":21,"axes.labelsize":19,
                     "xtick.labelsize":17,"ytick.labelsize":17,"legend.fontsize":16})
FS=(13.5,6.0)
NAVY="#0B2545"; RED="#C0392B"; GREEN="#1F9D62"; GRAY="#9AA5B1"; ORANGE="#E67E22"
BLUE="#2E6FD8"; PURPLE="#8E6FB0"; LGRAY="#BBBBBB"
J=lambda n: json.load(open(os.path.join(RES,n),encoding="utf-8"))


def barh(labels,vals,cols,xlabel,title,fname,note=None,hi=None,xmax=None):
    fig,ax=plt.subplots(figsize=FS)
    y=np.arange(len(vals))
    b=ax.barh(y,vals,color=cols,height=.68)
    if hi is not None:
        b[hi].set_edgecolor("#7a1f14"); b[hi].set_linewidth(2.5)
    for i,v in enumerate(vals):
        ax.text(v+max(vals)*0.012,i,f"{v:.3f}",va="center",fontsize=18,
                fontweight="bold" if i==hi else "normal")
    ax.set_yticks(y); ax.set_yticklabels(labels); ax.invert_yaxis()
    ax.set_xlabel(xlabel); ax.set_title(title,pad=14)
    ax.set_xlim(0,(xmax or max(vals)*1.18)); ax.grid(axis="x",alpha=.3)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    if note: fig.text(.99,.015,note,ha="right",fontsize=14,color="#54646f")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,fname),dpi=150); plt.close(fig)
    print("wrote",fname)


def main():
    # ---- ① ROM 真値の holdout（反り A−O）----
    d=J("disp_elsewhere_predict_AO.json")["configs"]
    keys=["変位なし（温度2点のみ）","A・O 自身を観測（循環・参考）","B・C：中段 z=75 mm の Uz",
          [k for k in d if k.startswith("選定2点")][0],"底面近く z=5 mm の Uz（悪い例）"]
    lab=["変位なし\n（温度2点のみ）","A・O 自身を測る\n（参考）","A・O の真下\n（B・C）",
         "選んだ別の2点\n（本研究）","底面近く\n（悪い例）"]
    barh(lab,[d[k]["AO_um"] for k in keys],[GRAY,LGRAY,ORANGE,RED,PURPLE],
         "観測していない反り A−O の誤差 [µm]",
         "A・O を一度も測らずに、反りを当てられるか",
         "slide_holdout_rom.png",note="真値＝ROM の双子実験・5 seed 平均",hi=3)

    # ---- ② 実ソルバ真値・3条件 ----
    c=J("holdout_AO_realsolver.json")["cases"]
    rows=["変位なし（温度2点のみ）","A・O 自身を観測（循環・参考）","B・C（中段 z=75 mm）","選定2点（上面）"]
    lab2=["変位なし","A・O 自身を測る","A・O の真下（B・C）","選んだ別の2点"]
    cases=[("learned","15 W"),("q25","25 W"),("intermittent","間欠加熱")]
    fig,ax=plt.subplots(figsize=FS)
    w=.26; y=np.arange(len(rows))
    for i,(ck,cl) in enumerate(cases):
        v=[c[ck][r]["AO_um"] for r in rows]
        ax.barh(y+(i-1)*w,v,height=w*.92,label=cl,
                color=[NAVY,BLUE,"#79a7e8"][i])
        for j,x in enumerate(v): ax.text(x+.02,y[j]+(i-1)*w,f"{x:.3f}",va="center",fontsize=14)
    ax.set_yticks(y); ax.set_yticklabels(lab2); ax.invert_yaxis()
    ax.set_xlabel("観測していない反り A−O の誤差 [µm]")
    ax.set_title("真値を OpenFOAM＋FrontISTR にしても、3条件とも同じ結論",pad=14)
    ax.legend(title="運転条件",loc="lower right"); ax.grid(axis="x",alpha=.3)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    fig.text(.99,.015,"5 seed 平均。A・O は一度も観測しない",ha="right",fontsize=14,color="#54646f")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"slide_holdout_real.png"),dpi=150); plt.close(fig)
    print("wrote slide_holdout_real.png")

    # ---- ③ 配置法の比較（15 W）----
    m=J("placement_method_comparison.json")
    al={"熱感度 |w| 最大（素朴）":"よく動く点に置く","熱感度の大きさ最大（素朴）":"よく動く点に置く"}
    mm={al.get(k,k):v for k,v in m["methods"].items()}
    order=[("変位なし（温度2点のみ）","変位なし",GRAY),
           ("__rand__","ランダム（中央値）",LGRAY),
           ("よく動く点に置く","よく動く点に置く",PURPLE),
           ("Q-DEIM（QR列ピボット）","Q-DEIM",BLUE),
           ("D最適（温度の det 最小）","D最適（温度場を狙う）",ORANGE),
           ("A最適（温度の trace 最小）","A最適（温度場を狙う）",GREEN),
           ("目的指向（本研究）","目的指向（反りを狙う）",RED)]
    vals=[m["random"]["learned"]["median"] if k=="__rand__" else mm[k]["AO_um"]["learned"] for k,_,_ in order]
    barh([l for _,l,_ in order],vals,[c for _,_,c in order],
         "観測していない反り A−O の誤差 [µm]（15 W）",
         "置き場所の決め方を、同じ土俵で比べる",
         "slide_placement.png",note="候補 14,754 通り・真値＝OpenFOAM＋FrontISTR・5 seed 平均",hi=6)

    # ---- ④ 従来の回帰補正との比較 ----
    b=J("baseline_regression_vs_da.json")["cases"]
    im=J("improve_heating_schedule.json")["cases"]
    labs=["15 W\n（係数を決めた条件）","25 W\n（決めていない条件）","間欠加熱\n（ON/OFF を与える）"]
    reg=[b["learned"]["P2+P4"]["regression_all"],b["q25"]["P2+P4"]["regression_all"],
         b["intermittent"]["P2+P4"]["regression_all"]]
    da =[b["learned"]["P2+P4"]["da_tempdisp_all"],b["q25"]["P2+P4"]["da_tempdisp_all"],
         im["intermittent"]["方法1：ON/OFF時刻を与える"]["AO_um"]]
    fig,ax=plt.subplots(figsize=FS)
    x=np.arange(3); w=.34
    ax.bar(x-w/2,reg,w,label="従来の回帰補正",color=GRAY)
    ax.bar(x+w/2,da ,w,label="データ同化（温度2点＋変位2点）",color=RED)
    for i in range(3):
        ax.text(x[i]-w/2,reg[i]+.008,f"{reg[i]:.2f}",ha="center",fontsize=17)
        ax.text(x[i]+w/2,da[i] +.008,f"{da[i]:.2f}",ha="center",fontsize=17,fontweight="bold")
        ax.annotate(f"回帰の 約 1/{reg[i]/da[i]:.1f}",xy=(x[i],max(reg[i],da[i])+.075),ha="center",
                    fontsize=18,color=RED,fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(labs); ax.set_ylabel("反り A−O の誤差 [µm]")
    ax.set_ylim(0,max(reg)*1.62); ax.legend(loc="upper right",framealpha=.95); ax.grid(axis="y",alpha=.3)
    ax.set_title("学習していない運転条件でも、回帰補正の約 1/3",pad=14)
    for s in ("top","right"): ax.spines[s].set_visible(False)
    fig.text(.99,.015,"温度センサ P2＋P4、同じ 0.3 K のノイズ、同じ 5 seed。真値＝OpenFOAM＋FrontISTR",
             ha="right",fontsize=14,color="#54646f")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"slide_regression.png"),dpi=150); plt.close(fig)
    print("wrote slide_regression.png")

    # ---- ⑤ 同じ場所に温度計か変位計か ----
    s=J("disp_vs_temp_same_place.json")["configs"]
    ks=list(s); lab5=["温度2点\nのみ","＋A・O に\n温度計2本","＋A・O に\n変位計2本","＋両方\n（6本）"]
    fig,axs=plt.subplots(1,2,figsize=(13.5,6.0))
    for ax,key,ttl,un in [(axs[0],"field_K","温度場全体","K"),(axs[1],"AO_um","反り A−O","µm")]:
        v=[s[k][key] for k in ks]; e=[s[k].get(key+"_sd",0) for k in ks]
        cols=[GRAY,ORANGE,RED,BLUE]
        ax.bar(range(4),v,color=cols,yerr=e,capsize=6,error_kw=dict(lw=1.8))
        for i,x in enumerate(v): ax.text(i,x+e[i]+max(v)*.04,f"{x:.3f}",ha="center",fontsize=16)
        ax.set_xticks(range(4)); ax.set_xticklabels(lab5,fontsize=14)
        ax.set_ylabel(f"誤差 [{un}]"); ax.set_title(ttl,pad=10)
        ax.set_ylim(0,(max(np.array(v)+np.array(e)))*1.22); ax.grid(axis="y",alpha=.3)
        for sp in ("top","right"): ax.spines[sp].set_visible(False)
    fig.suptitle("同じ場所に置くなら、温度計と変位計のどちらが良いか",fontsize=21)
    fig.text(.99,.015,"エラーバーは 5 seed のばらつき。真値＝ROM の双子実験",ha="right",fontsize=14,color="#54646f")
    fig.tight_layout(rect=(0,.03,1,.95)); fig.savefig(os.path.join(IMG,"slide_same_place.png"),dpi=150); plt.close(fig)
    print("wrote slide_same_place.png")


    # ---- ⑥ 変位が温度を直す（blog006 の4パネルから t=30 s の1枚だけ切り出す）----
    from PIL import Image
    im=Image.open(os.path.join(IMG,"blog006_disp_fixes_temp.png"))
    W,H=im.size
    im.crop((0,int(H*0.50),int(W*0.505),H)).save(os.path.join(IMG,"slide_disp_fixes_temp.png"))
    print("wrote slide_disp_fixes_temp.png")


if __name__=="__main__": main()
