"""センサの種類ごとに使う感度が違うことを1枚で示す.

検証（run/why_lowW_wins_check.py）の結論:
  温度センサ → 発熱感度 dT/dQ（温度RMSEとの相関 -0.835）
  変位センサ → 熱感度 W=Ks^-1 H_T の行ノルム（105で高W/低Wが3倍差）
W を温度センサに流用すると逆効果（相関 +0.405）。
「感度が高いところに置く」は正しく、感度の種類がセンサで違うだけ。

出力: docs/img/two_sensitivities.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_two_sensitivities_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")


def main():
    d=json.load(open(os.path.join(RES,"why_lowW_wins.json")))
    W=np.array(d["W"]); dq=np.array(d["dTdQ"]); rT=np.array(d["rmseT"])
    lab=[f"P{i}" for i in range(5)]

    fig=plt.figure(figsize=(17.0,7.0))
    gs=fig.add_gridspec(1,3,width_ratios=[1,1,1.12],wspace=.34)

    ax=fig.add_subplot(gs[0])
    ax.scatter(dq,rT,s=230,c="#1F4E9C",zorder=4)
    for i in range(5):
        ax.annotate(lab[i],(dq[i],rT[i]),textcoords="offset points",
                    xytext=(11,8),fontsize=13,weight="bold")
    z=np.polyfit(dq,rT,1); xx=np.linspace(dq.min()*.95,dq.max()*1.05,20)
    ax.plot(xx,np.polyval(z,xx),"--",color="#1F4E9C",lw=2,alpha=.6)
    ax.set_xlabel("発熱感度 $\\partial T_i/\\partial Q$",fontsize=13)
    ax.set_ylabel("温度RMSE [K]（小さいほど良い）",fontsize=13)
    ax.set_title(f"① 温度センサ ← 発熱感度 $\\partial T/\\partial Q$\n"
                 f"相関 {np.corrcoef(dq,rT)[0,1]:+.3f} ＝ 高いほど良い ✓",
                 fontsize=13.5,weight="bold",color="#2e7d32")
    ax.grid(alpha=.3)

    ax=fig.add_subplot(gs[1])
    ax.scatter(W,rT,s=230,c="#C0392B",zorder=4)
    for i in range(5):
        ax.annotate(lab[i],(W[i],rT[i]),textcoords="offset points",
                    xytext=(11,8),fontsize=13,weight="bold")
    z=np.polyfit(W,rT,1); xx=np.linspace(W.min()*.9,W.max()*1.05,20)
    ax.plot(xx,np.polyval(z,xx),"--",color="#C0392B",lw=2,alpha=.6)
    ax.set_xlabel("熱感度 $\\|W_{i,:}\\|$ [µm/K]",fontsize=13)
    ax.set_ylabel("温度RMSE [K]（小さいほど良い）",fontsize=13)
    ax.set_title(f"② 同じ温度センサを $W$ で選ぶと\n"
                 f"相関 {np.corrcoef(W,rT)[0,1]:+.3f} ＝ 高いほど悪い ✗",
                 fontsize=13.5,weight="bold",color="#C0392B")
    ax.grid(alpha=.3)

    ax=fig.add_subplot(gs[2]); ax.axis("off")
    ax.set_title("③ 使い分け ― どちらも「感度が高い点」で正しい",
                 fontsize=13.5,weight="bold")
    rows=[("温度センサ","発熱感度 $\\partial T/\\partial Q$",
           "「そこは発熱に敏感か」\n＝ 発熱量 $Q$ の情報が乗る点","#1F4E9C"),
          ("変位センサ","熱感度 $W=K_s^{-1}H_T$ の行ノルム",
           "「そこはよく動くか」\n＝ 固定端から遠い自由端","#C0392B")]
    for k,(who,what,why,col) in enumerate(rows):
        y=0.78-k*0.40
        ax.add_patch(plt.Rectangle((0.02,y-0.26),0.96,0.33,transform=ax.transAxes,
                     fc="#f4f7fb",ec=col,lw=2.2,clip_on=False))
        ax.text(0.06,y+0.02,who,transform=ax.transAxes,fontsize=15,weight="bold",color=col)
        ax.text(0.34,y+0.02,what,transform=ax.transAxes,fontsize=14)
        ax.text(0.06,y-0.16,why,transform=ax.transAxes,fontsize=12,color="#333")
    ax.text(0.02,0.06,"※ $W$ が高い点＝よく動く点＝固定端から遠い「場の端」。\n"
                      "　 端は場全体を代表しないので、温度センサには向かない。",
            transform=ax.transAxes,fontsize=12,color="#556")

    fig.suptitle("「感度が高いところにセンサを置く」は正しい ― ただし感度の種類がセンサで違う\n"
                 "5点それぞれを唯一の温度センサにして同化（60メンバー・5 seed平均・加熱期）",
                 fontsize=15.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.83])
    o=os.path.join(IMG,"two_sensitivities.png"); fig.savefig(o,dpi=135); plt.close(fig)
    print("wrote",o)


if __name__=="__main__":
    main()
