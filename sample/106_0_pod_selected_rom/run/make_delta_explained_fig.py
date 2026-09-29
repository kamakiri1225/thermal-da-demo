"""観測の価値 Δ を「分散が何%減るか」として説明し直す図.

Δ(X,y)=Cov(X,y)^2/(Var(y)+r) は「観測 y を1回使うと、推定したい X の
分散がどれだけ減るか」。生の Δ は X の分散と同じ単位なので行をまたいで
比べられない。そこで「事前分散の何%が消えるか」に直して読む。

① 分散が減る絵（事前→事後）と式の出どころ
② 分子・分母がそれぞれ何を意味するか
③ 推定対象×観測候補の表（%表記）

出力: docs/img/delta_explained.png, results/delta_explained.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_delta_explained_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")

B=np.array([[1.36,0.12,1.2],[0.12,1.04,0.4],[1.2,0.4,4.0]])
w=np.array([0.5,0.2,0.0]); rT=0.09; ru=0.01
cands=[("温度点1",np.array([1.,0,0]),rT),("温度点2",np.array([0,1.,0]),rT),("変位 u",w,ru)]
targets=[("$T_1$（点1の温度）",np.array([1.,0,0]),B[0,0]),
         ("$T_2$（点2の温度）",np.array([0,1.,0]),B[1,1]),
         ("発熱量 $Q$",np.array([0,0,1.]),B[2,2]),
         ("全温度場（$T_1$と$T_2$）",None,B[0,0]+B[1,1])]


def delta(cvec,h,r):
    var=h@B@h+r
    if cvec is None:
        return sum((np.eye(3)[i]@B@h)**2 for i in (0,1))/var
    return (cvec@B@h)**2/var


def main():
    M=np.array([[delta(cv,h,r) for _,h,r in cands] for _,cv,_ in targets])
    prior=np.array([p for _,_,p in targets])
    pct=M/prior[:,None]*100

    fig=plt.figure(figsize=(19.0,9.0))
    gs=fig.add_gridspec(2,3,height_ratios=[1,1.30],hspace=.34,wspace=.28)

    # ① 分散が減る絵
    ax=fig.add_subplot(gs[0,0])
    x=np.linspace(-5,5,400)
    v0=B[0,0]; d=M[0,0]; v1=v0-d
    ax.plot(x,np.exp(-x**2/(2*v0))/np.sqrt(2*np.pi*v0),lw=3,color="#7a8899",
            label=f"観測前  Var = {v0:.2f}")
    ax.fill_between(x,np.exp(-x**2/(2*v0))/np.sqrt(2*np.pi*v0),color="#7a8899",alpha=.16)
    ax.plot(x,np.exp(-x**2/(2*v1))/np.sqrt(2*np.pi*v1),lw=3,color="#C0392B",
            label=f"観測後  Var = {v1:.3f}")
    ax.fill_between(x,np.exp(-x**2/(2*v1))/np.sqrt(2*np.pi*v1),color="#C0392B",alpha=.16)
    ax.annotate("", xy=(1.1,0.30), xytext=(2.6,0.30),
                arrowprops=dict(arrowstyle="<->", color="#1F4E9C", lw=2.2))
    ax.text(1.85,0.34,f"$\\Delta$ = {d:.2f}\nが減った分",ha="center",fontsize=12,
            weight="bold",color="#1F4E9C")
    ax.set_xlabel("推定したい量 $X$ の誤差", fontsize=12)
    ax.set_ylabel("確からしさ", fontsize=12)
    ax.set_title("① $\\Delta$ とは「分散がどれだけ減るか」\n"
                 "（例：$T_1$ を温度点1で測ったとき）", fontsize=13, weight="bold")
    ax.legend(fontsize=11); ax.grid(alpha=.25)

    # ② 式の出どころ
    ax=fig.add_subplot(gs[0,1]); ax.axis("off")
    ax.set_title("② 式はカルマンフィルタの更新そのもの", fontsize=13, weight="bold")
    ax.text(0.02,0.80,"観測 $y$（ノイズ分散 $r$）を1回使うと、",fontsize=12.5,transform=ax.transAxes)
    ax.text(0.04,0.60,r"$\mathrm{Var}(X\,|\,y)=\mathrm{Var}(X)-"
                      r"\frac{\mathrm{Cov}(X,y)^2}{\mathrm{Var}(y)+r}$",
            fontsize=16,transform=ax.transAxes)
    ax.text(0.02,0.42,"この引かれる項が観測の価値：",fontsize=12.5,transform=ax.transAxes)
    ax.text(0.04,0.20,r"$\Delta(X,y)=\dfrac{\mathrm{Cov}(X,y)^2}{\mathrm{Var}(y)+r}$",
            fontsize=18,transform=ax.transAxes,color="#C0392B")
    ax.text(0.02,0.04,"（カルマンゲイン $K=\\mathrm{Cov}(X,y)/(\\mathrm{Var}(y)+r)$ を使えば "
                      "$\\Delta=K\\,\\mathrm{Cov}(X,y)$ ）",
            fontsize=11,transform=ax.transAxes,color="#556")

    # ③ 分子・分母の意味
    ax=fig.add_subplot(gs[0,2]); ax.axis("off")
    ax.set_title("③ 分子と分母が言っていること", fontsize=13, weight="bold")
    ax.text(0.02,0.86,"分子 $\\mathrm{Cov}(X,y)^2$",fontsize=13,weight="bold",
            color="#C0392B",transform=ax.transAxes)
    ax.text(0.04,0.50,"「測るもの $y$」と「知りたいもの $X$」が\n"
                      "連動しているほど価値が高い。\n"
                      "→ $X$ が変われば分子が変わる\n"
                      "＝ 推定対象ごとに最適な観測が違う",
            fontsize=12,transform=ax.transAxes)
    ax.text(0.02,0.32,"分母 $\\mathrm{Var}(y)+r$",fontsize=13,weight="bold",
            color="#1F4E9C",transform=ax.transAxes)
    ax.text(0.04,0.10,"その観測自体がばらつくほど、\n"
                      "またセンサのノイズ $r$ が大きいほど価値が下がる",
            fontsize=12,transform=ax.transAxes)

    # ④ 表（%表記）
    ax=fig.add_subplot(gs[1,:])
    ax.imshow(pct,cmap="RdYlGn",vmin=0,vmax=100,aspect="auto")
    ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in cands],fontsize=14,weight="bold")
    ax.set_yticks(range(4))
    ax.set_yticklabels([f"{t[0]}\n（観測前の分散 {t[2]:.2f}）" for t in targets],fontsize=12)
    for i in range(4):
        b=int(np.argmax(M[i]))
        for j in range(3):
            ax.text(j,i,f"{pct[i,j]:.0f} %\n"+("★この対象にはこれ" if j==b else f"($\\Delta$={M[i,j]:.3f})"),
                    ha="center",va="center",fontsize=13,
                    weight=("bold" if j==b else "normal"),
                    color=("#111" if 15<pct[i,j]<85 else "#111"))
            if j==b: ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,ec="#12308f",lw=4))
    ax.set_title("④ 推定したい量（行）× 観測の候補（列）：観測前の分散のうち何%が消えるか\n"
                 "★＝その行でいちばん価値が高い観測。★の列が行ごとに違うのがポイント",
                 fontsize=13.5,weight="bold")
    ax.set_xlabel("観測の候補（どこを何で測るか）", fontsize=12.5)

    fig.suptitle("観測の価値 $\\Delta$ ― 「その観測で、知りたい量の不確かさが何%消えるか」",
                 fontsize=16,weight="bold")
    fig.text(0.5,0.012,
             "※ 生の $\\Delta$ は対象 $X$ の分散と同じ単位なので、行をまたいで数値を比べても意味がない。"
             "必ず「同じ行の中で列どうし」を比べる（だから%に直してある）。3変数の小モデル "
             "$x=(T_1,T_2,Q)$、変位は $u=0.5T_1+0.2T_2$、$r_T=0.09$・$r_u=0.01$。",
             ha="center",fontsize=11.5,color="#444")
    fig.tight_layout(rect=[0,0.028,1,0.935])
    out=os.path.join(IMG,"delta_explained.png")
    fig.savefig(out,dpi=125); plt.close(fig)
    with open(os.path.join(RES,"delta_explained.json"),"w") as f:
        json.dump(dict(targets=[t[0] for t in targets],cands=[c[0] for c in cands],
                       delta=M.tolist(),prior_var=prior.tolist(),
                       pct_reduction=pct.tolist()),f,ensure_ascii=False,indent=2)
    print("Δ=\n",np.round(M,4),"\n%減=\n",np.round(pct,1))
    print("wrote",out)


if __name__=="__main__":
    main()
