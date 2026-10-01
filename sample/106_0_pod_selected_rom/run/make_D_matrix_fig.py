"""D行列とは何か・なぜ同化ループ中にFEMを呼ばずに済むのかを1枚で説明する.

① 熱弾性が線形なので重ね合わせが効く：u = W T̄ + Σ a_k (W φ_k)
   W φ_k は a に依存しない → 先に計算して保存できる
② D の作り方：FrontISTR を 6回（平均場 1 + モード 5）呼ぶだけ
③ 実際の D 行列（数値）
④ コスト：FrontISTR 1回＝分オーダー vs 行列積 4.4 µs

出力: docs/img/D_matrix_explained.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_D_matrix_fig.py
"""
from __future__ import annotations
import os, sys, time
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
BLUE="#1F4E9C"; RED="#C0392B"; GREEN="#2e7d32"


def main():
    d=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); D=d["D"]; um=d["u_mean"]
    rows=[(2,"A  $U_z$"),(5,"O  $U_z$"),(8,"B  $U_z$"),(11,"C  $U_z$")]
    a=np.random.randn(60,5); t0=time.perf_counter()
    for _ in range(2000): _=um+a@D.T
    us=(time.perf_counter()-t0)/2000*1e6

    fig=plt.figure(figsize=(19.0,10.2))
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.04],hspace=.30,wspace=.22)

    # ① なぜFEMを呼ばずに済むか
    ax=fig.add_subplot(gs[0,0]); ax.axis("off")
    ax.set_title("① なぜ同化ループ中にFEMを呼ばずに済むのか", fontsize=15,
                 weight="bold", color=BLUE, loc="left")
    ax.text(0.01,0.80,"熱弾性は温度について線形：",fontsize=13,transform=ax.transAxes)
    ax.text(0.05,0.645,r"$u = W\,(T-T_\mathrm{ref})$",fontsize=17,transform=ax.transAxes)
    ax.text(0.01,0.525,"温度はPODで「平均＋モードの足し算」：",fontsize=13,transform=ax.transAxes)
    ax.text(0.05,0.375,r"$T=\bar T+\sum_{k=1}^{5} a_k\,\varphi_k$",fontsize=17,transform=ax.transAxes)
    ax.text(0.01,0.255,"代入して、線形性で項ごとに分ける：",fontsize=13,transform=ax.transAxes)
    ax.text(0.03,0.075,r"$u = \underset{u_\mathrm{mean}}{W\bar T} \;+\; "
                        r"\sum_{k=1}^{5} a_k\;\underset{D[:,k]}{(W\varphi_k)}$",
            fontsize=19,transform=ax.transAxes,color=RED)
    ax.text(0.01,-0.075,"→ $W\\varphi_k$ は $a_k$ に依存しない。"
                        "だから先に1回だけ計算して保存しておけばよい",
            fontsize=13,transform=ax.transAxes,color=RED,weight="bold")

    # ② D の作り方
    ax=fig.add_subplot(gs[0,1]); ax.axis("off")
    ax.set_title("② $D$ の作り方 ― FrontISTR を 6回 呼ぶだけ", fontsize=15,
                 weight="bold", color=GREEN, loc="left")
    tb=[("1回目","$T=\\bar T$（平均場）","$u_\\mathrm{mean}=W\\bar T$"),
        ("2回目","$T=T_\\mathrm{ref}+\\varphi_1$","$D[:,1]=W\\varphi_1$"),
        ("3回目","$T=T_\\mathrm{ref}+\\varphi_2$","$D[:,2]=W\\varphi_2$"),
        ("4〜6回目","$T=T_\\mathrm{ref}+\\varphi_3,\\varphi_4,\\varphi_5$","$D[:,3..5]$")]
    ax.text(0.03,0.86,"回",fontsize=12.5,weight="bold",transform=ax.transAxes)
    ax.text(0.20,0.86,"FrontISTR に渡す温度場",fontsize=12.5,weight="bold",transform=ax.transAxes)
    ax.text(0.68,0.86,"得られるもの",fontsize=12.5,weight="bold",transform=ax.transAxes)
    for i,(a1,a2,a3) in enumerate(tb):
        y=0.71-i*0.155
        ax.add_patch(plt.Rectangle((0.01,y-0.055),0.97,0.125,transform=ax.transAxes,
                     fc="#eaf6ee",ec=GREEN,lw=1.5,clip_on=False))
        ax.text(0.03,y,a1,fontsize=12,transform=ax.transAxes,va="center")
        ax.text(0.20,y,a2,fontsize=13,transform=ax.transAxes,va="center")
        ax.text(0.68,y,a3,fontsize=13,transform=ax.transAxes,va="center",color=GREEN)
    ax.text(0.01,0.015,"これで終わり。以後の同化では $D$ を掛けるだけで、"
                       "FrontISTR は一度も呼ばない",
            fontsize=13,transform=ax.transAxes,color=GREEN,weight="bold")

    # ③ 実際の D
    ax=fig.add_subplot(gs[1,0])
    sub=np.array([D[i] for i,_ in rows])*1000      # µm→nm（モード1単位あたり）
    v=np.abs(sub).max()
    ax.imshow(sub,cmap="RdBu_r",vmin=-v,vmax=v,aspect="auto")
    for i in range(len(rows)):
        for k in range(5):
            ax.text(k,i,f"{sub[i,k]:+.2f}",ha="center",va="center",fontsize=12,
                    weight="bold",color=("white" if abs(sub[i,k])>0.66*v else "#1b2430"))
    ax.set_xticks(range(5)); ax.set_xticklabels([f"$\\varphi_{k+1}$" for k in range(5)],fontsize=14)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([l for _,l in rows],fontsize=13)
    ax.set_title("③ 実際の $D$（上面A/Oと中高さB/C の $U_z$ 行）単位 nm\n"
                 "$\\varphi_2$ の列は A が−、O が＋ ＝ 逆向き → 反りを生む",
                 fontsize=13.5,weight="bold")
    ax.set_xlabel("モード（列）",fontsize=12)

    # ④ コスト
    ax=fig.add_subplot(gs[1,1]); ax.axis("off")
    ax.set_title("④ だからこそ速い", fontsize=15, weight="bold", color=RED, loc="left")
    ax.text(0.02,0.80,"同化1サイクルで必要な計算（60メンバー分）",
            fontsize=13,transform=ax.transAxes,weight="bold")
    ax.text(0.05,0.655,r"$a=U_{r,P}^{+}(T_5-\bar T_P)$   … $5\times5$ の連立",
            fontsize=13.5,transform=ax.transAxes)
    ax.text(0.05,0.545,r"$u=u_\mathrm{mean}+D\,a$   … $12\times5$ の行列ベクトル積",
            fontsize=13.5,transform=ax.transAxes)
    ax.text(0.02,0.415,f"実測 {us:.1f} µs／サイクル（60メンバー込み）",
            fontsize=15,transform=ax.transAxes,color=RED,weight="bold")
    ax.add_patch(plt.Rectangle((0.02,0.10),0.95,0.235,transform=ax.transAxes,
                 fc="#fdecea",ec=RED,lw=2.2,clip_on=False))
    ax.text(0.06,0.275,"もし毎サイクルFrontISTRを呼んだら（見積もり。106では実施していない）",fontsize=13,
            transform=ax.transAxes,weight="bold")
    ax.text(0.06,0.185,"1回＝分オーダー × 60メンバー × 20サイクル = 数日",
            fontsize=13.5,transform=ax.transAxes)
    ax.text(0.06,0.125,f"→ $D$ を先に作っておくことで 約$10^7$〜$10^8$ 倍の差",
            fontsize=13,transform=ax.transAxes,color=RED,weight="bold")
    ax.text(0.02,0.015,"※ これが可能なのは熱弾性が線形だから。"
                       "塑性や接触があると $W$ が状態に依存し、この手は使えない",
            fontsize=11.5,transform=ax.transAxes,color="#556")

    fig.suptitle("$D$ 行列とは何か ― 「同化ループ中にFEMを一度も呼ばない」の中身",
                 fontsize=17,weight="bold",y=0.985)
    fig.tight_layout(rect=[0,0,1,0.945])
    o=os.path.join(IMG,"D_matrix_explained.png"); fig.savefig(o,dpi=125); plt.close(fig)
    print(f"行列積の実測 {us:.1f} µs")
    print("wrote",o)


if __name__=="__main__":
    main()
