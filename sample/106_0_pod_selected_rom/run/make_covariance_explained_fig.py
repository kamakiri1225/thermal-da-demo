"""「散布図が共分散」とはどういうことかを、実データの手計算で示す.

共分散は Cov(T,Q) = Σ_m (T^(m)-T̄)(Q^(m)-Q̄) / (N-1)。
各メンバーについて「温度の平均からのズレ」×「Qの平均からのズレ」を掛け、平均する。
散布図が斜めに伸びていれば積が同符号に偏り、Cov が非ゼロになる。
つまり散布図の形はそのまま共分散の計算を絵にしたものである。

出力: docs/img/covariance_explained.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_covariance_explained_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KC=273.15; NPT=5; IQ=5; IH=6; NAUG=7
DT=2.0; OBS_DT=30.0; T_END=330.0; N=60; SIG=0.30; P=2


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; K=rg.tri_to_matrix(d["K_upper"],NPT); h_t=float(d["h"]); heat=int(d["heat_node"])
    rng=np.random.default_rng(20260913); rng_o=np.random.default_rng(7)
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); tt=[0.0]; tr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        ts,x=rg.integrate_single(T,C,K,h_t,1.0,heat,a,b,DT)
        tt+=list(ts[1:]); tr+=list(x[1:]); T=x[-1]
    tt=np.array(tt); truth=np.array(tr)[:,P]-KC

    Z=np.zeros((N,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N),1e-3,0.1)
    snap0=(Z[:,P].copy()-KC, Z[:,IQ].copy()*15)
    tp=0.0
    for ci,tb in enumerate(cyc,1):
        Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,K,Z[:,IH],Z[:,IQ],heat,tp,tb,DT)
        yv=truth[np.argmin(abs(tt-tb))]+KC+rng_o.normal(0,SIG)
        Z=enkf_update(Z,np.array([yv]),None,np.diag([SIG**2]),rng,inflation=1.02,Yf=Z[:,[P]])
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        tp=tb
    snapN=(Z[:,P].copy()-KC, Z[:,IQ].copy()*15)

    def stats(sn):
        x,y=sn; dx,dy=x-x.mean(),y-y.mean(); pr=dx*dy
        return x,y,dx,dy,pr,float(pr.sum()/(N-1)),float(np.corrcoef(x,y)[0,1])

    fig=plt.figure(figsize=(19.0,9.6))
    gs=fig.add_gridspec(2,3,height_ratios=[1.12,1],hspace=.42,wspace=.30)

    for col,(sn,ttl) in enumerate([(snap0,"サイクル0（前進前）"),(snapN,"サイクル11（11回前進後）")]):
        x,y,dx,dy,pr,cov,cor=stats(sn)
        ax=fig.add_subplot(gs[0,col])
        # 4象限を塗り分け：積の符号がそのまま見える
        xm,ym=x.mean(),y.mean()
        xr=[x.min()-0.05*np.ptp(x),x.max()+0.05*np.ptp(x)]
        yr=[y.min()-0.05*np.ptp(y),y.max()+0.05*np.ptp(y)]
        ax.add_patch(plt.Rectangle((xm,ym),xr[1]-xm,yr[1]-ym,fc="#ffe9e4",zorder=0))
        ax.add_patch(plt.Rectangle((xr[0],yr[0]),xm-xr[0],ym-yr[0],fc="#ffe9e4",zorder=0))
        ax.add_patch(plt.Rectangle((xm,yr[0]),xr[1]-xm,ym-yr[0],fc="#e4ecff",zorder=0))
        ax.add_patch(plt.Rectangle((xr[0],ym),xm-xr[0],yr[1]-ym,fc="#e4ecff",zorder=0))
        ax.axvline(xm,color="0.4",lw=1.6,ls="--"); ax.axhline(ym,color="0.4",lw=1.6,ls="--")
        ax.scatter(x,y,s=34,c=np.where(pr>0,"#C0392B","#1F4E9C"),zorder=4)
        ax.set_xlim(xr); ax.set_ylim(yr)
        ax.text(xr[1],yr[1],"積が＋",ha="right",va="top",fontsize=12,color="#C0392B",weight="bold")
        ax.text(xr[0],yr[0],"積が＋",ha="left",va="bottom",fontsize=12,color="#C0392B",weight="bold")
        ax.text(xr[1],yr[0],"積が−",ha="right",va="bottom",fontsize=12,color="#1F4E9C",weight="bold")
        ax.text(xr[0],yr[1],"積が−",ha="left",va="top",fontsize=12,color="#1F4E9C",weight="bold")
        ax.set_xlabel("メンバーの温度 $T^{(m)}$ [℃]",fontsize=12.5)
        ax.set_ylabel("メンバーの発熱量 $Q^{(m)}$ [W]",fontsize=12.5)
        npos=int((pr>0).sum())
        ax.set_title(f"{ttl}\n＋が{npos}個 / −が{N-npos}個 → "
                     f"$\\mathrm{{Cov}}={cov:+.3f}$，相関 {cor:+.3f}",
                     fontsize=13.5,weight="bold",
                     color=("#556" if abs(cor)<0.15 else "#C0392B"))
        ax.grid(alpha=.25)

    # 右上：定義
    ax=fig.add_subplot(gs[0,2]); ax.axis("off")
    ax.set_title("共分散の定義 ― 散布図の形そのもの",fontsize=14.5,weight="bold")
    ax.text(0.02,0.86,"各メンバーについて",fontsize=13,transform=ax.transAxes)
    ax.text(0.05,0.70,r"$(T^{(m)}-\bar T)\times(Q^{(m)}-\bar Q)$",fontsize=16,transform=ax.transAxes)
    ax.text(0.02,0.58,"を計算し、全部足して $N-1$ で割る：",fontsize=13,transform=ax.transAxes)
    ax.text(0.03,0.44,r"$\mathrm{Cov}(T,Q)=\dfrac{1}{N-1}\sum_{m=1}^{N}"
                      r"(T^{(m)}-\bar T)(Q^{(m)}-\bar Q)$",
            fontsize=15,transform=ax.transAxes,color="#C0392B")
    ax.text(0.02,0.24,"・赤の象限（両方とも平均より上／両方とも下）は積が＋\n"
                      "・青の象限（片方だけ上）は積が−",
            fontsize=12.5,transform=ax.transAxes)
    ax.text(0.02,0.08,"丸い雲 → ＋と−が同数で打ち消し合う → $\\mathrm{Cov}\\simeq0$\n"
                      "斜めに伸びる → 片方の符号に偏る → $\\mathrm{Cov}\\neq0$",
            fontsize=12.5,transform=ax.transAxes,weight="bold")
    ax.text(0.02,-0.08,"つまり「散布図が斜めかどうか」を数値にしたものが共分散",
            fontsize=12.5,transform=ax.transAxes,color="#C0392B",weight="bold")

    # 下段：実際のメンバーで手計算
    ax=fig.add_subplot(gs[1,:]); ax.axis("off")
    x,y,dx,dy,pr,cov,cor=stats(snapN)
    ax.set_title("実際に手で計算してみる（サイクル11。60メンバーから4つ抜粋）",
                 fontsize=14.5,weight="bold")
    idx=list(np.argsort(dx)[:2])+list(np.argsort(dx)[-2:])
    hdr=["メンバー","温度 $T^{(m)}$ [℃]","ズレ $T-\\bar T$","発熱量 $Q^{(m)}$ [W]",
         "ズレ $Q-\\bar Q$","積"]
    xs=[0.05,0.19,0.34,0.47,0.62,0.76]
    for xx,hh in zip(xs,hdr):
        ax.text(xx,0.80,hh,fontsize=12.5,weight="bold",transform=ax.transAxes)
    for i,m in enumerate(idx):
        yy=0.64-i*0.115
        col="#C0392B" if pr[m]>0 else "#1F4E9C"
        for xx,v in zip(xs,[f"#{m}",f"{x[m]:.3f}",f"{dx[m]:+.3f}",f"{y[m]:.2f}",
                            f"{dy[m]:+.2f}",f"{pr[m]:+.4f}"]):
            ax.text(xx,yy,v,fontsize=12.5,transform=ax.transAxes,
                    color=(col if xx==xs[-1] else "#1b2430"),
                    weight=("bold" if xx==xs[-1] else "normal"))
    ax.text(0.05,0.14,"…（残り56メンバーも同様）",fontsize=12,transform=ax.transAxes,color="#556")
    ax.text(0.05,0.02,
            f"全部足すと $\\sum=${pr.sum():+.3f}　→　"
            f"$\\mathrm{{Cov}}={pr.sum():+.3f}/{N-1}={cov:+.4f}$ [℃·W]　"
            f"→　相関 $={cov:+.4f}/({dx.std(ddof=1):.4f}\\times{dy.std(ddof=1):.3f})={cor:+.3f}$",
            fontsize=13.5,transform=ax.transAxes,color="#C0392B",weight="bold")
    ax.text(0.05,-0.10,"この $\\mathrm{Cov}(T,Q)$ が非ゼロだから、"
                       "温度の残差がゲイン $K_Q=\\mathrm{Cov}(T,Q)/(\\mathrm{Var}(T)+r)$ を通って $Q$ を動かす",
            fontsize=12.5,transform=ax.transAxes,color="#556")

    fig.suptitle("「散布図が共分散」とはどういうことか ― 定義に戻って手で計算する",
                 fontsize=17,weight="bold",y=0.985)
    fig.tight_layout(rect=[0,0.03,1,0.94])
    o=os.path.join(IMG,"covariance_explained.png"); fig.savefig(o,dpi=125); plt.close(fig)
    print(f"サイクル0 相関 {stats(snap0)[6]:+.3f} / サイクル11 相関 {cor:+.3f}, Cov={cov:+.4f}")
    print("wrote",o)


if __name__=="__main__":
    main()
