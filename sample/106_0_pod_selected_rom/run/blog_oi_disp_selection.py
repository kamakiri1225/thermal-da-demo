"""blog_002用：高W/低W変位観測点の比較を「OI(固定B)」で回す.

run_disp_selection.py はEnKFだが、blog_002はOIの記事なので、同じ比較をOIで実施する。
状態 z=[T1..T5, Q, h](7次元)。固定背景共分散Bを初期アンサンブルの標本から一度だけ作って凍結し、
固定ゲイン K=B H^T (H B H^T + R)^{-1} を毎サイクル使う（＝OI）。
観測: 温度P2の1点に、①なし ②高W変位2点 ③低W変位2点 を足した3構成を比較。

出力: docs/img/blog_oi_disp_selection.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/blog_oi_disp_selection.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0; SIG_T=0.30; SIG_U=0.3
SEEDS=[20260913,20260914,20260915,20260916,20260917]
N_B=200   # 固定B作成用の初期サンプル数


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    heat=int(d["heat_node"]); C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])                     # (r,5)
    hi=np.load(os.path.join(RES,"dispop_hiW.npz")); lo=np.load(os.path.join(RES,"dispop_loW.npz"))
    # 変位観測演算子（T5→変位2点）: u = uz0 + ((T5-mean_p)@UPp) @ D^T
    def dispH(D):  return D@UPp                       # (2,5)  u ≈ uz0 + Hd@(T5-mean_p)
    Hd_hi=dispH(hi["D"]); Hd_lo=dispH(lo["D"])
    tP=2                                              # 温度観測点 P2

    # 真値トラジェクトリ
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def obsop(mode):
        """観測演算子 H (nobs x 7) と 参照変位オフセット uz0 を返す。"""
        rows=[np.eye(NAUG)[tP]]                       # 温度P2
        uz0=[]
        if mode in ("hi","lo"):
            Hd=Hd_hi if mode=="hi" else Hd_lo
            uz=hi["uz_mean"] if mode=="hi" else lo["uz_mean"]
            for j in range(2):
                row=np.zeros(NAUG); row[:NPT]=Hd[j]; rows.append(row); uz0.append(uz[j])
        return np.array(rows), np.array(uz0)

    def y_of(Tv, mode, H, uz0):
        """真値状態Tv(5)から観測ベクトルを作る。"""
        z=np.r_[Tv,0,0]; y=H@z
        if mode in ("hi","lo"): y[1:]=uz0+ (H[1:,:NPT]@(Tv-mean[pod]))  # 変位は基準+線形
        return y

    def fixedB(seed):
        """初期アンサンブルの標本からBを一度だけ作る（凍結）。"""
        rng=np.random.default_rng(seed)
        Z=np.zeros((N_B,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_B,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_B); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_B),1e-3,0.1)
        dz=Z-Z.mean(0); return dz.T@dz/(N_B-1)

    def run(mode, seed):
        H,uz0=obsop(mode); nobs=H.shape[0]
        R=np.diag([SIG_T**2]+([SIG_U**2]*2 if mode in ("hi","lo") else []))
        B=fixedB(seed)
        K=B@H.T@np.linalg.inv(H@B@H.T+R)              # 固定ゲイン(OI)
        rng=np.random.default_rng(seed+1); rng_o=np.random.default_rng(seed+7)
        z=np.zeros(NAUG); z[:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,NPT)
        z[IQ]=rng.uniform(0.3,1.8); z[IH]=np.clip(rng.normal(0.02,0.01),1e-3,0.1)
        rmse=[np.sqrt(((z[:NPT]-Ttr[0])**2).mean())]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            _,tr=rg.integrate_single(z[:NPT],C,Kmat,z[IH],z[IQ],heat,tp,tb,DT); z=z.copy(); z[:NPT]=tr[-1]; tp=tb
            yv=y_of(Ttr[ci],mode,H,uz0)+rng_o.normal(0,np.sqrt(np.diag(R)))
            ypred=y_of(z[:NPT],mode,H,uz0)
            z=z+K@(yv-ypred)
            z[IQ]=np.clip(z[IQ],0,3); z[IH]=np.clip(z[IH],1e-4,0.2)
            rmse.append(np.sqrt(((z[:NPT]-Ttr[ci])**2).mean()))
        return np.array(rmse)

    res={m:np.mean([run(m,s) for s in SEEDS],axis=0) for m in ["none","lo","hi"]}
    ht=(tg>0)&(tg<=300)
    e={m:res[m][ht].mean() for m in res}
    print("[oi-disp] 加熱期平均RMSE[K]:", {k:round(v,3) for k,v in e.items()})
    fig,ax=plt.subplots(figsize=(9,5.4)); ax.axvspan(0,300,color="orange",alpha=.06)
    ax.plot(tg,res["none"],":",color="tab:gray",lw=2.4,label=f"温度1点のみ（変位なし）  誤差{e['none']:.2f}K")
    ax.plot(tg,res["lo"],"--",color="tab:orange",lw=2.2,marker="s",ms=4,label=f"＋低W変位2点  誤差{e['lo']:.2f}K")
    ax.plot(tg,res["hi"],"--",color="tab:blue",lw=2.2,marker="o",ms=4,label=f"＋高W変位2点  誤差{e['hi']:.2f}K")
    ax.set_xlabel("time [s]"); ax.set_ylabel("全5点温度RMSE [K]")
    ax.set_title("OI（固定B）：高W／低Wの変位観測点で温度推定がどう変わるか\n"
                 "温度P2に変位2点を追加。高Wは速く下がり、低Wは変位なしとほぼ同じ（5seed平均）",fontsize=12)
    ax.grid(alpha=.3); ax.legend()
    fig.tight_layout(); out=os.path.join(IMG,"blog_oi_disp_selection.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[oi-disp] wrote",out)


if __name__=="__main__": main()
