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
    # 変位QoI（上面A/O差 Uz(A)-Uz(O)）を recT から出すための演算子
    do=np.load(os.path.join(RES,"disp_operator.npz")); uzAO=do["uz_mean"]; Dao=do["Dmode"]
    def qoiAO(T5):  a=(T5-mean[pod])@UPp.T; u=uzAO+a@Dao.T; return (u[...,0]-u[...,1])  # µm(Dmodeは既にµm)

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

    def run(mode, seed, da=True):
        """da=False なら観測補正を一切せず自由予測（＝データ同化なし）。"""
        H,uz0=obsop(mode); nobs=H.shape[0]
        R=np.diag([SIG_T**2]+([SIG_U**2]*2 if mode in ("hi","lo") else []))
        B=fixedB(seed)
        K=B@H.T@np.linalg.inv(H@B@H.T+R)              # 固定ゲイン(OI)
        rng=np.random.default_rng(seed+1); rng_o=np.random.default_rng(seed+7)
        z=np.zeros(NAUG); z[:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,NPT)
        z[IQ]=rng.uniform(0.3,1.8); z[IH]=np.clip(rng.normal(0.02,0.01),1e-3,0.1)
        rmse=[np.sqrt(((z[:NPT]-Ttr[0])**2).mean())]; recT=[z[:NPT].copy()]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            _,tr=rg.integrate_single(z[:NPT],C,Kmat,z[IH],z[IQ],heat,tp,tb,DT); z=z.copy(); z[:NPT]=tr[-1]; tp=tb
            if da:
                yv=y_of(Ttr[ci],mode,H,uz0)+rng_o.normal(0,np.sqrt(np.diag(R)))
                ypred=y_of(z[:NPT],mode,H,uz0)
                z=z+K@(yv-ypred)
                z[IQ]=np.clip(z[IQ],0,3); z[IH]=np.clip(z[IH],1e-4,0.2)
            rmse.append(np.sqrt(((z[:NPT]-Ttr[ci])**2).mean())); recT.append(z[:NPT].copy())
        return np.array(rmse), np.array(recT)

    out_rmse={}; out_recT={}
    # free=データ同化なし（自由予測）, none=温度のみ, lo/hi=温度+変位
    for m in ["free","none","lo","hi"]:
        obs_mode="none" if m=="free" else m
        rs=[run(obs_mode,s,da=(m!="free")) for s in SEEDS]
        out_rmse[m]=np.mean([r[0] for r in rs],axis=0)
        out_recT[m]=np.mean([r[1] for r in rs],axis=0)     # 5seed平均の解析温度
    ht=(tg>0)&(tg<=300)
    e={m:out_rmse[m][ht].mean() for m in out_rmse}
    # 変位QoI（真値・各構成）
    qoi_true=qoiAO(Ttr); qoi={m:qoiAO(out_recT[m]) for m in out_recT}
    print("[oi-disp] 加熱期平均RMSE[K]:", {k:round(v,3) for k,v in e.items()})

    fig,(ax0,ax1)=plt.subplots(1,2,figsize=(14.5,5.6))
    sty={"free":("-","tab:red","データ同化なし（自由予測）","x"),"none":(":","tab:gray","温度のみ（変位追加なし）","o"),"lo":("--","tab:orange","低W変位2点","s"),"hi":("--","tab:blue","高W変位2点","o")}
    # 左：温度（全5点RMSE）
    ax0.axvspan(0,300,color="orange",alpha=.06)
    for m,(ls,c,lab,mk) in sty.items():
        ax0.plot(tg,out_rmse[m],ls,color=c,lw=2.3,marker=mk,ms=4,label=f"{lab}  誤差{e[m]:.2f}K")
    ax0.set_xlabel("time [s]"); ax0.set_ylabel("全5点温度RMSE [K]"); ax0.grid(alpha=.3); ax0.legend()
    ax0.set_title("温度：推定の誤差（真値からのズレ）",fontsize=12)
    # 右：変位QoI Uz(A)-Uz(O)
    ax1.axvspan(0,300,color="orange",alpha=.06)
    ax1.plot(tg,qoi_true,"-",color="k",lw=3.6,alpha=.4,label="真値")
    for m,(ls,c,lab,mk) in sty.items():
        ax1.plot(tg,qoi[m],ls,color=c,lw=2.3,marker=mk,ms=4,label=lab)
    ax1.set_xlabel("time [s]"); ax1.set_ylabel("変位差 Uz(A)−Uz(O) [µm]"); ax1.grid(alpha=.3); ax1.legend()
    ax1.set_title("変位：同化後の温度から復元したA/O変位差",fontsize=12)
    fig.suptitle("OI（固定B・ROM）：観測構成で温度・変位の推定がどう変わるか（5seed平均）\n"
                 "赤＝データ同化なし（補正しないので真値から外れたまま）。温度P2で同化すると灰へ改善、さらに高W変位2点を足すと青が真値へ最速で追従",
                 fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.93]); out=os.path.join(IMG,"blog_oi_disp_selection.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[oi-disp] wrote",out)


if __name__=="__main__": main()
