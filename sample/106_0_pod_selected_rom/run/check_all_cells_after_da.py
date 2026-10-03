"""同化で復元した温度場が、代表点以外のセルでも合っているかを全20,696セルで確かめる.

真値は OpenFOAM の温度場そのもの（results/limit_truth_<case>.npz）。
観測は温度2点(P2+P4)＋変位2点(A/O)。EnKF・60メンバー・30秒ごと・5 seed。

調べること
  ① 全セルの誤差の分布（ヒストグラムと、半径・高さ別の分布図）
  ② 代表5点 と それ以外20,691セル で誤差に差があるか
  ③ 誤差が最大のセルはどこか
  ④ 誤差の時間変化（最大・95%・平均）

出力: results/check_all_cells_after_da.json, docs/img/all_cells_after_da.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/check_all_cells_after_da.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KC0=273.15
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TN=[2,4]
CASES=[("learned","15 W（ROM の係数を決めた計算）"),("q25","25 W（新しく計算した条件）"),
       ("intermittent","15 W 間欠加熱（新しく計算した条件）")]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:]); A=U@UPp
    op=np.load(os.path.join(RES,"disp_operator.npz")); ub=op["uz_mean"]; Dm=op["Dmode"]
    disp=lambda T5: ub+((T5-mean[pod])@UPp.T)@Dm.T
    field=lambda T5: mean+A@(T5-mean[pod])
    ncell=U.shape[0]; other=np.ones(ncell,bool); other[pod]=False
    out={"note":"真値＝OpenFOAMの温度場。観測＝温度2点(P2+P4)＋変位2点(A/O)。誤差は|推定−真値|、加熱期30〜300秒、5 seed平均","cases":{}}
    fig,axs=plt.subplots(4,len(CASES),figsize=(16.5,15.5),
                         gridspec_kw=dict(hspace=0.38,wspace=0.24,height_ratios=[1,1,1,0.85]))
    for j,(case,lab) in enumerate(CASES):
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,pod]; uz=z["uz"]
        cyc=t[1:]; heatw=(t[1:]<=300)
        errs=np.zeros((len(SEEDS),len(cyc),ncell)); ests=np.zeros_like(errs)
        for si,seed in enumerate(SEEDS):
            rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
            Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
            Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
            Rd=np.diag([SIG_T**2]*len(TN)+[SIG_U**2]*2); tp=0.0
            for ci,tb in enumerate(cyc,1):
                Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
                y=np.r_[T5[ci][TN],uz[ci]]+ro.normal(0,np.sqrt(np.diag(Rd)))
                Yf=np.column_stack([Z[:,TN],disp(Z[:,:NPT])])
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
                fe=field(Z[:,:NPT].mean(0))
                errs[si,ci-1]=np.abs(fe-Tf[ci]); ests[si,ci-1]=fe
        e=errs.mean(0)                       # (time, cell) 5 seed平均
        est=ests.mean(0)                     # 同化後の温度場（5 seed平均）
        eh=e[heatw].mean(0)                  # 加熱期平均（セルごと）
        imax=int(eh.argmax())
        rise=float(Tf.max()-rg.T_AIR_K)
        out["cases"][case]=dict(
            all_mean_K=float(eh.mean()),all_p95_K=float(np.percentile(eh,95)),all_max_K=float(eh.max()),
            pod5_mean_K=float(eh[pod].mean()),other_mean_K=float(eh[other].mean()),
            worst_cell=imax,worst_xyz_mm=np.round(Cc[imax]*1000,1).tolist(),
            max_rise_K=rise,max_err_pct_of_rise=float(100*eh.max()/rise),
            frac_below_noise=float((eh<SIG_T).mean()))
        print(f"{case}: 全セル平均 {eh.mean():.3f} K / 95% {np.percentile(eh,95):.3f} / 最大 {eh.max():.3f} K "
              f"(代表5点 {eh[pod].mean():.3f}, それ以外 {eh[other].mean():.3f})  ノイズ0.3K未満のセル {100*(eh<SIG_T).mean():.1f}%",flush=True)
        # ── 温度そのもの（主役）：誤差が最小・中央・最大のセル ──
        order=np.argsort(eh)
        picks=[(int(order[0]),"誤差が最小","#2E6FD8"),
               (int(order[len(order)//2]),"中央（50%点）","#E67E22"),
               (int(order[-1]),"誤差が最大","#C0392B")]
        for row,(c_,nm_,col_) in enumerate(picks):
            ax=axs[row,j]; ax.axvspan(0,300,color="#FDEBD0",alpha=.40)
            ax.plot(t,Tf[:,c_]-KC0,"o",color="k",ms=5,label="OpenFOAM（真値）")
            ax.plot(t[1:],est[:,c_]-KC0,"-",color=col_,lw=2.3,label="同化後の推定")
            ax.fill_between(t[1:],est[:,c_]-KC0-SIG_T,est[:,c_]-KC0+SIG_T,color=col_,alpha=.20,
                            label="±0.3 K（温度計のノイズ）")
            ax.grid(alpha=.3)
            xyz=np.round(Cc[c_]*1000,0).astype(int).tolist()
            ax.set_title(f"{nm_}のセル {xyz} mm　平均誤差 {eh[c_]:.3f} K",fontsize=11)
            if j==0: ax.set_ylabel("温度 [℃]")
            if row==0 and j==0: ax.legend(fontsize=9,loc="lower right")
            ax.tick_params(labelbottom=False)
            if row==0:
                ax.text(0.03,0.95,lab,transform=ax.transAxes,va="top",fontsize=11,weight="bold",
                        bbox=dict(fc="white",ec="#bbb"))
        # ── 補助：全セルの誤差の分布 ──
        ax=axs[3,j]
        ax.hist(eh[other],bins=60,color="#8FA8C8",alpha=.9,label=f"代表点以外 {other.sum():,} セル")
        for pi in pod: ax.axvline(eh[pi],color="#C0392B",lw=1.4)
        ax.axvline(SIG_T,color="k",ls="--",lw=2,label="温度計のノイズ 0.3 K")
        for c_,nm_,col_ in picks: ax.axvline(eh[c_],color=col_,lw=2.4,ls=":")
        ax.set_xlabel("加熱期の平均誤差 [K]"); ax.grid(alpha=.3)
        if j==0: ax.set_ylabel("セル数"); ax.legend(fontsize=9)
        ax.set_title(f"全セルの誤差　平均 {eh.mean():.3f} / 最大 {eh.max():.3f} K　0.3 K未満 {100*(eh<SIG_T).mean():.0f}%",
                     fontsize=10.5)
    json.dump(out,open(os.path.join(RES,"check_all_cells_after_da.json"),"w"),ensure_ascii=False,indent=1)
    fig.suptitle("同化で復元した温度は、代表5点以外のセルでも合っているか\n黒点＝OpenFOAM（真値）、色線＝同化後の推定、帯＝温度計のノイズ ±0.3 K。観測は温度2点＋変位2点、5 seed 平均",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.945)); fig.savefig(os.path.join(IMG,"all_cells_after_da.png"),dpi=150); plt.close(fig)
    print("wrote all_cells_after_da.png")


if __name__=="__main__": main()
