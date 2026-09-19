"""計算力学用：推定対象 × 配置方法 の最適センサ配置表（ROM＋拡大状態EnKF）.

主張:「推定したい物理量によって最適なセンサ配置は変わる」。
- 推定対象(行): 全温度場 / 変位差Uz(A)-Uz(O) / 熱量Q / 熱伝達率h
- 配置方法(列): ランダム / 等間隔 / POD-Q-DEIM / 感度最大(対象別)
- 各セル: その配置(温度2点)でEnKFを回し、その対象の推定誤差(5seed平均)
- 行ごとに最良の配置法をハイライト → 対象で最良法が変わることを示す

温度センサは全場20,696セルから自由に選べる(gappy-POD)。感度最大は対象別:
  全温度場→分散(モード)最大, Q→|∂T/∂Q|最大, h→|∂T/∂h|最大, 変位→温度→変位感度最大。

出力: docs/img/optimal_placement_table.png, results/optimal_placement.npz
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/optimal_placement_table.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from sensitivity_analysis import integrate_sensitivity
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]


def pick2_far(cells_rank, Cc):
    """感度順に並んだ候補から、互いに離れた2点を選ぶ（近接重複を避ける）。"""
    c1=cells_rank[0]
    for c in cells_rank[1:]:
        if np.linalg.norm(Cc[c]-Cc[c1])>0.02: return [int(c1),int(c)]
    return [int(cells_rank[0]),int(cells_rank[1])]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])
    do=np.load(os.path.join(RES,"disp_operator.npz")); uzAO=do["uz_mean"]; Dao=do["Dmode"]
    Q_true=rg.HEATER_RATED_W; N=U.shape[0]
    Gfull=UPp.T@(U.T@U)@UPp/N                     # 全場MSE = dT5^T Gfull dT5

    # 感度（全場ピーク）
    ts,T,SQ,Sh=integrate_sensitivity(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,Q_true,heat,0,T_END,DT)
    SQf_pk=np.abs(SQ@UPp.T@U.T).max(0); Shf_pk=np.abs(Sh@UPp.T@U.T).max(0)
    varf=np.abs(U[:,0])*1.0                       # モード1振幅=場の分散(温度場再構成の感度)

    # 温度→変位差(A-O) の感度: d(qoi)/d(cell温度) の大きさ
    w_ao=Dao[0]-Dao[1]                            # (5,) モード係数→A-O
    disp_sens=np.abs(U@ (UPp.T@w_ao))             # (N,) 各セル温度がA-O差にどれだけ効くか

    # --- 配置方法 → 温度2セル ---
    rng0=np.random.default_rng(2026)
    place={}
    place["等間隔"]=[int(np.argmax(Cc[:,0])),int(np.argmin(Cc[:,0]))]     # +X端と-X端
    place["POD/Q-DEIM"]=[int(pod[0]),int(pod[2])]                          # Q-DEIM点2つ
    sens_by_target={
        "全温度場": pick2_far(np.argsort(-varf),Cc),
        "変位差":   pick2_far(np.argsort(-disp_sens),Cc),
        "熱量Q":    pick2_far(np.argsort(-SQf_pk),Cc),
        "熱伝達率h":pick2_far(np.argsort(-Shf_pk),Cc),
    }
    RAND_SETS=[list(rng0.choice(N,2,replace=False)) for _ in range(5)]     # ランダムは5組平均

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]; ht=(tg>0)&(tg<=300)
    Tt=np.full(NPT,rg.T_AIR_K); Ttr=[Tt.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(Tt,C,Kmat,h_true,1.0,heat,a,b,DT); Tt=tr[-1]; Ttr.append(Tt.copy())
    Ttr=np.array(Ttr)
    def qoiAO(T5): a=(T5-mean[pod])@UPp.T; u=uzAO+a@Dao.T; return u[...,0]-u[...,1]
    qoi_true=qoiAO(Ttr)

    def run_enkf(cells, seed):
        cells=list(cells)
        Hrows=np.array([U[c,:]@UPp for c in cells])           # (nc,5)
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.eye(len(cells))*SIG_T**2
        recT=[Z[:,:NPT].mean(0).copy()]; Qh=[(Z[:,IQ].mean()*rg.HEATER_RATED_W,Z[:,IH].mean())]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            base=np.array([mean[c] for c in cells])
            Yf=base+ (Z[:,:NPT]-mean[pod])@Hrows.T
            yv=base+ (Ttr[ci]-mean[pod])@Hrows.T + rng_o.normal(0,SIG_T,len(cells))
            Z=enkf_update(Z,yv,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            recT.append(Z[:,:NPT].mean(0).copy()); Qh.append((Z[:,IQ].mean()*rg.HEATER_RATED_W,Z[:,IH].mean()))
        recT=np.array(recT); Qh=np.array(Qh)
        dT5=recT-Ttr
        eT=np.sqrt((np.einsum("ti,ij,tj->t",dT5,Gfull,dT5))[ht].mean())      # 全場RMSE
        eU=np.abs(qoiAO(recT)-qoi_true)[ht].mean()
        eQ=np.abs(Qh[-1,0]-Q_true); eH=np.abs(Qh[-1,1]-h_true)
        return eT,eU,eQ,eH

    def eval_place(cells):
        r=np.array([run_enkf(cells,s) for s in SEEDS]).mean(0)
        return dict(eT=r[0],eU=r[1],eQ=r[2],eH=r[3])

    # ランダム=5組×5seed平均
    rand=np.mean([[run_enkf(cs,s) for s in SEEDS] for cs in RAND_SETS],axis=(0,1))
    methods=["ランダム","等間隔","POD/Q-DEIM","感度最大(対象別)"]
    targets=["全温度場","変位差","熱量Q","熱伝達率h"]; units=["K","µm","W","W/K"]
    # 対象非依存な3法は1回のEnKFで4指標
    res_uniform=eval_place(place["等間隔"]); res_qdeim=eval_place(place["POD/Q-DEIM"])
    res_rand=dict(eT=rand[0],eU=rand[1],eQ=rand[2],eH=rand[3])
    # 感度最大は対象別
    res_sens={t:eval_place(sens_by_target[t]) for t in targets}

    keymap={"全温度場":"eT","変位差":"eU","熱量Q":"eQ","熱伝達率h":"eH"}
    M=np.zeros((4,4))
    for ti,t in enumerate(targets):
        k=keymap[t]
        M[ti,0]=res_rand[k]; M[ti,1]=res_uniform[k]; M[ti,2]=res_qdeim[k]; M[ti,3]=res_sens[t][k]
    for ti,t in enumerate(targets):
        print(f"[tbl] {t:8s}: "+"  ".join(f"{m}={M[ti,j]:.3f}" for j,m in enumerate(methods)))

    # --- 表(ヒートマップ)図: 行ごとに正規化して最良を強調 ---
    fig,ax=plt.subplots(figsize=(11,5.6))
    Mn=M/M.max(1,keepdims=True)                    # 行ごと正規化(1=最悪)
    im=ax.imshow(Mn,cmap="RdYlGn_r",vmin=0,vmax=1,aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(methods,fontsize=12)
    ax.set_yticks(range(4)); ax.set_yticklabels([f"{t}\n[{u}]" for t,u in zip(targets,units)],fontsize=12)
    for ti in range(4):
        best=int(np.argmin(M[ti]))
        for j in range(4):
            txt=f"{M[ti,j]:.3f}"
            ax.text(j,ti,txt+("\n★最良" if j==best else ""),ha="center",va="center",
                    fontsize=11,weight=("bold" if j==best else "normal"),
                    color=("black"))
            if j==best: ax.add_patch(plt.Rectangle((j-.5,ti-.5),1,1,fill=False,edgecolor="blue",lw=3))
    ax.set_title("推定対象 × 配置方法：最適センサ配置は『推定したい量』で変わる（ROM＋EnKF, 5seed平均）\n"
                 "各行＝推定誤差（小さいほど良い）。青枠＝その対象での最良配置。色は行内の相対（緑=良/赤=悪）",
                 fontsize=12,weight="bold")
    fig.tight_layout()
    out=os.path.join(IMG,"optimal_placement_table.png"); fig.savefig(out,dpi=140); plt.close(fig)
    print("[tbl] wrote",out)
    np.savez(os.path.join(RES,"optimal_placement.npz"),M=M,methods=methods,targets=targets,
             sens_cells={t:sens_by_target[t] for t in targets})


if __name__=="__main__": main()
