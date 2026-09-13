"""run1(実ソルバOI)の温度・変位時刻歴を「観測点＋未観測点」を場所ごとに分けて描く.

温度は run_fem_oi/{t}/solid/T（=OI同化後の全場）と truth/{t}/solid/T（真値）から任意セルを取り出す。
→ 観測点(hot/cold)だけでなく未観測点も、DAで真値へ寄るかを見られる。
free-run(同化なし run_fem_oi_free)が全時刻そろっていれば no-DA も重畳する。

出力: docs/img/oi_fullsolver_temp.png（106へもコピー）
再現: OMP_NUM_THREADS=4 python3 run/make_oi_run1_fig.py
"""
from __future__ import annotations
import os, sys, shutil
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from daof import of_case
from dacore import plots as _p
import matplotlib.pyplot as plt
K=273.15
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
ANA=os.path.join(ROOT,"openfoam","run_fem_oi")           # run1(σ_Q=6W) の同化後全場
ANA2=os.path.join(ROOT,"openfoam","run_fem_oi_b")        # run2(σ_Q=2W) の同化後全場
TRUTH=os.path.join(ROOT,"openfoam","run_fem_enkf","truth")
FREE=os.path.join(ROOT,"openfoam","run_fem_oi_free")     # 同化なし
# 温度プローブ [m]: 観測2点 + 未観測2点
PROBES=[("hot（観測・ヒータ側+X中央）",(0.032,0,0.05025),True),
        ("cold（観測・反対側−X中央）",(-0.032,0,0.05025),True),
        ("未観測：上部ヒータ側(+X,z=90)",(0.030,0,0.090),False),
        ("未観測：底部反対側(−X,z=15)",(-0.030,0,0.015),False)]
T_AIR=293.15; T0_OFF=5.0


def read_at(case, t, cells):
    try: return of_case.read_solid_T(case, t)[cells]
    except Exception: return None


def main():
    h=dict(np.load(os.path.join(RES,"oi_fullsolver_history.npz"), allow_pickle=True))
    times=h["times"]; cyc=[x for x in times if x>0]
    centres=of_case.solid_cell_centres()
    cells=of_case.nearest_cells([p[1] for p in PROBES], centres)
    # 温度: 各時刻で 同化後(ANA)・真値(TRUTH)・同化なし(FREE) を読む
    def series(case, t0val):
        arr=[np.full(len(cells), t0val)]
        for t in cyc:
            v=read_at(case, t, cells); arr.append(v if v is not None else np.full(len(cells),np.nan))
        return np.array(arr)-K
    Tana=series(ANA, T_AIR+T0_OFF); Ttru=series(TRUTH, T_AIR)
    Tana2=series(ANA2, T_AIR+T0_OFF)                      # run2(σ_Q=2W)
    Tfree=series(FREE, T_AIR+T0_OFF)
    free_valid=np.isfinite(Tfree).all(axis=1)
    free_end=float(times[free_valid][-1])

    fig,axes=plt.subplots(2,3,figsize=(16.5,9))
    # 温度4点
    for k,(title,_xyz,obs) in enumerate(PROBES):
        ax=axes[k//3][k%3] if k<3 else axes[1][k-3]
        ax=axes.ravel()[k]; ax.axvspan(0,300,color="orange",alpha=.06)
        ax.plot(times,Ttru[:,k],"-",color="k",lw=3.4,alpha=.45,label="正解(真値)")
        ax.plot(times,Tfree[:,k],":",color="tab:orange",lw=2.3,label=f"同化なし（{free_end:g}sまで）")
        ax.plot(times,Tana[:,k],color="tab:blue",lw=2.0,ls="--",dashes=(5,2.5),marker="o",ms=4,label="OI σ_Q=6W")
        ax.plot(times,Tana2[:,k],color="tab:green",lw=2.0,ls="--",dashes=(2,2),marker="s",ms=4,label="OI σ_Q=2W")
        c="tab:red" if obs else "dimgray"
        ax.set_title(title,fontsize=11,color=("black" if obs else "dimgray"))
        ax.set_ylabel("温度 [degC]"); ax.set_xlabel("time [s]"); ax.grid(alpha=.3); ax.legend(fontsize=8,handlelength=3)
    # 変位2点(A/O, 観測)。free-run(同化なし)の変位は free の history npz から（完了後に入る）
    def loadnpz(name):
        p=os.path.join(RES,name); return dict(np.load(p, allow_pickle=True)) if os.path.exists(p) else None
    fh=loadnpz("oi_fullsolver_free_history.npz")   # 同化なし
    bh=loadnpz("oi_fullsolver_b_history.npz")       # run2(σ_Q=2W)
    for k,(lab,j) in enumerate([("変位：A＝ヒータ側 上面Uz(+X, 観測)",0),("変位：O＝反対側 上面Uz(−X, 観測)",1)]):
        ax=axes[1][1+k]; ax.axvspan(0,300,color="orange",alpha=.06)
        ax.plot(times,h["u_truth"][:,j]*1000,"-",color="k",lw=3.4,alpha=.45,label="正解(真値)")
        if fh is not None:
            ax.plot(fh["times"],fh["u_analysis"][:,j]*1000,":",color="tab:orange",lw=2.3,label="同化なし(free-run)")
        ax.plot(times,h["u_analysis"][:,j]*1000,color="tab:blue",lw=2.0,ls="--",dashes=(5,2.5),marker="o",ms=4,label="OI σ_Q=6W")
        if bh is not None:
            ax.plot(bh["times"],bh["u_analysis"][:,j]*1000,color="tab:green",lw=2.0,ls="--",dashes=(2,2),marker="s",ms=4,label="OI σ_Q=2W")
        ax.set_title(lab,fontsize=11); ax.set_ylabel("Uz [µm]"); ax.set_xlabel("time [s]"); ax.grid(alpha=.3); ax.legend(fontsize=8,handlelength=3)
    note=f"（同化なし：{free_end:g}sまで計算済み）"
    fig.suptitle("実ソルバOI：観測点＋未観測点の温度／変位 時刻歴 "+note+"\n"
                 "黒=真値、橙点線=同化なし、青破線=OI(σ_Q=6W)、緑破線=OI(σ_Q=2W)。背景1本・60秒ごと更新。",
                 fontsize=12.5,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.93])
    out=os.path.join(IMG,"oi_fullsolver_temp.png"); fig.savefig(out,dpi=135); plt.close(fig)
    try: shutil.copy(out, os.path.join(ROOT,"..","106_0_pod_selected_rom","docs","img","oi_fullsolver_temp.png")); print("[run1] copied to 106")
    except Exception as e: print("[run1] copy skip:",e)
    print(f"[run1] wrote oi_fullsolver_temp.png (未観測点あり, free-run重畳={Tfree is not None})")


if __name__=="__main__": main()
