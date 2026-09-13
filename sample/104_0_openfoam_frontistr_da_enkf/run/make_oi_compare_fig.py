"""実ソルバOIを σ_Q を変えて2回回した結果を比較する図（＝OIは調整が必要で難しい）.

比較対象: σ_Q=0.5, 2, 6, 50 W
温度時刻歴（観測点）・変位時刻歴（上面2点）・発熱Qの推定 を並べ、
「決め打ちσ_Qを変えると結果（特にQの過補正・収束）が変わる」ことを示す。

出力: docs/img/oi_fullsolver_compare.png（106へもコピー）
再現: OMP_NUM_THREADS=4 python3 run/make_oi_compare_fig.py
"""
from __future__ import annotations
import os, sys, shutil
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
K=273.15
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
RUNS=[("sq0p5",0.5,"tab:purple"),("b",2.0,"tab:green"),("",6.0,"tab:blue"),("sq50",50.0,"tab:brown")]


def main():
    runs={}
    for suf, sig, color in RUNS:
        p=os.path.join(RES,f"oi_fullsolver{('_'+suf) if suf else ''}_history.npz")
        if not os.path.exists(p):
            print(f"[cmp] まだ無い: {p}（両方の run 完了後に実行）"); return
        runs[suf]=dict(np.load(p, allow_pickle=True))
    t=runs["" ]["times"]
    cols={suf:color for suf,sig,color in RUNS}
    lab={suf:f"OI σ_Q={sig:g}W" for suf,sig,color in RUNS}

    fig,axes=plt.subplots(1,3,figsize=(17,5.2))
    # 温度（hot観測点）
    ax=axes[0]; ax.axvspan(0,300,color="orange",alpha=.06)
    ax.plot(t, runs[""]["T_obs_truth"][:,0]-K,"-",color="k",lw=3.5,alpha=.4,label="真値(hot)")
    for suf, sig, color in RUNS:
        ax.plot(t, runs[suf]["T_obs_analysis"][:,0]-K,"--o",ms=3,color=cols[suf],label=lab[suf])
    ax.set_xlabel("time [s]"); ax.set_ylabel("hot観測点 温度 [degC]"); ax.grid(alpha=.3); ax.legend(fontsize=9)
    ax.set_title("温度時刻歴（観測点hot）")
    # 変位（ヒータ側上面）
    ax=axes[1]; ax.axvspan(0,300,color="orange",alpha=.06)
    ax.plot(t, runs[""]["u_truth"][:,0]*1000,"-",color="k",lw=3.5,alpha=.4,label="真値(A:ヒータ側Uz)")
    for suf, sig, color in RUNS:
        ax.plot(t, runs[suf]["u_analysis"][:,0]*1000,"--o",ms=3,color=cols[suf],label=lab[suf])
    ax.set_xlabel("time [s]"); ax.set_ylabel("A：ヒータ側上面 Uz [µm]"); ax.grid(alpha=.3); ax.legend(fontsize=9)
    ax.set_title("変位時刻歴（上面 A=ヒータ側+X）")
    # 発熱Q
    ax=axes[2]; ax.axvspan(0,300,color="orange",alpha=.06)
    ax.axhline(15.0,color="k",lw=1.5,label="真値 Q=15W")
    for suf, sig, color in RUNS:
        ax.plot(t, runs[suf]["q"],"--o",ms=3,color=cols[suf],label=lab[suf]+f"  最終{runs[suf]['q'][-1]:.1f}W")
    ax.set_xlabel("time [s]"); ax.set_ylabel("発熱量 Q の推定 [W]"); ax.grid(alpha=.3); ax.legend(fontsize=9)
    ax.set_title("発熱Qの推定：σ_Q次第で過補正・収束が変わる")
    fig.suptitle("実ソルバOI：決め打ち σ_Q を変えた4ケースの比較（真値 Q=15 W）",
                 fontsize=14,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.94])
    out='/mnt/d/work/002_cae/openfoam/20260505_datadoka/sample/104_0_openfoam_frontistr_da_enkf/docs/img/oi_fullsolver_compare.png'; fig.savefig(out,dpi=140); plt.close(fig)
    # 106のdocs/imgへもコピー（00から参照するため）
    dst='/mnt/d/work/002_cae/openfoam/20260505_datadoka/sample/106_0_pod_selected_rom/docs/img/oi_fullsolver_compare.png'
    try: shutil.copy(out, dst); print("[cmp] copied to 106/docs/img")
    except Exception as e: print("[cmp] copy skip:",e)
    for suf, sig, color in RUNS:
        q=runs[suf]["q"]; r=runs[suf]["rmse_field"]
        print(f"[cmp] σ_Q={sig:g}W: 場RMSE {r[0]:.2f}->{r[-1]:.3f}K, Q {q[0]:.0f}->{q[-1]:.2f}W")
    print("[cmp] wrote docs/img/oi_fullsolver_compare.png")


if __name__=="__main__": main()
