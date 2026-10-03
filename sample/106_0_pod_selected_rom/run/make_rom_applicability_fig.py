"""docs/22_rom_applicability.md 用：ROM 単体（同化なし）の予測を、OpenFOAM の結果と3条件で比べる.

ROM は 15 W・0〜300 s の CHT 1本から作ったもの（results/rom_calibrated_pod.npz）。
各条件で、発熱量とヒータの ON/OFF を正しく与えて ROM を 0→600 s 前進させ、5点の温度を OpenFOAM と比べる。
間欠加熱では「0〜300 s に一定と仮定」した場合も重ねる。

出力: docs/img/rom_applicability.png, results/rom_applicability.json
再現: OMP_NUM_THREADS=4 python3 run/make_rom_applicability_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from improve_heating_schedule import integrate, SCHED
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
CASES=[("learned",1.0,"15 W・0〜300 s（ROM を作った条件）"),("q25",25/15,"25 W・0〜300 s（使っていない）"),("intermittent",1.0,"15 W 間欠加熱（使っていない）")]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz")); C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])
    pod=np.load(os.path.join(RES,"qdeim_points.npz"))["cell_idx"]
    fig,axs=plt.subplots(1,3,figsize=(16,4.8),sharey=True); out={}; cols=["#1f77b4","#ff7f0e","#C0392B","#2ca02c","#9467bd"]
    for ax,(case,q,lab) in zip(axs,CASES):
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz")); t=z["times"]; T5=z["Tfield"][:,pod]
        def run(on):
            T=np.full((1,5),rg.T_AIR_K); tr=[T[0].copy()]
            for k in range(1,len(t)): T=integrate(T,C,Km,np.array([h]),np.array([q]),heat,on,t[k-1],t[k]); tr.append(T[0].copy())
            return np.array(tr)
        R=run(SCHED[case]); e=np.abs(R-T5)
        out[case]=dict(max_err_K=float(e.max()),mean_err_K=float(e[1:].mean()),max_rise_K=float(T5.max()-rg.T_AIR_K))
        for i in range(5):
            ax.plot(t,T5[:,i]-273.15,"o",color=cols[i],ms=4)
            ax.plot(t,R[:,i]-273.15,"-",color=cols[i],lw=1.8,label=f"P{i}")
        txt=f"ROM の誤差：最大 {e.max():.3f} K"
        if case=="intermittent":
            W=run(lambda tt:(tt<300.0)); ew=np.abs(W-T5); out[case]["fixed_schedule_max_err_K"]=float(ew.max())
            ax.plot(t,W[:,2]-273.15,"--",color="k",lw=1.5,label="P2：0〜300 s 一定と仮定")
            txt+=f"\n（0〜300 s 一定と仮定すると最大 {ew.max():.2f} K）"
        ax.text(0.97,0.04,txt,transform=ax.transAxes,ha="right",fontsize=10.5,bbox=dict(fc="white",ec="#bbb"))
        ax.set_title(lab,fontsize=12); ax.set_xlabel("時刻 [s]"); ax.grid(alpha=.3)
    axs[0].set_ylabel("温度 [℃]"); axs[0].legend(fontsize=9,ncol=2,loc="upper left"); axs[2].legend(fontsize=9,loc="upper left")
    fig.suptitle("ROM 単体の予測（線）と OpenFOAM（点）：発熱量と ON/OFF を正しく与え、同化なしで 0→600 s",fontsize=12.5)
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(IMG,"rom_applicability.png"),dpi=150); plt.close(fig)
    json.dump(out,open(os.path.join(RES,"rom_applicability.json"),"w"),ensure_ascii=False,indent=1); print(out)


if __name__=="__main__": main()
