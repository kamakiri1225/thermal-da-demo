"""従来の回帰補正とデータ同化を、新しく計算した条件運転条件で比べる.

従来法（工作機械の熱変位補正で一般的な重回帰）:
  反り U_z(A)−U_z(O) ＝ c0 + c1 ΔT_s1 + c2 ΔT_s2      （ΔT＝センサ温度−20 ℃）
  回帰式の係数を決めるデータ：15 W・0〜300 s の CHT（learned）の 0〜600 s、30 s ごと21時刻。
  入力は同化と同じ温度センサ2点（P2+P0 または P2+P4）。テストでは同化と同じく温度に 0.3 K のノイズ（同じ5 seed）。
データ同化：run/limit_realsolver_truth.py の結果（同じ真値・同じセンサ）。

真値はどれも OpenFOAM の温度場＋それを入れた FrontISTR の変位（results/limit_truth_<case>.npz）。

出力: results/baseline_regression_vs_da.json, docs/img/baseline_regression_vs_da.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/baseline_regression_vs_da.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
SEEDS=[20260913,20260914,20260915,20260916,20260917]; SIG_T=0.30; TAIR=293.15
CASES=[("learned","15 W・0〜300 s\n（ROM の係数を決めた計算）"),("q25","25 W・0〜300 s\n（新しく計算した条件）"),("intermittent","15 W 間欠加熱\n（新しく計算した条件）")]
SETS=[("P2+P0",[2,0],"温度2点 P2+P0（A側だけ）","温度 P2+P0＋変位 A/O"),("P2+P4",[2,4],"温度2点 P2+P4（両側）","温度 P2+P4＋変位 A/O")]


def load(case):
    f=os.path.join(RES,f"limit_truth_{case}.npz")
    if not os.path.exists(f): return None
    d=np.load(f); pod=np.load(os.path.join(RES,"qdeim_points.npz"))["cell_idx"]
    return d["times"],d["Tfield"][:,pod],d["uz"][:,0]-d["uz"][:,1]


def main():
    tr=load("learned"); t0,T0,ao0=tr
    out={"note":"反りA−Oの平均絶対誤差 [µm]。heat=30〜300 s、all=30〜600 s。5 seed 平均","cases":{}}
    for case,lab in CASES:
        L=load(case); J=os.path.join(RES,f"limit_realsolver_{case}.json")
        if L is None or not os.path.exists(J): print("skip",case); continue
        t,T,ao=L; da=json.load(open(J))["configs"]; heat=(t>0)&(t<=300); allw=t>0
        res={}
        for key,nodes,da_t,da_td in SETS:
            X0=np.column_stack([np.ones(len(t0)),T0[:,nodes]-TAIR]); c,*_=np.linalg.lstsq(X0,ao0,rcond=None)
            errs=[]
            for s in SEEDS:
                Tn=T[:,nodes]+np.random.default_rng(s+7).normal(0,SIG_T,(len(t),len(nodes)))
                pred=np.column_stack([np.ones(len(t)),Tn-TAIR])@c; e=np.abs(pred-ao)
                errs.append((e[heat].mean(),e[allw].mean()))
            errs=np.array(errs)
            res[key]=dict(coef=c.tolist(),regression_heat=float(errs[:,0].mean()),regression_all=float(errs[:,1].mean()),
                          da_temp_heat=da[da_t]["heat"]["AO_um"],da_temp_all=da[da_t]["all"]["AO_um"],
                          da_tempdisp_heat=da[da_td]["heat"]["AO_um"],da_tempdisp_all=da[da_td]["all"]["AO_um"])
        res["truth_AO_max_um"]=float(np.abs(ao).max()); res["no_DA_heat"]=da["同化なし"]["heat"]["AO_um"]
        out["cases"][case]=res
        print(case,json.dumps(res,ensure_ascii=False))
    json.dump(out,open(os.path.join(RES,"baseline_regression_vs_da.json"),"w"),ensure_ascii=False,indent=1)
    cs=[c for c,_ in CASES if c in out["cases"]]
    fig,axs=plt.subplots(1,len(cs),figsize=(5.2*len(cs),4.8),squeeze=False)
    for ax,case in zip(axs[0],cs):
        r=out["cases"][case]; k="P2+P4"
        names=["回帰式\n（温度2点）","データ同化\n（温度2点）","データ同化\n（温度2点＋変位2点）"]
        vals=[r[k]["regression_all"],r[k]["da_temp_all"],r[k]["da_tempdisp_all"]]
        ax.bar(range(3),vals,color=["#9AA5B1","#2E8B57","#C0392B"])
        for i,v in enumerate(vals): ax.text(i,v,f"{v:.2f}",ha="center",va="bottom",fontsize=11)
        ax.set_xticks(range(3)); ax.set_xticklabels(names,fontsize=10)
        ax.set_title(dict(CASES)[case],fontsize=12); ax.set_ylabel("反り A−O の誤差 [µm]（30〜600 s 平均）"); ax.grid(axis="y",alpha=.3)
    fig.suptitle("従来の回帰補正とデータ同化の比較（温度センサ P2+P4、真値＝OpenFOAM＋FrontISTR、5 seed）",fontsize=12.5)
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(IMG,"baseline_regression_vs_da.png"),dpi=150); plt.close(fig)


if __name__=="__main__": main()
