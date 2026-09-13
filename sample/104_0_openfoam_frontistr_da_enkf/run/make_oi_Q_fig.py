"""実ソルバOIの発熱量 Q の推定値 Q(t) を、σ_Q 違い＋同化なしで比較する図.

各runのログ（openfoam/run_fem_oi*.log）の "analysis t=.. Q=.." をパースして Q(t) を描く。
真値 Q=15W、でたらめ初期 Q0=22W を基準線に。

出力: docs/img/oi_Q_estimate.png（106へもコピー）
再現: OMP_NUM_THREADS=4 python3 run/make_oi_Q_fig.py
"""
from __future__ import annotations
import os, re, sys, shutil
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img"); OFDIR=os.path.join(ROOT,"openfoam")
Q0=22.0; QTRUE=15.0
# (ログ名, ラベル, 色, 線種)
RUNS=[("run_fem_oi_sq0p5.log","OI σ_Q=0.5W（小）","tab:purple","--"),
      ("run_fem_oi_b.log",    "OI σ_Q=2W",        "tab:green","--"),
      ("run_fem_oi.log",      "OI σ_Q=6W",        "tab:blue","--"),
      ("run_fem_oi_sq50.log", "OI σ_Q=50W（大）","tab:brown","--"),
      ("run_fem_oi_free.log", "同化なし(no-DA)",  "tab:orange",":")]
PAT=re.compile(r"analysis t=([\d.]+)s.*?Q=([-\d.]+)")


def parse(logname):
    p=os.path.join(OFDIR, logname)
    if not os.path.exists(p): return None
    t=[0.0]; q=[Q0]
    for line in open(p, errors="ignore"):
        m=PAT.search(line)
        if m: t.append(float(m.group(1))); q.append(float(m.group(2)))
    return (np.array(t), np.array(q)) if len(t)>1 else None


def main():
    fig,ax=plt.subplots(figsize=(11,6))
    ax.axvspan(0,300,color="orange",alpha=.06)
    ax.axhline(QTRUE,color="k",lw=2.2,label=f"真値 Q={QTRUE:.0f} W")
    ax.axhline(Q0,color="gray",lw=1.2,ls=":",alpha=.7)
    ax.text(605,Q0,"でたらめ初期 22W",fontsize=9,color="gray",va="center")
    nplot=0
    for logname,lab,col,ls in RUNS:
        r=parse(logname)
        if r is None: continue
        t,q=r; ax.plot(t,q,ls,color=col,lw=2.2,marker="o",ms=4,label=lab); nplot+=1
    ax.set_xlabel("time [s]"); ax.set_ylabel("発熱量 Q の推定 [W]")
    ax.set_ylim(0,24); ax.grid(alpha=.3); ax.legend(fontsize=10,loc="center right")
    ax.set_title("実ソルバOIの発熱量Qの推定：決め打ち σ_Q でどう変わるか\n"
                 "小さいσ_Qは補正が弱くQが動かない／大きいσ_Qは飽和して同じ値へ。"
                 "同化なしは22Wのまま。いずれも真値15Wには届きにくい（弱可同定）。",fontsize=12,weight="bold")
    fig.tight_layout()
    out=os.path.join(IMG,"oi_Q_estimate.png"); fig.savefig(out,dpi=140); plt.close(fig)
    try: shutil.copy(out, os.path.join(ROOT,"..","106_0_pod_selected_rom","docs","img","oi_Q_estimate.png")); print("[Q] copied to 106")
    except Exception as e: print("[Q] copy skip:",e)
    print(f"[Q] wrote oi_Q_estimate.png ({nplot} runs)")


if __name__=="__main__": main()
