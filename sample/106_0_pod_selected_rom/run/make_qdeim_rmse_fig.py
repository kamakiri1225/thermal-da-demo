"""Q-DEIMで点を1つずつ足したときの「温度[K]のRMSE」を示す.

選択の進行を測る残差ノルムは復元前のモード空間の無次元量なので、
実際の温度で何Kになるのかを別に出す。k点のときは k モードを使って
全20,696セルを復元し、OpenFOAM真値とのRMSE[K]を求める。
右はランダム5点（200通り）との比較。

出力: docs/img/qdeim_rmse_vs_npoints.png, results/qdeim_rmse_vs_npoints.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_qdeim_rmse_fig.py
"""
from __future__ import annotations
import os, sys, json, importlib.util
import numpy as np
from scipy.linalg import qr
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
spec=importlib.util.spec_from_file_location('q', os.path.join(HERE,'select_points_qdeim.py'))
q=importlib.util.module_from_spec(spec); spec.loader.exec_module(q)
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
KMAX=7; NRAND=200


def main():
    ts, Cc, X = q.load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,_=np.linalg.svd(Xc, full_matrices=False)

    def rmse(Um, pts):
        a=np.linalg.pinv(Um[pts,:])@Xc[pts,:]
        return float(np.sqrt(((Um@a-Xc)**2).mean()))

    ks=list(range(1,KMAX+1)); vals=[]
    for k in ks:
        Uk=U[:,:k]; pk=list(qr(Uk.T, pivoting=True)[2][:k])
        vals.append(rmse(Uk, pk))
        print(f"  {k}点(={k}モード): 全20,696セルの復元RMSE {vals[-1]:.5f} K")

    U5=U[:,:5]; p5=list(qr(U5.T, pivoting=True)[2][:5]); e5=rmse(U5,p5)
    rng=np.random.default_rng(0)
    rs=np.array([rmse(U5, list(rng.choice(len(Cc),5,replace=False))) for _ in range(NRAND)])

    fig,(a0,a1)=plt.subplots(1,2,figsize=(13.4,5.2))
    a0.plot(ks, vals, "-o", color="#1F4E9C", lw=2.6, ms=8)
    for k,v in zip(ks,vals):
        a0.annotate(f"{v*1000:.0f} mK" if v>=0.001 else f"{v*1000:.2f} mK",
                    (k,v), textcoords="offset points", xytext=(0,11),
                    ha="center", fontsize=10.5, weight="bold")
    a0.axhline(0.3, color="#2e9e5b", ls="--", lw=2.2)
    a0.text(KMAX, 0.34, "観測ノイズ 0.3 K", ha="right", fontsize=11,
            color="#2e9e5b", weight="bold")
    a0.scatter([5],[e5], s=260, marker="*", color="#C0392B", zorder=5,
               label=f"本研究の5点  {e5*1000:.2f} mK")
    a0.set_yscale("log"); a0.set_xticks(ks)
    a0.set_xlabel("測る点の数（＝使うモード数）", fontsize=12.5)
    a0.set_ylabel("全20,696セルの復元RMSE [K]（対数軸）", fontsize=12.5)
    a0.set_title(f"点を1つ足すごとに温度[K]で何が起きるか\n"
                 f"1点 {vals[0]:.3f} K → 5点 {vals[4]:.5f} K（{vals[0]/vals[4]:,.0f}分の1）",
                 fontsize=13, weight="bold")
    a0.legend(fontsize=11); a0.grid(alpha=.3, which="both")

    a1.hist(rs, bins=40, color="#e58f2a", alpha=.78, label=f"ランダム5点（{NRAND}通り）")
    a1.axvline(e5, color="#C0392B", lw=3.0, label=f"Q-DEIM 5点  {e5:.5f} K")
    a1.axvline(np.median(rs), color="#555", lw=2.2, ls="--",
               label=f"ランダム中央値  {np.median(rs):.5f} K")
    a1.axvline(rs.max(), color="#7a8899", lw=2.2, ls=":",
               label=f"ランダム最悪  {rs.max():.3f} K")
    a1.set_xscale("log")
    a1.set_xlabel("5点から全20,696セルを復元したRMSE [K]（対数軸）", fontsize=12.5)
    a1.set_ylabel("回数", fontsize=12.5)
    a1.set_title(f"同じ5点でも選び方で {rs.max()/e5:,.0f}倍 変わる\n"
                 "Q-DEIMの価値＝大外れを確実に避けること", fontsize=13, weight="bold")
    a1.legend(fontsize=10.5); a1.grid(alpha=.3)

    fig.suptitle("Q-DEIMの効果を温度[K]で見る（モード空間の残差ではなく、復元してからのRMSE）",
                 fontsize=14.5, weight="bold")
    fig.tight_layout(rect=[0,0,1,0.90])
    out=os.path.join(IMG,"qdeim_rmse_vs_npoints.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    with open(os.path.join(RES,"qdeim_rmse_vs_npoints.json"),"w") as f:
        json.dump(dict(k=ks, rmse_K=vals, qdeim5_K=e5,
                       random_median_K=float(np.median(rs)),
                       random_worst_K=float(rs.max())), f, indent=2)
    print("wrote", out)


if __name__=="__main__":
    main()
