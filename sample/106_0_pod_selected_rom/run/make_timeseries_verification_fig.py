"""「5点で本当に全体が合うのか」を棒グラフではなく温度の時刻歴で見せる.

3つの検証点（いずれも校正・選点に使っていない点）について
  黒  : OpenFOAM(CHT)の真値
  赤  : Q-DEIM5点の値から gappy-POD で復元した温度
  灰  : ランダム5点(200通り中の最悪)から復元した温度
を重ね、右端に全20,696セルRMSEの時刻歴を置く。

出力: docs/img/timeseries_verification.png, results/timeseries_verification.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_timeseries_verification_fig.py
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
KC=273.15; R=5; NRAND=200; SEED=20260930


def recon(U, mean, pts, X):
    return mean[:,None] + U @ (np.linalg.pinv(U[pts,:]) @ (X[pts,:]-mean[pts,None]))


def main():
    ts, Cc, X = q.load_snapshots()
    mean=X.mean(axis=1)
    U,_,_=np.linalg.svd(X-mean[:,None], full_matrices=False); Ur=U[:,:R]
    pod=list(qr(Ur.T, pivoting=True)[2][:R])

    Xq=recon(Ur, mean, pod, X)
    e_q=float(np.sqrt(((Xq-X)**2).mean()))

    rng=np.random.default_rng(SEED); best=None; worst=None; scores=[]
    for _ in range(NRAND):
        p=list(rng.choice(len(X), R, replace=False))
        if np.linalg.matrix_rank(Ur[p,:])<R: continue
        Xr=recon(Ur, mean, p, X); e=float(np.sqrt(((Xr-X)**2).mean()))
        scores.append(e)
        if worst is None or e>worst[0]: worst=(e,p,Xr)
    med=float(np.median(scores))
    print(f"Q-DEIM {e_q:.5f} K / ランダム中央 {med:.5f} K / ランダム最悪 {worst[0]:.4f} K")

    # 検証点：高さを上・中・下に分け、各帯でQ-DEIM5点から最も遠いセルを取る
    # （偏りを避けるため。単に「最遠」を並べると全部が上端に集まってしまう）
    dmin=np.linalg.norm(Cc[:,None,:]-Cc[pod][None,:,:],axis=2).min(axis=1)
    z=Cc[:,2]; picks=[]
    for lo,hi in [(0.067,0.1005),(0.033,0.067),(0.0,0.033)]:
        band=np.where((z>=lo)&(z<hi))[0]
        band=np.array([c for c in band if c not in pod])
        picks.append(int(band[np.argmax(dmin[band])]))

    fig,axes=plt.subplots(1,4,figsize=(19.2,5.4))
    for ax,c in zip(axes[:3], picks):
        ax.axvspan(0,300,color="orange",alpha=.07)
        ax.plot(ts, X[c]-KC, color="black", lw=5.4, alpha=.38, label="OpenFOAM（真値）")
        ax.plot(ts, Xq[c]-KC, "--", color="#C0392B", lw=2.4, dashes=(4,3),
                label="Q-DEIM 5点から復元")
        ax.plot(ts, worst[2][c]-KC, ":", color="#7a8899", lw=2.4,
                label="ランダム5点（最悪）から復元")
        ax.set_title(f"検証点 ({Cc[c][0]*1000:.0f}, {Cc[c][1]*1000:.0f}, {Cc[c][2]*1000:.0f}) mm\n"
                     f"5点から {dmin[c]*1000:.0f} mm　Q-DEIM誤差 "
                     f"{np.sqrt(((Xq[c]-X[c])**2).mean())*1000:.1f} mK",
                     fontsize=12.5, weight="bold")
        ax.set_xlabel("時刻 [s]", fontsize=11.5); ax.set_ylabel("温度 [℃]", fontsize=11.5)
        ax.grid(alpha=.3); ax.tick_params(labelsize=10.5)
    axes[0].legend(fontsize=10.5, loc="lower right")

    ax=axes[3]
    ax.axvspan(0,300,color="orange",alpha=.07)
    ax.plot(ts, np.sqrt(((Xq-X)**2).mean(axis=0)), "-o", color="#C0392B", ms=3.4, lw=2.4,
            label=f"Q-DEIM 5点（全時刻 {e_q*1000:.2f} mK）")
    ax.plot(ts, np.sqrt(((worst[2]-X)**2).mean(axis=0)), ":", color="#7a8899", lw=2.6,
            label=f"ランダム5点 最悪（{worst[0]:.3f} K）")
    ax.axhline(0.3, color="#2e9e5b", ls="--", lw=2)
    ax.text(600, 0.34, "観測ノイズ 0.3 K", fontsize=10.5, color="#2e9e5b",
            ha="right", weight="bold")
    ax.set_yscale("log"); ax.set_xlabel("時刻 [s]", fontsize=11.5)
    ax.set_ylabel("全20,696セルのRMSE [K]", fontsize=11.5)
    ax.set_title("全セルで見た復元誤差の時刻歴\nQ-DEIMは全時刻でノイズの2桁下",
                 fontsize=12.5, weight="bold")
    ax.legend(fontsize=10.5, loc="center right"); ax.grid(alpha=.3, which="both")
    ax.tick_params(labelsize=10.5)

    fig.suptitle("5点だけで全温度場が本当に合うのか ― 上部・中央・下部の検証点を時刻歴で確認する",
                 fontsize=15, weight="bold")
    fig.tight_layout(rect=[0,0,1,0.90])
    out=os.path.join(IMG,"timeseries_verification.png")
    fig.savefig(out, dpi=135); plt.close(fig)
    with open(os.path.join(RES,"timeseries_verification.json"),"w") as f:
        json.dump(dict(qdeim_rmse_K=e_q, random_median_K=med, random_worst_K=worst[0],
                       n_random=len(scores),
                       picks=[dict(cell=int(c), xyz_mm=[round(v*1000,1) for v in Cc[c]],
                                   dist_mm=round(float(dmin[c])*1000,1)) for c in picks]),
                  f, ensure_ascii=False, indent=2)
    print("wrote", out)


if __name__=="__main__":
    main()
