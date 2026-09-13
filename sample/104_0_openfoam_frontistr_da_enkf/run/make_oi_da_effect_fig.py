"""実ソルバOIの「データ同化の効果」図：温度・変位の時刻歴を
   正解 / 同化なし(free-run) / OI同化後 / 補正前の予報値 で比較する.

必要な入力（実ソルバ2本）:
  - OI本体（補正前予報も記録した拡張版）: results/oi_fullsolver{SUF}_history.npz
  - 同化なし free-run: results/oi_fullsolver_free_history.npz（OI_FREE=1で回したもの）

出力: docs/img/oi_da_effect.png（106へもコピー）
再現: OMP_NUM_THREADS=4 python3 run/make_oi_da_effect_fig.py [oi_tag]
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


def load(suf):
    p=os.path.join(RES, f"oi_fullsolver{('_'+suf) if suf else ''}_history.npz")
    return dict(np.load(p, allow_pickle=True)) if os.path.exists(p) else None


def main():
    oi_tag=sys.argv[1] if len(sys.argv)>1 else "b"
    oi=load(oi_tag); free=load("free")
    if oi is None or free is None:
        print(f"[da] まだ揃っていない: OI(tag={oi_tag})={oi is not None}, free={free is not None}"); return
    t=oi["times"]; tnames=["hot","cold"]; dnames=["A：ヒータ側上面 Uz(+X)","O：反対側上面 Uz(−X)"]

    fig,axes=plt.subplots(2,2,figsize=(14,9))
    # 温度2点
    for j in range(2):
        ax=axes[0][j]; ax.axvspan(0,300,color="orange",alpha=.06)
        ax.plot(t, oi["T_obs_truth"][:,j]-K,"-",color="k",lw=3.5,alpha=.45,label="正解(真値)")
        ax.plot(t, free["T_obs_analysis"][:,j]-K,":",color="tab:orange",lw=2.2,label="同化なし(free-run)")
        ax.plot(t, oi["T_obs_forecast"][:,j]-K,"-.",color="tab:gray",lw=1.6,label="補正前の予報値")
        ax.plot(t, oi["T_obs_analysis"][:,j]-K,"--o",color="tab:blue",ms=3,label="OI同化後")
        ax.set_title(f"温度時刻歴：観測点 {tnames[j]}"); ax.set_ylabel("温度 [degC]")
        ax.grid(alpha=.3); ax.legend(fontsize=8)
    # 変位2点
    for j in range(2):
        ax=axes[1][j]; ax.axvspan(0,300,color="orange",alpha=.06)
        ax.plot(t, oi["u_truth"][:,j]*1000,"-",color="k",lw=3.5,alpha=.45,label="正解(真値)")
        ax.plot(t, free["u_analysis"][:,j]*1000,":",color="tab:orange",lw=2.2,label="同化なし(free-run)")
        ax.plot(t, oi["u_forecast"][:,j]*1000,"-.",color="tab:gray",lw=1.6,label="補正前の予報値")
        ax.plot(t, oi["u_analysis"][:,j]*1000,"--o",color="tab:blue",ms=3,label="OI同化後")
        ax.set_title(f"変位時刻歴：{dnames[j]}"); ax.set_ylabel("Uz [µm]"); ax.set_xlabel("time [s]")
        ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.suptitle("実ソルバOIのデータ同化効果：正解／同化なし／補正前の予報値／OI同化後\n"
                 "（同化なしは初期の外れが残り、OIは各サイクルで予報→補正で真値へ）",fontsize=14,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.95])
    out=os.path.join(IMG,"oi_da_effect.png"); fig.savefig(out,dpi=140); plt.close(fig)
    dst=os.path.join(ROOT,"..","106_0_pod_selected_rom","docs","img","oi_da_effect.png")
    try: shutil.copy(out,dst); print("[da] copied to 106/docs/img")
    except Exception as e: print("[da] copy skip:",e)
    print("[da] wrote docs/img/oi_da_effect.png")


if __name__=="__main__": main()
