"""証拠図：温度がそれらしく合っても変位は合うとは限らない（W が温度誤差をµmに増幅）.

FrontISTRの変位オペレータ Dmode（温度→A/O上面Uz）から、QoI=Uz(A)-Uz(O) の
温度に対する感度 w[µm/K] を出し、
 (左) 各ROM点の温度が1KずれるとQoIが何µm動くか（感度ベクトル）
 (右) 温度誤差RMS[K] → QoI誤差RMS[µm] の関係（傾き=‖w‖）と、
      実測のデータ同化点（温度だけ/温度+変位）・QoI信号スケールを重ねる
を示す。

出力: docs/img/disp_amplification.png
再現: OMP_NUM_THREADS=4 python3 run/make_disp_amplification_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    do=np.load(os.path.join(RES,"disp_operator.npz"))
    Dmode=do["Dmode"]                                  # (2,r): A行, O行
    UPp=np.linalg.pinv(U[pod,:])                        # (r,5)
    w=(Dmode[0]-Dmode[1])@UPp                           # (5,) µm/K  QoIの点温度感度
    wn=float(np.linalg.norm(w))
    # 温度誤差RMS→QoI誤差RMS（独立誤差のモンテカルロ）
    rng=np.random.default_rng(0)
    sig=np.linspace(0,1.0,21); qoi_rms=[]
    for s in sig:
        dT=rng.normal(0,max(s,1e-9),(8000,5)); qoi_rms.append((dT@w).std())
    qoi_rms=np.array(qoi_rms)
    QOI_PEAK=2.7                                        # µm（02の変位追従図のピーク）
    # 実測のデータ同化点（02 §4-3 disp_sel_temp_vs_disp より、加熱期平均）
    da_pts={}  # 別観測点の平均絶対誤差は、このQoI標準偏差の図へ混ぜない。
    print("w[µm/K]=",np.round(w,3)," ‖w‖=%.3f µm/K  Σw=%.3f"%(wn,w.sum()))

    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(14,5.6))
    # 左：点ごとの感度
    ax1.bar(range(5),w,color=["tab:red" if x>0 else "tab:blue" for x in w])
    ax1.axhline(0,color="k",lw=1)
    ax1.set_xticks(range(5)); ax1.set_xticklabels([f"P{i}" for i in range(5)])
    ax1.set_ylabel("QoI感度 ∂(Uz_A−Uz_O)/∂T [µm/K]")
    ax1.set_title(f"ROM代表点の温度→全場復元→A/O変位差（‖w‖={wn:.2f} µm/K）",fontsize=12)
    ax1.text(0.03,0.04,"他のROM4点を固定し、1点だけ＋1 K\n節点1個の局所加熱とは異なる",
             transform=ax1.transAxes,fontsize=9.5,color="dimgray",va="bottom")
    ax1.grid(alpha=0.3,axis="y")
    # 右：温度誤差→変位誤差
    ax2.plot(sig,qoi_rms,"-",color="tab:gray",lw=2,label=f"独立温度誤差→QoI誤差（傾き{wn:.2f}µm/K）")
    ax2.axhline(QOI_PEAK,color="green",ls=":",lw=1.5,label=f"QoI信号のピーク ≈ {QOI_PEAK} µm")
    ax2.axvspan(0,0.30,color="orange",alpha=0.08)
    ax2.text(0.02,QOI_PEAK*0.80,"温度観測\nノイズ≈0.3K",fontsize=9,color="darkorange")
    for name,(t,u) in da_pts.items():
        ax2.scatter([t],[u],s=130,zorder=5)
        ax2.annotate(f"{name}\n(T={t}K, uz={u}µm)",(t,u),textcoords="offset points",
                     xytext=(-4,10),ha="right",fontsize=9)
    ax2.set_xlim(-0.03,1.05); ax2.set_ylim(-0.1,2.85)
    ax2.set_xlabel("各ROM点に仮定した独立誤差の標準偏差 [K]"); ax2.set_ylabel("A/O変位差の誤差の標準偏差 [µm]")
    ax2.set_title("独立・平均0・同じ分散の温度誤差を仮定（同化実験ではない）",fontsize=12)
    ax2.grid(alpha=0.3); ax2.legend(fontsize=9,loc="upper left")
    fig.suptitle("ROM温度5点からA/O変位差への感度：gappy-POD＋FrontISTRモード応答",
                 fontsize=14,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.93])
    fig.savefig(os.path.join(IMG,"disp_amplification.png"),dpi=140); plt.close(fig)
    print("[amp] wrote docs/img/disp_amplification.png")


if __name__=="__main__": main()
