"""blog_006 用：変位観測が「温度計のない点」の温度まで直す様子を示す.

本番計算（温度P2+P0＋変位A/O、seed 20260913）の1回目の更新を再現し、
  ① 変位1本が5点すべての温度に触れること（熱感度 W の行）
  ② 60メンバーの共分散 Cov(T_i, Uz) ＝ どの点がずれていると変位がずれるか
  ③ 平均の補正量を「温度計の分」と「変位計の分」に分けた内訳（t=30 s と t=150 s）
を1枚にする。数値は results/trace_displacement_and_Q.json から読む。

出力: docs/img/blog006_disp_fixes_temp.png
再現: OMP_NUM_THREADS=4 python3 run/make_disp_fixes_temp_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
OBS=[2,0]                       # 温度計を置いた点（P2, P0）


def main():
    o=json.load(open(os.path.join(RES,"trace_displacement_and_Q.json")))
    W=np.array(o["W_um_per_K"]); lab=[f"P{i}" for i in range(5)]
    fig=plt.figure(figsize=(15.5,8.4))
    gs=fig.add_gridspec(2,2,height_ratios=[1,1.05],hspace=0.42,wspace=0.22)
    # ① 熱感度 W の行
    ax=fig.add_subplot(gs[0,0]); x=np.arange(5); bw=0.38
    ax.bar(x-bw/2,W[0],bw,color="#C0392B",label="変位計 A が受ける影響")
    ax.bar(x+bw/2,W[1],bw,color="#2E6FD8",label="変位計 O が受ける影響")
    for i in range(5):
        ax.text(i-bw/2,W[0,i]+0.02,f"{W[0,i]:+.2f}",ha="center",fontsize=9.5)
        ax.text(i+bw/2,W[1,i]+0.02,f"{W[1,i]:+.2f}",ha="center",fontsize=9.5)
    ax.axhline(0,color="k",lw=.8); ax.set_xticks(x)
    ax.set_xticklabels([f"{l}\n{'温度計あり' if i in OBS else '温度計なし'}" for i,l in enumerate(lab)],fontsize=10)
    ax.set_ylabel("熱感度 W [µm/K]"); ax.set_ylim(-0.45,0.95); ax.legend(fontsize=9.5,loc="upper left")
    ax.set_title("① 変位計1本が、5点すべての温度に触れている\n（その点の温度が1 K上がると、変位が何µm動くか）",fontsize=12)
    ax.grid(axis="y",alpha=.3)
    # ② 共分散（1サイクル目）
    c0=o["cycles"][0]; cov=np.array(c0["Cov_q_Y"])
    ax=fig.add_subplot(gs[0,1])
    ax.text(0.5,0.5,"",transform=ax.transAxes)
    # 60メンバーの Cov(T_i, Uz(A)) は trace に無いので、W と温度のばらつきで説明する図にする
    sd=np.array(json.load(open(os.path.join(RES,"sensor_pair_goal_oriented.json")))["residual_150s"]["none"]["sd"])
    ax.bar(x,np.abs(W[0])*sd,color="#7D3C98")
    for i in range(5): ax.text(i,np.abs(W[0,i])*sd[i]+0.02,f"{np.abs(W[0,i])*sd[i]:.2f}",ha="center",fontsize=9.5)
    ax.set_xticks(x); ax.set_xticklabels(lab,fontsize=10.5)
    ax.set_ylabel("|W| × 温度のばらつき [µm]")
    ax.set_title("② どの点の温度がずれていると、変位がずれるか\n（熱感度 × まだ分かっていない温度の幅、150 s）",fontsize=12)
    ax.grid(axis="y",alpha=.3)
    # ③ 補正の内訳
    for j,ci in enumerate([0,1]):
        c=o["cycles"][ci]; K=np.array(c["K"]); r=np.array(c["residual"])
        t_part=(K[:,:2]@r[:2])[:5]; u_part=(K[:,2:]@r[2:])[:5]
        ax=fig.add_subplot(gs[1,j])
        ax.bar(x-bw/2,t_part,bw,color="#E67E22",label="温度計2本による補正")
        ax.bar(x+bw/2,u_part,bw,color="#2E8B57",label="変位計2本による補正")
        for i in range(5):
            ax.text(i-bw/2,t_part[i],f"{t_part[i]:+.2f}",ha="center",va="bottom" if t_part[i]>=0 else "top",fontsize=9)
            ax.text(i+bw/2,u_part[i],f"{u_part[i]:+.2f}",ha="center",va="bottom" if u_part[i]>=0 else "top",fontsize=9)
        ax.axhline(0,color="k",lw=.8); ax.set_xticks(x)
        ax.set_xticklabels([f"{l}\n{'温度計あり' if i in OBS else '温度計なし'}" for i,l in enumerate(lab)],fontsize=10)
        ax.set_ylabel("平均温度の補正量 [K]")
        sub="温度計のずれが大きい時期" if ci==0 else "温度計のずれがほぼ0になった後"
        ax.set_title(f"③ t={c['t_s']:.0f} s の補正の内訳（{sub}）",fontsize=12)
        ax.legend(fontsize=9.5); ax.grid(axis="y",alpha=.3)
        if ci==1: ax.set_ylim(-0.16,0.08)
    fig.suptitle("変位観測が温度を直す仕組み（本番計算・seed 20260913。温度計は P2・P0 の2本だけ）",fontsize=13.5,weight="bold")
    fig.savefig(os.path.join(IMG,"blog006_disp_fixes_temp.png"),dpi=150,bbox_inches="tight",pad_inches=0.15); plt.close(fig)
    print("wrote blog006_disp_fixes_temp.png")


if __name__=="__main__": main()
