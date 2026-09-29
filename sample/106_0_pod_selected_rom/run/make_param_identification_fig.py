"""発熱量Qと放熱hをデータ同化で同定する手順と、その結果を1枚で示す.

手順は「状態拡大」：推定したいパラメータを状態ベクトルに入れてしまう。
  z = (T1..T5, Q, h)
温度しか測らなくても、予報を通じて T と Q に相関が育つので、
観測行列 H の Q 成分がゼロでもゲインの Q 行が非ゼロになり Q が動く。

① 手順の流れ
② 温度とQの相関が育つ様子（アンサンブルの散布図）
③ Qの収束（平均±標準偏差、真値15 W）
④ hが決まらない理由（感度の時間形と、応答がノイズ以下であること）

出力: docs/img/param_identification.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_param_identification_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
Q_TRUE_W=15.0


def main():
    d=np.load(os.path.join(RES,"da_history.npz"), allow_pickle=True)
    t=d["time"]; q=d["q"]; h=d["h"]
    s=np.load(os.path.join(RES,"sensitivity_qh.npz"))
    ts=s["t"]; SQ=np.abs(s["S_Q"]).mean(1); Sh=np.abs(s["S_h"]).mean(1)
    hs=json.load(open(os.path.join(RES,"h_sensitivity_check.json")))
    w=np.load(os.path.join(RES,"qh_time_window.npz"))
    h_true=hs["h_true_W_K"]

    fig=plt.figure(figsize=(19.0,9.2))
    gs=fig.add_gridspec(2,3,height_ratios=[.92,1.0],hspace=.40,wspace=.30)

    # ① 手順
    ax=fig.add_subplot(gs[0,0]); ax.axis("off")
    ax.set_title("① 手順は「状態拡大」だけ", fontsize=13.5, weight="bold")
    steps=[("1","推定したいものを状態に入れる\n$z=(T_1,\\dots,T_5,\;Q,\;h)^{T}$（7次元）"),
           ("2","各メンバーに違う $Q,h$ を持たせて\nROMで前進（$Q,h$ は時間で動かさない）"),
           ("3","温度を観測し、EnKFで $z$ を丸ごと更新\n$z^a=z^b+K\\,(y-Hz^b)$"),
           ("4","$H$ の $Q$ 成分は 0 でも、$K$ の $Q$ 行は\n非ゼロ → $Q$ が動く")]
    for i,(n,txt) in enumerate(steps):
        y=0.80-i*0.215
        ax.add_patch(plt.Rectangle((0.01,y-0.10),0.97,0.155,transform=ax.transAxes,
                     fc="#eef3f9",ec="#2E75D4",lw=1.6,clip_on=False))
        ax.text(0.045,y+0.012,n,transform=ax.transAxes,fontsize=15,weight="bold",color="#1F4E9C")
        ax.text(0.13,y+0.012,txt,transform=ax.transAxes,fontsize=11.2,va="center")

    # ② なぜ温度だけでQが動くのか
    ax=fig.add_subplot(gs[0,1]); ax.axis("off")
    ax.set_title("② なぜ温度だけ測って $Q$ が当たるのか", fontsize=13.5, weight="bold")
    ax.text(0.02,0.86,"前進させると「$Q$ が大きいメンバーほど温度が高い」\n"
                      "という関係が自然にできる ＝ 共分散の非対角が育つ:",
            fontsize=11.8, transform=ax.transAxes)
    ax.text(0.05,0.63,r"$P\simeq\dfrac{dZ\,dZ^{T}}{N-1}\;,\qquad P_{TQ}\neq0$",
            fontsize=16, transform=ax.transAxes)
    ax.text(0.02,0.47,"ゲインの $Q$ 行は",fontsize=11.8,transform=ax.transAxes)
    ax.text(0.05,0.285,r"$K_Q=\dfrac{P_{T_1Q}}{P_{T_1T_1}+r}\;(\neq0)$",
            fontsize=16, transform=ax.transAxes, color="#C0392B")
    ax.text(0.02,0.045,"→ 温度の残差 $y-Hz^b$ が $Q$ を引っぱる。\n"
                      "　 OIはこの $P_{TQ}$ を決め打ちするので頭打ちになる（課題①）",
            fontsize=11.5, transform=ax.transAxes, color="#556")

    # ③ Qの収束
    ax=fig.add_subplot(gs[0,2])
    qm=q[:,0]*Q_TRUE_W; qs=q[:,1]*Q_TRUE_W
    ax.axvspan(0,300,color="orange",alpha=.07)
    ax.axhline(Q_TRUE_W,color="k",lw=2.6,ls="--",label=f"真値 {Q_TRUE_W:.0f} W")
    ax.plot(t,qm,"-o",color="#C0392B",lw=2.6,ms=5,label="推定値（アンサンブル平均）")
    ax.fill_between(t,qm-qs,qm+qs,color="#C0392B",alpha=.16,label="±1標準偏差（不確かさ）")
    ax.set_xlabel("時刻 [s]",fontsize=12); ax.set_ylabel("発熱量 $Q$ [W]",fontsize=12)
    ax.set_title(f"③ $Q$ は当たる：ばらつきが {qs[0]:.1f} → {qs[-1]:.1f} W に縮む\n"
                 "（温度1点のみ観測。2点＋変位なら 14.9 W）",fontsize=12.5,weight="bold")
    ax.legend(fontsize=10.5); ax.grid(alpha=.3)

    # ④ 感度の時間形
    ax=fig.add_subplot(gs[1,0])
    ax.axvspan(0,300,color="orange",alpha=.07)
    ax.plot(ts,SQ/SQ.max(),lw=2.8,color="#C0392B",label="$|\\partial T/\\partial Q|$（規格化）")
    ax.plot(ts,Sh/Sh.max(),lw=2.8,color="#1F4E9C",label="$|\\partial T/\\partial h|$（規格化）")
    ax.set_xlabel("時刻 [s]",fontsize=12); ax.set_ylabel("感度（各々の最大で規格化）",fontsize=12)
    ax.set_title("④ $Q$ と $h$ は「効く時間帯」が違う\n"
                 "$Q$＝加熱中に効く／$h$＝冷えていくほど効く",fontsize=12.5,weight="bold")
    ax.legend(fontsize=11); ax.grid(alpha=.3)

    # ⑤ 観測する時間帯とQ誤差
    ax=fig.add_subplot(gs[1,1])
    labs=["加熱期だけ\n(0–300 s)","冷却期だけ\n(300–600 s)","全区間\n(0–600 s)"]
    vals=[float(w["heat_eQ"]),float(w["cool_eQ"]),float(w["full_eQ"])]
    b=ax.bar(range(3),vals,color=["#e58f2a","#7a8899","#C0392B"])
    for r,v in zip(b,vals):
        ax.text(r.get_x()+r.get_width()/2,v,f"{v:.2f} W",ha="center",va="bottom",
                fontsize=12,weight="bold")
    ax.set_xticks(range(3)); ax.set_xticklabels(labs,fontsize=11)
    ax.set_ylabel("$Q$ の推定誤差 [W]",fontsize=12)
    ax.set_title("⑤ $Q$ は「いつ測るか」が効く\n加熱期を外すと 7倍 悪化する",
                 fontsize=12.5,weight="bold")
    ax.grid(alpha=.3,axis="y")

    # ⑥ hが決まらない理由
    ax=fig.add_subplot(gs[1,2])
    dT=hs["max_abs_temperature_difference_K"]; sig=hs["temperature_observation_sigma_K"]
    du=hs["max_abs_displacement_difference_um"]; sigu=hs["displacement_observation_sigma_um"]
    xs=np.arange(2)
    ax.bar(xs-0.19,[dT,du],width=.38,color="#1F4E9C",label="$h$ を40%変えたときの応答")
    ax.bar(xs+0.19,[sig,sigu],width=.38,color="#C0392B",label="観測ノイズ $\\sigma$")
    for x,(a,bb) in enumerate(zip([dT,du],[sig,sigu])):
        ax.text(x-0.19,a,f"{a:.3f}",ha="center",va="bottom",fontsize=11,weight="bold")
        ax.text(x+0.19,bb,f"{bb:.2f}",ha="center",va="bottom",fontsize=11,weight="bold")
    ax.set_xticks(xs); ax.set_xticklabels(["温度 [K]","変位 [µm]"],fontsize=12)
    ax.set_yscale("log")
    ax.set_title(f"⑥ $h$ は決まらない：応答がノイズの 1/6\n"
                 f"推定 {h[-1,0]:.4f} vs 真値 {h_true:.4f} W/K（ばらつきも縮まない）",
                 fontsize=12.5,weight="bold")
    ax.legend(fontsize=10.5); ax.grid(alpha=.3,axis="y",which="both")

    fig.suptitle("データ同化でパラメータを同定する ― 発熱量 $Q$ は当たり、放熱 $h$ は当たらない",
                 fontsize=16,weight="bold")
    fig.text(0.5,0.012,
             "※ 同じ枠組みで $h$ も状態に入れているが、$h$ を40%動かしても600秒での温度差は 0.049 K しかなく、"
             "観測ノイズ 0.30 K に埋もれる。\n"
             "  「同化の設定が悪い」のではなく「この実験条件では $h$ の情報が観測に現れない」という可同定性の問題。",
             ha="center",fontsize=11.5,color="#444")
    fig.tight_layout(rect=[0,0.030,1,0.935])
    out=os.path.join(IMG,"param_identification.png")
    fig.savefig(out,dpi=125); plt.close(fig)
    print(f"Q: {qm[-1]:.2f} W (真値15) / h: {h[-1,0]:.5f} (真値{h_true:.5f})")
    print("wrote",out)


if __name__=="__main__":
    main()
