"""熱感度 W = K_s^{-1} H_T がどこから来るのかを、導出の順に1枚で示す.

① 熱弾性のFEM式        K_s u = f_th
② 熱荷重は温度に線形    f_th = H_T (T - T_ref)
③ ゆえに               u = K_s^{-1} H_T (T - T_ref) = W (T - T_ref)
④ W の読み方           行＝どこがよく動くか / 列＝どの温度が効くか（実データの分布）
⑤ 本研究での実装        全W(5040x5040)を作らず、PODモード方向だけFrontISTRに通す

出力: docs/img/sensitivity_derivation.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_sensitivity_derivation_fig.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img")
S105=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results")


def main():
    kv=np.load(os.path.join(S105,"kinvh_sensitivity.npz"), allow_pickle=True)
    rs=kv["row_sens"]; valid=kv["valid"]; hi=kv["hi"]; lo=kv["lo"]
    rsv=rs[valid]
    print(f"行感度 |W_i| 有効{valid.sum()}節点: {rsv.min():.3f} 〜 {rsv.max():.3f} µm/K")
    print(f"  高W上位2 = {rs[hi]} / 低W下位2 = {rs[lo]}  比 {rs[hi][0]/rs[lo][0]:.0f}倍")

    fig=plt.figure(figsize=(18.4,9.6))
    gs=fig.add_gridspec(2,3,height_ratios=[1,1.02],hspace=.40,wspace=.32)

    # ── 上段：導出 ──
    ax=fig.add_subplot(gs[0,:]); ax.axis("off")
    ax.text(0.005,0.90,"① 熱弾性を有限要素で離散化すると、ふつうの構造解析と同じ形になる",
            fontsize=14, weight="bold", color="#1F4E9C", transform=ax.transAxes)
    ax.text(0.03,0.735, r"$K_s\,u \;=\; f_{\mathrm{th}}$", fontsize=19, transform=ax.transAxes)
    ax.text(0.20,0.745,
            "$K_s$＝剛性行列（底面固定などの拘束を含む）、$u$＝節点変位、"
            "$f_{\\mathrm{th}}$＝温度が生む等価節点力",
            fontsize=12.5, transform=ax.transAxes, color="#333")

    ax.text(0.005,0.575,"② その熱荷重は、温度について一次（線形）である",
            fontsize=14, weight="bold", color="#1F4E9C", transform=ax.transAxes)
    ax.text(0.03,0.395,
            r"$f_{\mathrm{th}}=\sum_e\int_{V_e} B_e^{T} D_e\,\alpha_e\,"
            r"(T-T_{\mathrm{ref}})\,m\;dV \;\equiv\; H_T\,(T-T_{\mathrm{ref}})$",
            fontsize=17, transform=ax.transAxes)
    ax.text(0.03,0.275,
            "熱ひずみ $\\varepsilon_{\\mathrm{th}}=\\alpha(T-T_{\\mathrm{ref}})\\,m$ "
            "（$m=(1,1,1,0,0,0)^{T}$）が温度に比例するので、"
            "$T$ を括り出して行列 $H_T$ にまとめられる",
            fontsize=12.5, transform=ax.transAxes, color="#333")

    ax.text(0.005,0.145,"③ 2つを合わせると、温度から変位への「掛け算1回」になる",
            fontsize=14, weight="bold", color="#C0392B", transform=ax.transAxes)
    ax.text(0.03,-0.025,
            r"$u \;=\; K_s^{-1}H_T\,(T-T_{\mathrm{ref}})\;\equiv\;W\,(T-T_{\mathrm{ref}})"
            r"\qquad\Longrightarrow\qquad W=\dfrac{\partial u}{\partial T}=K_s^{-1}H_T"
            r"\;\;[\mu\mathrm{m/K}]$",
            fontsize=17, transform=ax.transAxes, color="#C0392B")

    # ── 下段左：Wの読み方 ──
    ax=fig.add_subplot(gs[1,0]); ax.axis("off")
    ax.set_title("④ $W$ の読み方 ― 行と列で用途が違う", fontsize=13.5, weight="bold")
    ax.text(0.02,0.80,"行 $W_{i,:}$（1本の変位に着目）",
            fontsize=13, weight="bold", color="#C0392B", transform=ax.transAxes)
    ax.text(0.02,0.62,"「その点は全体の温度変化に対してどれだけ動くか」\n"
                      "→ 行ノルムが大きい点に 変位センサ を置く",
            fontsize=12, transform=ax.transAxes)
    ax.text(0.02,0.40,"列 $W_{:,j}$（1つの節点温度に着目）",
            fontsize=13, weight="bold", color="#1F4E9C", transform=ax.transAxes)
    ax.text(0.02,0.22,"「その場所が1 K上がると、全体がどう変形するか」\n"
                      "→ 変形を支配する温度が分かり、温度センサ の根拠になる",
            fontsize=12, transform=ax.transAxes)
    ax.text(0.02,0.03,"$K_s$ が拘束を含むので、固定面の近くは $K_s^{-1}$ が小さい\n"
                      "＝ 底面は測っても動かず、情報が少ない",
            fontsize=11.5, transform=ax.transAxes, color="#556")

    # ── 下段中：行感度の実データ分布 ──
    ax=fig.add_subplot(gs[1,1])
    ax.hist(rsv, bins=60, color="#3b74b8", alpha=.82)
    ax.axvline(rs[hi][0], color="#C0392B", lw=2.6, label=f"高W点 {rs[hi][0]:.2f} µm/K")
    ax.axvline(rs[lo][0], color="#7a8899", lw=2.6, ls="--",
               label=f"低W点 {rs[lo][0]:.2f} µm/K")
    ax.set_yscale("log")
    ax.set_xlabel("行感度 $\\|W_{i,:}\\|$ [µm/K]", fontsize=12)
    ax.set_ylabel("節点数", fontsize=12)
    ax.set_title(f"実データ：FrontISTRの $K,H$ から作った $W$\n"
                 f"高Wと低Wで {rs[hi][0]/rs[lo][0]:.0f}倍 の開き",
                 fontsize=13, weight="bold")
    ax.legend(fontsize=11); ax.grid(alpha=.3, which="both")

    # ── 下段右：実装 ──
    ax=fig.add_subplot(gs[1,2]); ax.axis("off")
    ax.set_title("⑤ 本研究での作り方（$W$ 全体は作らない）", fontsize=13.5, weight="bold")
    ax.text(0.02,0.84,"$W$ は 5,040×5,040。全部作るのは無駄。",
            fontsize=12.5, transform=ax.transAxes)
    ax.text(0.02,0.70,"温度場は $T=\\bar T+\\sum_k a_k\\varphi_k$ と5個の数で書けるので、",
            fontsize=12.5, transform=ax.transAxes)
    ax.text(0.02,0.515,
            r"$u \;=\; W\bar T\;+\;\sum_{k=1}^{5} a_k\,(W\varphi_k)$",
            fontsize=17, transform=ax.transAxes)
    ax.text(0.02,0.405,
            "$W\\bar T \\equiv u_{\\mathrm{mean}}$（平均場の変位）、"
            "$W\\varphi_k \\equiv D_{\\mathrm{mode}}[:,k]$（モードkの変位応答）",
            fontsize=11.5, transform=ax.transAxes, color="#556")
    ax.text(0.02,0.215,"→ FrontISTRを 1+5＝6回 呼んで $W\\varphi_k$ を作るだけ。\n"
                       "　 以後は $u=u_{\\mathrm{mean}}+D_{\\mathrm{mode}}\\,a$ の行列積1回",
            fontsize=12.5, transform=ax.transAxes, color="#C0392B", weight="bold")
    ax.text(0.02,0.03,"実装: run/run_da_compare.py の build_disp_operator()\n"
                      "同化のたびにFEMを解き直す必要がなくなる",
            fontsize=11.5, transform=ax.transAxes, color="#556")

    fig.suptitle("熱感度 $W=K_s^{-1}H_T$ はどこから来るのか ― 熱弾性の式からの導出",
                 fontsize=16, weight="bold")
    fig.tight_layout(rect=[0,0,1,0.935])
    out=os.path.join(IMG,"sensitivity_derivation.png")
    fig.savefig(out, dpi=125); plt.close(fig)
    print("wrote", out)


if __name__=="__main__":
    main()
