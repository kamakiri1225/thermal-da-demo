"""「PODでモード分解したこと」と「どの点を測るか」の関係を1枚で説明する.

POD は全場を5個の係数 a に集約する。だから測定の目的は a を当てることで、
1点測るごとに a についての式が1本立つ。5点で 5x5 の連立方程式になり、
その係数行列 U_{r,P}（モード行列から測る5行を抜いたもの）が
良条件かどうかで復元精度が決まる ― これが点の選び方の正体。

出力: docs/img/pod_to_points.png, results/pod_to_points.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_pod_to_points_fig.py
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
R=5; NRAND=200; KC=273.15


def matrix_panel(ax, A, title, sub, cmap="RdBu_r"):
    """5x5行列を数値つきで描く"""
    v=np.abs(A).max()
    ax.imshow(A, cmap=cmap, vmin=-v, vmax=v, aspect="equal")
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            ax.text(j, i, f"{A[i,j]*1000:+.1f}", ha="center", va="center",
                    fontsize=9.2, weight="bold",
                    color="white" if abs(A[i,j])>0.6*v else "#1b2430")
    ax.set_xticks(range(5)); ax.set_yticks(range(5))
    ax.set_xticklabels([f"$\\varphi_{k+1}$" for k in range(5)], fontsize=11)
    ax.set_yticklabels([f"{k+1}点目" for k in range(5)], fontsize=10)
    ax.set_title(title+"\n"+sub, fontsize=12, weight="bold")


def main():
    ts, Cc, X = q.load_snapshots()
    mean=X.mean(axis=1); Xc=X-mean[:,None]
    U,S,_=np.linalg.svd(Xc, full_matrices=False); Ur=U[:,:R]
    pod=list(qr(Ur.T, pivoting=True)[2][:R])

    def err(p):
        a=np.linalg.pinv(Ur[p,:])@Xc[p,:]
        return float(np.sqrt(((Ur@a-Xc)**2).mean()))
    Aq=Ur[pod,:]; cq=float(np.linalg.cond(Aq)); eq=err(pod)

    rng=np.random.default_rng(0); conds=[]; errs=[]; worst=None
    for _ in range(NRAND):
        p=list(rng.choice(len(Cc), R, replace=False))
        A=Ur[p,:]; c=float(np.linalg.cond(A)); e=err(p)
        conds.append(c); errs.append(e)
        if worst is None or e>worst[1]: worst=(p,e,c)
    pw,ew,cw=worst; Aw=Ur[pw,:]
    print(f"Q-DEIM cond={cq:.2f} err={eq:.5f} K / ランダム最悪 cond={cw:.0f} err={ew:.3f} K")

    fig=plt.figure(figsize=(18.6,9.4))
    gs=fig.add_gridspec(2,4,height_ratios=[1,1.12],hspace=.46,wspace=.40)

    # ① PODの結論
    ax=fig.add_subplot(gs[0,0:2]); ax.axis("off")
    ax.set_title("① PODが言っていること：全場は「5個の数」で決まる",
                 fontsize=13.5, weight="bold", loc="left")
    ax.text(0.01,0.70,
        r"$\hat T(x,t)\;=\;\bar T(x)\;+\;a_1(t)\,\varphi_1(x)+a_2(t)\,\varphi_2(x)"
        r"+\cdots+a_5(t)\,\varphi_5(x)$", fontsize=17, transform=ax.transAxes)
    ax.text(0.01,0.44,
        "形 $\\varphi_k$（モード）は計算済みで固定。\n時刻ごとに変わるのは $a_1\\ldots a_5$ の 5個だけ。",
        fontsize=13, transform=ax.transAxes)
    ax.text(0.01,0.235,
        "→ 20,696セルの温度を知ることは、5個の数 $a$ を当てることと同じになった。",
        fontsize=13.5, transform=ax.transAxes, color="#C0392B", weight="bold")
    ax.text(0.01,0.045,
        "（$\\varphi_k$＝モード形状、$a_k$＝その形がどれだけ混ざっているかの割合）",
        fontsize=11.5, transform=ax.transAxes, color="#556")

    # ② 5点＝5本の式
    ax=fig.add_subplot(gs[0,2:4]); ax.axis("off")
    ax.set_title("② だから「1点測る」＝「$a$ についての式が1本立つ」",
                 fontsize=13.5, weight="bold", loc="left")
    ax.text(0.01,0.70,
        r"$T(P_i)-\bar T(P_i)\;=\;a_1\varphi_1(P_i)+\cdots+a_5\varphi_5(P_i)$",
        fontsize=16, transform=ax.transAxes)
    ax.text(0.01,0.45,
        r"5点そろえると  $\;U_{r,P}\,a=T_P-\bar T_P\;$  という $5\times5$ の連立方程式",
        fontsize=14, transform=ax.transAxes)
    ax.text(0.01,0.235,
        "→ $U_{r,P}$ ＝ モード行列から「測る5行」を抜いたもの。\n"
        "　 どの行を抜くか ＝ どこに測定点を置くか そのもの。",
        fontsize=13.5, transform=ax.transAxes, color="#C0392B", weight="bold")
    ax.text(0.01,0.045,
        "（$U_r$ は 20,696行×5列。そこから5行だけ取り出す）",
        fontsize=11.5, transform=ax.transAxes, color="#556")

    # ③ 良い5点 / 悪い5点 の行列
    matrix_panel(fig.add_subplot(gs[1,0]), Aq,
                 "③ Q-DEIMが選んだ5点の $U_{r,P}$",
                 f"行が互いに向きを変えている\n条件数 {cq:.2f}　→ 復元誤差 {eq*1000:.2f} mK")
    matrix_panel(fig.add_subplot(gs[1,1]), Aw,
                 "④ ランダム5点（最悪）の $U_{r,P}$",
                 f"上3行がそっくり＝同じ式を3回書いている\n"
                 f"条件数 {cw:,.0f}　→ 復元誤差 {ew*1000:.0f} mK")

    # ⑤ 条件数と誤差の関係
    ax=fig.add_subplot(gs[1,2:4])
    ax.scatter(conds, np.array(errs)*1000, s=26, color="#e58f2a", alpha=.65,
               edgecolor="none", label=f"ランダム5点（{NRAND}通り）")
    ax.scatter([cq],[eq*1000], s=210, marker="*", color="#C0392B",
               zorder=5, label=f"Q-DEIM（条件数 {cq:.1f}）")
    ax.scatter([cw],[ew*1000], s=110, marker="X", color="#555",
               zorder=5, label="ランダム最悪")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("連立方程式の条件数 $\\mathrm{cond}(U_{r,P})$（対数軸）", fontsize=12.5)
    ax.set_ylabel("全20,696セルの復元誤差 [mK]（対数軸）", fontsize=12.5)
    ax.set_title("⑤ 点の選び方の良し悪し＝連立方程式の条件数\n"
                 "条件数が悪い点の組ほど復元が崩れる（きれいな比例関係）",
                 fontsize=13, weight="bold")
    ax.legend(fontsize=11, loc="upper left"); ax.grid(alpha=.3, which="both")
    r=np.corrcoef(np.log10(conds), np.log10(errs))[0,1]
    ax.text(.97,.06, f"対数どうしの相関 {r:+.3f}", transform=ax.transAxes,
            ha="right", fontsize=12, weight="bold", color="#1F4E9C",
            bbox=dict(fc="white", ec="#b9c6d6"))

    fig.suptitle("PODでモード分解したことと、どの点を測るかの関係\n"
                 "― 測定点選びは「$5\\times5$ の連立方程式を解きやすくする」問題になっている",
                 fontsize=15.5, weight="bold")
    fig.text(0.5,0.018,
        "Q-DEIM ＝ この連立が解きやすくなるように「測る5行」を選ぶ手続き（オフラインで1回。中身は枢軸付きQR分解）"
        "　／　gappy-POD ＝ 選んだ5点の測定値から $a$ を解いて全場に伸ばす計算（オンラインで毎時刻）",
        ha="center", fontsize=12.5, weight="bold", color="#1F4E9C",
        bbox=dict(fc="#f1f5fa", ec="#b9c6d6", pad=7))
    fig.tight_layout(rect=[0,0.052,1,0.905])
    out=os.path.join(IMG,"pod_to_points.png")
    fig.savefig(out, dpi=125); plt.close(fig)
    with open(os.path.join(RES,"pod_to_points.json"),"w") as f:
        json.dump(dict(qdeim_cells=[int(c) for c in pod], qdeim_cond=cq, qdeim_rmse_K=eq,
                       worst_cells=[int(c) for c in pw], worst_cond=cw, worst_rmse_K=ew,
                       log_corr=float(r), n_random=NRAND), f, ensure_ascii=False, indent=2)
    print("wrote", out)


if __name__=="__main__":
    main()
