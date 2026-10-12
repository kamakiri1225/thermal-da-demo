"""熱感度の重み × 温度の時刻歴 → 変位、を時間を追って見せるアニメーション（blog_008 §5-2）.

15 W、真値＝OpenFOAM（5点の温度）＋FrontISTR（全温度場を入れた変位）。
左：代表5点の温度（平均の温度分布からの上がり分）
中：U_z(A) の組み立て。平均の温度分布での変位 u[T̄] に、各点の「重み × 温度の上がり分」を積み上げる
右：時刻歴。重み × 温度（線）と、FrontISTR で全温度場を解いた変位（点）
出力: docs/img/weight_times_temp_anim.gif
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
COL = ["#C0392B", "#2E6FD8", "#E07A1F", "#2E8B57", "#8E44AD"]
plt.rcParams.update({"font.size": 15})


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); mean = kv["mean"].astype(float); pod = kv["cell_idx"]
    UPp = np.linalg.pinv(U[pod, :]); mp = mean[pod]
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]; um = op["u_mean"]
    W = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    iA, iO = nn((0.028, 0, 0.1005)), nn((-0.028, 0, 0.1005))
    hA, hO = W[iA, 2], W[iO, 2]
    z = np.load(os.path.join(RES, "limit_truth_learned.npz")); t = z["times"]; T5 = z["Tfield"][:, pod]
    ut = np.load(os.path.join(RES, "truth_disp_all_learned.npz"))["u"]
    dT = T5 - mp                                         # 平均の温度分布からの上がり分
    contA = dT * hA                                      # 各点の寄与
    estA = um[iA, 2] + contA.sum(1); estO = um[iO, 2] + dT @ hO
    trA = ut[:, iA, 2]; trO = ut[:, iO, 2]

    fig, axs = plt.subplots(1, 3, figsize=(21, 7.2), gridspec_kw=dict(width_ratios=[1, 1.25, 1.35]))
    fig.subplots_adjust(top=0.74, bottom=0.17, wspace=0.32)
    fig.suptitle("熱感度の重み × 温度の時刻歴 ＝ 変位　（15 W、温度は OpenFOAM、点は FrontISTR で全温度場を解いた変位）",
                 fontsize=19, weight="bold")
    sub = fig.text(0.5, 0.865, "", ha="center", fontsize=16, color="#333")
    ax0, ax1, ax2 = axs
    labels = ["P0", "P1", "P2", "P3", "P4"]
    ax2.plot(t, estA, color=COL[0], lw=2.5, label="U_z(A)：重み × 温度")
    ax2.plot(t, trA, "o", color=COL[0], ms=7, mfc="none", label="U_z(A)：FrontISTR")
    ax2.plot(t, estA - estO, color="#14213D", lw=2.5, label="反り A−O：重み × 温度")
    ax2.plot(t, trA - trO, "s", color="#14213D", ms=6, mfc="none", label="反り A−O：FrontISTR")
    ax2.axvspan(0, 300, color="#FDEBD0", alpha=.5)
    cur = ax2.axvline(0, color="#555", lw=2)
    ax2.set_xlabel("時刻 [s]"); ax2.set_ylabel("変位 [µm]"); ax2.grid(alpha=.3); ax2.legend(fontsize=12, loc="upper right", bbox_to_anchor=(1.0, 0.62))
    ax2.set_title("③ 時刻歴：線（掛け算）と点（FrontISTR）が重なる", fontsize=16)
    lo = min(contA.min(), 0) - 0.5; hi = max(estA.max(), um[iA, 2]) + 0.8

    def draw(k):
        ax0.clear(); ax1.clear()
        ax0.bar(labels, dT[k], color=COL)
        for i, v in enumerate(dT[k]): ax0.text(i, v + (0.15 if v >= 0 else -0.45), f"{v:+.2f}", ha="center", fontsize=13)
        ax0.set_ylim(min(dT.min(), 0) - 0.6, dT.max() + 0.6); ax0.axhline(0, color="#999", lw=1)
        ax0.set_ylabel("平均の温度分布からの上がり分 [K]"); ax0.grid(alpha=.3, axis="y")
        ax0.set_title("① 代表5点の温度（OpenFOAM）", fontsize=16)
        # 中：積み上げ（滝グラフ）
        base = um[iA, 2]; x = 0
        ax1.bar(x, base, color="#9AA5B1"); ax1.text(x, base + 0.15, f"{base:.2f}", ha="center", fontsize=12)
        lvl = base
        for i in range(5):
            c = contA[k, i]; x = i + 1
            ax1.bar(x, c, bottom=lvl, color=COL[i])
            ax1.text(x, max(lvl, lvl + c) + 0.12, f"{c:+.2f}", ha="center", fontsize=12)
            lvl += c
        ax1.bar(6, lvl, color="#C0392B", alpha=.85); ax1.text(6, lvl + 0.15, f"{lvl:.2f}", ha="center", fontsize=12, weight="bold")
        ax1.plot(6, trA[k], "o", ms=13, mfc="none", mec="#14213D", mew=2.5)
        ax1.set_xticks(range(7)); ax1.set_xticklabels(["u[T̄]"] + [f"P{i}\n{hA[i]:+.2f}\n×{dT[k, i]:+.2f}" for i in range(5)] + ["合計"], fontsize=11.5)
        ax1.set_ylim(lo, hi); ax1.axhline(0, color="#999", lw=1); ax1.grid(alpha=.3, axis="y")
        ax1.set_ylabel("U_z(A) [µm]")
        ax1.set_title("② U_z(A) の組み立て：重み × 上がり分 を積み上げる（○＝FrontISTR）", fontsize=14)
        cur.set_xdata([t[k], t[k]])
        sub.set_text(f"t = {t[k]:.0f} 秒　　U_z(A)：掛け算 {estA[k]:.3f} µm ／ FrontISTR {trA[k]:.3f} µm　　"
                     f"反り A−O：掛け算 {estA[k]-estO[k]:.3f} µm ／ FrontISTR {trA[k]-trO[k]:.3f} µm")
        return []

    anim = FuncAnimation(fig, draw, frames=len(t), interval=600)
    anim.save(os.path.join(IMG, "weight_times_temp_anim.gif"), writer=PillowWriter(fps=2), dpi=62)
    draw(10); fig.savefig(os.path.join(IMG, "weight_times_temp_t300.png"), dpi=80, facecolor="white")
    plt.close(fig); print("wrote weight_times_temp_anim.gif", os.path.getsize(os.path.join(IMG, "weight_times_temp_anim.gif")) // 1024, "KB")


if __name__ == "__main__":
    main()
