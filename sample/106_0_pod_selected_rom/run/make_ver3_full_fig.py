"""検証3（blog_008 §6）の詳細図を、観測していない温度・個別の変位まで広げて描く.

入力: results/unobserved_traj_v_<case>.npz, results/unobserved_check_v_<case>.json
      （run/selected_pair_unobserved_check.py v で作る）
出力: docs/img/ver3_detail_full.png
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
CASES = [("learned", "15 W"), ("q25", "25 W"), ("intermittent", "間欠加熱")]
# (種類, ラベル, 縦軸)
ROWS = [("R", "温度場全体（20,696 セル）の RMSE", "RMSE [K]"),
        ("T", "P4（代表点）", "温度 [°C]"),
        ("T", "反対側・上部（代表点以外）", "温度 [°C]"),
        ("T", "横（+y）・中段（代表点以外）", "温度 [°C]"),
        ("U", "反り A−O", "変位 [µm]"),
        ("U", "Uz(A)", "変位 [µm]"),
        ("U", "Uz(O)", "変位 [µm]"),
        ("U", "F（反対側・中段 z=50）", "変位 [µm]")]


def main():
    fig, axs = plt.subplots(len(ROWS), 3, figsize=(19, 4.0 * len(ROWS)))
    for j, (case, lab) in enumerate(CASES):
        z = np.load(os.path.join(RES, f"unobserved_traj_v_{case}.npz"))
        js = json.load(open(os.path.join(RES, f"unobserved_check_v_{case}.json"), encoding="utf-8"))["configs"]
        t = z["t"]; names = list(z["names"]); cols = list(z["colors"])
        tl = list(z["tlabels"]); ul = list(z["ulabels"])
        tl_find = lambda key: next(i for i, s in enumerate(tl) if s.startswith(key))
        for r, (kind, title, ylab) in enumerate(ROWS):
            ax = axs[r, j]
            for a0, a1 in ([(0, 150), (300, 450)] if case == "intermittent" else [(0, 300)]):
                ax.axvspan(a0, a1, color="#FDEBD0", alpha=.45)
            if kind == "R":
                for k, nm in enumerate(names): ax.plot(t, z["R"][k], color=cols[k], lw=2.3)
                ax.set_yscale("log"); errs = [js[nm]["field_rmse_K"] for nm in names]; unit = "K"
            else:
                if kind == "T":
                    i = tl_find(title.split("（")[0]); tv = z["Ttrue"][:, i] - 273.15
                    ys = [z["T"][k][:, i] - 273.15 for k in range(len(names))]
                    errs = [js[nm]["temperature_K"][tl[i]] for nm in names]; unit = "K"
                else:
                    i = ul.index(title); tv = z["Utrue"][:, i]; ys = [z["U"][k][:, i] for k in range(len(names))]
                    errs = [js[nm]["displacement_um"][title] for nm in names]; unit = "µm"
                ax.plot(t, tv, "o-", color="#111", lw=2.6, ms=4)
                for k in range(len(names)): ax.plot(t, ys[k], color=cols[k], lw=2.2)
            yy = 0.97
            for k, e in enumerate(errs):
                ref = (k == 1 and kind == "U" and title in ("反り A−O", "Uz(A)", "Uz(O)"))
                ax.text(0.98, yy, f"{e:.2f} {unit}" + ("（観測した点）" if ref else ""), transform=ax.transAxes,
                        ha="right", va="top", fontsize=12, color=cols[k], weight="bold",
                        bbox=dict(fc="white", ec="none", alpha=.8, pad=1)); yy -= 0.10
            ax.set_title(f"{lab}　{title}", fontsize=14.5); ax.grid(alpha=.3)
            if j == 0: ax.set_ylabel(ylab, fontsize=13)
            if r == len(ROWS) - 1: ax.set_xlabel("時刻 [s]")
    handles = [Line2D([], [], color="#111", lw=2.6, marker="o", ms=4, label="真値（OpenFOAM＋FrontISTR）")] + \
              [Line2D([], [], color=c, lw=2.5, label=n) for n, c in zip(names, cols)]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.983), ncol=4, fontsize=14, frameon=False)
    fig.suptitle("検証3　真値を OpenFOAM＋FrontISTR にして、観測していない温度・変位を3条件で確かめる（5 seed 平均）",
                 fontsize=18, weight="bold", y=0.997)
    fig.text(0.5, 0.962, "温度計は P2・P0（観測）。グラフの点は、青の「A・O 自身を測る」の A−O・Uz(A)・Uz(O) を除き、観測に使っていない。"
             "数字＝30〜300 秒の平均誤差（線の色と対応）。橙の帯＝ヒータ ON", ha="center", fontsize=13, color="#444")
    fig.tight_layout(rect=(0, 0, 1, 0.952))
    fig.savefig(os.path.join(IMG, "ver3_detail_full.png"), dpi=85, facecolor="white"); plt.close(fig)
    print("wrote ver3_detail_full.png")


if __name__ == "__main__":
    main()
