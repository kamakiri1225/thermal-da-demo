"""候補センサを「全体の温まりがどれだけ見えるか」「温度差がどれだけ見えるか」の2次元平面に描く（blog_009）.

温度のばらつき B（30〜300 秒の10時刻平均）を固有値分解すると、ほぼ2方向で説明できる：
  e1 ≈ (1,1,1,1,1)/√5  全体の温まり（約 96 %）
  e2                   ヒータ側と反対側の温度差
候補 c の重み h_c（温度 1 K あたりのセンサの値）をこの2方向に射影し、ノイズ 0.3 で割った値を座標にする：
  x_c = √λ1 (h_c·e1) / 0.3,  y_c = √λ2 (h_c·e2) / 0.3   （「信号／ノイズ」。大きいほどその方向がよく見える）
反りの重み w は x≈0（全体の温まりに反応しない）。
センサ1本の σ の減り方は、この平面で w と同じ向きにあり、原点から遠い候補ほど大きい。

出力: docs/img/two_axis_plane.png, results/two_axis_plane.json
再現: OMP_NUM_THREADS=4 python3 run/make_two_axis_plane_fig.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from select_nonrep_temp_disp import priors, NPT, EXCL, A_XYZ, O_XYZ

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
SIG = 0.3
plt.rcParams.update({"font.size": 15})


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); G = U @ UPp
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    nc = lambda x: int(np.linalg.norm(Cc - np.array(x), axis=1).argmin())
    iA, iO = nn(A_XYZ), nn(O_XYZ); w = Wall[iA, 2] - Wall[iO, 2]
    Bs = priors(); Bm = np.mean(Bs, axis=0)
    lam, E = np.linalg.eigh(Bm); order = np.argsort(lam)[::-1]; lam = lam[order]; E = E[:, order]
    e1, e2 = E[:, 0], E[:, 1]
    if e1.sum() < 0: e1 = -e1
    if e2[2] < 0: e2 = -e2                        # ヒータ側（P2）が＋になる向き
    share = lam / lam.sum()
    proj = lambda h: (np.sqrt(lam[0]) * (h @ e1) / SIG, np.sqrt(lam[1]) * (h @ e2) / SIG)
    # 変位計の候補（A・O の近くを除く）、温度計の候補（全セル）
    far = (np.linalg.norm(co - co[iA], axis=1) > EXCL) & (np.linalg.norm(co - co[iO], axis=1) > EXCL)
    Wf = Wall.reshape(-1, NPT); uidx = np.where(np.repeat(far, 3))[0]
    xu, yu = proj(Wf[uidx].T) if False else (np.sqrt(lam[0]) * (Wf[uidx] @ e1) / SIG, np.sqrt(lam[1]) * (Wf[uidx] @ e2) / SIG)
    xt, yt = np.sqrt(lam[0]) * (G @ e1) / SIG, np.sqrt(lam[1]) * (G @ e2) / SIG
    named = {"反り A−O の重み w": w, "変位 S1（blog_008、Ux）": Wall[nn((0.0097, -0.0362, 0.1005)), 0],
             "変位 S2（blog_008、Uz）": Wall[nn((-0.0346, 0.0144, 0.1005)), 2],
             "変位2本（blog_009、Ux、鏡写し）": Wall[nn((0.0188, 0.0325, 0.1005)), 0],
             "温度 P2（代表点）": np.eye(5)[2], "温度 P4（代表点）": np.eye(5)[4],
             "温度計1（blog_009、ヒータ側）": G[nc((0.0367, 0.0008, 0.0492))], "温度計2（blog_009、反対側）": G[nc((-0.0329, 0.0016, 0.0516))],
             "Uz(A)": Wall[iA, 2]}
    pts = {k: proj(v) for k, v in named.items()}
    out = {"note": "B の10時刻平均の固有値分解。x＝全体の温まりの信号/ノイズ、y＝温度差の信号/ノイズ（ノイズ 0.3）。",
           "eigen_share": share.round(4).tolist(), "e1": e1.round(3).tolist(), "e2": e2.round(3).tolist(),
           "points": {k: [float(a), float(b)] for k, (a, b) in pts.items()}}
    json.dump(out, open(os.path.join(RES, "two_axis_plane.json"), "w"), ensure_ascii=False, indent=1)
    print("share", share.round(3), "e1", e1.round(3), "e2", e2.round(3))
    for k, v in pts.items(): print(f"  {k:30s} x={v[0]:7.2f}  y={v[1]:7.2f}")

    fig, ax = plt.subplots(figsize=(13, 10))
    ax.scatter(xu, yu, s=6, color="#C7D2E0", alpha=.5, label="変位計の候補（14,754）")
    ax.scatter(xt, yt, s=6, color="#F2C4B0", alpha=.5, label="温度計の候補（20,696 セル）")
    style = {"反り A−O の重み w": ("*", "#14213D", 420), "変位 S1（blog_008、Ux）": ("s", "#2E6FD8", 160),
             "変位 S2（blog_008、Uz）": ("s", "#2E6FD8", 160), "変位2本（blog_009、Ux、鏡写し）": ("D", "#C0392B", 160),
             "温度 P2（代表点）": ("o", "#E07A1F", 160), "温度 P4（代表点）": ("o", "#E07A1F", 160),
             "温度計1（blog_009、ヒータ側）": ("^", "#2E8B57", 170), "温度計2（blog_009、反対側）": ("^", "#2E8B57", 170),
             "Uz(A)": ("P", "#8E44AD", 170)}
    offs = {"温度 P2（代表点）": (-150, 22), "温度計1（blog_009、ヒータ側）": (10, 10), "変位 S1（blog_008、Ux）": (-150, -26),
            "変位2本（blog_009、Ux、鏡写し）": (10, 10), "温度 P4（代表点）": (-120, -26), "温度計2（blog_009、反対側）": (10, 8),
            "変位 S2（blog_008、Uz）": (-60, 12), "Uz(A)": (10, 8), "反り A−O の重み w": (10, 8)}
    for k, (x, y) in pts.items():
        m, c, s = style[k]; ax.scatter(x, y, marker=m, s=s, color=c, edgecolor="k", zorder=5)
        ax.annotate(k, (x, y), xytext=offs.get(k, (8, 8)), textcoords="offset points", fontsize=12, color=c, weight="bold",
                    arrowprops=dict(arrowstyle="-", color=c, lw=0.8) if abs(offs.get(k, (8, 8))[0]) > 60 else None)
    ax.axhline(0, color="#999", lw=1); ax.axvline(0, color="#999", lw=1)
    ax.set_xlabel(f"全体の温まり（固有値の {share[0]*100:.0f} %）がどれだけ見えるか　＝ 信号 ÷ ノイズ")
    ax.set_ylabel(f"ヒータ側⇔反対側の温度差（{share[1]*100:.0f} %）がどれだけ見えるか　＝ 信号 ÷ ノイズ")
    ax.set_title("候補センサを「全体の温まり」と「温度差」の2軸で見る（30〜300 秒の平均）", fontsize=16)
    ax.legend(fontsize=12, loc="upper left"); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(IMG, "two_axis_plane.png"), dpi=100, facecolor="white"); plt.close(fig)
    print("wrote two_axis_plane.png")


if __name__ == "__main__":
    main()
