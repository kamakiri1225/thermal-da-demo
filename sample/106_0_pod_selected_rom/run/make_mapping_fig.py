"""OpenFOAM(セル中心) → FrontISTR(節点) の温度マッピングを図解する（blog_001 §2-3用）.

3枚組:
  A: 同じ断面にCFDセル中心とFEM節点を重ねる（位置が一致しないことを示す）
  B: 1節点を拡大し、最近傍k個のセルと逆距離重みを表示
  C: マッピング後のFEM節点温度
出力: docs/img/blog001_mapping.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_mapping_fig.py
"""
from __future__ import annotations
import os, sys, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
SAMPLE = os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE, "102_1_frontistr_hollow_cylinder_thermal_expansion", "python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
IMG = os.path.join(ROOT, "docs", "img")
KC = 273.15; ZSLICE = 0.050; TOL = 0.004; K = 8


def main():
    # --- CFD側: t=300s のセル中心と温度 ---
    Cc = m.read_foam_field(os.path.join(m.OF, "300", "solid", "C"))
    Tc = m.read_foam_field(os.path.join(m.OF, "300", "solid", "T"), len(Cc)) - KC
    # --- FEM側: 節点座標 ---
    mesh = cylinder_mesh.build_cylinder_mesh(4, 48, 20, 0.020, 0.0375, 0.1005)
    nodes = np.array([xyz for _n, xyz in mesh["nodes"]])
    print(f"[map] CFDセル {len(Cc)} 個 / FEM節点 {len(nodes)} 個")

    # --- IDWで全節点へ写像 ---
    Tn = np.zeros(len(nodes))
    for i, p in enumerate(nodes):
        d = np.linalg.norm(Cc - p, axis=1)
        nb = np.argsort(d)[:K]; dn = d[nb]
        if dn[0] < 1e-9:
            Tn[i] = Tc[nb[0]]
        else:
            w = 1.0 / dn
            Tn[i] = (w * Tc[nb]).sum() / w.sum()

    sc = np.abs(Cc[:, 2] - ZSLICE) < TOL          # CFDの断面
    sn = np.abs(nodes[:, 2] - ZSLICE) < TOL * 2   # FEMの断面
    vmin, vmax = Tc[sc].min(), Tc[sc].max()

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6))

    # A: 重ね合わせ
    ax = axes[0]
    ax.scatter(Cc[sc, 0]*1000, Cc[sc, 1]*1000, c=Tc[sc], cmap="turbo",
               s=14, vmin=vmin, vmax=vmax, label=f"CFDセル中心({sc.sum()}個)")
    ax.scatter(nodes[sn, 0]*1000, nodes[sn, 1]*1000, facecolors="none",
               edgecolors="k", s=42, lw=1.1, label=f"FEM節点({sn.sum()}個)")
    ax.set_title("A: 位置が一致しない\n（CFDはセル中心、FEMは節点）", fontsize=12.5, weight="bold")
    ax.legend(fontsize=8.5, loc="upper left", bbox_to_anchor=(0,-0.22), ncol=2, frameon=False); ax.set_aspect("equal")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")

    # B: 1節点の近傍k個と重み（温度差が大きい節点を選び、棒グラフで中身を見せる）
    ax = axes[1]
    cand = np.where(sn)[0]
    spread = []
    for i in cand:
        d = np.linalg.norm(Cc - nodes[i], axis=1); nb = np.argsort(d)[:K]
        spread.append(Tc[nb].max() - Tc[nb].min())
    tgt = int(cand[int(np.argmax(spread))])       # 近傍の温度差が最大の節点
    p = nodes[tgt]
    d = np.linalg.norm(Cc - p, axis=1); nb = np.argsort(d)[:K]
    dn = d[nb] * 1000; w = 1.0 / d[nb]; w /= w.sum()
    ypos = np.arange(K)[::-1]
    bars = ax.barh(ypos, w, color=plt.cm.turbo((Tc[nb] - vmin) / (vmax - vmin)),
                   edgecolor="k", lw=.5)
    for y, wj, dj, tj in zip(ypos, w, dn, Tc[nb]):
        ax.text(wj + 0.004, y, f"{tj:.2f}℃", va="center", fontsize=9.5, weight="bold")
    ax.set_yticks(ypos)
    ax.set_yticklabels([f"{i+1}番目  d={dj:.2f}mm" for i, dj in enumerate(dn)], fontsize=9)
    ax.set_xlabel("重み $w=1/d$ （合計1に正規化）")
    ax.set_xlim(0, w.max() * 1.45)
    ax.set_title(f"B: 近い{K}個のセルから逆距離加重平均\n"
                 f"→ この節点の温度 = {Tn[tgt]:.3f} ℃", fontsize=12.5, weight="bold")
    ax.axvline(0, color="k", lw=.8); ax.grid(alpha=.3, axis="x")
    ax.text(0.98, 0.03, f"近傍の温度差 {Tc[nb].max()-Tc[nb].min():.2f} ℃",
            transform=ax.transAxes, ha="right", fontsize=9.5, color="dimgray")

    # C: 写像後
    ax = axes[2]
    s2 = ax.scatter(nodes[sn, 0]*1000, nodes[sn, 1]*1000, c=Tn[sn], cmap="turbo",
                    s=46, vmin=vmin, vmax=vmax, edgecolors="k", lw=.3)
    ax.set_title("C: FEM節点に温度が乗った\n（これをFrontISTRへ渡す）", fontsize=12.5, weight="bold")
    ax.set_aspect("equal"); ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    fig.colorbar(s2, ax=ax, label="温度 [℃]")

    fig.suptitle(f"OpenFOAM → FrontISTR の温度マッピング（t=300 s、z={ZSLICE*1000:.0f} mm断面）"
                 f"：CFDセル{len(Cc):,}個 → FEM節点{len(nodes):,}個",
                 fontsize=13.5, weight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    out = os.path.join(IMG, "blog001_mapping.png")
    fig.savefig(out, dpi=140, bbox_inches="tight"); plt.close(fig)
    print("[map] wrote", out)
    print(f"[map] 節点温度の範囲: {Tn.min():.2f} 〜 {Tn.max():.2f} ℃"
          f"  (CFDセル: {Tc.min():.2f} 〜 {Tc.max():.2f} ℃)")


if __name__ == "__main__":
    main()
