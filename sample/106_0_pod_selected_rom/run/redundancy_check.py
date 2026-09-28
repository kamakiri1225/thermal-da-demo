"""「感度が高い点を並べる」となぜ悪いのかを、情報の重なりで定量化する（blog_005 §4-2の検証）.

同じ2点でも、近接させると観測が冗長になり連立が解けなくなる。
相関係数・観測行列の行列式・条件数で示す。

再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/redundancy_check.py
"""
from __future__ import annotations
import os, importlib.util
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)


def main():
    ts, Cc, X = m.load_snapshots()
    Xc = X - X.mean(axis=1)[:, None]
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    order = np.argsort(-np.abs(U[:, 0]))          # 場への効きが大きい順

    top2 = [int(order[0]), int(order[1])]          # 感度最大の上位2点
    p0 = int(order[0])
    far2 = None
    for c in order[1:]:
        if np.linalg.norm(Cc[c] - Cc[p0]) > 0.03:  # 30mm以上離す
            far2 = [p0, int(c)]; break

    for name, pts in [("感度最大の上位2点", top2), ("散らした2点", far2)]:
        d = np.linalg.norm(Cc[pts[0]] - Cc[pts[1]]) * 1000
        r = np.corrcoef(Xc[pts[0]], Xc[pts[1]])[0, 1]
        H = U[pts, :2]
        print(f"\n{name}: セル{pts}")
        print(f"  距離 {d:.1f} mm / 温度変動の相関 {r:+.5f}")
        print(f"  観測行列 |det| {abs(np.linalg.det(H)):.3e} / 条件数 {np.linalg.cond(H):.1f}")


if __name__ == '__main__':
    main()
