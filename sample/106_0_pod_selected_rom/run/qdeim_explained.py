"""Q-DEIMの貪欲選択を1ステップずつ再現し、枢軸QRと一致することを示す（blog_004 §4の検証）.

確認する内容:
  1. 「残差ノルム最大の点を選ぶ→その方向を差し引く」を手で回すと
     scipy.linalg.qr(pivoting=True) と同じ5点が出る
  2. 各ステップで残差ノルムがどれだけ減るか
  3. Q-DEIM点 vs ランダム5点 の復元精度・選点行列式の比較

再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/qdeim_explained.py
（102_0のOpenFOAM結果が必要）
"""
from __future__ import annotations
import os, sys, importlib.util
import numpy as np
from scipy.linalg import qr
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
R = 5


def main():
    np.set_printoptions(suppress=True, precision=4, floatmode='fixed')
    ts, C, X = m.load_snapshots()
    Xc = X - X.mean(axis=1)[:, None]
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    Ur = U[:, :R].copy()

    print("=" * 70)
    print("【貪欲選択を1ステップずつ】各セルは r=%d 次元のモード指紋（Uの行）を持つ" % R)
    W = Ur.copy(); picked = []
    for k in range(R):
        norms = np.linalg.norm(W, axis=1)
        p = int(np.argmax(norms)); picked.append(p)
        print(f"\n-- {k+1}点目: cell{p}  座標 {np.round(C[p]*1000,1)} mm  残差ノルム {norms[p]:.5f}")
        print(f"   指紋 = {Ur[p]}")
        v = W[p] / np.linalg.norm(W[p])
        W -= np.outer(W @ v, v)
        print(f"   この方向を全セルから差し引く → 残りの最大 {np.linalg.norm(W,axis=1).max():.5f}")

    _, _, piv = qr(Ur.T, pivoting=True)
    print("\n" + "=" * 70)
    print(f"手で追った貪欲選択: {picked}")
    print(f"scipy 枢軸QR      : {list(piv[:R])}")
    print(f"一致: {picked == list(piv[:R])}")

    def err(pts):
        a = np.linalg.pinv(Ur[pts, :]) @ Xc[pts, :]
        return float(np.sqrt(((Ur @ a - Xc) ** 2).mean()))

    rng = np.random.default_rng(0)
    rs = np.array([err(list(rng.choice(len(C), R, replace=False))) for _ in range(200)])
    print("\n" + "=" * 70)
    print("【5点から全場を復元したRMSE】")
    print(f"  Q-DEIM        : {err(picked):.5f} K")
    print(f"  ランダム200回 : 中央値 {np.median(rs):.5f} K / 最良 {rs.min():.5f} K / 最悪 {rs.max():.3f} K")
    d0 = abs(np.linalg.det(Ur[picked, :]))
    dr = [abs(np.linalg.det(Ur[list(rng.choice(len(C), R, replace=False)), :])) for _ in range(200)]
    print(f"  選点行列式 |det| : Q-DEIM {d0:.3e} / ランダム中央値 {np.median(dr):.3e}")


if __name__ == '__main__':
    main()
