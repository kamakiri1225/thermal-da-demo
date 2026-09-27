"""体積重み付きPODと等重みPODを比較する（blog_004 §2-3の限界検証）.

等重みSVDは各セルを同列に扱うため、細かいセルが密集した領域を実質的に重く見る。
厳密なPOD（連続体のエネルギー最大化）は体積重み sqrt(V_i) を掛けた行列を分解する。
本スクリプトは両者のモード・エネルギー・選点を比較し、影響の大きさを数値で示す。

前提: 102_0 で `postProcess -region solid -func writeCellVolumes -time 5` 実行済み
     （5/solid/V が存在すること）
出力: results/pod_volume_weight.json（比較結果）
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/pod_volume_weight_check.py
"""
from __future__ import annotations
import os, sys, json, importlib.util
import numpy as np
from scipy.linalg import qr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
RES = os.path.join(ROOT, 'results')
NPT = 5; RMODE = 5


def main():
    ts, C, X = m.load_snapshots()
    V = m.read_foam_field(os.path.join(m.OF, '5', 'solid', 'V'), len(C))
    print(f"[vol] セル体積 V: {V.shape}  min={V.min():.3e} max={V.max():.3e} m^3")
    print(f"[vol] 最大/最小の比 = {V.max()/V.min():.1f} 倍  （不均一の度合い）")
    print(f"[vol] 合計体積 = {V.sum()*1e6:.1f} cm^3")

    mean = X.mean(axis=1); Xc = X - mean[:, None]
    w = np.sqrt(V)                       # 体積重み

    U0, S0, _ = np.linalg.svd(Xc, full_matrices=False)          # 等重み（従来）
    Uw, Sw, _ = np.linalg.svd(w[:, None] * Xc, full_matrices=False)  # 体積重み
    Uw_phys = Uw / w[:, None]            # 物理空間のモードへ戻す

    e0 = S0**2 / (S0**2).sum(); ew = Sw**2 / (Sw**2).sum()
    print("\n[vol] 累積エネルギー[%]")
    print(f"   等重み  : {np.round(np.cumsum(e0[:4])*100, 4)}")
    print(f"   体積重み: {np.round(np.cumsum(ew[:4])*100, 4)}")

    # モードの向きの一致度（符号を除いた内積）
    print("\n[vol] モードの一致度（体積内積で正規化したcos）")
    cos = []
    for k in range(3):
        a = U0[:, k] * w; b = Uw[:, k]
        c = abs(a @ b) / (np.linalg.norm(a) * np.linalg.norm(b))
        cos.append(float(c))
        print(f"   モード{k+1}: cos = {c:.6f}  （1なら同じ向き）")

    # Q-DEIM選点が変わるか
    p0 = list(qr(U0[:, :RMODE].T, pivoting=True)[2][:NPT])
    pw = list(qr(Uw_phys[:, :RMODE].T, pivoting=True)[2][:NPT])
    print(f"\n[vol] Q-DEIM選点  等重み  : {p0}")
    print(f"[vol] Q-DEIM選点  体積重み: {pw}")
    print(f"[vol] 一致した点数: {len(set(p0) & set(pw))}/{NPT}")

    # 実害の評価: 5点から全場を復元したときの体積加重RMSE
    def recon_err(U, pts):
        UP = np.linalg.pinv(U[pts, :RMODE])
        a = UP @ Xc[pts, :]
        Xr = U[:, :RMODE] @ a
        d2 = (Xr - Xc)**2
        return float(np.sqrt((d2 * V[:, None]).sum() / (V.sum() * Xc.shape[1])))
    r00 = recon_err(U0, p0); rww = recon_err(Uw_phys, pw)
    print("\n[vol] 5点復元の体積加重RMSE [K]")
    print(f"   等重みPOD + 等重み選点  : {r00:.6f}")
    print(f"   体積重みPOD + 体積重み選点: {rww:.6f}")

    out = dict(nvol=int(len(V)), vmax_vmin=float(V.max()/V.min()),
               cum_uniform=[float(x) for x in np.cumsum(e0[:4])],
               cum_weighted=[float(x) for x in np.cumsum(ew[:4])],
               mode_cos=cos, pts_uniform=[int(x) for x in p0],
               pts_weighted=[int(x) for x in pw],
               rmse_uniform=r00, rmse_weighted=rww)
    with open(os.path.join(RES, 'pod_volume_weight.json'), 'w') as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    print("\n[vol] wrote results/pod_volume_weight.json")


if __name__ == '__main__':
    main()
