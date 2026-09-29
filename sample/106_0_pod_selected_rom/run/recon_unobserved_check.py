"""5点から復元した温度分布が「未観測セル」でも合っているかを検証する.

確認する3点:
  1. 未観測20,691セルだけでの復元誤差（観測点を除外して評価）
  2. 誤差の時間変化と空間分布（どこが苦手か）
  3. 訓練外での検証（前半でPOD/選点を作り、後半で評価＝外挿テスト）
出力: docs/img/recon_unobserved.png, results/recon_unobserved.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/recon_unobserved_check.py
"""
from __future__ import annotations
import os, sys, json, importlib.util
import numpy as np
from scipy.linalg import qr
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
spec = importlib.util.spec_from_file_location('q', os.path.join(HERE, 'select_points_qdeim.py'))
m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
IMG = os.path.join(ROOT, "docs", "img"); RES = os.path.join(ROOT, "results")
KC = 273.15; R = 5


def recon(U, mean, pts, X):
    """5点の値だけから全場を復元する（gappy-POD）"""
    a = np.linalg.pinv(U[pts, :]) @ (X[pts, :] - mean[pts, None])
    return mean[:, None] + U @ a


def main():
    ts, Cc, X = m.load_snapshots()
    mean = X.mean(axis=1); Xc = X - mean[:, None]
    U, S, _ = np.linalg.svd(Xc, full_matrices=False)
    Ur = U[:, :R]
    pts = list(qr(Ur.T, pivoting=True)[2][:R])
    unobs = np.setdiff1d(np.arange(len(X)), pts)
    print(f"観測 {len(pts)} セル / 未観測 {len(unobs):,} セル")

    # ── ① 訓練内（全121枚でPOD）──
    Xh = recon(Ur, mean, pts, X)
    err = Xh - X
    rmse_all = float(np.sqrt((err**2).mean()))
    rmse_un = float(np.sqrt((err[unobs]**2).mean()))
    rmse_ob = float(np.sqrt((err[pts]**2).mean()))
    rmse_t = np.sqrt((err[unobs]**2).mean(axis=0))
    max_t = np.abs(err[unobs]).max(axis=0)
    cell_rmse = np.sqrt((err[unobs]**2).mean(axis=1))
    print(f"\n【訓練内】observed {rmse_ob:.2e} K / 未観測 {rmse_un:.5f} K / 全体 {rmse_all:.5f} K")
    print(f"  未観測セルの最大誤差 {np.abs(err[unobs]).max():.4f} K"
          f"（温度変動幅 {X.max()-X.min():.2f} K に対し {np.abs(err[unobs]).max()/(X.max()-X.min())*100:.3f} %）")

    # ── ② 訓練外（前半でPOD・選点 → 後半で評価）──
    nh = X.shape[1] // 2
    Xtr, Xte = X[:, :nh], X[:, nh:]
    mtr = Xtr.mean(axis=1)
    Utr, _, _ = np.linalg.svd(Xtr - mtr[:, None], full_matrices=False)
    Utr = Utr[:, :R]
    ptr = list(qr(Utr.T, pivoting=True)[2][:R])
    un_tr = np.setdiff1d(np.arange(len(X)), ptr)
    Xte_h = recon(Utr, mtr, ptr, Xte)
    e_te = Xte_h - Xte
    rmse_te = float(np.sqrt((e_te[un_tr]**2).mean()))
    print(f"\n【訓練外】前半{nh}枚でPOD・選点 → 後半{X.shape[1]-nh}枚で評価")
    print(f"  未観測セルRMSE {rmse_te:.5f} K / 最大 {np.abs(e_te[un_tr]).max():.4f} K")

    # ── 図 ──
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(16.2, 5.6))
    a0.plot(ts, rmse_t, "-o", color="#c0392b", ms=4, lw=2.2, label="未観測セルのRMSE")
    a0.plot(ts, max_t, "--", color="#e8a99f", lw=2.0, label="未観測セルの最大誤差")
    a0.axvspan(0, 300, color="orange", alpha=.08)
    a0.set_yscale("log"); a0.set_xlabel("時刻 [s]"); a0.set_ylabel("復元誤差 [K]")
    a0.set_title(f"① 未観測{len(unobs):,}セルの誤差（時間変化）\n"
                 f"全時刻RMSE = {rmse_un:.5f} K", fontsize=12.5, weight="bold")
    a0.legend(fontsize=10); a0.grid(alpha=.3, which="both")

    a1.hist(cell_rmse * 1000, bins=60, color="#3b74b8", alpha=.8)
    a1.axvline(rmse_un * 1000, color="#c0392b", lw=2.5,
               label=f"全未観測セルのRMSE {rmse_un*1000:.2f} mK")
    a1.set_yscale("log"); a1.set_xlabel("セルごとの復元RMSE [mK]"); a1.set_ylabel("セル数")
    a1.set_title("② セルごとの誤差の分布\n"
                 f"最悪のセルでも {cell_rmse.max()*1000:.1f} mK", fontsize=12.5, weight="bold")
    a1.legend(fontsize=10); a1.grid(alpha=.3, which="both")

    lbl = ["観測点\n（5セル）", "未観測20,691\n（訓練内）", "未観測20,691\n（訓練外）"]
    val = [rmse_ob, rmse_un, rmse_te]
    FLOOR = 3e-5                       # 観測点は厳密に一致（0 K）なので対数軸では床値で描く
    b = a2.bar(range(3), [max(v, FLOOR) for v in val],
               color=["#8899aa", "#3b74b8", "#c0392b"])
    txt = ["0 K（定義上ぴったり一致）", f"{rmse_un*1000:.2f} mK", f"{rmse_te*1000:.2f} mK"]
    for r, v, tx in zip(b, val, txt):
        a2.text(r.get_x()+r.get_width()/2, max(v, FLOOR), tx,
                ha="center", va="bottom", fontsize=10.5, weight="bold")
    a2.set_ylim(FLOOR*0.7, 1.0)
    a2.axhline(0.30, color="#2e9e5b", ls="--", lw=2)
    a2.text(2.45, 0.32, "観測ノイズ 0.3 K", fontsize=10, color="#2e9e5b",
            ha="right", weight="bold")
    a2.set_yscale("log"); a2.set_xticks(range(3)); a2.set_xticklabels(lbl, fontsize=10.5)
    a2.set_ylabel("復元RMSE [K]")
    a2.set_title("③ 観測点 vs 未観測点 vs 訓練外\n未観測でも観測ノイズより3桁小さい",
                 fontsize=12.5, weight="bold")
    a2.grid(alpha=.3, axis="y", which="both")

    fig.suptitle("5点から復元した温度分布は、測っていないセルでも合っているか ― 検証結果",
                 fontsize=14.5, weight="bold", y=0.985)
    fig.tight_layout(rect=[0, 0, 1, 0.935])
    out = os.path.join(IMG, "recon_unobserved.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    with open(os.path.join(RES, "recon_unobserved.json"), "w") as f:
        json.dump(dict(n_obs=len(pts), n_unobs=int(len(unobs)), rmse_obs=rmse_ob,
                       rmse_unobs=rmse_un, rmse_unobs_max=float(np.abs(err[unobs]).max()),
                       rmse_unobs_holdout=rmse_te,
                       worst_cell_rmse=float(cell_rmse.max())), f, indent=2)
    print("\nwrote", out)


if __name__ == "__main__":
    main()
