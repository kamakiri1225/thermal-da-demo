"""ROM の係数合わせ（build_calibrate_rom.py と同じ最小二乗）を、ばらばらの初期値から何度もやり直す（多スタート）。
局所最適に落ちていないか、違う係数の組で同じくらい合うか（一意に決まらないか）を調べる。

  初期値：各係数を範囲の中から対数一様で選ぶ（C：10〜2000 J/K、K：0.01〜50 W/K、h：0.001〜2 W/K）。40 通り。
  あわせて、保存してある最適解でのヤコビアン（残差の係数についての微分）の特異値を求め、
  「係数を変えても残差がほとんど変わらない方向」があるかを見る。
  出力：results/calib_multistart.json、docs/img/calib_multistart.png
  再現：python3 run/calib_multistart.py
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.optimize import least_squares
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import rom_general as rg
from select_points_qdeim import load_snapshots

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
NPT = 5; DT = 2.0; NE = rg.n_edges(NPT); NSTART = 40
LB = np.concatenate([np.full(NPT, 1.0), np.full(NE, 1e-4), [1e-4]])
UB = np.concatenate([np.full(NPT, 5000.0), np.full(NE, 200.0), [50.0]])
G = {}


def init():
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); cells = kv["cell_idx"]
    ts, _, X = load_snapshots()
    G.update(ts=ts, Y=X[cells, :].T, heat=int(np.load(os.path.join(RES, "rom_calibrated_pod.npz"))["heat_node"]))


def forward(v):
    C = v[:NPT]; K = rg.tri_to_matrix(v[NPT:NPT + NE], NPT); h = v[NPT + NE]
    T0 = np.full(NPT, rg.T_AIR_K); Y = [T0.copy()]; T = T0.reshape(1, -1).copy()
    for a, b in zip(G["ts"][:-1], G["ts"][1:]):
        T = rg.integrate_ensemble(T, C, K, np.array([h]), np.array([1.0]), G["heat"], a, b, DT); Y.append(T[0].copy())
    return np.array(Y)


def resid(v):
    with np.errstate(all="ignore"):
        r = (forward(v) - G["Y"]).ravel()
    return np.clip(np.nan_to_num(r, nan=1e3, posinf=1e3, neginf=-1e3), -1e3, 1e3)   # 刻み 2 秒で発散する係数は大きな残差にする


def one(seed):
    rng = np.random.default_rng(seed)
    x0 = np.concatenate([np.exp(rng.uniform(np.log(10), np.log(2000), NPT)),
                         np.exp(rng.uniform(np.log(0.01), np.log(50), NE)),
                         [np.exp(rng.uniform(np.log(0.001), np.log(2)))]])
    try:
        sol = least_squares(resid, x0, bounds=(LB, UB), xtol=1e-12, ftol=1e-12)
    except Exception as e:
        return dict(seed=seed, rmse=float("nan"), x=[float("nan")] * len(x0), nfev=0, error=str(e))
    return dict(seed=seed, rmse=float(np.sqrt(np.mean(sol.fun ** 2))), x=sol.x.tolist(), nfev=int(sol.nfev))


def main():
    init(); t0 = time.time()
    with Pool(min(NSTART, 40), initializer=init) as p:
        runs = p.map(one, range(NSTART))
    best = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    xb = np.concatenate([best["C"], best["K_upper"], [best["h"]]])
    # 最適解でのヤコビアン（相対的な変化：∂r/∂ln θ）の特異値
    eps = 1e-4; r0 = resid(xb); J = np.zeros((r0.size, xb.size))
    for i in range(xb.size):
        xp = xb.copy(); xp[i] *= (1 + eps); J[:, i] = (resid(xp) - r0) / eps
    sv = np.linalg.svd(J, compute_uv=False)
    out = dict(note="多スタート（40 通り）。rmse は 5 点 × 121 時刻の温度の RMSE [K]。sv は ∂r/∂lnθ の特異値",
               time_s=time.time() - t0, best_saved=dict(rmse=float(np.sqrt(np.mean(r0 ** 2))), x=xb.tolist()),
               runs=sorted(runs, key=lambda r: r["rmse"]), singular_values=sv.tolist())
    json.dump(out, open(os.path.join(RES, "calib_multistart.json"), "w"), ensure_ascii=False, indent=1)
    rm = np.array([r["rmse"] for r in runs]); X = np.array([r["x"] for r in runs])
    print("失敗", int(np.isnan(rm).sum()), "/", len(rm))
    print(f"時間 {out['time_s']:.0f} s  RMSE 最小 {np.nanmin(rm):.4f} 中央 {np.nanmedian(rm):.4f} 最大 {np.nanmax(rm):.4f}  保存値 {out['best_saved']['rmse']:.4f}")
    good = np.nan_to_num(rm, nan=9) < 0.02
    print("RMSE<0.02 K の数", good.sum(), "/", len(rm))
    names = [f"C{i}" for i in range(5)] + [f"K{i}{j}" for i in range(5) for j in range(i + 1, 5)] + ["h"]
    for n, col in zip(names, X[good].T):
        print(f"{n:4s} 最小 {col.min():9.4f} 最大 {col.max():9.4f}")
    print("ΣC", np.round(X[good, :5].sum(1), 0))
    print("特異値（大→小）", np.round(sv, 4), "比", sv[0] / sv[-1])


if __name__ == "__main__":
    main()
