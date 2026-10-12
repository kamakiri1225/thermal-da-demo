"""4,000 本の総当たりではなく、理論から最適なセンサ位置を求める（blog_009）.

使うのは、式で計算した温度のばらつき B(t)（analytic_prior_covariance.py）と、FrontISTR の重みだけ。

① センサ1本：閉じた式で順位を付ける
     減る割合 = cos²θ_B(h, w) × SNR/(1+SNR),   cos θ_B = (wᵀBh)/√(wᵀBw · hᵀBh),  SNR = hᵀBh / r
   → 「反りと同じ向きか」×「ノイズに埋もれないか」。総当たりの σ と同じ順位になるはず（確かめる）。
② 2方向に落とす（ランク2近似）：B ≈ λ1 e1e1ᵀ + λ2 e2e2ᵀ（全体の温まり＋温度差）。
   候補は平面上の点 (x,y) になり、2本の組み合わせの σ は 2×2 の式で閉じる。
   → 全ペアを厳密に（貪欲法なしで）評価して最適な2本を求め、貪欲法・総当たりと比べる。
③ 目的を複数にする：反り・U_z(A)・平均温度の「残り割合²」の合計 J を2方向で最小にする2本。
④ 逐次（リッカチ）σ：30 秒ごとに補正を積み重ねた後の σ を、式で計算（測定値不要）。
   代表的な配置について、実際の EnKF 誤差との順位を比べる。

出力: results/theory_based_selection.json, docs/img/theory_pair_plane.png
再現: OMP_NUM_THREADS=4 python3 run/theory_based_selection.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p  # noqa: F401
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from analytic_prior_covariance import setup, analytic_B, TIMES, DT
from select_nonrep_temp_disp import EXCL, A_XYZ, O_XYZ

RES = os.path.join(ROOT, "results"); IMG = os.path.join(ROOT, "docs", "img")
R = 0.09; NPT = 5


def post_sigma(Bs, rows, targets):
    """rows: list of (h, r)。targets: dict name->v。10時刻平均の σ を各 target について返す"""
    out = {k: [] for k in targets}
    for B in Bs:
        P = B.copy()
        for h, r in rows:
            Ph = P @ h; P = P - np.outer(Ph, Ph) / (h @ Ph + r)
        for k, v in targets.items(): out[k].append(np.sqrt(max(v @ P @ v, 0)))
    return {k: float(np.mean(v)) for k, v in out.items()}


def main():
    C, Km, heat = setup()
    Bs = analytic_B(C, Km, heat); Bm = np.mean(Bs, axis=0)
    kv = np.load(os.path.join(RES, "qdeim_points.npz")); U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); G = U @ UPp
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]; Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin()); nc = lambda x: int(np.linalg.norm(Cc - np.array(x), axis=1).argmin())
    iA, iO = nn(A_XYZ), nn(O_XYZ); w = Wall[iA, 2] - Wall[iO, 2]; wA = Wall[iA, 2]; mean5 = np.ones(5) / 5
    far = (np.linalg.norm(co - co[iA], axis=1) > EXCL) & (np.linalg.norm(co - co[iO], axis=1) > EXCL)
    Wf = Wall.reshape(-1, NPT); uidx = np.where(np.repeat(far, 3))[0]; Hu = Wf[uidx]
    tmask = np.ones(len(Cc), bool); tmask[pod] = False; tidx = np.where(tmask)[0]; Ht = G[tidx]
    comp = "xyz"
    lab_u = lambda k: f"変位 {comp[int(uidx[k]) % 3].upper()} {np.round(co[int(uidx[k]) // 3] * 1000, 1).tolist()}"
    lab_t = lambda k: f"温度 {np.round(Cc[int(tidx[k])] * 1000, 1).tolist()}"
    out = {"note": "B は式で計算（共分散の伝播、30〜300 秒の10時刻）。ノイズ 0.3。候補：変位 14,754、温度 20,691（代表点以外）。"}

    # ---------- ① 1本：閉じた式 vs 総当たり σ ----------
    t0 = time.time()
    def closed_form(H):   # cos²θ_B × SNR/(1+SNR) から σ を組み立てる（総当たりと同じ量になるはず）
        s = np.zeros(len(H)); frac = np.zeros(len(H))
        for B in Bs:
            wBw = w @ B @ w; Bh = B @ H.T; wBh = w @ Bh; hBh = np.maximum(np.einsum("ij,ji->i", H, Bh), 1e-30)
            cos2 = wBh ** 2 / (wBw * hBh); snr = hBh / R
            red = cos2 * snr / (1 + snr)                     # 減る割合
            frac += red; s += np.sqrt(np.maximum(wBw * (1 - red), 0))
        return s / len(Bs), frac / len(Bs)
    def brute(H):
        s = np.zeros(len(H))
        for B in Bs:
            Bh = B @ H.T; s += np.sqrt(np.maximum(w @ B @ w - (w @ Bh) ** 2 / (np.einsum("ij,ji->i", H, Bh) + R), 0))
        return s / len(Bs)
    (cf_u, fr_u), br_u = closed_form(Hu), brute(Hu); (cf_t, fr_t), br_t = closed_form(Ht), brute(Ht)
    rk = lambda a: np.argsort(np.argsort(a))
    out["one_sensor"] = {"time_s": time.time() - t0,
        "disp_max_abs_diff_closed_vs_brute_um": float(np.abs(cf_u - br_u).max()),
        "temp_max_abs_diff_closed_vs_brute_um": float(np.abs(cf_t - br_t).max()),
        "disp_rank_corr_fraction_vs_brute": float(np.corrcoef(rk(-fr_u), rk(br_u))[0, 1]),
        "disp_best": lab_u(int(br_u.argmin())), "temp_best": lab_t(int(br_t.argmin())),
        "note": "σ = √(wᵀBw·(1−cos²θ_B·SNR/(1+SNR))) を時刻平均したものは総当たりの σ と一致（差は丸め）。「減る割合」を時刻平均して並べた順位は σ の順位と少し違う。"}
    print("① 1本", json.dumps(out["one_sensor"], ensure_ascii=False))

    FIG_ONLY = bool(os.environ.get("FIG_ONLY"))
    if FIG_ONLY:
        prev = json.load(open(os.path.join(RES, "theory_based_selection.json"), encoding="utf-8"))
    # ---------- ② 2本の全ペアを厳密に ----------
    lam, E = np.linalg.eigh(Bm); o = np.argsort(lam)[::-1]; lam, E = lam[o], E[:, o]
    e1, e2 = E[:, 0] * np.sign(E[:, 0].sum()), E[:, 1] * np.sign(E[2, 1])
    share = lam / lam.sum()
    def plane(H):   # (n,2)：信号/ノイズ（図のため。2方向近似）
        return np.stack([np.sqrt(lam[0]) * (H @ e1), np.sqrt(lam[1]) * (H @ e2)], 1) / np.sqrt(R)
    def target2(v): return np.array([np.sqrt(lam[0]) * (v @ e1), np.sqrt(lam[1]) * (v @ e2)])
    Zu, Zt = plane(Hu), plane(Ht)
    tw = target2(w); tA = target2(wA); tm = target2(mean5)
    # 厳密版：各時刻の B を B = L Lᵀ と分解し、座標 z = Lᵀh/√r（5次元）で同じ 2×2 の式を使う
    Ls = [np.linalg.cholesky(B + 1e-12 * np.eye(5)) for B in Bs]
    Zu5 = [(Hu @ L) / np.sqrt(R) for L in Ls]; Zt5 = [(Ht @ L) / np.sqrt(R) for L in Ls]
    tw5 = [L.T @ w for L in Ls]; tA5 = [L.T @ wA for L in Ls]; tm5 = [L.T @ mean5 for L in Ls]
    def pair_sigma2(Z1, Z2, tv):
        """2方向近似で、センサ1（Z1: (n,2) 固定1本ずつ）とセンサ2（Z2: (m,2)）の全組み合わせの σ²（目的 tv）を返す (n,m)"""
        # 事後：P' = I - Zᵀ(Z Zᵀ + I)^{-1} Z（単位化した座標）。 Z=(z1;z2) 2×2 なので閉じて書ける
        a11 = (Z1 ** 2).sum(1)[:, None] + 1; a22 = (Z2 ** 2).sum(1)[None, :] + 1; a12 = Z1 @ Z2.T
        det = a11 * a22 - a12 ** 2
        # P' tv = tv - Zᵀ S^{-1} Z tv ；S^{-1} = [[a22,-a12],[-a12,a11]]/det
        b1 = (Z1 @ tv)[:, None]; b2 = (Z2 @ tv)[None, :]
        q = (a22 * b1 ** 2 - 2 * a12 * b1 * b2 + a11 * b2 ** 2) / det
        return tv @ tv - q
    def best_pair(ZAs, ZBs, tvs, same, weights=None):
        """全ペアを厳密に評価。ZAs/ZBs/tvs は時刻ごとのリスト。weights があれば複数目的の J（残り割合²の和）"""
        best = (np.inf, -1, -1); t0 = time.time(); n = len(ZAs[0])
        for i in range(0, n, 300):
            acc = np.zeros((min(300, n - i), len(ZBs[0])))
            for k in range(len(ZAs)):
                if weights is None:
                    acc += np.sqrt(np.maximum(pair_sigma2(ZAs[k][i:i + 300], ZBs[k], tvs[k]), 0)) / len(ZAs)
                else:
                    for tv, s0 in zip(tvs, weights): acc += pair_sigma2(ZAs[k][i:i + 300], ZBs[k], tv[k]) / s0 ** 2 / len(ZAs)
            if same:
                for r_ in range(acc.shape[0]): acc[r_, :i + r_ + 1] = np.inf
            kk = int(np.argmin(acc)); r_, c_ = divmod(kk, acc.shape[1])
            if acc[r_, c_] < best[0]: best = (float(acc[r_, c_]), i + r_, int(c_))
        return best, time.time() - t0
    sig0_2 = float(np.sqrt(tw @ tw))
    res2 = {}
    for name, (ZA, la, ZB, lb, same) in {"変位2本": (Zu5, lab_u, Zu5, lab_u, True), "温度2本": (Zt5, lab_t, Zt5, lab_t, True),
                                       "温度1本＋変位1本": (Zt5, lab_t, Zu5, lab_u, False)}.items():
        if FIG_ONLY:
            res2[name] = prev["pair_exact"][name]; continue
        (sg, i, j), dt = best_pair(ZA, ZB, tw5, same)
        res2[name] = {"best": [la(i), lb(j)], "sigma_exact_pairs_um": sg, "time_s": dt}
        # 同じ2本を、フルの B で評価
        rows = [(Hu[i] if la is lab_u else Ht[i], R), (Hu[j] if lb is lab_u else Ht[j], R)]
        res2[name]["sigma_full_um"] = post_sigma(Bs, rows, {"AO": w, "UzA": wA, "meanT": mean5})
    out["pair_exact"] = {"eigen_share": share.round(4).tolist(), **res2,
                         "note": "全ペアを厳密に評価（5次元、10時刻の σ 平均）。貪欲法（blog_009）：変位2本 0.296、温度2本 0.365 µm"}
    print("② ペア（全ペア厳密）", json.dumps(res2, ensure_ascii=False, indent=1))

    # ---------- ③ 目的を複数にして、2方向で2本を選ぶ ----------
    none = post_sigma(Bs, [], {"AO": w, "UzA": wA, "meanT": mean5})
    wts = [none["AO"], none["UzA"], none["meanT"]]
    res3 = {}
    for name, (ZA, la, HA, ZB, lb, HB, same) in {"変位2本": (Zu5, lab_u, Hu, Zu5, lab_u, Hu, True), "温度2本": (Zt5, lab_t, Ht, Zt5, lab_t, Ht, True),
                                               "温度1本＋変位1本": (Zt5, lab_t, Ht, Zu5, lab_u, Hu, False)}.items():
        if FIG_ONLY:
            res3[name] = prev["multi_objective_exact"][name]; continue
        (J, i, j), dt = best_pair(ZA, ZB, [tw5, tA5, tm5], same, weights=wts)
        res3[name] = {"best": [la(i), lb(j)], "J_exact": J, "time_s": dt,
                      "sigma_full_um": post_sigma(Bs, [(HA[i], R), (HB[j], R)], {"AO": w, "UzA": wA, "meanT": mean5})}
    res3["参考：blog_009 の変位2本（反りだけで選んだ）"] = {"sigma_full_um": post_sigma(Bs, [(Wall[nn((0.0188, 0.0325, 0.1005)), 0], R), (Wall[nn((0.0188, -0.0325, 0.1005)), 0], R)], {"AO": w, "UzA": wA, "meanT": mean5})}
    res3["参考：何も測らない"] = {"sigma_full_um": post_sigma(Bs, [], {"AO": w, "UzA": wA, "meanT": mean5})}
    out["multi_objective_exact"] = res3
    print("③ 複数目的", json.dumps(res3, ensure_ascii=False, indent=1))

    # ---------- ④ 逐次（リッカチ）σ ----------
    n = 5; Ci = 1.0 / C; Tair = rg.T_AIR_K; Qr = rg.HEATER_RATED_W; mh = 0.02; mq = 1.05
    Ksum = Km.sum(1); A5 = (Km - np.diag(Ksum) - mh * np.eye(n)) * Ci[:, None]; g = np.zeros(n); g[heat] = Qr * Ci[heat]
    def riccati(rows):
        P = np.diag([15.0**2 / 12] * n + [1.5**2 / 12, 0.01**2]); T = np.full(n, Tair + 4.5); t = 0.0; sig = []
        H7 = np.array([np.r_[h, 0, 0] for h, _ in rows]) if rows else np.zeros((0, 7))
        for tb in TIMES:
            for _ in range(int(round((tb - t) / DT))):
                F = np.zeros((7, 7)); F[:5, :5] = A5; F[:5, 5] = g; F[:5, 6] = -Ci * (T - Tair)
                Phi = np.eye(7) + DT * F + 0.5 * DT**2 * F @ F
                P = Phi @ P @ Phi.T; T = T + DT * (A5 @ T + mh * Tair * Ci + g * mq)
            t = tb
            if len(rows):
                S = H7 @ P @ H7.T + R * np.eye(len(rows)); P = P - P @ H7.T @ np.linalg.solve(S, H7 @ P)
            sig.append(np.sqrt(max(w @ P[:5, :5] @ w, 0)))
        return float(np.mean(sig))
    S1 = Wall[nn((0.0097, -0.0362, 0.1005)), 0]; S2 = Wall[nn((-0.0346, 0.0144, 0.1005)), 2]
    U2a, U2b = Wall[nn((0.0188, 0.0325, 0.1005)), 0], Wall[nn((0.0188, -0.0325, 0.1005)), 0]
    T2a, T2b = G[nc((0.0367, 0.0008, 0.0492))], G[nc((-0.0329, 0.0016, 0.0516))]
    e = np.eye(5)
    cfgs = {"温度 P2・P0": [(e[2], R), (e[0], R)], "温度 P2・P4": [(e[2], R), (e[4], R)], "変位 S1・S2": [(S1, R), (S2, R)],
            "温度2本（blog_009）": [(T2a, R), (T2b, R)], "変位2本（blog_009）": [(U2a, R), (U2b, R)],
            "温度 P2・P0＋変位 S1・S2": [(e[2], R), (e[0], R), (S1, R), (S2, R)]}
    actual = {"温度 P2・P0": 0.769, "温度 P2・P4": 0.231, "変位 S1・S2": 0.164, "温度2本（blog_009）": 0.197, "変位2本（blog_009）": 0.120,
              "温度 P2・P0＋変位 S1・S2": 0.075}
    res4 = {}
    for k, rows in cfgs.items():
        res4[k] = {"sigma_single_um": post_sigma(Bs, rows, {"AO": w})["AO"], "sigma_sequential_um": riccati(rows), "actual_enkf_um": actual[k]}
    names = list(cfgs)
    rc = lambda a, b: float(np.corrcoef(rk(np.array(a)), rk(np.array(b)))[0, 1])
    out["sequential"] = {"configs": res4,
        "rank_corr_single_vs_actual": rc([res4[k]["sigma_single_um"] for k in names], [actual[k] for k in names]),
        "rank_corr_sequential_vs_actual": rc([res4[k]["sigma_sequential_um"] for k in names], [actual[k] for k in names]),
        "note": "actual は 15 W・実ソルバ真値・5 seed の反り A−O 誤差（blog_008 §6、blog_009 §4）。リッカチはモデル誤差 Q=0、平均軌道まわりの一次近似。"}
    print("④ 逐次", json.dumps(out["sequential"], ensure_ascii=False, indent=1))
    if not FIG_ONLY:
        json.dump(out, open(os.path.join(RES, "theory_based_selection.json"), "w"), ensure_ascii=False, indent=1, default=float)

    # ---------- 図：平面上に、②③で選ばれたペアを描く ----------
    fig, ax = plt.subplots(figsize=(13, 10))
    ax.scatter(Zu[:, 0], Zu[:, 1], s=5, color="#C7D2E0", alpha=.5, label="変位計の候補")
    ax.scatter(Zt[:, 0], Zt[:, 1], s=5, color="#F2C4B0", alpha=.5, label="温度計の候補")
    ax.annotate("反り A−O の向き", tw, xytext=(10, 10), textcoords="offset points", fontsize=12, color="#14213D", weight="bold")
    ax.arrow(0, 0, tw[0], tw[1], color="#14213D", width=0.05, length_includes_head=True)
    ax.arrow(0, 0, tA[0], tA[1], color="#8E44AD", width=0.05, length_includes_head=True)
    ax.annotate("U_z(A) の向き", tA, xytext=(10, -14), textcoords="offset points", fontsize=12, color="#8E44AD", weight="bold")
    for name, col, src in [("変位2本", "#C0392B", res2), ("温度1本＋変位1本", "#2E8B57", res2), ("温度2本", "#E07A1F", res2)]:
        pts = []
        for lab in src[name]["best"]:
            if lab.startswith("変位"): k = [i for i in range(len(uidx)) if lab_u(i) == lab][0]; pts.append(Zu[k])
            else: k = [i for i in range(len(tidx)) if lab_t(i) == lab][0]; pts.append(Zt[k])
        pts = np.array(pts); ax.scatter(pts[:, 0], pts[:, 1], s=180, marker="D", color=col, edgecolor="k", zorder=5, label=f"②目的＝反りだけ：{name}")
    for name, col in [("変位2本", "#C0392B"), ("温度1本＋変位1本", "#2E8B57")]:
        pts = []
        for lab in res3[name]["best"]:
            if lab.startswith("変位"): k = [i for i in range(len(uidx)) if lab_u(i) == lab][0]; pts.append(Zu[k])
            else: k = [i for i in range(len(tidx)) if lab_t(i) == lab][0]; pts.append(Zt[k])
        pts = np.array(pts); ax.scatter(pts[:, 0], pts[:, 1], s=220, marker="*", color=col, edgecolor="k", zorder=6, label=f"③目的＝反り＋Uz(A)＋平均温度：{name}")
    ax.axhline(0, color="#999", lw=1); ax.axvline(0, color="#999", lw=1); ax.grid(alpha=.3)
    ax.set_xlabel(f"全体の温まり（{share[0]*100:.0f} %）の信号÷ノイズ"); ax.set_ylabel(f"温度差（{share[1]*100:.0f} %）の信号÷ノイズ")
    ax.set_title("2本の全ペアを厳密に選んだ結果を、2方向の平面に描く", fontsize=16); ax.legend(fontsize=11, loc="upper left")
    fig.tight_layout(); fig.savefig(os.path.join(IMG, "theory_pair_plane.png"), dpi=100, facecolor="white"); plt.close(fig)
    print("wrote theory_pair_plane.png")


if __name__ == "__main__":
    main()
