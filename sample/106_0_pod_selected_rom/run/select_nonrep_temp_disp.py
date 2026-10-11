"""代表点以外の温度計2本＋変位計2本を、反り A−O の σ で選ぶ（blog_008 §5-2 の続き）.

blog_008 §5-2 では温度計を代表点 P2・P0 に固定し、変位計2本だけを選んだ。
ここでは温度計も「代表点以外の全セル」から選ぶ。
  温度計の候補：代表5点を除く全セル（20,691）。観測行列の行 g_j = ψ_j U_rP^{-1}（gappy-POD）
  変位計の候補：A・O から 12 mm 以内を除く全節点×3方向（14,754）。行 h_c（FrontISTR 応答）
  採点：σ_S(t_k) = sqrt(w^T P_S(t_k) w) を 30〜300 秒の10時刻で平均（§5-2 と同じ 4000 本の事前分布）
  選び方：何も測らない状態から1本ずつ、点数が最小の候補を足す（貪欲法）。温度計2本・変位計2本に達した種類は候補から外す。
センサを1本足したときの共分散の更新は、1行のカルマン更新（ランク1更新）で計算する。データ同化も測定値も使わない。

出力: results/nonrep_selection.json
再現: OMP_NUM_THREADS=4 python3 run/select_nonrep_temp_disp.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import rom_general as rg

RES = os.path.join(ROOT, "results")
NPT = 5; R_T = 0.30**2; R_U = 0.30**2; EXCL = 0.012
A_XYZ = (0.028, 0.0, 0.1005); O_XYZ = (-0.028, 0.0, 0.1005)
COMP = ["Ux", "Uy", "Uz"]


def priors():
    d = np.load(os.path.join(RES, "rom_calibrated_pod.npz"))
    C = d["C"]; Km = rg.tri_to_matrix(d["K_upper"], NPT); heat = int(d["heat_node"])
    rng = np.random.default_rng(20260930); N = 4000
    Tp = rng.uniform(rg.T_AIR_K - 3, rg.T_AIR_K + 12, (N, NPT))
    qp = rng.uniform(0.3, 1.8, N); hp = np.clip(rng.normal(0.02, 0.01, N), 1e-3, 0.1)
    Bs = []; tp = 0.0
    for tb in np.arange(30.0, 300.0 + 1e-9, 30.0):
        Tp = rg.integrate_ensemble(Tp, C, Km, hp, qp, heat, tp, tb, 2.0); tp = tb; Bs.append(np.cov(Tp.T))
    return Bs


def add_row(P, h, r):
    Ph = P @ h; return P - np.outer(Ph, Ph) / (h @ Ph + r)


def score(Ps, w, Hc, r):
    """候補行 Hc (n,5) それぞれを1本足したときの時刻平均 σ"""
    s = np.zeros(len(Hc))
    for P in Ps:
        HP = Hc @ P
        s += np.sqrt(np.maximum(w @ P @ w - (HP @ w) ** 2 / (np.einsum("ij,ij->i", HP, Hc) + r), 0))
    return s / len(Ps)


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); G = U @ UPp
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    w = Wall[nn(A_XYZ), 2] - Wall[nn(O_XYZ), 2]
    tmask = np.ones(len(Cc), bool); tmask[pod] = False
    tidx = np.where(tmask)[0]
    far = (np.linalg.norm(co - np.array(A_XYZ), axis=1) > EXCL) & (np.linalg.norm(co - np.array(O_XYZ), axis=1) > EXCL)
    Wf = Wall.reshape(-1, NPT); uidx = np.where(np.repeat(far, 3))[0]
    Bs = priors(); Ps = [B.copy() for B in Bs]
    sig0 = float(np.mean([np.sqrt(w @ B @ w) for B in Bs]))
    print(f"何も測らない σ = {sig0:.3f} µm   候補: 温度 {len(tidx)} セル、変位 {len(uidx)}", flush=True)
    picked = []; nT = nU = 0; hist = [{"step": 0, "sigma_um": sig0}]
    for step in range(4):
        best = None
        if nT < 2:
            sT = score(Ps, w, G[tidx], R_T); k = int(sT.argmin()); best = ("T", int(tidx[k]), float(sT[k]))
        if nU < 2:
            sU = score(Ps, w, Wf[uidx], R_U); k = int(sU.argmin())
            if best is None or sU[k] < best[2]: best = ("U", int(uidx[k]), float(sU[k]))
        kind, idx, s = best
        if kind == "T":
            h, r = G[idx], R_T; nT += 1; tidx = tidx[tidx != idx]
            info = {"kind": "温度", "cell": idx, "xyz_mm": np.round(Cc[idx] * 1000, 1).tolist()}
        else:
            h, r = Wf[idx], R_U; nU += 1; uidx = uidx[uidx != idx]; i, c = divmod(idx, 3)
            info = {"kind": "変位", "node": i, "comp": COMP[c], "xyz_mm": np.round(co[i] * 1000, 1).tolist()}
        Ps = [add_row(P, h, r) for P in Ps]
        info.update(step=step + 1, sigma_um=s, weights=np.round(h, 4).tolist()); picked.append(info); hist.append(info)
        print(f"  {step+1}本目: {info['kind']} {info['xyz_mm']} {info.get('comp','')}  σ = {s:.3f} µm", flush=True)

    # 比較用：今の配置（温度 P2・P0 ＋ S1・S2）と、温度だけの各配置の σ
    def sig_of(rows):
        out = []
        for B in Bs:
            P = B.copy()
            for h, r in rows: P = add_row(P, h, r)
            out.append(np.sqrt(max(w @ P @ w, 0)))
        return float(np.mean(out))
    e = np.eye(NPT)
    S1 = Wall[nn((0.0097, -0.0362, 0.1005)), 0]; S2 = Wall[nn((-0.0346, 0.0144, 0.1005)), 2]
    selT = [(G[p["cell"]], R_T) for p in picked if p["kind"] == "温度"]
    comp = {"温度 P2・P0 のみ": sig_of([(e[2], R_T), (e[0], R_T)]),
            "温度 P2・P4 のみ": sig_of([(e[2], R_T), (e[4], R_T)]),
            "温度 P2・P0 ＋ 変位 S1・S2（今の配置）": sig_of([(e[2], R_T), (e[0], R_T), (S1, R_U), (S2, R_U)]),
            "選んだ温度2点のみ（代表点以外）": sig_of(selT),
            "選んだ温度2点＋変位2点（代表点以外）": hist[-1]["sigma_um"]}
    for k, v in comp.items(): print(f"  σ  {k:32s} {v:.3f} µm")
    json.dump({"note": "反り A−O の σ（30〜300 秒の10時刻平均、µm）。4000本の事前分布、ノイズ 温度0.3 K・変位0.3 µm。データ同化はしていない。",
               "sigma_none_um": sig0, "picked": picked, "sigma_compare_um": comp},
              open(os.path.join(RES, "nonrep_selection.json"), "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
