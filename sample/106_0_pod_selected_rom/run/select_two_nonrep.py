"""代表点以外の場所から「2本だけ」選んで、反り A−O を当てる配置を σ で選ぶ（blog_009）.

候補
  温度計：代表5点を除く全セル（20,691）。行 g_j = ψ_j Φ_P^{-1}（gappy-POD）
  変位計：A・O から 12 mm 以内を除く全節点×3方向（14,754）。行 h_c（FrontISTR 応答）
選び方（3通り）
  T2：温度計を2本   U2：変位計を2本   ANY2：種類を問わず2本
  何も測らない状態から、σ（30〜300 秒の10時刻平均）が一番小さくなる候補を1本ずつ足す（貪欲法）。
比較用に、代表点の温度計 P2・P0（片側）と P2・P4（両側）の σ も計算する。データ同化も測定値も使わない。

出力: results/two_nonrep_selection.json
再現: OMP_NUM_THREADS=4 python3 run/select_two_nonrep.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from select_nonrep_temp_disp import priors, add_row, score, NPT, R_T, R_U, EXCL, A_XYZ, O_XYZ, COMP

RES = os.path.join(ROOT, "results")


def main():
    kv = np.load(os.path.join(RES, "qdeim_points.npz"))
    U = kv["pod_modes"].astype(float); pod = kv["cell_idx"]; Cc = kv["cell_centres"]
    UPp = np.linalg.pinv(U[pod, :]); G = U @ UPp
    op = np.load(os.path.join(RES, "dispop_allnodes.npz")); co = op["coords"]
    Wall = np.einsum("nck,kj->ncj", op["D"], UPp)
    nn = lambda x: int(np.linalg.norm(co - np.array(x), axis=1).argmin())
    w = Wall[nn(A_XYZ), 2] - Wall[nn(O_XYZ), 2]
    tmask = np.ones(len(Cc), bool); tmask[pod] = False; tall = np.where(tmask)[0]
    far = (np.linalg.norm(co - np.array(A_XYZ), axis=1) > EXCL) & (np.linalg.norm(co - np.array(O_XYZ), axis=1) > EXCL)
    Wf = Wall.reshape(-1, NPT); uall = np.where(np.repeat(far, 3))[0]
    Bs = priors()
    sig0 = float(np.mean([np.sqrt(w @ B @ w) for B in Bs]))
    out = {"note": "反り A−O の σ（30〜300 秒の10時刻平均、µm）。4000本の事前分布、ノイズ 温度0.3 K・変位0.3 µm。データ同化はしていない。",
           "sigma_none_um": sig0, "modes": {}}
    for mode, quota in [("T2", {"T": 2, "U": 0}), ("U2", {"T": 0, "U": 2}), ("ANY2", {"T": 2, "U": 2})]:
        Ps = [B.copy() for B in Bs]; tidx, uidx = tall.copy(), uall.copy(); n = {"T": 0, "U": 0}; picked = []
        for step in range(2):
            best = None
            if n["T"] < quota["T"]:
                s = score(Ps, w, G[tidx], R_T); k = int(s.argmin()); best = ("T", int(tidx[k]), float(s[k]))
            if n["U"] < quota["U"]:
                s = score(Ps, w, Wf[uidx], R_U); k = int(s.argmin())
                if best is None or s[k] < best[2]: best = ("U", int(uidx[k]), float(s[k]))
            kind, idx, sv = best; n[kind] += 1
            if kind == "T":
                h, r = G[idx], R_T; tidx = tidx[tidx != idx]
                info = {"kind": "温度", "cell": idx, "xyz_mm": np.round(Cc[idx] * 1000, 1).tolist()}
            else:
                h, r = Wf[idx], R_U; uidx = uidx[uidx != idx]; i, c = divmod(idx, 3)
                info = {"kind": "変位", "node": i, "comp": COMP[c], "xyz_mm": np.round(co[i] * 1000, 1).tolist()}
            Ps = [add_row(P, h, r) for P in Ps]
            info.update(step=step + 1, sigma_um=sv, weights=np.round(h, 4).tolist()); picked.append(info)
            print(f"[{mode}] {step+1}本目: {info['kind']} {info['xyz_mm']} {info.get('comp','')}  σ = {sv:.3f} µm", flush=True)
        out["modes"][mode] = picked

    def sig_of(rows):
        vals = []
        for B in Bs:
            P = B.copy()
            for h, r in rows: P = add_row(P, h, r)
            vals.append(np.sqrt(max(w @ P @ w, 0)))
        return float(np.mean(vals))
    e = np.eye(NPT)
    out["reference_um"] = {"温度 P2・P0（代表点・片側）": sig_of([(e[2], R_T), (e[0], R_T)]),
                           "温度 P2・P4（代表点・両側）": sig_of([(e[2], R_T), (e[4], R_T)])}
    for k, v in out["reference_um"].items(): print(f"  σ {k}: {v:.3f} µm")
    json.dump(out, open(os.path.join(RES, "two_nonrep_selection.json"), "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
