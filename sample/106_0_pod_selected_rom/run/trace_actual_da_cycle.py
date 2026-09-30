"""本研究の本番計算（run_da_compare.py の「温度2点+変位2点」、seed 20260913）を
そのまま再実行し、指定サイクルの更新で実際に使われた数値を記録する.

記録するもの: 観測 y、予報観測の平均、残差、Kalmanゲインの q 行、更新前後の q と Q。
説明用の仮の値は使わない。

出力: results/trace_actual_da_cycle.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/trace_actual_da_cycle.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02; SEED=20260913
TRACE_CYCLES=[1,5,10]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"disp_operator.npz")); uz_mean=op["uz_mean"]; Dm=op["Dmode"]
    def disp(T5): return uz_mean+((T5-mean[pod])@UPp.T)@Dm.T      # run_da_compare と同一
    b0=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.0,heat,0,300,DT)[1][-1]
    b1=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.1,heat,0,300,DT)[1][-1]
    hi2=list(np.argsort((b1-b0)/0.1)[::-1][:2])

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    rng=np.random.default_rng(SEED); rng_o=np.random.default_rng(SEED+7)
    Z=np.zeros((N_ENS,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2)
    out={"obs_temp_nodes":[f"P{i}" for i in hi2],"Q0_mean_W":float(Z[:,IQ].mean()*15),"cycles":[]}
    tp=0.0
    for ci,tb in enumerate(cyc,1):
        Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); Z=Z.copy(); tp=tb
        yv=list(Ttr[ci][hi2])+list(disp(Ttr[ci]))
        Yf=np.column_stack([Z[:,hi2],disp(Z[:,:NPT])])
        y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
        if ci in TRACE_CYCLES:
            # enkf_update と同じ膨張をかけてからゲインを計算（記録用。更新そのものは enkf_update が行う）
            zb=Z.mean(0); yb=Yf.mean(0)
            Zi=zb+INFL*(Z-zb); Yi=yb+INFL*(Yf-yb)
            dZ=Zi-Zi.mean(0); dY=Yi-Yi.mean(0)
            Czy=dZ.T@dY/(N_ENS-1); Cyy=dY.T@dY/(N_ENS-1)
            Kg=np.linalg.solve(Cyy+Rd,Czy.T).T
            q_before=float(Z[:,IQ].mean())
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        if ci in TRACE_CYCLES:
            res=y-Yf.mean(0)
            out["cycles"].append(dict(
                cycle=ci,t_s=float(tb),
                y=[float(v) for v in y], Yf_mean=[float(v) for v in Yf.mean(0)],
                residual=[float(v) for v in res],
                Cov_q_Y=[float(v) for v in Czy[IQ]],
                K_q=[float(v) for v in Kg[IQ]],
                K_q_times_res=[float(v) for v in Kg[IQ]*res],
                q_before=q_before, q_after=float(Z[:,IQ].mean()),
                Q_before_W=q_before*15, Q_after_W=float(Z[:,IQ].mean()*15)))
    out["Q_final_W"]=float(Z[:,IQ].mean()*15)
    json.dump(out,open(os.path.join(RES,"trace_actual_da_cycle.json"),"w"),ensure_ascii=False,indent=2)
    lab=["T(P%s)[K]"%hi2[0][-1] if False else f"T(P{hi2[0]})",f"T(P{hi2[1]})","Uz(A)","Uz(O)"]
    print("観測温度点:",out["obs_temp_nodes"],"  初期 Q平均 =",f"{out['Q0_mean_W']:.2f} W")
    for c in out["cycles"]:
        print(f"\n--- サイクル{c['cycle']}（t={c['t_s']:.0f} s） ---")
        print(f"{'観測':<8}{'y':>9}{'Yf平均':>9}{'残差':>9}{'Cov(q,·)':>10}{'K_q':>10}{'寄与':>10}")
        for i in range(4):
            print(f"{lab[i]:<8}{c['y'][i]:9.3f}{c['Yf_mean'][i]:9.3f}{c['residual'][i]:+9.3f}"
                  f"{c['Cov_q_Y'][i]:10.4f}{c['K_q'][i]:+10.4f}{c['K_q_times_res'][i]:+10.4f}")
        print(f"q: {c['q_before']:.4f} → {c['q_after']:.4f}   Q: {c['Q_before_W']:.2f} → {c['Q_after_W']:.2f} W"
              f"   （Σ寄与×15 = {sum(c['K_q_times_res'])*15:+.2f} W）")
    print(f"\n最終 Q = {out['Q_final_W']:.2f} W（真値15）")


if __name__=="__main__":
    main()
