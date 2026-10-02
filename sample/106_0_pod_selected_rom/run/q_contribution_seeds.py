"""発熱量Qの推定で、温度センサと変位計がそれぞれどれだけQを動かしたかを、5 seed で調べる.

run_da_compare.py の「温度2点(P2,P0)+変位2点(A,O)」と同じ設定。各サイクルの平均更新
  ΔQ = 15 Σ_j K_{q,j} d_j
を温度2列・変位2列に分け、seed ごとに合計する。あわせて、Qの誤差 |Q−15| がどちらの寄与で減ったかも見る。

出力: results/q_contribution_seeds.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/q_contribution_seeds.py
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
N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]; NODES=[2,0]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); op=np.load(os.path.join(RES,"disp_operator.npz")); ub=op["uz_mean"]; D=op["Dmode"]
    disp=lambda T5: ub+((T5-mean[pod])@UPp.T)@D.T
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2); out=[]
    for seed in SEEDS:
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Q0=15*Z[:,IQ].mean(); sT=sU=0.0; gT=gU=0.0; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); Z=Z.copy(); tp=tb
            Yf=np.column_stack([Z[:,NODES],disp(Z[:,:NPT])])
            y=np.array(list(Ttr[ci][NODES])+list(disp(Ttr[ci])))+ro.normal(0,np.sqrt(np.diag(Rd)))
            zb=Z.mean(0); yb=Yf.mean(0); Zi=zb+INFL*(Z-zb); Yi=yb+INFL*(Yf-yb)
            dZ=Zi-Zi.mean(0); dY=Yi-Yi.mean(0); Czy=dZ.T@dY/(N_ENS-1); Cyy=dY.T@dY/(N_ENS-1)
            K=np.linalg.solve(Cyy+Rd,Czy.T).T; c=15*K[IQ]*(y-yb)
            qb=15*Z[:,IQ].mean()
            # 「真値15 Wへ近づけた量」：補正の向きが 15−Q と同じなら正
            sgn=np.sign(15-qb)
            sT+=c[:2].sum(); sU+=c[2:].sum(); gT+=sgn*c[:2].sum(); gU+=sgn*c[2:].sum()
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        r=dict(seed=seed,Q0_W=float(Q0),Q_final_W=float(15*Z[:,IQ].mean()),sum_temp_W=float(sT),sum_disp_W=float(sU),
               toward_truth_temp_W=float(gT),toward_truth_disp_W=float(gU))
        out.append(r); print(r,flush=True)
    m={k:float(np.mean([o[k] for o in out])) for k in out[0] if k!="seed"}
    print("mean",m)
    json.dump(dict(note="ΔQの寄与 [W]。sum=符号つき合計、toward_truth=各サイクルで真値15 Wへ向かう向きを正とした合計",seeds=out,mean=m),
              open(os.path.join(RES,"q_contribution_seeds.json"),"w"),ensure_ascii=False,indent=1)


if __name__=="__main__": main()
