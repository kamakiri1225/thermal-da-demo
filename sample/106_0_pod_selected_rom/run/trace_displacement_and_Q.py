"""blog_006 用：変位の同化への入れ方と、発熱量Qの推定を、本番計算の実数値で記録する.

run_da_compare.py の「温度2点(P2,P0)+変位2点(A,O)」、seed 20260913 をそのまま再実行し、
  ・変位の写像に使う行列（平均場の5点温度、U_rP^+、D、平均場の変位、熱感度W）
  ・1メンバーを例に、温度 → POD係数 → 変位 の計算
  ・サイクル1・5・10の、カルマンゲイン全体（7×4）と、qの更新の内訳（温度由来・変位由来）
  ・全20サイクルの Q の推移と、温度由来・変位由来の寄与の累計
を記録する。説明用の仮の値は使わない。

出力: results/trace_displacement_and_Q.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/trace_displacement_and_Q.py
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
TRACE=[1,5,10]
NODES=[2,0]                      # run_da_compare の温度2点（P2, P0 の順）


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); Tbar=mean[pod]
    op=np.load(os.path.join(RES,"disp_operator.npz")); ubar=op["uz_mean"]; D=op["Dmode"]
    W=D@UPp                                                     # 2×5：温度→変位 [µm/K]
    disp=lambda T5: ubar+((T5-Tbar)@UPp.T)@D.T
    out=dict(Tbar_K=Tbar.tolist(),UPp=UPp.tolist(),D=D.tolist(),ubar_um=ubar.tolist(),W_um_per_K=W.tolist(),
             W_AO=(W[0]-W[1]).tolist(),C=C.tolist(),h_true=h_true,heat_node=heat)
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    rng=np.random.default_rng(SEED); ro=np.random.default_rng(SEED+7)
    Z=np.zeros((N_ENS,NAUG))
    Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2)
    out["Q0_mean_W"]=float(15*Z[:,IQ].mean()); out["Q0_sd_W"]=float(15*Z[:,IQ].std(ddof=1))
    hist=[]; cum_T=0.0; cum_U=0.0; tp=0.0; out["cycles"]=[]
    for ci,tb in enumerate(cyc,1):
        Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); Z=Z.copy(); tp=tb
        if ci==1:
            m0=Z[0,:NPT]; dl=m0-Tbar; a=UPp@dl
            out["member1_example"]=dict(T5_K=m0.tolist(),delta_K=dl.tolist(),a=a.tolist(),Uz_um=disp(m0).tolist(),q=float(Z[0,IQ]))
            out["truth_t30"]=dict(T5_K=Ttr[1].tolist(),Uz_um=disp(Ttr[1]).tolist())
        Yf=np.column_stack([Z[:,NODES],disp(Z[:,:NPT])])
        y=np.array(list(Ttr[ci][NODES])+list(disp(Ttr[ci])))+ro.normal(0,np.sqrt(np.diag(Rd)))
        zb=Z.mean(0); yb=Yf.mean(0); Zi=zb+INFL*(Z-zb); Yi=yb+INFL*(Yf-yb)
        dZ=Zi-Zi.mean(0); dY=Yi-Yi.mean(0); Czy=dZ.T@dY/(N_ENS-1); Cyy=dY.T@dY/(N_ENS-1)
        Kg=np.linalg.solve(Cyy+Rd,Czy.T).T; res=y-yb
        contrib=Kg[IQ]*res; cT=float(contrib[:2].sum()); cU=float(contrib[2:].sum())
        qb=float(Z[:,IQ].mean())
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        qa=float(Z[:,IQ].mean()); cum_T+=15*cT; cum_U+=15*cU
        hist.append(dict(t_s=float(tb),Q_before_W=15*qb,Q_after_W=15*qa,Q_sd_after_W=float(15*Z[:,IQ].std(ddof=1)),
                         dQ_temp_W=15*cT,dQ_disp_W=15*cU))
        if ci in TRACE:
            out["cycles"].append(dict(cycle=ci,t_s=float(tb),y=y.tolist(),Yf_mean=yb.tolist(),residual=res.tolist(),
                Cov_q_Y=Czy[IQ].tolist(),Var_Y=np.diag(Cyy).tolist(),Var_q=float(np.var(dZ[:,IQ],ddof=1)),
                K=Kg.tolist(),dz_mean=(Kg@res).tolist(),Q_before_W=15*qb,Q_after_W=15*qa,dQ_temp_W=15*cT,dQ_disp_W=15*cU))
    out["history"]=hist; out["cum_dQ_temp_W"]=cum_T; out["cum_dQ_disp_W"]=cum_U; out["Q_final_W"]=float(15*Z[:,IQ].mean())
    json.dump(out,open(os.path.join(RES,"trace_displacement_and_Q.json"),"w"),ensure_ascii=False,indent=1)
    np.set_printoptions(precision=4,suppress=True)
    print("Tbar",np.round(Tbar,3)); print("UPp\n",np.round(UPp,4)); print("D\n",np.round(D,5)); print("ubar",ubar); print("W\n",np.round(W,3))
    print("member1",{k:np.round(v,4) if isinstance(v,list) else v for k,v in out["member1_example"].items()})
    print("truth30",out["truth_t30"])
    for c in out["cycles"]:
        print(f"--- cycle {c['cycle']} t={c['t_s']}: Q {c['Q_before_W']:.2f}->{c['Q_after_W']:.2f}  dQ_temp {c['dQ_temp_W']:+.3f} dQ_disp {c['dQ_disp_W']:+.3f}")
        print(" y",np.round(c['y'],3)," Yf",np.round(c['Yf_mean'],3)); print(" K\n",np.round(np.array(c['K']),4)); print(" dz",np.round(c['dz_mean'],4))
        print(" Cov(q,Y)",np.round(c['Cov_q_Y'],4)," VarY",np.round(c['Var_Y'],4)," Var q",round(c['Var_q'],4))
    for hh in hist: print(f"t={hh['t_s']:5.0f} Q {hh['Q_before_W']:6.2f}->{hh['Q_after_W']:6.2f} sd {hh['Q_sd_after_W']:.2f}  temp {hh['dQ_temp_W']:+.3f} disp {hh['dQ_disp_W']:+.3f}")
    print("cum temp",cum_T,"cum disp",cum_U,"final",out["Q_final_W"])


if __name__=="__main__": main()
