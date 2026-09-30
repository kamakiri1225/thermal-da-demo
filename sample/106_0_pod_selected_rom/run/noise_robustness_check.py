"""観測ノイズを振って、本研究の結論が「その条件でしか成り立たない」かを検証する.

本研究はすべて σ_T=0.30 K, σ_u=0.30 µm の1条件で回している。
実機の熱電対は ±1 K 級もあるので、結論がノイズ水準に依存しないかを確かめる。

検証する結論:
  (1) 観測を増やすほど良くなる（温度1点 → 2点 → +変位2点）
  (2) 変位を足すと効く
  (3) 熱感度Wの高い温度点のほうが良い

出力: docs/img/noise_robustness.png, results/noise_robustness.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/noise_robustness_check.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; INFL=1.02; SEEDS=[20260913,20260914,20260915]
SIG_T_LIST=[0.1,0.3,1.0,3.0]      # K
SIG_U_LIST=[0.1,0.3,1.0]          # µm


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    iAz,iOz=2,5
    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T

    wAO=(D[iAz]-D[iOz])@UPp; sens=np.abs(wAO)
    hi=int(np.argmax(sens)); lo=int(np.argmin(sens)); hi2=list(np.argsort(sens)[::-1][:2])

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr)

    def run(nodes,use_disp,sT,sU,seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        if nodes is not None:
            Rd=np.diag([sT**2]*len(nodes)+([sU**2]*2 if use_disp else []))
        rec=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if nodes is not None:
                yv=list(Ttr[ci][nodes]); Yf=Z[:,nodes]
                if use_disp:
                    yv=yv+[Utr[ci][iAz],Utr[ci][iOz]]
                    Yf=np.column_stack([Yf,disp(Z[:,:NPT])[:,[iAz,iOz]]])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                rec.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean()))
        return float(np.mean(rec))

    cfgs=[("同化なし",None,False),(f"温度1点(W低 P{lo})",[lo],False),
          (f"温度1点(W高 P{hi})",[hi],False),(f"温度2点(P{hi2[0]}+P{hi2[1]})",hi2,False),
          ("温度2点+変位2点",hi2,True)]
    out={}
    print(f"{'σ_T [K]':>8}{'σ_u [µm]':>10}  " + "".join(f"{c[0]:>20}" for c in cfgs))
    for sT in SIG_T_LIST:
        for sU in SIG_U_LIST:
            if sU!=0.3 and sT!=0.3: continue      # 片方ずつ振る（計算量を抑える）
            vals=[]
            for name,nodes,ud in cfgs:
                v=np.mean([run(nodes,ud,sT,sU,s) for s in SEEDS])
                vals.append(float(v))
            out[f"sT={sT}_sU={sU}"]=dict(sig_T=sT,sig_u=sU,
                                         configs=[c[0] for c in cfgs],rmse_K=vals)
            print(f"{sT:8.1f}{sU:10.1f}  " + "".join(f"{v:20.4f}" for v in vals))

    # ── 図 ──
    fig,(a0,a1)=plt.subplots(1,2,figsize=(15.6,6.0))
    keysT=[k for k in out if out[k]["sig_u"]==0.3]
    xs=[out[k]["sig_T"] for k in keysT]
    cols=["0.45","tab:orange","tab:blue","tab:green","tab:red"]
    for i,(name,_n,_u) in enumerate(cfgs):
        a0.plot(xs,[out[k]["rmse_K"][i] for k in keysT],"-o",lw=2.4,ms=7,
                color=cols[i],label=name)
    a0.set_xscale("log"); a0.set_yscale("log")
    a0.set_xlabel("温度センサのノイズ $\\sigma_T$ [K]（対数軸）",fontsize=12.5)
    a0.set_ylabel("加熱期の全5点温度RMSE [K]",fontsize=12.5)
    a0.set_title("① 温度ノイズを振る（$\\sigma_u$=0.3 µm 固定）\n"
                 "順位が入れ替わらなければ結論は頑健",fontsize=13,weight="bold")
    a0.legend(fontsize=10.5); a0.grid(alpha=.3,which="both")

    keysU=[k for k in out if out[k]["sig_T"]==0.3]
    xu=[out[k]["sig_u"] for k in keysU]
    for i,(name,_n,_u) in enumerate(cfgs):
        a1.plot(xu,[out[k]["rmse_K"][i] for k in keysU],"-o",lw=2.4,ms=7,
                color=cols[i],label=name)
    a1.set_xscale("log"); a1.set_yscale("log")
    a1.set_xlabel("変位センサのノイズ $\\sigma_u$ [µm]（対数軸）",fontsize=12.5)
    a1.set_ylabel("加熱期の全5点温度RMSE [K]",fontsize=12.5)
    a1.set_title("② 変位ノイズを振る（$\\sigma_T$=0.3 K 固定）\n"
                 "赤が緑を下回り続けるなら「変位は効く」は頑健",fontsize=13,weight="bold")
    a1.legend(fontsize=10.5); a1.grid(alpha=.3,which="both")

    fig.suptitle("結論は「その条件でしか成り立たない」のか ― 観測ノイズを振って確かめる\n"
                 "本研究の標準条件は $\\sigma_T$=0.30 K, $\\sigma_u$=0.30 µm（3 seed平均）",
                 fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.885])
    o=os.path.join(IMG,"noise_robustness.png"); fig.savefig(o,dpi=130); plt.close(fig)
    json.dump(out,open(os.path.join(RES,"noise_robustness.json"),"w"),
              ensure_ascii=False,indent=2)
    print("\nwrote",o)


if __name__=="__main__":
    main()
