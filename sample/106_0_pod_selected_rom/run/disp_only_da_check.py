"""変位観測だけで温度場を直せるか（温度センサを1本も使わない同化）.

実機では内部に熱電対を入れられない場合がある。変位計だけで温度場と発熱量を
推定できるなら、センサ構成の選択肢が広がる。

比較する観測構成:
  ① 同化なし
  ② 温度2点のみ
  ③ 変位2点のみ（上面A/O）        ← 温度を1本も測らない
  ④ 変位2点のみ（中高さB/C）       ← 同上、別の場所
  ⑤ 変位4点（A/O + B/C）          ← 変位だけ増やす
  ⑥ 温度2点＋変位2点（本研究の標準）

評価: 全5点温度RMSE / 発熱量Q / 上面の反り

出力: docs/img/disp_only_da.png, results/disp_only_da.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/disp_only_da_check.py
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
N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
iAz,iOz,iBz,iCz=2,5,8,11


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T

    b0=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.0,heat,0,300,DT)[1][-1]
    b1=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.1,heat,0,300,DT)[1][-1]
    hi2=list(np.argsort((b1-b0)/0.1)[::-1][:2])

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr); QOI=Utr[:,iAz]-Utr[:,iOz]

    def run(tnodes, dcols, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        assim = (tnodes is not None) or (dcols is not None)
        if assim:
            Rd=np.diag([SIG_T**2]*(len(tnodes) if tnodes else 0)
                       +[SIG_U**2]*(len(dcols) if dcols else 0))
        rT=[]; rQ=[]; rB=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if assim:
                yv=[]; cols=[]
                if tnodes: yv+= list(Ttr[ci][tnodes]); cols.append(Z[:,tnodes])
                if dcols:  yv+= [Utr[ci][j] for j in dcols]; cols.append(disp(Z[:,:NPT])[:,dcols])
                Yf=np.column_stack(cols)
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            m=Z[:,:NPT].mean(0)
            if 60<=tb<=300:
                rT.append(np.sqrt(((m-Ttr[ci])**2).mean()))
                u=disp(m); rB.append(abs((u[iAz]-u[iOz])-QOI[ci]))
            rQ.append(Z[:,IQ].mean()*15)
        return float(np.mean(rT)), float(rQ[-1]), float(np.mean(rB))

    cfgs=[("① 同化なし",              None,  None,          "0.45"),
          ("② 温度2点のみ",           hi2,   None,          "tab:orange"),
          ("③ 変位2点のみ（上面A/O）",  None,  [iAz,iOz],     "tab:red"),
          ("④ 変位2点のみ（中高さB/C）",None,  [iBz,iCz],     "tab:purple"),
          ("⑤ 変位4点（A/O+B/C）",     None,  [iAz,iOz,iBz,iCz],"tab:brown"),
          ("⑥ 温度2点＋変位2点",       hi2,   [iAz,iOz],     "tab:blue")]
    rows=[]
    print(f"{'観測構成':<28}{'温度RMSE':>10}{'最終Q':>10}{'反りの誤差':>12}")
    for name,tn,dc,col in cfgs:
        r=[run(tn,dc,s) for s in SEEDS]
        T_=np.mean([x[0] for x in r]); Q_=np.mean([x[1] for x in r]); B_=np.mean([x[2] for x in r])
        rows.append(dict(name=name,color=col,rmseT=T_,Q=Q_,bow=B_,
                         n_temp=len(tn) if tn else 0, n_disp=len(dc) if dc else 0))
        print(f"{name:<28}{T_:8.3f} K{Q_:9.2f} W{B_:10.3f} µm")
    print(f"\n真値 Q = 15.00 W")

    fig,axes=plt.subplots(1,3,figsize=(18.0,5.8))
    names=[r["name"] for r in rows]; cols=[r["color"] for r in rows]
    for ax,(key,ttl,unit,ref) in zip(axes,
        [("rmseT","全5点の温度RMSE","[K]",None),
         ("Q","最終の発熱量 $Q$","[W]",15.0),
         ("bow","上面の反りの誤差","[µm]",None)]):
        v=[r[key] for r in rows]
        b=ax.bar(range(len(rows)),v,color=cols)
        for r_,x in zip(b,v):
            ax.text(r_.get_x()+r_.get_width()/2,x,f"{x:.3f}" if key!="Q" else f"{x:.2f}",
                    ha="center",va="bottom",fontsize=11,weight="bold")
        if ref: ax.axhline(ref,color="k",ls="--",lw=2.2,label=f"真値 {ref:.0f} W"); ax.legend(fontsize=10.5)
        else: ax.set_yscale("log")
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([n.split("（")[0].replace("①","①\n").replace("②","②\n")
                            .replace("③","③\n").replace("④","④\n").replace("⑤","⑤\n")
                            .replace("⑥","⑥\n") for n in names],fontsize=9.5)
        ax.set_ylabel(f"{ttl} {unit}",fontsize=12)
        ax.set_title(ttl,fontsize=13.5,weight="bold"); ax.grid(alpha=.3,axis="y",which="both")

    fig.suptitle("温度センサを1本も使わず、変位計だけで温度場を直せるか\n"
                 "③④＝変位だけ／②＝温度だけ／⑥＝両方（本研究の標準）。60メンバー・5 seed平均・加熱期60–300 s",
                 fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.86])
    o=os.path.join(IMG,"disp_only_da.png"); fig.savefig(o,dpi=130); plt.close(fig)
    json.dump(dict(Q_true=15.0,window="60-300s",rows=rows),
              open(os.path.join(RES,"disp_only_da.json"),"w"),ensure_ascii=False,indent=2)
    print("wrote",o)


if __name__=="__main__":
    main()
