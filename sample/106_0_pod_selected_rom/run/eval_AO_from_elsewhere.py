"""知りたいのは A/O の変位差。それを「A/Oを測らずに」当てられるかを検証する.

本当に知りたい量（評価点）: 上面の反り Uz(A) - Uz(O)   z=100.5 mm
  実機では工具先端側にあたり、加工中は変位計を置けない場所。

観測構成を4つ比べる:
  ① 同化なし
  ② 温度2点のみ                     ← 変位は一切測らない
  ③ 温度2点 + 変位B/C（中高さz=75mm） ← 別の場所の変位を測る。循環なし
  ④ 温度2点 + 変位A/O                ← 知りたい場所そのものを測る（参考・循環あり）

出力: docs/img/AO_from_elsewhere.png, results/AO_from_elsewhere.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/eval_AO_from_elsewhere.py
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
iAz,iOz,iBz,iCz = 2,5,8,11          # dispop_xyz_4pts の Uz 行


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz")); um=op["u_mean"]; D=op["D"]
    def disp(T5): return um+((T5-mean[pod])@UPp.T)@D.T

    # 温度観測点は発熱感度 dT/dQ の上位2点
    b0=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.0,heat,0,300,DT)[1][-1]
    b1=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Km,h_true,1.1,heat,0,300,DT)[1][-1]
    hi2=list(np.argsort((b1-b0)/0.1)[::-1][:2])
    print(f"温度観測点 = P{hi2}（発熱感度 上位2）")

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h_true,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp(Ttr)
    QOI_true=Utr[:,iAz]-Utr[:,iOz]          # 知りたい量：上面の反り

    def run(disp_obs, seed):
        """disp_obs: None / (iBz,iCz) / (iAz,iOz)"""
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        nobs=2+(2 if disp_obs else 0)
        Rd=np.diag([SIG_T**2]*2+([SIG_U**2]*2 if disp_obs else []))
        qoi=[]; uA=[]; uO=[]; tp=0.0
        u0=disp(Z[:,:NPT].mean(0)); qoi.append(u0[iAz]-u0[iOz]); uA.append(u0[iAz]); uO.append(u0[iOz])
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            yv=list(Ttr[ci][hi2]); Yf=Z[:,hi2]
            if disp_obs:
                i1,i2=disp_obs
                yv=yv+[Utr[ci][i1],Utr[ci][i2]]
                Yf=np.column_stack([Yf,disp(Z[:,:NPT])[:,[i1,i2]]])
            y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            u=disp(Z[:,:NPT].mean(0)); qoi.append(u[iAz]-u[iOz]); uA.append(u[iAz]); uO.append(u[iOz])
        return np.array(qoi),np.array(uA),np.array(uO)

    cfgs=[("同化なし",              None,        "0.45"),
          ("温度2点のみ",           None,        "tab:orange"),   # 同化はする
          ("温度2点＋変位B/C（中高さ）",(iBz,iCz), "tab:red"),
          ("温度2点＋変位A/O（参考）", (iAz,iOz),  "tab:blue")]

    res={}; heatw=(tg>=60)&(tg<=300)
    for k,(name,dobs,col) in enumerate(cfgs):
        if k==0:   # 同化なし＝自由予測
            def free(seed):
                rng=np.random.default_rng(seed)
                Z=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
                q=rng.uniform(0.3,1.8,N_ENS); hh=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
                out=[];uA=[];uO=[];tp=0.0
                u0=disp(Z.mean(0)); out.append(u0[iAz]-u0[iOz]); uA.append(u0[iAz]); uO.append(u0[iOz])
                for tb in cyc:
                    Z=rg.integrate_ensemble(Z,C,Km,hh,q,heat,tp,tb,DT); tp=tb
                    u=disp(Z.mean(0)); out.append(u[iAz]-u[iOz]); uA.append(u[iAz]); uO.append(u[iOz])
                return np.array(out),np.array(uA),np.array(uO)
            rs=[free(s) for s in SEEDS]
        else:
            rs=[run(dobs,s) for s in SEEDS]
        q=np.mean([r[0] for r in rs],axis=0)
        a=np.mean([r[1] for r in rs],axis=0); o=np.mean([r[2] for r in rs],axis=0)
        res[name]=dict(qoi=q.tolist(),uA=a.tolist(),uO=o.tolist(),color=col,
            mae_qoi=float(np.abs(q[heatw]-QOI_true[heatw]).mean()),
            mae_uA=float(np.abs(a[heatw]-Utr[heatw,iAz]).mean()),
            mae_uO=float(np.abs(o[heatw]-Utr[heatw,iOz]).mean()))

    print(f"\n知りたい量＝上面の反り Uz(A)-Uz(O)（真値の振幅 {np.abs(QOI_true).max():.3f} µm）")
    print(f"{'観測構成':<26}{'反りの誤差':>11}{'Uz(A)誤差':>11}{'Uz(O)誤差':>11}")
    for n,v in res.items():
        print(f"{n:<26}{v['mae_qoi']:9.3f} µm{v['mae_uA']:9.3f} µm{v['mae_uO']:9.3f} µm")

    fig,(a0,a1,a2)=plt.subplots(1,3,figsize=(18.6,5.8))
    for ax,(key,ttl,tru) in zip([a0,a1,a2],
        [("uA","$U_z(A)$  上面ヒータ側",Utr[:,iAz]),
         ("uO","$U_z(O)$  上面反対側",Utr[:,iOz]),
         ("qoi","反り $U_z(A)-U_z(O)$  ← 知りたい量",QOI_true)]):
        ax.axvspan(0,300,color="orange",alpha=.07)
        ax.plot(tg,tru,color="black",lw=5.0,alpha=.40,label="真値")
        for n,v in res.items():
            ls="--" if "B/C" in n else (":" if "A/O" in n else "-")
            ax.plot(tg,v[key],ls,color=v["color"],lw=2.4,label=n)
        ax.set_xlabel("時刻 [s]",fontsize=12); ax.set_ylabel("変位 [µm]",fontsize=12)
        ax.set_title(ttl,fontsize=13.5,weight="bold"); ax.grid(alpha=.3)
    a0.legend(fontsize=10,loc="upper left")

    fig.suptitle("知りたいのは上面の反り $U_z(A)-U_z(O)$ ― それを A/O を測らずに当てられるか\n"
                 f"評価点＝A/O（上面 z=100.5 mm）。"
                 f"赤破線＝別の場所（中高さ B/C）の変位を観測した構成。5 seed平均",
                 fontsize=15,weight="bold")
    fig.text(0.5,0.015,
        "反りの誤差（加熱期60–300 s のMAE）：" +
        " ／ ".join(f"{n} {v['mae_qoi']:.3f} µm" for n,v in res.items()),
        ha="center",fontsize=12,color="#333")
    fig.tight_layout(rect=[0,0.045,1,0.855])
    o=os.path.join(IMG,"AO_from_elsewhere.png"); fig.savefig(o,dpi=130); plt.close(fig)
    json.dump(dict(qoi="Uz(A)-Uz(O) at z=100.5mm",
                   qoi_amplitude_um=float(np.abs(QOI_true).max()),
                   window="60-300s MAE, 5 seed",
                   results={n:{k:v[k] for k in ("mae_qoi","mae_uA","mae_uO")} for n,v in res.items()}),
              open(os.path.join(RES,"AO_from_elsewhere.json"),"w"),ensure_ascii=False,indent=2)
    print("wrote",o)


if __name__=="__main__":
    main()
