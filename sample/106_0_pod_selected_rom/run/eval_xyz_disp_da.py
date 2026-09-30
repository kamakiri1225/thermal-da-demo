"""x・y・z の3成分すべてで、観測点と未観測点の変位がデータ同化でどうなるかを検証する.

これまでは z 成分（上下・軸方向）だけを見ていた。ここでは
  観測点   A(+28,0,100.5) / O(-28,0,100.5) mm   ← 同化に使う（ただし Uz のみ観測）
  未観測点 B(+28,0, 75.0) / C(-28,0, 75.0) mm   ← 同化に一切使わない
の Ux, Uy, Uz を、同化なし / 同化後 で真値と比べる。

狙い：「Uz しか観測していないのに、Ux・Uy も直るのか」「未観測点でも直るのか」を
データで示す。温度場が正しくなれば W を通じて全成分が従うはずだが、確認していなかった。

出力: docs/img/xyz_disp_da.png, results/xyz_disp_da.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/eval_xyz_disp_da.py
"""
from __future__ import annotations
import os, sys, json
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
# 観測点（上面）と未観測点（中高さ z=75mm）。どちらも +X / -X の対。
PTS=[("A",( 0.028,0,H    ),"観測"),("O",(-0.028,0,H    ),"観測"),
     ("B",( 0.028,0,0.075),"未観測"),("C",(-0.028,0,0.075),"未観測")]
COMP=["Ux","Uy","Uz"]


def build_op(U, mean, Cc):
    """全4点×3成分の変位演算子。 u = u_mean + D @ a   (u: 12成分)"""
    cache=os.path.join(RES,"dispop_xyz_4pts.npz")
    if os.path.exists(cache):
        d=np.load(cache); return d["u_mean"], d["D"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([p for _n,p in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    idx=[int(np.linalg.norm(coords-np.array(x),axis=1).argmin()) for _l,x,_r in PTS]
    for (lab,x,_r),j in zip(PTS,idx):
        print(f"  {lab} 指定{np.round(np.array(x)*1000,1)} → 最寄節点{np.round(coords[j]*1000,1)} mm")
    _,near=cKDTree(Cc).query(coords)
    work=os.path.join(ROOT,"openfoam","dispop_xyz"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))

    def run_field(Tfield):
        fistr_case.write_cnt(Path(work),node_ids,Tfield,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        return np.concatenate([np.array(disp[node_ids[j]])*1e6 for j in idx])   # (12,) µm

    r=U.shape[1]
    u_mean=run_field(mean[near]); D=np.zeros((12,r))
    for k in range(r):
        D[:,k]=run_field(Tref+U[near,k])-0.0
        print(f"  [disp] mode {k+1}/{r} FrontISTR 完了",flush=True)
    # モード応答は「平均場からの増分」ではなく絶対値なので、基準を引く
    zero=run_field(np.full(len(near),Tref))
    D=D-zero[:,None]; u_mean=u_mean-zero
    np.savez(cache,u_mean=u_mean,D=D,idx=np.array(idx))
    return u_mean, D


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"])
    heat_node=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float)
    Cc=kv["cell_centres"]; pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:])
    u_mean,D=build_op(U,mean,Cc)

    def disp_of_T5(T5):
        a=(T5-mean[pod])@UPp.T
        return u_mean + a@D.T                     # (...,12)

    # 熱感度Wで温度2点を選ぶ（本研究と同じ基準：A-O の z 変位差）
    iAz=2; iOz=5                                   # A の Uz / O の Uz
    wAO=(D[iAz]-D[iOz])@UPp
    hi2=list(np.argsort(np.abs(wAO))[::-1][:2])
    print(f"温度観測点 = P{hi2}（熱感度W上位2）")

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr); Utr=disp_of_T5(Ttr)         # 真値 (nt,12)

    def run(assimilate, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.diag([SIG_T**2]*2+[SIG_U**2]*2)
        out=[disp_of_T5(Z[:,:NPT].mean(0))]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy()
            Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
            tp=tb
            if assimilate:
                # 観測は「温度2点」＋「A/O の Uz だけ」。x,y は一切観測しない
                yv=list(Ttr[ci][hi2])+[Utr[ci][iAz],Utr[ci][iOz]]
                Yf=np.column_stack([Z[:,hi2],
                                    disp_of_T5(Z[:,:NPT])[:,[iAz,iOz]]])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            out.append(disp_of_T5(Z[:,:NPT].mean(0)))
        return np.array(out)

    free=np.mean([run(False,s) for s in SEEDS],axis=0)
    da  =np.mean([run(True ,s) for s in SEEDS],axis=0)

    heat=(tg>0)&(tg<=300)
    heat2=(tg>=60)&(tg<=300)      # 最初の1サイクル(t=30s)を除いた窓
    rows=[]
    for p,(lab,_x,role) in enumerate(PTS):
        for c in range(3):
            k=p*3+c
            e0=float(np.abs(free[heat,k]-Utr[heat,k]).mean())
            e1=float(np.abs(da  [heat,k]-Utr[heat,k]).mean())
            f0=float(np.abs(free[heat2,k]-Utr[heat2,k]).mean())
            f1=float(np.abs(da  [heat2,k]-Utr[heat2,k]).mean())
            rows.append(dict(point=lab,role=role,comp=COMP[c],
                             amp_um=float(np.abs(Utr[:,k]).max()),
                             mae_noDA_um=e0,mae_DA_um=e1,
                             improve=float(e0/e1) if e1>0 else float("inf"),
                             mae_noDA_after1cyc_um=f0,mae_DA_after1cyc_um=f1,
                             improve_after1cyc=float(f0/f1) if f1>0 else float("inf")))
    print(f"\n{'点':<4}{'区分':<6}{'成分':<5}{'振幅':>8}"
          f"{'--- 0-300s ---':>22}{'--- 60-300s(初回除く) ---':>28}")
    print(f"{'':<15}{'':>8}{'なし':>10}{'同化':>9}{'改善':>7}{'なし':>12}{'同化':>9}{'改善':>8}")
    for r in rows:
        print(f"{r['point']:<4}{r['role']:<6}{r['comp']:<5}{r['amp_um']:8.3f}"
              f"{r['mae_noDA_um']:10.3f}{r['mae_DA_um']:9.3f}{r['improve']:6.1f}x"
              f"{r['mae_noDA_after1cyc_um']:12.3f}{r['mae_DA_after1cyc_um']:9.3f}"
              f"{r['improve_after1cyc']:7.1f}x")

    # ── 図 ──
    fig,axes=plt.subplots(3,4,figsize=(19.2,10.4),sharex=True)
    for p,(lab,x,role) in enumerate(PTS):
        for c in range(3):
            ax=axes[c,p]; k=p*3+c
            ax.axvspan(0,300,color="orange",alpha=.07)
            ax.plot(tg,Utr[:,k],color="black",lw=4.6,alpha=.40,label="真値")
            ax.plot(tg,free[:,k],":",color="#7a8899",lw=2.2,label="同化なし")
            ax.plot(tg,da[:,k],"--",color="#C0392B",lw=2.3,dashes=(4,3),label="同化後")
            r=rows[k]
            ax.set_title(f"{lab}（{role}） {COMP[c]}\n"
                         f"MAE {r['mae_noDA_um']:.3f} → {r['mae_DA_um']:.3f} µm"
                         f"（{r['improve']:.0f}倍）",fontsize=11.5,
                         weight="bold",color=("#C0392B" if role=="未観測" else "#1F4E9C"))
            ax.grid(alpha=.3); ax.tick_params(labelsize=9.5)
            if c==2: ax.set_xlabel("時刻 [s]",fontsize=11)
            if p==0: ax.set_ylabel(f"{COMP[c]} [µm]",fontsize=12)
            if p==0 and c==0: ax.legend(fontsize=10,loc="upper left")
    fig.suptitle("観測は「温度2点＋A/O の $U_z$ だけ」― それでも $U_x,U_y$ と未観測点は直るか\n"
                 f"A/O＝上面 z=100.5 mm（観測に使用）　B/C＝中高さ z=75 mm（一切観測しない）　"
                 f"5 seed平均・加熱期 0–300 s のMAE",fontsize=15,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.915])
    out=os.path.join(IMG,"xyz_disp_da.png")
    fig.savefig(out,dpi=120); plt.close(fig)
    json.dump(dict(points=[dict(label=l,xyz_mm=[v*1000 for v in x],role=r) for l,x,r in PTS],
                   obs="温度2点(熱感度W上位2) + A/O の Uz のみ",
                   window="加熱期 0-300 s、5 seed平均", rows=rows),
              open(os.path.join(RES,"xyz_disp_da.json"),"w"),ensure_ascii=False,indent=2)
    print("\nwrote",out)


if __name__=="__main__":
    main()
