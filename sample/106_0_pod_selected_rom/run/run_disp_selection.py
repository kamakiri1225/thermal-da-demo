"""変位観測点の熱感度(W行)による選定を検証: 高W変位点 vs 低W変位点 vs 変位なし.

温度センサは高感度ノード(P2)に固定し、追加する変位2点を
  ・高W（105の W=K⁻¹H 行感度が高い上端2点）
  ・低W（行感度が低い底面2点）
で変えて同化。全5点温度の推定がどう違うかを時刻歴・per-node・棒グラフで示す。

出力(docs/img/): disp_sel_traj.png, disp_sel_pernode.png, disp_sel_rmse.png
再現: OMP_NUM_THREADS=4 python3 run/run_disp_selection.py
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from dacore import rom_general as rg
from dacore.enkf import enkf_update
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
K=273.15; NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2; DT=2.0; OBS_DT=30.0; T_END=600.0
N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02; SEEDS=[20260913,20260914,20260915,20260916,20260917]
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005; Tref=MATERIAL["reference_temperature_K"]
KINVH=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","kinvh_sensitivity.npz")


def build_disp_op(U, mean, Cc, fem_nodes, tag):
    """指定FEM節点(list)での uz を、PODモード応答 D で観測化する.
    戻り: uz_mean(m,), D(m,r) で uz = uz_mean + D@a."""
    cache=os.path.join(RES,f"dispop_{tag}.npz")
    if os.path.exists(cache):
        d=np.load(cache); return d["uz_mean"], d["D"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    node_ids=[n for n,_ in mesh["nodes"]]; coords=np.array([xyz for _n,xyz in mesh["nodes"]])
    _,near=cKDTree(Cc).query(coords)
    work=os.path.join(ROOT,"openfoam","dispop_"+tag); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    def run_field(Tfield):
        fistr_case.write_cnt(Path(work),node_ids,Tfield,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        return np.array([disp[node_ids[i]][2] for i in fem_nodes])*1e6
    r=U.shape[1]; m=len(fem_nodes)
    uz_mean=run_field(mean[near]); D=np.zeros((m,r))
    for k in range(r):
        D[:,k]=run_field(Tref+U[near,k]); print(f"[disp:{tag}] mode {k+1}/{r}",flush=True)
    np.savez(cache,uz_mean=uz_mean,D=D); return uz_mean, D


def main():
    import yaml
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Kmat=rg.tri_to_matrix(d["K_upper"],NPT); heat_node=int(d["heat_node"]); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]; pod_cells=kv["cell_idx"]
    UP_pinv=np.linalg.pinv(U[pod_cells,:])
    kh=np.load(KINVH); hiN=[int(i) for i in kh["hi"]]; loN=[int(i) for i in kh["lo"]]
    uz_hi,D_hi=build_disp_op(U,mean,Cc,hiN,"hiW")
    uz_lo,D_lo=build_disp_op(U,mean,Cc,loN,"loW")

    # dT/dQ で温度センサ(高感度)
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    Tsens=int(np.argmax((pert-base)/0.1))

    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tgrid=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    def dispA(T5,uzm,D): a=(T5-mean[pod_cells])@UP_pinv.T; return uzm+a@D.T

    def run(dispkind,seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        m=0 if dispkind is None else (uz_hi.size if dispkind=="hi" else uz_lo.size)
        Rd=np.diag([SIG_T**2]+[SIG_U**2]*m)
        rmse=[np.sqrt(((Z[:,:NPT].mean(0)-Ttr[0])**2).mean())]; recT=[Z[:,:NPT].mean(0).copy()]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Tn=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
            Z=Z.copy(); Z[:,:NPT]=Tn; tp=tb
            yv=[Ttr[ci][Tsens]]; Yf=Z[:,Tsens:Tsens+1]
            if dispkind=="hi":
                yv+=list(dispA(Ttr[ci],uz_hi,D_hi)); Yf=np.column_stack([Yf,dispA(Z[:,:NPT],uz_hi,D_hi)])
            elif dispkind=="lo":
                yv+=list(dispA(Ttr[ci],uz_lo,D_lo)); Yf=np.column_stack([Yf,dispA(Z[:,:NPT],uz_lo,D_lo)])
            y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            rmse.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean())); recT.append(Z[:,:NPT].mean(0).copy())
        return np.array(rmse),np.array(recT)

    cfgs=[("温度1点のみ(変位なし)",None,"tab:gray"),
          ("＋変位2点(低熱感度W)","lo","tab:orange"),
          ("＋変位2点(高熱感度W)","hi","tab:blue")]
    res={}
    for name,dk,col in cfgs:
        rs=[run(dk,s) for s in SEEDS]
        res[name]=(np.mean([r[0] for r in rs],axis=0), rs[0][1], col)
        print(f"[dsel] {name:22s} 過渡RMSE={res[name][0][(tgrid>0)&(tgrid<=300)].mean():.3f}K 最終={res[name][0][-1]:.4f}K")

    # 図1: RMSE時刻歴
    fig,ax=plt.subplots(figsize=(10,5.6)); ax.axvspan(0,300,color="orange",alpha=0.06)
    for name,(rmse,_,col) in res.items(): ax.semilogy(tgrid,rmse,"-o",ms=3,color=col,label=name)
    ax.set_xlabel("time [s]"); ax.set_ylabel("全5点温度の推定RMSE [K]")
    ax.set_title("変位観測点は熱感度 $W=K^{-1}H$ の高い場所を選ぶべき\n"
                 "高熱感度の変位点ほど温度推定が速く収束（温度センサは固定, 5seed平均）")
    ax.grid(alpha=0.3,which="both"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"disp_sel_traj.png"),dpi=140); plt.close(fig)

    # 図2: per-node温度時刻歴（正解／低熱感度／高熱感度の3本に絞り、高熱感度を強調）
    show=[("＋変位2点(低熱感度W)","tab:orange",1.8,(2,2)),
          ("＋変位2点(高熱感度W)","tab:blue",2.4,(5,2))]
    fig,axes=plt.subplots(2,3,figsize=(15,8.2),sharex=True)
    for i in range(NPT):
        ax=axes[i//3][i%3]; ax.axvspan(0,300,color="orange",alpha=0.06)
        ax.plot(tgrid,Ttr[:,i]-K,"-",color="k",lw=4.5,alpha=0.35)
        for name,col,lw,dash in show:
            ax.plot(tgrid,res[name][1][:,i]-K,"--",color=col,lw=lw,dashes=dash)
        role=" (温度センサ=観測)" if i==Tsens else " (未観測)"
        ax.set_title(f"P{i}{role}",fontsize=12); ax.grid(alpha=0.3)
        if i%3==0: ax.set_ylabel("温度 [degC]")
        if i//3==1: ax.set_xlabel("time [s]")
    st=[Line2D([0],[0],color="k",lw=4.5,alpha=0.35,label="正解(真値)"),
        Line2D([0],[0],color="tab:orange",lw=1.8,ls="--",label="低熱感度の変位点で同化"),
        Line2D([0],[0],color="tab:blue",lw=2.4,ls="--",label="高熱感度の変位点で同化")]
    axes[1][2].axis("off"); axes[1][2].legend(handles=st,fontsize=13,loc="center")
    fig.suptitle("各ノード温度: 変位観測点を高熱感度(W)/低熱感度(W)で変えた同化 vs 正解\n"
                 "高熱感度(青)ほど過渡期に真値へ速く追従＝同化精度が高い",fontsize=14,weight="bold")
    fig.tight_layout(rect=[0,0,1,0.95]); fig.savefig(os.path.join(IMG,"disp_sel_pernode.png"),dpi=130); plt.close(fig)

    # 図2c: 温度だけ同化 vs 温度+変位同化 — 温度は両者ほぼ合うが、変位には明確な差が出るか
    nmT="温度1点のみ(変位なし)"; nmTD="＋変位2点(高熱感度W)"
    uz_true=dispA(Ttr,uz_hi,D_hi)                       # 真値温度場での高W点変位 [nt,2]
    uz_T  =dispA(res[nmT][1],uz_hi,D_hi)
    uz_TD =dispA(res[nmTD][1],uz_hi,D_hi)
    heat=(tgrid>0)&(tgrid<=300)
    eT_T  =np.abs(res[nmT][1]-Ttr)[heat].mean();  eT_TD=np.abs(res[nmTD][1]-Ttr)[heat].mean()
    eU_T  =np.abs(uz_T-uz_true)[heat].mean();     eU_TD=np.abs(uz_TD-uz_true)[heat].mean()
    print(f"[dsel] 加熱期平均誤差: 温度のみ  T={eT_T:.3f}K uz={eU_T:.3f}µm")
    print(f"[dsel] 加熱期平均誤差: 温度+変位 T={eT_TD:.3f}K uz={eU_TD:.3f}µm")
    print(f"[dsel] 誤差比(温度のみ/温度+変位): 温度 {eT_T/eT_TD:.1f}倍  変位 {eU_T/eU_TD:.1f}倍")
    fig,axes=plt.subplots(2,3,figsize=(15,8.4),sharex=True)
    for i in range(NPT):
        ax=axes[i//3][i%3]; ax.axvspan(0,300,color="orange",alpha=0.06)
        ax.plot(tgrid,Ttr[:,i]-K,"-",color="k",lw=4.5,alpha=0.35)
        ax.plot(tgrid,res[nmT][1][:,i]-K,"--",color="tab:green",lw=1.8,dashes=(2,2))
        ax.plot(tgrid,res[nmTD][1][:,i]-K,"--",color="tab:blue",lw=2.2,dashes=(5,2))
        role=" (温度センサ=観測)" if i==Tsens else " (未観測)"
        ax.set_title(f"P{i} 温度{role}",fontsize=12); ax.grid(alpha=0.3)
        if i%3==0: ax.set_ylabel("温度 [degC]")
        if i//3==1: ax.set_xlabel("time [s]")
    axd=axes[1][2]; axd.axvspan(0,300,color="orange",alpha=0.06)
    axd.plot(tgrid,uz_true[:,0],"-",color="k",lw=4.5,alpha=0.35)
    axd.plot(tgrid,uz_T[:,0],"--",color="tab:green",lw=1.8,dashes=(2,2))
    axd.plot(tgrid,uz_TD[:,0],"--",color="tab:blue",lw=2.2,dashes=(5,2))
    axd.set_title("変位 uz(高W点)  ←ここに差が出る",fontsize=12,color="darkred")
    axd.set_ylabel("uz [µm]"); axd.set_xlabel("time [s]"); axd.grid(alpha=0.3)
    st=[Line2D([0],[0],color="k",lw=4.5,alpha=0.35,label="正解(真値)"),
        Line2D([0],[0],color="tab:green",lw=1.8,ls="--",label=f"温度P{Tsens}の1点のみ (T誤差{eT_T:.2f}K, uz誤差{eU_T:.2f}µm)"),
        Line2D([0],[0],color="tab:blue",lw=2.2,ls="--",label=f"同じ温度P{Tsens}＋高W変位2点 (T誤差{eT_TD:.2f}K, uz誤差{eU_TD:.2f}µm)")]
    fig.legend(handles=st,fontsize=12,loc="lower center",ncol=3,bbox_to_anchor=(0.5,-0.015))
    fig.suptitle("温度だけ同化 vs 温度+変位同化: 温度(P0-P4)はどちらも概ね真値に合うが、変位には明確な差が残る\n"
                 "→ 変位は温度のわずかな誤差を増幅して映す(感度Wが大きい)＝変位まで合わせたいなら変位観測が必要",
                 fontsize=13.5,weight="bold")
    fig.tight_layout(rect=[0,0.05,1,0.93]); fig.savefig(os.path.join(IMG,"disp_sel_temp_vs_disp.png"),dpi=130); plt.close(fig)

    # 図2b: PODで選んだ5点"以外"の点でも比較（gappy-POD全場復元→非ROM点で読む）
    # 非ROM点 = 102_0の元プローブ(hot/cold/mid/top)＝Q-DEIMの5点とは別の場所
    probes={"hot(32,0,50)":[0.032,0,0.05025],"cold(-32,0,50)":[-0.032,0,0.05025],
            "mid(0,24,50)":[0,0.024,0.05025],"top(25,0,95)":[0.025,0,0.095]}
    pcell=[int(np.linalg.norm(Cc-np.array(p),axis=1).argmin()) for p in probes.values()]
    def field_at(T5,cellidx):     # gappy-POD復元して指定セルの温度
        a=(T5-mean[pod_cells])@UP_pinv.T; return (mean[cellidx]+U[cellidx,:]@a)
    fig,axes=plt.subplots(2,2,figsize=(13,8))
    for ax,(name,pc) in zip(axes.ravel(),zip(probes.keys(),pcell)):
        ax.axvspan(0,300,color="orange",alpha=0.06)
        tr=np.array([field_at(Ttr[k],pc) for k in range(len(tgrid))])-K
        lo=np.array([field_at(res["＋変位2点(低熱感度W)"][1][k],pc) for k in range(len(tgrid))])-K
        hi=np.array([field_at(res["＋変位2点(高熱感度W)"][1][k],pc) for k in range(len(tgrid))])-K
        ax.plot(tgrid,tr,"-",color="k",lw=4.5,alpha=0.35)
        ax.plot(tgrid,lo,":",color="tab:orange",lw=1.8)
        ax.plot(tgrid,hi,"--",color="tab:blue",lw=2.4,dashes=(5,2))
        ax.set_title(f"{name}  ※非ROM点(gappy-POD復元)",fontsize=12); ax.grid(alpha=0.3)
        ax.set_ylabel("温度 [degC]"); ax.set_xlabel("time [s]")
    st=[Line2D([0],[0],color="k",lw=4.5,alpha=0.35,label="正解(真値)"),
        Line2D([0],[0],color="tab:orange",lw=1.8,ls=":",label="低熱感度の変位点で同化"),
        Line2D([0],[0],color="tab:blue",lw=2.4,ls="--",label="高熱感度の変位点で同化")]
    fig.suptitle("POD5点以外の場所(非ROM点)の温度も、高熱感度の変位点の方が真値に近い\n"
                 "（gappy-PODで全場復元→元プローブ hot/cold/mid/top で評価）",fontsize=14,weight="bold")
    fig.legend(handles=st,fontsize=12,loc="lower center",ncol=3,bbox_to_anchor=(0.5,-0.02))
    fig.tight_layout(rect=[0,0.05,1,0.94]); fig.savefig(os.path.join(IMG,"disp_sel_nonrom.png"),dpi=130); plt.close(fig)

    # 図3: 過渡期RMSE棒グラフ
    heat=(tgrid>0)&(tgrid<=300); names=[c[0] for c in cfgs]
    vals=[res[n][0][heat].mean() for n in names]; cols=[c[2] for c in cfgs]
    fig,ax=plt.subplots(figsize=(9,5.2)); b=ax.bar(range(len(names)),vals,color=cols)
    for bar,v in zip(b,vals): ax.text(bar.get_x()+bar.get_width()/2,v,f"{v:.3f}",ha="center",va="bottom",fontsize=11)
    ax.set_yscale("log"); ax.set_xticks(range(len(names))); ax.set_xticklabels(names,fontsize=11)
    ax.set_ylabel("加熱期(0-300s)平均RMSE [K]"); ax.grid(alpha=0.3,axis="y",which="both")
    ax.set_title("変位観測点は熱感度(W)の高い場所を選ぶべき")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"disp_sel_rmse.png"),dpi=140); plt.close(fig)
    yaml.safe_dump({n:{"transient_rmse":float(res[n][0][heat].mean())} for n in names},
                   open(os.path.join(RES,"disp_selection.yaml"),"w"),allow_unicode=True)
    print("[dsel] wrote disp_sel_traj/pernode/rmse.png")


if __name__=="__main__": main()
