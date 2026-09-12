"""POD選定5点ROMで、観測構成を変えてデータ同化の有用性を比較する（本研究の核）.

比較する構成:
  (0) 同化なし（free run: でたらめ初期のまま前進、補正しない）
  (1) 温度1点（低感度 dT/dQ 最小ノード）
  (2) 温度1点（高感度 dT/dQ 最大ノード）
  (3) 温度2点（高感度上位2ノード）
  (4) 温度2点 + 変位2点（上面Uz、FrontISTRのPODモード応答で観測化）

指標: 全5点温度の推定RMSEの時刻歴と最終値、発熱量Qの誤差。
変位観測演算子: 全温度場(gappy-POD復元) → FrontISTR → 上面Uz2点。
  uz = uz_mean + Dmode·a （a=PODモード振幅）。Dmode は各モードのFrontISTR応答。

出力(docs/img/): da_compare_rmse.png, da_compare_traj.png
出力: results/da_compare.yaml
再現: OMP_NUM_THREADS=4 python3 run/run_da_compare.py
"""
from __future__ import annotations
import os, sys
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
K=273.15; NPT=5; IQ=NPT; IH=NPT+1; NAUG=NPT+2
DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; SIG_U=0.3; INFL=1.02
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])


def build_disp_operator(U, mean, Cc):
    """gappy-POD場→FrontISTR→上面Uz2点 の演算子を作る.
    戻り: uz_mean(2,), Dmode(2,r)  で uz = uz_mean + Dmode @ a."""
    cache=os.path.join(RES,"disp_operator.npz")
    if os.path.exists(cache):
        d=np.load(cache); return d["uz_mean"], d["Dmode"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    qa=int(np.linalg.norm(coords-A_XYZ,axis=1).argmin())
    qo=int(np.linalg.norm(coords-O_XYZ,axis=1).argmin())
    _,near=cKDTree(Cc).query(coords)          # FEM節点→最寄りOpenFOAMセル
    work=os.path.join(ROOT,"openfoam","dispop"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    def run_field(Tfield):   # FEM節点温度場 -> uz(qa),uz(qo)
        fistr_case.write_cnt(Path(work),node_ids,Tfield,reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work)); disp=fistr_case.read_displacement(Path(work))
        return np.array([disp[node_ids[qa]][2], disp[node_ids[qo]][2]])*1e6  # µm
    r=U.shape[1]
    uz_mean=run_field(mean[near])                       # 平均場の変位
    Dmode=np.zeros((2,r))
    for k in range(r):
        Dmode[:,k]=run_field(Tref+U[near,k])            # モードkの変位応答
        print(f"[disp] mode {k+1}/{r} FrontISTR",flush=True)
    np.savez(cache,uz_mean=uz_mean,Dmode=Dmode,qa=qa,qo=qo)
    return uz_mean, Dmode


def main():
    import yaml
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    cells=d["cells"]; xyz=d["xyz"]; heat_node=int(d["heat_node"]); C=d["C"]
    Kmat=rg.tri_to_matrix(d["K_upper"],NPT); h_true=float(d["h"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    pod_cells=kv["cell_idx"]
    UP=U[pod_cells,:]; UP_pinv=np.linalg.pinv(UP)      # a = UP_pinv (T5 - mean[pod])

    # dT/dQ でノード感度
    base=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.0,heat_node,0,300,DT)[1][-1]
    pert=rg.integrate_single(np.full(NPT,rg.T_AIR_K),C,Kmat,h_true,1.1,heat_node,0,300,DT)[1][-1]
    sens=(pert-base)/0.1
    hi=int(np.argmax(sens)); lo=int(np.argmin(sens))
    hi2=list(np.argsort(sens)[::-1][:2])
    print(f"[cmp] dT/dQ={np.round(sens,2)}  高感度={hi} 低感度={lo} 上位2={hi2}")

    uz_mean, Dmode = build_disp_operator(U, mean, Cc)
    def disp_of_T5(T5):   # (…,5)->(…,2)
        a=(T5-mean[pod_cells])@UP_pinv.T
        return uz_mean + a@Dmode.T

    # 真値トラジェクトリ
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tgrid=np.r_[0,cyc]
    Ttr=[np.full(NPT,rg.T_AIR_K)]; T=np.full(NPT,rg.T_AIR_K)
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Kmat,h_true,1.0,heat_node,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run_cfg_seed(nodes, use_disp, seed):
        rng=np.random.default_rng(seed); rng_o=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        if nodes is not None:
            Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
        rec=[np.sqrt(((Z[:,:NPT].mean(0)-Ttr[0])**2).mean())]; qrec=[Z[:,IQ].mean()]
        urec=[disp_of_T5(Z[:,:NPT].mean(0))]
        hrec=[Z[:,IH].mean()]
        tp=0.0
        for ci,tb in enumerate(cyc,1):
            Tn=rg.integrate_ensemble(Z[:,:NPT],C,Kmat,Z[:,IH],Z[:,IQ],heat_node,tp,tb,DT)
            Z=Z.copy(); Z[:,:NPT]=Tn; tp=tb
            if nodes is not None:   # 同化する
                yv=list(Ttr[ci][nodes]); Yf=Z[:,nodes]
                if use_disp:
                    yv=yv+list(disp_of_T5(Ttr[ci])); Yf=np.column_stack([Yf,disp_of_T5(Z[:,:NPT])])
                y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            rec.append(np.sqrt(((Z[:,:NPT].mean(0)-Ttr[ci])**2).mean())); qrec.append(Z[:,IQ].mean())
            urec.append(disp_of_T5(Z[:,:NPT].mean(0)))
            hrec.append(Z[:,IH].mean())
        return np.array(rec), np.array(qrec), np.array(urec), np.array(hrec)

    SEEDS=[20260913,20260914,20260915,20260916,20260917]
    def run_cfg(nodes, use_disp):
        rs=[run_cfg_seed(nodes,use_disp,s) for s in SEEDS]
        return (np.mean([r[0] for r in rs],axis=0), np.mean([r[1] for r in rs],axis=0),
                np.array([r[2] for r in rs]), np.array([r[1] for r in rs])*15,
                np.array([r[3] for r in rs]))

    configs=[
        ("同化なし(free run)",   None,      False, "tab:gray"),
        ("温度1点(低感度)",      [lo],      False, "tab:orange"),
        ("温度1点(高感度)",      [hi],      False, "tab:blue"),
        ("温度2点(高感度)",      hi2,       False, "tab:green"),
        ("温度2点+変位2点",      hi2,       True,  "tab:red"),
    ]
    results={}
    disp_results={}
    parameter_results={}
    for name,nodes,ud,col in configs:
        rec,qrec,urec,qseeds,hseeds=run_cfg(nodes,ud)
        parameter_results[name]=(qseeds,hseeds)
        disp_results[name]=urec
        results[name]=(rec,qrec,col)
        print(f"[cmp] {name:18s} 最終温度RMSE={rec[-1]:.4f}K  最終Q={qrec[-1]*15:.2f}W")

    # 図1: 収束トラジェクトリ
    fig,ax=plt.subplots(figsize=(10,5.6)); ax.axvspan(0,300,color="orange",alpha=0.06)
    for name,(rec,qrec,col) in results.items():
        ax.semilogy(tgrid,rec,"-o",ms=3,color=col,label=name)
    ax.set_xlabel("time [s]"); ax.set_ylabel("全5点温度の推定RMSE [K]")
    ax.set_title("観測構成によるデータ同化精度の違い（POD選定5点ROM, 60メンバー）")
    ax.grid(alpha=0.3,which="both"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_compare_traj.png"),dpi=140); plt.close(fig)

    # 図2: 過渡期(加熱0-300s)の平均RMSE棒グラフ（最終値は皆収束するので過渡で比較）
    heat=(tgrid>0)&(tgrid<=300)
    names=[c[0] for c in configs]; vals=[results[n][0][heat].mean() for n in names]; cols=[c[3] for c in configs]
    fig,ax=plt.subplots(figsize=(10,5.4))
    b=ax.bar(range(len(names)),vals,color=cols)
    for bar,v in zip(b,vals): ax.text(bar.get_x()+bar.get_width()/2,v,f"{v:.3f}",ha="center",va="bottom",fontsize=11)
    ax.set_yscale("log"); ax.set_xticks(range(len(names))); ax.set_xticklabels(names,fontsize=11)
    ax.set_ylabel("加熱期(0-300s)の平均 全5点温度RMSE [K]"); ax.grid(alpha=0.3,axis="y",which="both")
    ax.set_title("本研究の有用性: 観測を選び・足すほど過渡期の推定が良くなる")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_compare_rmse.png"),dpi=140); plt.close(fig)

    # 全5構成の同じ同化実験から変位を評価。seed平均前に誤差を計算する。
    truth_u=disp_of_T5(Ttr)
    all_u=np.array([disp_results[n] for n in names])  # config, seed, time, A/O
    np.savez(os.path.join(RES,"da_compare_displacement.npz"),time=tgrid,
             names=np.array(names),seeds=SEEDS,truth_u_um=truth_u,estimate_u_um=all_u)
    labels=["Uz(A)：ヒータ側上面", "Uz(O)：反ヒータ側上面", "Uz(A) − Uz(O)：上面2点の差"]
    truth3=np.column_stack([truth_u,truth_u[:,0]-truth_u[:,1]])
    est3=np.concatenate([all_u,(all_u[:,:,:,0]-all_u[:,:,:,1])[...,None]],axis=-1)
    fig,axes=plt.subplots(3,2,figsize=(15,12),sharex=True)
    for j in range(3):
        ax,er=axes[j]
        ax.plot(tgrid,truth3[:,j],color="black",lw=3,label="ROM真値→同じ変位写像")
        for i,name in enumerate(names):
            ax.plot(tgrid,est3[i,:,:,j].mean(0),color=cols[i],ls="--",label=name)
            er.plot(tgrid,np.abs(est3[i,:,:,j]-truth3[:,j]).mean(0),color=cols[i],label=name)
        ax.set_title(labels[j]+"（5seedの推定平均）")
        er.set_title("各seedの絶対誤差を計算 → 5seed平均")
        for a in (ax,er):
            a.axvspan(0,300,color="orange",alpha=0.08); a.grid(alpha=0.3)
            a.set_xlabel("時刻 [s]"); a.set_ylabel("変位 / 絶対誤差 [µm]")
    axes[0,0].legend(fontsize=8); axes[0,1].legend(fontsize=8)
    fig.suptitle("温度RMSEと同じ5構成の変位比較：60メンバー × 5seed\n"
                 "最後の構成のみA/O変位を同化に使用。参照は実測ではなくROM双子実験。")
    fig.tight_layout(rect=[0,0,1,0.95])
    fig.savefig(os.path.join(IMG,"da_compare_displacement.png"),dpi=140); plt.close(fig)
    # 各seedで時間RMSEを求めてからseed平均。温度指標との定義の差を明記。
    delta=est3[:,:,heat,:]-truth3[heat,:]
    scores=np.sqrt(np.mean(delta**2,axis=2))   # config, seed, component
    pair_scores=np.sqrt(np.mean(delta[:,:,:,:2]**2,axis=(2,3)))
    fig,axes=plt.subplots(1,3,figsize=(16,5))
    for j,ax in enumerate(axes):
        v=scores[:,:,j].mean(1)
        ax.bar(range(5),v,color=cols)
        for i,x in enumerate(v): ax.text(i,x,f"{x:.3f}",ha="center",va="bottom")
        ax.set_xticks(range(5)); ax.set_xticklabels(names,rotation=25,ha="right",fontsize=8)
        ax.set_title(labels[j]); ax.set_ylabel("加熱期の時間RMSE [µm]（5seed平均）"); ax.grid(axis="y",alpha=0.3)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_compare_displacement_rmse.png"),dpi=140); plt.close(fig)
    rows=[]
    for i,name in enumerate(names):
        rows.append({"configuration":name,"temperature_mean_spatial_rmse_K":float(vals[i]),
                     "A_time_rmse_um":float(scores[i,:,0].mean()),
                     "O_time_rmse_um":float(scores[i,:,1].mean()),
                     "difference_time_rmse_um":float(scores[i,:,2].mean()),
                     "pair_time_space_rmse_um":float(pair_scores[i].mean())})
    import csv
    with open(os.path.join(RES,"da_compare_displacement.csv"),"w") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    lines=["| 観測構成 | 温度RMSE [K] | Uz(A) RMSE [µm] | Uz(O) RMSE [µm] | 差 Uz(A)−Uz(O) RMSE [µm] |",
           "|---|---:|---:|---:|---:|"]
    for r in rows:
        lines.append("| "+r["configuration"]+" | "+" | ".join(f"{r[k]:.3f}" for k in
             ["temperature_mean_spatial_rmse_K","A_time_rmse_um","O_time_rmse_um","difference_time_rmse_um"])+" |")
    Path(RES,"da_compare_displacement_table.md").write_text("\n".join(lines)+"\n")
    print("\n".join(lines))

    all_q=np.array([parameter_results[n][0] for n in names])
    all_h=np.array([parameter_results[n][1] for n in names])
    np.savez(os.path.join(RES,"da_compare_parameters.npz"),time=tgrid,seeds=SEEDS,
             names=np.array(names),Q_W=all_q,h_W_K=all_h,Q_true_W=15.,h_true_W_K=h_true)
    fig,axes=plt.subplots(1,2,figsize=(14,5))
    for ax,arr,truth,label in zip(axes,[all_q,all_h],[15.,h_true],["発熱量Q [W]","放熱係数h [W/K]"]):
        for i,name in enumerate(names):
            ax.plot(tgrid,arr[i].mean(0),color=cols[i],label=name)
            ax.fill_between(tgrid,arr[i].min(0),arr[i].max(0),color=cols[i],alpha=0.10)
        ax.axhline(truth,color="black",ls="--",label="真値")
        ax.axvspan(0,300,color="orange",alpha=0.06)
        ax.set_xlabel("時刻 [s]"); ax.set_ylabel(label); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=8)
    fig.suptitle("Q・hの推定：線は5seed平均、帯は5seedの最小〜最大（信頼区間ではない）")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"da_compare_parameters.png"),dpi=140); plt.close(fig)
    prows=[]
    for i,name in enumerate(names):
        for k in [10,20]:
            q=all_q[i,:,k]; h=all_h[i,:,k]
            prows.append(dict(configuration=name,time_s=float(tgrid[k]),
                Q_mean_W=float(q.mean()),Q_min_W=float(q.min()),Q_max_W=float(q.max()),
                Q_mean_abs_relative_error_pct=float(np.abs(q/15-1).mean()*100),
                h_mean_W_K=float(h.mean()),h_min_W_K=float(h.min()),h_max_W_K=float(h.max()),
                h_mean_abs_relative_error_pct=float(np.abs(h/h_true-1).mean()*100)))
    with open(os.path.join(RES,"da_compare_parameters.csv"),"w") as f:
        writer=csv.DictWriter(f,fieldnames=list(prows[0])); writer.writeheader(); writer.writerows(prows)
    print("PARAMETERS",prows)

    yaml.safe_dump({n:{"final_rmse_K":float(results[n][0][-1]),
        "final_Q_W":float(results[n][1][-1]*15)} for n in names}|
        {"dTdQ":[float(x) for x in sens],"hi_node":hi,"lo_node":lo},
        open(os.path.join(RES,"da_compare.yaml"),"w"),allow_unicode=True)
    print("[cmp] wrote da_compare_rmse.png, da_compare_traj.png, da_compare.yaml")


if __name__=="__main__": main()
