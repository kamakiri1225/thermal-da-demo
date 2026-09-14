"""Independent support experiment: T2 + high-W Uz2, held-out C/D histories.
Run: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 python3 run/run_holdout_displacement_support.py
"""
from pathlib import Path
import json, hashlib, csv
import numpy as np
import run_disp_selection as b

ROOT = Path(__file__).absolute().parent.parent
# Preserve writable spelling on case-insensitive mounted drives.
ROOT = Path(str(ROOT).replace("/002_CAE/", "/002_cae/"))
b.ROOT = str(ROOT); b.RES = str(ROOT / "results")
res = ROOT / "results"; img = ROOT / "docs/img"

def main():
    cal = np.load(res / "rom_calibrated_pod.npz")
    pod = np.load(res / "qdeim_points.npz")
    U, mean, cells = pod["pod_modes"].astype(float), pod["mean"], pod["cell_idx"]
    C, km, ht, hn = cal["C"], b.rg.tri_to_matrix(cal["K_upper"], 5), float(cal["h"]), int(cal["heat_node"])
    mesh = b.cylinder_mesh.build_cylinder_mesh(b.NR,b.NTH,b.NZ,b.R_IN,b.R_OUT,b.H)
    xyz = np.array([x for _,x in mesh["nodes"]]); ids = np.array([i for i,_ in mesh["nodes"]])
    high = np.load(b.KINVH)["hi"].astype(int)
    validation = np.array([np.linalg.norm(xyz-[x,0,b.H],axis=1).argmin() for x in [.028,-.028]])
    assert not set(high) & set(validation)
    nodes = np.r_[high,validation]
    digest = hashlib.sha256()
    for arr in [U,mean,pod["cell_centres"],xyz,nodes]: digest.update(np.ascontiguousarray(arr).tobytes())
    digest.update(json.dumps(b.MATERIAL,sort_keys=True).encode())
    tag = "holdout_support_" + digest.hexdigest()[:16]
    u0,D = b.build_disp_op(U,mean,pod["cell_centres"],nodes,tag)
    pinv = np.linalg.pinv(U[cells])
    def displacement(T): return u0 + ((T-mean[cells]) @ pinv.T) @ D.T
    t = np.arange(0,601,30)
    truth = [np.full(5,b.rg.T_AIR_K)]
    for ta,tb in zip(t[:-1],t[1:]):
        truth.append(b.rg.integrate_single(truth[-1],C,km,ht,1.,hn,ta,tb,2.)[1][-1])
    truth = np.array(truth); true_u = displacement(truth)
    temps = [2,0]; names = ["No DA","T2 only","T2 + high-W Uz2"]
    states, us = [], []
    for case in range(3):
        zs, uu = [], []
        for seed in b.SEEDS:
            initial = np.random.default_rng(seed)
            Z = np.zeros((60,7)); Z[:,:5] = initial.uniform(b.rg.T_AIR_K-3,b.rg.T_AIR_K+12,(60,5))
            Z[:,5] = initial.uniform(.3,1.8,60); Z[:,6] = np.clip(initial.normal(.02,.01,60),.001,.1)
            noiseT = np.random.default_rng(seed+100).normal(0,.3,(20,2))
            noiseU = np.random.default_rng(seed+200).normal(0,.3,(20,2))
            update_rng = np.random.default_rng(seed+300)
            zh = [Z.mean(0)]; uh = [displacement(Z[:,:5]).mean(0)]
            for k,(ta,tb) in enumerate(zip(t[:-1],t[1:])):
                Z[:,:5] = b.rg.integrate_ensemble(Z[:,:5],C,km,Z[:,6],Z[:,5],hn,ta,tb,2.)
                if case:
                    yf = Z[:,temps]; obs = truth[k+1,temps]+noiseT[k]
                    if case == 2:
                        yf = np.column_stack((yf,displacement(Z[:,:5])[:,:2]))
                        obs = np.r_[obs,true_u[k+1,:2]+noiseU[k]]
                    Z = b.enkf_update(Z,obs,None,np.eye(len(obs))*.3**2,update_rng,inflation=1.02,Yf=yf)
                    Z[:,5] = np.clip(Z[:,5],0,3); Z[:,6] = np.clip(Z[:,6],.0001,.2)
                zh.append(Z.mean(0)); uh.append(displacement(Z[:,:5]).mean(0))
            zs.append(zh); uu.append(uh)
        states.append(zs); us.append(uu)
    states, us = np.array(states),np.array(us)
    assert np.isfinite(states).all() and np.isfinite(us).all()
    assert np.allclose(states[:, :, 0], states[0, :, 0])
    diff = us[:,:,:,2]-us[:,:,:,3]; true_diff = true_u[:,2]-true_u[:,3]
    np.savez_compressed(res/"holdout_support_history.npz",time=t,states=states,u_um=us,truth_T=truth,truth_u_um=true_u,names=names,seeds=b.SEEDS,nodes=nodes,node_ids=ids[nodes],xyz_m=xyz[nodes],temperature_indices=temps)
    meta = dict(names=names,seeds=b.SEEDS,temperature_indices=temps,temperature_xyz_mm=(pod["cell_centres"][cells[temps]]*1000).tolist(),node_indices=nodes.tolist(),node_ids=ids[nodes].tolist(),xyz_mm=(xyz[nodes]*1000).tolist(),operator_tag=tag,ensemble=60,dt_s=2,obs_dt_s=30,sigma_T_K=.3,sigma_u_um=.3,inflation=1.02,h_true_W_per_K=ht,Q_true_W=15,truth="same ROM and same displacement operator; internal validation")
    (res/"holdout_support_conditions.json").write_text(json.dumps(meta,ensure_ascii=False,indent=2)+"\n")
    rows=[]
    for ci,name in enumerate(names):
        for phase,mask in [("heating",(t>0)&(t<=300)),("cooling",t>300)]:
            for label,pred,truthv in [("C",us[ci,:,:,2],true_u[:,2]),("D",us[ci,:,:,3],true_u[:,3]),("C-D",diff[ci],true_diff)]:
                perseed=np.sqrt(np.mean((pred[:,mask]-truthv[mask])**2,axis=1))
                rows.append([name,phase,label,float(perseed.mean()),float(perseed.std(ddof=1))])
    with (res/"holdout_support_rmse.csv").open("w") as f:
        w=csv.writer(f);w.writerow(["case","phase","target","mean_seed_RMSE_um","std_seed_RMSE_um"]);w.writerows(rows)
    colors=["tab:gray","tab:green","tab:blue"]
    fig,axs=b.plt.subplots(2,2,figsize=(13,9))
    for ax,label,values,tv in [(axs[0,0],"C: unobserved Uz",us[:,:,:,2],true_u[:,2]),(axs[0,1],"D: unobserved Uz",us[:,:,:,3],true_u[:,3]),(axs[1,0],"C - D: unobserved displacement difference",diff,true_diff)]:
        ax.plot(t,tv,color="black",lw=2.5,label="Truth (same ROM)")
        for i,name in enumerate(names): ax.plot(t,values[i].mean(0),"--o",ms=3,color=colors[i],label=name)
        ax.axvspan(0,300,color="orange",alpha=.07); ax.axvline(300,color="gray",ls=":")
        ax.set(title=label,xlabel="Time [s]",ylabel="Displacement [µm]");ax.grid(alpha=.25)
    ax=axs[1,1]
    for i,name in enumerate(names):
        err=diff[i]-true_diff
        ax.plot(t,err.mean(0),color=colors[i],label=name)
        ax.fill_between(t,err.min(0),err.max(0),color=colors[i],alpha=.12)
    ax.axhline(0,color="black",lw=1);ax.axvline(300,color="gray",ls=":")
    ax.set(title="C-D error: mean and range across 5 seeds",xlabel="Time [s]",ylabel="Estimate - truth [µm]");ax.grid(alpha=.25)
    handles,labels=axs[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc="lower center",ncol=4)
    fig.suptitle("106 ROM holdout validation | T: P2/P0; observed Uz: high-W S1/S2; held out: C/D\n60 members, analysis every 30 s, 5-seed means; heating 0-300 s, cooling 300-600 s",fontsize=12)
    fig.tight_layout(rect=[0,.05,1,.93]);fig.savefig(img/"holdout_support_displacement.png",dpi=160);b.plt.close(fig)
    fig,ax=b.plt.subplots(figsize=(8,7))
    theta=np.linspace(0,2*np.pi,300)
    for radius in [20,37.5]:ax.plot(radius*np.cos(theta),radius*np.sin(theta),color="gray")
    for j,label in enumerate(["S1 observed","S2 observed","C held out","D held out"]):
        xy=xyz[nodes[j],:2]*1000;ax.scatter(*xy,c="tab:blue" if j<2 else "tab:red",marker="o" if j<2 else "s",s=70);ax.annotate(label,xy,xytext=(5,8),textcoords="offset points")
    ax.set(xlabel="x [mm]",ylabel="y [mm]",title="Top surface z = 100.5 mm: displacement locations",aspect="equal",xlim=(-49,49),ylim=(-45,45));ax.grid(alpha=.2)
    fig.tight_layout();fig.savefig(img/"holdout_support_locations.png",dpi=140);b.plt.close(fig)
    print(json.dumps(rows,ensure_ascii=False,indent=2))

if __name__ == "__main__": main()
