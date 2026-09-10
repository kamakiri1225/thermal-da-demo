"""ROMで「温度のみ」vs「温度＋変位M」を比較（変位が効くことの実証）。

変位観測は、5点温度をIDWでメッシュ再構成→FrontISTRで上面変位、という物理チェーン。
FrontISTRは線形なので、5つの単位温度モードの応答を並べた行列 M（config/M_frontistr_operator.npy,
µm/K）で厳密に表せる。観測はアフィン uz=M(T-Tref) を各メンバーで評価して同化する。
出力: docs/img/rom_disp_comparison.png
"""
from __future__ import annotations
import os, sys
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore.calibrate import load_calibrated
from dacore.enkf import enkf_update
from dacore.ensemble import init_ensemble, clip_params, forecast, N_AUG
from dacore.observations import generate_truth, obs_node_indices
from dacore.cht_rom import N_NODES, T_AIR_K, NODE_NAMES
IMG=os.path.join(ROOT,"docs","img")
Mum=np.load(os.path.join(ROOT,"config","M_frontistr_operator.npy"))*1e6  # µm/K

def run(cfg,calib,tnodes,use_disp,sigU):
    sigT=cfg["observation"]["noise_C"]
    rng_o=np.random.default_rng(cfg["experiment"]["seed"]); rng=np.random.default_rng(cfg["experiment"]["seed"]+1)
    times,truth=generate_truth(calib,cfg); dt=cfg["experiment"]["obs_interval_s"]
    cyc=np.arange(dt,cfg["experiment"]["t_end_s"]+1e-9,dt); idx=np.clip(np.searchsorted(times,cyc),0,len(times)-1)
    nt=len(tnodes); no=nt+(2 if use_disp else 0)
    R=np.diag([sigT**2]*nt+([sigU**2]*2 if use_disp else []))
    Z=init_ensemble(cfg,rng); tp=0; rt=[0.0]; rr=[float(np.sqrt(((Z[:,:N_NODES].mean(0)-truth[0])**2).mean()))]
    for t1,ci in zip(cyc,idx):
        Z=forecast(Z,tp,t1,calib,cfg,rng); Tt=truth[ci]
        Yf=np.zeros((len(Z),no)); Yf[:,:nt]=Z[:,tnodes]
        if use_disp: Yf[:,nt:]=(Z[:,:N_NODES]-T_AIR_K)@Mum.T
        yv=list(Tt[tnodes])+(list(Mum@(Tt-T_AIR_K)) if use_disp else [])
        y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(R)))
        Z=enkf_update(Z,y,None,R,rng,inflation=cfg["filter"]["inflation"],Yf=Yf); clip_params(Z,cfg); tp=t1
        rt.append(float(t1)); rr.append(float(np.sqrt(((Z[:,:N_NODES].mean(0)-Tt)**2).mean())))
    return np.array(rt),np.array(rr)

def main():
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated(); hot=NODE_NAMES.index("hot"); cold=NODE_NAMES.index("cold")
    t0,r0=run(cfg,calib,[hot,cold],False,0)
    t1,r1=run(cfg,calib,[hot,cold],True,0.3)
    t2,r2=run(cfg,calib,[hot],False,0)
    t3,r3=run(cfg,calib,[hot],True,0.1)
    fig,ax=plt.subplots(figsize=(9.5,5.2)); ax.axvspan(0,300,color="orange",alpha=0.06)
    ax.plot(t0,r0,"-o",ms=3,color="tab:blue",label=f"温度2点のみ（最終 {r0[-1]:.3f} K）")
    ax.plot(t1,r1,"-s",ms=3,color="tab:green",label=f"温度2点＋変位M（最終 {r1[-1]:.3f} K）")
    ax.plot(t2,r2,"--o",ms=3,color="gray",label=f"温度1点のみ（最終 {r2[-1]:.3f} K）")
    ax.plot(t3,r3,"-^",ms=3,color="tab:red",label=f"温度1点＋変位M（最終 {r3[-1]:.3f} K）")
    ax.set_yscale("log"); ax.set_xlabel("time [s]"); ax.set_ylabel("温度RMSE(全5ノード) [K]")
    ax.set_title("ROM：変位観測（IDW→FrontISTR演算子M）を足すと改善する")
    ax.grid(alpha=0.3,which="both"); ax.legend(fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"rom_disp_comparison.png"),dpi=140)
    import subprocess; subprocess.run(["cp",os.path.join(IMG,"rom_disp_comparison.png"),
        os.path.join(ROOT,"presentation","assets","rom_disp_comparison.png")],check=False)
    print(f"温度2点:{r0[-1]:.3f} / +変位:{r1[-1]:.3f} / 温度1点:{r2[-1]:.3f} / 1点+変位:{r3[-1]:.3f} K")
    print("wrote docs/img/rom_disp_comparison.png")

if __name__=="__main__": main()
