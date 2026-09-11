"""温度センサをどこに置くとデータ同化の精度が出るか ― 感度分布と突き合わせる実験.

仮説: 「変位観測の熱感度が高い場所に温度センサを置くほど、(温度+変位2点の)
データ同化の精度が出やすい」。これを実際に計算して検証する。

実験:
  温度センサ1点を 5候補(hot/mid/cold/top/core) のそれぞれに置き、
  (a) 温度1点のみ同化   (b) 温度1点 + 変位2点(M演算子, 104の手法) を
  複数seedで回して最終温度RMSEを比較。
  各候補位置の「熱感度」は run/sensitivity_map.py の FrontISTR 感度分布
  S(x)=∂uz/∂T(x) から読む。

出力:
  results/sensor_placement.csv
  docs/img/sensor_placement_rmse.png (棒グラフ+感度注記)
  docs/img/sensitivity_vs_rmse.png  (感度 vs 精度の散布図)
"""
from __future__ import annotations
import os, sys
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
from dacore import plots as _p  # フォント
import matplotlib.pyplot as plt
from dacore.calibrate import load_calibrated
from dacore.enkf import enkf_update
from dacore.ensemble import init_ensemble, clip_params, forecast, N_AUG
from dacore.observations import generate_truth
from dacore.cht_rom import N_NODES, T_AIR_K, NODE_NAMES
from dacore.node_locations import NODE_XYZ, HEIGHT

IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
Mum=np.load(os.path.join(ROOT,"config","M_frontistr_operator.npy"))*1e6  # µm/K
SIG_U=0.3   # 変位ノイズ [µm]
SEEDS=[20260908,20260909,20260910,20260911,20260912]


def run_da(cfg,calib,sensor_idx,use_disp,seed):
    sigT=cfg["observation"]["noise_C"]
    rng_o=np.random.default_rng(seed); rng=np.random.default_rng(seed+1)
    times,truth=generate_truth(calib,cfg); dt=cfg["experiment"]["obs_interval_s"]
    cyc=np.arange(dt,cfg["experiment"]["t_end_s"]+1e-9,dt)
    idx=np.clip(np.searchsorted(times,cyc),0,len(times)-1)
    no=1+(2 if use_disp else 0)
    R=np.diag([sigT**2]+([SIG_U**2]*2 if use_disp else []))
    cfg2=dict(cfg); cfg2["experiment"]=dict(cfg["experiment"]); cfg2["experiment"]["seed"]=seed
    Z=init_ensemble(cfg2,rng); tp=0
    for t1,ci in zip(cyc,idx):
        Z=forecast(Z,tp,t1,calib,cfg,rng); Tt=truth[ci]
        Yf=np.zeros((len(Z),no)); Yf[:,0]=Z[:,sensor_idx]
        if use_disp: Yf[:,1:]=(Z[:,:N_NODES]-T_AIR_K)@Mum.T
        yv=[Tt[sensor_idx]]+(list(Mum@(Tt-T_AIR_K)) if use_disp else [])
        y=np.array(yv)+rng_o.normal(0,np.sqrt(np.diag(R)))
        Z=enkf_update(Z,y,None,R,rng,inflation=cfg["filter"]["inflation"],Yf=Yf)
        clip_params(Z,cfg); tp=t1
    m=Z[:,:N_NODES].mean(0)
    return float(np.sqrt(((m-truth[idx[-1]])**2).mean()))


def sensitivity_at_nodes():
    """FrontISTR感度分布から、各ROMノード位置の |∂uz/∂T| (2観測平均) を読む [µm/K]."""
    d=np.load(os.path.join(RES,"sensitivity_uz.npz"))
    S=d["S_um"]; NTH_P=int(d["NTH_P"]); NZ_P=int(d["NZ_P"])
    out={}
    for name,(x,y,z) in NODE_XYZ.items():
        th=np.arctan2(y,x)%(2*np.pi)
        a=min(int(th/(2*np.pi)*NTH_P),NTH_P-1); b=min(int(z/HEIGHT*NZ_P),NZ_P-1)
        out[name]=float(np.abs(S[:,a,b]).mean())
    return out


def main():
    os.makedirs(IMG,exist_ok=True); os.makedirs(RES,exist_ok=True)
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated()
    sens=sensitivity_at_nodes()
    print("[study] 各ノード位置の変位熱感度 [µm/K]:",
          {k:round(v,4) for k,v in sens.items()})

    rows=[]
    for i,name in enumerate(NODE_NAMES):
        r_t=[run_da(cfg,calib,i,False,s) for s in SEEDS]
        r_td=[run_da(cfg,calib,i,True,s) for s in SEEDS]
        rows.append((name,sens[name],np.mean(r_t),np.std(r_t),np.mean(r_td),np.std(r_td)))
        print(f"[study] sensor={name:5s} 感度={sens[name]:.4f}µm/K  "
              f"温度のみ {np.mean(r_t):.3f}±{np.std(r_t):.3f}K  "
              f"+変位2点 {np.mean(r_td):.4f}±{np.std(r_td):.4f}K",flush=True)

    with open(os.path.join(RES,"sensor_placement.csv"),"w") as f:
        f.write("sensor,sens_um_per_K,rmse_temponly_K,std_temponly,rmse_tempdisp_K,std_tempdisp\n")
        for r in rows: f.write(",".join(str(x) for x in r)+"\n")

    names=[r[0] for r in rows]; sv=[r[1] for r in rows]
    rt=[r[2] for r in rows]; rts=[r[3] for r in rows]
    rd=[r[4] for r in rows]; rds=[r[5] for r in rows]

    # 棒グラフ
    x=np.arange(len(names)); w=0.38
    fig,ax=plt.subplots(figsize=(10.5,5.6))
    ax.bar(x-w/2,rt,w,yerr=rts,capsize=4,color="tab:gray",label="温度1点のみ同化")
    ax.bar(x+w/2,rd,w,yerr=rds,capsize=4,color="tab:blue",label="温度1点＋変位2点同化(104の手法)")
    ax.set_yscale("log"); ax.set_xticks(x)
    ax.set_xticklabels([f"{n}\n感度{sens[n]:.3f}µm/K" for n in names],fontsize=13)
    ax.set_ylabel("最終温度RMSE(全5ノード) [K]"); ax.grid(alpha=0.3,axis="y",which="both")
    ax.set_title("温度センサの設置位置とデータ同化精度（5seed平均±σ）")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"sensor_placement_rmse.png"),dpi=140); plt.close(fig)

    # 散布図: 感度 vs RMSE
    fig,ax=plt.subplots(figsize=(8.5,5.4))
    ax.scatter(sv,rt,s=120,color="tab:gray",label="温度1点のみ")
    ax.scatter(sv,rd,s=120,color="tab:blue",marker="s",label="温度1点＋変位2点")
    for n,xx,y1,y2 in zip(names,sv,rt,rd):
        ax.annotate(n,(xx,y1),textcoords="offset points",xytext=(6,6),fontsize=13)
        ax.annotate(n,(xx,y2),textcoords="offset points",xytext=(6,-14),fontsize=13,color="tab:blue")
    ax.set_yscale("log"); ax.set_xlabel("センサ位置の変位熱感度 |∂uz/∂T| [µm/K] (FrontISTR実計算)")
    ax.set_ylabel("最終温度RMSE [K]"); ax.grid(alpha=0.3,which="both")
    ax.set_title("熱感度が高い場所ほどデータ同化の精度が出るか？")
    ax.legend(); fig.tight_layout()
    fig.savefig(os.path.join(IMG,"sensitivity_vs_rmse.png"),dpi=140); plt.close(fig)
    print("[study] wrote figs + results/sensor_placement.csv")

if __name__=="__main__": main()
