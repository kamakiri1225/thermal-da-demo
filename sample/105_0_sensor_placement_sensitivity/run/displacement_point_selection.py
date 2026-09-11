"""変位観測点そのものを熱感度で選ぶと同化精度が上がるか（ユーザー仮説の本命検証）.

別セッション(FrontISTR W=K^-1 H, 011_Tji_DUMPW)で「高感度変位2点でQoI誤差89%→8%」
が示された。同じ検証を本ケース(中空円筒×ROM×M演算子同化)で行う。

方法:
 1. 5つの単位温度モードをFrontISTRに通し、全5040節点のuz応答 M_all (n_nodes×5) を取得
    (5回の実行で全点の感度が得られる)
 2. 感度指標 s_i = ||M_all[i,:]|| で節点をランク付け。底面固定近傍(z<5mm)は除外
    (別セッションの教訓: 固定節点の感度はアーティファクト)
 3. 変位観測2点を {高感度2点 / 低感度2点 / 現行(上面ヒータ側+反対側) / ランダム2点}
    で選び、ROMデータ同化(温度hot1点+変位2点, 5seed)の精度を比較
出力:
  docs/img/disp_point_sensitivity_3d.png … 全節点の感度と選定点
  docs/img/disp_point_selection_rmse.png … 選び方별の精度
  results/disp_point_selection.csv
"""
from __future__ import annotations
import os, sys
from pathlib import Path
import numpy as np, yaml
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); sys.path.insert(0,ROOT)
sys.path.insert(0, os.path.abspath(os.path.join(
    ROOT,"..","102_1_frontistr_hollow_cylinder_thermal_expansion","python")))
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore.calibrate import load_calibrated
from dacore.enkf import enkf_update
from dacore.ensemble import init_ensemble, clip_params, forecast, N_AUG
from dacore.observations import generate_truth
from dacore.cht_rom import N_NODES, T_AIR_K, NODE_NAMES
from dacore.node_locations import NODE_XYZ, HEIGHT

NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
IMG=os.path.join(ROOT,"docs","img"); RES=os.path.join(ROOT,"results")
SIG_U=0.3; SEEDS=[20260908,20260909,20260910,20260911,20260912]
HOT=NODE_NAMES.index("hot")


def build_M_all():
    """5単位モード×FrontISTR5回で、全節点のuz応答 M_all[n_nodes,5] [µm/K] を得る."""
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([xyz for _n,xyz in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    romn=np.array(list(NODE_XYZ.values()))
    W=np.zeros((len(coords),5))
    for i,p in enumerate(coords):
        dd=np.linalg.norm(romn-p,axis=1)
        if dd.min()<1e-9: W[i,dd.argmin()]=1
        else: w=1/dd**2; W[i]=w/w.sum()
    work=os.path.join(ROOT,"openfoam","Mall_probe"); os.makedirs(work,exist_ok=True)
    fistr_case.write_mesh(Path(work),NR,NTH,NZ,R_IN,R_OUT,H,
        young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(Path(work))
    M=np.zeros((len(coords),5))
    for k in range(5):
        fistr_case.write_cnt(Path(work),node_ids,Tref+W@np.eye(5)[k],reference_temperature=Tref,
            young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
            thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(Path(work))
        disp=fistr_case.read_displacement(Path(work))
        M[:,k]=[disp[n][2] for n in node_ids]
        print(f"[Mall] mode {k+1}/5",flush=True)
    return coords, M*1e6   # µm/K


def run_da(cfg,calib,Mrows_um,seed):
    """観測= hot温度1点 + 指定2行の変位。最終温度RMSEを返す."""
    sigT=cfg["observation"]["noise_C"]
    rng_o=np.random.default_rng(seed); rng=np.random.default_rng(seed+1)
    times,truth=generate_truth(calib,cfg); dt=cfg["experiment"]["obs_interval_s"]
    cyc=np.arange(dt,cfg["experiment"]["t_end_s"]+1e-9,dt)
    idx=np.clip(np.searchsorted(times,cyc),0,len(times)-1)
    R=np.diag([sigT**2,SIG_U**2,SIG_U**2])
    cfg2=dict(cfg); cfg2["experiment"]=dict(cfg["experiment"]); cfg2["experiment"]["seed"]=seed
    Z=init_ensemble(cfg2,rng); tp=0
    for t1,ci in zip(cyc,idx):
        Z=forecast(Z,tp,t1,calib,cfg,rng); Tt=truth[ci]
        Yf=np.zeros((len(Z),3)); Yf[:,0]=Z[:,HOT]
        Yf[:,1:]=(Z[:,:N_NODES]-T_AIR_K)@Mrows_um.T
        y=np.concatenate([[Tt[HOT]],Mrows_um@(Tt-T_AIR_K)])+rng_o.normal(0,np.sqrt(np.diag(R)))
        Z=enkf_update(Z,y,None,R,rng,inflation=cfg["filter"]["inflation"],Yf=Yf)
        clip_params(Z,cfg); tp=t1
    m=Z[:,:N_NODES].mean(0)
    return float(np.sqrt(((m-truth[idx[-1]])**2).mean()))


def main():
    os.makedirs(IMG,exist_ok=True); os.makedirs(RES,exist_ok=True)
    cfg=yaml.safe_load(open(os.path.join(ROOT,"config","da_config.yaml")))
    calib=load_calibrated()
    coords,M=build_M_all()
    sens=np.abs(M).sum(axis=1)                      # 感度指標 [µm/K]
    free=coords[:,2]>0.005                          # 底面固定近傍を除外(教訓)
    order=np.argsort(sens)[::-1]
    order_free=[i for i in order if free[i]]
    hi=[order_free[0]]
    # 2点目は1点目と方向(モード応答)が異なる高感度点(冗長回避, 別セッションの知見)
    v0=M[hi[0]]/np.linalg.norm(M[hi[0]])
    for i in order_free[1:]:
        v=M[i]/max(np.linalg.norm(M[i]),1e-12)
        if abs(v@v0)<0.9: hi.append(i); break
    lo=[i for i in order[::-1] if free[i] and sens[i]>1e-4][:2]
    cur=[int(np.linalg.norm(coords-np.array(p),axis=1).argmin())
         for p in [[0.028,0,H],[-0.028,0,H]]]        # 現行(104)の2点
    rng=np.random.default_rng(0)
    rnd_sets=[list(rng.choice(np.where(free)[0],2,replace=False)) for _ in range(3)]

    cases={"高感度2点(方向分散)":hi,"現行(上面ヒータ側+反対側)":cur,"低感度2点":lo}
    rows=[]
    for name,pts in cases.items():
        r=[run_da(cfg,calib,M[pts],s) for s in SEEDS]
        rows.append((name,float(np.mean([sens[p] for p in pts])),np.mean(r),np.std(r)))
        print(f"[sel] {name:22s} 感度={rows[-1][1]:.3f}µm/K RMSE={np.mean(r):.4f}±{np.std(r):.4f}K",flush=True)
    rr=[np.mean([run_da(cfg,calib,M[p],s) for s in SEEDS]) for p in rnd_sets]
    rows.append(("ランダム2点(3組平均)",float(np.mean([sens[p].mean() for p in rnd_sets])),
                 float(np.mean(rr)),float(np.std(rr))))
    print(f"[sel] ランダム2点             RMSE={np.mean(rr):.4f}±{np.std(rr):.4f}K",flush=True)

    with open(os.path.join(RES,"disp_point_selection.csv"),"w") as f:
        f.write("case,sens_um_per_K,rmse_K,std_K\n")
        for r in rows: f.write(",".join(str(x) for x in r)+"\n")

    # 3D: 感度分布と選定点
    import pyvista as pv
    pv.OFF_SCREEN=True
    try: pv.start_xvfb()
    except Exception: pass
    pl=pv.Plotter(off_screen=True,window_size=(1000,820))
    cloud=pv.PolyData(coords); cloud["sens [um/K]"]=sens
    pl.add_mesh(cloud,scalars="sens [um/K]",cmap="viridis",point_size=6,
                render_points_as_spheres=True,
                scalar_bar_args={"title":"|dUz/dmode| sum [um/K]","title_font_size":18,"label_font_size":14})
    for p,c,lab in [(hi,"red","high"),(lo,"blue","low"),(cur,"orange","current")]:
        for i in p: pl.add_mesh(pv.Sphere(radius=0.0035,center=coords[i]),color=c)
    pl.add_text("displacement-point sensitivity\nred=high2  orange=current  blue=low2",
                font_size=16,color="black")
    pl.camera_position=[(0.24,-0.22,0.20),(0,0,0.05),(0,0,1)]; pl.set_background("white")
    pl.screenshot(os.path.join(IMG,"disp_point_sensitivity_3d.png")); pl.close()

    # 棒グラフ
    names=[r[0] for r in rows]; vals=[r[2] for r in rows]; errs=[r[3] for r in rows]
    fig,ax=plt.subplots(figsize=(10.5,5.4))
    ax.bar(range(len(rows)),vals,yerr=errs,capsize=4,
           color=["tab:red","tab:orange","tab:blue","tab:gray"])
    ax.set_xticks(range(len(rows))); ax.set_xticklabels(names,fontsize=12)
    ax.set_ylabel("最終温度RMSE [K]"); ax.set_yscale("log")
    ax.set_title("変位観測2点の選び方とデータ同化精度（温度hot1点+変位2点, 5seed）")
    ax.grid(alpha=0.3,axis="y",which="both")
    fig.tight_layout(); fig.savefig(os.path.join(IMG,"disp_point_selection_rmse.png"),dpi=140)
    print("[sel] wrote figs + results/disp_point_selection.csv")

if __name__=="__main__": main()
