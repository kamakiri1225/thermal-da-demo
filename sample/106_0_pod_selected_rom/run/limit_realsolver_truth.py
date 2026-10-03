"""限界の検証：真値を ROM ではなく実ソルバ（OpenFOAM の温度場＋FrontISTR の変位）にして同化する.

これまでの双子実験は、真値も同化モデルも同じ ROM で、変位の真値も同じ近似写像（5点→PODモード→D）で作っていた。
ここでは
  ・温度の真値     ＝ OpenFOAM(CHT) の固体温度場そのもの（5点はその値、全体は 20,696 セル）
  ・変位の真値     ＝ その温度場全体を FrontISTR に入れて直接計算した上面 A/O の Uz
  ・同化モデル     ＝ これまでと同じ ROM（発熱は 0〜300 s に一定、と仮定したまま）
とし、ROM の近似誤差・変位写像の近似誤差を含めて精度を評価する。

真値のケース（CASE）:
  learned      … 102_0 の CHT（15 W, 0〜300 s）。ROM の校正・POD に使った条件そのもの
  q25          … 25 W, 0〜300 s（新しく計算した条件発熱量）            openfoam/limit/q25
  intermittent … 15 W を 0〜150 s と 300〜450 s（新しく計算した条件加熱履歴。ROM の仮定と異なる）
                                                                   openfoam/limit/intermittent

出力: results/limit_realsolver_<CASE>.json, results/limit_truth_<CASE>.npz（FrontISTR 真値のキャッシュ）
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/limit_realsolver_truth.py learned
"""
from __future__ import annotations
import os, sys, json
from pathlib import Path
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from select_points_qdeim import read_foam_field
import cylinder_mesh, fistr_case
from fem.fem_obs import MATERIAL
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005
Tref=MATERIAL["reference_temperature_K"]
A_XYZ=np.array([0.028,0,H]); O_XYZ=np.array([-0.028,0,H])
CASES={"learned":(os.path.join(SAMPLE,"102_0_openfoam_hollow_cylinder_heat_transfer"),"15 W・0〜300 s（ROMの校正に使った条件）"),
       "q25":(os.path.join(ROOT,"openfoam","limit","q25"),"25 W・0〜300 s（新しく計算した条件発熱量）"),
       "intermittent":(os.path.join(ROOT,"openfoam","limit","intermittent"),"15 W・0〜150 s と 300〜450 s（新しく計算した条件加熱履歴）")}
CONFIGS=[("同化なし",None,False),("温度2点 P2+P0（A側だけ）",[2,0],False),("温度2点 P2+P4（両側）",[2,4],False),
         ("温度 P2+P0＋変位 A/O",[2,0],True),("温度 P2+P4＋変位 A/O",[2,4],True)]


def truth(case):
    cdir,_=CASES[case]; cache=os.path.join(RES,f"limit_truth_{case}.npz")
    kv=np.load(os.path.join(RES,"qdeim_points.npz")); Cc=kv["cell_centres"]; pod=kv["cell_idx"]
    times=np.r_[0,np.arange(OBS_DT,T_END+1e-9,OBS_DT)]
    if os.path.exists(cache):
        d=np.load(cache); return d["times"],d["Tfield"],d["uz"]
    Tf=[]
    for t in times:
        p=os.path.join(cdir,"0" if t==0 else f"{t:g}","solid","T")
        Tf.append(read_foam_field(p,len(Cc)))
    Tf=np.array(Tf)
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([x for _n,x in mesh["nodes"]]); node_ids=[n for n,_ in mesh["nodes"]]
    qa=int(np.linalg.norm(coords-A_XYZ,axis=1).argmin()); qo=int(np.linalg.norm(coords-O_XYZ,axis=1).argmin())
    _,near=cKDTree(Cc).query(coords)
    work=Path(ROOT,"openfoam","limit",f"fistr_{case}"); work.mkdir(parents=True,exist_ok=True)
    fistr_case.write_mesh(work,NR,NTH,NZ,R_IN,R_OUT,H,young_modulus=MATERIAL["young_modulus_Pa"],poisson_ratio=MATERIAL["poisson_ratio"],
        density=MATERIAL["density_kg_m3"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
    fistr_case.write_hecmw_ctrl(work)
    uz=[]
    for k,t in enumerate(times):
        fistr_case.write_cnt(work,node_ids,Tf[k][near],reference_temperature=Tref,young_modulus=MATERIAL["young_modulus_Pa"],
            poisson_ratio=MATERIAL["poisson_ratio"],thermal_expansion_coeff=MATERIAL["thermal_expansion_coeff_per_K"])
        fistr_case.run_fistr(work); disp=fistr_case.read_displacement(work)
        uz.append([disp[node_ids[qa]][2]*1e6,disp[node_ids[qo]][2]*1e6]); print(f"[truth:{case}] FrontISTR t={t:g}s",flush=True)
    uz=np.array(uz); np.savez(cache,times=times,Tfield=Tf,uz=uz); return times,Tf,uz


def main(case):
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz")); U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    UPp=np.linalg.pinv(U[pod,:]); op=np.load(os.path.join(RES,"disp_operator.npz")); uzm=op["uz_mean"]; Dm=op["Dmode"]
    disp=lambda T5: uzm+((T5-mean[pod])@UPp.T)@Dm.T          # 同化モデル側の変位（近似写像）
    field=lambda T5: mean+U@(UPp@(T5-mean[pod]))              # 5点 → 全体（gappy-POD）
    times,Tf,uz_true=truth(case); T5_true=Tf[:,pod]; cyc=times[1:]
    # 近似そのものの誤差（真の5点温度をそのまま入れた場合の上限性能）
    floor_field=np.mean([np.sqrt(((field(T5_true[k])-Tf[k])**2).mean()) for k in range(1,len(times))])
    floor_AO=np.mean([abs((disp(T5_true[k])[0]-disp(T5_true[k])[1])-(uz_true[k,0]-uz_true[k,1])) for k in range(1,len(times))])
    def run(nodes,use_disp,seed):
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG))
        Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        if nodes is not None: Rd=np.diag([SIG_T**2]*len(nodes)+([SIG_U**2]*2 if use_disp else []))
        e5=[];eF=[];eAO=[];tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if nodes is not None:
                yv=list(T5_true[ci][nodes]); Yf=Z[:,nodes]
                if use_disp: yv+=list(uz_true[ci]); Yf=np.column_stack([Yf,disp(Z[:,:NPT])])
                y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
                Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
                Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            m=Z[:,:NPT].mean(0); u=disp(m)
            e5.append(np.sqrt(((m-T5_true[ci])**2).mean())); eF.append(np.sqrt(((field(m)-Tf[ci])**2).mean()))
            eAO.append(abs((u[0]-u[1])-(uz_true[ci,0]-uz_true[ci,1])))
        return np.array(e5),np.array(eF),np.array(eAO),15*Z[:,IQ].mean()
    heatwin=cyc<=300; allwin=np.ones_like(cyc,bool)
    out=dict(case=case,label=CASES[case][1],
             note="誤差は5 seed平均。heat=30〜300 s、all=30〜600 s の平均。floor＝真の5点温度をそのまま入れたときの近似誤差",
             floor_field_K=float(floor_field),floor_AO_um=float(floor_AO),
             truth_AO_um=(uz_true[:,0]-uz_true[:,1]).tolist(),times=times.tolist(),configs={})
    print(f"== {case}: {CASES[case][1]}")
    print(f"  近似の誤差（真の5点温度を入れた場合）: 温度場 {floor_field:.3f} K, 反り {floor_AO:.3f} µm, 反りの最大 {np.abs(uz_true[:,0]-uz_true[:,1]).max():.2f} µm")
    for name,nodes,ud in CONFIGS:
        rs=[run(nodes,ud,s) for s in SEEDS]
        r={}
        for wn,w in [("heat",heatwin),("all",allwin)]:
            r[wn]=dict(T5_K=float(np.mean([x[0][w].mean() for x in rs])),field_K=float(np.mean([x[1][w].mean() for x in rs])),
                       AO_um=float(np.mean([x[2][w].mean() for x in rs])),AO_um_sd=float(np.std([x[2][w].mean() for x in rs],ddof=1)))
        r["Q_final_W"]=float(np.mean([x[3] for x in rs]))
        out["configs"][name]=r
        print(f"  {name:22s} 5点 {r['heat']['T5_K']:.3f}/{r['all']['T5_K']:.3f} K  温度場 {r['heat']['field_K']:.3f}/{r['all']['field_K']:.3f} K  反り {r['heat']['AO_um']:.3f}/{r['all']['AO_um']:.3f} µm  Q {r['Q_final_W']:.2f} W",flush=True)
    json.dump(out,open(os.path.join(RES,f"limit_realsolver_{case}.json"),"w"),ensure_ascii=False,indent=1)


if __name__=="__main__": main(sys.argv[1] if len(sys.argv)>1 else "learned")
