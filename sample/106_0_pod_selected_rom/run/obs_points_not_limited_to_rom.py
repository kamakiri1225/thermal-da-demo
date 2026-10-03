"""温度の観測点は ROM 代表5点に限られないことを確かめる.

ROM の状態は代表5点の温度だが、任意セル j の温度は POD 復元で
    T_j = mean_j + (U_r の j 行) @ U_rP^+ @ (T_5 - mean_P)
と 5 点の温度の一次式で書ける。これを観測演算子に使えば、代表点でないセルに
置いた温度計もそのまま同化に入れられる。

比較する構成（いずれも温度2点のみ、変位なし）:
  A: 代表点 P2+P0（本番と同じ。どちらもヒータ側）
  B: 代表点 P2+P4（A側とO側に分けた配置）
  C: 代表点でない任意の2セル（上面ヒータ側と中段の反対側）
  D: 代表点でない任意の2セル（どちらもヒータ側＝偏った配置）

真値も同じ ROM の双子実験。EnKF・60メンバー・30秒ごと20回・5 seed。

出力: results/obs_points_not_limited_to_rom.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/obs_points_not_limited_to_rom.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import rom_general as rg
from dacore.enkf import enkf_update
RES=os.path.join(ROOT,"results")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
# 代表点でない候補（座標 [m] で指定し、最寄りセルを使う）
FREE_SPLIT=[(0.034,0.0,0.095),(-0.034,0.0,0.030)]    # ヒータ側の上部 と 反対側の中段
FREE_SAME =[(0.034,0.0,0.095),(0.030,0.0,0.030)]     # どちらもヒータ側


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:])
    field=lambda T5: mean+U@(UPp@(T5-mean[pod]))
    near=lambda xyz: int(np.linalg.norm(Cc-np.array(xyz),axis=1).argmin())
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)

    def run(kind,idx,seed):
        """kind='pod' なら idx は代表点番号、'cell' なら idx はセル番号。"""
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.eye(len(idx))*SIG_T**2; e=[]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if kind=="pod":
                y=Ttr[ci][idx]; Yf=Z[:,idx]
            else:   # 任意セル：POD 復元を観測演算子にする
                y=field(Ttr[ci])[idx]
                Yf=np.column_stack([mean[c]+(U[c]@UPp)@(Z[:,:NPT]-mean[pod]).T for c in idx])
            y=y+ro.normal(0,SIG_T,len(idx))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            if tb<=300:
                m=Z[:,:NPT].mean(0); e.append(np.sqrt(((field(m)-field(Ttr[ci]))**2).mean()))
        return float(np.mean(e)), float(15*Z[:,IQ].mean())

    cfgs=[("代表点 P2+P0（どちらもヒータ側）","pod",[2,0]),
          ("代表点 P2+P4（A側とO側に分ける）","pod",[2,4]),
          ("代表点でない2セル（A側とO側に分ける）","cell",[near(x) for x in FREE_SPLIT]),
          ("代表点でない2セル（どちらもヒータ側）","cell",[near(x) for x in FREE_SAME])]
    out={"note":"温度2点のみ（変位なし）。温度場RMSEは全20,696セル・加熱期30〜300秒・5 seed平均","configs":{}}
    for name,kind,idx in cfgs:
        r=np.array([run(kind,idx,s) for s in SEEDS])
        pos=[np.round(Cc[pod[i]]*1000,1).tolist() if kind=="pod" else np.round(Cc[i]*1000,1).tolist() for i in idx]
        out["configs"][name]=dict(kind=kind,index=[int(i) for i in idx],xyz_mm=pos,
                                  field_K=float(r[:,0].mean()),field_K_sd=float(r[:,0].std(ddof=1)),Q_W=float(r[:,1].mean()))
        print(f"{name:34s} 座標{pos}  温度場 {r[:,0].mean():.3f}±{r[:,0].std(ddof=1):.3f} K  Q {r[:,1].mean():.2f} W",flush=True)
    json.dump(out,open(os.path.join(RES,"obs_points_not_limited_to_rom.json"),"w"),ensure_ascii=False,indent=1)


if __name__=="__main__": main()
