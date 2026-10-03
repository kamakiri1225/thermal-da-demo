"""blog_007 §8 用：代表点でないセルに温度計を置いた4配置の、位置図と温度の時刻歴.

run/obs_points_not_limited_to_rom.py と同じ4配置について
  左：どこに温度計を置いたか（真上から見た図＋高さ）
  右：代表的な未観測点の温度が同化でどう合うか（時刻歴）
を1枚にする。

出力: docs/img/obs_free_cells.png
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/make_obs_free_cells_fig.py
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
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60; SIG_T=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
FREE_SPLIT=[(0.034,0.0,0.095),(-0.034,0.0,0.030)]
FREE_SAME =[(0.034,0.0,0.095),(0.030,0.0,0.030)]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:]); A=U@UPp
    field=lambda T5: mean+A@(T5-mean[pod])
    near=lambda x: int(np.linalg.norm(Cc-np.array(x),axis=1).argmin())
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT); tg=np.r_[0,cyc]
    T=np.full(NPT,rg.T_AIR_K); Ttr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,DT); T=tr[-1]; Ttr.append(T.copy())
    Ttr=np.array(Ttr)
    # 評価用の未観測セル（どの配置でも観測していない）
    ev=[near((-0.0375,0.0,0.090)),near((0.0,0.0355,0.050))]
    def run(kind,idx,seed):
        rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
        Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
        Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
        Rd=np.eye(len(idx))*SIG_T**2; rec=[[mean[c]+A[c]@(Z[:,:NPT].mean(0)-mean[pod]) for c in ev]]; tp=0.0
        for ci,tb in enumerate(cyc,1):
            Z=Z.copy(); Z[:,:NPT]=rg.integrate_ensemble(Z[:,:NPT],C,Km,Z[:,IH],Z[:,IQ],heat,tp,tb,DT); tp=tb
            if kind=="pod": y=Ttr[ci][idx]; Yf=Z[:,idx]
            else:
                y=field(Ttr[ci])[idx]
                Yf=np.column_stack([mean[c]+(A[c]@(Z[:,:NPT]-mean[pod]).T) for c in idx])
            y=y+ro.normal(0,SIG_T,len(idx))
            Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
            Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
            mu=Z[:,:NPT].mean(0); rec.append([mean[c]+A[c]@(mu-mean[pod]) for c in ev])
        return np.array(rec)
    cfgs=[("代表点 P2+P0（どちらもヒータ側）","pod",[2,0],"#E67E22"),
          ("代表点 P2+P4（A側とO側に分ける）","pod",[2,4],"#2E8B57"),
          ("代表点でない2セル（A側とO側に分ける）","cell",[near(x) for x in FREE_SPLIT],"#2E6FD8"),
          ("代表点でない2セル（どちらもヒータ側）","cell",[near(x) for x in FREE_SAME],"#C0392B")]
    runs={nm:np.mean([run(k,i,s) for s in SEEDS],axis=0) for nm,k,i,_ in cfgs}
    truth=np.array([[mean[c]+A[c]@(t-mean[pod]) for c in ev] for t in Ttr])
    res=json.load(open(os.path.join(RES,"obs_points_not_limited_to_rom.json")))["configs"]
    fig=plt.figure(figsize=(16.5,5.4)); gs=fig.add_gridspec(1,3,width_ratios=[1,1.1,1.1],wspace=0.28)
    # 左：位置
    bx=fig.add_subplot(gs[0,0])
    # 横から見た図（x–z 断面）。高さが違う点を分離して描く
    bx.add_patch(plt.Rectangle((20,0),17.5,100.5,fc="#DCE3EA",ec="#7F8C8D"))
    bx.add_patch(plt.Rectangle((-37.5,0),17.5,100.5,fc="#DCE3EA",ec="#7F8C8D"))
    bx.add_patch(plt.Rectangle((38.2,25.25),2.6,50,fc="orangered",ec="none"))
    bx.text(42,50,"ヒータ",color="orangered",fontsize=10,va="center")
    bx.plot([-40,40],[0,0],color="k",lw=2.5); bx.text(0,-7,"底面（固定）",ha="center",fontsize=9.5)
    def put(x,z,mk,col,lab,off,fs=9.5,bold=False):
        bx.plot(x,z,mk,ms=11,color=col,mec="k")
        bx.annotate(lab,(x,z),textcoords="offset points",xytext=off,fontsize=fs,color=col,
                    weight=("bold" if bold else None))
    for i,lab,off in [(2,"P2",(7,-4)),(0,"P0",(-26,-4)),(4,"P4",(-24,2))]:
        p0=kv["xyz"][i]*1000; put(p0[0],p0[2],"o","#555",lab,off)
    put(Cc[near(FREE_SPLIT[0])][0]*1000,Cc[near(FREE_SPLIT[0])][2]*1000,"s","#2E6FD8","自由1",(7,-4))
    put(Cc[near(FREE_SPLIT[1])][0]*1000,Cc[near(FREE_SPLIT[1])][2]*1000,"s","#2E6FD8","自由2",(-32,-4))
    put(Cc[near(FREE_SAME[1])][0]*1000,Cc[near(FREE_SAME[1])][2]*1000,"s","#C0392B","自由3",(-32,8))
    put(Cc[ev[0]][0]*1000,Cc[ev[0]][2]*1000,"*","k","評価1",(-34,4),bold=True)
    put(0,Cc[ev[1]][2]*1000,"*","k","評価2\n(y=+35)",(-18,10),bold=True)
    bx.set_xlim(-58,62); bx.set_ylim(-14,112); bx.set_aspect("equal"); bx.axis("off")
    bx.set_title("温度計の位置（横から見た断面、寸法 mm）\n丸＝代表点、四角＝代表点でないセル、★＝評価点",fontsize=11.5)
    # 右2枚：時刻歴
    for k,(ax_i,ttl) in enumerate([(1,"評価点1（反ヒータ側 上部 z=90 mm）"),(2,"評価点2（+y 側 中段 z=50 mm）")]):
        ax=fig.add_subplot(gs[0,ax_i])
        ax.axvspan(0,300,color="#FDEBD0",alpha=.45)
        ax.plot(tg,truth[:,k]-273.15,"o-",color="k",lw=2.6,ms=4,label="真値")
        for nm,_kk,_ii,col in cfgs:
            ax.plot(tg,runs[nm][:,k]-273.15,"-",color=col,lw=1.8,label=f"{nm}（場 {res[nm]['field_K']:.3f} K）")
        ax.set_xlabel("時刻 [s]"); ax.set_ylabel("温度 [℃]"); ax.grid(alpha=.3); ax.set_title(ttl,fontsize=11.5)
        if k==0: ax.legend(fontsize=8.5,loc="lower right")
    fig.suptitle("代表点でないセルに温度計を置いても同化できる（観測は温度2点のみ、5 seed 平均）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"obs_free_cells.png"),dpi=150); plt.close(fig)
    print("wrote obs_free_cells.png")


if __name__=="__main__": main()
