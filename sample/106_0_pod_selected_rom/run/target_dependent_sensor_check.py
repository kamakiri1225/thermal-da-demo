"""「推定したい量が変われば最適センサも変わる」を、変位の推定位置でも確かめる.

これまでの Δ の表は推定対象が T1/T2/Q/全温度場 だった。ここでは
**推定したい変位の場所・方向を変えたとき**に、最適な観測が入れ替わるかを見る。
入れ替わらなければ主張は弱いので、正直に出す。

観測の価値:  Δ(X,y) = Cov(X,y)^2 / (Var(y)+r)
共分散は、でたらめ初期からのアンサンブル（事前分布）を前進させた実サンプルから作る。

出力: docs/img/target_dependent_sensors.png, results/target_dependent_sensors.json
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/target_dependent_sensor_check.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0
N=4000; SEED=20260930; SIG_T=0.30; SIG_U=0.30
T_EVAL=150.0          # 事前分布を評価する時刻（加熱中）


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]
    Cc=kv["cell_centres"]; UPp=np.linalg.pinv(U[pod,:])
    op=np.load(os.path.join(RES,"dispop_xyz_4pts.npz"))
    um=op["u_mean"]; D=op["D"]          # (12,5)  A,O,B,C × Ux,Uy,Uz

    # ── 事前アンサンブル（でたらめ初期＋不確かなQ,h）を T_EVAL まで前進 ──
    rng=np.random.default_rng(SEED)
    T=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N,NPT))
    q=rng.uniform(0.3,1.8,N); h=np.clip(rng.normal(0.02,0.01,N),1e-3,0.1)
    tp=0.0
    for tb in np.arange(OBS_DT,T_EVAL+1e-9,OBS_DT):
        T=rg.integrate_ensemble(T,C,Km,h,q,heat,tp,tb,DT); tp=tb
    a=(T-mean[pod])@UPp.T               # (N,5) モード係数
    u=um+a@D.T                          # (N,12) 変位

    iA=dict(Ux=0,Uy=1,Uz=2); iO=dict(Ux=3,Uy=4,Uz=5)
    iB=dict(Ux=6,Uy=7,Uz=8); iC=dict(Ux=9,Uy=10,Uz=11)

    # ── 推定したい量（すべて変位。場所と方向を変える）──
    targets=[
      ("上面の反り  $U_z(A)-U_z(O)$\n(z=100.5 mm)",        u[:,iA["Uz"]]-u[:,iO["Uz"]]),
      ("中高さの反り $U_z(B)-U_z(C)$\n(z=75 mm)",           u[:,iB["Uz"]]-u[:,iC["Uz"]]),
      ("上面の伸び  $U_z(A)+U_z(O)$\n(全体の熱膨張)",        u[:,iA["Uz"]]+u[:,iO["Uz"]]),
      ("反ヒータ側の半径変位 $U_x(O)$",                       u[:,iO["Ux"]]),
      ("発熱量 $Q$",                                        q*15.0),
      ("全温度場（5点の平均）",                               T.mean(1)),
    ]
    # ── 観測の候補 ──
    cands=[(f"温度 P{i}\n{np.round(Cc[pod[i]]*1000,0).astype(int)}", T[:,i], SIG_T**2) for i in range(NPT)]
    cands+= [("変位 $U_z(A)$\n上面+X", u[:,iA["Uz"]], SIG_U**2),
             ("変位 $U_z(O)$\n上面−X", u[:,iO["Uz"]], SIG_U**2),
             ("変位 $U_z(B)$\n中高さ+X", u[:,iB["Uz"]], SIG_U**2)]

    M=np.zeros((len(targets),len(cands))); prior=np.zeros(len(targets))
    for i,(_tl,X) in enumerate(targets):
        prior[i]=X.var()
        for j,(_cl,y,r) in enumerate(cands):
            M[i,j]=np.cov(X,y)[0,1]**2/(y.var()+r)
    pct=M/prior[:,None]*100

    best=[int(np.argmax(M[i])) for i in range(len(targets))]
    print(f"事前アンサンブル N={N}, 評価時刻 t={T_EVAL:.0f}s\n")
    print(f"{'推定したい量':<34}{'最良の観測':<18}{'分散減少':>9}")
    for i,(tl,_X) in enumerate(targets):
        print(f"{tl.split(chr(10))[0]:<34}{cands[best[i]][0].split(chr(10))[0]:<18}{pct[i,best[i]]:8.1f}%")
    uniq=sorted(set(best))
    print(f"\n→ 最良の観測は {len(uniq)} 種類に分かれた（候補{len(cands)}個、対象{len(targets)}個）")

    fig,ax=plt.subplots(figsize=(15.6,7.4))
    ax.imshow(pct,cmap="RdYlGn",vmin=0,vmax=100,aspect="auto")
    ax.set_xticks(range(len(cands)))
    ax.set_xticklabels([c[0] for c in cands],fontsize=10.5)
    ax.set_yticks(range(len(targets)))
    ax.set_yticklabels([t[0] for t in targets],fontsize=10.5)
    for i in range(len(targets)):
        for j in range(len(cands)):
            mk="\n★最適" if j==best[i] else ""
            ax.text(j,i,f"{pct[i,j]:.0f} %{mk}",ha="center",va="center",
                    fontsize=11,weight=("bold" if j==best[i] else "normal"))
            if j==best[i]:
                ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,ec="#12308f",lw=3.5))
    ax.set_xlabel("観測の候補（どこを何で測るか）",fontsize=12.5)
    ax.set_title("推定したい「変位の場所・方向」を変えると、最適なセンサは入れ替わるか\n"
                 f"数値＝その観測で対象の分散が何%消えるか（$\\Delta$／事前分散）。"
                 f"★＝その行の最良。事前アンサンブル {N}本・t={T_EVAL:.0f} s",
                 fontsize=13.5,weight="bold")
    fig.text(0.5,0.012,"※ 行ごとに★の列が変われば「対象が変われば最適センサも変わる」の証拠になる。"
             "変わらなければ主張は成り立たない。",ha="center",fontsize=11,color="#444")
    fig.tight_layout(rect=[0,0.035,1,1])
    out=os.path.join(IMG,"target_dependent_sensors.png")
    fig.savefig(out,dpi=135); plt.close(fig)
    json.dump(dict(n_ensemble=N,t_eval=T_EVAL,
                   targets=[t[0].replace("\n"," ") for t in targets],
                   candidates=[c[0].replace("\n"," ") for c in cands],
                   delta=M.tolist(),prior_var=prior.tolist(),pct=pct.tolist(),
                   best_index=best,n_distinct_best=len(uniq)),
              open(os.path.join(RES,"target_dependent_sensors.json"),"w"),
              ensure_ascii=False,indent=2)
    print("wrote",out)


if __name__=="__main__":
    main()
