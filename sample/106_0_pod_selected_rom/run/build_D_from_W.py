"""熱感度行列 W から D=WΦ を計算し、FrontISTR 6回で作った D と一致するか確かめる.

問い：D は「モードを入れて FrontISTR を6回解く」以外に、
      105 で構築済みの W = K_s^{-1} H_T を使って D = W Φ でも作れるのでは？

ここで
  W  … 105_0 の results/Wz_full.npy（Uz 行のみ、5,040節点 × 5,040温度自由度 [m/K]）
  Φ  … 106 の POD モードを FEM 節点温度に移したもの（5,040 × 5）
を使い、D = W Φ を計算して、FrontISTR 6回版（results/dispop_allnodes.npz の Uz 成分）と比べる。

出力: results/D_from_W_check.json, docs/img/D_from_W_check.png
再現: OMP_NUM_THREADS=4 python3 run/build_D_from_W.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
W105=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","Wz_full.npy")
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H=0.020,0.0375,0.1005


def main():
    if not os.path.exists(W105):
        print("W が無い（105_0 の run/build_W_from_dumps.py が必要）:",W105); return
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); Cc=kv["cell_centres"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H)
    coords=np.array([x for _n,x in mesh["nodes"]])
    _,near=cKDTree(Cc).query(coords)                 # FEM節点 → 最寄りOpenFOAMセル
    t0=time.time()
    Wz=np.load(W105)                                  # (5040節点, 5040温度) [m/K]
    print(f"W 読み込み {Wz.shape} {time.time()-t0:.1f}s",flush=True)
    Phi=U[near,:]                                     # (5040, 5) モードをFEM節点へ
    D_W=(Wz@Phi)*1e6                                  # (5040, 5) [µm]
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    D_fem=op["D"][:,2,:]                              # Uz 成分 (5040, 5) [µm]
    coords2=op["coords"]
    assert np.allclose(coords,coords2), "節点順が違う"
    # 比較
    num=np.abs(D_W-D_fem); den=np.abs(D_fem).max()
    out=dict(note="D=WΦ（105のW使用）と、FrontISTR 6回で作ったD（Uz成分）の比較 [µm]",
             shape=list(D_W.shape),
             max_abs_diff_um=float(num.max()),mean_abs_diff_um=float(num.mean()),
             max_rel_diff=float(num.max()/den),
             D_fem_absmax_um=float(den),
             corr=[float(np.corrcoef(D_W[:,k],D_fem[:,k])[0,1]) for k in range(5)])
    print(json.dumps(out,ensure_ascii=False,indent=1))
    json.dump(out,open(os.path.join(RES,"D_from_W_check.json"),"w"),ensure_ascii=False,indent=1)
    fig,axs=plt.subplots(1,3,figsize=(15.5,4.8))
    ax=axs[0]
    for k in range(5): ax.plot(D_fem[:,k],D_W[:,k],".",ms=2,label=f"モード{k+1}")
    lim=[D_fem.min()*1.05,D_fem.max()*1.05]; ax.plot(lim,lim,"k--",lw=1)
    ax.set_xlabel("FrontISTR 6回で作った $D$ [µm]"); ax.set_ylabel("$W\\Phi$ で作った $D$ [µm]")
    ax.set_title(f"全 5,040 節点 × 5 モード\n最大差 {num.max():.2e} µm",fontsize=11.5)
    ax.legend(fontsize=8.5,markerscale=4); ax.grid(alpha=.3)
    ax=axs[1]
    ax.hist(num.ravel(),bins=60,color="#2E6FD8")
    ax.set_xlabel("差の絶対値 [µm]"); ax.set_ylabel("個数"); ax.grid(alpha=.3)
    ax.set_title(f"差の分布（平均 {num.mean():.2e} µm）",fontsize=11.5)
    ax=axs[2]
    z=coords[:,2]*1000
    ax.plot(z,num.max(1),".",ms=3,color="#C0392B")
    ax.set_xlabel("高さ z [mm]"); ax.set_ylabel("差の最大 [µm]"); ax.grid(alpha=.3)
    ax.set_title("節点ごとの差",fontsize=11.5)
    fig.suptitle("$D=W\\Phi$ は FrontISTR 6回版と一致するか（Uz 成分、106 の POD モード）",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"D_from_W_check.png"),dpi=150); plt.close(fig)
    print("wrote D_from_W_check.png")


if __name__=="__main__": main()
