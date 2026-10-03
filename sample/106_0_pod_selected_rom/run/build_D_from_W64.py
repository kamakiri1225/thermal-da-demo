"""float64 で W を作り直し、D=WΦ が FrontISTR 6回版と一致するかを確かめる.

105 の Wz_full.npy は float32 で保存されている。D=WΦ に 8.6% のずれが出たため、
丸め誤差が原因かを確かめる目的で K, H のダンプから float64 のまま W を作り直して比較した。

【結論（2026-10-03）】float32 は原因ではなかった。3通り（float32のW / float64のW / WΦを5列だけ求解）
すべて最大差 1.881e-03 µm で同一。真の原因は要素タイプの違い：
  106 の FrontISTR 6回版 … TYPE=361（六面体1次、3,840要素）
  105 の W              … TYPE=341（四面体1次、六面体を6分割して23,040要素）
四面体1次要素は曲げに硬いため変形が小さく出る。高次モードほど差が大きい（傾き 0.996→0.893）。
CG の許容誤差（1e-8 → 1e-14）も原因ではないことを確認済み。
メモリ節約のため W 全体は保存せず、必要な WΦ（5,040×5）だけを作る。

  W Φ = K_s^{-1} H_T Φ   なので、(H_T Φ) の5列だけを右辺にして LU 求解すれば5回で済む。
  これは「モードを入れて FrontISTR を5回解く」のと同じ計算を、手元の行列でやっていることになる。

出力: results/D_from_W64_check.json, docs/img/D_from_W64_check.png
再現: OMP_NUM_THREADS=4 python3 run/build_D_from_W64.py
"""
from __future__ import annotations
import os, sys, json, time
import numpy as np
import scipy.io as sio
from scipy.sparse.linalg import splu
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(SAMPLE,"102_1_frontistr_hollow_cylinder_thermal_expansion","python"))
from dacore import plots as _p
import matplotlib.pyplot as plt
import cylinder_mesh
from scipy.spatial import cKDTree
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
CASE=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","openfoam","dumpw_case")
W32=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results","Wz_full.npy")
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H_=0.020,0.0375,0.1005


def main():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(np.float64); Cc=kv["cell_centres"]
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H_)
    coords=np.array([x for _n,x in mesh["nodes"]])
    _,near=cKDTree(Cc).query(coords)
    Phi=U[near,:]                                    # (5040,5) モードをFEM節点温度へ
    t0=time.time()
    K=sio.mmread(os.path.join(CASE,"dump_matrix_1_0.mm")).tocsc()
    Hm=sio.mmread(os.path.join(CASE,"H_matrix.mtx")).tocsc()
    print(f"読み込み {time.time()-t0:.1f}s  K{K.shape} H{Hm.shape}  dtype K={K.dtype} H={Hm.dtype}",flush=True)
    t0=time.time(); lu=splu(K.astype(np.float64)); print(f"LU分解 {time.time()-t0:.1f}s",flush=True)
    # ① WΦ を5列だけ解く（＝モードを入れて解くのと同じ）
    t0=time.time()
    rhs=np.asarray((Hm@Phi).todense(),dtype=np.float64) if hasattr(Hm@Phi,"todense") else np.asarray(Hm@Phi,dtype=np.float64)
    Xmode=lu.solve(rhs)                              # (15120,5)
    D_mode=Xmode[2::3,:]*1e6                         # Uz成分 [µm]
    t_mode=time.time()-t0
    print(f"WΦ を5列で求解 {t_mode:.1f}s",flush=True)
    # ② W を float64 で丸ごと作ってから掛ける
    t0=time.time()
    n_t=Hm.shape[1]; Wz=np.zeros((n_t,n_t),dtype=np.float64)
    for c0 in range(0,n_t,504):
        c1=min(c0+504,n_t)
        X=lu.solve(np.asarray(Hm[:,c0:c1].todense(),dtype=np.float64))
        Wz[:,c0:c1]=X[2::3,:]
    t_full=time.time()-t0
    D_W64=(Wz@Phi)*1e6
    print(f"W を丸ごと構築 {t_full:.1f}s",flush=True)
    # ③ float32 版（105 の保存値）
    D_W32=(np.load(W32).astype(np.float64)@Phi)*1e6
    D_fem=np.load(os.path.join(RES,"dispop_allnodes.npz"))["D"][:,2,:]
    valid=np.load(os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity","results",
                               "kinvh_sensitivity.npz"),allow_pickle=True)["valid"]
    out=dict(note="D=WΦ の3通りと FrontISTR 6回版の比較（Uz成分、[µm]）。valid＝拘束節点を除く4,800節点",
             time_mode_solve_s=round(t_mode,1),time_full_W_s=round(t_full,1),
             W64_MB=round(Wz.nbytes/1048576,1))
    for nm,Dx in [("WΦを5列で求解(float64)",D_mode),("Wをfloat64で構築",D_W64),("Wをfloat32で構築(105保存値)",D_W32)]:
        d=np.abs(Dx[valid]-D_fem[valid])
        out[nm]=dict(max_abs_diff_um=float(d.max()),mean_abs_diff_um=float(d.mean()),
                     max_rel_pct=float(100*d.max()/np.abs(D_fem[valid]).max()))
        print(f"{nm:30s} 最大差 {d.max():.3e} µm  平均 {d.mean():.3e}  最大相対 {100*d.max()/np.abs(D_fem[valid]).max():.3f}%",flush=True)
    json.dump(out,open(os.path.join(RES,"D_from_W64_check.json"),"w"),ensure_ascii=False,indent=1)
    # 図
    fig,axs=plt.subplots(1,2,figsize=(13.5,5.0))
    ax=axs[0]
    labs=[("float32 の W",D_W32,"#C0392B"),("float64 の W",D_W64,"#E67E22"),("5列だけ求解（float64）",D_mode,"#2E8B57")]
    for nm,Dx,c in labs:
        ax.plot(D_fem[valid].ravel(),Dx[valid].ravel(),".",ms=2,color=c,label=nm)
    lim=[D_fem[valid].min()*1.05,D_fem[valid].max()*1.05]; ax.plot(lim,lim,"k--",lw=1)
    ax.set_xlabel("FrontISTR 6回で作った $D$ [µm]"); ax.set_ylabel("$W\\Phi$ で作った $D$ [µm]")
    ax.set_title("拘束節点を除く 4,800 節点 × 5 モード",fontsize=11.5); ax.legend(fontsize=9,markerscale=5); ax.grid(alpha=.3)
    ax=axs[1]
    nms=[l[0] for l in labs]; vals=[out[k]["max_abs_diff_um"] for k in ["Wをfloat32で構築(105保存値)","Wをfloat64で構築","WΦを5列で求解(float64)"]]
    ax.bar(range(3),vals,color=[l[2] for l in labs])
    for i,v in enumerate(vals): ax.text(i,v,f"{v:.2e}",ha="center",va="bottom",fontsize=10)
    ax.set_yscale("log"); ax.set_xticks(range(3)); ax.set_xticklabels(nms,fontsize=9.5)
    ax.set_ylabel("FrontISTR 6回版との最大差 [µm]"); ax.grid(axis="y",alpha=.3)
    ax.set_title("精度の比較（対数目盛）",fontsize=11.5)
    fig.suptitle("$D=W\\Phi$ は FrontISTR 6回版と一致するか ― float32 と float64 の違い",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"D_from_W64_check.png"),dpi=150); plt.close(fig)
    print("wrote D_from_W64_check.png")


if __name__=="__main__": main()
