"""W の「必要な行だけ」を随伴法（単位荷重法）で速く求める.

W は 15,120行(変位自由度) × 5,040列(温度自由度)。
  ・列ごと（順方向）：節点 j の温度を 1 K にして解く → 5,040 回の求解が必要
  ・行ごと（随伴）  ：知りたい変位自由度 i に単位荷重を掛けて解く → 1 回の求解で W の i 行が丸ごと出る

K_s が対称なので
    W[i,:] = e_i^T K_s^{-1} H_T = (K_s^{-1} e_i)^T H_T = λ_i^T H_T
    （λ_i は「点 i に単位荷重を掛けたときの変位場」＝随伴場）
が成り立つ。つまり FrontISTR で単位荷重を掛けて1回解けば、その点の
「全 5,040 節点の温度に対する感度」が一度に得られる。

ここでは 105 のダンプ（K, H）を使って随伴法を実装し、
全体を構築した W の該当行と一致するか・どれだけ速いかを確かめる。

出力: results/W_rows_adjoint.json, docs/img/W_rows_adjoint.png
再現: OMP_NUM_THREADS=4 python3 run/W_rows_by_adjoint.py
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
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
P105=os.path.join(SAMPLE,"105_0_sensor_placement_sensitivity")
CASE=os.path.join(P105,"openfoam","dumpw_case")
NR,NTH,NZ=4,48,20; R_IN,R_OUT,H_=0.020,0.0375,0.1005
A_XYZ=(0.028,0.0,H_); O_XYZ=(-0.028,0.0,H_)


def main():
    mesh=cylinder_mesh.build_cylinder_mesh(NR,NTH,NZ,R_IN,R_OUT,H_)
    coords=np.array([x for _n,x in mesh["nodes"]])
    iA=int(np.linalg.norm(coords-np.array(A_XYZ),axis=1).argmin())
    iO=int(np.linalg.norm(coords-np.array(O_XYZ),axis=1).argmin())
    t0=time.time()
    K=sio.mmread(os.path.join(CASE,"dump_matrix_1_0.mm")).tocsc().astype(np.float64)
    Hm=sio.mmread(os.path.join(CASE,"H_matrix.mtx")).tocsc().astype(np.float64)
    t_read=time.time()-t0
    t0=time.time(); lu=splu(K); t_lu=time.time()-t0
    print(f"読み込み {t_read:.1f}s  LU分解 {t_lu:.1f}s  K{K.shape} H{Hm.shape}",flush=True)
    # ── 随伴法：知りたい2点（A・O の Uz）だけ ──
    dof=[3*iA+2,3*iO+2]
    t0=time.time()
    E=np.zeros((K.shape[0],len(dof)))
    for c,d in enumerate(dof): E[d,c]=1.0          # 単位荷重
    Lam=lu.solve(E)                                 # (15120, 2) 随伴場
    W_rows=(Lam.T@Hm)                               # (2, 5040) [m/K]
    t_adj=time.time()-t0
    print(f"随伴法（2点）: {t_adj:.3f}s  求解 {len(dof)} 回",flush=True)
    # ── 比較：全体を構築した W の該当行 ──
    W32=np.load(os.path.join(P105,"results","Wz_full.npy"))   # (5040節点, 5040温度) Uz行のみ
    ref=np.vstack([W32[iA],W32[iO]]).astype(np.float64)
    d=np.abs(W_rows-ref)
    valid=np.load(os.path.join(P105,"results","kinvh_sensitivity.npz"),allow_pickle=True)["valid"]
    out=dict(note="W の行を随伴法（単位荷重1回／点）で求め、全体構築版（105のWz_full）の同じ行と比較",
             node_A=iA,node_O=iO,
             time_read_s=round(t_read,1),time_lu_s=round(t_lu,1),
             time_adjoint_2pts_s=round(t_adj,3),time_full_W_s=54.4,
             max_abs_diff=float(d.max()),rel_to_rowmax=float(d.max()/np.abs(ref).max()),
             rowsum_adjoint_um_per_K=[float(np.abs(W_rows[k][valid]).sum()*1e6) for k in range(2)],
             rowsum_full_um_per_K=[float(np.abs(ref[k][valid]).sum()*1e6) for k in range(2)])
    print(json.dumps(out,ensure_ascii=False,indent=1))
    json.dump(out,open(os.path.join(RES,"W_rows_adjoint.json"),"w"),ensure_ascii=False,indent=1)
    # 図
    fig,axs=plt.subplots(1,2,figsize=(14,5.0))
    ax=axs[0]
    for k,nm,c in [(0,"A（+X 上面）","#C0392B"),(1,"O（−X 上面）","#2E6FD8")]:
        ax.plot(ref[k][valid]*1e6,W_rows[k][valid]*1e6,".",ms=3,color=c,label=nm)
    lim=[ref[:,valid].min()*1e6,ref[:,valid].max()*1e6]; ax.plot(lim,lim,"k--",lw=1)
    ax.set_xlabel("全体を構築した $W$ の行 [µm/K]"); ax.set_ylabel("随伴法で求めた行 [µm/K]")
    ax.set_title(f"一致の確認（拘束節点を除く 4,800 列）\n最大差 {d.max()*1e6:.2e} µm/K",fontsize=11.5)
    ax.legend(fontsize=10,markerscale=4); ax.grid(alpha=.3)
    ax=axs[1]
    nm=["随伴法\n（2点ぶん）","$W\\Phi$ を5列\n（モード5つ）","$W$ を全部\n（5,040列）"]
    vals=[t_adj,0.1,54.4]; cols=["#2E8B57","#E67E22","#C0392B"]
    ax.bar(range(3),vals,color=cols)
    for i,v in enumerate(vals): ax.text(i,v,f"{v:.2f} s",ha="center",va="bottom",fontsize=11)
    ax.set_yscale("log"); ax.set_xticks(range(3)); ax.set_xticklabels(nm,fontsize=10)
    ax.set_ylabel("求解にかかった時間 [s]（対数目盛）"); ax.grid(axis="y",alpha=.3)
    ax.set_title("求める範囲を絞ると速い（LU分解 2.4 s は共通）",fontsize=11.5)
    fig.suptitle("$W$ の必要な行だけを随伴法（単位荷重）で求める",fontsize=13)
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"W_rows_adjoint.png"),dpi=150); plt.close(fig)
    print("wrote W_rows_adjoint.png")


if __name__=="__main__": main()
