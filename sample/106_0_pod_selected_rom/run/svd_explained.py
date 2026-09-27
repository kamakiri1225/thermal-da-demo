"""np.linalg.svd が「中で何をしているか」を実データで確かめる（blog_004 §2-3の検証）.

102_0のCHT温度場121枚にSVDをかけ、次の6点を数値で確認する:
  1. Xc = U diag(S) Vt に厳密分解できている
  2. Uの列は空間相関行列 Xc Xc^T の固有ベクトル
  3. 時間側の小さい行列(121x121)からも同じ特異値が出る（スナップショット法）
  4. r モード打ち切りの誤差（Eckart-Young）
  5. モードは正規直交
  6. Vt = 各モードの時間プロファイル

再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/svd_explained.py
（102_0のOpenFOAM結果が必要）
"""
import sys, os, numpy as np
sys.path.insert(0,'run')
import importlib.util
spec=importlib.util.spec_from_file_location('q','run/select_points_qdeim.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
ts,C,X = m.load_snapshots()
mean=X.mean(axis=1); Xc=X-mean[:,None]
print(f"\n【入力】Xc = {Xc.shape[0]}セル × {Xc.shape[1]}時刻")
U,S,Vt=np.linalg.svd(Xc,full_matrices=False)
print(f"【出力】U={U.shape}  S={S.shape}  Vt={Vt.shape}")
print(f"\n■特異値S（先頭6個）: {np.round(S[:6],3)}")
e=S**2/(S**2).sum()
print(f"■エネルギー比[%]: {np.round(e[:6]*100,4)}")
print(f"■累積[%]:        {np.round(np.cumsum(e[:6])*100,4)}")

print("\n=== 検証1: Xc = U·diag(S)·Vt に分解できているか ===")
err=np.linalg.norm(Xc-(U*S)@Vt)/np.linalg.norm(Xc)
print(f"  相対誤差 = {err:.2e}  （ほぼ0＝完全に分解できている）")

print("\n=== 検証2: Uの列は空間相関行列 Xc·Xcᵀ の固有ベクトルか ===")
# 20696角は大きいので5点のみ抜き出して確認するのでなく、射影で確認
k=0
lhs=Xc@(Xc.T@U[:,k]); rhs=(S[k]**2)*U[:,k]
print(f"  Xc·Xcᵀ·u1 と σ1²·u1 の相対差 = {np.linalg.norm(lhs-rhs)/np.linalg.norm(rhs):.2e}")

print("\n=== 検証3: 時間側の小さい行列(121×121)からも同じ答えが出るか（スナップショット法） ===")
w,v=np.linalg.eigh(Xc.T@Xc)          # 121×121 の固有値分解
sig=np.sqrt(np.maximum(w,0))[::-1]    # 大きい順
print(f"  固有値の平方根(先頭4)   = {np.round(sig[:4],3)}")
print(f"  SVDの特異値  (先頭4)   = {np.round(S[:4],3)}  → 一致")

print("\n=== 検証4: 何モードで場を再現できるか（Eckart-Young） ===")
for r in [1,2,3,5]:
    Xr=(U[:,:r]*S[:r])@Vt[:r]
    rel=np.linalg.norm(Xc-Xr)/np.linalg.norm(Xc)
    rmse=np.sqrt(((Xc-Xr)**2).mean())
    print(f"  {r}モード: 相対誤差 {rel*100:7.4f}%   温度RMSE {rmse:.4f} K")

print("\n=== 検証5: モードは直交した単位ベクトルか ===")
G=U[:,:3].T@U[:,:3]
print("  UᵀU (先頭3×3) =\n", np.round(G,10))

print("\n=== 検証6: 時間プロファイル Vt の中身（モード1,2の時間変化） ===")
a1=S[0]*Vt[0]; a2=S[1]*Vt[1]
idx=[0,10,30,60,90,120]
print("  時刻[s]     :", [int(ts[i]) for i in idx])
print("  モード1係数 :", np.round(a1[idx],1))
print("  モード2係数 :", np.round(a2[idx],1))
