"""バグ修正: OI各runの u_analysis を「保存済みの同化後温度場」から再計算して npz を直す.

現行ドライバの旧不具合で u_analysis に u_forecast（補正前の場の変位）が入っていた。
同化後の温度場は run_fem_oi_<tag>/{t}/solid/T に正しく保存されているので、
そこへ FrontISTR(displacement_obs) を掛け直して u_analysis を上書きする（OpenFOAMの再計算は不要）。

対象: run_fem_oi_b(σ_Q=2), run_fem_oi_sq0p5(σ_Q=0.5), run_fem_oi_sq50(σ_Q=50)
   free(同化なし)は u=予報のままが正しいので対象外。run1(σ_Q=6)は旧正しいコードなので対象外。
再現: OMP_NUM_THREADS=4 python3 run/fix_oi_disp.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from fem.fem_obs import displacement_obs
RES=os.path.join(ROOT,"results"); OF=os.path.join(ROOT,"openfoam")
TARGETS=[("b","run_fem_oi_b"),("sq0p5","run_fem_oi_sq0p5"),("sq50","run_fem_oi_sq50")]


def _tname(t): return str(int(t)) if float(t)==int(t) else ("%g"%float(t))


def main():
    for tag,case in TARGETS:
        npz=os.path.join(RES,f"oi_fullsolver_{tag}_history.npz")
        casedir=os.path.join(OF,case)
        if not os.path.exists(npz): print(f"[{tag}] npzなし→skip"); continue
        h=dict(np.load(npz,allow_pickle=True))
        t=h["times"]; ua=h["u_analysis"].copy()
        changed=0
        for i,ti in enumerate(t):
            if ti<=0: continue                                   # t=0 は 0 のまま
            Tf=os.path.join(casedir,_tname(ti),"solid","T")
            if not os.path.exists(Tf): print(f"[{tag}] t={ti:g} 場なし→そのまま"); continue
            u_new=displacement_obs(casedir,_tname(ti),os.path.join(casedir,f"fem_afix_t{ti:g}"))
            old=ua[i].copy(); ua[i]=u_new
            if not np.allclose(old,u_new,atol=1e-9): changed+=1
            print(f"[{tag}] t={ti:g}s  A: {old[0]*1000:6.2f}→{u_new[0]*1000:6.2f}  O: {old[1]*1000:6.2f}→{u_new[1]*1000:6.2f} um")
        h["u_analysis"]=ua
        np.savez(npz,**h)
        print(f"[{tag}] DONE 上書き（{changed}/{len(t)-1}点を修正） -> {os.path.basename(npz)}\n")


if __name__=="__main__": main()
