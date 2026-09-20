"""blog_005用：観測の価値 Δ(X,y)=Cov²/(Var+r) の4×3表をヒートマップで可視化.

3変数モデル x=(T1,T2,Q), B=[[1.36,0.12,1.2],[0.12,1.04,0.4],[1.2,0.4,4]],
候補 y∈{温度点1, 温度点2, 変位u(w=(0.5,0.2))}, r_T=0.09, r_u=0.01。
行=推定対象、★=行内最大（その対象に最適な観測）→ 行ごとに★の列が変わる。
出力: docs/img/blog005_delta_table.png
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
IMG=os.path.join(ROOT,"docs","img")

B=np.array([[1.36,0.12,1.2],[0.12,1.04,0.4],[1.2,0.4,4.0]])
w=np.array([0.5,0.2,0.0]); rT=0.09; ru=0.01
cands=[("温度点1",np.array([1,0,0]),rT),("温度点2",np.array([0,1,0]),rT),("変位 u",w,ru)]
targets=[("T1（点1の温度）",np.array([1,0,0])),
         ("T2（点2の温度）",np.array([0,1,0])),
         ("発熱量 Q",np.array([0,0,1])),
         ("全温度場（T1+T2）",None)]

def delta(cvec,h,r):
    var=h@B@h+r
    if cvec is None:   # 全場: T1とT2の分散減少の和
        return sum((np.eye(3)[i]@B@h)**2 for i in (0,1))/var
    return (cvec@B@h)**2/var

M=np.array([[delta(cv,h,r) for _,h,r in cands] for _,cv in targets])
fig,ax=plt.subplots(figsize=(10.5,5.6))
Mn=M/M.max(1,keepdims=True)
ax.imshow(Mn,cmap="RdYlGn",vmin=0,vmax=1,aspect="auto")
ax.set_xticks(range(3)); ax.set_xticklabels([c[0] for c in cands],fontsize=13,weight="bold")
ax.set_yticks(range(4)); ax.set_yticklabels([t[0] for t in targets],fontsize=12)
for i in range(4):
    b=int(np.argmax(M[i]))
    for j in range(3):
        ax.text(j,i,f"{M[i,j]:.3f}"+("\n★最適" if j==b else ""),ha="center",va="center",
                fontsize=12,weight=("bold" if j==b else "normal"))
        if j==b: ax.add_patch(plt.Rectangle((j-.5,i-.5),1,1,fill=False,ec="blue",lw=3))
ax.set_title("観測の価値 $\\Delta$（＝その観測で対象の分散がどれだけ減るか）\n"
             "行＝推定したい量、列＝観測の候補。★の列が行ごとに違う＝最適センサは対象で変わる",
             fontsize=12.5,weight="bold")
fig.tight_layout()
out=os.path.join(IMG,"blog005_delta_table.png"); fig.savefig(out,dpi=150,bbox_inches="tight")
print("wrote",out); print(np.round(M,4))
