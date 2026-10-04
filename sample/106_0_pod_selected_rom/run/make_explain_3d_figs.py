"""解説用の図（3D の位置図＋グラフ）を作る.

数値の表だけでは「どこに何があり、何が起きているか」が分からないので、
pyvista（ParaView と同じ VTK）で中空円筒・代表点・センサを3Dで描き、
matplotlib のグラフと組み合わせる。3D 側の文字は記号（P0, S1 など）だけにし、
日本語の説明は matplotlib 側に書く（VTK は日本語を描けないため）。

  explain_free_cell_weights.png  … blog_006 §7-0：代表点以外の温度計は、どの代表点にどれだけ乗るか
  explain_weighted_sum.png       … blog_006 §3-3：「ずれ × 重み → 足す」でセンサの読みを予想する
  explain_correction_by_sensor.png … blog_006 §3-4：各センサのずれが、代表5点をどれだけ直したか
  explain_after_correction.png   … blog_006 §3-4 段階4：補正の前後と、そこから計算した次の 30 秒

数値は blog_006 と同じ計算（seed 20260913、t=30 秒の1回目の補正、真値＝ROM の双子実験）。
再現: OMP_NUM_THREADS=4 python3 run/make_explain_3d_figs.py
"""
from __future__ import annotations
import os, sys
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
import pyvista as pv
pv.OFF_SCREEN=True
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
R_IN,R_OUT,H=20.0,37.5,100.5
RED="#C0392B"; GREEN="#1F9D62"; BLUE="#2E6FD8"; ORANGE="#E67E22"; GRAY="#8696a7"; NAVY="#0B2545"
plt.rcParams.update({"font.size":15,"axes.titlesize":17,"axes.labelsize":15})


def base_scene(p):
    """中空円筒（半透明）とヒータ面（赤）"""
    tube=pv.CylinderStructured(radius=np.linspace(R_IN,R_OUT,2),height=H,center=(0,0,H/2),
                               direction=(0,0,1),theta_resolution=96,z_resolution=2)
    p.add_mesh(tube.extract_surface(algorithm="dataset_surface"),color="#d8dee6",opacity=0.28,show_edges=False)
    p.add_mesh(tube.extract_feature_edges(feature_angle=60),color="#7c8a99",line_width=1.5)
    th=np.radians(np.linspace(-76,76,40)); z=np.linspace(25,75,2)
    TH,Z=np.meshgrid(th,z); X=(R_OUT+0.3)*np.cos(TH); Y=(R_OUT+0.3)*np.sin(TH)
    p.add_mesh(pv.StructuredGrid(X,Y,Z),color="#e74c3c",opacity=0.30)
    p.add_point_labels([(R_OUT+14,0,50)],["HEATER"],font_size=22,text_color="#c0392b",
                       shape=None,show_points=False,always_visible=True)


def shot(p,name):
    p.camera_position=[(185,-175,125),(0,0,50),(0,0,1)]
    p.camera.zoom(0.85)
    f=os.path.join(ROOT,"results",f"_{name}.png"); p.screenshot(f,transparent_background=False); p.close()
    return plt.imread(f)


def load():
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]*1000
    return U,mean,pod,Cc


def fig_free_weights():
    U,mean,pod,Cc=load(); inv=np.linalg.inv(U[pod])
    P=Cc[pod]
    free={"F1":(35.1,-0.8,94.6),"F2":(-32.9,-1.6,29.7)}
    G={}
    for k,x in free.items():
        j=int(np.linalg.norm(Cc-np.array(x),axis=1).argmin()); G[k]=(Cc[j],U[j]@inv)
    p=pv.Plotter(off_screen=True,window_size=(1100,1000)); p.set_background("white"); base_scene(p)
    for i,x in enumerate(P):
        p.add_mesh(pv.Sphere(radius=3.2,center=x),color="#c0392b")
    p.add_point_labels(P+np.array([0,0,7]),[f"P{i}" for i in range(5)],font_size=26,text_color="#7a1f14",
                       shape=None,show_points=False,always_visible=True,bold=True)
    for k,(x,g) in G.items():
        p.add_mesh(pv.Cube(center=x,x_length=6,y_length=6,z_length=6),color="#2e6fd8")
        p.add_point_labels([x+np.array([0,0,8])],[k],font_size=28,text_color="#1f4e9c",shape=None,
                           show_points=False,always_visible=True,bold=True)
        for i in range(5):
            if abs(g[i])>=0.1:
                p.add_mesh(pv.Line(x,P[i]),color="#1f4e9c",line_width=4+14*abs(g[i]))
                mid=(x+P[i])/2
                p.add_point_labels([mid],[f"{g[i]:.2f}"],font_size=24,text_color="#1f4e9c",shape="rounded_rect",
                                   shape_color="white",shape_opacity=0.85,show_points=False,always_visible=True)
    img=shot(p,"free")
    fig=plt.figure(figsize=(16,7.8))
    ax=fig.add_axes([0.0,0.02,0.52,0.86]); ax.imshow(img); ax.axis("off")
    ax.set_title("● 赤＝代表点 P0〜P4　■ 青＝代表点以外の温度計 F1・F2\n線＝重み（太いほど強く乗る。0.1 未満は省略）",fontsize=14)
    ax2=fig.add_axes([0.58,0.14,0.40,0.66])
    names=["P0","P1","P2","P3","P4"]; x=np.arange(5); w=0.26
    unit=np.eye(5)[2]
    ax2.bar(x-w,unit,w,color=RED,label="代表点 P2 の温度計")
    ax2.bar(x,G["F1"][1],w,color=BLUE,label="F1（ヒータ側の上部）")
    ax2.bar(x+w,G["F2"][1],w,color="#79a7e8",label="F2（反対側の下部）")
    ax2.axhline(0,color="k",lw=1); ax2.set_xticks(x); ax2.set_xticklabels(names)
    ax2.set_ylabel("重み g（その代表点に何割乗るか）"); ax2.set_ylim(-0.2,1.15); ax2.grid(axis="y",alpha=.3)
    ax2.legend(loc="upper left",fontsize=12.5)
    ax2.set_title("温度計の読み ＝ 5点の温度 × 重み の合計",fontsize=15)
    fig.suptitle("代表点以外に置いた温度計も、代表5点の温度で書ける（gappy-POD の1行＝重み）",fontsize=18,weight="bold")
    fig.text(0.58,0.03,"代表点 P2 の温度計は「P2 だけ 1」。\nF1 はほぼ P3 に、F2 は P4 と P1 に分かれて乗る",fontsize=13,color=NAVY)
    fig.savefig(os.path.join(IMG,"explain_free_cell_weights.png"),dpi=130,facecolor="white"); plt.close(fig)
    print("wrote explain_free_cell_weights.png")


def fig_weighted_sum():
    T=np.array([297.366,295.304,298.191,296.583,293.997]); mp=np.array([296.659,295.698,296.938,296.146,295.361])
    d=T-mp; w=np.array([-0.532,0.596,-0.763,0.286,0.290]); pr=d*w
    names=["P0","P1","P2","P3","P4"]; x=np.arange(5)
    fig,axs=plt.subplots(1,4,figsize=(17,5.6),gridspec_kw=dict(width_ratios=[1,1,1,0.8]))
    for ax,v,c,t,u in [(axs[0],d,ORANGE,"① 5点の温度のずれ\n（予報 − 平均）","K"),
                       (axs[1],w,GREEN,"② 変位計 S1 の重み\n（1 K あたり何 µm）","µm/K"),
                       (axs[2],pr,BLUE,"③ ① × ②","µm")]:
        ax.bar(x,v,color=c); ax.axhline(0,color="k",lw=1)
        for i,y in enumerate(v): ax.text(i,y+(0.06 if y>=0 else -0.06)*max(abs(v)),f"{y:+.2f}",ha="center",
                                         va="bottom" if y>=0 else "top",fontsize=12.5)
        ax.set_xticks(x); ax.set_xticklabels(names); ax.set_title(t,fontsize=15); ax.set_ylabel(u)
        ax.set_ylim(-1.35*max(abs(v)),1.35*max(abs(v))); ax.grid(axis="y",alpha=.3)
    ax=axs[3]; ax.axis("off")
    ax.text(0.5,0.80,"④ ③ を全部足す",ha="center",fontsize=16,weight="bold",transform=ax.transAxes)
    ax.text(0.5,0.62,f"{pr.sum():+.2f} µm",ha="center",fontsize=22,color=BLUE,weight="bold",transform=ax.transAxes)
    ax.text(0.5,0.46,"＋ 平均温度のときの変位\n−1.74 µm",ha="center",fontsize=14,transform=ax.transAxes)
    ax.text(0.5,0.18,"S1 は −3.58 µm を\n示すはず",ha="center",fontsize=19,color=RED,weight="bold",transform=ax.transAxes)
    fig.suptitle("センサの読みの予想 ＝「5点の温度のずれ × 重み」を足すだけ（メンバー1、t=30 秒、代表点ではない変位計 S1）",
                 fontsize=16,weight="bold")
    fig.tight_layout(rect=(0,0,1,0.92)); fig.savefig(os.path.join(IMG,"explain_weighted_sum.png"),dpi=130,facecolor="white"); plt.close(fig)
    print("wrote explain_weighted_sum.png")


def kalman_example():
    """blog_006 §3-4 と同じ計算（seed 20260913、t=30 秒）"""
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])
    U,mean,pod,Cc=load(); UPp=np.linalg.pinv(U[pod]); mp=mean[pod]
    op=np.load(os.path.join(RES,"dispop_allnodes.npz")); um=op["u_mean"]; Dall=op["D"]; co=op["coords"]
    Wall=np.einsum("nck,kj->ncj",Dall,UPp)
    near=lambda x:int(np.linalg.norm(co-np.array(x),axis=1).argmin())
    sel=[(near((0.0097,-0.0362,0.1005)),0),(near((-0.0346,0.0144,0.1005)),2)]
    rng=np.random.default_rng(20260913)
    Z=np.zeros((60,7)); Z[:,:5]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(60,5))
    Z[:,5]=rng.uniform(0.3,1.8,60); Z[:,6]=np.clip(rng.normal(0.02,0.01,60),1e-3,0.1)
    Z[:,:5]=rg.integrate_ensemble(Z[:,:5],C,Km,Z[:,6],Z[:,5],heat,0.0,30.0,2.0)
    T=np.full(5,rg.T_AIR_K); _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,0.0,30.0,2.0); Tt=tr[-1]
    ws=[Wall[i,c] for i,c in sel]; u0=[um[i,c] for i,c in sel]
    Yf=np.column_stack([Z[:,2],Z[:,0]]+[b+(Z[:,:5]-mp)@w for w,b in zip(ws,u0)])
    yt=np.r_[Tt[2],Tt[0],[b+w@(Tt-mp) for w,b in zip(ws,u0)]]
    zb=Z.mean(0); Zi=zb+1.02*(Z-zb); yb=Yf.mean(0); Yi=yb+1.02*(Yf-yb)
    dZ=Zi-Zi.mean(0); dY=Yi-Yi.mean(0)
    K=(dZ.T@dY/59)@np.linalg.inv(dY.T@dY/59+np.eye(4)*0.09)
    dy=yt-Yf.mean(0)
    return K,dy,Cc[pod],[co[i]*1000 for i,_ in sel]


def fig_correction():
    K,dy,P,S=kalman_example()
    ct=K[:5,:2]@dy[:2]; c1=K[:5,2]*dy[2]; c2=K[:5,3]*dy[3]; tot=ct+c1+c2
    p=pv.Plotter(off_screen=True,window_size=(1100,1000)); p.set_background("white"); base_scene(p)
    for i,x in enumerate(P):
        col="#c0392b" if i in (0,2) else "#95a5a6"
        p.add_mesh(pv.Sphere(radius=3.4,center=x),color=col)
    p.add_point_labels(P+np.array([0,0,7]),[f"P{i}" for i in range(5)],font_size=26,text_color="#333333",
                       shape=None,show_points=False,always_visible=True,bold=True)
    for k,(x,dirv) in enumerate(zip(S,[(1,0,0),(0,0,1)])):
        p.add_mesh(pv.Cube(center=x,x_length=6,y_length=6,z_length=6),color="#1f9d62")
        p.add_mesh(pv.Arrow(start=x,direction=dirv,scale=16),color="#1f9d62")
        out=np.array([x[0],x[1],0.0]); out=out/np.linalg.norm(out)
        p.add_point_labels([x+14*out+np.array([0,0,4])],[f"S{k+1}"],font_size=28,text_color="#1f7a4c",
                           shape=None,show_points=False,always_visible=True,bold=True)
    img=shot(p,"corr")
    fig=plt.figure(figsize=(17,7.6))
    ax=fig.add_axes([0.0,0.02,0.46,0.86]); ax.imshow(img); ax.axis("off")
    ax.set_title("● 赤＝温度計（P2・P0、代表点）　● 灰＝温度計のない代表点\n■ 緑＝変位計 S1（Ux）・S2（Uz）＝代表点ではない上面の点",fontsize=13.5)
    ax2=fig.add_axes([0.52,0.13,0.46,0.66])
    x=np.arange(5); w=0.2
    ax2.bar(x-1.5*w,ct,w,color=RED,label="温度計2本から")
    ax2.bar(x-0.5*w,c1,w,color="#7fc8a0",label="変位計 S1 から")
    ax2.bar(x+0.5*w,c2,w,color=GREEN,label="変位計 S2 から")
    ax2.bar(x+1.5*w,tot,w,color=NAVY,label="合計")
    ax2.axhline(0,color="k",lw=1); ax2.set_xticks(x); ax2.set_xticklabels([f"P{i}" for i in range(5)])
    ax2.set_ylabel("直す量 [K]"); ax2.grid(axis="y",alpha=.3); ax2.legend(loc="upper left",fontsize=12,ncol=2)
    ax2.annotate(f"S2 のずれ {dy[3]:+.2f} µm が\nP4 を {c2[4]:+.1f} K 直す",xy=(4+0.5*w,c2[4]),xytext=(2.6,c2[4]+1.0),
                 fontsize=13.5,color=GREEN,weight="bold",arrowprops=dict(arrowstyle="->",color=GREEN,lw=2))
    ax2.set_title("1回目（t=30 秒）の補正：どのセンサが、どの代表点をどれだけ直したか",fontsize=14.5)
    fig.suptitle("代表点ではない変位計のずれが、カルマンゲインを通じて代表5点の温度を直す",fontsize=18,weight="bold")
    fig.text(0.52,0.02,"温度計だけだと P4 は +6 K 動くが、S2 が打ち消して、\n5点とも −4.3〜−5.0 K にそろう",fontsize=13,color=NAVY)
    fig.savefig(os.path.join(IMG,"explain_correction_by_sensor.png"),dpi=130,facecolor="white"); plt.close(fig)
    print("wrote explain_correction_by_sensor.png")
    print("check: dy",np.round(dy,3)," c2",np.round(c2,2)," tot",np.round(tot,2))


def fig_after_correction():
    """§3-4 段階4：補正の前後と、そこから計算した次の 30 秒（t=30 → 60 秒）"""
    from dacore.enkf import enkf_update
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],5); h=float(d["h"]); heat=int(d["heat_node"])
    U,mean,pod,Cc=load(); UPp=np.linalg.pinv(U[pod]); mp=mean[pod]
    op=np.load(os.path.join(RES,"dispop_allnodes.npz")); um=op["u_mean"]; Dall=op["D"]; co=op["coords"]
    Wall=np.einsum("nck,kj->ncj",Dall,UPp)
    near=lambda x:int(np.linalg.norm(co-np.array(x),axis=1).argmin())
    sel=[(near((0.0097,-0.0362,0.1005)),0),(near((-0.0346,0.0144,0.1005)),2)]
    ws=[Wall[i,c] for i,c in sel]; u0=[um[i,c] for i,c in sel]
    rng=np.random.default_rng(20260913)
    Z=np.zeros((60,7)); Z[:,:5]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(60,5))
    Z[:,5]=rng.uniform(0.3,1.8,60); Z[:,6]=np.clip(rng.normal(0.02,0.01,60),1e-3,0.1)
    Z[:,:5]=rg.integrate_ensemble(Z[:,:5],C,Km,Z[:,6],Z[:,5],heat,0,30,2.0)
    T=np.full(5,rg.T_AIR_K); Tt=[T]
    for a,b in [(0,30),(30,60)]:
        _,tr=rg.integrate_single(T,C,Km,h,1.0,heat,a,b,2.0); T=tr[-1]; Tt.append(T)
    Yf=np.column_stack([Z[:,2],Z[:,0]]+[b+(Z[:,:5]-mp)@w for w,b in zip(ws,u0)])
    y=np.r_[Tt[1][2],Tt[1][0],[b+w@(Tt[1]-mp) for w,b in zip(ws,u0)]]   # §3-4 と同じく観測ノイズなし
    Za=enkf_update(Z.copy(),y,None,np.diag([0.09]*4),np.random.default_rng(1),inflation=1.02,Yf=Yf)
    Za[:,5]=np.clip(Za[:,5],0,3); Za[:,6]=np.clip(Za[:,6],1e-4,0.2)
    Zn=Za.copy(); Zn[:,:5]=rg.integrate_ensemble(Za[:,:5],C,Km,Za[:,6],Za[:,5],heat,30,60,2.0)
    Zo=Z.copy();  Zo[:,:5]=rg.integrate_ensemble(Z[:,:5],C,Km,Z[:,6],Z[:,5],heat,30,60,2.0)
    zb=Z.mean(0); za=Za.mean(0)
    fig,axs=plt.subplots(1,2,figsize=(16,6.2))
    x=np.arange(5); lab=[f"P{i}" for i in range(5)]
    ax=axs[0]
    ax.plot(x,zb[:5]-273.15,"s",color=GRAY,ms=13,label="補正前（予報の平均）")
    ax.plot(x,za[:5]-273.15,"o",color=RED,ms=13,label="補正後（平均）")
    ax.plot(x,Tt[1]-273.15,"_",color="k",ms=34,mew=3.5,label="真値")
    for i in range(5): ax.annotate("",xy=(i,za[i]-273.15),xytext=(i,zb[i]-273.15),arrowprops=dict(arrowstyle="->",color=RED,lw=2))
    ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel("温度 [℃]"); ax.grid(alpha=.3)
    ax.set_title(f"t = 30 秒：補正で5点の温度が下がる\n発熱量 Q：{15*zb[5]:.1f} W → {15*za[5]:.1f} W（真値 15 W）",fontsize=15)
    ax.legend(loc="upper right",fontsize=12.5)
    ax=axs[1]
    ax.plot(x,Zo[:,:5].mean(0)-273.15,"s",color=GRAY,ms=13,label="補正しなかった場合の予報")
    ax.plot(x,Zn[:,:5].mean(0)-273.15,"o",color=RED,ms=13,label="補正後の値から予報")
    ax.plot(x,Tt[2]-273.15,"_",color="k",ms=34,mew=3.5,label="真値")
    ax.set_xticks(x); ax.set_xticklabels(lab); ax.set_ylabel("温度 [℃]"); ax.grid(alpha=.3)
    e1=np.abs(Zn[:,:5].mean(0)-Tt[2]).mean(); e0=np.abs(Zo[:,:5].mean(0)-Tt[2]).mean()
    ax.set_title(f"t = 60 秒：直した値から ROM で 30 秒進めた結果\n真値とのずれ：補正なし {e0:.2f} K → 補正あり {e1:.2f} K",fontsize=15)
    ax.legend(loc="upper right",fontsize=12.5)
    fig.suptitle("段階4：直した5点温度・Q・h を初期値にして、次の30秒を ROM で計算する（1回目、60メンバーの平均）",fontsize=16,weight="bold")
    fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig(os.path.join(IMG,"explain_after_correction.png"),dpi=130,facecolor="white"); plt.close(fig)
    print("wrote explain_after_correction.png")
    return zb,za,Tt,Zo[:,:5].mean(0),Zn[:,:5].mean(0)


if __name__=="__main__":
    fig_free_weights(); fig_weighted_sum(); fig_correction(); fig_after_correction()
