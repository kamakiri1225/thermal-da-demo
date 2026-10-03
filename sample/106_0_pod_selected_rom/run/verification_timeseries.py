"""各検証について「観測点・温度の時刻歴・変位の時刻歴」をそろえて示す図を作る.

これまでは検証ごとに「誤差の棒グラフ」しか出しておらず、
  ・どこを観測して、どこを評価したのか
  ・温度はどう合ったのか
  ・変位（反り）はどう合ったのか
が分からなかった。4つの検証すべてで同じ形の図を作る。

  検証1（§4）同じ場所に温度計か変位計か        真値＝ROM の双子実験
  検証2（§5）別の場所の変位で A・O を当てる     真値＝ROM の双子実験
  検証3（§6）真値を OpenFOAM＋FrontISTR に      3条件
  検証4（§7）配置の決め方を既存法と比べる        真値＝OpenFOAM＋FrontISTR（15 W）

出力: results/verification_timeseries.json
      docs/img/ver{1..4}_detail.png（ブログ用・図の見出しつき）
      docs/img/ver{1..4}_slide.png（スライド用・見出しなしで大きく見せる）
再現: OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 python3 run/verification_timeseries.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE); SAMPLE=os.path.dirname(ROOT)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from dacore.enkf import enkf_update
from improve_heating_schedule import integrate as integ_on
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; IQ=5; IH=6; NAUG=7; DT=2.0; OBS_DT=30.0; T_END=600.0; N_ENS=60
SIG_T=0.30; SIG_U=0.30; INFL=1.02
SEEDS=[20260913,20260914,20260915,20260916,20260917]
TN=[2,0]; COMP=["Ux","Uy","Uz"]
A_XYZ=(0.028,0.0,0.1005); O_XYZ=(-0.028,0.0,0.1005)
SEL=[((0.0097,-0.0362,0.1005),0),((-0.0346,0.0144,0.1005),2)]
BC =[((0.028,0.0,0.075),2),((-0.028,0.0,0.075),2)]
LOW=[((0.0375,0.0,0.005),2),((-0.0375,0.0,0.005),2)]
SCHED={"learned":lambda t:(t<300.0),"q25":lambda t:(t<300.0),
       "intermittent":lambda t:(t<150.0)|((t>=300.0)&(t<450.0))}
CASES=[("learned","15 W"),("q25","25 W"),("intermittent","間欠加熱")]
GRAY="#9AA5B1"; RED="#C0392B"; GREEN="#1F9D62"; ORANGE="#E67E22"; BLUE="#2E6FD8"; PURPLE="#8E6FB0"
plt.rcParams.update({"font.size":15,"axes.titlesize":16,"axes.labelsize":15,
                     "xtick.labelsize":13,"ytick.labelsize":13,"legend.fontsize":12})


# ---------------- 共通：上面を真上から見た「観測点・評価点」の図 ----------------
def draw_points(ax,temp_pts,disp_pts,title,note=""):
    th=np.linspace(0,2*np.pi,241)
    ax.plot(37.5*np.cos(th),37.5*np.sin(th),color="#8696a7",lw=1.5)
    ax.plot(20.0*np.cos(th),20.0*np.sin(th),color="#8696a7",lw=1.5)
    hs=np.linspace(-np.radians(76),np.radians(76),80)
    ax.plot(38.8*np.cos(hs),38.8*np.sin(hs),color=RED,lw=6,solid_capstyle="butt")
    ax.text(46,0,"ヒータ側",color=RED,fontsize=12,ha="center",va="center",rotation=-90)
    ax.plot([28.7,-28.7],[0,0],"o",color="k",ms=12,zorder=6)
    ax.text(28.7,7,"A",fontsize=14,ha="center",fontweight="bold")
    ax.text(-28.7,7,"O",fontsize=14,ha="center",fontweight="bold")
    ax.text(0,43,"●黒＝評価点 A・O（測らない）",fontsize=12,ha="center",color="k")
    for (x,y),lab in temp_pts:
        ax.plot(x,y,"o",color=RED,ms=12,mec="w",mew=1.4,zorder=5)
        ax.text(x-(7 if x>25 else 0),y+6.5,lab,color=RED,fontsize=12,
                ha="right" if x>25 else "center",fontweight="bold")
    for (x,y),c,lab in disp_pts:
        ax.plot(x,y,"s",color=GREEN,ms=13,mec="w",mew=1.5,zorder=5)
        d={"Ux":(1,0),"Uy":(0,1),"Uz":(0,0)}[c]
        if d!=(0,0):
            ax.annotate("",xy=(x+13*d[0],y+13*d[1]),xytext=(x,y),
                        arrowprops=dict(arrowstyle="-|>",color=GREEN,lw=2.4))
        ax.text(x,y-8,f"{lab} {c}",color=GREEN,fontsize=12,ha="center",va="top",fontweight="bold")
    ax.set_aspect("equal"); ax.set_xlim(-58,58); ax.set_ylim(-50,50); ax.axis("off")
    ax.set_title(title+("\n"+note if note else ""),fontsize=14)


def graphs_only(name,cyc,truth,res,cols,aotrue,node,t_title,a_title):
    """スライド用：観測点の図を外し、温度と変位の2枚だけを大きく出す"""
    fig,axs=plt.subplots(1,2,figsize=(15.5,6.2))
    panel_T(axs[0],cyc,truth,res,cols,node,"P4（温度計なし）の温度 [K]")
    axs[0].set_title(t_title); axs[0].legend(loc="lower right",fontsize=12)
    panel_AO(axs[1],cyc,aotrue,res,cols)
    axs[1].set_title(a_title); axs[1].legend(loc="lower right",fontsize=12)
    fig.tight_layout(); fig.savefig(os.path.join(IMG,f"{name}_slide.png"),dpi=150); plt.close(fig)
    print(f"wrote {name}_slide.png（2枚組）")


def save2(fig,name,sup,rect=(0,0,1,0.93)):
    """ブログ用（図の見出しつき）とスライド用（見出しなし・余白最小）を両方書き出す"""
    fig.suptitle(sup,fontsize=17); fig.tight_layout(rect=rect)
    fig.savefig(os.path.join(IMG,f"{name}_detail.png"),dpi=150)
    plt.close(fig); print(f"wrote {name}_detail.png")


def load_common():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    kv=np.load(os.path.join(RES,"qdeim_points.npz"))
    U=kv["pod_modes"].astype(float); mean=kv["mean"].astype(float); pod=kv["cell_idx"]; Cc=kv["cell_centres"]
    UPp=np.linalg.pinv(U[pod,:]); Af=U@UPp
    op=np.load(os.path.join(RES,"dispop_allnodes.npz"))
    um=op["u_mean"]; Dall=op["D"]; coords=op["coords"]
    Wall=np.einsum("nck,kj->ncj",Dall,UPp)
    near=lambda x:int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    return dict(C=C,Km=Km,h=h,heat=heat,U=U,mean=mean,pod=pod,Cc=Cc,UPp=UPp,Af=Af,
                um=um,coords=coords,Wall=Wall,near=near,mp=mean[pod])


def rom_truth(cm):
    cyc=np.arange(OBS_DT,T_END+1e-9,OBS_DT)
    T=np.full(NPT,rg.T_AIR_K); tr=[T.copy()]
    for a,b in zip(np.r_[0,cyc[:-1]],cyc):
        _,x=rg.integrate_single(T,cm["C"],cm["Km"],cm["h"],1.0,cm["heat"],a,b,DT); T=x[-1]; tr.append(T.copy())
    return cyc,np.array(tr)


def enkf_traj(cm,obs_disp,obs_temp_cells,T5true,utrue,cyc,on,seed):
    """1 seed 分を回し、各サイクルの (5点温度の推定, 反りA−Oの推定) を返す"""
    Wall=cm["Wall"]; um=cm["um"]; mp=cm["mp"]; Af=cm["Af"]; mean=cm["mean"]
    iA,iO=cm["iA"],cm["iO"]
    ws=[Wall[i,c] for i,c in obs_disp]; u0=[um[i,c] for i,c in obs_disp]
    wt=[Af[c] for c in obs_temp_cells]; t0=[mean[c] for c in obs_temp_cells]
    rng=np.random.default_rng(seed); ro=np.random.default_rng(seed+7)
    Z=np.zeros((N_ENS,NAUG)); Z[:,:NPT]=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(N_ENS,NPT))
    Z[:,IQ]=rng.uniform(0.3,1.8,N_ENS); Z[:,IH]=np.clip(rng.normal(0.02,0.01,N_ENS),1e-3,0.1)
    Rd=np.diag([SIG_T**2]*(len(TN)+len(wt))+[SIG_U**2]*len(ws))
    wA=Wall[iA,2]; wO=Wall[iO,2]; aA=um[iA,2]; aO=um[iO,2]
    Tm=[]; ao=[]; tp=0.0
    for ci,tb in enumerate(cyc,1):
        Z=Z.copy()
        Z[:,:NPT]=(integ_on(Z[:,:NPT],cm["C"],cm["Km"],Z[:,IH],Z[:,IQ],cm["heat"],on,tp,tb) if on is not None
                   else rg.integrate_ensemble(Z[:,:NPT],cm["C"],cm["Km"],Z[:,IH],Z[:,IQ],cm["heat"],tp,tb,DT))
        tp=tb
        yv=list(T5true[ci][TN]); Yf=Z[:,TN]
        for w,b0 in zip(wt,t0):                              # 代表点でないセルの温度観測
            yv.append(b0+w@(T5true[ci]-mp)); Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
        for (i,c),w,b0 in zip(obs_disp,ws,u0):               # 変位観測
            yv.append(utrue[ci,i,c] if utrue is not None else b0+w@(T5true[ci]-mp))
            Yf=np.column_stack([Yf,b0+(Z[:,:NPT]-mp)@w])
        y=np.array(yv)+ro.normal(0,np.sqrt(np.diag(Rd)))
        Z=enkf_update(Z,y,None,Rd,rng,inflation=INFL,Yf=Yf)
        Z[:,IQ]=np.clip(Z[:,IQ],0,3); Z[:,IH]=np.clip(Z[:,IH],1e-4,0.2)
        m=Z[:,:NPT].mean(0); dm=m-mp
        Tm.append(m); ao.append((aA+wA@dm)-(aO+wO@dm))
    return np.array(Tm),np.array(ao)


def series(cm,cfgs,T5true,utrue,cyc,on):
    out={}
    for name,obs_d,obs_t in cfgs:
        r=[enkf_traj(cm,obs_d,obs_t,T5true,utrue,cyc,on,s) for s in SEEDS]
        out[name]=dict(T=np.mean([x[0] for x in r],axis=0),AO=np.mean([x[1] for x in r],axis=0))
    return out


def panel_T(ax,cyc,truth,res,cols,node,ylab):
    ax.plot(cyc,truth[1:,node],"o-",color="k",lw=2.4,ms=5,label="真値",zorder=5)
    for (nm,c) in cols: ax.plot(cyc,res[nm]["T"][:,node],"-",color=c,lw=2.2,label=nm)
    ax.axvspan(0,300,color="#FDEBD0",alpha=.4)
    ax.set_xlabel("時刻 [s]"); ax.set_ylabel(ylab); ax.grid(alpha=.3)


def panel_AO(ax,cyc,aotrue,res,cols):
    ax.plot(cyc,aotrue,"o-",color="k",lw=2.4,ms=5,label="真値",zorder=5)
    for (nm,c) in cols: ax.plot(cyc,res[nm]["AO"],"-",color=c,lw=2.2,label=nm)
    ax.axvspan(0,300,color="#FDEBD0",alpha=.4)
    ax.set_xlabel("時刻 [s]"); ax.set_ylabel("反り A−O [µm]"); ax.grid(alpha=.3)


def main():
    cm=load_common(); near=cm["near"]; coords=cm["coords"]; Cc=cm["Cc"]
    cm["iA"]=near(A_XYZ); cm["iO"]=near(O_XYZ)
    iA,iO=cm["iA"],cm["iO"]; Wall=cm["Wall"]; um=cm["um"]
    sel=[(near(x),c) for x,c in SEL]; bc=[(near(x),c) for x,c in BC]; low=[(near(x),c) for x,c in LOW]
    tcell=[int(np.linalg.norm(Cc-np.array(x),axis=1).argmin()) for x in (A_XYZ,O_XYZ)]
    xy=lambda i:(coords[i][0]*1000,coords[i][1]*1000)
    P=lambda k:(cm["Cc"][cm["pod"][k]][0]*1000,cm["Cc"][cm["pod"][k]][1]*1000)
    tpts=[(P(2),"P2"),(P(0),"P0")]
    cyc,Ttr=rom_truth(cm)
    aotrue_rom=np.array([ (um[iA,2]+Wall[iA,2]@(Ttr[c]-cm["mp"])) - (um[iO,2]+Wall[iO,2]@(Ttr[c]-cm["mp"]))
                          for c in range(1,len(Ttr))])
    NODE=4   # 温度計を置いていない代表点 P4（ヒータの反対側）
    store={}

    # ================= 検証1：同じ場所に温度計か変位計か =================
    cfg1=[("温度2点のみ",[],[]),
          ("＋A・O に温度計2本",[],tcell),
          ("＋A・O に変位計2本",[(iA,2),(iO,2)],[])]
    r1=series(cm,cfg1,Ttr,None,cyc,None)
    fig,axs=plt.subplots(1,3,figsize=(17.5,5.6))
    draw_points(axs[0],tpts,[(xy(iA),"Uz","A"),(xy(iO),"Uz","O")],
                "検証1：同じ場所で測る","●赤＝温度計（P2・P0）／■緑＝A・O に温度計 or 変位計")
    cols=[("温度2点のみ",GRAY),("＋A・O に温度計2本",ORANGE),("＋A・O に変位計2本",RED)]
    panel_T(axs[1],cyc,Ttr,r1,cols,NODE,"P4（温度計なし）の温度 [K]")
    axs[1].set_title("温度：温度計を置いていない P4"); axs[1].legend(loc="lower right")
    panel_AO(axs[2],cyc,aotrue_rom,r1,cols)
    axs[2].set_title("変位：反り A−O"); axs[2].legend(loc="lower right")
    save2(fig,"ver1","検証1　同じ場所に温度計を置くか、変位計を置くか（真値＝ROM の双子実験・5 seed 平均）")
    graphs_only("ver1",cyc,Ttr,r1,cols,aotrue_rom,NODE,"温度：温度計を置いていない P4","変位：反り A−O")
    store["ver1"]={k:{"T_P4":v["T"][:,NODE].tolist(),"AO":v["AO"].tolist()} for k,v in r1.items()}

    # ================= 検証2：別の場所の変位で A・O を当てる =================
    cfg2=[("変位なし（温度2点のみ）",[],[]),
          ("A・O 自身を測る（参考）",[(iA,2),(iO,2)],[]),
          ("A・O の真下（B・C）",bc,[]),
          ("選んだ別の2点（本研究）",sel,[]),
          ("底面近く（悪い例）",low,[])]
    r2=series(cm,cfg2,Ttr,None,cyc,None)
    fig,axs=plt.subplots(1,3,figsize=(17.5,5.6))
    draw_points(axs[0],tpts,[(xy(i),COMP[c],"選定") for i,c in sel],
                "検証2：A・O は一度も測らない","●赤＝温度計／■緑＝選んだ変位2点")
    cols=[("変位なし（温度2点のみ）",GRAY),("A・O 自身を測る（参考）",BLUE),
          ("A・O の真下（B・C）",ORANGE),("選んだ別の2点（本研究）",RED),("底面近く（悪い例）",PURPLE)]
    panel_T(axs[1],cyc,Ttr,r2,cols,NODE,"P4（温度計なし）の温度 [K]")
    axs[1].set_title("温度：温度計を置いていない P4"); axs[1].legend(loc="lower right",fontsize=11)
    panel_AO(axs[2],cyc,aotrue_rom,r2,cols)
    axs[2].set_title("変位：観測していない反り A−O"); axs[2].legend(loc="lower right",fontsize=11)
    save2(fig,"ver2","検証2　別の場所の変位2点で、一度も測らない A・O の反りを当てる（真値＝ROM の双子実験・5 seed 平均）")
    graphs_only("ver2",cyc,Ttr,r2,cols,aotrue_rom,NODE,"温度：温度計を置いていない P4","変位：観測していない反り A−O")
    store["ver2"]={k:{"T_P4":v["T"][:,NODE].tolist(),"AO":v["AO"].tolist()} for k,v in r2.items()}

    # ================= 検証3：真値を実ソルバにする（3条件）=================
    cfg3=[("変位なし（温度2点のみ）",[],[]),
          ("A・O 自身を測る（参考）",[(iA,2),(iO,2)],[]),
          ("選んだ別の2点（本研究）",sel,[])]
    cols3=[("変位なし（温度2点のみ）",GRAY),("A・O 自身を測る（参考）",BLUE),("選んだ別の2点（本研究）",RED)]
    fig,axs=plt.subplots(2,3,figsize=(18.0,8.2))
    figT,axT=plt.subplots(1,3,figsize=(17.0,5.6))   # スライド用（温度だけ3条件）
    figU,axU=plt.subplots(1,3,figsize=(17.0,5.6))   # スライド用（変位だけ3条件）
    store["ver3"]={}
    for j,(case,lab) in enumerate(CASES):
        z=np.load(os.path.join(RES,f"limit_truth_{case}.npz"))
        t=z["times"]; Tf=z["Tfield"]; T5=Tf[:,cm["pod"]]
        ut=np.load(os.path.join(RES,f"truth_disp_all_{case}.npz"))["u"]
        aot=ut[:,iA,2]-ut[:,iO,2]
        r=series(cm,cfg3,T5,ut,t[1:],SCHED[case])
        panel_T(axs[0,j],t[1:],T5,r,cols3,NODE,"P4 の温度 [K]" if j==0 else "")
        axs[0,j].set_title(f"{lab}　温度（P4：温度計なし）")
        if j==0: axs[0,j].legend(loc="lower right",fontsize=11)
        panel_AO(axs[1,j],t[1:],aot[1:],r,cols3)
        axs[1,j].set_title(f"{lab}　変位（観測していない反り A−O）")
        if j>0: axs[1,j].set_ylabel("")
        panel_T(axT[j],t[1:],T5,r,cols3,NODE,"P4 の温度 [K]" if j==0 else "")
        axT[j].set_title(lab)
        if j==0: axT[j].legend(loc="lower right",fontsize=12)
        panel_AO(axU[j],t[1:],aot[1:],r,cols3)
        axU[j].set_title(lab)
        if j>0: axU[j].set_ylabel("")
        if j==0: axU[j].legend(loc="upper left",fontsize=12)
        store["ver3"][case]={k:{"T_P4":v["T"][:,NODE].tolist(),"AO":v["AO"].tolist()} for k,v in r.items()}
    fig.suptitle("検証3　真値を OpenFOAM＋FrontISTR にして、3条件で確かめ直す（A・O は一度も観測しない・5 seed 平均）",fontsize=17)
    fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(os.path.join(IMG,"ver3_detail.png"),dpi=150)
    fig.suptitle(""); fig.tight_layout(rect=(0,0,1,1))
    fig.savefig(os.path.join(IMG,"ver3_slide.png"),dpi=150); plt.close(fig)
    print("wrote ver3_detail.png / ver3_slide.png")
    figT.tight_layout(); figT.savefig(os.path.join(IMG,"ver3_slide_T.png"),dpi=150); plt.close(figT)
    figU.tight_layout(); figU.savefig(os.path.join(IMG,"ver3_slide_U.png"),dpi=150); plt.close(figU)
    print("wrote ver3_slide_T.png / ver3_slide_U.png")

    # ================= 検証4：配置の決め方を比べる（実ソルバ 15 W）=================
    pm=json.load(open(os.path.join(RES,"placement_method_comparison.json"),encoding="utf-8"))
    al={"熱感度 |w| 最大（素朴）":"よく動く点に置く","熱感度の大きさ最大（素朴）":"よく動く点に置く"}
    mm={al.get(k,k):v for k,v in pm["methods"].items()}
    pick=lambda nm:[(near(np.array(p)/1000.0),COMP.index(c)) for p,c in zip(mm[nm]["xyz_mm"],mm[nm]["comp"])]
    cfg4=[("変位なし（温度2点のみ）",[],[]),
          ("よく動く点に置く",pick("よく動く点に置く"),[]),
          ("A最適（温度場を狙う）",pick("A最適（温度の trace 最小）"),[]),
          ("目的指向（反りを狙う・本研究）",pick("目的指向（本研究）"),[])]
    cols4=[("変位なし（温度2点のみ）",GRAY),("よく動く点に置く",PURPLE),
           ("A最適（温度場を狙う）",GREEN),("目的指向（反りを狙う・本研究）",RED)]
    z=np.load(os.path.join(RES,"limit_truth_learned.npz"))
    t=z["times"]; T5=z["Tfield"][:,cm["pod"]]
    ut=np.load(os.path.join(RES,"truth_disp_all_learned.npz"))["u"]
    aot=ut[:,iA,2]-ut[:,iO,2]
    r4=series(cm,cfg4,T5,ut,t[1:],SCHED["learned"])
    fig,axs=plt.subplots(1,3,figsize=(17.5,5.6))
    dp=[(xy(i),COMP[c],"") for i,c in pick("目的指向（本研究）")]
    draw_points(axs[0],tpts,dp,"検証4：決め方だけを変える","■緑＝目的指向が選んだ2点（他の決め方は別の場所）")
    panel_T(axs[1],t[1:],T5,r4,cols4,NODE,"P4（温度計なし）の温度 [K]")
    axs[1].set_title("温度：温度計を置いていない P4"); axs[1].legend(loc="lower right",fontsize=11)
    panel_AO(axs[2],t[1:],aot[1:],r4,cols4)
    axs[2].set_title("変位：観測していない反り A−O"); axs[2].legend(loc="lower right",fontsize=11)
    save2(fig,"ver4","検証4　置き場所の決め方を変えて比べる（真値＝OpenFOAM＋FrontISTR・15 W・5 seed 平均）")
    graphs_only("ver4",t[1:],T5,r4,cols4,aot[1:],NODE,"温度：温度計を置いていない P4","変位：観測していない反り A−O")
    store["ver4"]={k:{"T_P4":v["T"][:,NODE].tolist(),"AO":v["AO"].tolist()} for k,v in r4.items()}

    json.dump({"note":"各検証の時刻歴（5 seed 平均）。T_P4＝温度計を置いていない代表点P4の推定温度、AO＝反りA−Oの推定",
               "cycles_s":cyc.tolist(),"data":store},
              open(os.path.join(RES,"verification_timeseries.json"),"w"),ensure_ascii=False,indent=1)
    print("wrote results/verification_timeseries.json")


if __name__=="__main__": main()
