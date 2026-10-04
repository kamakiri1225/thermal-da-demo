"""ROM そのもの（同化なし）が間欠運転を再現できるかを、OpenFOAM の温度と比べて示す.

同化は一切しない。発熱量とヒータの ON/OFF を正しく与えて ROM を 0→600 秒まで一気に積分し、
OpenFOAM の温度（results/limit_truth_intermittent.npz）と重ねる。
比較のため「0〜300 秒に一定と仮定した（決め打ちのままの）ROM」も描く。

出力: docs/img/rom_alone_intermittent.png, results/rom_alone_intermittent.json
再現: OMP_NUM_THREADS=4 python3 run/rom_alone_intermittent.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); ROOT=os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
from dacore import plots as _p
import matplotlib.pyplot as plt
from dacore import rom_general as rg
from improve_heating_schedule import integrate
RES=os.path.join(ROOT,"results"); IMG=os.path.join(ROOT,"docs","img")
NPT=5; DT=2.0; SIG_T=0.30
ON_TRUE=lambda t:(t<150.0)|((t>=300.0)&(t<450.0))     # 実際の加熱（0〜150秒、300〜450秒）
ON_FIXED=lambda t:(t<300.0)                            # 決め打ちのまま
COLS=["#1f77b4","#ff7f0e","#C0392B","#2ca02c","#9467bd"]


def main():
    d=np.load(os.path.join(RES,"rom_calibrated_pod.npz"))
    C=d["C"]; Km=rg.tri_to_matrix(d["K_upper"],NPT); h=float(d["h"]); heat=int(d["heat_node"])
    pod=np.load(os.path.join(RES,"qdeim_points.npz"))["cell_idx"]
    z=np.load(os.path.join(RES,"limit_truth_intermittent.npz")); t=z["times"]; T5=z["Tfield"][:,pod]
    def run(on):
        T=np.full((1,5),rg.T_AIR_K); tr=[T[0].copy()]
        for k in range(1,len(t)):
            T=integrate(T,C,Km,np.array([h]),np.array([1.0]),heat,on,t[k-1],t[k]); tr.append(T[0].copy())
        return np.array(tr)
    Rt=run(ON_TRUE); Rf=run(ON_FIXED)
    et=Rt-T5; ef=Rf-T5
    out=dict(note="同化なし。ROMを0→600秒まで積分し、OpenFOAMの5点温度と比較",
             correct_max_K=float(np.abs(et).max()),correct_mean_K=float(np.abs(et)[1:].mean()),
             fixed_max_K=float(np.abs(ef).max()),fixed_mean_K=float(np.abs(ef)[1:].mean()),
             max_rise_K=float(T5.max()-rg.T_AIR_K))
    json.dump(out,open(os.path.join(RES,"rom_alone_intermittent.json"),"w"),ensure_ascii=False,indent=1)
    print(out)
    fig=plt.figure(figsize=(16,8.6))
    gs=fig.add_gridspec(2,2,height_ratios=[1.35,1],hspace=0.30,wspace=0.20)
    # 上段：温度の時刻歴
    for j,(R,lab,sub) in enumerate([(Rt,"【条件を合わせた ROM】→ 合う","ON/OFF を OpenFOAM と同じ 0〜150 秒・300〜450 秒に"),
                                     (Rf,"【比較用：わざと条件をずらした ROM】→ 外れる","ON/OFF を与えず 0〜300 秒ずっと ON のまま")]):
        ax=fig.add_subplot(gs[0,j])
        for k in range(600+1):
            pass
        tt=np.linspace(0,600,1201)
        ax.fill_between(tt,0,40,where=[ON_TRUE(x) for x in tt],color="#FDEBD0",alpha=.7,label="実際にヒータON")
        if j==1:
            ax.fill_between(tt,0,40,where=[ON_FIXED(x) for x in tt],color="none",hatch="///",edgecolor="#999",
                            linewidth=0.0,label="ROMが仮定したON")
        for i in range(5):   # OpenFOAM＝太い薄い線、ROM＝細い破線（重なっても両方見えるように）
            ax.plot(t,T5[:,i]-273.15,"-",color=COLS[i],lw=8,alpha=.28,solid_capstyle="round")
            ax.plot(t,R[:,i]-273.15,"--",color=COLS[i],lw=1.8,label=f"P{i}")
        from matplotlib.lines import Line2D
        style=[Line2D([0],[0],color="#555",lw=8,alpha=.28,label="OpenFOAM（太い薄い線）"),
               Line2D([0],[0],color="#555",lw=1.8,ls="--",label="ROM（細い破線）")]
        ax.set_ylim(19,27); ax.set_ylabel("温度 [℃]"); ax.grid(alpha=.3)
        ax.set_title(f"{lab}\n{sub}",fontsize=13,color=("#1F9D62" if j==0 else "#C0392B"),weight="bold")
        ax.tick_params(labelbottom=False)
        if j==0:
            l1=ax.legend(fontsize=9,ncol=3,loc="upper left"); ax.add_artist(l1)
        ax.legend(handles=style,fontsize=10,loc="lower right",framealpha=.95)
    # 下段：誤差
    for j,(e,lab) in enumerate([(et,"差（条件を合わせた ROM）"),(ef,"差（わざと条件をずらした ROM）")]):
        ax=fig.add_subplot(gs[1,j])
        ax.axhspan(-SIG_T,SIG_T,color="#F6C85F",alpha=.35,zorder=0,label="温度計のノイズ ±0.3 K")
        for i in range(5): ax.plot(t,e[:,i],"-",color=COLS[i],lw=1.8)
        ax.axhline(0,color="k",lw=.8); ax.set_xlabel("時刻 [s]"); ax.set_ylabel("ROM − OpenFOAM [K]")
        ax.grid(alpha=.3); ax.set_ylim(-4.2,4.2)
        ax.text(0.98,0.93,f"最大 {np.abs(e).max():.3f} K",transform=ax.transAxes,ha="right",va="top",
                fontsize=11,bbox=dict(fc="white",ec="#bbb"))
        ax.set_title(lab,fontsize=11.5)
        if j==0: ax.legend(fontsize=9.5,loc="lower right")
    fig.suptitle("ROM そのもの（データ同化なし）は間欠運転を再現できるか ― 15 W 間欠加熱、OpenFOAM との比較",fontsize=13.5)
    fig.tight_layout(rect=(0,0,1,0.95)); fig.savefig(os.path.join(IMG,"rom_alone_intermittent.png"),dpi=150); plt.close(fig)
    print("wrote rom_alone_intermittent.png")


if __name__=="__main__": main()
