"""保存済みIDW比較結果を、主要8パネルと全14パネルに分けて描く。"""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from dacore import plots as _p
from dacore.temperature_mapping import MAPPING_METHOD
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
RES=ROOT/'results';IMG=ROOT/'docs/img'
SHORT=['同化なし','温度P2・P0のみ','温度P2・P0＋変位S1・S2','代表点以外の温度2点','代表点以外の温度2点＋変位2点']

def draw(case, setname="n"):
    ids=[0,1,2,4] if setname=="s" else list(range(5))
    z=np.load(RES/f'unobserved_traj_{setname}_{case}.npz');tt=z['t'];assert len(z['names'])==5
    handles=[Line2D([],[],color='black',lw=2.7,label='真値（OpenFOAM＋FrontISTR）')]
    handles += [Line2D([],[],color=str(z['colors'][i]),lw=2,label=(SHORT[i] if setname=="n" else str(z["names"][i]))) for i in ids]
    overview=[('全20,696セルの温度場RMSE','R',None),('反り Uz(A)−Uz(O)（未観測）','U',0),
        ('未観測温度 P1（横側・中段）','T',0),('未観測温度 P3（ヒータ側・下部）','T',1),
        ('未観測温度 P4（反対側・下部）','T',2),('未観測 Uz(A)（上面・ヒータ側）','U',1),
        ('未観測 Uz(O)（上面・反対側）','U',2),('未観測 Uz(E)（ヒータ側・中段）','U',4),
        ('未観測 Uz(F)（反対側・中段）','U',5)]
    allpanels=[]
    for j in range(7):
        allpanels.append((str(z['tlabels'][j]),'T',j))
        allpanels.append((str(z['ulabels'][j]),'U',j) if j<6 else ('全温度場のRMSE','R',None))
    for full,panels in [(False,overview),(True,allpanels)]:
        rows=(len(panels)+1)//2;fig,axs=plt.subplots(rows,2,figsize=(12,4*rows+1.2))
        for ax,(title,key,j) in zip(axs.ravel(),panels):
            unit='温度 [°C]' if key=='T' else '変位 [µm]' if key=='U' else 'RMSE [K]'
            if key!='R':
                true=z['Ttrue'][:,j]-273.15 if key=='T' else z['Utrue'][:,j]
                ax.plot(tt,true,color='black',lw=2.7,zorder=5)
            else:ax.axhline(0,color='black',lw=1)
            for i in ids:
                v=z['R'][i] if key=='R' else z[key][i,:,j]-(273.15 if key=='T' else 0)
                ax.plot(tt,v,color=str(z['colors'][i]),lw=2,ls='--' if i in [1,2,3] else '-')
            ax.set_title(title,fontsize=15,pad=8);ax.set_xlabel('時刻 [s]',fontsize=13);ax.set_ylabel(unit,fontsize=13)
            ax.tick_params(labelsize=12);ax.grid(alpha=.25);ax.set_xlim(30,600)
            if case=='intermittent':
                ax.axvspan(30,150,color='#FDEBD0',alpha=.3);ax.axvspan(300,450,color='#FDEBD0',alpha=.3)
            else:ax.axvspan(30,300,color='#FDEBD0',alpha=.3)
        for ax in axs.ravel()[len(panels):]:ax.axis('off')
        fig.suptitle({'learned':'15 W','q25':'25 W','intermittent':'15 W 間欠加熱'}[case]+('：代表点以外から選ぶ配置の比較' if setname=='n' else '：温度P2・P0固定、変位配置の比較'),y=.995,fontsize=18,weight='bold')
        fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.5,.973),ncol=2,fontsize=12,frameon=False)
        fig.subplots_adjust(top=1-2.25/(4*rows+1.2),bottom=.04,left=.09,right=.985,hspace=.40,wspace=.28)
        fig.text(.5,.008,'真値＝保存済みOpenFOAM温度＋IDWで写像して解き直したFrontISTR。EnKF 60メンバー、30秒ごと、5 seed平均。',ha='center',fontsize=11)
        dest=IMG/f'unobserved_{setname}_{case}{"_all" if full else ""}.png'
        fig.savefig(dest,dpi=140,facecolor='white')
        if setname=='s' and not full:fig.savefig(IMG/f'blog008_unobserved_{case}.png',dpi=140,facecolor='white')
        plt.close(fig)
if __name__=='__main__':
    for setname in ['n','s']:
        for case in ['learned','q25','intermittent']:draw(case,setname)
