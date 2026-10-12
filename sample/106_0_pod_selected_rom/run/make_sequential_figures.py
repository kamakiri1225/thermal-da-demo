from pathlib import Path
import sys,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from dacore import plots
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
IMG=ROOT/'docs/img'
for case in ['learned','q25','intermittent']:
 d=np.load(ROOT/f'results/unobserved_traj_seq_{case}.npz');t=d['t'];fig,axs=plt.subplots(4,2,figsize=(12,18))
 for ax in axs.flat:ax.grid(alpha=.25);ax.set_xlabel('時刻 [s]');ax.axvspan(0,150 if case=='intermittent' else 300,color='orange',alpha=.10)
 if case=='intermittent':
  for ax in axs.flat:ax.axvspan(300,450,color='orange',alpha=.10)
 panels=[('T',0,'未観測温度 P1（横側・中段）','℃'),
         ('T',1,'未観測温度 P3（ヒータ側・下部）','℃'),
         ('T',2,'未観測温度 P4（反対側・下部）','℃'),
         ('U',0,'未観測の反り Uz(A)−Uz(O)','µm'),
         ('U',1,'未観測の Uz(A)（上面・ヒータ側）','µm'),
         ('U',2,'未観測の Uz(O)（上面・反対側）','µm'),
         ('U',4,'未観測の Uz(E)（ヒータ側・中段）','µm'),
         ('U',5,'未観測の Uz(F)（反対側・中段）','µm')]
 for ax,(key,j,title,unit) in zip(axs.flat,panels):
  if key!='R':ax.plot(t,d[key+'true'][:,j]-(273.15 if key=='T' else 0),color='black',lw=2.4,label='正解（OpenFOAM＋FrontISTR）')
  for i in range(5):
   v=d[key][i] if j is None else d[key][i,:,j]-(273.15 if key=='T' else 0)
   ax.plot(t,v,color=str(d['colors'][i]),ls='--' if i in [0,1,3] else '-',lw=1.8,label=str(d['names'][i]))
  ax.set_title(title,fontsize=15);ax.set_ylabel(unit)
 handles,labels=axs[0,1].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.965),ncol=2,fontsize=11)
 fig.suptitle({'learned':'15 W・0〜300秒加熱','q25':'25 W・0〜300秒加熱','intermittent':'15 W・間欠加熱'}[case]+'：温度配置→変位配置の順で改善を検証',fontsize=16)
 fig.subplots_adjust(top=.89,bottom=.04,hspace=.38,wspace=.25);fig.savefig(IMG/f'sequential_history_{case}.png',dpi=140);plt.close(fig)
fig,axs=plt.subplots(1,3,figsize=(12,4.5));axs[0].bar(['OI'],[.09],color='#0072B2');axs[0].set_ylabel('最終温度場RMSE [K]');axs[0].set_title('温度場（OI）');axs[0].text(0,.09,'0.09 K',ha='center',va='bottom');axs[0].set_ylim(0,.13)
axs[1].bar(['OI','EnKF'],[8.9,14.93],color=['#0072B2','#CC3311']);axs[1].axhline(15,color='black',ls='--',label='正解 15 W');axs[1].legend(loc='upper center');axs[1].set_ylabel('最終発熱量 [W]');axs[1].set_ylim(0,20)
for i,v in enumerate([8.9,14.93]):axs[1].text(i,v,f'{v:g} W',ha='center',va='bottom')
axs[2].bar(['EnKF'],[13],color='#CC3311');axs[2].set_ylim(0,17);axs[2].set_ylabel('一連の計算時間 [時間]');axs[2].text(0,13,'約13時間',ha='center',va='bottom');axs[2].set_title('104：5メンバー・10更新')
fig.suptitle('高忠実度のデータ同化：既存結果の要約（異なる計算条件の結果）',fontsize=14);fig.text(.5,.015,'対象現象は600秒。EnKFの時間は1更新分ではない。OIの計算時間・EnKFの温度RMSEは本図では比較しない。',ha='center',fontsize=10);fig.tight_layout(rect=[0,.05,1,.93]);fig.savefig(IMG/'oi_enkf_summary.png',dpi=150);plt.close(fig)
