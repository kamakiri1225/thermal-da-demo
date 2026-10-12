"""Plot the four selected configurations from cached EnKF evaluation metrics."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from dacore import plots
import matplotlib.pyplot as plt
r=json.loads((ROOT/'results/unobserved_check_seq_learned.json').read_text())['configs']
keys=['手動の温度2点 P2・P0','選定した温度2点','選定温度2点＋底面変位2点','選定温度2点＋選定変位2点']
labels=['手動温度2点\nP2・P0','選定温度2点\nT1・T2','T1・T2\n＋底面変位2点','T1・T2\n＋選定変位2点']
colors=['#E69F00','#0072B2','#009E73','#CC3311']
fig,axs=plt.subplots(2,2,figsize=(12,9))
for ax,quantity,title,unit in zip(axs.flat,['field_rmse_K','反り A−O','Uz(A)','Uz(O)'],['全温度場の誤差','未観測の反り Uz(A)−Uz(O)','未観測の Uz(A)','未観測の Uz(O)'],['K','µm','µm','µm']):
 vals=[r[k]['field_rmse_K'] if quantity=='field_rmse_K' else r[k]['displacement_um'][quantity] for k in keys]
 bars=ax.bar(np.arange(4),vals,color=colors,width=.65)
 ax.bar_label(bars,labels=[f'{v:.3f}' for v in vals],padding=5,fontsize=12)
 ax.set_xticks(np.arange(4),labels,fontsize=11);ax.set_ylim(0,max(vals)*1.22)
 ax.set_ylabel(('全セルRMSEの時間平均' if quantity=='field_rmse_K' else '平均絶対誤差')+f' [{unit}]',fontsize=12)
 ax.set_title(title,fontsize=15);ax.grid(axis='y',alpha=.25);ax.set_axisbelow(True)
fig.suptitle('15 W：温度2点を選び、変位2点を追加すると誤差はどう変わるか',fontsize=17,y=.985)
fig.text(.5,.935,'30〜300秒の10時刻・5試行の平均／60メンバー／正解＝OpenFOAM＋FrontISTR',ha='center',fontsize=12)
fig.text(.5,.025,'すべて小さいほど良い。反りはA・Oの差、Uz(A)・Uz(O)は各点単体。誤差棒は表示していません。',ha='center',fontsize=11)
fig.subplots_adjust(top=.86,bottom=.13,left=.08,right=.98,hspace=.48,wspace=.28)
fig.savefig(ROOT/'docs/img/sequential_error_bars_learned.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,5.7))
for ax,case,title in zip(axs,['q25','intermittent'],['25 W：0〜300秒ON','15 W：0〜150秒・300〜450秒ON']):
 data=json.loads((ROOT/f'results/unobserved_check_seq_{case}.json').read_text())['configs'];idx=[0,1,3];vals=[data[keys[i]]['displacement_um']['反り A−O'] for i in idx]
 bars=ax.bar(np.arange(3),vals,color=[colors[i] for i in idx],width=.65);ax.bar_label(bars,labels=[f'{v:.3f}' for v in vals],padding=5,fontsize=13)
 ax.set_xticks(np.arange(3),[labels[i] for i in idx],fontsize=11);ax.set_ylim(0,max(vals)*1.23);ax.set_ylabel('未観測反りの平均絶対誤差 [µm]');ax.set_title(title,fontsize=14);ax.grid(axis='y',alpha=.25);ax.set_axisbelow(True)
fig.suptitle('加熱条件を変えた場合も、温度選定 → 変位追加で反りが改善した',fontsize=16,y=.99)
fig.text(.5,.9,'誤差の評価区間は両条件とも30〜300秒／60メンバー・5試行平均',ha='center',fontsize=12)
fig.subplots_adjust(top=.78,bottom=.19,left=.08,right=.98,wspace=.3)
fig.savefig(ROOT/'docs/img/sequential_error_bars_other_conditions.png',dpi=150);plt.close(fig)
