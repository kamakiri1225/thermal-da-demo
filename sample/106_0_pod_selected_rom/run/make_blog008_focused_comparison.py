"""blog008の比較図を、保存済み結果からB/C参考条件を除いて描く。再解析はしない。"""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dacore import plots as _p
import matplotlib.pyplot as plt
RES=ROOT/'results'; IMG=ROOT/'docs/img'

def main():
    fig,axs=plt.subplots(5,3,figsize=(16,21))
    for col,(ax,case,label) in enumerate(zip(axs[-1],['learned','q25','intermittent'],['15 W','25 W','間欠加熱'])):
        z=np.load(RES/f'unobserved_traj_s_{case}.npz')
        assert str(z['ulabels'][0]) == '反り A−O'
        assert str(z['ulabels'][1]) == 'Uz(A)' and str(z['ulabels'][2]) == 'Uz(O)'
        tt=z['t']; true=z['Utrue'][:,0]
        ax.plot(tt,true,color='black',lw=2.5,label='真値（OpenFOAM＋FrontISTR）')
        for i in [0,1,2,4]:
            pred=z['U'][i,:,0]
            ax.plot(tt,pred,color=str(z['colors'][i]),lw=1.8,label=str(z['names'][i]))
        ax.set_title(label);ax.set_xlabel('時刻 [s]');ax.grid(alpha=.25)
        ax.axvspan(0,300,color='#FDEBD0',alpha=.2)
        for row in range(4):
            a=axs[row,col]
            if row==0:
                for i in [0,1,2,4]:a.plot(tt,z['R'][i],color=str(z['colors'][i]),lw=1.8)
                a.set_ylabel('温度場RMSE [K]');a.set_title(label+'：全20,696セル')
            else:
                j=row-1
                a.plot(tt,z['Ttrue'][:,j]-273.15,color='black',lw=2.5)
                for i in [0,1,2,4]:a.plot(tt,z['T'][i,:,j]-273.15,color=str(z['colors'][i]),lw=1.8)
                a.set_ylabel('温度 [℃]');a.set_title('未観測温度 '+str(z['tlabels'][j]))
            a.set_xlabel('時刻 [s]');a.grid(alpha=.25)
        for a in axs[:,col]:
            # The intermittent heater is ON during 0–150 and 300–450 s.
            for patch in list(a.patches):patch.remove()
            a.axvspan(0,150 if case=='intermittent' else 300,color='#FDEBD0',alpha=.2)
            if case=='intermittent':a.axvspan(300,450,color='#FDEBD0',alpha=.2)

        f,aa=plt.subplots(3,2,figsize=(11,13.5))
        targets=[('温度場RMSE [K]',None,None),('未観測温度 '+str(z['tlabels'][0])+' [K]','T',0),
                 ('未観測の反り A−O [µm]','d',None),('未観測 Uz(A) [µm]','U',1),
                 ('未観測 Uz(O) [µm]','U',2),('未観測 '+str(z['ulabels'][-1])+' [µm]','U',-1)]
        for a,(title,key,j) in zip(aa.ravel(),targets):
            if key=='T':a.plot(tt,z['Ttrue'][:,j],color='black',lw=2.4)
            elif key=='U':a.plot(tt,z['Utrue'][:,j],color='black',lw=2.4)
            elif key=='d':a.plot(tt,true,color='black',lw=2.4)
            for i in [0,1,2,4]:
                v=z['R'][i] if key is None else z['T'][i,:,j] if key=='T' else z['U'][i,:,j] if key=='U' else z['U'][i,:,0]
                a.plot(tt,v,color=str(z['colors'][i]),lw=1.8)
            a.set_title(title);a.set_xlabel('時刻 [s]');a.grid(alpha=.25);a.axvspan(0,300,color='#FDEBD0',alpha=.2)
        handles,labels=ax.get_legend_handles_labels()
        f.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.96),ncol=2,fontsize=11)
        f.suptitle(label+'：良い配置・悪い配置の未観測量（IDW、5 seed平均）',y=.995,fontsize=16)
        f.subplots_adjust(top=.855,bottom=.055,left=.09,right=.98,hspace=.34,wspace=.25);f.savefig(IMG/f'blog008_unobserved_{case}.png',dpi=120);plt.close(f)
    for a in axs[-1]:a.set_ylabel('反り A−O [µm]')
    fig.suptitle('実ソルバ真値での主要配置比較：温度場 → 未観測温度3点 → 未観測の反り',fontsize=19,y=.995)
    fig.legend(*axs[-1,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.976),ncol=2,fontsize=12)
    fig.text(.5,.915,'温度P2・P0固定／60メンバー／30秒ごと20更新／5試行平均／正解＝OpenFOAM＋FrontISTR',ha='center',fontsize=13)
    fig.tight_layout(rect=(0,.01,1,.89));fig.savefig(IMG/'blog008_AO_focused.png',dpi=130);plt.close(fig)
    out=json.loads((RES/'disp_elsewhere_predict_AO.json').read_text())['configs']
    names=[n for n in out if 'B・C' not in n and '自身' not in n]
    fig,axs=plt.subplots(1,3,figsize=(17,5))
    for ax,key,title,unit in zip(axs,['AO_um','A_um','field_K'],['反り A−O','Uz(A)','全温度場'],['µm','µm','K']):
        values=[out[n][key] for n in names];ax.barh(range(len(names)),values,color='#4786b8')
        ax.set_yticks(range(len(names)),names if ax is axs[0] else ['']*len(names),fontsize=8)
        for i,v in enumerate(values):ax.text(v,i,f' {v:.3f}',va='center',fontsize=9)
        ax.invert_yaxis();ax.set_title(title);ax.set_xlabel('平均誤差 ['+unit+']');ax.grid(axis='x',alpha=.2)
    fig.suptitle('ROM双子実験：主要配置の比較（温度P2＋P0は共通、5 seed平均）')
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(IMG/'blog008_ROM_focused.png',dpi=130);plt.close(fig)
    method=json.loads((RES/'placement_method_comparison.json').read_text())['methods']
    names=[n for n in method if '自身' not in n]
    fig,axs=plt.subplots(1,3,figsize=(16,5))
    for ax,case,label in zip(axs,['learned','q25','intermittent'],['15 W','25 W','間欠加熱']):
        values=[method[n]['AO_um'][case] for n in names]
        ax.barh(range(len(names)),values,color=['#C0392B' if '目的指向' in n else '#638dae' for n in names])
        for i,v in enumerate(values):ax.text(v,i,f' {v:.3f}',va='center',fontsize=9)
        ax.set_yticks(range(len(names)),names if ax is axs[0] else ['']*len(names),fontsize=9)
        ax.invert_yaxis();ax.set_title(label);ax.set_xlabel('未観測A−Oの平均絶対誤差 [µm]');ax.grid(axis='x',alpha=.25)
    fig.suptitle('配置の選び方の比較：A・O自身を観測する参考条件を除外')
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(IMG/'blog008_methods_focused.png',dpi=130);plt.close(fig)
    points=[n for n in names if method[n]['xyz_mm']]
    fig,axs=plt.subplots(2,3,figsize=(13,8))
    theta=np.linspace(0,2*np.pi,180)
    for ax,n in zip(axs.ravel(),points):
        ax.plot(37.5*np.cos(theta),37.5*np.sin(theta),color='#71808d')
        ax.plot(20*np.cos(theta),20*np.sin(theta),color='#71808d')
        for xyz,comp in zip(method[n]['xyz_mm'],method[n]['comp']):
            ax.plot(xyz[0],xyz[1],'s',color='#2E8B57',ms=9)
            ax.annotate(comp,(xyz[0],xyz[1]),xytext=(5,6),textcoords='offset points',fontsize=10)
        ax.plot([28.75,-28.75],[0,0],'k*',ms=10)
        ax.set_title(n,fontsize=11);ax.set_aspect('equal');ax.set_xlim(-47,47);ax.set_ylim(-47,47);ax.set_xlabel('x [mm]');ax.set_ylabel('y [mm]')
    for ax in axs.ravel()[len(points):]:ax.axis('off')
    fig.suptitle('上から見た変位観測位置（■選定点、★未観測の評価点A・O）')
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(IMG/'blog008_method_points_focused.png',dpi=130);plt.close(fig)
    print('wrote focused figures from cached results')
if __name__=='__main__':main()
