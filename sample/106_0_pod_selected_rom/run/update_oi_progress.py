"""完了済みnpzと各runの解析完了ログからOIの進捗図を作る。"""
from pathlib import Path
import re
import json
import sys
import numpy as np
ROOT=Path(__file__).absolute().parents[1]
S104=ROOT.parent/'104_0_openfoam_frontistr_da_enkf'
sys.path.insert(0,str(ROOT))
from dacore import plots
import matplotlib.pyplot as plt


def main():
    entries=[]
    fig,axes=plt.subplots(1,2,figsize=(14,5.5))
    configs=[('', 'run1（σQ=6 W）','tab:blue'), ('_b','run2（σQ=2 W）','tab:green'),
             ('_free','同化なし（Q固定22 W）','tab:orange'),
             ('_sq0p5','σQ=0.5 W','tab:purple'),('_sq50','σQ=50 W','tab:red')]
    for suffix,label,color in configs:
        log=S104/f'openfoam/run_fem_oi{suffix}.log'
        if not log.exists():continue
        content=log.read_text()
        matches=re.findall(r'analysis t=([\d.]+)s RMSE\(field\)=([\d.]+)K Q=([\d.]+)',content)
        saved=S104/f'results/oi_fullsolver{suffix}_history.npz'
        complete=saved.exists() and '[oi] DONE.' in content
        if complete:
            d=np.load(saved);t=np.array(d['times']);q=np.array(d['q']);err=np.array(d['rmse_field'])
        elif matches:
            arr=np.array(matches,float);t=np.r_[0,arr[:,0]];q=np.r_[22,arr[:,2]];err=np.r_[5,arr[:,1]]
        else:
            entries.append(dict(run=label,status='解析更新の記録なし',latest_time_s=None))
            continue
        state='600秒完了' if complete else f'{t[-1]:g}秒までの途中結果'
        entries.append(dict(run=label,status=state,latest_time_s=float(t[-1]),Q_W=float(q[-1]),
                            absolute_relative_Q_error_pct=float(abs(q[-1]/15-1)*100),temperature_rmse_K=float(err[-1])))
        axes[0].plot(t,q,'o--',color=color,ms=4,label=f'{label}／{state}')
        axes[1].semilogy(t,err,'o--',color=color,ms=4,label=label)
    axes[0].axhline(15,color='black',lw=2,label='真値Q=15 W')
    axes[0].set_ylabel('発熱量パラメータ Q [W]');axes[1].set_ylabel('全20,696セルの温度RMSE [K]')
    for ax in axes:
        ax.axvspan(0,300,color='orange',alpha=.07);ax.set_xlabel('時刻 [s]');ax.grid(alpha=.3);ax.legend(fontsize=8)
    fig.suptitle('実ソルバOI：温度場とQの推定は別々に評価する\n'
                 'OpenFOAM＋FrontISTR／温度hot・cold＋上面2領域Uz／60秒ごと更新／hは推定対象外\n'
                 '途中結果は記録済み時刻で止める。300秒以降の実ヒータ入力は0 W。',fontsize=11)
    fig.tight_layout(rect=[0,0,1,.89]);fig.savefig(ROOT/'docs/img/oi_parameter_progress.png',dpi=140);plt.close(fig)
    (ROOT/'results/oi_parameter_progress.json').write_text(json.dumps(entries,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(entries,ensure_ascii=False,indent=2))


if __name__=='__main__':main()
