"""Plot single-update variance reduction and sequential MAE from cached results.

No ROM or finite-element simulations are run. Source: results/delta_vs_da.json.
"""
from pathlib import Path
import json
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dacore import plots  # Configure the same Japanese font as other figures.
import matplotlib.pyplot as plt

def main():
    rows = json.loads((ROOT/'results/delta_vs_da.json').read_text())['rows']
    predicted = np.array([row['pred_pct'] for row in rows])
    actual = np.array([row['act_pct'] for row in rows])
    mae = np.array([row['mae_um'] for row in rows])
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for axis, values, title, ylabel in [
        (axes[0],actual,'単回の分散減少との整合','1回目の実測分散減少 [%]'),
        (axes[1],mae,'逐次更新後の誤差とは順位が異なる','加熱期の変位差MAE [µm]')]:
        for index, row in enumerate(rows):
            axis.scatter(predicted[index], values[index], color=plt.cm.tab10(index),
                         marker='o' if row['label'].startswith('温度') else 'D',
                         s=65, label=row['label'])
        axis.set_xlabel('最初の予報分布で計算したΔ / 事前分散 [%]')
        axis.set_ylabel(ylabel)
        axis.grid(alpha=.25)
        axis.set_title(title+'\nPearson r = '+
                       f'{np.corrcoef(predicted,values)[0,1]:+.3f}')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.91),
               ncol=4,frameon=False,fontsize=11)
    fig.suptitle('単回の分散減少と時間履歴の真値誤差は別の評価量',fontsize=14,y=.99)
    fig.text(.5,.015,'既存delta_vs_da.json：ROM双子実験／60メンバー／5seed／温度または変位1点の観測／対象＝Uz(A)−Uz(O)',
             ha='center',fontsize=10)
    fig.tight_layout(rect=[0,.05,1,.78])
    fig.savefig(ROOT/'docs/img/review_delta_vs_final_error.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':
    main()
