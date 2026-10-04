"""blog_005の条件付き感度配置の証明を、説明用の数値で図示する。
新しいCAE・同化実験ではない。python3 run/make_sensitivity_design_example.py
"""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).absolute().parent.parent
sys.path.insert(0,str(ROOT))
from dacore import plots as _p
import matplotlib.pyplot as plt

def main():
    fig,axs=plt.subplots(1,2,figsize=(13,5.4))
    sig=np.array([1.,2.,3.]); post=1/(1+sig**2)
    axs[0].bar(['感度1','感度2','感度3'],post,color=['#a6bddb','#3690c0','#045a8d'])
    for i,v in enumerate(post):axs[0].text(i,v+.02,f'{v:.2f}',ha='center',fontsize=13)
    axs[0].set_ylim(0,.65);axs[0].set_ylabel('Qの事後分散（説明用の規格化値）')
    axs[0].set_title('Qだけが未知：大きなQ感度ほど分散が小さい\n同じ事前分散1・観測ノイズ分散1・1時刻')
    A=np.array([[10.,10.],[10.,10.]]);B=np.array([[3.,0.],[0.,3.]])
    eig=[np.linalg.eigvalsh(x.T@x)[0] for x in [A,B]]
    axs[1].bar(['A：大感度だが比例','B：小感度でも独立'],eig,color=['#d95f0e','#238b45'])
    for i,v in enumerate(eig):axs[1].text(i,v+.25,f'{v:.0f}',ha='center',fontsize=13)
    axs[1].set_ylim(0,12);axs[1].set_ylabel('Q/hを区別する情報（最小固有値）', fontsize=13)
    axs[1].set_title('Q・hが未知：感度の大きさだけでは選べない\nAの感度行列=[10,10;10,10] ／ B=[3,0;0,3]')
    for ax in axs:
        ax.title.set_fontsize(14)
        ax.tick_params(labelsize=12)
    fig.suptitle('感度から何が言えるか：条件を分けて考える',fontsize=18)
    fig.text(.5,.02,'説明用の線形・ガウス模型。106のCAE結果やEnKF比較結果ではありません。',ha='center',fontsize=12)
    fig.tight_layout(rect=[0,.07,1,.91])
    fig.savefig(ROOT/'docs/img/sensitivity_design_proof.png',dpi=150)
    plt.close(fig)
if __name__=='__main__':main()
