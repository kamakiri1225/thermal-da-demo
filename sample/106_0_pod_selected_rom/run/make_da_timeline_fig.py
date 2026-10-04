"""ROM代表点と観測位置の違い、観測から状態へ戻るEnKF更新を図解する。

再現: python3 run/make_da_timeline_fig.py
出力: docs/img/da_cycle_timeline.png / .svg
画像だけを更新し、ブログ本文や解析結果には触れない。
"""
from pathlib import Path
import sys
ROOT = Path(__file__).absolute().parent.parent
sys.path.insert(0, str(ROOT))
from dacore import plots as _p
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

NAVY = '#17324f'
BLUE = '#2463b4'
GREEN = '#19764f'
ORANGE = '#bd5b14'
RED = '#b73434'


def panel(ax, x, y, w, h, title, lines, color, fill):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=0.12',
                              facecolor=fill,edgecolor=color,linewidth=2))
    ax.text(x+.18,y+h-.22,title,ha='left',va='top',fontsize=18,color=color,weight='bold')
    for offset,text,size in lines:
        ax.text(x+.18,y+h-offset,text,ha='left',va='top',fontsize=size,color=NAVY,linespacing=1.5)


def arrow(ax, start, end, color=NAVY):
    ax.add_patch(FancyArrowPatch(start,end,arrowstyle='-|>',mutation_scale=20,
                               linewidth=2.2,color=color,shrinkA=0,shrinkB=0))


def main():
    fig,ax = plt.subplots(figsize=(20,12))
    fig.subplots_adjust(left=.02,right=.98,bottom=.02,top=.98)
    ax.set_xlim(0,20); ax.set_ylim(0,12); ax.axis('off')
    ax.text(.2,11.8,'代表点以外の観測を、どうROMの補正に使うか',fontsize=27,weight='bold',color=NAVY,va='top')
    ax.text(.2,11.13,'観測点の値をROMへ代入するのではなく、観測との差から「代表5点・Q・h」を更新する。',fontsize=18,color=NAVY)

    # 上段: 予報と補正を同じ物理時刻として区別する。
    ax.text(.2,10.55,'時間：ROMは2秒刻みで進める ／ 観測による補正は30秒ごと',fontsize=18,weight='bold',color=NAVY)
    xs=[.9,4.0,7.1,10.2,13.3,16.4,19.3]
    ax.plot([xs[0],xs[-1]],[9.88,9.88],color='#8696a7',lw=3)
    for x,label in zip(xs,['0秒','30秒','60秒','90秒','…','300秒','600秒']):
        ax.plot([x,x],[9.71,10.05],color=RED if label!='0秒' else NAVY,lw=2.5)
        ax.text(x,9.47,label,ha='center',fontsize=16,color=NAVY)
    ax.text(5.55,10.08,'補正 → 予報 → 補正',ha='center',fontsize=14,color=RED)
    ax.text(.9,9.02,'0〜300秒：加熱　／　300〜600秒：冷却　／　30, 60, …, 600秒で計20回補正',fontsize=16,color=NAVY)

    # センサ観測は計算とは別の入力。60秒の更新へ明示してつなぐ。
    panel(ax,10.4,7.70,4.35,.98,'60秒の観測値',[(.50,'温度センサ・変位センサの値',16)],RED,'#fff1ee')
    arrow(ax,(12.57,7.68),(12.57,6.91),RED)
    ax.text(13.02,7.23,'同じ60秒で比較',fontsize=14,color=RED,va='center')

    y=3.43; h=3.45; w=4.35
    panel(ax,.2,y,w,h,'① ROMで30 → 60秒を予報',[
        (.78,'30秒で補正した状態から出発',16),
        (1.37,'各メンバーのQ・hを使い、\n代表5点の温度を時間積分',17),
        (2.63,'60メンバーそれぞれで計算',15)],BLUE,'#edf4ff')
    panel(ax,5.3,y,w,h,'② 60秒のセンサ値を予測',[
        (.78,'温度：gappy-PODで\n観測セルの温度を求める',17),
        (1.81,'変位：FrontISTRで事前計算\nした応答から求める',17),
        (2.90,'観測位置は代表5点以外も可',15)],GREEN,'#eef8f2')
    panel(ax,10.4,y,w,h,'③ 60秒の状態をEnKFで補正',[
        (.78,'観測 − 予測 ＝ 観測との差',17),
        (1.42,'60メンバーの共分散から\nカルマンゲインKを求める',17),
        (2.49,'K × 観測との差で\n代表5点・Q・hを更新',17)],ORANGE,'#fff4e6')
    panel(ax,15.5,y,w,h,'④ 補正後から次の予報へ',[
        (.78,'60秒の補正後の状態を保存',16),
        (1.45,'代表5点温度 → 次の初期温度\nQ・h → 次の計算のパラメータ',16),
        (2.58,'ROMで60 → 90秒を予報\n90秒の観測で再び補正',17)],BLUE,'#edf4ff')
    for x in [4.57,9.67,14.77]:
        arrow(ax,(x,5.14),(x+.68,5.14))
    ax.text(.2,7.20,'1サイクルの例：30秒の補正後 → 60秒の補正 → 次の予報',fontsize=19,weight='bold',color=NAVY)

    # 一番の疑問: センサのずれが代表点に戻る経路を独立表示。
    panel(ax,.2,.98,19.65,2.0,'なぜ、代表点以外の観測で代表5点を直せるのか？',[
        (.58,'各メンバーで「代表点の温度・Q・h」と「センサ位置の予測値」を対応させ、共分散を計算する。',17),
        (1.02,'その関係から、各センサのずれを7変数の補正量へ変換する係数K（7行 × 観測数の列）を作る。',17),
        (1.47,'補正後の状態 ＝ 補正前の状態 ＋ K ×（観測 ＋ メンバーごとの観測摂動 − 予測）',17)],ORANGE,'#fff9ef')
    ax.text(.2,.52,'更新する状態：温度5個・q_scale・h_ROM の7変数（Q = 15 q_scale）。h_ROMは放熱コンダクタンス［W/K］。',fontsize=14,color=NAVY)
    ax.text(.2,.17,'注：全温度場や変位を直接更新するのではない。補正は同じ時刻で行い、物理時間を進めない。',fontsize=14,color=NAVY)
    img=ROOT/'docs'/'img'
    fig.savefig(img/'da_cycle_timeline.png',dpi=160,facecolor='white')
    fig.savefig(img/'da_cycle_timeline.svg',facecolor='white')
    plt.close(fig)
    print('Updated docs/img/da_cycle_timeline.png and .svg')


if __name__ == '__main__':
    main()
