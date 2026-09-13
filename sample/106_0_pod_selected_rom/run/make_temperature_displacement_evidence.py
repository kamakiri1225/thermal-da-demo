"""同一の温度P2同化結果で、温度5点とA/O変位差の追従を比較する。"""
from pathlib import Path
import sys
import json
import numpy as np

ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))
from dacore import plots
import matplotlib.pyplot as plt


def main():
    res = ROOT / 'results'
    d = np.load(res / 'da_history.npz')
    pod = np.load(res / 'qdeim_points.npz')
    op = np.load(res / 'disp_operator.npz')
    time = d['time']
    Ta, Tt = d['analysis_T'], d['truth_T']
    inv = np.linalg.pinv(pod['pod_modes'][pod['cell_idx']].astype(float))

    def displacement(T):
        a = (T - pod['mean'][pod['cell_idx']]) @ inv.T
        u = op['uz_mean'] + a @ op['Dmode'].T
        return u[:, 0] - u[:, 1]

    da, dt = displacement(Ta), displacement(Tt)
    heat = (time > 0) & (time <= 300)
    metrics = {'temperature_time_rmse_by_point_K': np.sqrt(np.mean((Ta[heat]-Tt[heat])**2, axis=0)).tolist(),
               'difference_time_rmse_um': float(np.sqrt(np.mean((da[heat]-dt[heat])**2))),
               'time_s': time.tolist(), 'difference_truth_um': dt.tolist(), 'difference_estimate_um': da.tolist()}
    xyz = pod['xyz'] * 1000
    locations = ['ヒータ側・内周近く／中間高さ', '+Y側の外周近く／中間高さ',
                 'ヒータ側・外周近く／中間高さ', 'ヒータ側・底面近く', '反ヒータ側・底面近く']
    fig = plt.figure(figsize=(15, 17), layout='constrained')
    layout = fig.add_gridspec(4, 2, height_ratios=[1.3, 1, 1, 1])
    top, side = fig.add_subplot(layout[0, 0]), fig.add_subplot(layout[0, 1])
    from matplotlib.patches import Circle, Rectangle
    for radius in [20, 37.5]:
        top.add_patch(Circle((0, 0), radius, fill=False, color='gray'))
    top.set(xlim=(-46, 55), ylim=(-46, 48), xlabel='x [mm]', ylabel='y [mm]')
    top.set_aspect('equal'); top.set_title('上から見た位置（高さは右図で確認）')
    side.add_patch(Rectangle((-37.5, 0), 75, 100.5, fill=False, color='gray'))
    # 中空円筒の内面（内径40 mm）は、くり抜き部であることが分かるよう点線で示す。
    side.add_patch(Rectangle((-20, 0), 40, 100.5, fill=False,
                             edgecolor='gray', linestyle='--', linewidth=1.4))
    side.text(0, 104, '中空部（内径40 mm）', ha='center', va='bottom',
              fontsize=9, color='dimgray')
    side.set(xlim=(-46, 56), ylim=(-12, 120), xlabel='x [mm]', ylabel='高さ z [mm]')
    side.set_title('横から見た位置（底面 z=0／上面 z=100.5 mm）')
    for i, (x,y,z) in enumerate(xyz):
        color = 'crimson' if i == 2 else 'tab:blue'
        marker = '*' if i == 2 else 'o'
        for ax, point in [(top,(x,y)), (side,(x,z))]:
            ax.scatter(*point, color=color, marker=marker, s=150 if i==2 else 65, zorder=3)
            ax.annotate(f'P{i}'+(' 温度観測' if i==2 else ''), point,
                        xytext=(4,9 if i!=0 else -15), textcoords='offset points', fontsize=9, color=color)
    for name,x in [('A',28),('O',-28)]:
        side.scatter(x,100.5,marker='s',color='purple',s=65)
        side.annotate(name+' 変位評価', (x,100.5),xytext=(-18,9),textcoords='offset points',fontsize=9,color='purple')
    top.text(39,-35,'ヒータ側\n＋X',color='darkred',fontsize=10)
    side.text(39,18,'ヒータ側\n＋X',color='darkred',fontsize=10)
    for ax in (top,side): ax.grid(alpha=.2)
    axes = np.array([[fig.add_subplot(layout[row, col]) for col in range(2)] for row in range(1,4)])
    for i, ax in enumerate(axes.flat):
        ax.axvspan(0, 300, color='orange', alpha=.08)
        if i < 5:
            ax.plot(time, Tt[:, i]-273.15, 'k-', lw=2.5, label='ROM真値')
            ax.plot(time, Ta[:, i]-273.15, '--o', color='tab:green', ms=3, label='温度P2のみ同化')
            x,y,z=xyz[i]
            ax.set_title(f'P{i}：{locations[i]}\n'+
                         f'({x:.1f}, {y:.1f}, {z:.1f}) mm／'+('観測あり' if i==2 else '観測なし・推定を評価')+
                         f"／RMSE {metrics['temperature_time_rmse_by_point_K'][i]:.2f} K", fontsize=10)
            ax.set_ylabel('温度 [℃]')
        else:
            ax.plot(time, dt, 'k-', lw=2.5, label='ROM真値→変位写像')
            ax.plot(time, da, '--o', color='tab:green', ms=3, label='同じ推定温度→変位写像')
            ax.set_title('重要：2点間変位差 Uz(A)−Uz(O)\n'+
                         f"加熱期RMSE {metrics['difference_time_rmse_um']:.2f} µm", color='darkred')
            ax.set_ylabel('2点間変位差 [µm]')
        ax.grid(alpha=.3)
        ax.legend(fontsize=8)
        ax.set_xlabel('時刻 [s]')
    fig.suptitle('温度と2点間変位差は、同じように追従しているか？\n'
                 '106のROM双子実験：温度P2だけを同化。P0/P1/P3/P4は未観測、A/O変位も未観測\n'
                 '60メンバー・seed=20260913・30秒ごと20回更新／真値Q=15 W、h=0.015351 W/K\n'
                 '黒＝ROM真値、緑＝同化後平均。橙背景＝加熱0〜300秒、以降は冷却。RMSEは30〜300秒。', fontsize=12)
    # Windows側の大文字パスは実行環境によって読み取り専用に見えるため、
    # 成果物はワークスペースの書き込み可能な表記で保存する。
    out_root = Path('/mnt/d/work/002_cae/openfoam/20260505_datadoka/sample/106_0_pod_selected_rom')
    fig.savefig(out_root/'docs/img/temperature_displacement_evidence.png', dpi=140)
    plt.close(fig)
    (out_root/'results/temperature_displacement_evidence.json').write_text(json.dumps(metrics, ensure_ascii=False, indent=2)+'\n')
    print(metrics)


if __name__ == '__main__':
    main()
