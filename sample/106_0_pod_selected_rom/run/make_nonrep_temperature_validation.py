"""blog009/010: 同化に使わない代表点以外の4温度点を、保存済み同一EnKF試験から比較する。"""
from pathlib import Path
import sys,json,shutil
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run')]
from dacore import plots
from selected_pair_unobserved_check import FREE_T
from make_blog001_figs import heater_shell
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge,Circle
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
pv.OFF_SCREEN=True
kv=np.load(ROOT/'results/qdeim_points.npz');cc=kv['cell_centres'];pod=kv['cell_idx']
indices=[int(np.linalg.norm(cc-np.array(x),axis=1).argmin()) for _,x in FREE_T]
assert not set(indices)&set(pod)
coords=cc[indices]*1000

for number,sub in [('009',ROOT/'results'),('010',ROOT/'results/blog010')]:
 img=ROOT/'docs/img'/('blog010' if number=='010' else '');img.mkdir(exist_ok=True)
 sel=json.loads((sub/'sequential_selection.json').read_text());ts=[p['index'] for p in sel['temperature']]
 assert not set(indices)&set(ts), '検証点が観測点と重複しています'
 metadata={'note':'4 points outside ROM representative points and outside temperature observation points','points':[{'id':i+1,'label':FREE_T[i][0],'cell_index':idx,'xyz_mm':coords[i].tolist()} for i,idx in enumerate(indices)],'temperature_observation_cells':ts,'ROM_representative_cells':pod.tolist(),'cases':{}}
 # Grey model with heater and separate markers. All markers projected so back points remain visible.
 mesh=pv.read(ROOT/'results/vtk/W_AO.vtu');fig,axs=plt.subplots(1,2,figsize=(16,9))
 pl=pv.Plotter(off_screen=True,window_size=(1000,1100));pl.set_background('white');pl.add_mesh(mesh,color='#d4dce2')
 sh=heater_shell();sh.points*=1000;pl.add_mesh(sh,color='#ff3500',lighting=False)
 pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.18);axs[0].imshow(pl.screenshot(return_img=True))
 for i,pt in enumerate(coords):
  c=vtkCoordinate();c.SetCoordinateSystemToWorld();c.SetValue(*pt);x,y=c.GetComputedDoubleDisplayValue(pl.renderer)
  axs[0].scatter(x,1100-y,s=170,c='#0066ee',edgecolors='white',lw=1.2,zorder=10)
  axs[0].annotate(f'検証点{i+1}',(x,1100-y),xytext=(10,10),textcoords='offset points',fontsize=18,weight='bold',bbox=dict(fc='white',ec='none',alpha=.9),zorder=11)
 pl.close();axs[0].axis('off');axs[0].set_title('青：温度を確かめる未観測点',fontsize=22)
 ax=axs[1]
 for radius in [20,37.5]:ax.add_patch(Circle((0,0),radius,fill=False,color='#777'))
 ax.add_patch(Wedge((0,0),38,-80,80,width=2,fc='#ff3500',alpha=.9))
 ax.scatter(cc[pod,0]*1000,cc[pod,1]*1000,s=50,c='#777',marker='x',label='ROM代表点（灰）')
 ax.scatter(cc[ts,0]*1000,cc[ts,1]*1000,s=130,c='#dd00aa',marker='s',label='観測温度T1・T2（紫）')
 for i,pt in enumerate(coords):
  ax.scatter(pt[0],pt[1],s=140,c='#0066ee',edgecolors='white',zorder=8,label='未観測の検証点（青）' if i==0 else None)
  offset=[(8,10),(12,25),(8,8),(12,-35)][i]
  ax.annotate(f'検証点{i+1}',pt[:2],xytext=offset,textcoords='offset points',fontsize=17,weight='bold',bbox=dict(fc='white',ec='none',alpha=.85),zorder=9)
 ax.set_xlim(-57,57);ax.set_ylim(-50,55);ax.set_aspect('equal');ax.grid(alpha=.2);ax.set_xlabel('x [mm]',fontsize=17);ax.set_ylabel('y [mm]',fontsize=17);ax.set_title('上から見る：赤はヒータ側（＋x）',fontsize=22)
 ax.legend(loc='upper center',bbox_to_anchor=(.5,1.0),fontsize=13)
 fig.suptitle(f'blog{number}：代表点でも観測点でもない4点で温度を確認',fontsize=25)
 fig.text(.5,.035,'斜め図は背面の点も投影して表示。点の高さと正確なセル中心座標は記事の表に記載。',ha='center',fontsize=17)
 fig.tight_layout(rect=[0,.13,1,.94]);fig.savefig(img/'nonrep_temperature_positions.png',dpi=150);plt.close(fig)
 for case,title in [('learned','15 W・0〜300秒ON'),('q25','25 W・0〜300秒ON'),('intermittent','15 W・間欠加熱')]:
  d=np.load(sub/f'unobserved_traj_seq_{case}.npz');met=json.loads((sub/f'unobserved_check_seq_{case}.json').read_text())['configs'];t=d['t']
  fig,axs=plt.subplots(2,2,figsize=(16,12))
  for i,ax in enumerate(axs.flat):
   col=i+3;lab=str(d['tlabels'][col]);ax.plot(t,d['Ttrue'][:,col]-273.15,color='black',lw=2.8,label='正解（OpenFOAM）')
   for k in range(len(d['names'])):
    ax.plot(t,d['T'][k,:,col]-273.15,color=str(d['colors'][k]),ls='--' if k in [0,1,3] else '-',lw=2,label=str(d['names'][k]))
   ax.axvspan(0,150 if case=='intermittent' else 300,color='orange',alpha=.10)
   if case=='intermittent':ax.axvspan(300,450,color='orange',alpha=.10)
   ax.set_title(f'検証点{i+1}：{FREE_T[i][0]}\n({coords[i,0]:.1f}, {coords[i,1]:.1f}, {coords[i,2]:.1f}) mm',fontsize=20)
   ax.set_xlabel('時刻 [s]',fontsize=16);ax.set_ylabel('温度 [℃]',fontsize=16);ax.tick_params(labelsize=14);ax.grid(alpha=.25)
  handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.947),ncol=2,fontsize=15)
  fig.suptitle(f'blog{number}・{title}：未観測・非代表点の温度は合うか',fontsize=25)
  fig.text(.5,.025,'60メンバー・30秒ごと更新・5seed平均。線は各試行の60メンバー平均をさらに5試行で平均。\n黒＝OpenFOAM正解。橙帯＝加熱中。4点とも観測に使用せず、gappy-PODで復元した温度。',ha='center',fontsize=16)
  fig.subplots_adjust(top=.79,bottom=.11,left=.07,right=.98,hspace=.40,wspace=.23);fig.savefig(img/f'nonrep_temperature_{case}.png',dpi=150);plt.close(fig)
  metadata['cases'][case]={name:[vals['temperature_K'][str(d['tlabels'][i+3])] for i in range(4)] for name,vals in met.items()}
 (sub/'nonrep_temperature_validation.json').write_text(json.dumps(metadata,ensure_ascii=False,indent=2))
 print('verified no overlap and plotted blog',number,flush=True)
