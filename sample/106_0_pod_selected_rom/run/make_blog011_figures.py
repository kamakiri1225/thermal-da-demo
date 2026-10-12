"""Figures for blog011; use saved calculations, no solvers rerun."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run'),str(ROOT.parent/'102_1_frontistr_hollow_cylinder_thermal_expansion/python')]
from dacore import plots
from dacore.temperature_mapping import build_idw_mapping,map_temperature
from analytic_covariance_sensor_study import operators,A_XYZ,O_XYZ,EXCL
import matplotlib.pyplot as plt
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
import cylinder_mesh
RES=ROOT/'results/blog011';IMG=ROOT/'docs/img/blog011'
info=json.loads((RES/'study.json').read_text());z=np.load(RES/'covariances.npz');t=z['times']
cc,coords,G,H,w,tm,um=operators();sel=info['selections']['quad32'];sc=np.load(RES/'scores_quad32.npz');mc=np.load(RES/'scores_MC4000.npz')
fig,ax=plt.subplots(1,3,figsize=(19,5.8));fig.subplots_adjust(top=.77,bottom=.16,wspace=.32)
for key,col,ls in [('MC4000','#0072B2','--'),('quad32','#CC3311','-'),('tangent','#009E73',':')]:
 B=z[key];ax[0].plot(t,np.sqrt(np.einsum('i,tij,j->t',w,B,w)),ls,color=col,lw=2.5,label=key)
 for i in range(5):ax[1].plot(t,np.sqrt(B[:,i,i]),ls,color=col,alpha=.65)
ax[0].set_ylabel('観測前の反りの標準偏差 [µm]');ax[1].set_ylabel('代表5点それぞれの標準偏差 [K]')
for a in ax[:2]:a.set_xlabel('時刻 [s]');a.grid(alpha=.3)
Bq=z['quad32'][4];delta=z['MC4000'][4]-Bq;v=abs(delta).max();im=ax[2].imshow(delta,cmap='coolwarm',vmin=-v,vmax=v);ax[2].set_xticks(range(5),[f'P{i}' for i in range(5)]);ax[2].set_yticks(range(5),[f'P{i}' for i in range(5)]);ax[2].set_title('150秒の共分散差：MC4000 − quad32');fig.colorbar(im,ax=ax[2],label='共分散の差 [K²]')
fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.91),ncol=3);fig.suptitle('熱方程式から求めた共分散と、4,000ケースの共分散',fontsize=21);fig.savefig(IMG/'covariance_comparison.png',dpi=150);plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(14,5.8));fig.subplots_adjust(top=.80,bottom=.17,wspace=.3)
for stage,lab in [('temperature_1','温度1本目'),('temperature_2','温度2本目'),('displacement_1','変位1本目'),('displacement_2','変位2本目')]:
 mask=tm if stage.startswith('temperature') else um
 ax[0].scatter(mc[stage][mask][::30],sc[stage][mask][::30],s=8,alpha=.5,label=lab)
lo,hi=ax[0].get_xlim();ax[0].plot([lo,hi],[lo,hi],'k--',lw=1);ax[0].set_xlabel('4,000ケース法の $J_c$ [µm]');ax[0].set_ylabel('熱方程式＋数値積分の $J_c$ [µm]');ax[0].grid(alpha=.3)
keys=['MC4000','quadrature_32','tangent'];vals=[info['timing_s'][k] for k in keys];bars=ax[1].bar(['4,000ケース','数値積分32点＋両端','局所線形化'],vals,color=['#0072B2','#CC3311','#009E73'])
for b,v in zip(bars,vals):ax[1].text(b.get_x()+b.get_width()/2,v,f'{v:.3f} s',ha='center',va='bottom')
ax[1].set_ylim(0,max(vals)*1.3);ax[1].set_ylabel('共分散を作る計算時間 [s]');ax[1].grid(axis='y',alpha=.3)
fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.3,.90),ncol=2);fig.suptitle('候補の採点値と共分散計算時間の比較',fontsize=21);fig.text(.5,.03,'時間＝30〜300秒の共分散作成のみ、1回の測定。候補探索・描画・EnKF・高忠実度解析は含まない。',ha='center');fig.savefig(IMG/'scores_and_cost.png',dpi=150);plt.close(fig)
# Mesh and all three displacement-direction maps, same score calculation as blog009.
M=build_idw_mapping(cc,coords);mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005);ids={n:i for i,(n,_) in enumerate(mesh['nodes'])};cells=[]
for _,conn in mesh['elements']:cells.extend([8,*[ids[n] for n in conn]])
grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh['elements']),12,np.uint8),coords*1000);pv.OFF_SCREEN=True
for stage in range(4):
 kind='temperature' if stage<2 else 'displacement';k=stage%2;key=f'{kind}_{k+1}';n=1 if stage<2 else 3
 values=map_temperature(M,sc[key])[:,None] if n==1 else sc[key].reshape(-1,3)
 valid=np.ones(len(coords),bool) if n==1 else um.reshape(-1,3).any(1)
 lim=(float(values[valid].min()),float(values[valid].max()))
 fig,axs=plt.subplots(1,n,figsize=(10 if n==1 else 24,10),squeeze=False);fig.subplots_adjust(top=.80,bottom=.15,left=.01,right=.97,wspace=.04)
 for c,ax in enumerate(axs[0]):
  gg=grid.copy();gg['J']=np.where(valid,values[:,c],np.nan)
  pl=pv.Plotter(off_screen=True,window_size=(850,950));pl.set_background('white');pl.add_mesh(gg,scalars='J',cmap='RdYlBu',clim=lim,show_scalar_bar=False,nan_color='#cccccc')
  records=[(f'T{i+1}',p,'#ff00cc') for i,p in enumerate(sel['temperature'][:min(stage+1,2)])]
  if stage>=2:records += [(f'S{i+1}',p,'#00b86b') for i,p in enumerate(sel['displacement'][:k+1])]
  for name,rec,col in records:
   if 'component' in rec and rec['component']==c:pl.add_mesh(pv.Arrow(start=np.array(rec['xyz_m'])*1000,direction=np.eye(3)[c],scale=25),color='black',lighting=False)
  pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.22);ax.imshow(pl.screenshot(return_img=True))
  for name,rec,col in records:
   if 'component' in rec and rec['component']!=c:continue
   coo=vtkCoordinate();coo.SetCoordinateSystemToWorld();coo.SetValue(*(np.array(rec['xyz_m'])*1000));x,y=coo.GetComputedDoubleDisplayValue(pl.renderer)
   ax.scatter(x,950-y,s=120,color=col,edgecolor='black',zorder=10);ax.annotate(name,(x,950-y),xytext=(8,10),textcoords='offset points',fontsize=18,bbox=dict(fc='white',ec='none',alpha=.85))
  pl.close();ax.axis('off');ax.set_title('温度を測る候補' if n==1 else ['Ux','Uy','Uz'][c],fontsize=23)
  if n==1:grid.point_data[key+'_um']=values[:,0]
  else:grid.point_data[key+'_'+['Ux','Uy','Uz'][c]+'_um']=np.where(valid,values[:,c],np.nan)
 fig.suptitle(['① 温度1本目','② 温度2本目（T1追加後）','③ 変位1本目（T1・T2追加後）','④ 変位2本目（T1・T2・S1追加後）'][stage],fontsize=25)
 fig.text(.5,.91,'熱方程式＋hの数値積分で求めた $J_c$ 分布／4,000ケース法と選択点が一致',ha='center',fontsize=17)
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*lim),cmap='RdYlBu'),ax=axs.ravel().tolist(),fraction=.025,pad=.015);cb.set_label('$J_c$ [µm]：小さいほど良い',fontsize=15)
 fig.text(.5,.045,'紫＝温度センサ、緑＝変位センサ、黒矢印＝測定方向。30〜300秒で採点。\n温度の色はOpenFOAMセル値をIDWで表示メッシュへ写したもの。変位の灰色領域はA/O近傍12 mmの除外範囲。',ha='center',fontsize=13)
 fig.savefig(IMG/f'{key}.png',dpi=130);plt.close(fig)
grid.save(RES/'vtk/analytic_sensor_scores.vtu')
# Clear history plots: identical MC/tangent curves verified, displayed once as quad32.
for case,title in [('learned','15 W・0〜300秒ON'),('q25','25 W・0〜300秒ON'),('intermittent','15 W・間欠加熱')]:
 d=np.load(RES/f'unobserved_traj_analytic_{case}.npz');tt=d['t']
 for arr in ['T','U','R']:
  assert np.array_equal(d[arr][3],d[arr][5]) and np.array_equal(d[arr][3],d[arr][7]),(case,arr)
  assert np.array_equal(d[arr][2],d[arr][4]) and np.array_equal(d[arr][2],d[arr][6]),(case,arr)
 fig,axs=plt.subplots(2,4,figsize=(23,11));fig.subplots_adjust(top=.80,bottom=.1,hspace=.38,wspace=.25)
 lines=[]
 for row in range(2):
  for k,ax in enumerate(axs[row]):
   a='T' if row==0 else 'U';j=k+3 if row==0 else k
   true=d[a+'true'][:,j]-(273.15 if row==0 else 0);ax.plot(tt,true,color='black',lw=5,alpha=.30,label='真値：OpenFOAM＋FrontISTR')
   for idx,col,ls,lab in [(1,'#E69F00','--','手動温度 P2・P0'),(4,'#0072B2','--','方程式から選定：温度2点'),(5,'#CC3311','-','方程式から選定：温度2点＋変位2点')]:
    ax.plot(tt,d[a][idx,:,j]-(273.15 if row==0 else 0),color=col,ls=ls,lw=2.3,label=lab)
   for lo,hi in ([(0,150),(300,450)] if case=='intermittent' else [(0,300)]):ax.axvspan(lo,hi,color='#ffe6bf',alpha=.3)
   ax.set_title(str(d['tlabels' if row==0 else 'ulabels'][j]),fontsize=15);ax.set_xlabel('時刻 [s]');ax.set_ylabel('温度 [°C]' if row==0 else '変位 [µm]');ax.grid(alpha=.3)
 fig.legend(*axs[0,0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.92),ncol=2,fontsize=14);fig.suptitle('観測に使わない温度4点と変位：'+title,fontsize=24)
 fig.text(.5,.035,'60メンバーEnKF・30秒ごと20更新・5 seed平均。4,000ケース法・局所線形化も同じ配置となり、計算された時刻歴は一致（重複線は省略）。',ha='center',fontsize=13)
 fig.savefig(IMG/f'holdout_history_{case}.png',dpi=145);plt.close(fig)
# Multi-case bar plot, identical three selection methods displayed once.
fig,axs=plt.subplots(1,3,figsize=(19,6));fig.subplots_adjust(top=.78,bottom=.21,wspace=.3)
for ax,metric,title in zip(axs,['field_rmse_K','反り A−O','Uz(O)'],['全温度場RMSE [K]','反りA−OのMAE [µm]','Uz(O)のMAE [µm]']):
 for k,(nm,lab,col) in enumerate([('手動温度 P2・P0','手動温度2点','#E69F00'),('quad32：温度2点','選定温度2点','#0072B2'),('quad32：温度＋変位4本','＋選定変位2点','#CC3311')]):
  vals=[]
  for case in ['learned','q25','intermittent']:
   v=json.loads((RES/f'unobserved_check_analytic_{case}.json').read_text())['configs'][nm];vals.append(v[metric] if metric=='field_rmse_K' else v['displacement_um'][metric])
  ax.bar(np.arange(3)+(k-1)*.24,vals,width=.24,color=col,label=lab)
 ax.set_xticks(range(3),['15 W','25 W','15 W間欠']);ax.set_title(title,fontsize=18);ax.grid(axis='y',alpha=.3)
fig.legend(*axs[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.91),ncol=3,fontsize=14);fig.suptitle('方程式に基づいて選んだセンサで、未観測の推定精度を確認',fontsize=22)
fig.text(.5,.07,'30〜300秒・5 seed平均。温度場は時刻ごとの全セルRMSEの平均、変位は平均絶対誤差。\n3つの配置選定法は同じ点を選んだため、同化結果も同一。',ha='center',fontsize=13)
fig.savefig(IMG/'holdout_error_bars.png',dpi=150);plt.close(fig)
print('blog011 figures and exact trajectory consistency checks complete',flush=True)
