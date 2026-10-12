"""blog010: rotated targets, cached scores/EnKF histories and FrontISTR static W rows."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run')]
from dacore import plots
from dacore.temperature_mapping import build_idw_mapping,map_temperature
from make_blog001_figs import heater_shell
from select_nonrep_temp_disp import add_row
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
import scipy.io as sio
from scipy.sparse.linalg import splu
RES=ROOT/'results/blog010';IMG=ROOT/'docs/img/blog010'
info=json.loads((RES/'conditions.json').read_text());sel=json.loads((RES/'sequential_selection.json').read_text())
kv=np.load(RES/'qdeim_points.npz');op=np.load(RES/'dispop_allnodes.npz');co=op['coords'];cc=kv['cell_centres']
A=np.array(info['A_actual_m'])*1000;O=np.array(info['O_actual_m'])*1000
base=pv.read(ROOT/'results/vtk/W_AO.vtu');pv.OFF_SCREEN=True
compnames=['Ux','Uy','Uz']

def projection(pl,pt):
 c=vtkCoordinate();c.SetCoordinateSystemToWorld();c.SetValue(*pt);x,y=c.GetComputedDoubleDisplayValue(pl.renderer);return x,1000-y

def render(ax,grid,markers,rgb=False,scalars=None,clim=None,cmap='RdYlBu',heater=False):
 pl=pv.Plotter(off_screen=True,window_size=(850,1000));pl.set_background('white')
 pl.add_mesh(grid,scalars=scalars,rgb=rgb,clim=clim,cmap=cmap,show_scalar_bar=False,nan_color='#cccccc',color='#cbd6df' if scalars is None else None)
 if heater:
  sh=heater_shell();sh.points*=1000;pl.add_mesh(sh,color='#ff2500',lighting=False)
 for name,pt,color,c in markers:
  if c is not None:pl.add_mesh(pv.Arrow(start=pt,direction=np.eye(3)[c],scale=25,tip_radius=.13,shaft_radius=.05),color='black',lighting=False)
 pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.20);ax.imshow(pl.screenshot(return_img=True))
 for name,pt,color,c in markers:
  x,y=projection(pl,pt);ax.scatter(x,y,s=130,c=color,edgecolors='black',zorder=10)
  offset={'A（新）':(10,25),'O（新）':(-75,20),'A（旧）':(15,-25),'O（旧）':(-65,20)}.get(name,(8,12))
  ax.annotate(name,(x,y),xytext=offset,textcoords='offset points',fontsize=17,weight='bold',bbox=dict(fc='white',ec='none',alpha=.93),zorder=11)
 pl.close();ax.axis('off')

# Location diagram, same camera as W.
markers=[('A（新）',A,'black',None),('O（新）',O,'black',None),('A（旧）',np.array([28.75,0,100.5]),'#888888',None),('O（旧）',np.array([-28.75,0,100.5]),'#888888',None)]
fig,ax=plt.subplots(figsize=(11,12));render(ax,base,markers,heater=True)
fig.suptitle('ヒータ固定：評価点A・Oだけをz軸まわり＋90°回転',fontsize=24)
fig.text(.5,.045,'赤＝ヒータ（＋x側の側面）／黒＝新A・O／灰＝blog009の旧A・O\n＋zから見て反時計回り。実際の評価節点：A=(0, +28.75, 100.5)、O=(0, −28.75, 100.5) mm',ha='center',fontsize=16)
fig.subplots_adjust(top=.94,bottom=.13);fig.savefig(IMG/'model_rotated.png',dpi=150,bbox_inches='tight');plt.close(fig)

# Four selection stages, no duplicated VTK/Matplotlib labels.
scores=np.load(RES/'sequential_scores.npz');M=build_idw_mapping(cc,co)
fig,axs=plt.subplots(1,4,figsize=(26,12));fig.subplots_adjust(left=.015,right=.985,top=.80,bottom=.18,wspace=.13)
for stage,ax in enumerate(axs):
 kind='temperature' if stage<2 else 'displacement';k=stage%2;rec=sel[kind][k];c=rec.get('component',0)
 vals=map_temperature(M,scores[kind][k]) if stage<2 else scores[kind][k].reshape(-1,3)[:,c]
 eligible=np.ones(len(co),bool) if stage<2 else ((np.linalg.norm(co-np.array(info['A_requested_m']),axis=1)>.012)&(np.linalg.norm(co-np.array(info['O_requested_m']),axis=1)>.012))
 arr=vals if stage<2 else scores[kind][k].reshape(-1,3)[eligible];lim=(float(arr.min()),float(arr.max()))
 g=base.copy();g['J']=np.where(eligible,vals,np.nan)
 marks=[(('T' if stage<2 else 'S')+str(j+1),np.array(p['xyz_m'])*1000,'#ff00cc',None if stage<2 else p['component']) for j,p in enumerate(sel[kind][:k+1])]
 render(ax,g,marks,scalars='J',clim=lim)
 title=[f'① 温度T1',f'② 温度T2',f'③ 変位S1：{compnames[c]}',f'④ 変位S2：{compnames[c]}'][stage]
 cond=['まだ観測なし','T1を使う条件','T1・T2を使う条件','T1・T2・S1を使う条件'][stage]
 ax.set_title(title+'\n'+cond,fontsize=21,pad=12)
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*lim),cmap='RdYlBu'),ax=ax,orientation='horizontal',fraction=.045,pad=.02);cb.set_label('$J_c$ [µm]：小さいほど良い',fontsize=18);cb.ax.tick_params(labelsize=16)
fig.suptitle('新A−Oを対象に、温度2点 → 変位2点を選び直す',fontsize=29,y=.985)
fig.text(.5,.92,'4,000ケース／30〜300秒の標準偏差の平均／各段階で共分散を更新して再採点',ha='center',fontsize=21)
fig.text(.5,.065,'紫点＝選択位置、矢印＝測定方向。各段階で色範囲は異なります。\nUx・Uy・Uzを全て比較し、選ばれた方向の分布を表示。新A・Oの近傍12 mmは変位候補から除外。',ha='center',fontsize=18)
fig.savefig(IMG/'selection_overview.png',dpi=150);plt.close(fig)

# W rows for the new target. Native Wdiff in blog009 belongs to OLD A/O: never reuse it.
case=ROOT.parent/'105_0_sensor_placement_sensitivity/openfoam/dumpw_case'
K=sio.mmread(case/'dump_matrix_1_0.mm').tocsc().astype(float);H=sio.mmread(case/'H_matrix.mtx').tocsc().astype(float)
E=np.zeros((K.shape[0],2));E[3*info['A_node_index']+2,0]=1;E[3*info['O_node_index']+2,1]=1
adj=splu(K).solve(E,trans='T');Wr=(adj.T@H)*1e6
residual=float(np.linalg.norm(K.T@adj-E)/np.linalg.norm(E))
g=base.copy()
for field,v in zip(['W_A_um_per_K','W_O_um_per_K','W_AminusO_um_per_K'],[Wr[0],Wr[1],Wr[0]-Wr[1]]):g[field]=v
(RES/'vtk').mkdir(exist_ok=True);g.save(RES/'vtk/W_AO_rotated.vtu')
limit=float(np.max(np.abs(Wr)));limit=max(limit,float(np.max(np.abs(Wr[0]-Wr[1]))));norm=SymLogNorm(linthresh=3e-4,vmin=-limit,vmax=limit,base=10);cmap=plt.get_cmap('jet')
marks=[('A（新）',A,'black',None),('O（新）',O,'black',None)]+[(f'T{i+1}',np.array(p['xyz_m'])*1000,'#ff00cc',None) for i,p in enumerate(sel['temperature'])]+[(f'S{i+1} {compnames[p["component"]]}',np.array(p['xyz_m'])*1000,'#00cc55',p['component']) for i,p in enumerate(sel['displacement'])]
fig,axs=plt.subplots(1,3,figsize=(24,12));fig.subplots_adjust(left=.015,right=.94,top=.79,bottom=.12,wspace=.035)
for ax,f,title in zip(axs,['W_A_um_per_K','W_O_um_per_K','W_AminusO_um_per_K'],['新Uz(A)への熱感度','新Uz(O)への熱感度','新反りA−Oへの熱感度']):
 gg=g.copy();gg['rgb']=(cmap(norm(gg[f]))[:,:3]*255).astype(np.uint8);render(ax,gg,marks,rgb=True,scalars='rgb');ax.set_title(title,fontsize=25,pad=14)
cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=fig.add_axes([.95,.25,.015,.47]));cb.set_label('熱感度 [µm/K]：符号付き・対称対数',fontsize=19);cb.ax.tick_params(labelsize=17)
fig.suptitle('評価点を90°回した場合のFrontISTR熱感度と選定位置',fontsize=29,y=.985)
fig.text(.5,.915,'色＝FEM節点を1 K温めたときの対象変位の変化／Wは時間不変',ha='center',fontsize=21)
fig.text(.5,.87,'紫＝温度T1・T2　緑＝変位S1・S2　黒矢印＝選定された測定方向Uy',ha='center',fontsize=20)
fig.text(.5,.045,'FrontISTRのダンプK・Hから新しい評価自由度のW行を再計算。旧A/OのWdiffは使用していません。\n背面の点も投影して表示。3図共通色範囲、ゼロ近傍のみ線形・その他は対数。',ha='center',fontsize=17)
fig.savefig(IMG/'W_selected_overlay.png',dpi=150,bbox_inches='tight');plt.close(fig)
points=pv.PolyData(np.array([m[1] for m in marks]));points['sensor_code']=np.arange(len(marks));points['direction']=np.array([np.eye(3)[m[3]] if m[3] is not None else np.zeros(3) for m in marks]);points.save(RES/'vtk/selected_points.vtp')
(RES/'W_metadata.json').write_text(json.dumps({'source_K':str(case/'dump_matrix_1_0.mm'),'source_H':str(case/'H_matrix.mtx'),'adjoint_relative_residual':residual,'W_limit_um_per_K':limit,'sensor_codes':[{'code':i,'name':m[0],'xyz_mm':m[1].tolist()} for i,m in enumerate(marks)]},ensure_ascii=False,indent=2))
print('W adjoint residual',residual,flush=True)

# Identical eight unobserved histories and the same metrics/configurations.
truth_amplitudes={}
for case in ['learned','q25','intermittent']:
 d=np.load(RES/f'unobserved_traj_seq_{case}.npz');t=d['t'];fig,axs=plt.subplots(4,2,figsize=(14,20))
 panels=[('T',0,'未観測温度 P1（横側・中段）','℃'),('T',1,'未観測温度 P3（ヒータ側・下部）','℃'),('T',2,'未観測温度 P4（反対側・下部）','℃'),('U',0,'未観測の新反り Uz(A)−Uz(O)','µm'),('U',1,'新Uz(A)：＋y側・上面','µm'),('U',2,'新Uz(O)：−y側・上面','µm'),('U',4,'未観測Uz(E)：ヒータ側・中段','µm'),('U',5,'未観測Uz(F)：反対側・中段','µm')]
 for ax,(key,j,title,unit) in zip(axs.flat,panels):
  ax.grid(alpha=.25);ax.set_xlabel('時刻 [s]');ax.axvspan(0,150 if case=='intermittent' else 300,color='orange',alpha=.10)
  if case=='intermittent':ax.axvspan(300,450,color='orange',alpha=.10)
  ax.plot(t,d[key+'true'][:,j]-(273.15 if key=='T' else 0),color='black',lw=2.4,label='正解（OpenFOAM＋FrontISTR）')
  for i in range(5):ax.plot(t,d[key][i,:,j]-(273.15 if key=='T' else 0),color=str(d['colors'][i]),ls='--' if i in [0,1,3] else '-',lw=1.8,label=str(d['names'][i]))
  ax.set_title(title,fontsize=17);ax.set_ylabel(unit)
 handles,labels=axs[0,1].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.965),ncol=2,fontsize=13)
 title={'learned':'15 W・0〜300秒ON','q25':'25 W・0〜300秒ON','intermittent':'15 W・0〜150秒と300〜450秒ON'}[case]
 fig.suptitle(title+'：評価点を90°回転した同条件比較',fontsize=20);fig.subplots_adjust(top=.89,bottom=.04,hspace=.40,wspace=.25);fig.savefig(IMG/f'history_{case}.png',dpi=140);plt.close(fig)
 truth_amplitudes[case]={'max_abs_truth_new_difference_um':float(np.abs(d['Utrue'][:,0]).max()),'peak_to_peak_truth_new_difference_um':float(np.ptp(d['Utrue'][:,0]))}
keys=['手動の温度2点 P2・P0','選定した温度2点','選定温度2点＋底面変位2点','選定温度2点＋選定変位2点'];labels=['手動温度2点\nP2・P0','選定温度2点\nT1・T2','T1・T2\n＋底面変位2点','T1・T2\n＋選定変位2点'];colors=['#E69F00','#0072B2','#009E73','#CC3311']
metric=json.loads((RES/'unobserved_check_seq_learned.json').read_text())['configs']
fig,axs=plt.subplots(2,2,figsize=(14,10))
for ax,key,title,unit in zip(axs.flat,['field_rmse_K','反り A−O','Uz(A)','Uz(O)'],['全温度場','新反りA−O','新Uz(A)','新Uz(O)'],['K','µm','µm','µm']):
 vals=[metric[k]['field_rmse_K'] if key=='field_rmse_K' else metric[k]['displacement_um'][key] for k in keys];bars=ax.bar(np.arange(4),vals,color=colors);ax.bar_label(bars,fmt='%.3f',padding=5,fontsize=14);ax.set_xticks(np.arange(4),labels,fontsize=12);ax.set_ylim(0,max(vals)*1.25);ax.set_title(title,fontsize=18);ax.set_ylabel(('RMSEの時間平均' if key=='field_rmse_K' else '平均絶対誤差')+f' [{unit}]');ax.grid(axis='y',alpha=.25)
fig.suptitle('15 W：新しい評価点での配置比較（30〜300秒・5seed平均）',fontsize=22);fig.tight_layout(rect=[0,0,1,.94]);fig.savefig(IMG/'error_bars_learned.png',dpi=150);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(14,6))
for ax,case,title in zip(axs,['q25','intermittent'],['25 W・0〜300秒ON','15 W・間欠加熱']):
 met=json.loads((RES/f'unobserved_check_seq_{case}.json').read_text())['configs'];ii=[0,1,3];vals=[met[keys[i]]['displacement_um']['反り A−O'] for i in ii];bars=ax.bar(range(3),vals,color=[colors[i] for i in ii]);ax.bar_label(bars,fmt='%.3f',padding=5,fontsize=15);ax.set_xticks(range(3),[labels[i] for i in ii],fontsize=13);ax.set_ylim(0,max(vals)*1.25);ax.set_title(title,fontsize=19);ax.set_ylabel('新反りA−Oの平均絶対誤差 [µm]');ax.grid(axis='y',alpha=.25)
fig.suptitle('加熱条件を変えた場合の比較（30〜300秒・5seed平均）',fontsize=22);fig.tight_layout(rect=[0,0,1,.92]);fig.savefig(IMG/'error_bars_other_conditions.png',dpi=150);plt.close(fig)

# t150 remaining covariance, same re-evaluation as blog009 §3-4.
B=np.load(RES/'prior_covariances.npz')['B'][4];inv=np.linalg.pinv(kv['pod_modes'][kv['cell_idx']]);G=kv['pod_modes']@inv;Wall=np.einsum('nck,kj->ncj',op['D'],inv);w=Wall[info['A_node_index'],2]-Wall[info['O_node_index'],2]
rows=[G[p['index']] for p in sel['temperature']]+[Wall[p['index']//3,p['component']] for p in sel['displacement']];stages=[];P=B.copy()
for k in range(5):
 if k:P=add_row(P,rows[k-1],.09)
 if k:
  HH=np.array(rows[:k]);AA=B@HH.T;direct=B-AA@np.linalg.solve(HH@AA+.09*np.eye(k),AA.T);assert np.allclose(P,direct,atol=1e-10)
 stages.append({'n_obs':k,'P_K2':P.tolist(),'sigma2_um2':float(w@P@w),'sigma_um':float(np.sqrt(max(w@P@w,0)))})
fig,axs=plt.subplots(1,2,figsize=(18,7));lab=['測る前','T1','T1・T2','T1・T2\n＋S1','T1・T2\n＋S1・S2'];vals=[st['sigma2_um2'] for st in stages]
bars=axs[0].bar(range(5),vals,color=['#888','#0072b2','#56b4e9','#e69f00','#cc3311'])
for k,st in enumerate(stages):axs[0].text(k,vals[k],f"分散 {vals[k]:.4f}\nσ={st['sigma_um']:.3f} µm",ha='center',va='bottom',fontsize=14)
axs[0].set_ylim(0,max(vals)*1.35);axs[0].set_ylabel('新反りA−Oの分散 [µm²]');axs[0].set_title('対象変位の残存分散',fontsize=19)
for j in range(5):axs[1].plot(range(5),[st['P_K2'][j][j] for st in stages],'-o',label=f'P{j}')
axs[1].set_ylabel('代表5点温度の残存分散 [K²]');axs[1].legend(loc='upper right',ncol=3);axs[1].set_title('温度誤差の共分散Pの対角成分',fontsize=19)
for ax in axs:ax.set_xticks(range(5),lab,fontsize=15);ax.grid(axis='y',alpha=.25)
fig.suptitle('t=150秒：新しいA−Oに対して残る不確かさ',fontsize=23);fig.text(.5,.035,'4,000ケース・温度0.3 K／変位0.3 µmの観測ノイズ。実測値は使わない線形推定の予測。',ha='center',fontsize=16);fig.tight_layout(rect=[0,.10,1,.93]);fig.savefig(IMG/'posterior_t150.png',dpi=150);plt.close(fig)
(RES/'posterior_t150.json').write_text(json.dumps({'t_s':150,'B_K2':B.tolist(),'stages':stages},ensure_ascii=False,indent=2))
(RES/'truth_amplitudes.json').write_text(json.dumps(truth_amplitudes,ensure_ascii=False,indent=2))
print('figures complete',truth_amplitudes,flush=True)
