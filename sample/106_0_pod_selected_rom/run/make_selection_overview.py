"""Four successive selection stages, from cached J scores (no DA rerun)."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
sys.path.insert(0,str(ROOT.parent/'102_1_frontistr_hollow_cylinder_thermal_expansion/python'))
from dacore import plots
from dacore.temperature_mapping import build_idw_mapping,map_temperature
import cylinder_mesh,pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
import matplotlib.pyplot as plt
RES=ROOT/'results';IMG=ROOT/'docs/img'
kv=np.load(RES/'qdeim_points.npz');op=np.load(RES/'dispop_allnodes.npz');scores=np.load(RES/'sequential_scores.npz');sel=json.loads((RES/'sequential_selection.json').read_text())
co=op['coords'];M=build_idw_mapping(kv['cell_centres'],co)
mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005);ids={n:i for i,(n,_) in enumerate(mesh['nodes'])};cells=[]
for _,conn in mesh['elements']:cells.extend([8,*[ids[n] for n in conn]])
grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh['elements']),12,np.uint8),co*1000)
pv.OFF_SCREEN=True
fig,axs=plt.subplots(1,4,figsize=(26,12));fig.subplots_adjust(left=.015,right=.985,top=.81,bottom=.18,wspace=.13)
for stage,ax in enumerate(axs):
 kind='temperature' if stage<2 else 'displacement';k=stage%2
 rec=sel[kind][k];comp=rec.get('component',0)
 values=map_temperature(M,scores[kind][k]) if stage<2 else scores[kind][k].reshape(-1,3)[:,comp]
 eligible=np.ones(len(co),bool) if stage<2 else (np.linalg.norm(co-[.028,0,.1005],axis=1)>.012)&(np.linalg.norm(co-[-.028,0,.1005],axis=1)>.012)
 # Same scale as the corresponding original stage (all 3 displacement directions).
 arr=values if stage<2 else scores[kind][k].reshape(-1,3)[eligible]
 lim=(float(np.min(arr)),float(np.max(arr)));gg=grid.copy();gg['J']=np.where(eligible,values,np.nan)
 pl=pv.Plotter(off_screen=True,window_size=(650,850));pl.set_background('white');pl.add_mesh(gg,scalars='J',cmap='RdYlBu',clim=lim,show_scalar_bar=False,nan_color='#cccccc')
 for name,rec0 in [(('T' if kind=='temperature' else 'S')+str(j+1),v) for j,v in enumerate(sel[kind][:k+1])]:
  pt=np.array(rec0['xyz_m'])*1000
  if stage>=2:pl.add_mesh(pv.Arrow(start=pt,direction=np.eye(3)[rec0['component']],scale=28,tip_length=.32,tip_radius=.13,shaft_radius=.055),color='black',lighting=False)
 pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.25);ax.imshow(pl.screenshot(return_img=True))
 # Overlay the projected markers so that selected points on the back remain visible.
 for j,rec0 in enumerate(sel[kind][:k+1]):
  coord=vtkCoordinate();coord.SetCoordinateSystemToWorld();coord.SetValue(*(np.array(rec0['xyz_m'])*1000))
  xx,yy=coord.GetComputedDoubleDisplayValue(pl.renderer)
  ax.scatter([xx],[850-yy],s=120,c='#ff00cc',edgecolors='black',linewidths=1.4,zorder=10)
  name=('T' if stage<2 else 'S')+str(j+1)
  if stage==1 and j==1:name+='（背面）'
  ax.annotate(name,(xx,850-yy),xytext=(8,12),textcoords='offset points',fontsize=20,weight='bold',bbox=dict(facecolor='white',alpha=.9,edgecolor='none'),zorder=11)
 pl.close();ax.axis('off')
 title=['① 温度T1を選ぶ','② 温度T2を選ぶ','③ 変位S1を選ぶ（Ux）','④ 変位S2を選ぶ（Ux）'][stage]
 condition=['まだ何も観測しない','T1を観測する条件','T1・T2を観測する条件','T1・T2・S1を観測する条件'][stage]
 ax.set_title(title+'\n'+condition,fontsize=21,pad=12)
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*lim),cmap='RdYlBu'),ax=ax,orientation='horizontal',fraction=.045,pad=.02);cb.set_label('$J_c$ [µm]：小さいほど良い',fontsize=18);cb.ax.tick_params(labelsize=16)
fig.suptitle('温度2点 → 変位2点：観測を1本ずつ加えて採点し直す',fontsize=30,y=.985)
fig.text(.5,.92,'目的：未観測の反り Uz(A)−Uz(O)／4,000ケース／30〜300秒の標準偏差を平均して採点',ha='center',fontsize=21)
fig.text(.5,.065,'紫球＝選択点、太い矢印＝測定方向。カラーバーは各段階で異なる。\n変位はUx・Uy・Uzの全候補から選び、ここでは選ばれた方向Uxの分布を表示。A/O近傍12 mmは除外。',ha='center',fontsize=19)
fig.savefig(IMG/'sequential_selection_overview.png',dpi=150);plt.close(fig)
