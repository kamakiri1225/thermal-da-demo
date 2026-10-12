"""FrontISTR-derived static W rows with sequentially selected sensor positions."""
from pathlib import Path
import sys,json,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from dacore import plots
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.lines import Line2D
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
RES=ROOT/'results';IMG=ROOT/'docs/img';grid=pv.read(RES/'vtk/W_AO.vtu')
sel=json.loads((RES/'sequential_selection.json').read_text());pv.OFF_SCREEN=True
fields=['W_A_um_per_K','W_O_um_per_K','W_AminusO_um_per_K']
limit=max(float(np.max(np.abs(grid[f]))) for f in fields);norm=SymLogNorm(linthresh=3e-4,vmin=-limit,vmax=limit,base=10);cmap=plt.get_cmap('jet')
fig,axs=plt.subplots(1,3,figsize=(24,12));fig.subplots_adjust(left=.015,right=.94,top=.79,bottom=.12,wspace=.035)
markers=[]
for i,p in enumerate(sel['temperature']):markers.append((f'T{i+1}'+('（背面）' if i==1 else ''),np.array(p['xyz_m'])*1000,'#ff00cc',None))
for i,p in enumerate(sel['displacement']):markers.append((f'S{i+1} Ux',np.array(p['xyz_m'])*1000,'#00cc55',p['component']))
for lab,pt in [('A',[28.75,0,100.5]),('O',[-28.75,0,100.5])]:markers.append((lab,np.array(pt),'black',None))
for ax,field,title in zip(axs,fields,['Uz(A)の熱感度 $W_A$','Uz(O)の熱感度 $W_O$','反りの熱感度 $W_A-W_O$']):
 g=grid.copy();g['display_rgb']=(cmap(norm(np.asarray(g[field])))[:,:3]*255).astype(np.uint8)
 pl=pv.Plotter(off_screen=True,window_size=(850,1000));pl.set_background('white');pl.add_mesh(g,scalars='display_rgb',rgb=True,show_scalar_bar=False)
 for name,pt,color,comp in markers:
  if comp is not None:pl.add_mesh(pv.Arrow(start=pt,direction=np.eye(3)[comp],scale=27,tip_length=.32,tip_radius=.13,shaft_radius=.055),color='black',lighting=False)
 pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.23);ax.imshow(pl.screenshot(return_img=True))
 offsets={'A':(-40,22),'O':(10,22),'S1 Ux':(10,-25),'S2 Ux':(-55,28)}
 for name,pt,color,comp in markers:
  coord=vtkCoordinate();coord.SetCoordinateSystemToWorld();coord.SetValue(*pt);x,y=coord.GetComputedDoubleDisplayValue(pl.renderer)
  ax.scatter([x],[1000-y],s=100 if name in ['A','O'] else 150,c=color,edgecolors='white' if name in ['A','O'] else 'black',linewidths=1.3,zorder=10)
  ax.annotate(name,(x,1000-y),xytext=offsets.get(name,(8,12)),textcoords='offset points',fontsize=19,weight='bold',bbox=dict(facecolor='white',alpha=.94,edgecolor='none'),arrowprops=dict(arrowstyle='-',color='black',lw=.9),zorder=11)
 pl.close();ax.axis('off');ax.set_title(title,fontsize=25,pad=14)
cbax=fig.add_axes([.95,.25,.015,.47]);cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),cax=cbax,ticks=[-.05,-.01,-.001,0,.001,.01,.05]);cb.set_label('熱感度 [µm/K]（符号付き・対称対数）',fontsize=19);cb.ax.tick_params(labelsize=17)
fig.suptitle('FrontISTRの熱感度分布と、今回選んだ観測位置を重ねる',fontsize=30,y=.985)
fig.text(.5,.915,'色＝そのFEM節点の温度を1 K変えたとき、A・OのUzまたはその差がいくら変わるか',ha='center',fontsize=21)
fig.text(.5,.87,'紫＝温度T1・T2　緑＝変位S1・S2　太い黒矢印＝測定方向Ux　黒点＝評価位置A・O',ha='center',fontsize=20)
fig.text(.5,.045,'Wは時間不変。3図は共通の色範囲。ゼロ付近は線形、それ以外は対数で表示。\n背面の選択点も投影位置に重ねて表示。温度センサはOpenFOAMセル中心、変位センサはFEM節点。',ha='center',fontsize=17)
fig.savefig(IMG/'W_AO_selected_sensors_oblique.png',dpi=160,bbox_inches='tight');plt.close(fig)
# ParaView-native files: W grid (already in mm), sensor coordinates and direction vectors.
grid.save(RES/'vtk/W_AO_selected_sensors.vtu')
pts=pv.PolyData(np.array([m[1] for m in markers]));pts['sensor_code']=np.arange(len(markers));pts['direction']=np.array([np.eye(3)[m[3]] if m[3] is not None else np.zeros(3) for m in markers]);pts.save(RES/'vtk/W_AO_selected_sensor_points.vtp')
(RES/'W_AO_selected_overlay.json').write_text(json.dumps({'source':'vtk/W_AO.vtu','note':'A/O rows from FrontISTR dumped K,H; difference verified against native FrontISTR Wdiff','coordinate_unit':'mm','sensor_codes':[{'code':i,'name':m[0],'xyz_mm':m[1].tolist()} for i,m in enumerate(markers)],'color_scale':'symmetric log, linthresh=0.0003 um/K','W_range_um_per_K':[-limit,limit]},ensure_ascii=False,indent=2))
