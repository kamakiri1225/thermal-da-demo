"""W重ね図と同じ視点で、実モデルのヒータ領域を赤く表示する。"""
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'run'))
from dacore import plots
from make_blog001_figs import heater_shell
import matplotlib.pyplot as plt
import pyvista as pv
from vtkmodules.vtkRenderingCore import vtkCoordinate
pv.OFF_SCREEN=True
mesh=pv.read(ROOT/'results/vtk/W_AO.vtu')
heater=heater_shell();heater.points*=1000
pl=pv.Plotter(off_screen=True,window_size=(1100,1000));pl.set_background('white')
pl.add_mesh(mesh,color='#cbd6df',show_scalar_bar=False)
pl.add_mesh(heater,color='#ff2500',lighting=False)
pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.16)
fig,ax=plt.subplots(figsize=(12,10));ax.imshow(pl.screenshot(return_img=True));ax.axis('off')
for name,point,offset,color in [
 ('ヒータ（赤い領域）',[37.6,0,50],(25,0),'#c52000'),
 ('A：ヒータ側',[28.75,0,100.5],(30,30),'black'),
 ('O：反対側',[-28.75,0,100.5],(-90,25),'black')]:
 co=vtkCoordinate();co.SetCoordinateSystemToWorld();co.SetValue(*point);x,y=co.GetComputedDoubleDisplayValue(pl.renderer)
 ax.annotate(name,(x,1000-y),xytext=offset,textcoords='offset points',fontsize=21,color=color,weight='bold',bbox=dict(fc='white',ec='none',alpha=.94),arrowprops=dict(arrowstyle='->',color=color,lw=2))
pl.close()
fig.suptitle('ヒータの位置：熱感度図と同じ向きから見る',fontsize=25)
fig.text(.5,.055,'外周の＋x側、周方向約±80°、高さ24.75〜75.75 mmを加熱\n赤はヒータ領域を表し、温度や熱感度の色ではありません。',ha='center',fontsize=19)
fig.subplots_adjust(top=.94,bottom=.12)
fig.savefig(ROOT/'docs/img/heater_location_oblique.png',dpi=150,bbox_inches='tight');plt.close(fig)
