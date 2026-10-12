"""Temperature-first greedy selection, followed by displacement selection; IDW operators."""
from pathlib import Path
import sys,json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'run')]
from select_nonrep_temp_disp import priors,score,add_row,R_T,R_U,A_XYZ,O_XYZ,EXCL
from dacore.temperature_mapping import build_idw_mapping,map_temperature
from dacore import plots
import matplotlib.pyplot as plt
import pyvista as pv
sys.path.insert(0,str(ROOT.parent/'102_1_frontistr_hollow_cylinder_thermal_expansion/python'))
import cylinder_mesh
RES=ROOT/'results';IMG=ROOT/'docs/img'

def main():
    kv=np.load(RES/'qdeim_points.npz');co=np.load(RES/'dispop_allnodes.npz');coords=co['coords'];cc=kv['cell_centres'];pod=kv['cell_idx']
    inv=np.linalg.pinv(kv['pod_modes'][pod]);G=kv['pod_modes']@inv
    W=np.einsum('nck,kj->ncj',co['D'],inv);nn=lambda p:int(np.linalg.norm(coords-p,axis=1).argmin())
    w=W[nn(A_XYZ),2]-W[nn(O_XYZ),2];Ps=priors();out={'temperature':[],'displacement':[], 'criterion':'mean posterior SD of A-O at 30,...,300 s', 'prior_cases':4000}
    fields=[];validT=np.ones(len(cc),bool);validT[pod]=False
    validU=np.repeat((np.linalg.norm(coords-A_XYZ,axis=1)>EXCL)&(np.linalg.norm(coords-O_XYZ,axis=1)>EXCL),3)
    for kind,H,r,mask in [('temperature',G,R_T,validT),('displacement',W.reshape(-1,5),R_U,validU)]:
        for k in range(2):
            sc=score(Ps,w,H,r);idx=int(np.where(mask,sc,np.inf).argmin());mask[idx]=False
            xyz=cc[idx] if kind=='temperature' else coords[idx//3]
            rec={'index':idx,'xyz_m':xyz.tolist(),'J_um':float(sc[idx])}
            if kind=='displacement':rec['component']=idx%3
            out[kind].append(rec);fields.append((kind,k,sc.copy(),idx));Ps=[add_row(P,H[idx],r) for P in Ps]
            print(kind,k+1,rec,flush=True)
    (RES/'sequential_selection.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
    np.savez_compressed(RES/'sequential_scores.npz',temperature=np.array([f[2] for f in fields[:2]]),displacement=np.array([f[2] for f in fields[2:]]))
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005);ids={n:i for i,(n,_) in enumerate(mesh['nodes'])};cells=[]
    for _,conn in mesh['elements']:cells.extend([8,*[ids[n] for n in conn]])
    grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh['elements']),12,np.uint8),coords*1000)
    M=build_idw_mapping(cc,coords);pv.OFF_SCREEN=True
    for kind,k,sc,idx in fields:
        nc=1 if kind=='temperature' else 3
        fig,axes=plt.subplots(1,nc,figsize=(6*nc,8),squeeze=False)
        vals=map_temperature(M,sc)[:,None] if nc==1 else sc.reshape(-1,3)
        valid=np.ones(len(coords),bool) if nc==1 else validU.reshape(-1,3).any(1)
        lim=(float(vals[valid].min()),float(vals[valid].max()))
        for c,ax in enumerate(axes[0]):
            gg=grid.copy();gg['J']=np.where(valid,vals[:,c],np.nan)
            pl=pv.Plotter(off_screen=True,window_size=(800,900));pl.set_background('white')
            pl.add_mesh(gg,scalars='J',cmap='RdYlBu',clim=lim,nan_color='#cccccc',show_scalar_bar=False)
            for q,rec in enumerate(out[kind][:k+1]):
                if nc==3 and rec['component']!=c:continue
                pt=np.array(rec['xyz_m'])*1000
                pl.add_mesh(pv.Sphere(radius=4,center=pt),color='#ff00cc',lighting=False)
                if nc==3:
                    pl.add_mesh(pv.Arrow(start=pt,direction=np.eye(3)[c],scale=28,tip_length=.32,tip_radius=.13,shaft_radius=.055),color='black',lighting=False)
                    pl.add_mesh(pv.Arrow(start=pt,direction=np.eye(3)[c],scale=28,tip_length=.30,tip_radius=.10,shaft_radius=.04),color='#ff00cc',lighting=False)
                pl.add_point_labels([pt+[0,0,15]],[('T' if nc==1 else 'S')+str(q+1)],font_size=35,text_color='black',shape_color='white',shape_opacity=1,always_visible=True)
            pl.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];pl.camera.zoom(1.2)
            ax.imshow(pl.screenshot(return_img=True));pl.close();ax.axis('off');ax.set_title('温度' if nc==1 else ['Ux','Uy','Uz'][c],fontsize=22)
        fig.suptitle(f'測定候補：'+('温度' if nc==1 else '変位')+f'{k+1}本目の選択',fontsize=22)
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*lim),cmap='RdYlBu'),ax=axes.ravel().tolist(),fraction=.035,pad=.02);cb.set_label('$J_c$ [µm]：小さいほど良い',fontsize=14)
        fig.text(.5,.035,'紫の球＝選択点、矢印＝測定方向。'+('温度候補は全セル（代表5点を除外）。' if nc==1 else '選定温度T1・T2を固定。A/O近傍12 mmは除外。'),ha='center',fontsize=12)
        fig.savefig(IMG/f'sequential_{kind}_{k+1}.png',dpi=130);plt.close(fig)
    for kind,k,sc,idx in fields:
        if kind=='temperature':grid.point_data[f'J_temperature_{k+1}_um']=map_temperature(M,sc)
        else:
            for c,name in enumerate(['Ux','Uy','Uz']):grid.point_data[f'J_displacement_{k+1}_{name}_um']=sc.reshape(-1,3)[:,c]
    grid.save(RES/'vtk/sequential_sensor_scores.vtu')
if __name__=='__main__':main()
