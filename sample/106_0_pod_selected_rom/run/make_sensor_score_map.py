"""IDW operator: J_c maps and greedy selected displacement sensors, oblique view.

J_c is mean posterior standard deviation of A-O at 30,...,300 s [um].
Temperature P2/P0 are fixed; exclude nodes within 12 mm of A/O.
"""
from pathlib import Path
import sys, json, csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(ROOT/'run'))
sys.path.insert(0,str(ROOT.parent/'102_1_frontistr_hollow_cylinder_thermal_expansion/python'))
from dacore import plots as _p
from dacore.temperature_mapping import mapping_cache_valid, MAPPING_METHOD
from select_nonrep_temp_disp import priors, add_row, score, R_T, R_U, A_XYZ, O_XYZ, EXCL, COMP
import cylinder_mesh
import matplotlib.pyplot as plt
import pyvista as pv
RES=ROOT/'results'; IMG=ROOT/'docs/img'


def main():
    assert mapping_cache_valid(RES/'dispop_allnodes.npz'), 'Rebuild IDW operator first'
    kv=np.load(RES/'qdeim_points.npz'); U=kv['pod_modes'].astype(float); pod=kv['cell_idx']
    UPp=np.linalg.pinv(U[pod]); op=np.load(RES/'dispop_allnodes.npz')
    coords=op['coords']; Wall=np.einsum('nck,kj->ncj',op['D'],UPp); Wf=Wall.reshape(-1,5)
    nn=lambda x:int(np.linalg.norm(coords-np.array(x),axis=1).argmin())
    iA,iO=nn(A_XYZ),nn(O_XYZ); w=Wall[iA,2]-Wall[iO,2]
    valid=(np.linalg.norm(coords-A_XYZ,axis=1)>EXCL)&(np.linalg.norm(coords-O_XYZ,axis=1)>EXCL)
    available=np.repeat(valid,3); Bs=priors()
    Ps=[]
    for B in Bs:
        P=B.copy()
        for j in [2,0]:P=add_row(P,np.eye(5)[j],R_T)
        Ps.append(P)
    before=[float(np.mean([np.sqrt(max(w@P@w,0)) for P in Ps]))]
    records=[]; maps=[]; selected=[]; decompositions=[]
    for step in range(2):
        values=score(Ps,w,Wf,R_U); masked=np.where(available,values,np.inf)
        pick=int(masked.argmin()); selected.append(pick); n,c=divmod(pick,3)
        maps.append(np.where(np.repeat(valid,3),values,np.nan).reshape(-1,3))
        records.append(dict(step=step+1,node_index=n,node_id=int(op['node_ids'][n]),
            xyz_mm=(coords[n]*1000).tolist(),component=COMP[c],J_c_um=float(values[pick]),
            J_before_um=before[-1],weights=Wf[pick].tolist()))
        # Exact score explanation at t=150 s (index 4), conditional on previous sensors.
        P=Ps[4]; h=Wf[pick]; vd=float(w@P@w); cov=float(w@P@h); vy=float(h@P@h)
        decompositions.append(dict(step=step+1,time_s=150,variance_before_um2=vd,
            covariance=cov,signal_variance=vy,noise_variance=R_U,
            variance_reduction_um2=cov**2/(vy+R_U),sigma_after_um=float(np.sqrt(max(vd-cov**2/(vy+R_U),0)))))
        Ps=[add_row(P,Wf[pick],R_U) for P in Ps]; available[pick]=False
        before.append(float(values[pick]))
        print(records[-1],flush=True)
    out=dict(mapping_method=MAPPING_METHOD,definition='J_c=mean sigma(A-O) at 30:30:300 s; smaller is better',
        fixed_temperature_cells=[int(x) for x in pod[[2,0]]],prior_cases=4000,prior_seed=20260930,
        temperature_noise_K=.3,displacement_noise_um=.3,exclusion_radius_mm=12,
        selected=records,explanation_t150=decompositions)
    (RES/'sensor_score_idw.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
    np.savez_compressed(RES/'sensor_score_idw.npz',mapping_method=MAPPING_METHOD,
        scores=np.array(maps),selected=np.array(selected),coords=coords,prior_covariances=np.array(Bs),w=w)
    with (RES/'sensor_score_idw.csv').open('w',newline='') as f:
        writer=csv.writer(f);writer.writerow(['node_id','x_mm','y_mm','z_mm','direction','eligible','J_step1_um','J_step2_um'])
        for n in range(len(coords)):
            for c in range(3):writer.writerow([int(op['node_ids'][n]),*(coords[n]*1000),COMP[c],bool(valid[n]),maps[0][n,c],maps[1][n,c]])
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005)
    ids={nid:i for i,(nid,_) in enumerate(mesh['nodes'])};cells=[]
    for _,conn in mesh['elements']:cells.extend([8,*[ids[nid] for nid in conn]])
    grid=pv.UnstructuredGrid(np.array(cells),np.full(len(mesh['elements']),12,dtype=np.uint8),coords*1000)
    for st in range(2):
        for c in range(3):grid.point_data[f'J_step{st+1}_{COMP[c]}_um']=maps[st][:,c]
    (RES/'vtk').mkdir(exist_ok=True);grid.save(RES/'vtk/sensor_score_idw.vtu')
    pv.OFF_SCREEN=True
    fig,axs=plt.subplots(2,3,figsize=(16,12));fig.subplots_adjust(top=.87,bottom=.10,hspace=.19,wspace=.03,right=.89)
    for st in range(2):
        finite=maps[st][np.isfinite(maps[st])];clim=[float(finite.min()),float(finite.max())]
        for c,ax in enumerate(axs[st]):
            plot=pv.Plotter(off_screen=True,window_size=(700,640));plot.set_background('white')
            plot.add_mesh(grid,scalars=f'J_step{st+1}_{COMP[c]}_um',cmap='RdYlBu',clim=clim,nan_color='#cccccc',show_scalar_bar=False)
            plot.add_point_labels(coords[[iA,iO]]*1000+[0,0,7],['A','O'],font_size=28,shape=None,point_color='black',text_color='black',always_visible=True)
            for k,pick in enumerate(selected[:st+1]):
                ni,ci=divmod(pick,3)
                if ci!=c:continue
                pt=coords[ni]*1000;plot.add_mesh(pv.Sphere(radius=4.0,center=pt),color='#cc00cc',show_edges=False)
                direction=np.eye(3)[ci]
                plot.add_mesh(pv.Arrow(start=pt,direction=direction,scale=27, tip_length=.32, tip_radius=.13, shaft_radius=.055),color='#ff00cc',lighting=False)
                plot.add_point_labels([pt+[0,0,14]],[f'S{k+1} {COMP[ci]}'],font_size=25,text_color='#990099',shape=None,always_visible=True)
            plot.camera_position=[(180,-220,200),(0,0,50),(0,0,1)];plot.camera.zoom(1.30)
            ax.imshow(plot.screenshot(return_img=True));plot.close();ax.axis('off')
            ax.set_title(COMP[c],fontsize=22)
        barax=fig.add_axes([.91,.55 if st==0 else .13,.018,.29])
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(*clim),cmap='RdYlBu'),cax=barax)
        cb.set_label('$J_c$ [µm]：小さいほど良い',fontsize=15);cb.ax.tick_params(labelsize=14)
    fig.suptitle('反り A−O の不確かさを減らす観測位置：候補の採点と実際の選択',fontsize=22,weight='bold')
    fig.text(.5,.945,'温度P2・P0は固定／4000本のROM事前分布／変位応答は8セルIDW＋FrontISTR',ha='center',fontsize=16)
    fig.text(.5,.91,'測定候補',ha='center',fontsize=22,weight='bold')
    fig.text(.02,.87,'① 変位1本目を選ぶ前',fontsize=18,weight='bold')
    fig.text(.02,.48,'② S1を観測に加え、2本目を選ぶ前',fontsize=18,weight='bold')
    fig.text(.5,.045,'紫の球・矢印＝選択した位置・測定方向。灰＝A/O近傍12 mmの除外領域。\n上段と下段は別の色範囲（各段の3方向は共通）。これは熱感度Wそのものではなく、観測後に残る反りの標準偏差。',ha='center',fontsize=14)
    fig.savefig(IMG/'sensor_score_idw_oblique.png',dpi=160,facecolor='white');plt.close(fig)

if __name__=='__main__':main()
