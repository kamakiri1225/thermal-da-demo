"""104の保存済み同化後場から温度履歴・FrontISTR変位差・温度GIFを生成。"""
from pathlib import Path
import sys, json
import numpy as np
ROOT=Path(__file__).absolute().parents[1]; SAMPLE=ROOT.parent
S104=SAMPLE/'104_0_openfoam_frontistr_da_enkf'
sys.path.insert(0,str(S104));sys.path.insert(0,str(SAMPLE/'102_1_frontistr_hollow_cylinder_thermal_expansion/python'))
from dacore import plots
from daof import of_case
from fem.fem_obs import MATERIAL
import cylinder_mesh, fistr_case
from openfoam_temperature import align_cell_centers_to_node_mesh, interpolate_to_nodes
import matplotlib.pyplot as plt
from PIL import Image
import pyvista as pv

def main():
    res=ROOT/'results';img=ROOT/'docs/img';work=ROOT/'openfoam/fullsolver_evidence';work.mkdir(parents=True,exist_ok=True)
    src=S104/'openfoam/run_fem_enkf';fld=src/'fields'
    times=np.arange(0,601,60)
    centres=np.load(fld/'cell_centres.npy')
    analysis=np.array([np.load(fld/f'ensmean_t{t}.npy') for t in times])
    truth=np.array([np.full(len(centres),293.15) if t==0 else of_case.read_solid_T(str(src/'truth'),float(t)) for t in times])
    probes=np.array([[.032,0,.05025],[-.032,0,.05025]])
    cells=((centres[:,None,:]-probes[None,:,:])**2).sum(2).argmin(0)
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005)
    ids=[n for n,_ in mesh['nodes']];coords=np.array([v for _,v in mesh['nodes']])
    aligned,_=align_cell_centers_to_node_mesh(centres,coords)
    groups=[mesh['top_heater_nodes'],mesh['top_opposite_nodes']]
    from scipy.spatial import cKDTree
    distances,near=cKDTree(aligned).query(coords,k=8)
    weights=1/np.maximum(distances,1e-30);weights/=weights.sum(1,keepdims=True)
    exact=distances[:,0]<1e-9;weights[exact]=0;weights[exact,0]=1
    ta=np.sum(analysis[:,near]*weights[None,:,:],axis=2)
    tt=np.sum(truth[:,near]*weights[None,:,:],axis=2)
    ua=[];ut=[]
    for k,t in enumerate(times):
        case=work/f'analysis_t{t}'
        if not (case/'hollow_cylinder_thermal_expansion.res.0.1').exists():
            case.mkdir(exist_ok=True)
            fistr_case.write_mesh(case,4,48,20,.020,.0375,.1005,young_modulus=MATERIAL['young_modulus_Pa'],poisson_ratio=MATERIAL['poisson_ratio'],density=MATERIAL['density_kg_m3'],thermal_expansion_coeff=MATERIAL['thermal_expansion_coeff_per_K'])
            fistr_case.write_hecmw_ctrl(case)
            fistr_case.write_cnt(case,ids,ta[k],reference_temperature=MATERIAL['reference_temperature_K'],young_modulus=MATERIAL['young_modulus_Pa'],poisson_ratio=MATERIAL['poisson_ratio'],thermal_expansion_coeff=MATERIAL['thermal_expansion_coeff_per_K'])
            fistr_case.run_fistr(case)
        disp=fistr_case.read_displacement(case)
        ua.append([np.mean([disp[n][2] for n in g])*1e6 for g in groups])
        if t==0:ut.append([0.,0.])
        else:
            disp=fistr_case.read_displacement(src/'truth'/f'fem_t{t}')
            ut.append([np.mean([disp[n][2] for n in g])*1e6 for g in groups])
        print('FrontISTR analysis',t,flush=True)
    ua=np.array(ua);ut=np.array(ut)
    np.savez(res/'fullsolver_evidence.npz',time=times,temperature_analysis_K=analysis[:,cells],temperature_truth_K=truth[:,cells],probe_cells=cells,probe_xyz=centres[cells],u_analysis_um=ua,u_truth_um=ut)
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for j,ax in enumerate(axes):
        ax.plot(times,truth[:,cells[j]]-273.15,'k-',label='OpenFOAM真値')
        ax.plot(times,analysis[:,cells[j]]-273.15,'o--',label='同化後平均（5メンバー）',ms=4)
        ax.set_title(['hot：ヒータ側','cold：反ヒータ側'][j]);ax.set_ylabel('温度 [℃]');ax.set_xlabel('時刻 [s]');ax.grid(alpha=.3);ax.legend()
    fig.tight_layout();fig.savefig(img/'fullsolver_temperature_history.png',dpi=140);plt.close(fig)
    da=ua[:,0]-ua[:,1];dt=ut[:,0]-ut[:,1]
    fig,axes=plt.subplots(2,1,figsize=(10,8),sharex=True)
    axes[0].plot(times,dt,'k-',lw=2,label='実ソルバ真値')
    axes[0].plot(times,da,'o--',ms=4,label='同化後温度→FrontISTR再計算')
    axes[0].set_ylabel('ヒータ側 − 反ヒータ側のUz [µm]');axes[0].legend()
    axes[0].set_title('104：上面2測定領域の変位差（各領域の節点平均の差）')
    axes[1].plot(times,da-dt,'o-');axes[1].axhline(0,color='k',lw=1)
    axes[1].set_ylabel('変位差の推定誤差 [µm]');axes[1].set_xlabel('時刻 [s]')
    for ax in axes:ax.axvspan(0,300,color='orange',alpha=.08);ax.grid(alpha=.3)
    fig.tight_layout();fig.savefig(img/'fullsolver_displacement_difference.png',dpi=140);plt.close(fig)
    pv.OFF_SCREEN=True
    index={n:i for i,n in enumerate(ids)};conn=[]
    for _,ns in mesh['elements']:conn.extend([8]+[index[n] for n in ns])
    grid=pv.UnstructuredGrid(np.array(conn),np.full(len(mesh['elements']),12,np.uint8),coords)
    clim=(float(min(ta.min(),tt.min())-273.15),float(max(ta.max(),tt.max())-273.15))
    frames=[]
    for k,t in enumerate(times):
        pl=pv.Plotter(shape=(1,2),off_screen=True,window_size=(1000,600))
        for j,temps in enumerate([ta[k],tt[k]]):
            pl.subplot(0,j);g=grid.copy();g['T']=temps-273.15
            pl.add_mesh(g,scalars='T',cmap='turbo',clim=clim,n_colors=16,scalar_bar_args={'title':'Temperature [C]'})
            pl.add_text(('Analysis mean' if j==0 else 'Truth')+f' / t={t}s',font_size=13)
            pl.camera_position=[(.24,-.22,.20),(0,0,.05),(0,0,1)];pl.set_background('white')
        frames.append(Image.fromarray(pl.screenshot()));pl.close()
    frames[0].save(img/'fullsolver_temperature_comparison.gif',save_all=True,append_images=frames[1:],duration=600,loop=0)
    (res/'fullsolver_evidence_summary.json').write_text(json.dumps({'difference_error_final_um':float(da[-1]-dt[-1]),'difference_max_abs_error_after_update_um':float(abs(da[1:]-dt[1:]).max()),'probe_xyz_m':centres[cells].tolist(),'temperature_clim_C':clim},indent=2)+'\n')
    print('done',da[-1]-dt[-1],flush=True)

if __name__=='__main__':main()
