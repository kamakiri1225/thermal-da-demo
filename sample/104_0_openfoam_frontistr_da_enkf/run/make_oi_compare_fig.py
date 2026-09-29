"""Consolidated OI sigma_Q comparison. Recompute analysis Uz from saved T fields.
Run from 104: OMP_NUM_THREADS=4 python3 run/make_oi_compare_fig.py
"""
from pathlib import Path
import sys, hashlib, json, shutil
import numpy as np
ROOT=Path(str(Path(__file__).absolute().parent.parent).replace('/002_CAE/','/002_cae/'))
sys.path.insert(0,str(ROOT))
from daof import of_case
from fem.fem_obs import displacement_obs, rte
from dacore import plots as _plots
import matplotlib.pyplot as plt
RUNS=[('sq0p5',0.5,'tab:purple','o'),('b',2,'tab:green','s'),('',6,'tab:blue','^'),('sq50',50,'tab:brown','x')]
PROBES=[('hot（観測）',(0.032,0,.05025)),('cold（観測）',(-.032,0,.05025)),('上部ヒータ側（未観測）',(.030,0,.090)),('底部反対側（未観測）',(-.030,0,.015))]

def main():
    res=ROOT/'results'; out=ROOT/'docs/img/oi_fullsolver_compare.png'
    runs={s:dict(np.load(res/f"oi_fullsolver{'_'+s if s else ''}_history.npz")) for s,_,_,_ in RUNS}
    free=dict(np.load(res/'oi_fullsolver_free_history.npz'))
    ts=runs['']['times']; cells=of_case.nearest_cells([p for _,p in PROBES],of_case.solid_cell_centres())
    truth=ROOT/'openfoam/run_fem_enkf/truth'
    def temperatures(case,initial):
        return np.array([np.full(4,initial) if t==0 else of_case.read_solid_T(str(case),t)[cells] for t in ts])-273.15
    tt=temperatures(truth,293.15); tf=temperatures(ROOT/'openfoam/run_fem_oi_free',298.15)
    temp={}; disp={}
    for s,sig,_,_ in RUNS:
        case=ROOT/f"openfoam/run_fem_oi{'_'+s if s else ''}"
        temp[s]=temperatures(case,298.15)
        # Re-evaluate all cases consistently; legacy u_analysis could contain forecast Uz.
        digest=hashlib.sha256()
        for t in ts[1:]: digest.update((case/f'{t:g}/solid/T').read_bytes())
        digest.update((ROOT/'fem/fem_obs.py').read_bytes())
        cache=res/f"oi_sigma_unified_uz_{s or '6'}_{digest.hexdigest()[:12]}.npz"
        if cache.exists(): u=np.load(cache)['u_mm']
        else:
            u=np.full((len(ts),2),np.nan) # no physical initial FEM result in legacy history
            for i,t in enumerate(ts[1:],1):
                u[i]=displacement_obs(str(case),f'{t:g}',str(ROOT/f'openfoam/oi_sigma_unified_fem/{s or "6"}/{t:g}'))
                print(f'sigma={sig} analysis FEM t={t:g}',flush=True)
            np.savez_compressed(cache,times=ts,u_mm=u)
        disp[s]=u*1000
    # Holdout points fixed before viewing results: middle-height +/-X, not upper observed groups.
    import cylinder_mesh
    mesh=cylinder_mesh.build_cylinder_mesh(4,48,20,.020,.0375,.1005)
    xyz=np.array([x for _,x in mesh['nodes']]); ids=np.array([i for i,_ in mesh['nodes']])
    vn=[int(np.linalg.norm(xyz-[x,0,.075],axis=1).argmin()) for x in [.028,-.028]]
    v_ids=ids[vn]
    assert not set(v_ids)&(set(mesh['top_heater_nodes'])|set(mesh['top_opposite_nodes']))
    def read_validation(work):
        dis=rte.fistr_case.read_displacement(work)
        return np.array([dis[int(i)][2] for i in v_ids])*1e6
    vrun={s:np.array([read_validation(ROOT/f'openfoam/oi_sigma_unified_fem/{s or "6"}/{t:g}') for t in ts[1:]]) for s,_,_,_ in RUNS}
    vfree=np.array([read_validation(ROOT/f'openfoam/run_fem_oi_free/fem_t{t:g}') for t in ts[1:]])
    vtruth=[]
    for t in ts[1:]:
        work=ROOT/f'openfoam/oi_sigma_unified_fem/truth/{t:g}'
        if not (work/'hollow_cylinder_thermal_expansion.res.0.1').exists():
            displacement_obs(str(truth),f'{t:g}',str(work))
            print(f'holdout truth FEM t={t:g}',flush=True)
        vtruth.append(read_validation(work))
    vtruth=np.array(vtruth)
    np.savez_compressed(res/'oi_sigma_unified_validation.npz',times=ts[1:],node_ids=v_ids,xyz_m=xyz[vn],truth_um=vtruth,free_um=vfree,**{f'analysis_{s or "6"}_um':v for s,v in vrun.items()})
    fig,axes=plt.subplots(3,4,figsize=(21,11.5))
    labels={s:rf'$\sigma_Q$={sig:g} W' for s,sig,_,_ in RUNS}
    for k in range(4):
        ax=axes.flat[k];ax.plot(ts,tt[:,k],color='black',lw=2.5,label='真値')
        ax.plot(ts,tf[:,k],':',color='tab:orange',lw=2,label='同化なし')
        for s,sig,c,m in RUNS: ax.plot(ts,temp[s][:,k],ls='--',marker=m,ms=4,color=c,label=labels[s])
        ax.set_title(PROBES[k][0]+'：同化後温度');ax.set_ylabel('温度 [℃]')
    for k,title in enumerate(['A：ヒータ側上面 平均Uz','O：反対側上面 平均Uz','A−O：2点間変位差']):
        ax=axes.flat[4+k]
        def component(u): return u[:,k] if k<2 else u[:,0]-u[:,1]
        # Exclude artificial initial zeros stored by the old driver.
        ax.plot(ts[1:],component(runs['']['u_truth']*1000)[1:],color='black',lw=2.5)
        ax.plot(free['times'][1:],component(free['u_analysis']*1000)[1:],':',color='tab:orange',lw=2)
        for s,sig,c,m in RUNS: ax.plot(ts,component(disp[s]),ls='--',marker=m,ms=4,color=c)
        ax.set_title(title+'\n同化後温度から再計算');ax.set_ylabel('変位 [µm]')
    for k,title in enumerate(['C：未観測・中高さヒータ側Uz','D：未観測・中高さ反対側Uz','C−D：未観測の2点間変位差']):
        ax=axes.flat[7+k]
        def component_v(u): return u[:,k] if k<2 else u[:,0]-u[:,1]
        ax.plot(ts[1:],component_v(vtruth),color='black',lw=2.5)
        ax.plot(ts[1:],component_v(vfree),':',color='tab:orange',lw=2)
        for s,sig,c,m in RUNS:ax.plot(ts[1:],component_v(vrun[s]),ls='--',marker=m,ms=4,color=c)
        ax.set_title(title);ax.set_ylabel('変位 [µm]')
    ax=axes.flat[10];ax.axhline(15,color='black',lw=2)
    ax.plot(free['times'],free['q'],':',color='tab:orange',lw=2)
    for s,sig,c,m in RUNS:ax.plot(ts,runs[s]['q'],'--',marker=m,ms=4,color=c)
    ax.set_title('発熱パラメータQ\n300秒以降の実入力は0 W');ax.set_ylabel('Q [W]')
    ax=axes.flat[11]
    ax.plot(free['times'],free['rmse_field'],':',color='tab:orange',lw=2)
    for s,sig,c,m in RUNS:ax.plot(ts,runs[s]['rmse_field'],'--',marker=m,ms=4,color=c)
    ax.set_title('全20,696セルの同化後温度RMSE');ax.set_ylabel('RMSE [K]');ax.set_yscale('log')
    for ax in axes.flat:
        ax.title.set_fontsize(12)
        ax.axvspan(0,300,color='orange',alpha=.05);ax.axvline(300,color='gray',ls=':',lw=1);ax.grid(alpha=.25);ax.set_xlabel('時刻 [s]')
    handles,labs=axes.flat[0].get_legend_handles_labels()
    fig.legend(handles,labs,loc='upper center',bbox_to_anchor=(.5,.95),ncol=6,fontsize=12,borderaxespad=0)
    fig.suptitle('実ソルバOI：σQの4条件を温度・変位・変位差・Qで同時比較\n背景1本／60秒ごと10回更新／加熱0–300秒・冷却300–600秒／変位は保存温度場からFrontISTRで再評価',fontsize=17)
    fig.subplots_adjust(left=.05,right=.99,bottom=.06,top=.855,wspace=.26,hspace=.48);fig.savefig(out,dpi=140,bbox_inches='tight',pad_inches=.12);plt.close(fig)
    dst=ROOT.parent/'106_0_pod_selected_rom/docs/img/oi_fullsolver_compare.png';shutil.copyfile(out,dst)
    metrics={}
    for s,sig,_,_ in RUNS:
        e=disp[s]-runs[s]['u_truth']*1000;heat=(ts>0)&(ts<=300)
        metrics[str(sig)]={'Q_final_W':float(runs[s]['q'][-1]),'T_field_final_RMSE_K':float(runs[s]['rmse_field'][-1]),'heating_Uz_RMSE_um':np.sqrt(np.mean(e[heat]**2,axis=0)).tolist(),'heating_difference_RMSE_um':float(np.sqrt(np.mean((e[heat,0]-e[heat,1])**2)))}
    for s,sig,_,_ in RUNS:
        err=vrun[s]-vtruth;heat=ts[1:]<=300
        metrics[str(sig)]['unobserved_heating_Uz_RMSE_um']=np.sqrt(np.mean(err[heat]**2,axis=0)).tolist()
        metrics[str(sig)]['unobserved_heating_difference_RMSE_um']=float(np.sqrt(np.mean((err[heat,0]-err[heat,1])**2)))
    (res/'oi_sigma_unified_metrics.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics,indent=2))
if __name__=='__main__':main()
