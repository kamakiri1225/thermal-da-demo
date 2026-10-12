"""blog011: conditional linear covariance + 1D quadrature vs MC4000 and tangent.
Existing ROM/EnKF unchanged. Outputs isolated in results/blog011, docs/img/blog011.
Run: OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python run/analytic_covariance_sensor_study.py
"""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.special import ndtr
from numpy.polynomial.legendre import leggauss
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run')]
from dacore import rom_general as rg
from select_nonrep_temp_disp import priors,score,add_row,R_T,R_U,A_XYZ,O_XYZ,EXCL
import selected_pair_unobserved_check as evaluation
RES=ROOT/'results/blog011';IMG=ROOT/'docs/img/blog011';TIMES=np.arange(30.,301.,30.)

def h_quadrature(order):
    # clip(N(.02,.01^2), .001,.1): continuous density plus atoms at endpoints.
    x,w=leggauss(order);lo,hi=.001,.1
    h=(hi+lo)/2+(hi-lo)/2*x
    weights=w*(hi-lo)/2*np.exp(-.5*((h-.02)/.01)**2)/(.01*np.sqrt(2*np.pi))
    h=np.r_[lo,h,hi];weights=np.r_[ndtr((lo-.02)/.01),weights,ndtr((.02-hi)/.01)]
    return h,weights

def model():
    d=np.load(RES/'rom_calibrated_pod.npz');return d['C'],rg.tri_to_matrix(d['K_upper'],5),int(d['heat_node'])

def conditional_covariance(order):
    C,K,heat=model();hs,weights=h_quadrature(order);n=len(hs)
    # Seven deterministic responses per h: baseline, five initial-temperature units, q unit.
    T=np.full((n,7,5),rg.T_AIR_K);T[:,1:6,:]+=np.eye(5)[None,:,:]
    qs=np.zeros((n,7));qs[:,6]=1
    mu=np.r_[np.full(5,4.5),1.05];V=np.diag(np.r_[np.full(5,15.**2/12),1.5**2/12])
    Bs=[];means=[];tp=0.
    for tb in TIMES:
        T=rg.integrate_ensemble(T.reshape(-1,5),C,K,np.repeat(hs,7),qs.ravel(),heat,tp,tb,2.).reshape(n,7,5);tp=tb
        F=(T[:,1:,:]-T[:,0:1,:]).transpose(0,2,1) # (h,5,6)
        mh=T[:,0,:]+F@mu
        m=weights@mh
        # Central form avoids cancellation of ~293 K offsets.
        cov=np.einsum('nik,kl,njl->nij',F,V,F)
        dm=mh-m
        B=np.einsum('n,nij->ij',weights,cov)+np.einsum('n,ni,nj->ij',weights,dm,dm)
        Bs.append(B);means.append(m)
    return np.array(Bs),np.array(means),float(weights.sum())

def tangent_covariance():
    C,K,heat=model();hs,p=h_quadrature(64);hm=p@hs;hv=p@((hs-hm)**2)
    T=np.full(5,rg.T_AIR_K+4.5);S=np.zeros((5,7));S[:,:5]=np.eye(5)
    V=np.diag(np.r_[np.full(5,15.**2/12),1.5**2/12,hv]);L=K-np.diag(K.sum(1));A=(L-hm*np.eye(5))/C[:,None]
    def rhs(T,S,t):
        b=np.zeros(5);b[heat]=15.*float(t<300)/C[heat]
        f=A@(T-rg.T_AIR_K)+b*1.05
        dS=A@S;dS[:,5]+=b;dS[:,6]-=(T-rg.T_AIR_K)/C
        return f,dS
    Bs=[]
    for step in range(150):
        t=2.*step;dt=2.
        a,b=rhs(T,S,t);c,d=rhs(T+a,S+b,t+1);e,f=rhs(T+c,S+d,t+1);g,h=rhs(T+2*e,S+2*f,t+2)
        T+=dt/6*(a+2*c+2*e+g);S+=dt/6*(b+2*d+2*f+h)
        if (step+1)%15==0:Bs.append(S@V@S.T)
    return np.array(Bs)

def operators():
    kv=np.load(RES/'qdeim_points.npz');op=np.load(RES/'dispop_allnodes.npz')
    cc=kv['cell_centres'];coords=op['coords'];inv=np.linalg.pinv(kv['pod_modes'][kv['cell_idx']]);G=kv['pod_modes']@inv
    H=np.einsum('nck,kj->ncj',op['D'],inv).reshape(-1,5)
    nn=lambda x:int(np.linalg.norm(coords-x,axis=1).argmin())
    w=H[3*nn(A_XYZ)+2]-H[3*nn(O_XYZ)+2]
    tm=np.ones(len(cc),bool);tm[kv['cell_idx']]=False
    um=np.repeat((np.linalg.norm(coords-A_XYZ,axis=1)>EXCL)&(np.linalg.norm(coords-O_XYZ,axis=1)>EXCL),3)
    return cc,coords,G,H,w,tm,um

def select(Bs):
    cc,coords,G,H,w,tm,um=operators();Ps=list(Bs.copy());out={'temperature':[],'displacement':[]};fields={}
    for kind,rows,r,mask in [('temperature',G,R_T,tm),('displacement',H,R_U,um)]:
        for k in range(2):
            sc=score(Ps,w,rows,r);idx=int(np.where(mask,sc,np.inf).argmin());mask[idx]=False
            rec={'index':idx,'xyz_m':(cc[idx] if kind=='temperature' else coords[idx//3]).tolist(),'J_um':float(sc[idx])}
            if kind=='displacement':rec['component']=idx%3
            out[kind].append(rec);fields[f'{kind}_{k+1}']=sc
            Ps=[add_row(P,rows[idx],r) for P in Ps]
    return out,fields

def main():
    RES.mkdir(exist_ok=True);IMG.mkdir(exist_ok=True);(RES/'vtk').mkdir(exist_ok=True)
    for name in ['rom_calibrated_pod.npz','qdeim_points.npz','dispop_allnodes.npz']+[f'{p}_{c}.npz' for c in evaluation.CASES for p in ['limit_truth','truth_disp_all']]:
        dst=RES/name
        if not dst.exists():dst.symlink_to(Path('..')/name)
    timing={};Bs={};convergence={}
    for order in (8,16,32,64):
        t=time.perf_counter();B,mean,mass=conditional_covariance(order);timing[f'quadrature_{order}']=time.perf_counter()-t
        convergence[str(order)]={'weight_sum':mass};Bs[f'quad{order}']=B
    assert abs(convergence['32']['weight_sum']-1)<1e-10
    for order in (8,16,32):
        convergence[str(order)]['relative_B_vs_64']=float(np.linalg.norm(Bs[f'quad{order}']-Bs['quad64'])/np.linalg.norm(Bs['quad64']))
    assert convergence['32']['relative_B_vs_64']<1e-8
    t=time.perf_counter();Bs['MC4000']=np.array(priors());timing['MC4000']=time.perf_counter()-t
    t=time.perf_counter();Bs['tangent']=tangent_covariance();timing['tangent']=time.perf_counter()-t
    # Validate the affine response against independent trajectories at an off-grid h.
    C,K,heat=model();hh=.027;T0=rg.T_AIR_K+np.array([1.,-2.,4.,8.,3.]);q=1.37
    base=np.full((7,5),rg.T_AIR_K);base[1:6]+=np.eye(5);qq=np.r_[np.zeros(6),1.]
    batch=rg.integrate_ensemble(base,C,K,np.full(7,hh),qq,heat,0.,150.,2.)
    predicted=batch[0]+(batch[1:]-batch[0]).T@np.r_[T0-rg.T_AIR_K,q]
    direct=rg.integrate_ensemble(T0[None,:],C,K,np.array([hh]),np.array([q]),heat,0.,150.,2.)[0]
    affine_error=float(np.max(np.abs(predicted-direct)));assert affine_error<1e-9
    sels={};scores={}
    for key in ('MC4000','quad32','tangent'):
        t=time.perf_counter();sels[key],scores[key]=select(Bs[key]);timing[f'selection_{key}']=time.perf_counter()-t
        print(key,json.dumps(sels[key]),flush=True)
        np.savez_compressed(RES/f'scores_{key}.npz',**scores[key])
        (RES/f'selection_{key}.json').write_text(json.dumps(sels[key],ensure_ascii=False,indent=2))
        for B in Bs[key]:
            assert np.allclose(B,B.T,atol=1e-12)
            assert np.linalg.eigvalsh(B).min()>-1e-9
    np.savez_compressed(RES/'covariances.npz',times=TIMES,**Bs)
    info={'timing_s':timing,'quadrature_convergence':convergence,'affine_validation_max_K':affine_error,
          'relative_B_MC_vs_quad32':float(np.linalg.norm(Bs['MC4000']-Bs['quad32'])/np.linalg.norm(Bs['quad32'])),
          'relative_B_tangent_vs_quad32':float(np.linalg.norm(Bs['tangent']-Bs['quad32'])/np.linalg.norm(Bs['quad32'])),
          'selections':sels,'prior_seed':20260930,'EnKF_seeds':evaluation.SEEDS,'EnKF_members':60,
          'note':'Same priors/candidates/target as blog009. Quadrature includes clipped-normal endpoint masses. RK4 dt=2 s, same 300 s endpoint rule.'}
    (RES/'study.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
    evaluation.RES=str(RES);evaluation.IMG=str(IMG)
    def configs(_):
        cfg=[('同化なし',None,[],'#999999','--'),('手動温度 P2・P0',evaluation.P2P0,[],'#E69F00','--')]
        for key,col in [('MC4000','#0072B2'),('quad32','#CC3311'),('tangent','#009E73')]:
            s=sels[key];ts=[('xyz',tuple(p['xyz_m'])) for p in s['temperature']];us=[(tuple(p['xyz_m']),p['component']) for p in s['displacement']]
            cfg.append((key+'：温度2点',ts,[],col,'--'))
            cfg.append((key+'：温度＋変位4本',ts,us,col,'-'))
        return '4,000ケースと熱方程式に基づく観測配置の比較',cfg
    evaluation.configs=configs
    for case in evaluation.CASES:
        t=time.perf_counter();evaluation.main('analytic',case);info['timing_s'][f'EnKF_{case}']=time.perf_counter()-t
    (RES/'study.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
    print('blog011 calculation and checks complete',json.dumps(info['timing_s']),flush=True)
if __name__=='__main__':main()
