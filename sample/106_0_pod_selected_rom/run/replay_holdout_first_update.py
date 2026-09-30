"""Replay one existing holdout experiment update using cached operators.

No OpenFOAM or FrontISTR runs. The replay must agree with the saved history.
Run: python3 run/replay_holdout_first_update.py
"""
from pathlib import Path
import sys
import json
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from dacore import rom_general as rg

def main():
    results = ROOT/'results'
    cal = np.load(results/'rom_calibrated_pod.npz')
    pod = np.load(results/'qdeim_points.npz')
    history = np.load(results/'holdout_support_history.npz')
    op = np.load(results/'dispop_holdout_support_a4dadbb851dca4d7.npz')
    cells = pod['cell_idx']
    pinv = np.linalg.pinv(pod['pod_modes'][cells].astype(float))
    mean = pod['mean'][cells]
    def displacement(T):
        return op['uz_mean'] + ((T-mean)@pinv.T)@op['D'].T
    seed = 20260913
    initial = np.random.default_rng(seed)
    Z = np.zeros((60,7))
    Z[:,:5] = initial.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(60,5))
    Z[:,5] = initial.uniform(.3,1.8,60)
    Z[:,6] = np.clip(initial.normal(.02,.01,60),.001,.1)
    Z[:,:5] = rg.integrate_ensemble(Z[:,:5],cal['C'],rg.tri_to_matrix(cal['K_upper'],5),
                                  Z[:,6],Z[:,5],int(cal['heat_node']),0,30,2)
    Y = np.column_stack((Z[:,[2,0]],displacement(Z[:,:5])[:,:2]))
    obs = np.r_[history['truth_T'][1,[2,0]]+
                np.random.default_rng(seed+100).normal(0,.3,(20,2))[0],
                history['truth_u_um'][1,:2]+
                np.random.default_rng(seed+200).normal(0,.3,(20,2))[0]]
    zm,ym = Z.mean(0),Y.mean(0)
    Zi,Yi = zm+1.02*(Z-zm),ym+1.02*(Y-ym)
    dZ,dY = Zi-Zi.mean(0),Yi-Yi.mean(0)
    Czy,Cyy = dZ.T@dY/59,dY.T@dY/59
    R = np.eye(4)*.09
    gain = np.linalg.solve(Cyy+R,Czy.T).T
    eps = np.random.default_rng(seed+300).multivariate_normal(np.zeros(4),R,size=60)
    Za = Zi+(obs+eps-Yi)@gain.T
    raw_Q = 15*Za[:,5].mean()
    clipped_q = int(np.count_nonzero((Za[:,5]<0)|(Za[:,5]>3)))
    Za[:,5] = np.clip(Za[:,5],0,3)
    Za[:,6] = np.clip(Za[:,6],.0001,.2)
    np.testing.assert_allclose(Za.mean(0),history['states'][2,0,1],rtol=0,atol=1e-10)
    report = dict(seed=seed,time_s=30,matched_saved_history=True,
                  observation=obs.tolist(),forecast_mean=ym.tolist(),
                  mean_observation_perturbation=eps.mean(0).tolist(),gain_q=gain[5].tolist(),
                  Q_contributions_W=(15*gain[5]*(obs+eps.mean(0)-ym)).tolist(),
                  Q_before_W=float(15*zm[5]),Q_after_before_clipping_W=float(raw_Q),
                  Q_after_W=float(15*Za[:,5].mean()),q_clipped_members=clipped_q)
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':
    main()
