"""理論解説の行列恒等式・数値例を検算する（解析ケースは変更しない）。"""
from pathlib import Path
import json
import numpy as np

def main():
    B=np.array([[4.,1.,.5],[1.,3.,.2],[.5,.2,2.]])
    H=np.array([[1.,0.,.3],[0.,1.,0.]])
    R=np.diag([1.,.5]);xb=np.array([20.,18.,2.]);y=np.array([22.,19.])
    Bi=np.linalg.inv(B);Ri=np.linalg.inv(R)
    A=Bi+H.T@Ri@H;S=H@B@H.T+R
    K=np.linalg.solve(S,H@B).T
    xa=xb+K@(y-H@xb)
    xa_normal=np.linalg.solve(A,Bi@xb+H.T@Ri@y)
    np.testing.assert_allclose(xa,xa_normal,atol=1e-12)
    np.testing.assert_allclose(A@K,H.T@Ri,atol=1e-12)
    np.testing.assert_allclose(Bi@(xa-xb)-H.T@Ri@(y-H@xa),0,atol=1e-12)
    J=np.eye(3)-K@H
    Pa=J@B@J.T+K@R@K.T
    np.testing.assert_allclose(Pa,B-B@H.T@np.linalg.solve(S,H@B),atol=1e-12)
    np.testing.assert_allclose(Pa,np.linalg.inv(A),atol=1e-12)
    assert np.linalg.eigvalsh(Pa).min()>0
    mean=20.;P=4.;kf=[]
    for yy in [22.,21.]:
        mf=.9*mean+2;Pf=.9**2*P+.16;gain=Pf/(Pf+1)
        mean=mf+gain*(yy-mf);P=(1-gain)**2*Pf+gain**2
        kf.append(dict(forecast_mean=mf,forecast_variance=Pf,gain=gain,analysis_mean=mean,analysis_variance=P))
    np.testing.assert_allclose([r['analysis_mean'] for r in kf],[21.545454545454547,21.218885212522274])
    particles=np.arange(18.,23.);logL=-.5*(21.5-particles)**2/.25
    w=np.exp(logL-logL.max());w/=w.sum();ess=1/(w@w)
    np.testing.assert_allclose(w@particles,21.486380358560382)
    np.testing.assert_allclose(ess,2.0364698303474653)
    resampled=particles[np.searchsorted(np.cumsum(w),.08+np.arange(5)/5)]
    np.testing.assert_array_equal(resampled,[21,21,21,22,22])
    out=dict(oi_normal_equation_and_gain_equivalent=True,oi_gradient_zero=True,oi_joseph_and_inverse_hessian_equal=True,kf=kf,pf=dict(weights=w.tolist(),mean=float(w@particles),ess=float(ess),resampled=resampled.tolist()))
    Path(__file__).absolute().with_name('theory_validation.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
    print('Verified: OI normal equations/gain/Joseph covariance; KF two steps; PF likelihood/ESS/resampling.')
if __name__=='__main__':main()
