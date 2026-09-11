"""スライドの説明用仮想例。104の保存結果の再現ではない。
温度2列＋Qに縮小。平均・共分散の手計算を追うためinflation=1。
実ソルバや既存ケースは実行・変更しない。
"""
from pathlib import Path
import json
import numpy as np

def calculate():
    Zf = np.array([[290,296,10], [292,294,14], [294,298,12],
                   [296,300,18], [298,292,16]], dtype=float)
    Yf = np.column_stack((Zf[:,:2], np.arange(1,6)*.001, np.arange(-2,3)*.001))
    zbar, ybar = Zf.mean(axis=0), Yf.mean(axis=0)
    dZ, dY = Zf-zbar, Yf-ybar
    C_zy, C_yy = dZ.T@dY/4, dY.T@dY/4
    R = np.diag([.30**2,.30**2,1e-4**2,1e-4**2])
    S = C_yy+R
    gain_T = np.linalg.solve(S,C_zy.T)
    K = gain_T.T
    y = np.array([293,295,.0025,-.0005])
    eps = np.random.default_rng(104).multivariate_normal(np.zeros(4),R,size=5)
    innov = y[None,:]+eps-Yf
    delta = innov@gain_T
    Za = Zf+delta
    return {name: value for name,value in locals().items() if isinstance(value,np.ndarray)}

if __name__=='__main__':
    data=calculate()
    for name,value in data.items():print(name,'=',value)
    out=Path(__file__).absolute().with_suffix('.json')
    out.write_text(json.dumps({k:v.tolist() for k,v in data.items()},ensure_ascii=False,indent=2))
