"""実際のROM座標・校正済みDを用いた説明用温度ベクトルの計算例。
入力温度と補間位置は説明用。特定時刻の解析結果ではない。
"""
from pathlib import Path
import json,sys
import numpy as np
import yaml
HERE=Path(__file__).absolute().parent
sys.path.insert(0,str(HERE.parent))
from dacore.node_locations import NODE_XYZ
from dacore.displacement import displacement

def calculate():
    xyz=np.array(list(NODE_XYZ.values()))
    targets=np.array([[.030,0,.05025],[-.030,0,.05025],[.025,0,.095]])
    d=np.linalg.norm(targets[:,None,:]-xyz[None,:,:],axis=2)
    W=np.zeros_like(d)
    for i,row in enumerate(d):
        if row.min()<1e-9:W[i,row.argmin()]=1
        else:
            w=1/row**2;W[i]=w/w.sum()
    T_C=np.array([25,23,21,24,22.])
    op=yaml.safe_load((HERE.parent/'config/displacement_operator.yaml').read_text())
    D_mm=np.array(op['D'])
    delta_T=T_C+273.15-op['T_ref_K']
    u_mm=displacement(T_C+273.15,op)
    return dict(xyz=xyz,targets=targets,distance_m=d,W=W,T_C=T_C,Tmesh_C=W@T_C,
                D_mm=D_mm,delta_T=delta_T,contribution_um=D_mm*1000*delta_T,u_um=u_mm*1000)
if __name__=='__main__':
    v=calculate();HERE.joinpath('rom_numeric_example.json').write_text(json.dumps({k:a.tolist() for k,a in v.items()},indent=2))
    for k,a in v.items():print(k,'=',a)
