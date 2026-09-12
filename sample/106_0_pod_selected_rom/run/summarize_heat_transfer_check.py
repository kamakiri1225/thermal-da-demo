"""OpenFOAM出力の壁面hと熱流束からの逆算を同じ定義で比較する。"""
from pathlib import Path
import json
import numpy as np
import pyvista as pv
ROOT=Path(__file__).absolute().parents[1]
r=pv.POpenFOAMReader(str(ROOT.parent/'102_0_openfoam_hollow_cylinder_heat_transfer/hollowCylinder.foam'))
r.enable_all_cell_arrays();out=[]
for t in [300,600]:
    r.set_active_time_value(t)
    patch=r.read()['solid']['boundary']['solid_to_fluid']
    T=np.asarray(patch.cell_data['T']);q=np.asarray(patch.cell_data['wallHeatFlux'])
    h=np.asarray(patch.cell_data['heatTransferCoeff(T)']);area=patch.compute_cell_sizes()['Area']
    calc=-q/(T-293.15)
    row=dict(time_s=t,area_m2=float(area.sum()),outward_heat_W=float((-q*area).sum()),
        wall_mean_T_K=float(np.average(T,weights=area)),
        effective_htc_W_m2_K=float((-q*area).sum()/((T-293.15)*area).sum()),
        area_mean_local_htc_W_m2_K=float(np.average(h,weights=area)),
        max_abs_output_vs_formula_W_m2_K=float(abs(h-calc).max()))
    out.append(row)
(ROOT/'results/heat_transfer_check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
