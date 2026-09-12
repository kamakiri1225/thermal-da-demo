"""hだけを40%増やす決定論的比較。EnKFを使わず、他の条件は真値に固定する。"""
from pathlib import Path
import sys
import json
import numpy as np
ROOT = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(ROOT))
from dacore import rom_general as rg


def main():
    res = ROOT / 'results'
    d = np.load(res / 'rom_calibrated_pod.npz')
    h = float(d['h'])
    C = d['C']
    km = rg.tri_to_matrix(d['K_upper'], 5)
    args = (np.full(5, rg.T_AIR_K), C, km)
    t, baseline = rg.integrate_single(*args, h, 1., int(d['heat_node']), 0, 600, 2)
    _, perturbed = rg.integrate_single(*args, 1.4*h, 1., int(d['heat_node']), 0, 600, 2)
    pod = np.load(res / 'qdeim_points.npz')
    op = np.load(res / 'disp_operator.npz')
    inverse = np.linalg.pinv(pod['pod_modes'][pod['cell_idx'], :].astype(float))
    delta_T = perturbed-baseline
    delta_u = delta_T @ inverse.T @ op['Dmode'].T
    summary = {
        'h_true_W_K': h, 'h_perturbed_W_K': 1.4*h,
        'max_abs_temperature_difference_K': float(abs(delta_T).max()),
        'max_abs_displacement_difference_um': float(abs(delta_u).max()),
        'approx_uniform_cooling_time_constant_s': float(C.sum()/(5*h)),
        'temperature_observation_sigma_K': .3,
        'displacement_observation_sigma_um': .3,
        'evaluation_interval_s': 2, 'end_time_s': 600,
        'note': 'No assimilation; identical initial temperature and Q. Displacements use cached FrontISTR POD responses.'
    }
    np.savez(res/'h_sensitivity_check.npz', time=t, baseline_T=baseline,
             perturbed_T=perturbed, delta_T_K=delta_T, delta_u_um=delta_u)
    (res/'h_sensitivity_check.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
