"""blog010: A/Oだけをz軸まわり+90度回転し、blog009と同一の選定・EnKFを再実行。
再現: OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 python run/rotated_target_sensor_study.py
既存の選定・EnKF関数を利用。旧結果は変更せず results/blog010/ と docs/img/blog010/ へ出力。
"""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run')]
import sequential_sensor_study as selection
import selected_pair_unobserved_check as evaluation
from select_nonrep_temp_disp import priors

RES=ROOT/'results/blog010';IMG=ROOT/'docs/img/blog010'
A=(0.,.028,.1005);O=(0.,-.028,.1005)

def main():
 RES.mkdir(exist_ok=True);(RES/'vtk').mkdir(exist_ok=True);IMG.mkdir(exist_ok=True)
 # Input arrays are shared unchanged with blog009; only output folder and A/O differ.
 inputs=['qdeim_points.npz','dispop_allnodes.npz','rom_calibrated_pod.npz']+[f'{prefix}_{case}.npz' for case in evaluation.CASES for prefix in ['limit_truth','truth_disp_all']]
 for name in inputs:
  dest=RES/name
  if not dest.exists():dest.symlink_to(Path('..')/name)
 co=np.load(RES/'dispop_allnodes.npz')['coords'];nn=lambda p:int(np.linalg.norm(co-p,axis=1).argmin())
 info={'rotation_deg':90,'rotation_convention':'counterclockwise viewed from +z (A:+x -> +y, O:-x -> -y)',
       'A_requested_m':A,'O_requested_m':O,'A_node_index':nn(A),'O_node_index':nn(O),
       'A_actual_m':co[nn(A)].tolist(),'O_actual_m':co[nn(O)].tolist(),
       'heater':'unchanged, +x side','prior_cases':4000,'prior_seed':20260930,'EnKF_members':60,'EnKF_seeds':evaluation.SEEDS,
       'observation_interval_s':30,'inflation':1.02,'sigma_T_K':.3,'sigma_u_um':.3,'evaluation_interval_s':[30,300],
       'displacement_exclusion_radius_m':.012,'cases':list(evaluation.CASES),'shared_inputs':inputs}
 (RES/'conditions.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
 start=time.perf_counter();Bs=priors();np.savez_compressed(RES/'prior_covariances.npz',B=np.array(Bs),times=np.arange(30,301,30))
 selection.RES=RES;selection.IMG=IMG;selection.A_XYZ=A;selection.O_XYZ=O
 selection.priors=lambda: [B.copy() for B in Bs]
 selection.main();info['selection_elapsed_s']=time.perf_counter()-start
 evaluation.RES=str(RES);evaluation.IMG=str(IMG);evaluation.A_XYZ=A;evaluation.O_XYZ=O
 info['EnKF_elapsed_s']={}
 for case in evaluation.CASES:
  start=time.perf_counter();evaluation.main('seq',case);info['EnKF_elapsed_s'][case]=time.perf_counter()-start
 (RES/'conditions.json').write_text(json.dumps(info,ensure_ascii=False,indent=2))
 print('blog010 calculation complete',json.dumps(info['EnKF_elapsed_s']),flush=True)
if __name__=='__main__':main()
