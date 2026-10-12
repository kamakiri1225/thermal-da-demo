"""150秒の事前共分散と候補追加後の共分散を再計算（blog008/009）。"""
from pathlib import Path
import sys,json,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT),str(ROOT/'run')]
from dacore import rom_general as rg,plots
from select_nonrep_temp_disp import add_row
import matplotlib.pyplot as plt
RES=ROOT/'results';IMG=ROOT/'docs/img'

def main():
 start=time.perf_counter()
 d=np.load(RES/'rom_calibrated_pod.npz');kv=np.load(RES/'qdeim_points.npz');op=np.load(RES/'dispop_allnodes.npz')
 rng=np.random.default_rng(20260930);n=4000
 T=rng.uniform(rg.T_AIR_K-3,rg.T_AIR_K+12,(n,5));q=rng.uniform(.3,1.8,n);h=np.clip(rng.normal(.02,.01,n),.001,.1)
 T=rg.integrate_ensemble(T,d['C'],rg.tri_to_matrix(d['K_upper'],5),h,q,int(d['heat_node']),0,150,2)
 B=np.cov(T.T);inv=np.linalg.pinv(kv['pod_modes'][kv['cell_idx']]);G=kv['pod_modes']@inv
 Wall=np.einsum('nck,kj->ncj',op['D'],inv);co=op['coords'];nn=lambda p:int(np.linalg.norm(co-np.array(p),axis=1).argmin())
 ia=nn((.028,0,.1005));io=nn((-.028,0,.1005));w=Wall[ia,2]-Wall[io,2]
 js=json.loads((RES/'sequential_selection.json').read_text())
 fixed=[('温度P2',np.eye(5)[2],.09),('温度P0',np.eye(5)[0],.09),('旧S1：Ux',Wall[nn((.0097,-.0362,.1005)),0],.09),('旧S2：Uz',Wall[nn((-.0346,.0144,.1005)),2],.09)]
 seq=[(f'温度T{i+1}',G[p['index']],.09) for i,p in enumerate(js['temperature'])]+[(f'変位S{i+1}：Ux',Wall[p['index']//3,p['component']],.09) for i,p in enumerate(js['displacement'])]
 out={'t_s':150,'prior_cases':n,'seed':20260930,'dt_s':2,'mapping_method':str(op['mapping_method']) if 'mapping_method' in op else 'not recorded','B_K2':B.tolist(),'w_um_per_K':w.tolist(),'studies':{}}
 arrays={'T150':T,'B':B,'w':w}
 for kind,rows,steps in [('fixed',fixed,[0,2,3,4]),('seq',seq,[0,1,2,3,4])]:
  stages=[]
  for k in steps:
   H=np.array([r[1] for r in rows[:k]]).reshape(k,5);R=np.diag([r[2] for r in rows[:k]])
   if k:
    A=B@H.T;S=H@A+R;P=B-A@np.linalg.solve(S,A.T)
   else:P=B.copy()
   Pi=B.copy()
   for _,hc,rc in rows[:k]:Pi=add_row(Pi,hc,rc)
   assert np.allclose(P,Pi,rtol=1e-9,atol=1e-10)
   assert np.linalg.eigvalsh(P).min()>-1e-9
   variance=float(w@P@w);sig=np.sqrt(max(variance,0));label='測る前' if not k else ('温度P2・P0' if kind=='fixed' else '＋'.join(r[0] for r in rows[:min(k,2)]))
   if k>2:label+=f'\n＋変位{k-2}点'
   stages.append({'label':label,'n_obs':k,'H':H.tolist(),'R':R.tolist(),'P_K2':P.tolist(),'temperature_variance_K2':np.diag(P).tolist(),'sigma2_um2':variance,'sigma_um':sig})
  out['studies'][kind]=stages
  arrays[f'P_{kind}']=np.array([st['P_K2'] for st in stages])
  fig,axs=plt.subplots(1,2,figsize=(18,7));x=np.arange(len(stages));cols=['#888888','#0072b2','#56b4e9','#e69f00','#cc3311'][:len(stages)]
  plotlabels=['測る前','P2・P0','P2・P0\n＋旧S1','P2・P0\n＋旧S1・旧S2'] if kind=='fixed' else ['測る前','T1','T1・T2','T1・T2\n＋S1','T1・T2\n＋S1・S2']
  vals=[st['sigma2_um2'] for st in stages];axs[0].bar(x,vals,color=cols)
  for i,st in enumerate(stages):axs[0].text(i,vals[i],f"分散 {vals[i]:.3f}\nσ={st['sigma_um']:.3f} µm",ha='center',va='bottom',fontsize=15)
  axs[0].set_ylim(0,max(vals)*1.35);axs[0].set_xticks(x,plotlabels,fontsize=16);axs[0].set_ylabel('反りの分散 σ² [µm²]',fontsize=18);axs[0].set_title('対象：未観測の反り Uz(A)−Uz(O)',fontsize=20)
  for j in range(5):axs[1].plot(x,[st['temperature_variance_K2'][j] for st in stages],'-o',lw=2.3,label=f'P{j}')
  axs[1].set_xticks(x,plotlabels,fontsize=16);axs[1].set_ylabel('代表点温度の分散：Pの対角成分 [K²]',fontsize=18);axs[1].set_title('測定後に残る5点温度の不確かさ',fontsize=20);axs[1].legend(loc='upper right',ncol=3,fontsize=15)
  for ax in axs:ax.grid(axis='y',alpha=.25);ax.tick_params(axis='y',labelsize=15)
  title='blog008：温度P2・P0を固定して変位を追加' if kind=='fixed' else 'blog009：選定温度T1・T2の後に選定変位S1・S2を追加'
  fig.suptitle(title+'（t=150秒で再計算）',fontsize=23)
  fig.text(.5,.035,'4,000ケース・seed=20260930・RK4 Δt=2秒・15 W加熱中／観測ノイズ：温度0.3 K、変位0.3 µm\n線形推定で残る分散の予測。実測値は使わず、実際のEnKF誤差ではありません。',ha='center',fontsize=15)
  fig.tight_layout(rect=[0,.12,1,.93]);fig.savefig(IMG/f'posterior_t150_{kind}.png',dpi=150);plt.close(fig)
  print(kind,[(st['n_obs'],round(st['sigma_um'],5)) for st in stages],flush=True)
 out['elapsed_s']=time.perf_counter()-start
 np.savez_compressed(RES/'posterior_t150.npz',**arrays)
 (RES/'posterior_t150.json').write_text(json.dumps(out,ensure_ascii=False,indent=2))
 print('elapsed',out['elapsed_s'])
if __name__=='__main__':main()
