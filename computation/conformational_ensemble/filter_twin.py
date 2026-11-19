from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
from scipy.linalg import expm

def traj(tq,tb,tx):
 A=np.array([[-1/tq,0,0],[1/tb,-1/tb,0],[0,1/tx,-1/tx]],float)
 y=np.array([.5,.5,.5],float);saved=[y.copy()]
 for t in range(300):
  f=1/(1+np.exp(-(1.5 if 20<=t<80 else 0)))
  aug=np.zeros((4,4));aug[:3,:3]=A;aug[0,3]=f/tq
  y=(expm(aug)@np.r_[y,1.])[:3];saved.append(y.copy())
 return np.array(saved)
a=traj(8,25,20);b=traj(25,8,20);r={'type':'posthoc_constructed_analytic_twin','max_deformation_abs_difference':float(np.max(np.abs(.5*a[:,2]-.5*b[:,2]))),'gate_abs_difference_t35':float(abs(a[35,0]-b[35,0])),'gate_35_values':[float(a[35,0]),float(b[35,0])],'n_timepoints':len(a),'limitations':'R=0 constant H; dimensionless output; no biological measurements; exact matrix exponential.'}
if r['max_deformation_abs_difference']>1e-12 or r['gate_abs_difference_t35']<=.05:raise ValueError('predeclared criterion failed')
(output_directory(Path(__file__)) / 'results_v1').joinpath('filter_twin.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
