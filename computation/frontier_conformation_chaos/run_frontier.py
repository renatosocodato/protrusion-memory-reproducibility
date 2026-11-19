#!/usr/bin/env python3
"""Frozen independent synthetic perturbation and optical-reference study."""
import ast, json, math
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
from scipy.linalg import expm
ROOT=output_directory(Path(__file__))
OUT=ROOT/'results_v1'
if OUT.exists(): raise FileExistsError('Consumed run preserved')
OUT.mkdir()
source=Path(__file__).resolve().parent.parent/'joint_fret_actin_mechanics'/'run_joint_observation.py'
tree=ast.parse(source.read_text())
names={'streams','draw','rhs','sim','observations','feature','fit','predict','wilson','rep'}
functions=ast.Module(body=[x for x in tree.body if isinstance(x,ast.FunctionDef) and x.name in names],type_ignores=[])
T=np.array([0,20,35,50,65,80,95,110,140,155,170,185,200,215,230,260,300]); TX=np.array([0,60,120,180,240,300]); TALL=sorted(set(T)|set(TX))
exec(compile(functions,str(source),'exec'),globals())
# Exact reconstruction of the previous development cohort; no holdout is reused.
trrng,_=streams(202611180101); ptrain=draw(trrng[0],300)
trainstates,_=sim(ptrain,'both'); train,_=observations(ptrain,trainstates,trrng[3])
models={kind:fit(feature(train,kind)) for kind in ('joint','x')}
seed=202611180301
ss=np.random.SeedSequence(seed); r=[np.random.default_rng(s) for s in ss.spawn(10)]
p=draw(r[0],200)
wide={k:v.copy() for k,v in p.items()}
for key,a,b in [('alpha',-.12,.12),('beta',-.6,.6),('chi',-.12,.12),('bleed',0,.20),('k',0,.15),('j',0,.15)]: wide[key]=r[1].uniform(a,b,(200,1))

def simulate(p,reverse,mix):
 n=len(p['N']); y=np.full((n,2,5),.5); out={0:y.copy()}; dt=.25; ids={round(t/dt):t for t in TALL}
 def deriv(y,u,v):
  result=rhs(y,p,u,v,False)
  q,a,_,s,_=np.moveaxis(y,-1,0)
  result[:,:,3]=((1-mix)*q+mix*a-s)/p['ts']
  return result
 for k in range(1200):
  t=(k+.5)*dt
  first=float(20<=t<80); second=float(140<=t<200)
  u,v=(second,first) if reverse else (first,second)
  k1=deriv(y,u,v); k2=deriv(y+.5*dt*k1,u,v); k3=deriv(y+.5*dt*k2,u,v); k4=deriv(y+dt*k3,u,v)
  y=y+dt*(k1+2*k2+2*k3+k4)/6
  if not np.isfinite(y).all() or (y[:,:,0]<0).any() or (y[:,:,0]>1).any() or (y[:,:,1:3]<0).any() or (y[:,:,3:]<0).any(): raise ValueError('nonphysical state')
  if k+1 in ids:out[ids[k+1]]=y.copy()
 return out
conditions=[('same',p,False,0),('reverse',p,True,0),('wide',wide,False,.5),('reverse_wide',wide,True,.5)]
result={}; predictions={}; floors={}
for i,(name,params,reverse,mix) in enumerate(conditions):
 state=simulate(params,reverse,mix); obs,fl=observations(params,state,r[i+2]);floors[name]=fl
 if any(fl.values()): raise ValueError('denominator floor')
 predictions[name]={k:predict(m,feature(obs,k)) for k,m in models.items()}
 result[name]={k:rep(v) for k,v in predictions[name].items()}
idx=r[6].integers(0,200,(5000,200)); labs=np.array([0,1]); comparisons={}
for name in ('reverse','wide','reverse_wide'):
 comparisons[name]={}
 for kind in models:
  a=(predictions[name][kind]==labs).mean(axis=1); b=(predictions['same'][kind]==labs).mean(axis=1);d=a-b
  comparisons[name][kind]={'difference_percentage_points':float(d.mean()*100),'parent_bootstrap95_percentage_points':[float(v*100) for v in np.quantile(d[idx].mean(axis=1),[.025,.975])]}
# Idealized calibration, with standards and heldout mixture generated independently.
M=np.array([[1.,.15],[.25,1.]])
levels=r[7].uniform(500,3000,(2,60)); cols=[]
for c in range(2):
 y=r[7].poisson(40+levels[c,:,None]*M[:,c][None,:]); cols.append(((levels[c,:,None]*(y-40)).sum(axis=0)/(levels[c]@levels[c])))
Mhat=np.column_stack(cols)
xtrue=r[8].uniform(500,3000,(100,2)); y=r[8].poisson(40+xtrue@M.T); recovered=np.linalg.solve(Mhat,(y-40).T).T
if (recovered<=0).any() or np.linalg.cond(Mhat)>100:raise ValueError('reference inversion failed')
true=xtrue[:,1]/xtrue[:,0]; naive=(y[:,1]-40)/(y[:,0]-40);corrected=recovered[:,1]/recovered[:,0]
optical={'true_matrix':M.tolist(),'fitted_matrix':Mhat.tolist(),'training_standards_per_component':60,'heldout_mixture_standards':100,'median_absolute_relative_ratio_error_naive':float(np.median(abs(naive/true-1))),'median_absolute_relative_ratio_error_corrected':float(np.median(abs(corrected/true-1))),'known_blank_and_known_component_levels':'assumed, synthetic'}
# Analytic conditional spectrum check, exact linear Jacobian and e^At.
tau=np.array([8.,25.,20.,4.,5.]); A=np.diag(-1/tau); A[1,0]=1/tau[1];A[2,1]=1/tau[2];A[3,0]=1/tau[3];A[4,1]=1/tau[4]
eigs=np.linalg.eigvals(A);diag=-1/tau
if not np.allclose(np.sort(eigs),np.sort(diag)): raise ValueError('eigenvalues')
P=expm(A*300); analytic={'taus_assigned':tau.tolist(),'eigenvalues':eigs.tolist(),'largest_conditional_exponent':float(max(eigs)),'perturbation_propagator_frobenius_norm_t300':float(np.linalg.norm(P)),'scope':'original feed-forward cascade under identical prescribed inputs'}
output={'synthetic_only':True,'training_seed_reconstructed':202611180101,'independent_root_seed':seed,'n_parents':200,'parent_paired_two_families':True,'conditions':result,'relative_to_same_schedule':comparisons,'optical_reference_calibration':optical,'conditional_stability':analytic,'denominator_floors':floors}
(OUT/'results.json').write_text(json.dumps(output,indent=2)+'\n')
np.savez_compressed(OUT/'predictions.npz',**{f'{n}_{k}':v for n,by in predictions.items() for k,v in by.items()})
print(json.dumps(output,indent=2))
