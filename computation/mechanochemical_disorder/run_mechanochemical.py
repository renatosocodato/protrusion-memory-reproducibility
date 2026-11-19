#!/usr/bin/env python3
import json
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import kendalltau
P=output_directory(Path(__file__));O=P/'synthetic_v1'
if O.exists():raise FileExistsError('Preserve prior run')
O.mkdir();k0=.5;kp=1.;kd=.4;rho=.5;zeta=2.;Fp=1.;k=1.;alpha=.4;b=.2;E0=.26;g=.08;taus=1.
def equilibrium_res(a):
 q=expit(b-alpha*Fp*a);return kp*q*(1-a)*(1+rho*a)-kd*a
astar=brentq(equilibrium_res,1e-10,1-1e-10);qstar=expit(b-alpha*Fp*astar);xstar=Fp*astar/k
B=b-alpha*k*xstar;plus=k0*np.exp(B/2);minus=k0*np.exp(-B/2);h=kp*(1-astar)*(1+rho*astar)
J=np.array([[-(plus+minus),0,-alpha*k/2*(plus*(1-qstar)+minus*qstar)],[h,kp*qstar*(rho-1-2*rho*astar)-kd,0],[0,Fp/zeta,-k/zeta]])
eig=np.linalg.eigvals(J)
def rhs(y,u):
 q,a,x,s=y;bias=b+u-alpha*k*x;on=k0*np.exp(bias/2);off=k0*np.exp(-bias/2)
 return np.array([on*(1-q)-off*q,kp*q*(1-a)*(1+rho*a)-kd*a,(Fp*a-k*x)/zeta,(E0-g*k*x-s)/taus])
def run(arm,dt,seed=None):
 rng=np.random.default_rng(seed);y=np.array([.3,.1,.1,.25]);ou=0.;n=round(400/dt);sample=round(2/dt);t=[0.];states=[y.copy()];uvals=[0.]
 for i in range(n):
  mid=(i+.5)*dt
  u=0. if arm=='autonomous' else (.4*np.sin(2*np.pi*mid/30) if arm=='periodic' else ou)
  k1=rhs(y,u);k2=rhs(y+dt*k1/2,u);k3=rhs(y+dt*k2/2,u);k4=rhs(y+dt*k3,u);y=y+dt*(k1+2*k2+2*k3+k4)/6
  if not np.isfinite(y).all() or not (0<=y[0]<=1 and 0<=y[1]<=1 and 0<=y[2]<=max(.1,Fp/k)+1e-12):raise ValueError('state boundary')
  if arm=='stochastic':ou=np.exp(-dt/4)*ou+.4*np.sqrt(1-np.exp(-2*dt/4))*rng.normal()
  if (i+1)%sample==0:t.append((i+1)*dt);states.append(y.copy());uvals.append(u)
 return np.array(t),np.array(states),np.array(uvals)
result={};arrays={}
for arm in ['autonomous','periodic','stochastic']:
 t,y,u=run(arm,.02,202611180401 if arm=='stochastic' else None);mask=t>=200;v=y[mask];F=k*v[:,2];s=v[:,3]
 result[arm]={'post_transient_sample_count':len(v),'state_mean':v.mean(0).tolist(),'state_sd':v.std(0,ddof=1).tolist(),'state_min':v.min(0).tolist(),'state_max':v.max(0).tolist(),'force_sensor_kendall_tau':float(kendalltau(F,s).statistic) if np.std(F)>1e-12 else None,'input_sd':float(u[mask].std(ddof=1))};arrays[arm]=y;arrays['t']=t;arrays[arm+'_input']=u
_,reference,_=run('periodic',.01);err=float(np.max(abs(reference-arrays['periodic'])));assert err<1e-4
out={'kind':'synthetic assumed non-Lorenz model','equilibrium':{'q':qstar,'a':astar,'x':xstar,'residual':float(equilibrium_res(astar)),'eigenvalues':[[float(z.real),float(z.imag)] for z in eig],'largest_real':float(eig.real.max())},'small_gain_margin_assumed':kd-alpha*Fp*kp/2,'periodic_step_halving_max_difference':err,'seed':202611180401,'cases':result,'limitations':'No fitted biological parameters; FRET gate q is distinct from mechanically responsive s. No identification of stochastic vs driven from real data.'}
(O/'results.json').write_text(json.dumps(out,indent=2)+'\n');np.savez_compressed(O/'trajectories.npz',**arrays);print(json.dumps(out,indent=2))
