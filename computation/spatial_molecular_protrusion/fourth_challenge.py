from pathlib import Path
from reproduction_io import output_directory
import json,time
import numpy as np
import pandas as pd
from run_spatial import draw,h,wilson,L,NTEST
P=output_directory(Path(__file__))/'results_v1'
if (P/'fourth_challenge.json').exists():raise SystemExit('Preserve previous result')
AMPS=np.array([.12,.28,.4]);N=64

def simulate(par,phase,beta,m,kind,n=64,dt=.05):
 batch=len(par); x=np.arange(n)*L/n;k=np.arange(n//2+1);eig=4*np.sin(np.pi*k/n)**2/(L/n)**2
 M,kap,_,_,_,tb,tx,_=par.T;gamma=M*beta*eig[6]
 phi=AMPS[None,:,None]*np.cos(2*np.pi*2*x[None,None,:]/L+phase[:,None,None]); phi0=phi.copy()
 b=np.full((batch,3,n),.7*h(m));pos=np.full((batch,3),.7*h(m));F0=b.mean(axis=2)
 e_b=np.exp(-dt/tb)[:,None,None];e_x=np.exp(-dt/tx)[:,None]
 minc=10.;mass=0.;maxE=-np.inf
 def energy(q):
  grad=(np.roll(q,-1,axis=2)-q)/(L/n)
  return (beta[:,None,None]*q**4/4-q*q/2+kap[:,None,None]*grad**2/2).sum(axis=2)*L/n
 if kind=='D':Eold=energy(phi)
 for _ in range(round(5/dt)):
  cube=phi**3;z=np.fft.rfft(phi,axis=2);zc=np.fft.rfft(cube,axis=2)
  if kind=='D':z=(z-dt*M[:,None,None]*beta[:,None,None]*eig*zc)/(1+dt*M[:,None,None]*(kap[:,None,None]*eig**2-eig))
  else:
   cube-=cube.mean(axis=2,keepdims=True);zc=np.fft.rfft(cube,axis=2)
   z=(z-dt*gamma[:,None,None]*zc)/(1+dt*M[:,None,None]*(kap[:,None,None]*eig**2-eig))
  phi=np.fft.irfft(z,n=n,axis=2);c=m*(1+.6*phi);b=b*e_b+.7*h(c)*(1-e_b);F=b.mean(axis=2);pos=pos*e_x+F*(1-e_x)
  if not np.all(np.isfinite(phi)):raise RuntimeError('nonfinite')
  minc=min(minc,float(np.min(c)));mass=max(mass,float(np.max(abs(phi.mean(axis=2)))))
  if kind=='D':
   E=energy(phi);maxE=max(maxE,float(np.max(E-Eold)));Eold=E
 if minc<=0 or mass>1e-10 or (kind=='D' and maxE>1e-9):raise RuntimeError('invariant failure')
 return {'phi0':phi0,'phi5':phi,'F0':F0,'F5':F,'x5':pos,'mass':mass,'minc':minc,'max_energy_inc':maxE}

def observe(s,phase,rng):
 ob0=s['phi0']+rng.normal(0,.002,s['phi0'].shape);ob5=s['phi5']+rng.normal(0,.002,s['phi5'].shape)
 z0=np.fft.rfft(ob0,axis=2);z5=np.fft.rfft(ob5,axis=2);g=np.log(abs(z5[:,:,2])/abs(z0[:,:,2]))/5
 w=AMPS**2-(AMPS**2).mean();C1=-(g*w).sum(axis=1)/(w*w).sum()
 a3=2*np.real((z5[:,2,6]-z0[:,2,6])*np.exp(-3j*phase))/N
 C3=-a3/(5*AMPS[2]**3); Q=C3/C1
 return C1,C3,Q

q=np.array([[.2,1.,.15,.15,.065,15.,20.,.03]]);phase=np.array([.3]);beta=np.array([1.])
checks={}
for kind in ('D','H'):
 a=simulate(q,phase,beta,.3,kind);b=simulate(q,phase,beta,.3,kind,dt=.025);c=simulate(q,phase,beta,.3,kind,n=128)
 difft=float(np.max(abs(a['phi5']-b['phi5'])));diffs=float(np.max(abs(a['phi5']-c['phi5'][:,:,::2])));ft=float(np.max(abs(a['F5']-b['F5'])));fs=float(np.max(abs(a['F5']-c['F5'])))
 checks[kind]={'dt_phi':difft,'grid_phi':diffs,'dt_force':ft,'grid_force':fs}
(P/'fourth_checks.json').write_text(json.dumps(checks,indent=2))
if any(v['dt_phi']>.02 or v['grid_phi']>.02 or v['dt_force']>.002 or v['grid_force']>.002 for v in checks.values()):raise SystemExit('Reference comparison failed; preserved and stopped')
out={};rows=[]
for m,seed in ((.3,30112026),(1.,30112027)):
 rng=np.random.default_rng(seed);par,ph=draw(rng,NTEST);phase=ph[:,0];beta=rng.uniform(.7,1.3,NTEST);reg={}
 for kind in ('D','H'):
  s=simulate(par,phase,beta,m,kind);C1,C3,Q=observe(s,phase,rng);pred=(C1>0)&(Q>1.2);k=int(pred.sum());
  reg[kind]={'called_demixing':k,'n':NTEST,'fraction':k/NTEST,'wilson95':wilson(k,NTEST),'invalid_C1':int((C1<=0).sum()),'Q_median':float(np.median(Q)),'Q_quantiles_025_975':np.quantile(Q,[.025,.975]).tolist(),'force_change_median_largest_amplitude':float(np.median(s['F5'][:,2]-s['F0'][:,2])),'min_concentration':s['minc'],'max_mass_error':s['mass'],'max_energy_increment':s['max_energy_inc'] if kind=='D' else None}
  for i in range(NTEST):rows.append({'m':m,'parent':i,'model':kind,'C1':C1[i],'C3':C3[i],'Q':Q[i],'called_demixing':bool(pred[i]),'force_change':s['F5'][i,2]-s['F0'][i,2]})
 out[str(m)]=reg
pd.DataFrame(rows).to_csv(P/'fourth_challenge_parents.csv',index=False)
(P/'fourth_challenge.json').write_text(json.dumps(out,indent=2));print(json.dumps({'checks':checks,'results':out},indent=2))
