from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
import pandas as pd
from run_spatial import draw,initial,h,wilson,L,MODES,NTEST
P=output_directory(Path(__file__))/'results_v1'
if (P/'third_challenge.json').exists():raise SystemExit('Preserve existing challenge')
def run(par,phases,gamma,m,n=64,dt=.25):
 eig=4*np.sin(np.pi*np.arange(n//2+1)/n)**2/(L/n)**2
 M,kap,_,_,_,tb,tx,_=par.T
 phi=initial(n,par,phases);b=np.full((len(par),n),.7*h(m));pos=np.full(len(par),.7*h(m));
 phis=[phi.copy()];Fs=[b.mean(axis=1)];Xs=[pos.copy()];mass=0.;minc=10.
 for j in range(round(60/dt)):
  cube=phi**3;cube-=cube.mean(axis=1,keepdims=True)
  z=(np.fft.rfft(phi,axis=1)-dt*gamma[:,None]*np.fft.rfft(cube,axis=1))/(1+dt*M[:,None]*(kap[:,None]*eig**2-eig))
  phi=np.fft.irfft(z,n=n,axis=1); c=m*(1+.6*phi)
  eb=np.exp(-dt/tb);ex=np.exp(-dt/tx)
  b=b*eb[:,None]+.7*h(c)*(1-eb[:,None]);F=b.mean(axis=1);pos=pos*ex+F*(1-ex)
  mass=max(mass,float(np.max(np.abs(phi.mean(axis=1)))));minc=min(minc,float(np.min(c)))
  if not np.all(np.isfinite(phi)) or minc<=0 or mass>1e-10:raise RuntimeError('inadmissible')
  t=(j+1)*dt
  if t in (15.,30.,60.):phis.append(phi.copy());Fs.append(F.copy());Xs.append(pos.copy())
 return {'phi':np.stack(phis,axis=1),'F':np.stack(Fs,axis=1),'x':np.stack(Xs,axis=1),'mass_error':mass,'min_c':minc}
q=np.array([[.2,1.,.15,.15,.065,15.,20.,.03]]);ph=np.array([[.2,1.1,2.2]]);gamma=np.array([.1])
a=run(q,ph,gamma,.3);b=run(q,ph,gamma,.3,dt=.125);c=run(q,ph,gamma,.3,n=128)
timephi=float(np.max(abs(a['phi']-b['phi'])));gridphi=float(np.max(abs(a['phi']-c['phi'][:,:,::2])))
timeforce=max(float(np.max(abs(a[k]-b[k]))) for k in ('F','x'));gridforce=max(float(np.max(abs(a[k]-c[k]))) for k in ('F','x'))
checks={'dt_phi':timephi,'grid_phi':gridphi,'dt_force_deformation':timeforce,'grid_force_deformation':gridforce}
(P/'third_checks.json').write_text(json.dumps(checks,indent=2))
if timephi>.02 or gridphi>.02 or timeforce>.002 or gridforce>.002:raise SystemExit('New-operator check failed; stop')
out={};rows=[]
for m,seed in ((.3,29112026),(1.,29112027)):
 rng=np.random.default_rng(seed);par,ph=draw(rng,NTEST);gamma=rng.uniform(.05,.15,NTEST);s=run(par,ph,gamma,m)
 ob=s['phi']+rng.normal(0,.002,s['phi'].shape);ff=np.fft.rfft(ob,axis=2);v=np.abs(ff[:,:,MODES]);g=np.log(v[:,1,:]/v[:,0,:])/15;S=g[:,1]-g[:,0]
 fo=s['F'][:,[0,3]]+rng.normal(0,.0002,(NTEST,2));dF=fo[:,1]-fo[:,0];one=S>0;two=one&((dF>0) if m==.3 else (dF<0))
 def rec(mask):
  k=int(mask.sum());return {'false_positive':k,'n':NTEST,'fraction':k/NTEST,'wilson95':wilson(k,NTEST)}
 out[str(m)]={'spectral':rec(one),'spectral_plus_force':rec(two),'S_median':float(np.median(S)),'S_quantiles_025_975':np.quantile(S,[.025,.975]).tolist(),'force_change_measured_median':float(np.median(dF)),'mass_error':s['mass_error'],'min_c':s['min_c']}
 for i in range(NTEST):rows.append({'m':m,'parent':i,'S':S[i],'force_change':dF[i],'spectral_positive':bool(one[i]),'joint_positive':bool(two[i])})
pd.DataFrame(rows).to_csv(P/'third_challenge_parents.csv',index=False)
(P/'third_challenge.json').write_text(json.dumps(out,indent=2));print(json.dumps({'checks':checks,'results':out},indent=2))
