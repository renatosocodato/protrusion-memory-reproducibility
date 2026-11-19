from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
import pandas as pd
from run_spatial import draw,simulate,wilson,MODES,NTEST
P=output_directory(Path(__file__))/'results_v1'
if (P/'second_challenge.json').exists():raise SystemExit('Preserve previously generated challenge')
rows=[];out={}
for m,seed in ((.3,28112026),(1.,28112027)):
 rng=np.random.default_rng(seed); par,ph=draw(rng,NTEST); regime={}
 for kind in ('D','R','A'):
  s=simulate(kind,par,ph,m)
  if s['mass_err']>1e-10 or s['min_c']<=0 or (kind=='D' and s['max_energy_increase']>1e-9):raise RuntimeError('invalid physical state')
  obs=s['phi']+rng.normal(0,.002,s['phi'].shape)
  ft=np.fft.rfft(obs,axis=2); amps=np.abs(ft[:,:,MODES]);g=np.log(amps[:,1,:]/amps[:,0,:])/15
  S=g[:,1]-g[:,0]
  forceobs=s['F'][:,[0,3]]+rng.normal(0,.0002,(NTEST,2)); dF=forceobs[:,1]-forceobs[:,0]
  spectral=S>0;joint=spectral&((dF>0) if m==.3 else (dF<0))
  def entry(bits):
   k=int(bits.sum()); return {'predicted_demixing':k,'n':NTEST,'fraction':k/NTEST,'wilson95':wilson(k,NTEST)}
  regime[kind]={'spectral':entry(spectral),'spectral_plus_force':entry(joint),'S_median':float(np.median(S)),'S_quantiles_025_975':np.quantile(S,[.025,.975]).tolist(),'force_change_measured_median':float(np.median(dF)),'force_change_measured_quantiles_025_975':np.quantile(dF,[.025,.975]).tolist()}
  for i in range(NTEST):rows.append({'m':m,'parent':i,'model':kind,'S':S[i],'delta_force_measured':dF[i],'spectral_positive':bool(spectral[i]),'joint_positive':bool(joint[i])})
 out[str(m)]=regime
pd.DataFrame(rows).to_csv(P/'second_challenge_parents.csv',index=False)
(P/'second_challenge.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
