from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
import pandas as pd
p=output_directory(Path(__file__))/'results_v1'
d=pd.read_csv(p/'heldout_predictions.csv')
rng=np.random.default_rng(241117)
out={}
for m in (.3,1.):
 a=d[d.mean_concentration==m]
 xs=[]
 for mode in ('position','position_plus_organization'):
  x=a[a.channel==mode].pivot(index='parent_id',columns='true_class',values='predicted_class').sort_index()
  x=x[['uniform','finite_recruitment','demixing_analogue']]
  xs.append(x.to_numpy())
 assert all(len(x)==300 for x in xs)
 ix=rng.integers(0,300,(5000,300))
 true=np.array(['uniform','finite_recruitment','demixing_analogue'])
 boot_acc=[]; boot_fpr=[]
 for x in xs:
  sampled=x[ix]; correct=(sampled==true)
  boot_acc.append(correct.mean(axis=(1,2)))
  boot_fpr.append((sampled[:,:,:2]=='demixing_analogue').mean(axis=(1,2)))
 gain=boot_acc[1]-boot_acc[0]
 out[str(m)]={'accuracy_gain_organization_minus_position':{'estimate':float(np.mean(xs[1]==true)-np.mean(xs[0]==true)),'cluster_bootstrap_95':[float(z) for z in np.quantile(gain,[.025,.975])]},'false_C_among_A_B':{mode:{'estimate':float(np.mean(xs[i][:,:2]=='demixing_analogue')),'cluster_bootstrap_95':[float(z) for z in np.quantile(boot_fpr[i],[.025,.975])]} for i,mode in enumerate(('position','position_plus_organization'))}}
(p/'paired_uncertainty.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
