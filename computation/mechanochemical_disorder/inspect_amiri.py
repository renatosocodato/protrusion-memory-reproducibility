#!/usr/bin/env python3
from pathlib import Path
from reproduction_io import amiri_workbook, output_directory
import hashlib,json,numbers
import numpy as np,openpyxl
from scipy.stats import kendalltau
P=output_directory(Path(__file__));src=amiri_workbook();out=P/'amiri_summary.json'
if out.exists():raise FileExistsError('Preserve analysis')
hash=hashlib.sha256(src.read_bytes()).hexdigest();assert hash=='dcf21574fbcd5318d1581b3e1c80bd77893c9214dd871a71f174bda49f23b907'
sheet=openpyxl.load_workbook(src,read_only=True,data_only=True)['Fig 2'];rows=list(sheet.iter_rows(values_only=True));outcomes={}
for name,j in [('sensor_calyculin',0),('sensor_blebbistatin',5),('control_calyculin',10)]:
 arr=np.array([[float(v) for v in row[j:j+3]] for row in rows[2:] if all(isinstance(v,numbers.Real) and not isinstance(v,bool) for v in row[j:j+3])]);assert len(arr)>=3 and np.isfinite(arr).all() and np.all(np.diff(arr[:,0])>0)
 t,T,E=arr.T;k=kendalltau(T,E).statistic
 outcomes[name]={'n_timepoints_one_cell':len(t),'time_minutes':t.tolist(),'gaps_minutes':np.diff(t).tolist(),'first':{'time_min':float(t[0]),'traction_pa':float(T[0]),'fret_efficiency':float(E[0])},'last':{'time_min':float(t[-1]),'traction_pa':float(T[-1]),'fret_efficiency':float(E[-1])},'first_last_change':{'traction_pa':float(T[-1]-T[0]),'fret_efficiency':float(E[-1]-E[0])},'first_three_medians':{'traction_pa':float(np.median(T[:3])),'fret_efficiency':float(np.median(E[:3]))},'last_three_medians':{'traction_pa':float(np.median(T[-3:])),'fret_efficiency':float(np.median(E[-3:]))},'within_cell_kendall_tau':float(k),'values':arr.tolist()}
result={'source_doi':'10.1038/s41467-023-43612-5','source_sha256':hash,'cell_type':'U2OS, nonmicroglial','source':'public Source Data, Fig 2 representative traces c,d,f','post_hoc_exploratory':True,'n_independent_cells_per_curve':1,'conditions':outcomes,'interpretation':'FRET sensor changes and substrate traction are paired within each example; no N-WASP gate, actin-assembly time series, or local cellular force calibration.'}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
