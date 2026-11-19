import json
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
out=output_directory(Path(__file__))/'feedback_boundary_results.json'
if out.exists():raise FileExistsError('Do not overwrite a consumed case')
g=40.;rows=[]
for q in (.1,.5):
 z=float(np.log(q/(1-q))+g*q); k=g*q*(1-q)
 J=np.array([[-1.,0.,-k],[1.,-1.,0.],[0.,1.,-1.]])
 eig=np.linalg.eigvals(J); spectral=float(np.max(eig.real)); theory=float(-1+k**(1/3)/2)
 assert abs(spectral-theory)<1e-12
 rows.append({'q':q,'z0':z,'feedback_gain':k,'eigenvalues':[[float(v.real),float(v.imag)] for v in eig],'largest_real_part':spectral,'analytic_largest_real_part':theory,'linear_status':'stable' if spectral<0 else 'unstable'})
result={'assumed_synthetic_feedback':True,'gain_g':g,'equal_unit_time_constants':True,'linear_instability_threshold_K':8.,'cases':rows}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
