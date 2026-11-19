from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from run_comparison import rates, simulate
P=output_directory(Path(__file__))
q=np.array([[9.,8.,14.,18.,24.,.6]])
records={}
for m in (.3,1.):
    _, states=simulate(q,m,'demixing_analogue',dt=.5,full=True)
    states=states[0]; t=np.arange(0,300.01,.5)
    p,a,r,b,x=states.T
    pR=.6
    target=r-r*(r*r-pR*pR*p)/(pR*pR)
    spline=CubicSpline(t,target)
    def rhs(t,y):
        pp,aa,rr,bb,xx=y
        op=rates(np.array([y]),t,q,m,'uniform',False)[0]
        op[2]=-(1-rr*rr)*(rr-spline(t))/14.
        return op
    y0=states[0]
    sol=solve_ivp(rhs,(0,300),y0,rtol=2e-9,atol=1e-11,t_eval=t,max_step=.25)
    if not sol.success: raise RuntimeError(sol.message)
    diff=np.abs(sol.y.T-states)
    records[str(m)]={'target_range':[float(np.min(target)),float(np.max(target))], 'max_abs_all_times_r_b_x':float(diff[:,2:].max()),'max_abs_at_minute_r_b_x':float(diff[::120,2:].max()),'initial_zero_noise_state':[float(i) for i in y0]}
(P/'results_v1'/'exploratory_matched_recruitment_v2.json').write_text(json.dumps(records,indent=2))
print(json.dumps(records,indent=2))
