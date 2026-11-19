#!/usr/bin/env python3
import json
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
from scipy.integrate import solve_ivp
p=output_directory(Path(__file__))/'dop853_check.json'
if p.exists():raise FileExistsError('Preserve check')
def rhs(t,z):
 U,V,W,x=z[:4]; T=z[4:].reshape(4,4);r=28.;eta=.5
 m=1/(1+np.exp(-(W-25)/30)); y=np.array([10*(V-U)+eta*np.tanh(U/10),U*(r-W)-V,U*V-8*W/3,(m-x)/2])
 J=np.array([[-10+eta*(1-np.tanh(U/10)**2)/10,10,0,0],[r-W,-1,-U,0],[V,U,-8/3,0],[0,0,m*(1-m)/60,-.5]])
 return np.r_[y,(J@T).ravel()]
z=np.r_[[1.,1.,1.,.5],np.eye(4).ravel()];kw={'method':'DOP853','rtol':1e-9,'atol':1e-11,'max_step':.1}
t=solve_ivp(rhs,(0,200),z,**kw)
if not t.success:raise ValueError(t.message)
z=np.r_[t.y[:4,-1],np.eye(4).ravel()];tot=np.zeros(4);tr=0.
for i in range(2000):
 # Integral of trace with midpoint quadrature at 0.5 unit intervals, while integrating tangent accurately.
 sol=solve_ivp(rhs,(0,.5),z,t_eval=[.25,.5],**kw)
 if not sol.success:raise ValueError(sol.message)
 mid=sol.y[:,0];tr+=.5*(-10+.05*(1-np.tanh(mid[0]/10)**2)-1-8/3-.5)
 z=sol.y[:,-1];Q,R=np.linalg.qr(z[4:].reshape(4,4));tot+=np.log(np.abs(np.diag(R)));z[4:]=Q.ravel()
ex=tot/1000;trace=tr/1000;ref=.903587316347384
out={'method':'DOP853','rtol':1e-9,'atol':1e-11,'max_step':.1,'qr_interval':.5,'measurement_time':1000,'transient':200,'spectrum':ex.tolist(),'mean_trace_midpoint_approx':float(trace),'trace_difference':float(ex.sum()-trace),'rk4_reference_maximum_exponent':ref,'largest_exponent_difference':float(abs(ex[0]-ref)),'check_pass':bool(abs(ex[0]-ref)<.05 and abs(ex.sum()-trace)<.02)}
p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
