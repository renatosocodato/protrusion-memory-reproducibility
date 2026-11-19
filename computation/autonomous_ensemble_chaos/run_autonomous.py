#!/usr/bin/env python3
"""Frozen autonomous benchmark: RK4 states/variational QR, no biological fitting."""
import json, math
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
from numba import njit
ROOT=output_directory(Path(__file__));OUT=ROOT/'results_v1'
if OUT.exists(): raise FileExistsError('Preserve original outcomes')
OUT.mkdir()
@njit
def fun(y,r,eta):
 U,V,W,x=y; m=1/(1+np.exp(-(W-25)/30))
 return np.array([10*(V-U)+eta*np.tanh(U/10),U*(r-W)-V,U*V-(8/3)*W,(m-x)/2])
@njit
def jac(y,r,eta):
 U,V,W,x=y; m=1/(1+np.exp(-(W-25)/30));th=np.tanh(U/10)
 J=np.zeros((4,4));J[0,0]=-10+eta*(1-th*th)/10;J[0,1]=10;J[1,0]=r-W;J[1,1]=-1;J[1,2]=-U;J[2,0]=V;J[2,1]=U;J[2,2]=-8/3;J[3,2]=m*(1-m)/60;J[3,3]=-.5
 return J
@njit
def state_step(y,h,r,eta):
 k1=fun(y,r,eta);k2=fun(y+h*k1/2,r,eta);k3=fun(y+h*k2/2,r,eta);k4=fun(y+h*k3,r,eta)
 return y+h*(k1+2*k2+2*k3+k4)/6
@njit
def joint_step(y,Q,h,r,eta):
 a=fun(y,r,eta);A=jac(y,r,eta)@Q
 y2=y+h*a/2;Q2=Q+h*A/2;b=fun(y2,r,eta);B=jac(y2,r,eta)@Q2
 y3=y+h*b/2;Q3=Q+h*B/2;c=fun(y3,r,eta);C=jac(y3,r,eta)@Q3
 y4=y+h*c;Q4=Q+h*C;d=fun(y4,r,eta);D=jac(y4,r,eta)@Q4
 return y+h*(a+2*b+2*c+d)/6, Q+h*(A+2*B+2*C+D)/6
@njit
def calculate(init,r,eta,h,record):
 y=init.copy();transient=round(200/h);n=round(1000/h);freq=round(.1/h);block=round(100/h)
 for k in range(transient):y=state_step(y,h,r,eta)
 Q=np.eye(4);cum=np.zeros(4);bcum=np.zeros(4);tracesum=0.;blocks=np.zeros((10,4));samples=np.zeros((n//freq+1,4)) if record else np.zeros((1,4));ylo=np.full(4,np.inf);yhi=np.full(4,-np.inf);si=0
 for k in range(n):
  y,Q=joint_step(y,Q,h,r,eta)
  if not np.all(np.isfinite(y)) or not (0<y[3]<1):raise ValueError('finite state or bounded output failed')
  for j in range(4):ylo[j]=min(ylo[j],y[j]);yhi[j]=max(yhi[j],y[j])
  tracesum+=np.trace(jac(y,r,eta))*h
  if (k+1)%freq==0:
   Q,R=np.linalg.qr(Q);d=np.log(np.abs(np.diag(R)));cum+=d;bcum+=d
   if record:samples[si]=y;si+=1
  if (k+1)%block==0:
   blocks[(k+1)//block-1]=bcum/100;bcum=np.zeros(4)
 return cum/1000,tracesum/1000,blocks,ylo,yhi,samples[:si] if record else samples
cases=[('canonical',28.,0.,[1,1,1,.5],.01),('feedback',28.,.5,[1,1,1,.5],.01),('initial_b',28.,.5,[-8,8,27,.5],.01),('initial_c',28.,.5,[5,-3,20,.5],.01),('half_step',28.,.5,[1,1,1,.5],.005),('low_control',10.,.5,[1,1,1,.5],.01)]
results={}
for name,r,eta,initial,h in cases:
 ex,tr,blocks,lo,hi,samples=calculate(np.array(initial,dtype=float),r,eta,h,name=='feedback')
 if abs(ex.sum()-tr)>.02:raise ValueError('volume balance mismatch: '+name)
 results[name]={'r':r,'eta':eta,'initial':initial,'dt':h,'spectrum':ex.tolist(),'mean_jacobian_trace':float(tr),'spectrum_sum_minus_mean_trace':float(ex.sum()-tr),'block_spectra':blocks.tolist(),'state_min':lo.tolist(),'state_max':hi.tolist()}
 if name=='feedback':np.savez_compressed(OUT/'representative_trace.npz',t=np.arange(1,len(samples)+1)*.1,states=samples)
base=results['feedback']['spectrum'];check={'positive_all_three':all(results[n]['spectrum'][0]>.1 for n in ['feedback','initial_b','initial_c']),'near_zero_second_all_three':all(abs(results[n]['spectrum'][1])<.05 for n in ['feedback','initial_b','initial_c']),'negative_remaining_all_three':all(max(results[n]['spectrum'][2:])<0 for n in ['feedback','initial_b','initial_c']),'step_difference_lt_0.05':abs(base[0]-results['half_step']['spectrum'][0])<.05}
res={'synthetic_only':True,'autonomous':True,'canonical_lorenz_derived':True,'assigned_units':True,'cases':results,'frozen_numerical_checks':check,'all_checks_pass':all(check.values()),'cautions':'Finite-time numerical estimate under a constructed system, not proof for cells or a mathematical proof of an attractor.'}
(OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps(res,indent=2))
