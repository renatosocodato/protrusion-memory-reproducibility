#!/usr/bin/env python3
"""Frozen paired synthetic gate challenge. No empirical fitting."""
from pathlib import Path
from reproduction_io import output_directory
import json, math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=output_directory(Path(__file__))
OUT=ROOT/'results_v1';OUT.mkdir(exist_ok=True)
N=200

def draw(seed,n):
 r=np.random.default_rng(seed)
 return np.column_stack([r.uniform(a,b,n) for a,b in [(6,12),(8,20),(10,25),(12,35),(.45,.75),(5,10),(35,70)]])

def rhs(y,par,m,U,z,tau):
 p,r,q,b,x=np.moveaxis(y,-1,0);tp,tr,tb,tx,R=par[:,:5].T
 h=lambda c:c*c/(1+c*c)
 H=(h(m*(1-r))+h(m*(1+r)))/2
 qe=1/(1+np.exp(-z))
 return np.stack(((U-p)/tp, (1-r*r)*(R*p-r)/tr,(qe-q)/tau,(q*H-b)/tb,(b-x)/tx),axis=-1)

def simulate(par,m,tau,dt=.25,history=False,offset=0.):
 n=len(par);h=m*m/(1+m*m);y=np.column_stack((np.zeros(n),np.zeros(n),np.full(n,.5),np.full(n,.5*h),np.full(n,.5*h)))
 times=(0,35,60,65,120,180,240,300); indices={round(t/dt):t for t in times};out={0:y.copy()};hist=[y.copy()] if history else None
 for k in range(round(300/dt)):
  t=(k+.5)*dt;U=float(20<=t<200);u=1.5*float(20<=t<80)
  # Offset is added to BOTH sigma and epsilon to test exact gauge.
  z=(offset-offset)+u
  k1=rhs(y,par,m,U,z,tau);k2=rhs(y+.5*dt*k1,par,m,U,z,tau);k3=rhs(y+.5*dt*k2,par,m,U,z,tau);k4=rhs(y+dt*k3,par,m,U,z,tau)
  y=y+dt*(k1+2*k2+2*k3+k4)/6
  if not np.isfinite(y).all() or (y[:,:3]<-1e-11).any() or (y[:,:3]>1+1e-11).any() or (y[:,3:]<-1e-11).any() or (m*(1-y[:,1])<=0).any(): raise ValueError('violated physical or numerical bounds at step '+str(k))
  if k+1 in indices:out[indices[k+1]]=y.copy()
  if history:hist.append(y.copy())
 return out, np.stack(hist) if history else None

def features(fast,slow,rng):
 classes=(fast,slow);x=np.stack([np.stack([a[t][:,4] for t in (0,60,120,180,240,300)],axis=1) for a in classes],axis=1)
 q=np.stack([np.stack([a[t][:,2] for t in (35,65)],axis=1) for a in classes],axis=1)
 x=x+rng.normal(0,.005,x.shape);q=q+rng.normal(0,.03,q.shape)
 return (x[:,:,1:]-x[:,:,[0]]),q

def fit(X):
 x=X.reshape(-1,X.shape[-1]);labels=np.tile([0,1],X.shape[0]);mean=x.mean(0);scale=np.maximum(x.std(0),1e-12);z=(x-mean)/scale
 centers=np.stack([z[labels==i].mean(0) for i in range(2)]);delta=z-centers[labels];cov=delta.T@delta/len(labels)+.01*np.eye(X.shape[-1]);return mean,scale,centers,np.linalg.inv(cov)

def predict(model,X):
 mean,scale,centers,inv=model;z=(X-mean)/scale;d=z[:,:,None,:]-centers[None,None,:,:];return np.einsum('...i,ij,...j->...',d,inv,d).argmin(-1)

def wilson(k,n):
 z=1.959963984540054;p=k/n;d=1+z*z/n;mid=(p+z*z/(2*n))/d;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;return [mid-h,mid+h]

def report(pred):
 correct=pred==np.array([0,1]);k=correct.sum(axis=0);return {'n_parents':len(pred),'correct_fast':int(k[0]),'correct_slow':int(k[1]),'accuracy':float(correct.mean()),'class_wilson95':[wilson(int(k[i]),len(pred)) for i in range(2)],'calls_fast_by_true_class':[int((pred[:,i]==0).sum()) for i in range(2)]}

nom=np.array([[9.,14.,17.5,23.5,.6,7.5,52.5]])
checks={}
for col,label in [(5,'fast'),(6,'slow')]:
 a,ha=simulate(nom,.3,nom[:,col],.25,True);b,hb=simulate(nom,.3,nom[:,col],.125,True);err=float(np.max(np.abs(ha-hb[::2])));checks['dt_'+label]=err
 if err>=1e-4: raise ValueError('prespecified step comparison failed')
 a,_=simulate(nom,.3,nom[:,col],.25,False,0);b,_=simulate(nom,.3,nom[:,col],.25,False,math.log(2));err=max(float(np.max(np.abs(a[t]-b[t]))) for t in a);checks['entropy_energy_twin_'+label]=err
 if err>1e-12:raise ValueError('gauge identity failed')

output={'status':'synthetic_only','n':N,'bounds':'checked during every integration step','checks':checks,'m':{}}
rows=[];fig,axs=plt.subplots(2,3,figsize=(12,6.9))
for row,m in enumerate([.3,1.]):
 pars_train=draw(2026111801+row*10,N);pars_test=draw(2026111802+row*10,N);rngtr=np.random.default_rng(2026111811+row);rngte=np.random.default_rng(2026111821+row);rngnull=np.random.default_rng(2026111831+row)
 trf,_=simulate(pars_train,m,pars_train[:,5]);trs,_=simulate(pars_train,m,pars_train[:,6]);tef,_=simulate(pars_test,m,pars_test[:,5]);tes,_=simulate(pars_test,m,pars_test[:,6]);Xtr,Qtr=features(trf,trs,rngtr);Xte,Qte=features(tef,tes,rngte)
 modelx=fit(Xtr);modelq=fit(np.concatenate((Xtr,Qtr),axis=2));px=predict(modelx,Xte);pq=predict(modelq,np.concatenate((Xte,Qte),axis=2))
 cx=px==np.array([0,1]);cq=pq==np.array([0,1]);gain=cq.mean()-cx.mean();brng=np.random.default_rng(2026111803+row);inds=brng.integers(0,N,(5000,N));g=(cq[inds].mean((1,2))-cx[inds].mean((1,2)))
 null1,_=simulate(pars_test,m,pars_test[:,5]);null2,_=simulate(pars_test,m,pars_test[:,5]);Xnull,Qnull=features(null1,null2,rngnull);nullx=predict(modelx,Xnull);nullq=predict(modelq,np.concatenate((Xnull,Qnull),axis=2))
 # Median trajectories are paired parameter-wise; they are not data.
 output['m'][str(m)]={'position_only':report(px),'position_plus_gate':report(pq),'paired_gain':float(gain),'paired_gain_bootstrap95':[float(x) for x in np.quantile(g,[.025,.975])],'matched_null_position':report(nullx),'matched_null_position_plus_gate':report(nullq),'gate_q35_median_fast':float(np.median(tef[35][:,2])),'gate_q35_median_slow':float(np.median(tes[35][:,2])),'x60_median_fast':float(np.median(tef[60][:,4])),'x60_median_slow':float(np.median(tes[60][:,4]))}
 for i in range(N):rows.append([m,i,int(px[i,0]),int(px[i,1]),int(pq[i,0]),int(pq[i,1])])
 ts=np.array([0,35,60,65,120,180,240,300]);
 for c,lab in [(tef,'Fast gate'),(tes,'Slow gate')]:
  for ax,state in [(axs[row,0],2),(axs[row,1],4)]:ax.plot(ts,[np.median(c[t][:,state]) for t in ts],marker='o',ms=3,label=lab)
 ax=axs[row,2];v=[cx.mean(),cq.mean()];ax.bar(['Position','Position + gate'],v,color=['#6b7b8b','#237a64']);ax.axhline(.5,color='black',linestyle=':',lw=1);ax.set_ylim(0,1.05);ax.set_ylabel('Held-out accuracy');axs[row,0].set_title(f'm = {m}: accessible fraction');axs[row,1].set_title(f'm = {m}: effective deformation')
 for ax in axs[row,:2]:ax.set_xlabel('Assigned time (s)');ax.set_ylabel('Dimensionless');ax.axvspan(20,80,color='orange',alpha=.12)
axs[0,0].legend(frameon=False,fontsize=8);fig.suptitle('Synthetic conformational gate: conditional observation test',fontsize=13);fig.tight_layout();fig.savefig(OUT/'ensemble_gate.png',dpi=180);fig.savefig(OUT/'ensemble_gate.svg');plt.close(fig)
(OUT/'results.json').write_text(json.dumps(output,indent=2)+'\n');np.savetxt(OUT/'heldout_predictions.csv',np.array(rows),delimiter=',',header='m,parent,pos_fast,pos_slow,posq_fast,posq_slow',comments='',fmt=['%.1f','%d','%d','%d','%d','%d']);print(json.dumps(output,indent=2))
