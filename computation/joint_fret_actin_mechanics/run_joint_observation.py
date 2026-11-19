#!/usr/bin/env python3
"""Frozen joint photon/readout synthetic challenge; no biological data fitting."""
from pathlib import Path
from reproduction_io import output_directory
import json, math
import numpy as np
from protrusion_memory.model import simulate_pair
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
T=np.array([0,20,35,50,65,80,95,110,140,155,170,185,200,215,230,260,300]);TX=np.array([0,60,120,180,240,300]);TALL=sorted(set(T)|set(TX))

def streams(seed):
 ss=np.random.SeedSequence(seed);children=ss.spawn(7);return [np.random.default_rng(c) for c in children], [list(c.spawn_key) for c in children]

def draw(rng,n):
 def U(a,b):return rng.uniform(a,b,(n,1))
 p={k:U(*bounds) for k,bounds in {'fast':(6,12),'slow':(30,50),'tx':(12,30),'ts':(2,6),'tl':(2,6),'alpha':(-.03,.03),'beta':(-.2,.2),'chi':(-.03,.03),'N':(1500,3000),'gd':(.85,1.15),'ga':(.85,1.15),'gl':(.85,1.15),'bd':(0,.00025),'ba':(0,.00025),'bl':(0,.00025),'bleed':(0,.04),'k':(0,.03),'j':(0,.03)}.items()};return p

def rhs(y,p,u,v,equal):
 q,a,x,s,l=np.moveaxis(y,-1,0);fast=p['fast'];slow=p['slow'];tq=np.concatenate((fast,slow),axis=1);ta=np.concatenate((slow,fast),axis=1)
 if equal:tq=ta=np.repeat(np.sqrt(fast*slow),2,axis=1)
 target=1/(1+np.exp(-(u+p['beta']*v)))
 return np.stack(((target-q)/tq,(q+p['alpha']*u+.20*v-a)/ta,(a+p['chi']*v-x)/p['tx'],(q-s)/p['ts'],(a-l)/p['tl']),axis=-1)

def sim(p,arm,dt=.25,equal=False,history=False):
 n=len(p['N']);y=np.full((n,2,5),.5);out={0:y.copy()};hist=[y.copy()] if history else None;ids={round(t/dt):t for t in TALL}
 for k in range(round(300/dt)):
  t=(k+.5)*dt;u=float(arm in ('u','both') and 20<=t<80);v=float(arm in ('v','both') and 140<=t<200)
  k1=rhs(y,p,u,v,equal);k2=rhs(y+.5*dt*k1,p,u,v,equal);k3=rhs(y+.5*dt*k2,p,u,v,equal);k4=rhs(y+dt*k3,p,u,v,equal);y=y+(dt/6)*(k1+2*k2+2*k3+k4)
  if not np.isfinite(y).all() or (y[:,:,0]<-1e-12).any() or (y[:,:,0]>1+1e-12).any() or (y[:,:,1:3]<-1e-12).any() or (y[:,:,3]<-1e-12).any() or (y[:,:,3]>1+1e-12).any() or (y[:,:,4]<-1e-12).any():raise ValueError('bound failure step '+str(k))
  if k+1 in ids:out[ids[k+1]]=y.copy()
  if history:hist.append(y.copy())
 return out,np.stack(hist) if history else None

def observations(p,series,rng):
 states=np.stack([series[int(t)] for t in T],axis=2);s=states[:,:,:,3];l=states[:,:,:,4]; tt=T[None,None,:]
 N=p['N'][:,:,None];E=.45-.20*np.tanh(2*(s-.5));dl=1-E
 def pp(k):return p[k][:,:,None]
 muD=40+pp('gd')*N*np.exp(-pp('bd')*tt)*dl;muA=40+pp('ga')*N*np.exp(-pp('ba')*tt)*(E+pp('bleed')*dl+pp('k')*l);muL=40+pp('gl')*N*np.exp(-pp('bl')*tt)*(l+pp('j')*E)
 if min(muD.min(),muA.min(),muL.min())<0 or not all(np.isfinite(z).all() for z in [muD,muA,muL]):raise ValueError('photon mean failure')
 d=rng.poisson(muD)+rng.normal(0,2,muD.shape);a=rng.poisson(muA)+rng.normal(0,2,muA.shape);lobs=rng.poisson(muL)+rng.normal(0,2,muL.shape)
 denom=d-40;l0=lobs[:,:,[0]]-40;floorD=int((denom<1).sum());floorL=int((l0<1).sum());ratio=(a-40)/np.maximum(denom,1);act=(lobs-40)/np.maximum(l0,1)
 x=np.stack([series[int(t)][:,:,2] for t in TX],axis=2);x+=rng.normal(0,.005,x.shape)
 return {'FRET':ratio[:,:,1:]-ratio[:,:,[0]],'actin':act[:,:,1:]-1,'x':x[:,:,1:]-x[:,:,[0]]}, {'D':floorD,'L0':floorL}

def feature(d,kind):
 return np.concatenate([d[k] for k in (('x','FRET','actin') if kind=='joint' else (kind,))],axis=2)

def fit(X):
 x=X.reshape(-1,X.shape[-1]);lab=np.tile([0,1],X.shape[0]);mean=x.mean(0);sd=np.maximum(x.std(0),1e-12);z=(x-mean)/sd;cen=np.stack([z[lab==i].mean(0) for i in range(2)]);e=z-cen[lab];inv=np.linalg.inv(e.T@e/len(lab)+.1*np.eye(X.shape[-1]));return mean,sd,cen,inv

def predict(M,X):
 mean,sd,cen,inv=M;d=(X-mean)/sd;d=d[:,:,None,:]-cen[None,None,:,:];return np.einsum('...i,ij,...j->...',d,inv,d).argmin(-1)

def wilson(k,n):
 z=1.959963984540054;p=k/n;D=1+z*z/n;m=(p+z*z/(2*n))/D;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/D;return [m-h,m+h]

def rep(pr):
 ok=pr==np.array([0,1]);k=ok.sum(0);return {'accuracy':float(ok.mean()),'correct':[int(v) for v in k],'class_wilson95':[wilson(int(v),len(pr)) for v in k],'predicted_fast':[int((pr[:,i]==0).sum()) for i in range(2)]}

def twin():
 A,B,f=simulate_pair();fp=B[:,0];err={'FRET_preoptical':float(np.max(np.abs(f-fp))),'actin':float(np.max(np.abs(A[:,1]-B[:,1]))),'deformation':float(np.max(np.abs(A[:,2]-B[:,2]))),'latent_gate_max_difference':float(np.max(np.abs(A[:,0]-B[:,0])))}
 if max(err[k] for k in ('FRET_preoptical','actin','deformation'))>1e-12 or err['latent_gate_max_difference']<.01:raise ValueError('twin check failed')
 return err,A,B,f

def main():
    out = output_directory(Path(__file__)) / "results_v1"
    out.mkdir()
    NTRAIN = 300
    NTEST = 200
    rngsT,keysT=streams(202611180101);rngsE,keysE=streams(202611180102);ptr=draw(rngsT[0],NTRAIN);pte=draw(rngsE[0],NTEST)
    # Independent implementation of a nominal step check with matched family parameters.
    nom={k:v[:1].copy() for k,v in ptr.items()};_,h1=sim(nom,'both',.25,False,True);_,h2=sim(nom,'both',.125,False,True);err=float(np.max(np.abs(h1-h2[::2])))
    if err>1e-4:raise ValueError('step check failed')
    train={};test={};floors={'train':{},'test':{}}
    for i,arm in enumerate(['u','v','both','sham'],1):
     tr,_=sim(ptr,arm);te,_=sim(pte,arm);train[arm],floors['train'][arm]=observations(ptr,tr,rngsT[i]);test[arm],floors['test'][arm]=observations(pte,te,rngsE[i])
    models={};result={};preds={}
    for arm in ('u','v','both'):
     result[arm]={};models[arm]={};preds[arm]={}
     for kind in ('x','FRET','actin','joint'):
      model=fit(feature(train[arm],kind));pr=predict(model,feature(test[arm],kind));models[arm][kind]=model;preds[arm][kind]=pr;result[arm][kind]=rep(pr)
    cx=preds['both']['x']==np.array([0,1]);cj=preds['both']['joint']==np.array([0,1]);idx=rngsE[6].integers(0,NTEST,(5000,NTEST));gain=cj.mean()-cx.mean();gs=cj[idx].mean((1,2))-cx[idx].mean((1,2))
    # Outcome of a no-input, held-out sham detector, not class-order false positives.
    null_train=train['sham']['FRET'][:,0,:];mean=null_train.mean(0);sd=np.maximum(null_train.std(0,ddof=1),1e-9)
    score=lambda f:np.max(np.abs((f-mean)/sd),axis=-1)
    threshold=float(np.quantile(score(null_train),.95));false=score(test['sham']['FRET'][:,0,:])>threshold;detect=score(test['both']['FRET'])>threshold
    # Equal-time cases do not have correct fast/slow labels.
    eq,_=sim(pte,'both',equal=True);eqobs,eqfloor=observations(pte,eq,rngsE[5]);eqpred=predict(models['both']['joint'],feature(eqobs,'joint'))
    twinerr,A,B,f=twin()
    result_out={'kind':'synthetic_only','training_parents':NTRAIN,'heldout_parents':NTEST,'seeds':{'train_root':202611180101,'train_spawn_keys':keysT,'test_root':202611180102,'test_spawn_keys':keysE},'dt_reference_max_state_difference':err,'classification':result,'primary_gain_joint_minus_x':float(gain),'primary_gain_bootstrap95':[float(v) for v in np.quantile(gs,[.025,.975])], 'sham_detection':{'training_threshold':threshold,'false_positive_count':int(false.sum()),'false_positive_n':NTEST,'false_positive_wilson95':wilson(int(false.sum()),NTEST),'both_detected_by_family':[int(detect[:,i].sum()) for i in range(2)]},'equal_time_null':{'n_arbitrary_repeats':NTEST,'called_fast_by_repeat':[int((eqpred[:,i]==0).sum()) for i in range(2)],'interpretation':'no true fast/slow labels'},'denominator_floors':floors,'equal_time_denominator_floors':eqfloor,'adversarial_twin':twinerr,'limitations':'Uncalibrated fixed input schedules and restricted optical/nuisance ranges. No biological data, fitted protein or force.'}
    (out/'results.json').write_text(json.dumps(result_out,indent=2)+'\n')
    np.savez_compressed(out/'heldout_predictions.npz',**{arm+'_'+k:v for arm in preds for k,v in preds[arm].items()},sham_false=false,equal_time_predictions=eqpred)
    fig,ax=plt.subplots(1,3,figsize=(13,4));x=np.arange(3);w=.18;colors=['#617387','#9a6b43','#648a74','#1f6c5b']
    for j,(kind,label) in enumerate([('x','Deformation'),('FRET','FRET'),('actin','Actin'),('joint','Joint')]):ax[0].bar(x+(j-1.5)*w,[result[a][kind]['accuracy']*100 for a in ('u','v','both')],w,label=label,color=colors[j])
    ax[0].set_xticks(x,['u only','v only','u + v']);ax[0].set_ylim(0,107);ax[0].axhline(50,color='black',ls=':',lw=.8);ax[0].set_ylabel('Held-out accuracy (%)');ax[0].set_title('Restricted model families')
    t=np.arange(301);ax[1].plot(t,A[:,0],label='q: model A',color='#245c91');ax[1].plot(t,B[:,0],label="q′: model B",color='#bf6d20',ls='--');ax[1].set_title('Adversarial latent states');ax[1].set_xlabel('Assigned time (s)');ax[1].set_ylabel('Deviation from baseline');ax[1].legend(frameon=False,fontsize=8)
    ax[2].plot(t,A[:,1],label='actin A / B',color='#648a74');ax[2].plot(t,A[:,2],label='deformation A / B',color='#617387');ax[2].plot(t,f,label='FRET pre-optical A / B',color='#1f6c5b',ls='--');ax[2].set_title('Identical predicted readouts');ax[2].set_xlabel('Assigned time (s)');ax[2].legend(frameon=False,fontsize=8)
    for z in ax[1:]:z.axvspan(20,80,color='#f8cf88',alpha=.3);z.axvspan(140,200,color='#9bc4ec',alpha=.25)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',ncol=4,frameon=False);fig.suptitle('Joint observation and an exact out-of-family ambiguity',fontsize=13);fig.tight_layout(rect=(0,.10,1,.93));fig.savefig(out/'joint_observation.png',dpi=180);fig.savefig(out/'joint_observation.svg');plt.close(fig)
    print(json.dumps(result_out,indent=2))


if __name__ == "__main__":
    main()
