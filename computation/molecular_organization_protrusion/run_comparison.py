#!/usr/bin/env python3
"""Frozen synthetic molecular organization -> actin -> deformation comparison.
See PROTOCOL_FROZEN.md. No empirical observations are used by this simulator.
"""
from pathlib import Path
from reproduction_io import output_directory
import json
import math
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = output_directory(Path(__file__))
OUT = ROOT / 'results_v1'
CLASSES = ('uniform', 'finite_recruitment', 'demixing_analogue')
MVALUES = (0.3, 1.0)
TIMES = np.arange(0, 301, 60)
N = 300


def parameters(rng, n):
    ranges = ((6,12), (5,12), (8,20), (10,25), (12,35), (.45,.75))
    return np.column_stack([rng.uniform(a,b,n) for a,b in ranges])


def signal(c):
    return c*c / (1.0+c*c)


def rates(y, t, par, m, kind, null):
    p, a, r, b, x = y.T
    tp, ta, tr, tb, tx, R = par.T
    u = float(20 <= t < 200)
    dp = (u-p)/tp
    da = (.5+.4*p-a)/ta
    mob = 1-r*r
    if kind == 'uniform':
        dr = np.zeros_like(r)
    elif kind == 'finite_recruitment':
        dr = -mob*(r-R*p)/tr
    elif kind == 'demixing_analogue':
        dr = -mob*r*(r*r-R*R*p)/(tr*R*R)
    elif kind == 'driven_recruitment':
        dr = -mob*(r-R*p*(.65+.35*np.sin(2*np.pi*t/75)))/tr
    else:
        raise ValueError(kind)
    H = signal(m) if null else .5*(signal(m*(1-r))+signal(m*(1+r)))
    db = (a*H-b)/tb
    dx = (b-x)/tx
    return np.column_stack((dp, da, dr, db, dx))


def simulate(par, m, kind, dt=.5, null=False, full=False):
    n=len(par)
    p=np.zeros(n); a=np.full(n,.5)
    r=np.zeros(n) if kind=='uniform' else np.full(n,.03)
    H=np.full(n,signal(m)) if null else .5*(signal(m*(1-r))+signal(m*(1+r)))
    b=a*H; x=b.copy()
    y=np.column_stack((p,a,r,b,x))
    steps=round(300/dt); stride=round(60/dt)
    sampled=[y.copy()]; hist=[y.copy()] if full else None
    for j in range(steps):
        t=j*dt
        k1=rates(y,t,par,m,kind,null)
        y=y+dt*rates(y+.5*dt*k1,t+.5*dt,par,m,kind,null)
        if not np.all(np.isfinite(y)) or np.any(y[:,2]<-1e-10) or np.any(y[:,2]>1+1e-10) or np.any(y[:,0:2]<-1e-10) or np.any(y[:,3:]<-1e-10):
            raise ValueError(f'physical bounds or finite check at {t}, {kind}')
        if (j+1)%stride==0: sampled.append(y.copy())
        if full: hist.append(y.copy())
    result=np.stack(sampled,axis=1)
    return (result,np.stack(hist,axis=1)) if full else result


def observe(trajectory, rng):
    n=trajectory.shape[0]
    noisy_x=trajectory[:,:,4]+rng.normal(0,.005,(n,6))
    noisy_r=trajectory[:,:,2]+rng.normal(0,.03,(n,6))
    return noisy_x[:,1:]-noisy_x[:,[0]], noisy_r


def fit(X, labels):
    mean=X.mean(axis=0); scale=X.std(axis=0); scale=np.maximum(scale,1e-12)
    z=(X-mean)/scale
    means=np.array([z[labels==k].mean(axis=0) for k in range(3)])
    res=z-means[labels]
    covariance=res.T@res/len(labels)+.01*np.eye(X.shape[1])
    return mean,scale,means,np.linalg.inv(covariance)


def predict(model,X):
    mean,scale,centers,inv=model
    v=(X-mean)/scale
    delta=v[:,None,:]-centers[None,:,:]
    dist=np.einsum('nki,ij,nkj->nk',delta,inv,delta)
    return np.argmin(dist,axis=1)


def wilson(k,n):
    z=1.95996398454
    p=k/n; d=1+z*z/n
    center=(p+z*z/(2*n))/d
    half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return [center-half,center+half]


def counts_and_stats(y, yhat):
    cm=np.array([[int(np.sum((y==i)&(yhat==j))) for j in range(3)] for i in range(3)])
    sens=[{'success':int(cm[i,i]),'n':int(cm[i].sum()),'estimate':float(cm[i,i]/cm[i].sum()),'wilson95':wilson(int(cm[i,i]),int(cm[i].sum()))} for i in range(3)]
    fp=int(cm[0,2]+cm[1,2]); fn=int(cm[0].sum()+cm[1].sum())
    return {'confusion':cm.tolist(),'sensitivity':sens,'overall_accuracy':float(np.trace(cm)/cm.sum()),'phase_label_fpr':{'success':fp,'n':fn,'estimate':fp/fn,'wilson95':wilson(fp,fn)}}


def dataset(m,seed,n=300,variant='ordinary'):
    rng=np.random.default_rng(seed)
    par=parameters(rng,n)
    xsets=[]; rsets=[]; tr=[]
    kinds=CLASSES if variant in ('ordinary','null') else ('driven_recruitment',)
    for kind in kinds:
        raw=simulate(par,m,kind,null=(variant=='null'))
        xx,rr=observe(raw,rng)
        tr.append(raw);xsets.append(xx);rsets.append(rr)
    labels=np.repeat(np.arange(len(kinds)),n)
    return np.concatenate(xsets),np.concatenate(rsets),labels,tr


def main():
    if OUT.exists(): raise SystemExit('Output directory already exists; preserve frozen result.')
    OUT.mkdir()
    start=time.monotonic()
    results={'design':{'n_per_class_per_regime':N,'training_seed':7717,'test_seed':9917,'out_of_family_seed':11917,'null_seed':13917,'noise_x':.005,'noise_r':.03,'solver':'RK2 dt 0.5 s','estimated_position_units':'dimensionless','biological_calibration':False},'regimes':{},'checks':{}}
    # Algebraic check independent of implementation of the signal function.
    for m in MVALUES:
        d=.8*m
        observed=.5*(signal(m-d)+signal(m+d))-signal(m)
        expected=d*d*(1+d*d-3*m*m)/((1+m*m)*(1+(m-d)**2)*(1+(m+d)**2))
        assert abs(observed-expected)<1e-14
    # Same parameters and initial conditions, exact invariant of null model.
    par=np.array([[9.,8.,14.,18.,24.,.6]])
    ref_data={}; maxdiff=0.; min_diss=0.
    for m in MVALUES:
        for name in CLASSES:
            slow=simulate(par,m,name,dt=.5)
            fine=simulate(par,m,name,dt=.25)
            err=float(np.max(np.abs(slow[:,:,2:]-fine[:,:,2:])))
            maxdiff=max(maxdiff,err)
            ref_data[(m,name)]=simulate(par,m,name,dt=.5,full=True)[1][0]
        nulls=[simulate(par,m,name,null=True)[:,:,3:] for name in CLASSES]
        assert all(np.array_equal(nulls[0],v) for v in nulls[1:])
        # Fixed signal landscapes: gradient times dynamics <=0, including entire r grid.
        rr=np.linspace(0,1,101); pp=.7; R=.6
        for g in (rr-R*pp, rr*(rr*rr-R*R*pp)/(R*R)):
            diss=-(1-rr*rr)*g*g/14
            assert np.max(diss)<=1e-14
            min_diss=min(min_diss,float(np.min(diss)))
    results['checks']={'analytic_identity':'passed','null_position_identical':'passed','fixed_p_gradient_dissipation':'passed','maximum_abs_dt_half_at_observation_times':maxdiff,'tolerance':.001,'timestep_check_passed':maxdiff<=.001,'bounds':'checked at every integration step'}
    if maxdiff>.001:
        (OUT/'results.json').write_text(json.dumps(results,indent=2)); raise SystemExit('Frozen numerical acceptance criterion failed.')
    data=[]; predstore={}
    for idx,m in enumerate(MVALUES):
        trainx,trainr,yt,_=dataset(m,7717+idx)
        testx,testr,y,raw=dataset(m,9917+idx)
        outx,outr,_,_=dataset(m,11917+idx,variant='out')
        nullx,nullr,yn,_=dataset(m,13917+idx,variant='null')
        item={}
        for mode in ('position','position_plus_organization'):
            X=trainx if mode=='position' else np.column_stack((trainx,trainr))
            Z=testx if mode=='position' else np.column_stack((testx,testr))
            O=outx if mode=='position' else np.column_stack((outx,outr))
            model=fit(X,yt); pred=predict(model,Z); outpred=predict(model,O)
            stats=counts_and_stats(y,pred)
            stats['out_of_family_prediction_fractions']={k:float(np.mean(outpred==i)) for i,k in enumerate(CLASSES)}
            if mode=='position':
                npred=predict(model,nullx)
                stats['null_confusion']=counts_and_stats(yn,npred)['confusion']
                stats['null_C_label_fraction']=float(np.mean(npred==2))
            item[mode]=stats
            predstore[(m,mode)]=pred
            for ii in range(len(y)):
                data.append({'mean_concentration':m,'channel':mode,'parent_id':ii%N,'true_class':CLASSES[y[ii]],'predicted_class':CLASSES[pred[ii]]})
        results['regimes'][str(m)]=item
    pd.DataFrame(data).to_csv(OUT/'heldout_predictions.csv',index=False)
    results['elapsed_seconds']=round(time.monotonic()-start,3)
    (OUT/'results.json').write_text(json.dumps(results,indent=2))
    # Two simple scientific figures, no figure-numbered manuscript export.
    fig,ax=plt.subplots(2,3,figsize=(12,6.5),layout='constrained')
    tt=np.arange(0,300.01,.5)
    colors=('#35608c','#cd8529','#9b4775')
    for row,m in enumerate(MVALUES):
        for i,name in enumerate(CLASSES):
            a=ref_data[(m,name)]
            ax[row,0].plot(tt,a[:,2],color=colors[i],label=name.replace('_',' '))
            ax[row,1].plot(tt,a[:,3],color=colors[i])
            ax[row,2].plot(tt,a[:,4]-a[0,4],color=colors[i])
            ax[row,2].scatter(TIMES,(a[::120,4]-a[0,4]),color=colors[i],s=12)
        ax[row,0].set_ylabel(f'm = {m}\nEffective contrast')
        ax[row,1].set_ylabel('Actin output (a.u.)')
        ax[row,2].set_ylabel('Displacement change (a.u.)')
        for a in ax[row]:
            a.axvspan(20,200,color='#888888',alpha=.08)
            a.set_xlabel('Time (s)')
            a.spines[['top','right']].set_visible(False)
    ax[0,0].legend(frameon=False,fontsize=8,loc='upper left')
    for a,title in zip(ax[0],('Molecular organization','Actin output','Position and 60-s samples')):a.set_title(title)
    fig.suptitle('Constructed two-domain systems; parameters assumed',fontsize=12)
    fig.savefig(OUT/'trajectories.png',dpi=190);fig.savefig(OUT/'trajectories.svg');plt.close(fig)
    fig,ax=plt.subplots(2,2,figsize=(8,7),layout='constrained')
    for row,m in enumerate(MVALUES):
        for col,mode in enumerate(('position','position_plus_organization')):
            cm=np.array(results['regimes'][str(m)][mode]['confusion'])/N
            ax[row,col].imshow(cm, vmin=0,vmax=1,cmap='Blues')
            for i in range(3):
                for j in range(3):ax[row,col].text(j,i,f'{cm[i,j]:.2f}',ha='center',va='center',color='white' if cm[i,j]>.55 else 'black')
            ax[row,col].set_xticks(range(3),['uniform','finite','demixing'],rotation=30,ha='right'); ax[row,col].set_yticks(range(3),['uniform','finite','demixing'])
            ax[row,col].set_xlabel('Predicted synthetic class');ax[row,col].set_ylabel('Generating class')
            ax[row,col].set_title(f'm = {m}; '+('position' if col==0 else 'position + contrast'))
    fig.suptitle('Independent synthetic challenge; fractions per class',fontsize=12)
    fig.savefig(OUT/'heldout_confusion.png',dpi=190);fig.savefig(OUT/'heldout_confusion.svg');plt.close(fig)
    print(json.dumps(results,indent=2))

if __name__=='__main__':
    main()
