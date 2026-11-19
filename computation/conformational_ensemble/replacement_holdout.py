#!/usr/bin/env python3
"""Reconstruct frozen m=1 training and apply it to new, independent parents."""
from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np

from run_challenge import draw, features, fit, predict, report, simulate

root = output_directory(Path(__file__))
n=200;m=1.
train=draw(2026111811,n);test=draw(2026111842,n)
trf,_=simulate(train,m,train[:,5]);trs,_=simulate(train,m,train[:,6]);tef,_=simulate(test,m,test[:,5]);tes,_=simulate(test,m,test[:,6])
xt,qt=features(trf,trs,np.random.default_rng(2026111812));xe,qe=features(tef,tes,np.random.default_rng(2026111852))
mx=fit(xt);mq=fit(np.concatenate((xt,qt),axis=2));px=predict(mx,xe);pq=predict(mq,np.concatenate((xe,qe),axis=2));cx=px==np.array([0,1]);cq=pq==np.array([0,1])
idx=np.random.default_rng(2026111872).integers(0,n,(5000,n));g=cq[idx].mean((1,2))-cx[idx].mean((1,2))
na,_=simulate(test,m,test[:,5]);nb,_=simulate(test,m,test[:,5]);nx,nq=features(na,nb,np.random.default_rng(2026111862));nullx=predict(mx,nx);nullq=predict(mq,np.concatenate((nx,nq),axis=2))
r={'m':m,'source':'new reserved synthetic parent set','seed_roles':{'training_parameters':2026111811,'training_noise':2026111812,'test_parameters':2026111842,'test_noise':2026111852,'null_noise':2026111862,'bootstrap':2026111872},'position_only':report(px),'position_plus_gate':report(pq),'paired_gain':float(cq.mean()-cx.mean()),'paired_gain_bootstrap95':[float(x) for x in np.quantile(g,[.025,.975])],'matched_null_position':report(nullx),'matched_null_position_plus_gate':report(nullq),'gate_q35_median_fast':float(np.median(tef[35][:,2])),'gate_q35_median_slow':float(np.median(tes[35][:,2])),'x60_median_fast':float(np.median(tef[60][:,4])),'x60_median_slow':float(np.median(tes[60][:,4])),'supersedes':'m=1.0 results_v1 held-out independence claim; earlier results retained'}
out=root/'results_v2_fresh';out.mkdir(exist_ok=False);(out/'results.json').write_text(json.dumps(r,indent=2)+'\n');np.savetxt(out/'heldout_predictions.csv',np.column_stack((np.arange(n),px,pq)),delimiter=',',header='parent,pos_fast,pos_slow,posq_fast,posq_slow',comments='',fmt='%d');print(json.dumps(r,indent=2))
