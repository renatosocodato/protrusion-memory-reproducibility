from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=output_directory(Path(__file__))/'results_v1'
a=pd.read_csv(p/'second_challenge_parents.csv');b=pd.read_csv(p/'third_challenge_parents.csv');b=b.rename(columns={'force_change':'delta_force_measured'});b['model']='H'
d=pd.concat([a,b],ignore_index=True)
colors=['#a0476c','#337596','#c28224','#696493'];labels=['Demixing','Passive\nrecruit.','Active\nrecruit.','Active\nlookalike'];ks=['D','R','A','H']
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
for row,m in enumerate((.3,1.)):
 sub=d[d.m==m]
 for col,key in enumerate(('S','delta_force_measured')):
  values=[sub[sub.model==k][key].to_numpy() for k in ks]
  bp=axes[row,col].boxplot(values,patch_artist=True,showfliers=False,whis=(2.5,97.5),tick_labels=labels)
  for patch,c in zip(bp['boxes'],colors):patch.set_facecolor(c);patch.set_alpha(.7)
  axes[row,col].axhline(0,color='#777777',lw=.8)
  axes[row,col].set_ylabel('S (per assigned second)' if col==0 else 'Force change (dimensionless)')
  axes[row,col].tick_params(axis='x',labelsize=8)
 ax=axes[row,2];barx=np.arange(4);w=.34
 vals1=[sub[sub.model==k].spectral_positive.mean() for k in ks]
 vals2=[sub[sub.model==k].joint_positive.mean() for k in ks]
 ax.bar(barx-w/2,vals1,w,label='Spectrum',color='#4478a3');ax.bar(barx+w/2,vals2,w,label='Spectrum + force',color='#d5983b')
 ax.set_xticks(barx,labels,fontsize=8);ax.set_ylim(0,1.12);ax.set_ylabel('Fraction labeled demixing');ax.spines[['top','right']].set_visible(False)
 for i,(v1,v2) in enumerate(zip(vals1,vals2)):
  ax.text(i-w/2,min(v1+.025,1.04),f'{v1:.2f}',ha='center',fontsize=8);ax.text(i+w/2,min(v2+.025,1.04),f'{v2:.2f}',ha='center',fontsize=8)
 for c in range(3):
  axes[row,c].spines[['top','right']].set_visible(False)
  if c==0:axes[row,c].set_title(f'm={m}: 15-s spectral response')
  if c==1:axes[row,c].set_title(f'm={m}: 60-s force response')
  if c==2:axes[row,c].set_title(f'm={m}: classification')
axes[0,2].legend(frameon=False,fontsize=8,loc='upper center')
fig.suptitle('Independent synthetic challenges; active lookalike was held out',fontsize=13)
fig.savefig(p/'spectral_mechanical_challenge.png',dpi=190);fig.savefig(p/'spectral_mechanical_challenge.svg');plt.close(fig)
