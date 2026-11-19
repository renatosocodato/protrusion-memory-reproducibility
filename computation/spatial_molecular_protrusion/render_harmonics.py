from pathlib import Path
from reproduction_io import output_directory
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=output_directory(Path(__file__))/'results_v1';d=pd.read_csv(p/'fourth_challenge_parents.csv')
fig,axes=plt.subplots(2,2,figsize=(9,7),layout='constrained')
colors=['#a0476c','#696493']; labs=['Demixing','Active lookalike']
for j,m in enumerate((.3,1.)):
 a=d[d.m==m]
 for col,key in enumerate(('Q','force_change')):
  values=[a[a.model==k][key].to_numpy() for k in ('D','H')]
  bp=axes[j,col].boxplot(values,patch_artist=True,showfliers=True,flierprops={'marker':'.','markersize':3,'alpha':0.7},whis=(2.5,97.5),tick_labels=labs)
  for patch,c in zip(bp['boxes'],colors):patch.set_facecolor(c);patch.set_alpha(.75)
  axes[j,col].spines[['top','right']].set_visible(False)
  axes[j,col].set_title(f'm = {m}')
  axes[j,col].set_ylabel('Nonlinear ratio Q' if col==0 else 'Effective force change (dimensionless)')
  if col==0:axes[j,col].axhline(1.2,color='#555555',linestyle='--',lw=1,label='Frozen threshold');axes[j,col].legend(frameon=False,fontsize=8)
  else:axes[j,col].axhline(0,color='#777777',lw=.8)
fig.suptitle('Independent amplitude challenge: model-specific nonlinear response',fontsize=12)
fig.savefig(p/'nonlinear_harmonic_challenge.png',dpi=190);fig.savefig(p/'nonlinear_harmonic_challenge.svg');plt.close(fig)
