from pathlib import Path
from reproduction_io import output_directory
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=output_directory(Path(__file__));j=json.loads((p/'amiri_summary.json').read_text());d=np.load(p/'synthetic_v1/trajectories.npz');
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,3,figsize=(12,3.7),gridspec_kw={'width_ratios':[1,1,1]});styles={'sensor_calyculin':('#b9632f','TS + calyculin'),'sensor_blebbistatin':('#367b95','TS + blebbistatin'),'control_calyculin':('#737373','CTS + calyculin')}
for name,(color,label) in styles.items():
 a=np.array(j['conditions'][name]['values']);breaks=np.flatnonzero(np.diff(a[:,0])>5)+1
 for n,seg in enumerate(np.split(a,breaks)):
  axs[0].plot(seg[:,0],seg[:,1],'.-',ms=2.4,lw=1,color=color,label=label if n==0 else None)
  axs[1].plot(seg[:,0],100*seg[:,2],'.-',ms=2.4,lw=1,color=color)
axs[0].set(xlabel='Time (min)',ylabel='Mean traction magnitude (Pa)',title='Published examples: one cell each');axs[1].set(xlabel='Time (min)',ylabel='FRET efficiency (%)',title='Paired sensor signal')
t=d['t'];mask=(t>=200)&(t<=320)
for key,color,label in [('autonomous','#777777','autonomous'),('periodic','#c06b37','periodic drive'),('stochastic','#327c94','stochastic drive')]:axs[2].plot(t[mask]-200,d[key][mask,2],lw=1,color=color,label=label)
axs[2].set(xlabel='Assigned model time',ylabel='Model load, dimensionless',title='Separate synthetic trajectories');axs[2].legend(frameon=False,fontsize=8)
axs[0].legend(frameon=False,fontsize=7,loc='upper left');
for i,ax in enumerate(axs):ax.text(-.14,1.05,chr(65+i),transform=ax.transAxes,fontsize=11,fontweight='bold');ax.tick_params(length=3)
fig.tight_layout();out=p/'synthetic_v1';fig.savefig(out/'data_and_synthetic.png',dpi=200,bbox_inches='tight');fig.savefig(out/'data_and_synthetic.svg',bbox_inches='tight');plt.close(fig)
