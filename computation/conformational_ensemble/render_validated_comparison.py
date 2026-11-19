from pathlib import Path
from reproduction_io import output_directory
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
r=output_directory(Path(__file__))
lo=json.loads((r/'results_v1/results.json').read_text())['m']['0.3'];hi=json.loads((r/'results_v2_fresh/results.json').read_text());data=[lo,hi]
fig,axes=plt.subplots(1,2,figsize=(10.5,4.1),gridspec_kw={'width_ratios':[1.2,1]})
x=np.arange(2);width=.32
for label,key,offset,color in [('Position','position_only',-width/2,'#6b7b8b'),('Position + gate','position_plus_gate',width/2,'#237a64')]:
 vals=[100*d[key]['accuracy'] for d in data];bars=axes[0].bar(x+offset,vals,width,label=label,color=color)
 for b,v in zip(bars,vals):axes[0].text(b.get_x()+b.get_width()/2,v+1,f'{v:g}%',ha='center',fontsize=9)
axes[0].axhline(50,color='black',ls=':',lw=1);axes[0].set_ylim(0,112);axes[0].set_ylabel('Held-out accuracy (%)');axes[0].set_xticks(x,['m = 0.3','m = 1.0']);fig.legend(*axes[0].get_legend_handles_labels(),frameon=False,loc='lower center',ncol=2)
for i,d in enumerate(data):
 est=100*d['paired_gain'];lo,hi=(100*v for v in d['paired_gain_bootstrap95']);axes[1].errorbar(est,i,xerr=[[est-lo],[hi-est]],fmt='o',capsize=4,color='#237a64');axes[1].text(hi+.4,i,f'{est:g} pp',va='center',fontsize=9)
axes[1].axvline(0,color='black',lw=.8);axes[1].set_yticks(x,['m = 0.3','m = 1.0']);axes[1].invert_yaxis();axes[1].set_xlim(-1,24);axes[1].set_xlabel('Paired accuracy gain (percentage points)');axes[1].set_title('Parent bootstrap 95% interval')
fig.suptitle('Synthetic conformational-gate challenge',fontsize=13);fig.tight_layout(rect=(0,.11,1,.94));out=r/'results_v2_fresh';fig.savefig(out/'ensemble_gate.png',dpi=180);fig.savefig(out/'ensemble_gate.svg');plt.close(fig)
