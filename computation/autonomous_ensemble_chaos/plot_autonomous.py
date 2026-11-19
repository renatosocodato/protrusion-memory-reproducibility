from pathlib import Path
from reproduction_io import output_directory
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=output_directory(Path(__file__))/'results_v1';d=json.loads((p/'results.json').read_text());data=np.load(p/'representative_trace.npz');t=data['t'];u,v,w,x=data['states'].T;q=1/(1+np.exp(-u/20));a=1/(1+np.exp(-v/30));m=1/(1+np.exp(-(w-25)/30))
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axs=plt.subplots(1,3,figsize=(11.5,3.65),gridspec_kw={'width_ratios':[1.25,1,1]})
sel=t<=60;axs[0].plot(t[sel],q[sel],color='#385f88',lw=1,label='q: gate proxy');axs[0].plot(t[sel],a[sel],color='#a66b27',lw=.85,alpha=.8,label='a: actin proxy');axs[0].plot(t[sel],x[sel],color='#278267',lw=1.2,label='x: displacement');axs[0].set(xlabel='Assigned time',ylabel='Normalized proxy',ylim=(.2,.85),title='Autonomous output');axs[0].legend(frameon=False,fontsize=7,loc='upper right')
sel=t<=200;axs[1].plot(q[sel],m[sel],color='#326f86',alpha=.35,lw=.45);axs[1].set(xlabel='q: gate proxy',ylabel='m: feedback proxy',title='Bounded state projection')
keys=['low_control','canonical','feedback','initial_b','initial_c','half_step'];labels=['r=10','Lorenz','r=28, A','r=28, B','r=28, C','r=28, h/2'];vals=[d['cases'][z]['spectrum'][0] for z in keys];y=np.arange(len(keys));axs[2].axvline(0,color='#606970',lw=.8);axs[2].scatter(vals,y,c=['#777777','#7794ac','#bd6a38','#bd6a38','#bd6a38','#bd6a38'],s=28,zorder=3);axs[2].set(yticks=y,yticklabels=labels,xlabel='Largest finite-time exponent / model time',title='Full-state tangent estimate',xlim=(-.75,1.15));axs[2].invert_yaxis()
for i,ax in enumerate(axs):ax.text(-.15,1.05,chr(65+i),transform=ax.transAxes,fontweight='bold',fontsize=11);ax.tick_params(length=3)
fig.suptitle('Autonomous dynamics in a constructed feedback benchmark',fontsize=12,y=1.02);fig.tight_layout();fig.savefig(p/'autonomous_chaos.png',dpi=200,bbox_inches='tight');fig.savefig(p/'autonomous_chaos.svg',bbox_inches='tight');plt.close(fig)
