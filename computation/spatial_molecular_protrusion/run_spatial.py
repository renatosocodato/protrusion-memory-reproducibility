#!/usr/bin/env python3
"""Frozen spatial perturbation test: see PROTOCOL_FROZEN.md."""
from pathlib import Path
from reproduction_io import output_directory
import json, math, time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=output_directory(Path(__file__))
OUT=HERE/'results_v1'
L=32.;T=60.;NTEST=200
MODES=np.array([2,3,4]); SAMPLE=(0.,15.,30.,60.)


def h(c): return c*c/(1+c*c)


def draw(rng,n):
    z=np.column_stack([rng.uniform(.16,.24,n),rng.uniform(.9,1.1,n),rng.uniform(.1,.2,n),rng.uniform(.1,.2,n),rng.uniform(.055,.075,n),rng.uniform(10,20,n),rng.uniform(15,25,n),rng.uniform(.025,.035,n)])
    phases=rng.uniform(0,2*np.pi,(n,3))
    return z,phases


def initial(n,pars,phases,amplitude_override=None):
    x=np.arange(n)*L/n
    waves=np.cos(2*np.pi*MODES[None,:,None]*x[None,None,:]/L+phases[:,:,None]).sum(axis=1)
    waves-=waves.mean(axis=1,keepdims=True)
    sig=waves.std(axis=1,keepdims=True)
    amp=pars[:,[7]] if amplitude_override is None else np.full((len(pars),1),amplitude_override)
    return waves/sig*amp


def energy(phi,kappa,dx):
    grad=(np.roll(phi,-1,axis=1)-phi)/dx
    return (0.25*(phi*phi-1)**2+0.5*kappa[:,None]*grad*grad).sum(axis=1)*dx


def simulate(kind,pars,phases,m,n=64,dt=.25,amplitude_override=None,all_checks=True):
    batch=len(pars); dx=L/n
    k=np.arange(n//2+1)
    eig=4*np.sin(np.pi*k/n)**2/dx**2
    M,kap,D,Da,lam,tb,tx,_=pars.T
    ph=initial(n,pars,phases,amplitude_override)
    baseline=.7*h(m); b=np.full((batch,n),baseline); pos=np.full(batch,baseline)
    bdecay=np.exp(-dt/tb); xdecay=np.exp(-dt/tx)
    sample=[ph.copy()];fs=[b.mean(axis=1)];xs=[pos.copy()]
    out_times=[0.]; maxmass=float(np.max(np.abs(ph.mean(axis=1)))); minconc=float(np.min(m*(1+.6*ph))); max_energy_inc=-np.inf
    Eprev=energy(ph,kap,dx) if kind=='D' else None
    nums=round(T/dt)
    for j in range(nums):
        z=np.fft.rfft(ph,axis=1)
        if kind=='D':
            cubic=np.fft.rfft(ph**3,axis=1)
            z=(z-dt*M[:,None]*eig*cubic)/(1+dt*M[:,None]*(kap[:,None]*eig*eig-eig))
        elif kind=='R': z=z*np.exp(-dt*D[:,None]*eig)
        elif kind=='A':
            react=ph-ph**3
            react-=react.mean(axis=1,keepdims=True)
            z=(z+dt*lam[:,None]*np.fft.rfft(react,axis=1))/(1+dt*Da[:,None]*eig)
        else: raise ValueError(kind)
        ph=np.fft.irfft(z,n=n,axis=1)
        c=m*(1+.6*ph)
        b=b*bdecay[:,None]+(.7*h(c))*(1-bdecay[:,None])
        F=b.mean(axis=1); pos=pos*xdecay+F*(1-xdecay)
        if all_checks:
            if not (np.all(np.isfinite(ph)) and np.all(np.isfinite(b)) and np.all(np.isfinite(pos))): raise FloatingPointError('nonfinite state')
            maxmass=max(maxmass,float(np.max(np.abs(ph.mean(axis=1)))))
            minconc=min(minconc,float(np.min(c)))
            if kind=='D':
                enow=energy(ph,kap,dx); max_energy_inc=max(max_energy_inc,float(np.max(enow-Eprev))); Eprev=enow
        t=(j+1)*dt
        if any(abs(t-u)<1e-8 for u in SAMPLE[1:]):
            sample.append(ph.copy());fs.append(F.copy());xs.append(pos.copy());out_times.append(t)
    assert out_times==list(SAMPLE),out_times
    return {'phi':np.stack(sample,axis=1),'F':np.stack(fs,axis=1),'x':np.stack(xs,axis=1),'mass_err':maxmass,'min_c':minconc,'max_energy_increase':max_energy_inc if kind=='D' else None}


def wilson(k,n):
    z=1.95996398454;p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den;half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(mid-half),float(mid+half)]


def check_reference():
    pars=np.array([[.2,1.,.15,.15,.065,15.,20.,.03]])
    phases=np.array([[.2,1.1,2.2]])
    delta_time=0.;delta_space=0.;time_fx=0.;space_fx=0.;energymax=-np.inf;mass=0.;minc=10.
    for kind in ('D','R','A'):
        base=simulate(kind,pars,phases,.3)
        half=simulate(kind,pars,phases,.3,dt=.125)
        fine=simulate(kind,pars,phases,.3,n=128)
        delta_time=max(delta_time,float(np.max(np.abs(base['phi']-half['phi']))))
        time_fx=max(time_fx,float(np.max(np.abs(base['F']-half['F']))),float(np.max(np.abs(base['x']-half['x']))))
        delta_space=max(delta_space,float(np.max(np.abs(base['phi']-fine['phi'][:,:,::2]))))
        space_fx=max(space_fx,float(np.max(np.abs(base['F']-fine['F']))),float(np.max(np.abs(base['x']-fine['x']))))
        for s in (base,half,fine):
            mass=max(mass,s['mass_err']);minc=min(minc,s['min_c'])
            if s['max_energy_increase'] is not None:energymax=max(energymax,s['max_energy_increase'])
    # Separate small-amplitude linear-rate test over five seconds.
    modes={}
    eig3=4*np.sin(np.pi*3/64)**2/(L/64)**2
    ph0=initial(64,pars,phases,amplitude_override=.0001)
    for kind in ('D','R'):
        ph=ph0.copy();dt=.25
        for j in range(round(5/dt)):
            z=np.fft.rfft(ph,axis=1)
            eig=4*np.sin(np.pi*np.arange(33)/64)**2/(L/64)**2
            if kind=='D':
                z=(z-dt*pars[:,[0]]*eig*np.fft.rfft(ph**3,axis=1))/(1+dt*pars[:,[0]]*(pars[:,[1]]*eig**2-eig))
            else:
                z=z*np.exp(-dt*pars[:,[2]]*eig)
            ph=np.fft.irfft(z,n=64,axis=1)
        observed=float(np.log(np.abs(np.fft.rfft(ph,axis=1)[0,3])/np.abs(np.fft.rfft(ph0,axis=1)[0,3]))/5)
        expected=float(.2*eig3*(1-eig3) if kind=='D' else -.15*eig3)
        modes[kind]={'observed':observed,'semidiscrete_continuous_expected':expected,'error':abs(observed-expected)}
    valid=(delta_time<=.02 and delta_space<=.02 and time_fx<=.002 and space_fx<=.002 and max(v['error'] for v in modes.values())<=.002 and mass<=1e-10 and minc>0 and energymax<=1e-9)
    return {'max_reference_dt_difference_phi':delta_time,'max_reference_grid_difference_phi':delta_space,'max_reference_dt_difference_F_x':time_fx,'max_reference_grid_difference_F_x':space_fx,'max_mass_error':mass,'minimum_concentration':minc,'max_energy_increment':energymax,'growth_rate_check':modes,'passed':bool(valid)}


def main():
    if OUT.exists(): raise SystemExit('Refusing to overwrite existing frozen output')
    OUT.mkdir();started=time.monotonic()
    checks=check_reference()
    results={'reference_checks':checks,'design':{'parents_per_regime':NTEST,'seed':27112026,'paired_models':True,'threshold_G60':1.3,'grid':64,'domain':32,'dt':.25,'T':60,'noise_sd_phi':.002,'synthetic_only':True},'regimes':{}}
    (OUT/'checks.json').write_text(json.dumps(checks,indent=2))
    if not checks['passed']: raise SystemExit('Reference acceptance failed: preserve checks and stop before reserved test')
    rng=np.random.default_rng(27112026); rows=[]; exemplars={}
    for m in (.3,1.):
        pars,phases=draw(rng,NTEST)
        per={}; perrep={}
        for kind in ('D','R','A'):
            s=simulate(kind,pars,phases,m)
            if s['mass_err']>1e-10 or s['min_c']<=0 or (kind=='D' and s['max_energy_increase']>1e-9): raise RuntimeError('Reserved ensemble violated preset invariants')
            measured=s['phi']+rng.normal(0,.002,s['phi'].shape)
            sigma=measured.std(axis=2)
            g=sigma[:,-1]/sigma[:,0]
            ff=np.fft.rfft(measured,axis=2)
            amps=np.abs(ff[:,:,MODES])
            rates15=np.log(amps[:,1,:]/amps[:,0,:])/15
            rates60=np.log(amps[:,3,:]/amps[:,0,:])/60
            pred=g>1.3; k=int(pred.sum()); n=len(pred)
            forcechange=s['F'][:,-1]-s['F'][:,0]; dx=s['x'][:,-1]-s['x'][:,0]
            per[kind]={'G60_median':float(np.median(g)),'G60_quantiles_025_975':np.quantile(g,[.025,.975]).tolist(),'called_demixing':k,'n':n,'called_demixing_fraction':k/n,'wilson95':wilson(k,n),'force_change_median':float(np.median(forcechange)),'force_change_quantiles_025_975':np.quantile(forcechange,[.025,.975]).tolist(),'position_change_median':float(np.median(dx)),'position_change_quantiles_025_975':np.quantile(dx,[.025,.975]).tolist(),'mass_error':s['mass_err'],'minimum_c':s['min_c'],'max_energy_increase':s['max_energy_increase']}
            perrep[kind]=s
            for i in range(NTEST):
                rows.append({'m':m,'parent':i,'model':kind,'G60':g[i],'prediction_D':bool(pred[i]),'force_change':forcechange[i],'position_change':dx[i],'g2_15':rates15[i,0],'g3_15':rates15[i,1],'g4_15':rates15[i,2],'g2_60':rates60[i,0],'g3_60':rates60[i,1],'g4_60':rates60[i,2]})
        # Deterministic identical initial condition and same-data zero separation.
        assert all(np.array_equal(perrep['D']['phi'][:,0],perrep[k]['phi'][:,0]) for k in ('R','A'))
        assert all(np.array_equal(perrep['D']['F'][:,0],perrep[k]['F'][:,0]) for k in ('R','A'))
        exemplars[m]=perrep;results['regimes'][str(m)]=per
    pd.DataFrame(rows).to_csv(OUT/'heldout_parents.csv',index=False)
    results['elapsed_seconds_before_figures']=round(time.monotonic()-started,3)
    (OUT/'results.json').write_text(json.dumps(results,indent=2))
    fig,axes=plt.subplots(2,3,figsize=(12,6.2),layout='constrained'); colors={'D':'#a0476c','R':'#337596','A':'#c28224'}
    for j,(m,per) in enumerate(exemplars.items()):
        coords=np.arange(64)*L/64
        for kind in ('D','R','A'):
            s=per[kind]; label={'D':'Demixing model','R':'Passive recruitment','A':'Active recruitment'}[kind]
            axes[j,0].plot(coords,s['phi'][0,3],label=label,color=colors[kind])
            axes[j,1].plot(SAMPLE,s['F'][0]-s['F'][0,0],color=colors[kind],marker='o',ms=3)
            axes[j,2].plot(SAMPLE,s['x'][0]-s['x'][0,0],color=colors[kind],marker='o',ms=3)
        axes[j,0].set_ylabel(f'm={m}\ncontrast φ')
        axes[j,1].set_ylabel('Change in actin/load (a.u.)')
        axes[j,2].set_ylabel('Deformation change (a.u.)')
        axes[j,0].set_xlabel('Position (dimensionless)')
        for a in axes[j,1:]:a.set_xlabel('Assigned time (s)')
        for a in axes[j]:a.spines[['top','right']].set_visible(False)
    axes[0,0].legend(frameon=False,fontsize=8);axes[0,0].set_title('Spatial profile at 60 s');axes[0,1].set_title('Actin / effective force');axes[0,2].set_title('Protrusion displacement')
    fig.suptitle('Constructed spatial models; one matched synthetic parent',fontsize=12)
    fig.savefig(OUT/'spatial_representative.png',dpi=190);fig.savefig(OUT/'spatial_representative.svg');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,4.2),layout='constrained')
    for j,m in enumerate((.3,1.)):
        vals=[results['regimes'][str(m)][kind]['called_demixing_fraction'] for kind in ('D','R','A')]
        axes[j].bar(['Demixing','Passive\nrecruitment','Active\nrecruitment'],vals,color=[colors[k] for k in ('D','R','A')])
        for i,v in enumerate(vals): axes[j].text(i,min(v+.025,1.02),f'{v:.1%}',ha='center',fontsize=10)
        axes[j].set_ylim(0,1.15); axes[j].set_ylabel('Fraction labeled demixing');axes[j].set_title(f'm={m}; 200 parent cases');axes[j].spines[['top','right']].set_visible(False)
    fig.suptitle('Predeclared 60-s spatial-contrast rule',fontsize=12)
    fig.savefig(OUT/'spatial_challenge.png',dpi=190);fig.savefig(OUT/'spatial_challenge.svg');plt.close(fig)
    print(json.dumps(results,indent=2))

if __name__=='__main__':main()
