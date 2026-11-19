#!/usr/bin/env python3
"""Finite illustrative dynamics after stationary sensitivity; no cell fitting."""
import json
from pathlib import Path
from reproduction_io import output_directory
import numpy as np

HERE = output_directory(Path(__file__))
M = np.array([0.73, 1.0, 1.46])
SEED = 202611180401


def rhs(y, drive):
    q, a, x, s = y.T
    bias = .2 + drive - .4*x
    on, off = .5*np.exp(bias/2), .5*np.exp(-bias/2)
    z = M*q
    return np.stack((on*(1-q)-off*q,
                     z/(1+z)*(1-a)*(1+.5*a)-.4*a,
                     (a-x)/2, .26-.08*x-s), axis=1)


def run(arm, dt):
    rng = np.random.default_rng(SEED)
    y = np.tile([.3, .1, .1, .25], (len(M), 1))
    ou = 0.
    samples = [y.copy()]
    for i in range(round(400/dt)):
        tmid = (i+.5)*dt
        u = 0. if arm == 'autonomous' else (.4*np.sin(2*np.pi*tmid/30) if arm == 'periodic' else ou)
        k1 = rhs(y, u)
        k2 = rhs(y+dt*k1/2, u)
        k3 = rhs(y+dt*k2/2, u)
        k4 = rhs(y+dt*k3, u)
        y = y+dt*(k1+2*k2+2*k3+k4)/6
        if not np.isfinite(y).all() or np.any((y[:, :3] < 0) | (y[:, :3] > 1+1e-12)):
            raise ValueError('Nonfinite or out-of-domain state')
        if arm == 'stochastic':
            ou = np.exp(-dt/4)*ou + .4*np.sqrt(1-np.exp(-2*dt/4))*rng.normal()
        if (i+1) % round(2/dt) == 0:
            samples.append(y.copy())
    return np.asarray(samples)


def main():
    arrays = {}
    cases = {}
    for arm in ('autonomous', 'periodic', 'stochastic'):
        arr = run(arm, .02)
        arrays[arm] = arr
        post = arr[100:]
        mu = post.mean(axis=0)
        sd = post.std(axis=0, ddof=1)
        cases[arm] = [{'m': float(m), 'q_mean': float(mu[i, 0]), 'q_sd': float(sd[i, 0]),
                       'load_mean': float(mu[i, 2]), 'load_sd': float(sd[i, 2]),
                       'delta_load_mean_from_m_1': float(mu[i, 2]-mu[1, 2])}
                      for i, m in enumerate(M)]
    err = float(np.max(np.abs(run('periodic', .01)-arrays['periodic'])))
    assert err < 1e-4
    result = {'purpose': 'Prior-informed scenarios with assumed RNA-to-protein transfer, not biological calibration',
              'times': '0..400 with samples every 2; summaries use 200..400 inclusive',
              'seed_ou': SEED, 'periodic_step_halving_max_error': err, 'cases': cases}
    np.savez_compressed(HERE/'dynamic_trajectories.npz', **arrays, t=np.arange(201)*2, m=M)
    (HERE/'dynamic_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
