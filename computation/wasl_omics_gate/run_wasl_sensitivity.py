#!/usr/bin/env python3
"""Prespecified prior-predictive illustration, not cellular data fitting."""
import json
from pathlib import Path
from reproduction_io import output_directory

import numpy as np
from scipy.optimize import brentq

HERE = output_directory(Path(__file__))
RATIO = 1.46
SEED = 202611180501
N = 20000


def state(m):
    """Unique stationary root for the frozen parameters, vectorized bisection."""
    m = np.asarray(m, dtype=float)
    lo = np.zeros_like(m)
    hi = np.ones_like(m)
    for _ in range(60):
        a = (lo + hi) / 2
        q = 1 / (1 + np.exp(-(0.2 - 0.4 * a)))
        z = m * q
        f = z / (1 + z) * (1 - a) * (1 + 0.5 * a) - 0.4 * a
        lo = np.where(f > 0, a, lo)
        hi = np.where(f > 0, hi, a)
    a = (lo + hi) / 2
    q = 1 / (1 + np.exp(-(0.2 - 0.4 * a)))
    z = m * q
    f = z / (1 + z) * (1 - a) * (1 + 0.5 * a) - 0.4 * a
    return q, a, f


def main():
    baseline_q, baseline_a, baseline_resid = [float(x) for x in state(1.0)]
    scenarios = []
    for beta in (0.0, 0.5, 1.0):
        for xi in (0.5, 1.0, 2.0):
            m1 = RATIO**beta * xi
            q, a, resid = [float(x) for x in state(m1)]
            scenarios.append({'beta': beta, 'xi': xi, 'm1': m1, 'q1': q, 'a1': a,
                              'F1': a, 'delta_F': a-baseline_a, 'root_residual': resid})
    rng = np.random.default_rng(SEED)
    prior = []
    max_residual = abs(baseline_resid)
    for width in (1, 2):
        for beta_mode in ('fixed_1', 'uniform_0_1'):
            m0 = 10 ** rng.uniform(-1, 1, N)
            xi = 10 ** rng.uniform(-width, width, N)
            beta = np.ones(N) if beta_mode == 'fixed_1' else rng.uniform(0, 1, N)
            m1 = m0 * RATIO**beta * xi
            q0, a0, r0 = state(m0)
            q1, a1, r1 = state(m1)
            if not (np.all((q0 >= 0) & (q0 <= 1) & (q1 >= 0) & (q1 <= 1)) and
                    np.all((a0 >= 0) & (a0 <= 1) & (a1 >= 0) & (a1 <= 1))):
                raise ValueError('State outside physical interval')
            delta = a1-a0
            prob = float(np.mean(delta > 0))
            prior.append({'xi_log10_half_width': width, 'beta_mode': beta_mode,
                          'delta_F_quantiles_025_50_975': np.quantile(delta, [0.025, 0.5, 0.975]).tolist(),
                          'fraction_delta_F_positive': prob, 'monte_carlo_se_fraction': float(np.sqrt(prob*(1-prob)/N))})
            max_residual = max(max_residual, float(np.max(np.abs(r0))), float(np.max(np.abs(r1))))
    checks = []
    for m in (0.1, 1.0, 10.0):
        _, a, _ = state(m)
        def residual(y):
            q = 1/(1+np.exp(-(0.2-0.4*y)))
            z = m*q
            return z/(1+z)*(1-y)*(1+.5*y)-.4*y
        checks.append(abs(float(a)-brentq(residual, 0, 1, xtol=1e-14)))
    assert max_residual < 1e-12 and max(checks) < 1e-12
    # Equal accessible pool z at different gate states and partition factors.
    z = 0.25
    equilibrium_a = brentq(lambda a: z/(1+z)*(1-a)*(1+.5*a)-.4*a, 0, 1)
    twins = []
    for p, q in ((.5, .5), (1., .25)):
        bias = np.log(q/(1-q)) + .4*equilibrium_a
        twins.append({'P': 1., 'p': p, 'q': q, 'b': float(bias), 'z': p*q,
                      'a': equilibrium_a, 'F': equilibrium_a, 's': .26-.08*equilibrium_a})
    output = {'description': 'Assumption-sensitive prior predictive, not a posterior or biological prediction',
              'seed': SEED, 'draws_per_combination': N, 'rna_ratio_scenario': RATIO,
              'parameters': {'kp': 1., 'kd': .4, 'rho': .5, 'alpha': .4, 'b': .2, 'Fp': 1., 'k': 1.},
              'baseline_m_1': {'q': baseline_q, 'a': baseline_a, 'F': baseline_a},
              'fixed_scenarios': scenarios, 'prior_predictive': prior,
              'stationary_nonidentifiability_example': twins,
              'checks': {'max_bisection_residual': max_residual, 'max_difference_from_brentq': max(checks)}}
    (HERE/'results.json').write_text(json.dumps(output, indent=2)+'\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
