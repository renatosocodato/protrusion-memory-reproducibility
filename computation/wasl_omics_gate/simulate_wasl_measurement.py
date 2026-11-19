#!/usr/bin/env python3
"""Synthetic relative RNA-to-protein sensitivity; the numerical atlas reference cancels."""
import json
from pathlib import Path
from reproduction_io import output_directory
import numpy as np
from run_wasl_sensitivity import state

HERE = output_directory(Path(__file__))
SEED = 202611180502
N = 20000
R_ATLAS = 34.7


def main():
    rng = np.random.default_rng(SEED)
    cases = []
    for geometric_sd in (2, 10):
        r = R_ATLAS*np.exp(np.log(geometric_sd)*rng.normal(size=N))
        pref = 10**rng.uniform(-1, 1, N)
        beta = rng.uniform(0, 1, N)
        ptotal = pref*(r/R_ATLAS)**beta*np.exp(np.log(2)*rng.normal(size=N))
        partition = 10**rng.uniform(-2, 0, N)
        m = ptotal*partition/.1
        y = ptotal*np.exp(np.log(1.2)*rng.normal(size=N))
        q, actin, residual = state(m)
        assert np.isfinite(np.stack((r, ptotal, partition, m, y, q, actin))).all()
        assert np.max(np.abs(residual)) < 1e-12
        assert np.all((q >= 0) & (q <= 1) & (actin >= 0) & (actin <= 1))
        values = {'normalized_rna_ntpm': r, 'relative_total_protein': ptotal,
                  'simulated_relative_protein_readout': y, 'effective_local_pool_m': m,
                  'gate_fraction_q': q, 'actin_a_and_load_F': actin}
        cases.append({'assumed_rna_geometric_sd': geometric_sd,
                      'quantiles_025_50_975': {key: np.quantile(val, [.025, .5, .975]).tolist()
                                                for key, val in values.items()},
                      'max_residual': float(np.max(np.abs(residual)))})
    result = {'kind': 'Synthetic prior-predictive protein readout; not an actual protein observation or posterior',
              'source_human_atlas_ntpm': R_ATLAS, 'seed': SEED, 'draws_per_case': N,
              'cases': cases}
    (HERE/'measurement_results.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
