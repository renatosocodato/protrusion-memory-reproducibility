"""Quantify recovery in fresh independent synthetic cohorts at the frozen design."""

from collections import Counter
import json
from pathlib import Path
import tempfile

import numpy as np
from scipy.stats import binomtest

from challenge import CASES, VERSION, known_response, plan_at, write_gains, write_responses
from protrusion_memory import mode_test

N_COHORTS = 100


def study():
    roots = np.random.SeedSequence(181009).spawn(len(CASES))
    results = {}
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        for (name, (strength, coupling, mismatch, leak, bound)), seed in zip(CASES.items(), roots):
            plan = plan_at(root/"plan.toml", bound)
            template = known_response(strength, coupling)
            decisions = Counter()
            for cohort in seed.spawn(N_COHORTS):
                gain_rng, test_rng = map(np.random.default_rng, cohort.spawn(2))
                write_gains(root/"gains.csv", gain_rng)
                write_responses(root/"reserve.csv", test_rng, template, mismatch, leak)
                frozen = mode_test.calibrate(plan, root/"gains.csv", VERSION)
                (root/"gains.json").write_text(json.dumps(frozen))
                report = mode_test.score(plan, root/"reserve.csv", root/"gains.json", VERSION)
                decisions[report["decision"]] += 1
            n = decisions["deformation_beyond_leakage"]
            ci = binomtest(n, N_COHORTS).proportion_ci(method="wilson")
            results[name] = {"decisions": dict(decisions),
                             "beyond_leakage_95pct_mc_interval": [ci.low, ci.high]}
    return {"data_kind": "synthetic cohorts", "cohorts_per_case": N_COHORTS, "cases": results}


if __name__ == "__main__":
    print(json.dumps(study(), indent=2, allow_nan=False))
