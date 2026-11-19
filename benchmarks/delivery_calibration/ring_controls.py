"""Simulate separate standards, mechanical controls and reserved fixed-drive ring cohorts."""

import argparse
from collections import Counter
import json
from pathlib import Path
import sys
import tempfile

import numpy as np

from protrusion_memory import mode_test

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "mode_intervention"))
from challenge import CALIBRATION, RESERVED, known_response, write_gains, write_responses

VERSION = "2.0.0"
MECHANICAL = [f"blocked-{index}" for index in range(10)]
CASES = {
    # Coupling, mechanical effect in control, effect in reserved, residual actin scale.
    "mediated": (.07, .004, .004, 0),
    "direct_only": (0, .004, .004, 0),
    "transferred_leakage": (0, .24, .24, 0),
    "incomplete_block": (.07, .004, .004, .50),
    "nontransfer": (0, .004, .24, 0),
}


def write_plan(path, mechanical=MECHANICAL):
    path.write_text(
        'schema = "mode-test-v2"\nmode = 2\naxis = [1, 0]\nwindow_s = [10, 30]\n'
        f'[split]\ncalibration = {json.dumps(CALIBRATION)}\nreserved = {json.dumps(RESERVED)}\n'
        f'mechanical = {json.dumps(mechanical)}\n'
        '[criteria]\nactin_min = 0.15\nshape_min = 0.08\ndrive_max = 0.075\ndrive_gain = 0.55\n'
        'blocked_actin_max = 0.15\nblocked_actin_gain = 0.45\ntransfer_margin = 0.020\n'
        'control_shape_error = 0.0\n[resampling]\nseed = 49681\ndraws = 1200\n'
    )
    return mode_test.Plan.load(path)


def example(directory, case, response, blocked, seed):
    plan = write_plan(directory / "plan.toml")
    gain, control, reserve = map(np.random.default_rng, seed.spawn(3))
    gains, controls, reserved = (directory / f"{name}.csv" for name in ("gains", "controls", "reserved"))
    frozen, delivery = (directory / f"{name}.json" for name in ("angular", "mechanical"))
    write_gains(gains, gain)
    write_responses(controls, control, blocked, .012, case[1], MECHANICAL)
    write_responses(reserved, reserve, response, .012, case[2])
    frozen.write_text(json.dumps(mode_test.calibrate(plan, gains, VERSION)))
    delivery.write_text(json.dumps(mode_test.calibrate_delivery(plan, controls, frozen, VERSION)))
    return mode_test.score(plan, reserved, frozen, VERSION, delivery)


def benchmark(directory, cohorts):
    result = {"data_kind": "synthetic", "cohorts_per_case": cohorts,
              "preparations_per_cohort": {"angular": 8, "mechanical": 10, "reserved": 12}, "cases": {}}
    rep_seeds = np.random.SeedSequence(449162).spawn(len(CASES))
    trial_seeds = np.random.SeedSequence(449163).spawn(len(CASES))
    for (name, case), rep_seed, trial_seed in zip(CASES.items(), rep_seeds, trial_seeds):
        response = known_response(1, case[0])
        blocked = known_response(case[3], case[0])
        files = directory / name
        files.mkdir()
        representative = example(files, case, response, blocked, rep_seed)
        counts, interval_covers, qualified = Counter(), 0, 0
        with tempfile.TemporaryDirectory() as scratch:
            for seed in trial_seed.spawn(cohorts):
                report = example(Path(scratch), case, response, blocked, seed)
                counts[report["decision"]] += 1
                interval = report["shape_interval_with_leakage"]
                if interval is not None:
                    qualified += 1
                    interval_covers += int(interval[0] <= response[1] <= interval[1])
        result["cases"][name] = {"planted_actin_deformation_um": float(response[1]),
                                 "mechanical_odd_control_um": case[1], "mechanical_odd_reserved_um": case[2],
                                 "counts": dict(counts), "interval_covers_latent": interval_covers,
                                 "qualified_intervals": qualified, "representative": representative}
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, required=True, help="new directory for illustrative inputs")
    parser.add_argument("--cohorts", type=int, default=60)
    args = parser.parse_args()
    if args.cohorts < 1:
        parser.error("cohorts must be positive")
    args.datasets.mkdir()
    print(json.dumps(benchmark(args.datasets, args.cohorts), indent=2, allow_nan=False))
