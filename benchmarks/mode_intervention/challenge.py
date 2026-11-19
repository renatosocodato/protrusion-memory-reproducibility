"""Synthetic independent gain calibration and signed actin intervention on ring mode 2."""

import argparse
import csv
import json
from pathlib import Path
import tempfile

import numpy as np

from protrusion_memory import mode_test
from protrusion_memory.spatial import simulate

VERSION = "1.9.0"
CALIBRATION = [f"gain-{i}" for i in range(8)]
RESERVED = [f"reserved-{i}" for i in range(12)]
ARM_ORDER = mode_test.ARMS
GAINS = np.array([[[1.12,-.17],[.15,.92]],
                  [[.97,.09],[-.12,1.08]],
                  [[1.03,-.08],[.06,.89]]])
CASES = {
    # actin change scale, mechanical coupling, upstream mismatch, delivery effect, assumed bound
    "mediated": (1, .07, .012, .004, .025),
    "direct_only": (1, 0, .012, .004, .025),
    "bounded_delivery": (1, 0, .012, .12, .16),
    "drive_mismatch": (1, 0, .13, .004, .025),
    "weak_actin": (.04, .07, .012, .004, .025),
    "hidden_odd_leakage": (1, 0, .012, .24, .025),
}


def plan_at(path, limit):
    path.write_text(
        'schema = "mode-test-v1"\nmode = 2\naxis = [1, 0]\nwindow_s = [10, 30]\n'
        f'[split]\ncalibration = {json.dumps(CALIBRATION)}\nreserved = {json.dumps(RESERVED)}\n'
        f'[criteria]\nactin_min = 0.15\nshape_min = 0.08\ndrive_max = 0.075\n'
        f'drive_gain = 0.55\ndelivery_bound = {limit}\n'
        '[resampling]\nseed = 99685\ndraws = 1500\n'
    )
    return mode_test.Plan.load(path)


def known_response(strength, coupling):
    def one(second_input):
        config = {
            "schema": "ring-v1", "ring": {"nodes": 64, "modes": [2]},
            "time": {"stop_s": 32, "sample_s": 1},
            "parameters": {"u": 0.8, "r": 0, "g": coupling},
            "modulations": [{"parameter": "r", "start_s": 10, "stop_s": 15,
                             "operation": "set", "value": second_input}],
            "inputs": {"stimulus": {"parameter": "u", "cos": {"2": 1}},
                       "actin_pulse": {"parameter": "r", "cos": {"2": 1}}},
            "fields": {
                "actin": {"decay_s_inv": .125, "diffusion_rad2_per_s": .015,
                          "drives": {"stimulus": 1, "actin_pulse": 1}},
                "shape": {"decay_s_inv": .05, "diffusion_rad2_per_s": .002,
                          "drives": {"stimulus": .03}, "coupling": {"actin": "g"}},
            },
            "reporters": {"actin": {"field": "actin", "mode": 2, "response_s": 4, "initial": 0},
                          "shape": {"field": "shape", "mode": 2}},
            "observations": {
                name: {"unit": unit, "mix": {name: 1}, "exposure_s": .7, "delay_s": .6}
                for name, unit in (("actin", "au"), ("shape", "um"))
            },
        }
        readouts = simulate(config)["observations"]
        return np.asarray([np.mean([y for t,y in zip(item["time_s"],item["values"])
                                   if 10 <= t < 30]) for item in readouts.values()])
    return one(.35*strength)-one(0)


def write_gains(path, rng):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["preparation","mode","channel","known_cos","known_sin","read_cos","read_sin"])
        for prep in CALIBRATION:
            for channel, gain in zip(mode_test.CHANNELS, GAINS):
                actual = gain + rng.normal(0,.025,(2,2))
                for basis in ([1,0],[-1,0],[0,1],[0,-1],[.6,.8],[-.8,.6]):
                    measured = actual @ basis + rng.normal(0,.01,2)
                    writer.writerow([prep,2,channel,*basis,*measured])


def write_responses(path, rng, response, mismatch, leakage, preparations=RESERVED):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["preparation","mode","arm",*(f"{c}_{a}" for c in mode_test.CHANNELS
                                                       for a in ("cos","sin"))])
        for prep in preparations:
            individual_gain = GAINS + rng.normal(0,.025,GAINS.shape)
            base = rng.normal(0,.1,(3,2))
            base[2,0] += .8
            effect = response * rng.normal(1,.08)
            for arm in ARM_ORDER:
                nominal = 1 if arm.endswith("plus") else -1
                values = base.copy()
                values[1] += nominal*np.array([.025,.009])  # shared delivery offset
                if arm.startswith("active"):
                    values[0,0] += nominal*effect[0]
                    values[1,0] += nominal*(effect[1]+.4*mismatch+leakage)
                    values[2,0] += nominal*mismatch
                observed = (individual_gain @ values[...,None])[...,0] + rng.normal(0,[.028,.025,.004],(2,3)).T
                writer.writerow([prep,2,arm,*observed.ravel()])


def benchmark(output=None):
    out = {"data_kind": "synthetic calibration and independently generated reserved cohorts", "cases": {}}
    root_sequence = np.random.SeedSequence(91723)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) if output is None else output
        for case, (scale, coupling, mismatch, leak, bound) in CASES.items():
            path = root/case
            path.mkdir()
            plan = plan_at(path/"plan.toml", bound)
            response = known_response(scale, coupling)
            gain_rng, test_rng = map(np.random.default_rng, root_sequence.spawn(2))
            write_gains(path/"gains.csv", gain_rng)
            write_responses(path/"reserved.csv", test_rng, response, mismatch, leak)
            frozen = path/"calibration.json"
            frozen.write_text(json.dumps(mode_test.calibrate(plan,path/"gains.csv",VERSION)))
            report = mode_test.score(plan,path/"reserved.csv",frozen,VERSION)
            out["cases"][case] = {"latent_response_to_actin": response.tolist(),
                                 "assumed_delivery_bound": bound, **report}
    return out


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", type=Path, help="write synthetic input CSVs in a new directory")
    args = parser.parse_args()
    if args.datasets:
        args.datasets.mkdir()
    print(json.dumps(benchmark(args.datasets), indent=2, allow_nan=False))
