"""Calibrate and score two known synthetic edge mechanisms in fresh preparations."""

import argparse
import csv
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from protrusion_memory.forecast import Plan, fit, score

VERSION = "1.5.0"
REPLICATES = 40
PRACTICAL_GAIN = 0.05  # um/min of RMSE, specific to these simulations
TRAIN = tuple(f"cal-{i}" for i in range(6))
TEST = tuple(f"reserved-{i}" for i in range(8))


def schedule(preparation, i):
    if preparation in TRAIN:
        return 10 + i % 3, 24 + i % 3
    return 13 + i % 4, 32 + i % 4


def write_plan(path):
    text = """schema = "edge-forecast-v1"
sample_s = 5
horizon_frames = 1
prefix_frames = 6
motion_delay_frames = 0
actin_delay_frames = 0
motion_unit = "um/min"

[split]
"""
    text += f"calibration = {json.dumps(TRAIN)}\nheldout = {json.dumps(TEST)}\n"
    for split in (TRAIN, TEST):
        for i, preparation in enumerate(split):
            onset, offset = schedule(preparation, i)
            text += f"\n[perturbations.{preparation}]\nonset_frame = {onset}\noffset_frame = {offset}\n"
    path.write_text(text)


def write_split(path, preparations, rng, coupling):
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["preparation", "unit", "region", "frame", "actin", "protrusion"])
        for i, preparation in enumerate(preparations):
            onset, offset = schedule(preparation, i)
            common_noise = rng.normal(0, 0.08, 52)
            for cell in range(2):
                for region in range(2):
                    log_signal = 0.0
                    prev_velocity = 0.0
                    velocity = rng.normal(0, 0.16)
                    edge_offset = rng.normal(0, 0.12)
                    for frame in range(52):
                        pulse = float(onset <= frame < offset)
                        previous_pulse = float(onset <= frame - 1 < offset)
                        log_signal = 0.48 * log_signal + 0.48 * pulse + rng.normal(0, 0.36)
                        writer.writerow([preparation, f"cell-{cell}", f"region-{region}", frame,
                                         math.exp(log_signal), velocity + edge_offset])
                        prev_velocity, velocity = velocity, (
                            0.48 * velocity + 0.17 * prev_velocity + 0.44 * pulse
                            - 0.19 * previous_pulse + coupling * log_signal
                            + common_noise[frame] + rng.normal(0, 0.19)
                        )


def benchmark():
    results = {"data_kind": "synthetic", "replicates_per_mechanism": REPLICATES,
               "practical_gain_um_per_min": PRACTICAL_GAIN, "horizon_s": 5, "cases": {}}
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        write_plan(root / "plan.toml")
        plan = Plan(root / "plan.toml")
        for case, coupling in (("pulse_and_motion", 0.0), ("local_actin_increment", 1.1)):
            increments = []
            for i in range(REPLICATES):
                train_seed, test_seed = np.random.SeedSequence(96021 + i).spawn(2)
                write_split(root / "calibration.csv", TRAIN, np.random.default_rng(train_seed), coupling)
                write_split(root / "heldout.csv", TEST, np.random.default_rng(test_seed), coupling)
                frozen = fit(plan, root / "calibration.csv", VERSION)
                (root / "fit.json").write_text(json.dumps(frozen))
                report = score(plan, root / "heldout.csv", root / "fit.json", VERSION)
                increments.append(report["actin_modulated_minus"]["motion_only"]["rmse_difference"])
                if i == 0:
                    first = report
            results["cases"][case] = {
                "true_coupling": coupling,
                "rmse_difference_um_per_min": increments,
                "median_rmse_difference_um_per_min": float(np.median(increments)),
                "gains_exceeding_practical_threshold": sum(d < -PRACTICAL_GAIN for d in increments),
                "first_replicate": first,
            }
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(benchmark(), stream, indent=2, allow_nan=False)
        stream.write("\n")
