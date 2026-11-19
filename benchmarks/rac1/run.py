"""Evaluate the portable solver against the external fibroblast velocity recordings."""

import argparse
from collections import defaultdict
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import tomllib
from urllib.request import urlopen

import numpy as np

from protrusion_memory.portable import simulate

ROOT = Path(__file__).resolve().parents[2]
COMMIT = "a78f2d8d46d8e03a85cf8d568b5ace10fbaa8e72"
SOURCE = f"https://raw.githubusercontent.com/DanuserLab/Protrusion_Retraction_Database/{COMMIT}/rawData"
INPUTS = {
    "GithubExampleData.txt": "2e47b30bfe8026bcede354c9f2cac2e9fa6a295ee7966704e101b569dd91a4b1",
    "PA_Info.txt": "b9e31c0d7448bdea1e621f74d1d2d793004a7f9dce64c991f29938c478114aec",
}
GROUPS = {"PA_Rac1": set(range(1, 8)), "Control_Rac1": set(range(1, 7))}


def read(name, directory):
    if directory is None:
        with urlopen(f"{SOURCE}/{name}", timeout=45) as stream:
            raw = stream.read(15_000_000)
    else:
        raw = (directory / name).read_bytes()
    if sha256(raw).hexdigest() != INPUTS[name]:
        raise ValueError(f"source checksum mismatch for {name}")
    return csv.DictReader(io.StringIO(raw.decode("utf-8")))


def traces(directory):
    illuminated = {(row["Label"], row["cellIndex"], row["windowIndex"])
                   for row in read("PA_Info.txt", directory) if row["PA_window"] == "1"}
    sums = defaultdict(float)
    counts = defaultdict(int)
    for row in read("GithubExampleData.txt", directory):
        frame = int(row["frameIndex"])
        window = (row["Label"], row["cellIndex"], row["windowIndex"])
        if row["layerIndex"] != "1" or frame <= 1 or window not in illuminated or row["Ch0"] == "NA":
            continue
        key = (row["Label"], int(row["cellIndex"]), frame)
        sums[key] += float(row["Ch0"])
        counts[key] += 1
    by_cell = defaultdict(dict)
    for (label, cell, frame), total in sorted(sums.items()):
        by_cell[label, cell][frame] = total / counts[label, cell, frame]
    for group, cells in GROUPS.items():
        for epoch in ("Before", "After"):
            observed = {cell for label, cell in by_cell if label == f"{group}_{epoch}"}
            if observed != cells:
                raise ValueError(f"{group} {epoch} cell identifiers do not match the source design")
    return by_cell


def run(directory):
    with (ROOT / "examples/rac1-edge.toml").open("rb") as stream:
        config = tomllib.load(stream)
    response = simulate(config)["observations"]["edge_velocity"]
    times = np.array(response["time_s"])
    template = np.array(response["values"])
    tau = config["parameters"]["relaxation_s"]
    step = config["time"]["sample_s"]
    target = 1 + (tau / step) * (np.exp(-times / tau) - np.exp(-(times - step) / tau))
    error = float(np.max(abs(template - target)))
    if error > 1e-7:
        raise ValueError("portable solver disagrees with the independent frame-average calculation")

    observed = traces(directory)

    def align(group, cell):
        before = np.mean(list(observed[f"{group}_Before", cell].values()))
        after = observed[f"{group}_After", cell]
        frames = list(after)
        if not set(frames) <= set(range(2, len(times) + 2)):
            raise ValueError("frames exceed the configured activation epoch")
        return template[np.array(frames) - 2], np.array(list(after.values())) - before

    calibration_cells = [1, 3, 5, 7]
    held_out_cells = [2, 4, 6]
    fitted = []
    for cell in calibration_cells:
        shape, velocity = align("PA_Rac1", cell)
        fitted.append(float(shape @ velocity / (shape @ shape)))
    amplitude = float(np.mean(fitted))

    # Invoke the same declarative platform with the actual calibrated parameter.
    config["parameters"]["gain_nm_s"] = amplitude
    calibrated = np.array(simulate(config)["observations"]["edge_velocity"]["values"])
    if not np.allclose(calibrated, amplitude * template, rtol=1e-8, atol=1e-8):
        raise ValueError("the calibrated portable model failed its linear scaling check")
    results = []
    for cell in held_out_cells:
        _, velocity = align("PA_Rac1", cell)
        frames = list(observed["PA_Rac1_After", cell])
        predicted = calibrated[np.array(frames) - 2]
        results.append({
            "cell": cell, "frames": len(velocity), "observed_mean_change_nm_s": float(velocity.mean()),
            "model_mse_nm2_s2": float(np.mean((velocity - predicted) ** 2)),
            "stationary_mse_nm2_s2": float(np.mean(velocity ** 2)),
        })
    controls = {cell: float(align("Control_Rac1", cell)[1].mean()) for cell in held_out_cells}
    model_rmse = float(np.sqrt(np.mean([row["model_mse_nm2_s2"] for row in results])))
    stationary_rmse = float(np.sqrt(np.mean([row["stationary_mse_nm2_s2"] for row in results])))
    return {
        "study_doi": "10.7554/eLife.57681", "source_commit": COMMIT, "source_sha256": INPUTS,
        "maximum_frame_average_solver_error_nm_s": error, "calibration_pa_cells": calibration_cells,
        "calibration_gains_nm_s": fitted, "fitted_gain_nm_s": amplitude,
        "held_out_pa_cells": results, "held_out_light_control_mean_changes_nm_s": controls,
        "held_out_model_rmse_nm_s": model_rmse, "held_out_stationary_rmse_nm_s": stationary_rmse,
        "model_rmse_minus_stationary_nm_s": model_rmse - stationary_rmse,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, help="directory containing the pinned original public files")
    parser.add_argument("--output", type=Path, help="new report (default: stdout)")
    args = parser.parse_args()
    report = json.dumps(run(args.data_dir), indent=2, allow_nan=False) + "\n"
    if args.output is None:
        print(report, end="")
    else:
        with args.output.open("x", encoding="utf-8") as stream:
            stream.write(report)


if __name__ == "__main__":
    main()
