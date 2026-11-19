"""Project recomputed results onto the published figure inputs and check parity."""
import argparse
import json
import math
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "figures/data"


def get_json(path):
    return json.loads(path.read_text())


def select(observed, template):
    """Keep precisely the keys needed by the artwork."""
    if isinstance(template, dict):
        return {key: select(observed[key], value) for key, value in template.items()}
    if isinstance(template, list):
        return [select(actual, wanted) for actual, wanted in zip(observed, template, strict=True)]
    return observed


def agree(actual, reference, path, atol):
    if isinstance(reference, dict):
        for key in reference:
            agree(actual[key], reference[key], f"{path}/{key}", atol)
    elif isinstance(reference, list):
        for i, (a, b) in enumerate(zip(actual, reference, strict=True)):
            agree(a, b, f"{path}/{i}", atol)
    elif isinstance(reference, (int, float)) and not isinstance(reference, bool):
        if not math.isclose(actual, reference, abs_tol=atol, rel_tol=0):
            raise ValueError(f"{path}: recomputed {actual}, expected {reference} within {atol}")
    elif actual != reference:
        raise ValueError(f"{path}: recomputed {actual!r}, expected {reference!r}")


def write_json(output, name, observed, atol=1e-8):
    reference = get_json(REFERENCE / name)
    result = select(observed, reference)
    agree(result, reference, name, atol)
    (output / name).write_text(json.dumps(result, indent=2) + "\n")


def write_table(output, name, observed):
    reference = pd.read_csv(REFERENCE / name)
    pd.testing.assert_frame_equal(
        observed, reference, check_dtype=False, check_exact=False, rtol=1e-9, atol=1e-10
    )
    observed.to_csv(output / name, index=False)


def check_chaotic_trace(output, trace):
    name = "S02_state.csv"
    reference = pd.read_csv(REFERENCE / name)
    np.testing.assert_allclose(trace["time"], reference["time"], rtol=0, atol=1e-10)
    early = trace["time"] <= 5
    np.testing.assert_allclose(trace.loc[early], reference.loc[early], rtol=1e-9, atol=1e-10)
    states = trace[["q", "a", "m", "x"]].to_numpy()
    if not (np.isfinite(states).all() and (states >= 0).all() and (states <= 1).all()):
        raise ValueError(f"{name}: recomputed bounded projections left their stated domain.")
    for lower, upper in ((5, 100), (100, 200)):
        late = (trace["time"] > lower) & (trace["time"] <= upper)
        columns = ["q", "a", "m", "x"]
        np.testing.assert_allclose(trace.loc[late, columns].std(), reference.loc[late, columns].std(),
                                   rtol=.35, atol=0, err_msg=f"Late variability on ({lower}, {upper}]")
    # Chaotic trajectories can diverge pointwise across numeric environments.
    # Rebuilt figures show the new trajectory on the same times, with exponents checked separately.
    trace.to_csv(output / name, index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("computations", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    calculated, output = args.computations.resolve(), args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = lambda key: get_json(calculated / key)

    joint = "joint_fret_actin_mechanics"
    write_json(output, "classification.json", source(f"{joint}/results_v1/results.json"))
    reporter = source(f"{joint}/independent_gate_readout.json")
    key = lambda row: tuple(row[k] for k in ("tau", "cadence", "bound"))
    by_setting = {key(row): row for row in reporter["rows"]}
    desired = get_json(REFERENCE / "independent_gate_readout.json")["rows"]
    write_json(output, "independent_gate_readout.json",
               {"rows": [by_setting[key(row)] for row in desired]})

    frontier = source("frontier_conformation_chaos/results_v1/results.json")
    write_json(output, "optical_reference.json", frontier["optical_reference_calibration"])
    order = pd.DataFrame(
        [(label, sum(frontier["conditions"][k]["joint"]["correct"]), 2 * frontier["n_parents"])
         for label, k in [("Gate then actin", "same"), ("Actin then gate", "reverse")]],
        columns=["order", "correct", "total"],
    )
    write_table(output, "pulse_order_counts.csv", order)
    chaos = "autonomous_ensemble_chaos/results_v1"
    spectrum = source(f"{chaos}/results.json")["cases"]
    name = "chaos_exponents.json"
    expected = get_json(REFERENCE / name)
    exponents = select({"cases": {
        label: {"spectrum": row["spectrum"][:1]} for label, row in spectrum.items()
    }}, expected)
    for label, row in exponents["cases"].items():
        agree(row, expected["cases"][label], label, .002 if label == "low_control" else .05)
    (output / name).write_text(json.dumps(exponents, indent=2) + "\n")
    with np.load(calculated / chaos / "representative_trace.npz") as raw:
        chosen = raw["t"] <= 200
        u, v, w, x = raw["states"][chosen].T
        sigmoid = lambda z: 1 / (1 + np.exp(-z))
        trace = pd.DataFrame({"time": raw["t"][chosen], "q": sigmoid(u / 20),
                              "a": sigmoid(v / 30), "m": sigmoid((w - 25) / 30), "x": x})
    check_chaotic_trace(output, trace)

    mechanical = "mechanochemical_disorder"
    write_json(output, "S01_amiri_traces.json", source(f"{mechanical}/amiri_summary.json"))
    with np.load(calculated / mechanical / "synthetic_v1/trajectories.npz") as raw:
        chosen = (raw["t"] >= 200) & (raw["t"] <= 320)
        load = pd.concat([
            pd.DataFrame({"drive": drive, "time": raw["t"][chosen] - 200,
                          "load": raw[drive][chosen, 2]})
            for drive in ("autonomous", "periodic", "stochastic")
        ], ignore_index=True)
    write_table(output, "S01_load.csv", load)

    organization = "molecular_organization_protrusion"
    os.environ["PROTRUSION_REPRO_OUTPUT"] = str(calculated)
    os.environ["MPLCONFIGDIR"] = str(calculated / ".matplotlib")
    sys.path[:0] = [str(ROOT / "computation"), str(ROOT / "computation" / organization)]
    from run_comparison import CLASSES, MVALUES, simulate
    rows = []
    for m in MVALUES:
        for label in CLASSES:
            _, history = simulate(np.array([[9., 8., 14., 18., 24., .6]]),
                                  m, label, dt=.5, full=True)
            h = history[0]
            rows.extend((m, label, .5*i, row[2], row[3], row[4]-h[0,4]) for i, row in enumerate(h))
    write_table(output, "S03_trajectory.csv", pd.DataFrame(
        rows, columns=["m", "model", "time", "contrast", "actin", "displacement"]
    ))

    spatial = "spatial_molecular_protrusion/results_v1"
    for name, src in {
        "S04_heldout_predictions.csv": f"{organization}/results_v1/heldout_predictions.csv",
        "S05_second_challenge_parents.csv": f"{spatial}/second_challenge_parents.csv",
        "S05_third_challenge_parents.csv": f"{spatial}/third_challenge_parents.csv",
        "S06_fourth_challenge_parents.csv": f"{spatial}/fourth_challenge_parents.csv",
    }.items():
        columns = pd.read_csv(REFERENCE / name, nrows=0).columns.tolist()
        write_table(output, name, pd.read_csv(calculated / src)[columns])

    ensemble = "conformational_ensemble"
    write_json(output, "S07_conformation.json", {
        "0.3": source(f"{ensemble}/results_v1/results.json")["m"]["0.3"],
        "1.0": source(f"{ensemble}/results_v2_fresh/results.json"),
    })
    print(f"All 14 figure input sets regenerated and checked in {output}")


if __name__ == "__main__":
    main()
