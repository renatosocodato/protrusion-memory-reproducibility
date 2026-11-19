"""Stress the acquisition map on one public paired actin/membrane recording."""

import argparse
from collections import defaultdict
import csv
from decimal import Decimal
import json
from pathlib import Path
import tempfile

import numpy as np

from protrusion_memory.portable import simulate

from guv_g01 import extract, SAMPLE_S

CASES = {
    "point": (0, 0, 0, 0),
    "reporter": (1, 0, 0, 0),
    "exposure": (0, 1, 0, 0),
    "delay": (0, 0, 1, 0),
    "mixing": (0, 0, 0, 0.1),
    "combined": (1, 1, 1, 0.1),
}


def replay(segment, design, baseline):
    tau, exposure, delay, mixing = design
    # A sample becomes available at its frame and is held until the next frame.
    step = Decimal(str(SAMPLE_S))
    times = [float(i * step) for i in range(len(segment) + 1)]
    signals = np.asarray([[float(row[k]) for k in ("actin_au", "membrane_au")]
                          for row in segment]) / baseline - 1
    modulations = [{"parameter": channel, "start_s": times[i], "stop_s": times[i+1],
                    "operation": "set", "value": float(value)}
                   for i, values in enumerate(signals)
                   for channel, value in zip(("a", "m"), values)]
    config = {
        "schema": "ode-v2",
        "time": {"stop_s": times[-1], "sample_s": SAMPLE_S},
        "states": {"unused": 0},
        "parameters": {"a": 0.0, "m": 0.0},
        "equations": {"unused": "0"},
        "modulations": modulations,
        "reporters": {"actin": {"expression": "a"}, "membrane": {"expression": "m"}},
        "observations": {"actin": {"unit": "fraction", "mix": {"actin": 1, "membrane": mixing},
                                   "exposure_s": exposure * SAMPLE_S, "delay_s": delay * SAMPLE_S}},
    }
    if tau:
        config["reporters"]["actin"].update(response_s=tau * SAMPLE_S, initial=float(signals[0, 0]))
    obs = simulate(config)["observations"]["actin"]
    return {int(segment[0]["frame"]) + round(t / SAMPLE_S): x
            for t, x in zip(obs["time_s"], obs["values"]) if t < times[-1]}


def contiguous(rows):
    segments = []
    for row in rows:
        if not row["actin_au"] or not row["membrane_au"]:
            continue
        if not segments or int(row["frame"]) != int(segments[-1][-1]["frame"]) + 1:
            segments.append([])
        segments[-1].append(row)
    return segments


def summarize(rows):
    by_region = defaultdict(list)
    for row in rows:
        by_region[row["region"]].append(row)
    sector = []
    for region, series in sorted(by_region.items()):
        prefix = [r for r in series if int(r["frame"]) < 4 and r["actin_au"] and r["membrane_au"]]
        if len(prefix) != 4:
            continue
        baseline = np.mean([[float(r[key]) for key in ("actin_au", "membrane_au")] for r in prefix], axis=0)
        acquired = {name: {} for name in CASES}
        for segment in contiguous(series):
            if len(segment) < 4:
                continue
            for name, design in CASES.items():
                acquired[name].update(replay(segment, design, baseline))
        velocity = {int(r["frame"]): float(r["radial_velocity_px_per_min"]) for r in series
                    if r["radial_velocity_px_per_min"]}
        common = sorted(set.intersection(*(set(x) for x in acquired.values()))
                        & {t - 1 for t in velocity if t - 1 >= 7})
        if len(common) < 3:
            continue
        target = [velocity[t + 1] for t in common]
        if np.std(target) == 0:
            continue
        correlations = {}
        for name, traces in acquired.items():
            signal = [traces[t] for t in common]
            correlations[name] = float(np.corrcoef(signal, target)[0, 1]) if np.std(signal) else None
        sector.append({"region": region, "origins": len(common), "correlations": correlations})
    median = {name: float(np.median([s["correlations"][name] for s in sector
                                    if s["correlations"][name] is not None])) for name in CASES}
    return {"data_kind": "public one-vesicle data + assigned observation operators", "sample_s": SAMPLE_S,
            "preparations": 1, "sector_summaries": sector, "median_sector_correlation": median}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("membrane", type=Path)
    parser.add_argument("actin", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "paired.csv"
        extract(args.membrane, args.actin, path, include_membrane=True)
        with path.open() as stream:
            result = summarize(list(csv.DictReader(stream)))
    print(json.dumps(result, indent=2, allow_nan=False))
