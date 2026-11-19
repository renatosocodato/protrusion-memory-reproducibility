"""Re-score the existing H1299 date split after a matched acquisition replay."""

import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import tempfile

from protrusion_memory.forecast import Plan, fit, score
from protrusion_memory.portable import simulate

VERSION = "1.6.0"


def acquired_rows(path, sample_s):
    # Use the recorded positive signal as a causal step input. A new segment starts at each gap.
    groups = defaultdict(list)
    with path.open() as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        for row in reader:
            groups[row["preparation"], row["unit"]].append(row)
    pairs = []
    for rows in groups.values():
        rows.sort(key=lambda r: int(r["frame"]))
        segments = []
        for r in rows:
            if not r["actin"]:
                continue
            if not segments or int(r["frame"]) != int(segments[-1][-1]["frame"]) + 1:
                segments.append([])
            segments[-1].append(r)
        for segment in segments:
            if len(segment) < 3:
                continue
            signal = [float(r["actin"]) for r in segment]
            config = {
                "schema": "ode-v2",
                "time": {"stop_s": len(signal)*sample_s, "sample_s": sample_s},
                "states": {"control": 0},
                "equations": {"control": "0"},
                "parameters": {"signal": signal[0]},
                "modulations": [{"parameter": "signal", "operation": "set", "value": value,
                                 "start_s": i*sample_s, "stop_s": (i+1)*sample_s}
                                for i, value in enumerate(signal)],
                "reporters": {"lifeact": {"expression": "signal", "response_s": sample_s,
                                           "initial": signal[0]}},
                "observations": {"lifeact": {"unit": "au", "mix": {"lifeact": 1},
                                             "delay_s": sample_s, "exposure_s": sample_s}},
            }
            obs = simulate(config)["observations"]["lifeact"]
            for t, value in zip(obs["time_s"], obs["values"]):
                index = round(t / sample_s)
                if index < len(segment):
                    point = segment[index]
                    pairs.append((point, {**point, "actin": str(value)}))
    return fields, pairs


def run(plan, calibration, reserved):
    train = acquired_rows(calibration, plan.sample_s)
    test = acquired_rows(reserved, plan.sample_s)
    results = {}
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        for variant, column in (("matched_point", 0), ("coarse_acquisition_stress", 1)):
            for name, (fields, rows) in (("calibration", train), ("reserved", test)):
                with (root/f"{name}.csv").open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=fields)
                    writer.writeheader()
                    writer.writerows(item[column] for item in rows)
            frozen = fit(plan, root/"calibration.csv", VERSION)
            (root/"fit.json").write_text(json.dumps(frozen))
            report = score(plan, root/"reserved.csv", root/"fit.json", VERSION)
            results[variant] = {
                "counts": report["counts"],
                "rmse": {kind: report["candidates"][kind]["scaled_rmse"]
                         for kind in ("motion_only", "actin_modulated")},
                "actin_increment": report["actin_modulated_minus"]["motion_only"],
            }
    assert results["matched_point"]["counts"] == results["coarse_acquisition_stress"]["counts"]
    assert results["matched_point"]["rmse"]["motion_only"] == results["coarse_acquisition_stress"]["rmse"]["motion_only"]
    return {"data_kind": "public H1299 data + assigned observation operators", "sample_s": plan.sample_s,
            "results": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("calibration", type=Path)
    parser.add_argument("reserved", type=Path)
    args = parser.parse_args()
    print(json.dumps(run(Plan(args.plan), args.calibration, args.reserved), indent=2, allow_nan=False))
