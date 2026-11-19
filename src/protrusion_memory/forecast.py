"""Causal-prefix forecasts for paired actin and protrusion time series."""

from collections import defaultdict
import csv
from dataclasses import dataclass
import io
import json
import math

import numpy as np

from .calibration import digest, load_table, metrics, names, weights
from .portable import number, table

SCHEMA = "forecast-v1"
MODELS = ("stationary", "persistence", "motion_only", "actin_modulated")
COLUMNS = {"preparation", "unit", "frame", "actin", "protrusion"}


@dataclass(frozen=True)
class Origin:
    preparation: str
    unit: str
    frame: int
    target: float
    baseline: float
    motion: float
    actin: float
    channel = "protrusion"


def integer(value, context, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{context} must be an integer >= {minimum}")
    return value


class Plan:
    def __init__(self, path):
        config, self.hash = load_table(path)
        table(config, "forecast plan", fields={"schema", "sample_s", "horizon_frames", "prefix_frames",
                                              "motion_delay_frames", "motion_unit", "split"})
        if config.get("schema") != SCHEMA:
            raise ValueError(f"expected schema {SCHEMA}")
        self.sample_s = number(config["sample_s"], "sample_s")
        if self.sample_s <= 0:
            raise ValueError("sample_s must be positive")
        self.horizon = integer(config["horizon_frames"], "horizon_frames", 1)
        self.prefix = integer(config["prefix_frames"], "prefix_frames", 2)
        self.motion_delay = integer(config["motion_delay_frames"], "motion_delay_frames", 0)
        self.unit = config["motion_unit"]
        if not isinstance(self.unit, str) or not self.unit.strip():
            raise ValueError("motion_unit must name the units of the protrusion column")
        split = table(config["split"], "split", fields={"calibration", "heldout"})
        self.training = names(split["calibration"], "calibration")
        self.heldout = names(split["heldout"], "heldout")
        if set(self.training) & set(self.heldout):
            raise ValueError("forecast preparation splits overlap")

    @property
    def clock(self):
        return {"sample_s": self.sample_s, "horizon_s": self.sample_s * self.horizon,
                "prefix_frames": self.prefix, "motion_delay_frames": self.motion_delay,
                "motion_unit": self.unit}


def observations(path, expected):
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    if set(reader.fieldnames or ()) != COLUMNS or len(reader.fieldnames) != len(COLUMNS):
        raise ValueError("CSV header must contain exactly " + ", ".join(sorted(COLUMNS)))
    by_cell = defaultdict(dict)
    for row in reader:
        if row.keys() != COLUMNS or any(value is None for value in row.values()):
            raise ValueError(f"malformed CSV row {reader.line_num}")
        preparation, unit = row["preparation"], row["unit"]
        if preparation not in expected or not unit or unit != unit.strip():
            raise ValueError(f"unexpected preparation or invalid unit on row {reader.line_num}")
        frame = integer(int(row["frame"]), "frame", 0)
        series = by_cell[preparation, unit]
        if frame in series:
            raise ValueError(f"duplicate frame for {preparation}/{unit}: {frame}")
        actin = number(float(row["actin"]), "actin") if row["actin"] else None
        motion = number(float(row["protrusion"]), "protrusion") if row["protrusion"] else None
        if actin is not None and actin <= 0:
            raise ValueError("recorded actin must be positive for the log transform")
        if motion is not None and motion < 0:
            raise ValueError("recorded protrusion must be nonnegative")
        series[frame] = (math.log(actin) if actin is not None else None, motion)
    if {p for p, _ in by_cell} != set(expected):
        raise ValueError("input must contain all declared preparation IDs for this split")
    return by_cell, digest(raw)


def eligible_origins(by_cell, plan):
    result = []
    for (preparation, unit), series in sorted(by_cell.items()):
        baseline_frames = []
        baseline_end = None
        for frame, (actin, motion) in sorted(series.items()):
            if actin is None or motion is None:
                baseline_frames.clear()
                continue
            if baseline_frames and frame != baseline_frames[-1] + 1:
                baseline_frames.clear()
            baseline_frames.append(frame)
            if len(baseline_frames) == plan.prefix:
                baseline_end = frame
                break
        if baseline_end is None:
            continue
        prefix_actin, prefix_motion = np.mean([series[t] for t in baseline_frames], axis=0)
        for frame, (actin, _) in sorted(series.items()):
            lag, future = frame - plan.motion_delay, frame + plan.horizon
            if (frame <= baseline_end or lag < baseline_end or actin is None
                    or lag not in series or future not in series
                    or not all(t in series for t in range(lag, future + 1))):
                continue
            motion, target = series[lag][1], series[future][1]
            if motion is not None and target is not None:
                result.append(Origin(preparation, unit, frame, target, float(prefix_motion),
                                     float(motion - prefix_motion), float(actin - prefix_actin)))
    absent = {p for p, _ in by_cell} - {o.preparation for o in result}
    if absent:
        raise ValueError("no eligible forecasts in preparations: " + ", ".join(sorted(absent)))
    return result


def _matrix(origins):
    return np.array([(o.motion, o.actin) for o in origins])


def _counts(origins):
    return {"preparations": len({o.preparation for o in origins}),
            "cells": len({(o.preparation, o.unit) for o in origins}), "origins": len(origins)}


def fit(plan, data, version):
    records, source_hash = observations(data, plan.training)
    origins = eligible_origins(records, plan)
    x = _matrix(origins)
    target = np.array([o.target - o.baseline for o in origins])
    weight = np.sqrt(weights(origins))
    fitted, ranks = {}, {}
    for candidate, size in (("motion_only", 1), ("actin_modulated", 2)):
        coef, _, rank, _ = np.linalg.lstsq(x[:, :size] * weight[:, None], target * weight, rcond=None)
        fitted[candidate] = list(map(float, coef))
        ranks[candidate] = int(rank)
    return {"schema": SCHEMA, "software_version": version, "plan_sha256": plan.hash,
            "calibration_sha256": source_hash, "calibration_preparations": plan.training,
            "heldout_preparations_reserved": plan.heldout, "clock": plan.clock,
            "counts": _counts(origins), "coefficients": fitted, "rank": ranks}


def score(plan, data, fit_file, version):
    raw = fit_file.read_bytes()
    frozen = json.loads(raw)
    if (not isinstance(frozen, dict) or frozen.get("schema") != SCHEMA
            or frozen.get("plan_sha256") != plan.hash or frozen.get("software_version") != version
            or frozen.get("calibration_preparations") != plan.training
            or frozen.get("heldout_preparations_reserved") != plan.heldout or frozen.get("clock") != plan.clock):
        raise ValueError("frozen forecast fit does not match the plan, split and software version")
    records, source_hash = observations(data, plan.heldout)
    if source_hash == frozen.get("calibration_sha256"):
        raise ValueError("held-out CSV equals the calibration CSV")
    origins = eligible_origins(records, plan)
    baseline = np.array([o.baseline for o in origins])
    target = np.array([o.target for o in origins])
    x = _matrix(origins)
    predictions = {"stationary": baseline, "persistence": baseline + x[:, 0]}
    for candidate, size in (("motion_only", 1), ("actin_modulated", 2)):
        coefficients = frozen["coefficients"][candidate]
        if not isinstance(coefficients, list) or len(coefficients) != size:
            raise ValueError(f"invalid frozen {candidate} coefficient vector")
        calibrated = [number(value, candidate) for value in coefficients]
        predictions[candidate] = baseline + x[:, :size] @ calibrated
    rows = {}
    totals = {}
    for candidate in MODELS:
        residual = predictions[candidate] - target
        totals[candidate] = metrics(origins, residual, np.ones(len(origins)))
        rows[candidate] = {}
        for p in plan.heldout:
            idx = [i for i, o in enumerate(origins) if o.preparation == p]
            subset = [origins[i] for i in idx]
            summary = metrics(subset, residual[idx], np.ones(len(idx)))
            rows[candidate][p] = {"rmse": summary["scaled_rmse"], **_counts(subset)}
    contrasts = {}
    for reference in ("stationary", "motion_only"):
        paired = {p: rows["actin_modulated"][p]["rmse"] ** 2 - rows[reference][p]["rmse"] ** 2
                  for p in plan.heldout}
        leave_out = []
        if len(plan.heldout) > 1:
            for omit in plan.heldout:
                errors = [math.sqrt(float(np.mean([r["rmse"] ** 2 for p, r in rows[m].items() if p != omit])))
                          for m in ("actin_modulated", reference)]
                leave_out.append(errors[0] - errors[1])
        contrasts[reference] = {
            "rmse_difference": totals["actin_modulated"]["scaled_rmse"] - totals[reference]["scaled_rmse"],
            "paired_mse_difference_by_preparation": paired,
            "leave_one_preparation_out_rmse_difference_range": [min(leave_out), max(leave_out)] if leave_out else None,
        }
    return {"schema": SCHEMA, "software_version": version, "fit_sha256": digest(raw),
            "heldout_sha256": source_hash, "heldout_preparations": plan.heldout, "clock": plan.clock,
            "counts": _counts(origins), "candidates": totals, "by_preparation": rows,
            "actin_modulated_minus": contrasts}
