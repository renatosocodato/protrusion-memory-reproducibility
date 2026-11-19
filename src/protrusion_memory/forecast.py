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
EDGE_SCHEMA = "edge-forecast-v1"
MODELS = ("stationary", "persistence", "motion_only", "actin_modulated")
COLUMNS = {"preparation", "unit", "frame", "actin", "protrusion"}
EDGE_COLUMNS = COLUMNS | {"region"}


@dataclass(frozen=True)
class Origin:
    preparation: str
    unit: str
    frame: int
    target: float
    baseline: float
    motion: float
    actin: float
    channel: str = "protrusion"
    context: tuple = ()


def integer(value, context, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{context} must be an integer >= {minimum}")
    return value


class Plan:
    def __init__(self, path):
        config, self.hash = load_table(path)
        self.schema = config.get("schema")
        self.edge = self.schema == EDGE_SCHEMA
        fields = {"schema", "sample_s", "horizon_frames", "prefix_frames",
                  "motion_delay_frames", "motion_unit", "split"}
        if self.edge:
            fields |= {"actin_delay_frames", "perturbations"}
        table(config, "forecast plan", fields=fields)
        if self.schema not in (SCHEMA, EDGE_SCHEMA):
            raise ValueError(f"expected schema {SCHEMA} or {EDGE_SCHEMA}")
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

        if self.edge:
            if min(len(self.training), len(self.heldout)) < 2:
                raise ValueError("edge forecasting requires at least two preparations in each split")
            self.actin_delay = integer(config["actin_delay_frames"], "actin_delay_frames", 0)
            perturbations = config["perturbations"]
            table(perturbations, "perturbations", fields=set(self.training + self.heldout))
            self.perturbations = {}
            for name in self.training + self.heldout:
                pulse = table(perturbations[name], name, fields={"onset_frame", "offset_frame"})
                onset = integer(pulse["onset_frame"], f"{name} onset_frame", self.prefix)
                offset = integer(pulse["offset_frame"], f"{name} offset_frame", onset + 1)
                self.perturbations[name] = onset, offset

    @property
    def clock(self):
        clock = {"sample_s": self.sample_s, "horizon_s": self.sample_s * self.horizon,
                 "prefix_frames": self.prefix, "motion_delay_frames": self.motion_delay,
                 "motion_unit": self.unit}
        if self.edge:
            clock["actin_delay_frames"] = self.actin_delay
        return clock

    @property
    def features(self):
        if self.edge:
            return ("last_motion", "previous_motion", "pulse_active", "pulse_previous", "local_log_actin")
        return ("last_motion", "log_actin")


def observations(path, expected, edge=False):
    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    columns = EDGE_COLUMNS if edge else COLUMNS
    if set(reader.fieldnames or ()) != columns or len(reader.fieldnames) != len(columns):
        raise ValueError("CSV header must contain exactly " + ", ".join(sorted(columns)))
    by_cell = defaultdict(dict)
    for row in reader:
        if row.keys() != columns or any(value is None for value in row.values()):
            raise ValueError(f"malformed CSV row {reader.line_num}")
        preparation, unit = row["preparation"], row["unit"]
        if preparation not in expected or not unit or unit != unit.strip():
            raise ValueError(f"unexpected preparation or invalid unit on row {reader.line_num}")
        frame = integer(int(row["frame"]), "frame", 0)
        key = (preparation, unit)
        if edge:
            region = row["region"]
            if not region or region != region.strip():
                raise ValueError(f"invalid edge region on row {reader.line_num}")
            key += (region,)
        series = by_cell[key]
        if frame in series:
            raise ValueError(f"duplicate frame for {'/'.join(key)}: {frame}")
        actin = number(float(row["actin"]), "actin") if row["actin"] else None
        motion = number(float(row["protrusion"]), "protrusion") if row["protrusion"] else None
        if actin is not None and actin <= 0:
            raise ValueError("recorded actin must be positive for the log transform")
        if not edge and motion is not None and motion < 0:
            raise ValueError("recorded protrusion must be nonnegative")
        series[frame] = (math.log(actin) if actin is not None else None, motion)
    if {key[0] for key in by_cell} != set(expected):
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
    absent = {key[0] for key in by_cell} - {o.preparation for o in result}
    if absent:
        raise ValueError("no eligible forecasts in preparations: " + ", ".join(sorted(absent)))
    return result


def edge_origins(records, plan):
    origins = []
    for (preparation, unit, region), series in sorted(records.items()):
        onset, offset = plan.perturbations[preparation]
        prefix = [series.get(frame, (None, None)) for frame in range(plan.prefix)]
        if any(None in values for values in prefix):
            continue
        base_actin, base_motion = map(float, np.mean(prefix, axis=0))
        for frame in sorted(series):
            lag = frame - plan.motion_delay
            actin_frame = frame - plan.actin_delay
            future = frame + plan.horizon
            first = min(lag - 1, actin_frame)
            if (frame < onset or first < plan.prefix - 1
                    or any(t not in series for t in range(first, future + 1))):
                continue
            actin = series[actin_frame][0]
            recent, previous, target = (series[t][1] for t in (lag, lag - 1, future))
            if None in (actin, recent, previous, target):
                continue
            context = (previous - base_motion, float(onset <= frame < offset),
                       float(onset <= frame - 1 < offset))
            origins.append(Origin(preparation, unit, frame, target, base_motion,
                                  recent - base_motion, actin - base_actin, region, context))
    absent = {key[0] for key in records} - {o.preparation for o in origins}
    if absent:
        raise ValueError("no eligible edge forecasts in preparations: " + ", ".join(sorted(absent)))
    return origins


def _edge_coverage(records, origins):
    return {"recorded_regions": len(records),
            "excluded_regions": len(records) - _counts(origins)["regions"]}


def _matrix(origins):
    return np.array([(o.motion, *o.context, o.actin) for o in origins])


def _counts(origins):
    return {"preparations": len({o.preparation for o in origins}),
            "cells": len({(o.preparation, o.unit) for o in origins}), "origins": len(origins),
            **({"regions": len({(o.preparation, o.unit, o.channel) for o in origins})}
               if origins[0].context else {})}


def fit(plan, data, version):
    records, source_hash = observations(data, plan.training, plan.edge)
    origins = edge_origins(records, plan) if plan.edge else eligible_origins(records, plan)
    x = _matrix(origins)
    target = np.array([o.target - o.baseline for o in origins])
    weight = np.sqrt(weights(origins))
    fitted, ranks = {}, {}
    for candidate, size in (("motion_only", len(plan.features) - 1), ("actin_modulated", len(plan.features))):
        coef, _, rank, _ = np.linalg.lstsq(x[:, :size] * weight[:, None], target * weight, rcond=None)
        if plan.edge and rank != size:
            raise ValueError(f"{candidate}: calibration cannot resolve all declared features")
        fitted[candidate] = list(map(float, coef))
        ranks[candidate] = int(rank)
    return {"schema": plan.schema, "software_version": version, "plan_sha256": plan.hash,
            "calibration_sha256": source_hash, "calibration_preparations": plan.training,
            "heldout_preparations_reserved": plan.heldout, "clock": plan.clock,
            "counts": _counts(origins), "coefficients": fitted, "rank": ranks,
            **({"features": list(plan.features), "coverage": _edge_coverage(records, origins)} if plan.edge else {})}


def score(plan, data, fit_file, version):
    raw = fit_file.read_bytes()
    frozen = json.loads(raw)
    if (not isinstance(frozen, dict) or frozen.get("schema") != plan.schema
            or frozen.get("plan_sha256") != plan.hash or frozen.get("software_version") != version
            or frozen.get("calibration_preparations") != plan.training
            or frozen.get("heldout_preparations_reserved") != plan.heldout or frozen.get("clock") != plan.clock):
        raise ValueError("frozen forecast fit does not match the plan, split and software version")
    if plan.edge and frozen.get("features") != list(plan.features):
        raise ValueError("frozen edge forecast feature order does not match the plan")
    records, source_hash = observations(data, plan.heldout, plan.edge)
    if source_hash == frozen.get("calibration_sha256"):
        raise ValueError("held-out CSV equals the calibration CSV")
    origins = edge_origins(records, plan) if plan.edge else eligible_origins(records, plan)
    baseline = np.array([o.baseline for o in origins])
    target = np.array([o.target for o in origins])
    x = _matrix(origins)
    predictions = {"stationary": baseline, "persistence": baseline + x[:, 0]}
    for candidate, size in (("motion_only", len(plan.features) - 1), ("actin_modulated", len(plan.features))):
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
    references = ("motion_only",) if plan.edge else ("stationary", "motion_only")
    for reference in references:
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
    return {"schema": plan.schema, "software_version": version, "fit_sha256": digest(raw),
            "heldout_sha256": source_hash, "heldout_preparations": plan.heldout, "clock": plan.clock,
            "counts": _counts(origins), "candidates": totals, "by_preparation": rows,
            "actin_modulated_minus": contrasts,
            **({"coverage": _edge_coverage(records, origins)} if plan.edge else {})}
