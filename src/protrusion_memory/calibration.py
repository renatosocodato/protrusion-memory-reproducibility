"""Fit declarative candidates and score them on disjoint preparations."""

from collections import Counter, defaultdict
from copy import deepcopy
import csv
from dataclasses import dataclass
from hashlib import sha256
import io
import json
import math
import tomllib

import numpy as np
from scipy.optimize import least_squares

from .portable import ODE_SCHEMAS, number, simulate, table

SCHEMA = "calibration-v1"
MODELS = ("stationary", "modulated")
CSV_FIELDS = {"preparation", "unit", "scenario", "channel", "time_s", "value", "measurement_unit"}


@dataclass(frozen=True, order=True)
class Row:
    preparation: str
    unit: str
    scenario: str
    channel: str
    time_s: float
    value: float


def digest(raw):
    return sha256(raw).hexdigest()


def load_table(path):
    raw = path.read_bytes()
    return tomllib.loads(raw.decode("utf-8")), digest(raw)


def names(items, context):
    if (not isinstance(items, list) or not items
            or any(not isinstance(item, str) or not item.strip() or item != item.strip() for item in items)):
        raise ValueError(f"{context} must be a nonempty list of preparation IDs")
    if len(items) != len(set(items)):
        raise ValueError(f"{context} contains duplicate preparation IDs")
    return sorted(items)


class Plan:
    def __init__(self, path):
        config, plan_hash = load_table(path)
        table(config, "plan", fields={"schema", "split", "scenarios", "channels", "models"})
        if config.get("schema") != SCHEMA:
            raise ValueError(f"expected schema {SCHEMA}")
        self.identity = {"plan_sha256": plan_hash, "model_sha256": {}}
        split = table(config["split"], "split", fields={"calibration", "heldout"})
        self.training = names(split["calibration"], "calibration")
        self.heldout = names(split["heldout"], "heldout")
        if set(self.training) & set(self.heldout):
            raise ValueError("calibration and held-out preparation IDs overlap")

        self.channels = {}
        for channel, item in table(config["channels"], "channels").items():
            table(item, channel, fields={"unit", "scale"})
            scale = number(item["scale"], f"{channel} scale")
            unit = item["unit"]
            if not isinstance(unit, str) or not unit.strip() or scale <= 0:
                raise ValueError(f"{channel} requires a unit and a positive fixed scale")
            self.channels[channel] = (unit, scale)

        self.scenarios = {}
        for scenario, item in table(config["scenarios"], "scenarios").items():
            if not isinstance(item, dict):
                raise ValueError(f"scenario {scenario} must be a table of parameter overrides")
            self.scenarios[scenario] = {k: number(v, f"{scenario}.{k}") for k, v in item.items()}

        models = table(config["models"], "models", fields=set(MODELS))
        self.models = {}
        self.specs = {}
        for candidate in MODELS:
            item = table(models[candidate], candidate, fields={"path", "fit"})
            model, model_hash = load_table(path.parent / item["path"])
            if model.get("schema") not in ODE_SCHEMAS:
                raise ValueError(f"{candidate} must use the ode-v1 or ode-v2 schema")
            self.identity["model_sha256"][candidate] = model_hash
            self.models[candidate] = model
            specs = item.get("fit", {})
            parameters = model.get("parameters", {})
            if not isinstance(specs, dict) or specs.keys() - parameters.keys():
                raise ValueError(f"{candidate} fitted parameters must be declared in its model")
            limits = {}
            for key, bound in sorted(specs.items()):
                table(bound, key, fields={"initial", "lower", "upper"})
                initial, lower, upper = (number(bound[k], f"{key}.{k}") for k in ("initial", "lower", "upper"))
                if not lower < upper or not lower <= initial <= upper:
                    raise ValueError(f"{key} requires lower < upper and lower <= initial <= upper")
                limits[key] = (initial, lower, upper)
            self.specs[candidate] = limits
            for scenario, overrides in self.scenarios.items():
                if overrides.keys() - parameters.keys():
                    raise ValueError(f"{candidate}: unknown parameter in scenario {scenario}")
                if overrides.keys() & limits.keys():
                    raise ValueError(f"{candidate}: scenario {scenario} overrides a fitted parameter")
        for channel, (unit, _) in self.channels.items():
            laws = [self.models[candidate].get("observations", {}).get(channel, {}) for candidate in MODELS]
            if any(law.get("unit") != unit for law in laws):
                raise ValueError(f"{channel}: both models must declare the expected measurement unit")
            if laws[0].get("operator", "point") != laws[1].get("operator", "point"):
                raise ValueError(f"{channel}: both models must use the same observation operator")


def load_rows(path, plan, expected):
    raw = path.read_bytes()
    records = csv.DictReader(io.StringIO(raw.decode("utf-8-sig"), newline=""))
    if set(records.fieldnames or ()) != CSV_FIELDS or len(records.fieldnames) != len(CSV_FIELDS):
        raise ValueError("CSV header must contain exactly " + ", ".join(sorted(CSV_FIELDS)))
    rows, seen = [], set()
    for record in records:
        if set(record) != CSV_FIELDS or any(value is None or not value.strip() for value in record.values()):
            raise ValueError(f"missing or extra CSV value on line {records.line_num}")
        p, u, s, c = (record[key] for key in ("preparation", "unit", "scenario", "channel"))
        if p not in expected:
            raise ValueError(f"unexpected preparation in this split: {p}")
        if s not in plan.scenarios or c not in plan.channels:
            raise ValueError(f"unknown scenario or channel on line {records.line_num}")
        if record["measurement_unit"] != plan.channels[c][0]:
            raise ValueError(f"measurement unit mismatch for channel {c}")
        time = number(float(record["time_s"]), "time_s")
        value = number(float(record["value"]), "value")
        key = (p, u, s, c, time)
        if key in seen:
            raise ValueError(f"duplicate observation: {key}")
        seen.add(key)
        rows.append(Row(p, u, s, c, time, value))
    if {r.preparation for r in rows} != set(expected):
        raise ValueError("input must contain every preparation declared for this split")
    return sorted(rows), digest(raw)


def weights(rows):
    """Equal preparations, units per preparation, and channels per unit."""
    preparations = defaultdict(set)
    units = defaultdict(set)
    counts = Counter()
    for row in rows:
        preparations[row.preparation].add(row.unit)
        units[row.preparation, row.unit].add(row.channel)
        counts[row.preparation, row.unit, row.channel] += 1
    return np.array([
        1 / (len(preparations) * len(preparations[r.preparation]) * len(units[r.preparation, r.unit])
             * counts[r.preparation, r.unit, r.channel]) for r in rows
    ])


def predict(plan, candidate, rows, fitted):
    readouts = {}
    for scenario in sorted({r.scenario for r in rows}):
        cfg = deepcopy(plan.models[candidate])
        cfg["parameters"] = {**cfg.get("parameters", {}), **plan.scenarios[scenario], **fitted}
        observations = simulate(cfg)["observations"]
        for channel in {r.channel for r in rows if r.scenario == scenario}:
            if channel not in observations or observations[channel]["unit"] != plan.channels[channel][0]:
                raise ValueError(f"{candidate}/{scenario}: missing observation or unit mismatch for {channel}")
            readout = observations[channel]
            readouts[scenario, channel] = dict(zip(readout["time_s"], readout["values"]))
    predicted = []
    for row in rows:
        readout = readouts[row.scenario, row.channel]
        if row.time_s not in readout:
            raise ValueError(f"{candidate}/{row.scenario}/{row.channel}: time {row.time_s} is not on the readout grid")
        predicted.append(readout[row.time_s])
    return np.array(predicted)


def metrics(rows, error, scales):
    weighted = weights(rows)
    normalized_mse = float(weighted @ (error / scales) ** 2)
    by_channel, signed = {}, {}
    for channel in sorted({r.channel for r in rows}):
        idx = [i for i, r in enumerate(rows) if r.channel == channel]
        channel_weights = weights([rows[i] for i in idx])
        by_channel[channel] = math.sqrt(float(channel_weights @ error[idx] ** 2))
        signed[channel] = float(channel_weights @ error[idx])
    return {"scaled_mse": normalized_mse, "scaled_rmse": math.sqrt(normalized_mse),
            "channel_rmse": by_channel, "channel_mean_signed_error": signed}


def fit(plan, data, version):
    rows, train_hash = load_rows(data, plan, plan.training)
    targets = np.array([r.value for r in rows])
    scales = np.array([plan.channels[r.channel][1] for r in rows])
    weighted = np.sqrt(weights(rows)) / scales
    candidates = {}
    for candidate in MODELS:
        specs = plan.specs[candidate]
        estimate = {}
        details = {"evaluations": 0, "jacobian_rank": 0, "active_bounds": []}
        if specs:
            initial, lower, upper = zip(*specs.values())

            def residual(parameters):
                fitted = dict(zip(specs, map(float, parameters)))
                return weighted * (predict(plan, candidate, rows, fitted) - targets)

            result = least_squares(residual, initial, bounds=(lower, upper), method="dogbox", max_nfev=300)
            if not result.success:
                raise ValueError(f"{candidate}: fit failed: {result.message}")
            estimate = dict(zip(specs, map(float, result.x)))
            details = {"evaluations": int(result.nfev), "jacobian_rank": int(np.linalg.matrix_rank(result.jac)),
                       "optimality": float(result.optimality), "active_bounds": [key for key, flag in zip(specs, result.active_mask) if flag]}
        error = predict(plan, candidate, rows, estimate) - targets
        candidates[candidate] = {"parameters": estimate, "diagnostics": details,
                                 "calibration": metrics(rows, error, scales)}
    return {"schema": SCHEMA, "software_version": version, "identity": plan.identity,
            "calibration_sha256": train_hash, "calibration_preparations": plan.training,
            "heldout_preparations_reserved": plan.heldout, "candidates": candidates}


def score(plan, data, fit_file, version):
    raw = fit_file.read_bytes()
    frozen = json.loads(raw)
    if (not isinstance(frozen, dict) or frozen.get("schema") != SCHEMA
            or frozen.get("software_version") != version or frozen.get("identity") != plan.identity
            or frozen.get("calibration_preparations") != plan.training
            or frozen.get("heldout_preparations_reserved") != plan.heldout):
        raise ValueError("frozen fit does not match the software, plan, candidate models and split")
    rows, target_hash = load_rows(data, plan, plan.heldout)
    if target_hash == frozen.get("calibration_sha256"):
        raise ValueError("held-out CSV equals the calibration CSV")
    targets = np.array([r.value for r in rows])
    scales = np.array([plan.channels[r.channel][1] for r in rows])
    totals, individual = {}, {}
    for candidate in MODELS:
        parameters = frozen["candidates"][candidate]["parameters"]
        if not isinstance(parameters, dict) or parameters.keys() != plan.specs[candidate].keys():
            raise ValueError(f"{candidate}: fitted parameter names do not match the plan")
        for key, value in parameters.items():
            _, lower, upper = plan.specs[candidate][key]
            if not lower <= number(value, key) <= upper:
                raise ValueError(f"{candidate}: fitted {key} is outside the bounds")
        error = predict(plan, candidate, rows, parameters) - targets
        totals[candidate] = metrics(rows, error, scales)
        individual[candidate] = {}
        for p in plan.heldout:
            idx = [i for i, row in enumerate(rows) if row.preparation == p]
            individual[candidate][p] = metrics([rows[i] for i in idx], error[idx], scales[idx])
    paired = {p: individual["modulated"][p]["scaled_mse"] - individual["stationary"][p]["scaled_mse"]
              for p in plan.heldout}
    leave_one_out = []
    for p in plan.heldout if len(plan.heldout) > 1 else ():
        rmse = {label: math.sqrt(float(np.mean([v["scaled_mse"] for key, v in individual[label].items() if key != p])))
                for label in MODELS}
        leave_one_out.append(rmse["modulated"] - rmse["stationary"])
    return {"schema": SCHEMA, "software_version": version, "fit_sha256": digest(raw),
            "heldout_sha256": target_hash, "heldout_preparations": plan.heldout,
            "channel_units": {channel: spec[0] for channel, spec in plan.channels.items()},
            "candidates": totals, "by_preparation": individual,
            "modulated_minus_stationary_scaled_rmse": totals["modulated"]["scaled_rmse"] - totals["stationary"]["scaled_rmse"],
            "modulated_minus_stationary_scaled_mse_by_preparation": paired,
            "leave_one_preparation_out_rmse_difference_range":
                [min(leave_one_out), max(leave_one_out)] if leave_one_out else None}
