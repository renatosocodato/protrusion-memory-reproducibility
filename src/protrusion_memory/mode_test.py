"""Calibrate angular readouts and test a reserved, signed intervention on one mode."""

import csv
from dataclasses import dataclass
import io
import json

import numpy as np

from .calibration import digest, load_table, names
from .portable import number, table
from .spatial import integer

SCHEMA = "mode-test-v1"
MECHANICAL_SCHEMA = "mode-test-v2"
MECHANICAL_CRITERIA = ("blocked_actin_max", "blocked_actin_gain", "transfer_margin", "control_shape_error")
CHANNELS = ("actin", "shape", "drive")
ARMS = ("active_plus", "delivery_plus", "active_minus", "delivery_minus")
GAIN_COLUMNS = {"preparation", "mode", "channel", "known_cos", "known_sin", "read_cos", "read_sin"}
RESPONSE_COLUMNS = {"preparation", "mode", "arm", *(f"{c}_{a}" for c in CHANNELS for a in ("cos", "sin"))}


@dataclass(frozen=True)
class Plan:
    digest: str
    mode: int
    axis: np.ndarray
    window: tuple
    calibration: list
    reserved: list
    actin_min: float
    shape_min: float
    drive_max: float
    drive_gain: float
    delivery_bound: float
    drive_error: float
    seed: int
    bootstrap: int
    schema: str
    mechanical: list
    blocked_actin_max: float
    blocked_actin_gain: float
    transfer_margin: float
    control_shape_error: float

    @classmethod
    def load(cls, path):
        raw, fingerprint = load_table(path)
        table(raw, "mode plan", fields={"schema", "mode", "axis", "window_s", "split", "criteria", "resampling"})
        schema = raw.get("schema")
        if schema not in (SCHEMA, MECHANICAL_SCHEMA):
            raise ValueError(f"expected {SCHEMA} or {MECHANICAL_SCHEMA}")
        measured = schema == MECHANICAL_SCHEMA
        angle = raw["axis"]
        if not isinstance(angle, list) or len(angle) != 2:
            raise ValueError("axis must give two cosine/sine weights")
        axis = np.asarray([number(x, "axis") for x in angle])
        if not np.isclose(np.linalg.norm(axis), 1, atol=1e-10, rtol=0):
            raise ValueError("axis must have unit norm")
        window = raw["window_s"]
        if not isinstance(window, list) or len(window) != 2:
            raise ValueError("window_s must give two finite endpoints")
        window = tuple(number(x, "window_s") for x in window)
        if not 0 <= window[0] < window[1]:
            raise ValueError("require 0 <= window start < stop")
        split_names = ("calibration", "reserved", "mechanical") if measured else ("calibration", "reserved")
        split = table(raw["split"], "split", fields=set(split_names))
        groups = [names(split.get(k), k) for k in split_names]
        if min(map(len, groups)) < 3 or len(set().union(*groups)) != sum(map(len, groups)):
            raise ValueError("require at least three disjoint preparations per split")
        calibration, reserved = groups[:2]
        keys = ("actin_min", "shape_min", "drive_max", "drive_gain")
        extra = MECHANICAL_CRITERIA if measured else ("delivery_bound",)
        criteria = table(raw["criteria"], "criteria", fields=set(keys + extra) | {"drive_error"})
        values = [number(criteria.get(k), k) for k in keys]
        mechanical_values = [number(criteria.get(k), k) for k in extra] if measured else [0] * 4
        values.append(0 if measured else number(criteria.get("delivery_bound"), "delivery_bound"))
        values.append(number(criteria.get("drive_error", 0), "drive_error"))
        if min(values + mechanical_values) < 0 or min(values[:2]) <= 0:
            raise ValueError("physical limits must be nonnegative and meaningful effects positive")
        rs = table(raw["resampling"], "resampling", fields={"seed", "draws"})
        draws = integer(rs["draws"], "draws", 200, 20000)
        if draws * max(map(len, groups)) > 2_000_000:
            raise ValueError("mode test exceeds the resampling budget")
        return cls(fingerprint, integer(raw["mode"], "mode", 2, 64), axis, window, calibration, reserved,
                   *values, integer(rs["seed"], "seed", 0, 2**63-1), draws,
                   schema, groups[2] if measured else [], *mechanical_values)


def rows(path, columns, preparations, mode):
    data = path.read_bytes()
    reader = csv.DictReader(io.StringIO(data.decode("utf-8-sig"), newline=""))
    if set(reader.fieldnames or ()) != columns or len(reader.fieldnames) != len(columns):
        raise ValueError("wrong CSV columns; see the mode-test specification")
    records = []
    for row in reader:
        if row.keys() != columns or any(v is None for v in row.values()):
            raise ValueError(f"invalid CSV row {reader.line_num}")
        if row["preparation"] not in preparations or int(row["mode"]) != mode:
            raise ValueError("CSV contains undeclared preparation or mode")
        records.append(row)
    if {r["preparation"] for r in records} != set(preparations):
        raise ValueError("CSV is missing a declared preparation")
    return records, digest(data)


def coordinates(row, channel):
    return [number(float(row[f"{channel}_{axis}"]), f"{channel}_{axis}") for axis in ("cos", "sin")]


def angular_matrices(gram, cross):
    gains = np.swapaxes(np.linalg.solve(gram, cross), -2, -1)
    if not np.isfinite(gains).all() or np.any(np.linalg.cond(gains) > 100):
        raise ValueError("angular gains are ill-conditioned")
    return gains


def calibrate(plan, path, version):
    records, fingerprint = rows(path, GAIN_COLUMNS, plan.calibration, plan.mode)
    groups = {(p, c): [] for p in plan.calibration for c in CHANNELS}
    for row in records:
        key = row["preparation"], row["channel"]
        if key not in groups:
            raise ValueError("unknown angular calibration channel")
        groups[key].append((coordinates(row, "known"), coordinates(row, "read")))
    gram, cross = [], []
    for p in plan.calibration:
        design, response = [], []
        for c in CHANNELS:
            values = np.asarray(groups[p, c])
            if len(values) < 3 or np.linalg.matrix_rank(values[:, 0, :]) != 2:
                raise ValueError(f"{p}/{c}: calibration must span both angular axes")
            x, y = values[:, 0, :], values[:, 1, :]
            design.append(x.T @ x / len(x))
            response.append(x.T @ y / len(y))
        gram.append(design)
        cross.append(response)
    gram, cross = np.array(gram), np.array(cross)
    angular_matrices(gram, cross)
    gains = angular_matrices(gram.mean(axis=0), cross.mean(axis=0))
    return {"schema": plan.schema, "software_version": version, "plan_sha256": plan.digest,
            "calibration_sha256": fingerprint, "preparations": plan.calibration,
            "channel_order": CHANNELS, "gain_matrices": dict(zip(CHANNELS, gains.tolist())),
            "gram": gram.tolist(), "cross": cross.tolist()}


def paired_responses(plan, path, preparations=None):
    preparations = plan.reserved if preparations is None else preparations
    records, fingerprint = rows(path, RESPONSE_COLUMNS, preparations, plan.mode)
    grouped = {}
    for row in records:
        key = row["preparation"], row["arm"]
        if key[1] not in ARMS or key in grouped:
            raise ValueError("unknown or duplicate intervention arm")
        grouped[key] = [coordinates(row, c) for c in CHANNELS]
    if any((p, a) not in grouped for p in preparations for a in ARMS):
        raise ValueError("each preparation needs all four intervention arms")
    raw = np.asarray([[grouped[p, a] for a in ARMS] for p in preparations])
    return raw[:, ::2]-raw[:, 1::2], fingerprint


def decode(raw, gain):
    # Axes: bootstrap draw, preparation, sign, channel, cosine/sine.
    return np.linalg.solve(gain[:, None, None, :, :, :], raw[..., None])[..., 0]


def effects(raw, gain, plan):
    latent = decode(raw, gain)
    projected = latent @ plan.axis
    actin = np.stack((projected[:, :, 0, 0], -projected[:, :, 1, 0]), axis=-1)
    shape = .5*(projected[:, :, 0, 1]-projected[:, :, 1, 1])
    residual_drive = np.linalg.norm(latent[:, :, :, 2, :], axis=-1)
    perpendicular = np.array([-plan.axis[1], plan.axis[0]])
    off_axis = .5*(latent[:, :, 0, 1, :]-latent[:, :, 1, 1, :]) @ perpendicular
    return actin, shape, residual_drive, off_axis


def angular_draws(plan, frozen_path, version):
    frozen_raw = frozen_path.read_bytes()
    frozen = json.loads(frozen_raw)
    if (not isinstance(frozen, dict) or frozen.get("schema") != plan.schema
            or frozen.get("software_version") != version or frozen.get("plan_sha256") != plan.digest
            or frozen.get("preparations") != plan.calibration or frozen.get("channel_order") != list(CHANNELS)):
        raise ValueError("frozen angular calibration does not match this plan and version")
    gram, cross = (np.asarray(frozen[k], dtype=float) for k in ("gram", "cross"))
    if any(m.shape != (len(plan.calibration), 3, 2, 2) or not np.isfinite(m).all() for m in (gram, cross)):
        raise ValueError("invalid frozen angular calibration matrices")
    a_rng, reserve_rng, control_rng = map(np.random.default_rng, np.random.SeedSequence(plan.seed).spawn(3))
    a = a_rng.integers(0, len(gram), size=(plan.bootstrap, len(gram)))
    gains = angular_matrices(gram[a].mean(axis=1), cross[a].mean(axis=1))
    pooled = angular_matrices(gram.mean(axis=0), cross.mean(axis=0))[None, ...]
    return digest(frozen_raw), gains, pooled, reserve_rng, control_rng


def control_summary(plan, raw, gains, pooled, rng):
    def components(latent):
        shape = .5 * ((latent[:, :, 0, 1] - latent[:, :, 1, 1]) @ plan.axis)
        actin = np.linalg.norm(latent[:, :, :, 0], axis=-1)
        drive = np.linalg.norm(latent[:, :, :, 2], axis=-1)
        even = .5 * ((latent[:, :, 0, 1] + latent[:, :, 1, 1]) @ plan.axis)
        return shape, actin, drive, even

    def allowance(shape, actin, drive, _even):
        return (np.abs(shape.mean(axis=1)) + plan.blocked_actin_gain * actin.mean(axis=(1, 2))
                + plan.drive_gain * (drive.mean(axis=(1, 2)) + plan.drive_error)
                + plan.control_shape_error + plan.transfer_margin)

    observed = components(decode(raw[None, ...], pooled))
    # Keep every observed preparation in the qualification guard while varying gain uncertainty.
    _, actin_guard, drive_guard, _ = components(decode(raw[None, ...], gains))
    actin_upper = float(np.quantile(actin_guard.max(axis=(1, 2)), .975))
    drive_upper = float(np.quantile(drive_guard.max(axis=(1, 2)), .975)) + plan.drive_error
    indices = rng.integers(0, len(raw), size=(plan.bootstrap, len(raw)))
    bounds = allowance(*components(decode(raw[indices], gains)))
    summary = {"status": "qualified" if actin_upper <= plan.blocked_actin_max and drive_upper <= plan.drive_max
               else "control_unqualified",
               "control_actin_upper": actin_upper, "control_drive_upper": drive_upper,
               "control_odd_mean": float(observed[0].mean()),
               "control_even_mean": float(observed[3].mean()),
               "mechanical_allowance": float(allowance(*observed)[0]),
               "mechanical_allowance_upper": float(np.quantile(bounds, .975))}
    return summary, bounds


def calibrate_delivery(plan, path, gain_path, version):
    if plan.schema != MECHANICAL_SCHEMA:
        raise ValueError("mechanical calibration requires a mode-test-v2 plan")
    gain_sha, gains, pooled, _, rng = angular_draws(plan, gain_path, version)
    raw, fingerprint = paired_responses(plan, path, plan.mechanical)
    summary, _ = control_summary(plan, raw, gains, pooled, rng)
    return {"schema": "mode-delivery-calibration-v1", "software_version": version, "plan_sha256": plan.digest,
            "angular_calibration_sha256": gain_sha, "mechanical_csv_sha256": fingerprint,
            "mechanical_preparations": plan.mechanical, "raw_contrasts": raw.tolist(), "summary": summary}


def score(plan, path, frozen_path, version, delivery_path=None):
    measured = plan.schema == MECHANICAL_SCHEMA
    if measured != (delivery_path is not None):
        raise ValueError("mode-test-v2 requires a separate mechanical calibration; mode-test-v1 uses its fixed bound")
    gain_sha, gains, pooled, r_rng, control_rng = angular_draws(plan, frozen_path, version)
    delivery_bounds = plan.delivery_bound
    mechanical = None
    if measured:
        frozen_delivery = delivery_path.read_bytes()
        control = json.loads(frozen_delivery)
        if (not isinstance(control, dict) or control.get("schema") != "mode-delivery-calibration-v1"
                or control.get("software_version") != version or control.get("plan_sha256") != plan.digest
                or control.get("angular_calibration_sha256") != gain_sha
                or control.get("mechanical_preparations") != plan.mechanical):
            raise ValueError("mechanical calibration does not match this plan, version and angular calibration")
        contrasts = np.asarray(control.get("raw_contrasts"), dtype=float)
        if contrasts.shape != (len(plan.mechanical), 2, 3, 2) or not np.isfinite(contrasts).all():
            raise ValueError("invalid mechanical calibration contrasts")
        mechanical, delivery_bounds = control_summary(plan, contrasts, gains, pooled, control_rng)
        mechanical.update(calibration_sha256=digest(frozen_delivery), preparations=plan.mechanical)
    raw, fingerprint = paired_responses(plan, path)
    r = r_rng.integers(0, len(raw), size=(plan.bootstrap, len(raw)))
    actin, shape, drive, _ = effects(raw[r], gains, plan)
    first_stage = actin.mean(axis=1).min(axis=1)
    shape_mean = shape.mean(axis=1)
    nuisance = plan.drive_gain*(drive.mean(axis=(1,2))+plan.drive_error) + delivery_bounds
    robust = [float(np.quantile(shape_mean-nuisance, .025)), float(np.quantile(shape_mean+nuisance, .975))]
    if mechanical is not None and mechanical["status"] == "qualified":
        mechanical["transfer_margin_limit_for_signed_effect"] = float(max(0, plan.transfer_margin
            + max(robust[0], -robust[1]) - plan.shape_min))
    actin_interval = np.quantile(first_stage, [.025, .975]).tolist()
    # Use simultaneous drive max across both signs and across the observed reserve preparations.
    actin_point, shape_point, drive_point, off_axis = (x[0] for x in effects(raw[None, ...], pooled, plan))
    drive_guard = effects(np.broadcast_to(raw, (len(gains), *raw.shape)), gains, plan)[2]
    drive_upper = float(np.quantile(drive_guard.max(axis=(1, 2)), .975))+plan.drive_error
    if mechanical is not None and mechanical["status"] != "qualified":
        decision = "mechanical_control_unqualified"
    elif drive_upper > plan.drive_max:
        decision = "drive_not_matched"
    elif actin_interval[0] < plan.actin_min:
        decision = "actin_unresolved"
    elif robust[0] > plan.shape_min or robust[1] < -plan.shape_min:
        decision = "deformation_beyond_leakage"
    elif -plan.shape_min <= robust[0] <= robust[1] <= plan.shape_min:
        decision = "deformation_bounded_small"
    else:
        decision = "deformation_unresolved"
    limits = ("actin_min", "shape_min", "drive_max", "drive_gain", "drive_error")
    limits += MECHANICAL_CRITERIA if measured else ("delivery_bound",)
    return {"schema": plan.schema, "software_version": version, "mode": plan.mode,
            "window_s": plan.window, "axis": plan.axis.tolist(), "decision": decision,
            "limits": {name: getattr(plan, name) for name in limits},
            **({"mechanical_calibration": mechanical} if measured else {}),
            "calibration_sha256": gain_sha, "reserved_sha256": fingerprint,
            "reserved_preparations": plan.reserved, "drive_upper": drive_upper,
            "actin_both_signs_interval": actin_interval,
            "shape_odd_mean": float(shape_point.mean()),
            "shape_odd_interval": np.quantile(shape_mean, [.025, .975]).tolist(),
            "shape_interval_with_leakage": None if decision == "mechanical_control_unqualified" else robust,
            "mean_leakage_allowance": (
                None if decision == "mechanical_control_unqualified" else
                float(plan.drive_gain*(drive_point.mean()+plan.drive_error)
                      + (mechanical["mechanical_allowance"] if measured else plan.delivery_bound))),
            "by_preparation": {p: {"actin_signed_by_arm": actin_point[i].tolist(),
                                   "shape_odd": float(shape_point[i]),
                                   "shape_perpendicular_odd": float(off_axis[i]),
                                   "residual_drive_norms": drive_point[i].tolist()}
                               for i, p in enumerate(plan.reserved)}}
