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

    @classmethod
    def load(cls, path):
        raw, fingerprint = load_table(path)
        table(raw, "mode plan", fields={"schema", "mode", "axis", "window_s", "split", "criteria", "resampling"})
        if raw.get("schema") != SCHEMA:
            raise ValueError(f"expected {SCHEMA}")
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
        split = table(raw["split"], "split", fields={"calibration", "reserved"})
        calibration, reserved = (names(split[k], k) for k in ("calibration", "reserved"))
        if min(len(calibration), len(reserved)) < 3 or set(calibration) & set(reserved):
            raise ValueError("require at least three disjoint preparations per split")
        keys = ("actin_min", "shape_min", "drive_max", "drive_gain", "delivery_bound")
        criteria = table(raw["criteria"], "criteria", fields=set(keys)|{"drive_error"})
        values = [number(criteria[k], k) for k in keys]
        values.append(number(criteria.get("drive_error", 0), "drive_error"))
        if min(values) < 0 or min(values[:2]) <= 0:
            raise ValueError("physical limits must be nonnegative and meaningful effects positive")
        rs = table(raw["resampling"], "resampling", fields={"seed", "draws"})
        draws = integer(rs["draws"], "draws", 200, 20000)
        if draws * max(len(calibration), len(reserved)) > 2_000_000:
            raise ValueError("mode test exceeds the resampling budget")
        return cls(fingerprint, integer(raw["mode"], "mode", 2, 64), axis, window, calibration, reserved,
                   *values, integer(rs["seed"], "seed", 0, 2**63-1), draws)


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
    return {"schema": SCHEMA, "software_version": version, "plan_sha256": plan.digest,
            "calibration_sha256": fingerprint, "preparations": plan.calibration,
            "channel_order": CHANNELS, "gain_matrices": dict(zip(CHANNELS, gains.tolist())),
            "gram": gram.tolist(), "cross": cross.tolist()}


def paired_responses(plan, path):
    records, fingerprint = rows(path, RESPONSE_COLUMNS, plan.reserved, plan.mode)
    grouped = {}
    for row in records:
        key = row["preparation"], row["arm"]
        if key[1] not in ARMS or key in grouped:
            raise ValueError("unknown or duplicate intervention arm")
        grouped[key] = [coordinates(row, c) for c in CHANNELS]
    if any((p, a) not in grouped for p in plan.reserved for a in ARMS):
        raise ValueError("each reserved preparation needs all four intervention arms")
    raw = np.asarray([[grouped[p, a] for a in ARMS] for p in plan.reserved])
    return raw[:, ::2]-raw[:, 1::2], fingerprint


def effects(raw, gain, plan):
    # Axes: bootstrap draw, preparation, sign, channel, cosine/sine.
    latent = np.linalg.solve(gain[:, None, None, :, :, :], raw[..., None])[..., 0]
    projected = latent @ plan.axis
    actin = np.stack((projected[:, :, 0, 0], -projected[:, :, 1, 0]), axis=-1)
    shape = .5*(projected[:, :, 0, 1]-projected[:, :, 1, 1])
    residual_drive = np.linalg.norm(latent[:, :, :, 2, :], axis=-1)
    perpendicular = np.array([-plan.axis[1], plan.axis[0]])
    off_axis = .5*(latent[:, :, 0, 1, :]-latent[:, :, 1, 1, :]) @ perpendicular
    return actin, shape, residual_drive, off_axis


def score(plan, path, frozen_path, version):
    frozen_raw = frozen_path.read_bytes()
    frozen = json.loads(frozen_raw)
    if (not isinstance(frozen, dict) or frozen.get("schema") != SCHEMA
            or frozen.get("software_version") != version or frozen.get("plan_sha256") != plan.digest
            or frozen.get("preparations") != plan.calibration or frozen.get("channel_order") != list(CHANNELS)):
        raise ValueError("frozen angular calibration does not match this plan and version")
    raw, fingerprint = paired_responses(plan, path)
    gram, cross = (np.asarray(frozen[k], dtype=float) for k in ("gram", "cross"))
    if any(m.shape != (len(plan.calibration), 3, 2, 2) or not np.isfinite(m).all() for m in (gram, cross)):
        raise ValueError("invalid frozen angular calibration matrices")
    root = np.random.SeedSequence(plan.seed)
    a_rng, r_rng = map(np.random.default_rng, root.spawn(2))
    a = a_rng.integers(0, len(gram), size=(plan.bootstrap,len(gram)))
    r = r_rng.integers(0, len(raw), size=(plan.bootstrap,len(raw)))
    gains = angular_matrices(gram[a].mean(axis=1), cross[a].mean(axis=1))
    actin, shape, drive, _ = effects(raw[r], gains, plan)
    first_stage = actin.mean(axis=1).min(axis=1)
    shape_mean = shape.mean(axis=1)
    nuisance = plan.drive_gain*(drive.mean(axis=(1,2))+plan.drive_error) + plan.delivery_bound
    robust = [float(np.quantile(shape_mean-nuisance, .025)), float(np.quantile(shape_mean+nuisance, .975))]
    actin_interval = np.quantile(first_stage, [.025, .975]).tolist()
    # Use simultaneous drive max across both signs and across the observed reserve preparations.
    pooled = angular_matrices(gram.mean(axis=0), cross.mean(axis=0))[None, ...]
    actin_point, shape_point, drive_point, off_axis = (x[0] for x in effects(raw[None, ...], pooled, plan))
    drive_guard = effects(np.broadcast_to(raw, (len(gains), *raw.shape)), gains, plan)[2]
    drive_upper = float(np.quantile(drive_guard.max(axis=(1, 2)), .975))+plan.drive_error
    if drive_upper > plan.drive_max:
        decision = "drive_not_matched"
    elif actin_interval[0] < plan.actin_min:
        decision = "actin_unresolved"
    elif robust[0] > plan.shape_min or robust[1] < -plan.shape_min:
        decision = "deformation_beyond_leakage"
    elif -plan.shape_min <= robust[0] <= robust[1] <= plan.shape_min:
        decision = "deformation_bounded_small"
    else:
        decision = "deformation_unresolved"
    return {"schema": SCHEMA, "software_version": version, "mode": plan.mode,
            "window_s": plan.window, "axis": plan.axis.tolist(), "decision": decision,
            "limits": {name: getattr(plan, name) for name in ("actin_min", "shape_min", "drive_max",
                                                               "drive_gain", "delivery_bound", "drive_error")},
            "calibration_sha256": digest(frozen_raw), "reserved_sha256": fingerprint,
            "reserved_preparations": plan.reserved, "drive_upper": drive_upper,
            "actin_both_signs_interval": actin_interval,
            "shape_odd_mean": float(shape_point.mean()),
            "shape_odd_interval": np.quantile(shape_mean, [.025, .975]).tolist(),
            "shape_interval_with_leakage": robust,
            "mean_leakage_allowance": float(plan.drive_gain*(drive_point.mean()+plan.drive_error)+plan.delivery_bound),
            "by_preparation": {p: {"actin_signed_by_arm": actin_point[i].tolist(),
                                   "shape_odd": float(shape_point[i]),
                                   "shape_perpendicular_odd": float(off_axis[i]),
                                   "residual_drive_norms": drive_point[i].tolist()}
                               for i, p in enumerate(plan.reserved)}}
