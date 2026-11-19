"""Separate image-plane centroid motion and angular shape in one public vesicle."""

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np
from PIL import Image, ImageSequence
from scipy.ndimage import label, map_coordinates, shift

MASK_HASH = "df0a0cccd13d721bee54fb237dd9047785bf70257910cc5850f07fee131ec767"
ACTIN_HASH = "53e9be6425736e7b8a3f84ff8d45c5d36b9c1e97629290eb855c075453ef1e52"
PIXELS_PER_UM = 2.276379
SAMPLE_S = 65.58921813964844
MODES = (1, 2, 3, 4)
THETA = (np.arange(96) + .5) * (2*np.pi/96)
DIRECTIONS = np.array([np.sin(THETA), np.cos(THETA)])


def radial(mask, center):
    if not map_coordinates(mask, np.asarray(center)[:, None], order=1)[0] >= .5:
        return None
    radii = np.arange(0., np.linalg.norm(mask.shape)+1., .25)
    coordinates = center[:, None, None] + DIRECTIONS[:, :, None]*radii
    samples = map_coordinates(mask, coordinates, order=1, mode="constant")
    inside = samples >= .5
    transitions = np.diff(inside.astype(np.int8), axis=1)
    if np.any((transitions != 0).sum(axis=1) != 1):
        return None
    indices = np.argmax(transitions == -1, axis=1)
    first = samples[np.arange(len(THETA)), indices]
    last = samples[np.arange(len(THETA)), indices+1]
    return radii[indices] + .25*(first-.5)/(first-last)


def projections(values):
    return {str(n): [float(2*np.mean(values*np.cos(n*THETA))),
                     float(2*np.mean(values*np.sin(n*THETA)))] for n in MODES}


def describe(mask, actin, reference=None):
    occupied = mask >= .5
    if label(occupied)[1] != 1:
        raise ValueError("mask must contain one connected object")
    if np.any(occupied[[0,-1]]) or np.any(occupied[:, [0,-1]]):
        raise ValueError("mask touches image border")
    center = np.asarray(np.nonzero(occupied)).mean(axis=1)
    centered = radial(mask, center)
    laboratory = radial(mask, center if reference is None else reference)
    record = {
        "centroid_xy_um": (center[::-1]/PIXELS_PER_UM).tolist(),
        "area_equivalent_radius_um": float(np.sqrt(occupied.sum()/np.pi)/PIXELS_PER_UM),
        "lab_radius_modes_um": None if laboratory is None else projections(laboratory/PIXELS_PER_UM),
        "centered_radius_modes_um": None,
        "cortical_actin_modes_relative": None,
    }
    if centered is not None:
        band = center[:,None,None] + DIRECTIONS[:,:,None]*(centered[:,None]-np.array([1.,1.5,2.]))
        intensity = map_coordinates(actin, band, order=1).mean(axis=1)
        record["centered_radius_modes_um"] = projections(centered/PIXELS_PER_UM)
        if intensity.mean() > 0:
            record["cortical_actin_modes_relative"] = projections(intensity/intensity.mean()-1)
    return center, record


def geometry_controls(mask, center, original):
    moved = shift(mask, (0.75, -1.25), order=1)
    translated, record = describe(moved, np.ones_like(moved))
    known_shift_error = float(np.linalg.norm(translated-center-[0.75,-1.25])/PIXELS_PER_UM)
    if record["centered_radius_modes_um"] is None:
        # A raster staircase can produce multiple ray crossings.
        n2_error = None
    else:
        n2_error = float(np.linalg.norm(np.asarray(record["centered_radius_modes_um"]["2"])
                                        - original["centered_radius_modes_um"]["2"]))
    yy, xx = np.indices((121,121), dtype=float)
    r = np.hypot(yy-60, xx-60)
    phi = np.arctan2(yy-60, xx-60)
    bumped = (r <= 30+6*np.exp(-phi**2/.16)).astype(float)
    plus, _ = describe(bumped, np.ones_like(bumped))
    return {"subpixel_shift": {"known_shift_yx_px": [.75,-1.25], "centroid_error_um": known_shift_error,
                               "centered_mode2_change_um": n2_error},
            "stationary_center_bump": {"projected_centroid_shift_um": float(
                np.linalg.norm(plus-[60,60])/PIXELS_PER_UM)}}


def analyze(mask_path, actin_path):
    for path, expected in ((mask_path, MASK_HASH), (actin_path, ACTIN_HASH)):
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"source checksum mismatch for {path.name}")
    records = []
    with Image.open(mask_path) as mask, Image.open(actin_path) as actin:
        if mask.n_frames != 28 or actin.n_frames < mask.n_frames or mask.size != actin.size:
            raise ValueError("unexpected G01 image geometry")
        reference = None
        for i, (m, a) in enumerate(zip(ImageSequence.Iterator(mask), ImageSequence.Iterator(actin))):
            binary = (np.asarray(m) > 0).astype(float)
            center, frame = describe(binary, np.asarray(a, dtype=float), reference)
            if reference is None:
                reference = center
                checks = geometry_controls(binary, center, frame)
            records.append({"frame": i, "time_s": i*SAMPLE_S, **frame})
    centered = [r for r in records if r["centered_radius_modes_um"] is not None]
    both = [r for r in centered if r["lab_radius_modes_um"] is not None]
    shifted = np.array([r["centroid_xy_um"] for r in records])
    dist = np.linalg.norm(shifted-shifted[0], axis=1)
    baseline = [r["centered_radius_modes_um"] for r in records[:6]
                if r["centered_radius_modes_um"] is not None]
    if len(baseline) < 4:
        raise ValueError("too few valid centered baseline frames")
    summaries = {}
    for mode in ("2", "3", "4"):
        norm = lambda modes: float(np.linalg.norm(modes[mode]))
        offsets = np.asarray([b[mode] for b in baseline])
        baseline_mean = offsets.mean(axis=0)
        summaries[mode] = {
            "baseline_component_sd_um": np.std(offsets, axis=0, ddof=1).tolist(),
            "max_centered_amplitude_um": max(norm(r["centered_radius_modes_um"]) for r in centered),
            "max_distance_from_baseline_mode_um": max(float(np.linalg.norm(
                np.asarray(r["centered_radius_modes_um"][mode])-baseline_mean)) for r in centered),
        }
    return {"data_kind": "one G01 vesicle, membrane-defined contour", "vesicles": 1,
            "sample_s": SAMPLE_S, "pixels_per_um": PIXELS_PER_UM, "frames": records,
            "geometry_controls": checks,
            "summary": {"valid_centered_frames": len(centered), "common_radial_frames": len(both),
                        "valid_baseline_frames": len(baseline),
                        "max_centroid_shift_um": float(dist.max()),
                        "area_radius_first_last_um": [records[0]["area_equivalent_radius_um"],
                                                     records[-1]["area_equivalent_radius_um"]],
                        "median_lab_mode1_um": float(np.median([np.linalg.norm(r["lab_radius_modes_um"]["1"]) for r in both])),
                        "median_centered_mode1_um": float(np.median([np.linalg.norm(r["centered_radius_modes_um"]["1"]) for r in both])),
                        "higher_modes": summaries}}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mask", type=Path)
    parser.add_argument("actin", type=Path)
    args = parser.parse_args()
    print(json.dumps(analyze(args.mask, args.actin), indent=2, allow_nan=False))
