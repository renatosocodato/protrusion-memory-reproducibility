"""Extract paired angular readouts from an example GUV, without fitting a model."""

import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
from PIL import Image, ImageSequence
from scipy.ndimage import gaussian_filter, map_coordinates

CHECKSUMS = (
    "d8ccb159075a0ee3361f0fb4f6bdd4ac1b85814563aee8182c1c594181dd9c4e",
    "53e9be6425736e7b8a3f84ff8d45c5d36b9c1e97629290eb855c075453ef1e52",
)
SAMPLE_S = 65.58921813964844
FRAMES = 28


def extract(membrane_path, actin_path, output, include_membrane=False):
    for path, expected in zip((membrane_path, actin_path), CHECKSUMS):
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError(f"source checksum mismatch for {path.name}")
    with Image.open(membrane_path) as membrane, Image.open(actin_path) as actin:
        first = gaussian_filter(np.asarray(membrane, dtype=float), 1)
        ring = first > np.quantile(first, 0.90)
        center = np.array(np.nonzero(ring)).mean(axis=1)
        radii = np.arange(10.0, 58.5, 0.5)
        theta = (np.arange(96) + 0.5) * (2 * np.pi / 96)
        sampled = center[:, None, None] + np.array([np.sin(theta), np.cos(theta)])[:, :, None] * radii
        last = None
        previous_ridge = None
        with output.open("x", newline="") as stream:
            writer = csv.writer(stream)
            columns = ["frame", "region", "radius_px", "radial_velocity_px_per_min", "actin_au"]
            if include_membrane:
                columns.append("membrane_au")
            writer.writerow(columns)
            mem_pages = ImageSequence.Iterator(membrane)
            act_pages = ImageSequence.Iterator(actin)
            for frame, (mem, act) in enumerate(zip(mem_pages, act_pages)):
                if frame >= FRAMES:
                    break
                profile = map_coordinates(gaussian_filter(np.asarray(mem, dtype=float), 1),
                                          sampled, order=1)
                if previous_ridge is not None:
                    allowed = np.abs(np.arange(len(radii)) - previous_ridge[:, None]) <= 12
                    profile = np.where(allowed, profile, -np.inf)
                ridge = np.argmax(profile, axis=1)
                valid = (ridge >= 4) & (ridge < len(radii) - 2)
                if previous_ridge is not None:
                    valid &= np.abs(ridge - previous_ridge) < 12
                previous_ridge = np.where(valid, ridge, previous_ridge) if previous_ridge is not None else ridge
                # Measure the cytoplasmic side only. The membrane channel defines the ridge.
                band = np.asarray(act, dtype=float)
                act_profile = map_coordinates(band, sampled, order=1)
                ray_indices = np.maximum(0, ridge[:, None] - np.arange(2, 5))
                local = np.take_along_axis(act_profile, ray_indices, axis=1).mean(axis=1)
                # Keep like-for-like rays for the sector averages and differences.
                good_sector = valid.reshape(8, 12).all(axis=1)
                intensity = local.reshape(8, 12).mean(axis=1)
                position = np.where(good_sector, radii[ridge].reshape(8, 12).mean(axis=1), np.nan)
                velocity = np.full(8, np.nan) if last is None else (position - last) * 60 / SAMPLE_S
                if include_membrane:
                    ridge_signal = profile[np.arange(len(ridge)), ridge].reshape(8, 12).mean(axis=1)
                for region in range(8):
                    row = [frame, region, position[region] if good_sector[region] else "",
                           velocity[region] if np.isfinite(velocity[region]) else "",
                           intensity[region] if good_sector[region] else ""]
                    if include_membrane:
                        row.append(ridge_signal[region] if good_sector[region] else "")
                    writer.writerow(row)
                last = position


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("membrane", type=Path)
    parser.add_argument("actin", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--include-membrane", action="store_true", help="also export local mCherry")
    args = parser.parse_args()
    extract(args.membrane, args.actin, args.output, args.include_membrane)
