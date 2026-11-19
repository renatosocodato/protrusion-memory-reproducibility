"""Synthetic calibrated discrimination at cadences from the acquisition studies."""

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import ndtr

from challenge import basis, pulse
from guv_g01 import SAMPLE_S

PIXELS_PER_UM = 2.276379
SCALE_MULTIPLIERS = (.5, 1., 2.)
PHASES = (0., .25, .5, .75, .99)
LABELS = ("unresolved", "mediated", "direct")
N_CAL = 6000
N_TEST = 4000


def prefix_scale(path):
    sector = defaultdict(list)
    with path.open() as stream:
        for row in csv.DictReader(stream):
            if int(row["frame"]) < 6 and row["radius_px"]:
                sector[row["region"]].append(float(row["radius_px"]) / PIXELS_PER_UM)
    sd = {k: float(np.std(values, ddof=1)) for k, values in sorted(sector.items()) if len(values) == 6}
    if len(sd) != 8:
        raise ValueError("the baseline requires six valid frames in each of the eight sectors")
    return {"sectors": len(sd), "vesicles": 1, "baseline_frames": 6, "sector_sd_um": sd,
            "median_sector_sd_um": float(np.median(list(sd.values())))}


def covariance(times, scale, kind, exposure_s=.7, tau=20):
    if kind == "white":
        return scale**2 * np.eye(len(times))
    # Covariance of a stationary OU process averaged in disjoint equal windows.
    lag = abs(times[:, None] - times)
    fraction = -math.expm1(-exposure_s / tau)
    var = 2*tau / exposure_s**2 * (exposure_s - tau*fraction)
    result = ((tau*fraction/exposure_s)**2 *
              np.exp(-np.maximum(0, lag-exposure_s)/tau)) / var
    np.fill_diagonal(result, 1)
    return scale**2 * result


def interval(wins, n):
    p = wins/n
    z = 1.959963984540054
    middle = (p + z*z/(2*n))/(1+z*z/n)
    margin = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/(1+z*z/n)
    return [max(0, middle-margin), min(1, middle+margin)]


def classify(template, cov, streams):
    lower = np.linalg.cholesky(cov)
    whitened = solve_triangular(lower, template, lower=True)
    offset = solve_triangular(lower, np.ones(len(template)), lower=True)
    offset /= np.linalg.norm(offset)
    centered = whitened - np.outer(offset, offset @ whitened)
    norms = np.linalg.norm(centered, axis=0)
    rays = centered / np.where(norms > 1e-12, norms, np.inf)
    # Two projections of the full Gaussian observation retain all information
    # needed by the two likelihood scores. This avoids generating long traces.
    correlation = rays.T @ rays
    eig, vectors = np.linalg.eigh(correlation)
    transform = vectors * np.sqrt(np.maximum(eig, 0))
    cal = np.random.default_rng(streams[0]).normal(size=(N_CAL, 2)) @ transform.T
    cutoff = float(np.quantile((.5*np.maximum(cal, 0)**2).max(axis=1), .95, method="higher"))
    means = np.column_stack((np.zeros(2), rays.T @ centered))
    counts = {}
    for i, label in enumerate(("null", "mediated", "direct")):
        noise = np.random.default_rng(streams[i+1]).normal(size=(N_TEST, 2)) @ transform.T
        scores = .5*np.maximum(noise+means[:, i], 0)**2
        best = scores.argmax(axis=1)
        detectable = scores.max(axis=1) > cutoff
        decided = detectable & (abs(scores[:, 0] - scores[:, 1]) >= 2)
        tags = np.where(decided, best+1, 0)
        counts[label] = dict(zip(LABELS, np.bincount(tags, minlength=3).tolist()))
        counts[label]["detected_before_margin"] = int(detectable.sum())
    dprime = float(np.linalg.norm(centered[:, 1]-centered[:, 0]))
    return {"dprime_known_amplitudes": dprime, "ideal_balanced_accuracy": float(ndtr(dprime/2)),
            "confusion": counts, "null_detection_interval": interval(counts["null"]["detected_before_margin"], N_TEST),
            "correct_interval": {k: interval(counts[k][k], N_TEST) for k in LABELS[1:]}}


def study(path):
    empirical = prefix_scale(path)
    scale = empirical["median_sector_sd_um"]
    designs = [(name, step, phase) for name, step in (("1s", 1.), ("G01", SAMPLE_S), ("H1299_clock", 300.))
               for phase in PHASES]
    designs.append(("G01_near_frame", SAMPLE_S, None))
    generator = np.random.SeedSequence(980218)
    results = []
    for clock, sampling, phase in designs:
        if phase is None:
            onset = 130.
        else:
            onset = 130 + phase if clock == "1s" else (2+phase)*sampling
        config = pulse(sampling, .7, .6, onset)
        if clock == "H1299_clock":
            config["time"]["stop_s"] = 1200
        times, template = basis(config)
        for kind in ("white", "OU"):
            for multiplier in SCALE_MULTIPLIERS:
                result = classify(template, covariance(times, scale*multiplier, kind), generator.spawn(4))
                results.append({"clock": clock, "sample_s": sampling, "phase": phase, "pulse_start_s": onset,
                                "noise": kind, "scale_multiplier": multiplier, **result})
    return {"data_kind": "synthetic, anchored to one-vesicle baseline fluctuation",
            "empirical_scale": empirical, "calibration_nulls_per_setting": N_CAL,
            "test_paths_per_class": N_TEST, "results": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("guv_csv", type=Path)
    args = parser.parse_args()
    print(json.dumps(study(args.guv_csv), indent=2, allow_nan=False))
