"""Compare mode information when angular gains are known up to a shared scale or free."""

import argparse
from copy import deepcopy
import json
import math
from pathlib import Path

import numpy as np

from protrusion_memory.spatial import simulate

from guv_geometry import SAMPLE_S

MODES = (2, 3, 4)
N_CAL, N_TEST = 6000, 4000


def ring(cadence, onset, diffusion):
    config = {
        "schema": "ring-v1",
        "ring": {"nodes": 64, "modes": list(MODES)},
        "time": {"stop_s": 320, "sample_s": cadence},
        "parameters": {"pulse": 0, "g": .05, "direct": 0},
        "modulations": [{"parameter": "pulse", "start_s": onset, "stop_s": onset+5,
                         "operation": "set", "value": 1}],
        "inputs": {"drive": {"parameter": "pulse", "cos": {str(n): 1 for n in MODES}}},
        "fields": {
            "actin": {"decay_s_inv": .125, "diffusion_rad2_per_s": diffusion, "drives": {"drive": 1}},
            "shape": {"decay_s_inv": .05, "diffusion_rad2_per_s": .002,
                      "coupling": {"actin": "g"}, "drives": {"drive": "direct"}},
        },
        "reporters": {f"mode{n}": {"field": "shape", "mode": n} for n in MODES},
        "observations": {f"mode{n}": {"unit": "um", "mix": {f"mode{n}": 1},
                                     "delay_s": .6, "exposure_s": .7} for n in MODES},
    }
    return config


def templates(cadence, onset, diffusivity):
    mediated = ring(cadence, onset, diffusivity)
    direct = deepcopy(mediated)
    discrete_k2 = 4*math.sin(math.pi*2/64)**2 / (2*math.pi/64)**2
    direct["parameters"].update(g=0, direct=.05/(.125+diffusivity*discrete_k2))
    acquired = [simulate(c)["observations"] for c in (mediated, direct)]
    return np.asarray([[run[f"mode{n}"]["values"] for run in acquired] for n in MODES]).transpose(0, 2, 1)


def scores(z, grams, shared):
    if shared:
        return .5*np.maximum(z.sum(axis=1), 0)**2 / np.maximum(np.diagonal(grams.sum(axis=0)), 1e-24)
    return (.5*np.maximum(z, 0)**2 / np.maximum(np.diagonal(grams, axis1=1, axis2=2), 1e-24)).sum(axis=1)


def evaluate(basis, shared, streams):
    grams = np.einsum("mti,mtj->mij", basis, basis)
    eig, vec = np.linalg.eigh(grams)
    transform = vec * np.sqrt(np.maximum(0, eig))[:,None,:]

    def draw(path, n, truth=None):
        values = np.einsum("mij,pmj->pmi", transform,
                           np.random.default_rng(path).normal(size=(n,len(basis),2)))
        return values if truth is None else values+grams[:,:,truth]

    threshold = float(np.quantile(scores(draw(streams[0], N_CAL), grams, shared).max(axis=1),
                                  .95, method="higher"))
    matrix = {}
    for i, truth in enumerate(("null", "mediated", "direct")):
        ll = scores(draw(streams[i+1], N_TEST, None if i==0 else i-1), grams, shared)
        detected = ll.max(axis=1) > threshold
        chosen = detected & (abs(ll[:,0]-ll[:,1]) >= 2)
        counts = np.bincount(np.where(chosen,ll.argmax(axis=1)+1,0), minlength=3)
        matrix[truth] = {**dict(zip(("unresolved", "mediated", "direct"), counts.tolist())),
                         "detected": int(detected.sum())}
    return matrix


def study(geometry):
    sigmas = np.asarray([np.sqrt(np.mean(np.square(geometry["summary"]["higher_modes"][str(n)]
                                                   ["baseline_component_sd_um"]))) for n in MODES])
    designs = [(clock, cadence, (130+p if clock=="1s" else cadence*(2+p)), p)
               for clock, cadence in (("1s", 1.), ("G01", SAMPLE_S))
               for p in (0.,.25,.50,.75,.99)]
    designs.append(("G01_near_frame", SAMPLE_S, 130., None))
    rng = np.random.SeedSequence(917201)
    rows = []
    for clock, cadence, onset, phase in designs:
        for diffusion in (0., .015):
            raw = templates(cadence, onset, diffusion)
            base = (raw - raw.mean(axis=1, keepdims=True)) / sigmas[:,None,None]
            for label, shaped, paired in (("mode2_only", base[:1], True), ("common_gain", base, True),
                                          ("mode_specific_gains", base, False)):
                rows.append({"clock": clock, "pulse_s": onset, "phase": phase,
                             "actin_diffusivity": diffusion, "model": label,
                             "confusion": evaluate(shaped, paired, rng.spawn(4))})
    return {"data_kind": "synthetic spatial pulse, radial fluctuation scale from one vesicle",
            "mode_order": list(MODES), "baseline_sd_um": sigmas.tolist(),
            "calibration_nulls": N_CAL, "test_paths_per_class": N_TEST, "results": rows}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("geometry", type=Path)
    args = parser.parse_args()
    print(json.dumps(study(json.loads(args.geometry.read_text())), indent=2, allow_nan=False))
