"""Check whether camera acquisition preserves the direct-versus-mediated pulse signature."""

import json
from copy import deepcopy

import numpy as np

from protrusion_memory.portable import simulate


def pulse(cadence, exposure, delay, onset=130):
    return {
        "schema": "ode-v2",
        "time": {"stop_s": 320, "sample_s": cadence},
        "states": {"a": 0, "h": 0},
        "parameters": {"u": 0, "g": 0, "direct": 0},
        "equations": {"a": "u - a/8", "h": "g*a + direct*u - h/20"},
        "modulations": [{"parameter": "u", "start_s": onset, "stop_s": onset + 5, "operation": "set", "value": 1}],
        "reporters": {
            "actin": {"expression": "a", "response_s": 4, "initial": 0},
            "membrane": {"expression": "h"},
            "actuator_signal": {"expression": "u"},
        },
        "observations": {
            "actin": {"unit": "au", "mix": {"actin": 1, "actuator_signal": 0.1},
                      "exposure_s": exposure, "delay_s": delay},
            "deformation": {"unit": "um", "mix": {"membrane": 1},
                            "exposure_s": exposure, "delay_s": delay},
        },
    }


def basis(config):
    columns = []
    for parameter, coefficient in (("g", 0.05), ("direct", 0.4)):
        case = deepcopy(config)
        case["parameters"][parameter] = coefficient
        readout = simulate(case)["observations"]["deformation"]
        columns.append(readout["values"])
    return np.asarray(readout["time_s"]), np.array(columns).T


def recover(config):
    time, acquired = basis(config)
    naive_config = deepcopy(config)
    for item in naive_config["observations"].values():
        item.update(exposure_s=0, delay_s=0)
    naive_time, naive = basis(naive_config)
    naive = naive[np.searchsorted(naive_time, time)]
    synth = np.array([[1., 0.], [0., 1.], [1., 1.]])
    recovered = {}
    for name, matrix in (("acquisition_aware", acquired), ("instantaneous", naive)):
        coef, _, rank, sv = np.linalg.lstsq(matrix, acquired @ synth.T, rcond=None)
        recovered[name] = {"rank": int(rank), "relative_smallest_singular_value": float(sv[-1] / sv[0]),
                           "estimated_components": coef.T.tolist()}
    separation = np.linalg.norm(acquired[:, 0] - acquired[:, 1])
    if recovered["acquisition_aware"]["rank"] == 2:
        np.testing.assert_allclose(recovered["acquisition_aware"]["estimated_components"], synth, atol=1e-9)
    return {"samples": len(time), "equal_area_case_l2_shape_difference_um": float(separation),
            "acquisition": recovered}


if __name__ == "__main__":
    designs = {"fast_1s": pulse(1, 0.7, 0.6), "guv_clock": pulse(65.58921813964844, 0.7, 0.6)}
    cadence = 65.58921813964844
    phase_scan = []
    for fraction in (0.0, 0.25, 0.50, 0.75, 0.99):
        onset = cadence * (2 + fraction)
        result = recover(pulse(cadence, 0.7, 0.6, onset=onset))
        phase_scan.append({"onset_fraction_of_frame": fraction,
                           "relative_smallest_singular_value":
                           result["acquisition"]["acquisition_aware"]["relative_smallest_singular_value"],
                           "l2_shape_difference_um": result["equal_area_case_l2_shape_difference_um"]})
    print(json.dumps({"data_kind": "synthetic assigned parameters",
                      "coefficient_order": ["actin_mediated", "direct"],
                      "challenge_order": ["mediated_only", "direct_only", "both"],
                      "designs": {name: recover(config) for name, config in designs.items()},
                      "guv_clock_phase_scan": phase_scan},
                     indent=2, allow_nan=False))
