"""Numerical and input contracts for portable models."""

from copy import deepcopy
import math
from pathlib import Path
import tomllib
import unittest

import numpy as np

from protrusion_memory.portable import compare, simulate

ROOT = Path(__file__).resolve().parents[1]


def sample():
    return {
        "schema": "ode-v1", "time": {"sample_s": 1, "stop_s": 3},
        "states": {"x": 0}, "parameters": {"u": 0},
        "equations": {"x": "u"},
        "modulations": [{"parameter": "u", "start_s": 0.5, "stop_s": 2,
                         "operation": "set", "value": 2}],
        "observations": {
            "slope": {"expression": "x", "operator": "interval_rate", "unit": "nm/s", "noise_sd": 2},
            "drive": {"expression": "u", "unit": "assigned"},
        },
    }


class PortableTests(unittest.TestCase):
    def test_off_grid_intervention_integrates_continuously(self):
        result = simulate(sample())
        np.testing.assert_allclose(result["states"]["x"], [0, 1, 3, 3], atol=1e-10)
        self.assertEqual(result["observations"]["drive"]["values"], [0, 2, 0, 0])
        np.testing.assert_allclose(result["observations"]["slope"]["values"], [1, 2, 0], atol=1e-10)
        self.assertEqual(result["observations"]["slope"]["time_s"], [1, 2, 3])

    def test_fractional_clock_samples_right_side_of_boundary(self):
        config = sample()
        config["time"] = {"sample_s": 0.3, "stop_s": 1.8}
        config["modulations"][0].update(start_s=0.9, stop_s=1.5)
        result = simulate(config)
        self.assertEqual(result["time_s"], [0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8])
        self.assertEqual(result["observations"]["drive"]["values"], [0, 0, 0, 2, 2, 0, 0])

    def test_modulations_are_relative_to_baseline(self):
        config = sample()
        config["parameters"]["u"] = 3
        config["modulations"][0].update(operation="scale")
        config["modulations"].append({"parameter": "u", "start_s": 2, "stop_s": 3,
                                      "operation": "add", "value": 1})
        result = simulate(config)
        np.testing.assert_allclose(result["states"]["x"], [0, 4.5, 10.5, 14.5], atol=1e-10)
        self.assertEqual(result["observations"]["drive"]["values"], [3, 6, 4, 3])

    def test_compare_requires_commensurate_readouts(self):
        config = sample()
        other = deepcopy(config)
        other.pop("modulations")
        result = compare(config, other)
        self.assertAlmostEqual(result["channels"]["slope"]["rmse"], math.sqrt(5/3))
        self.assertAlmostEqual(result["channels"]["slope"]["rmse_in_noise_sd"], math.sqrt(5/3) / 2)
        other["observations"]["slope"]["unit"] = "um/s"
        with self.assertRaisesRegex(ValueError, "units"):
            compare(config, other)

    def test_rejects_unsupported_expressions_and_conflicting_modulations(self):
        config = sample()
        config["equations"]["x"] = "x.real"
        with self.assertRaisesRegex(ValueError, "unsupported expression"):
            simulate(config)
        config = sample()
        config["modulations"] *= 2
        with self.assertRaisesRegex(ValueError, "overlap"):
            simulate(config)

    def test_benchmark_matches_independent_interval_average(self):
        with (ROOT / "examples/rac1-edge.toml").open("rb") as stream:
            config = tomllib.load(stream)
        channel = simulate(config)["observations"]["edge_velocity"]
        t = np.array(channel["time_s"])
        reference = 1 + 6 * (np.exp(-t / 30) - np.exp(-(t - 5) / 30))
        np.testing.assert_allclose(channel["values"], reference, atol=1e-7, rtol=0)


if __name__ == "__main__":
    unittest.main()
