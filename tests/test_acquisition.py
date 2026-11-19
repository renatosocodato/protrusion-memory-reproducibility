"""Analytical checks on shutter timing, acquisition and legacy compatibility."""

from copy import deepcopy
import unittest

import numpy as np

from protrusion_memory.portable import compare, simulate


def model():
    return {"schema": "ode-v2",
            "time": {"stop_s": 3, "sample_s": .5},
            "states": {"x": 0},
            "equations": {"x": "0"},
            "parameters": {"pulse": 0},
            "modulations": [{"parameter": "pulse", "start_s": .75, "stop_s": 2.25,
                             "operation": "set", "value": 1}],
            "reporters": {"source": {"expression": "pulse"}, "other": {"expression": "2"}},
            "observations": {"a": {"unit": "au", "mix": {"source": 1}}}}


class AcquisitionTests(unittest.TestCase):
    def test_exposure_delay_and_mixing_analytical_rectangles(self):
        config = model()
        config["observations"]["a"].update(mix={"source": 1, "other": .1},
                                          exposure_s=.5, delay_s=.25)
        out = simulate(config)["observations"]["a"]
        self.assertEqual(out["time_s"], [1, 1.5, 2, 2.5, 3])
        np.testing.assert_allclose(out["values"], [.2, 1.2, 1.2, 1.2, .2], atol=1e-10)

    def test_reporter_matches_step_and_washout(self):
        config = model()
        config["reporters"]["source"].update(response_s=.2, initial=0)
        obs = simulate(config)["observations"]["a"]
        t = np.asarray(obs["time_s"])
        expected = (1-np.exp(-np.clip(t-.75, 0, 1.5)/.2))*np.exp(-np.maximum(t-2.25, 0)/.2)
        np.testing.assert_allclose(obs["values"], expected, atol=1e-8)
        self.assertEqual(obs["values"][:2], [0, 0])

    def test_acquisition_cannot_see_future_pulse(self):
        config = model()
        config["reporters"]["source"].update(response_s=.4, initial=.1)
        config["modulations"][0]["start_s"] = 1.25
        config["observations"]["a"].update(exposure_s=.5, delay_s=.5, operator="interval_rate")
        no_pulse = deepcopy(config)
        no_pulse["modulations"] = []
        a, b = [simulate(c)["observations"]["a"] for c in (config, no_pulse)]
        self.assertEqual(a["time_s"][0], 1.5)
        np.testing.assert_allclose(a["values"][:1], b["values"][:1], atol=1e-10)
        self.assertGreater(abs(a["values"][1]-b["values"][1]), .1)

    def test_instantaneous_identity_matches_v1(self):
        config = model()
        legacy = deepcopy(config)
        legacy["schema"] = "ode-v1"
        del legacy["reporters"]
        legacy["observations"]["a"] = {"unit": "au", "expression": "pulse"}
        self.assertEqual(compare(legacy, config)["channels"]["a"]["rmse"], 0)

    def test_v2_rejects_empty_support_and_missing_initial_state(self):
        config = model()
        config["observations"]["a"]["exposure_s"] = 4
        with self.assertRaisesRegex(ValueError, "leaves no samples"):
            simulate(config)
        config["reporters"]["source"]["response_s"] = 1
        with self.assertRaisesRegex(ValueError, "explicit initial"):
            simulate(config)
