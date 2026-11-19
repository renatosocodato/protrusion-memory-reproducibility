"""Seeded process, camera and event-aligned shutter checks."""

from copy import deepcopy
import unittest

import numpy as np

from protrusion_memory.portable import compare, ensemble, simulate
from protrusion_memory.stochastic import integrate


def model():
    return {
        "schema": "sde-v1",
        "time": {"stop_s": 1, "sample_s": .25},
        "states": {"x": 0, "y": 0},
        "equations": {"x": "0", "y": "0"},
        "diffusion": {"shared": {"x": "1", "y": "2"}},
        "stochastic": {"seed": 18439, "max_step_s": .05},
        "reporters": {"x": {"expression": "x"}, "y": {"expression": "y"}},
        "observations": {"camera": {"unit": "au", "mix": {"x": 1, "y": .5},
                                     "delay_s": .07, "exposure_s": .12, "noise_sd": .2}},
    }


class StochasticTests(unittest.TestCase):
    def test_wiener_driver_is_shared_and_replay_is_exact(self):
        config = model()
        a = simulate(config)
        self.assertEqual(a, simulate(config))
        np.testing.assert_allclose(a["states"]["y"], 2*np.asarray(a["states"]["x"]), atol=0)
        self.assertEqual(a["observations"]["camera"]["time_s"], [.25, .5, .75, 1.])
        config["stochastic"]["seed"] += 1
        self.assertNotEqual(a["states"]["x"], simulate(config)["states"]["x"])

    def test_explicit_detector_noise_keeps_process_replay(self):
        config = model()
        clean = simulate(config)
        config["stochastic"]["measurement_noise"] = True
        noisy = simulate(config)
        self.assertEqual(clean["states"], noisy["states"])
        self.assertNotEqual(clean["observations"], noisy["observations"])

    def test_shutter_ends_after_noninteger_delay_and_respects_modulation(self):
        config = model()
        config["diffusion"]["shared"] = {"x": "0"}
        config["parameters"] = {"u": 0}
        config["reporters"]["x"] = {"expression": "u"}
        config["modulations"] = [{"parameter": "u", "start_s": .15, "stop_s": .20,
                                  "operation": "set", "value": 1}]
        recorded = simulate(config)["observations"]["camera"]["values"]
        np.testing.assert_allclose(recorded, [.03/.12, 0, 0, 0], atol=1e-10)
        legacy = deepcopy(config)
        legacy["schema"] = "ode-v2"
        del legacy["stochastic"], legacy["diffusion"]
        np.testing.assert_allclose(recorded, simulate(legacy)["observations"]["camera"]["values"], atol=1e-10)

    def test_ensemble_matches_ou_end_moments(self):
        config = model()
        config["equations"]["x"] = "-x"
        config["observations"]["camera"] = {"unit": "au", "mix": {"x": 1}}
        config["stochastic"]["max_step_s"] = .05
        out = ensemble(config, 256)
        camera = out["observations"]["camera"]
        self.assertEqual(len(set(out["path_seeds"])), 256)
        exact_sd = np.sqrt((1-np.exp(-2))/2)
        self.assertLess(abs(camera["mean"][-1]), .15)
        self.assertLess(abs(camera["sd"][-1] - exact_sd), .10)
        with self.assertRaisesRegex(ValueError, "deterministic"):
            compare(config, config)


    def test_ou_refinement_with_shared_wiener_path(self):
        rng = np.random.default_rng(9465)
        increments = rng.normal(0, 1/8, size=(64, 16))

        class FixedIncrements:
            def __init__(self, drivers):
                self.drivers = iter(drivers)

            def standard_normal(self, count):
                return next(self.drivers)

        def sample(n):
            dt = 1/n
            blocks = increments.reshape(n, 64//n, 16).sum(axis=1)
            returns = integrate(lambda t, y: -y, lambda t, y: np.eye(16),
                                np.zeros(16), np.array([0., 1.]), dt,
                                FixedIncrements(blocks/np.sqrt(dt)))
            return returns.path[-1]

        reference = sample(64)
        self.assertLess(np.linalg.norm(sample(16)-reference),
                        .6*np.linalg.norm(sample(8)-reference))
