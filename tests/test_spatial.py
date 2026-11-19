"""Analytical periodic-grid transport and causal projected readouts."""

from copy import deepcopy
import math
import unittest

import numpy as np

from protrusion_memory.spatial import compare, simulate


def config(nodes=32):
    return {
        "schema": "ring-v1",
        "ring": {"nodes": nodes, "modes": [0, 2, 3]},
        "time": {"stop_s": 2, "sample_s": .25},
        "fields": {"a": {"decay_s_inv": 0, "diffusion_rad2_per_s": .2,
                         "initial_cos": {"0": 1, "2": 2},
                         "initial_sin": {"3": 3}}},
        "reporters": {"shape": {"field": "a", "mode": 2}},
        "observations": {"camera": {"unit": "um", "mix": {"shape": 1}}},
    }


class RingTests(unittest.TestCase):
    def test_discrete_diffusion_and_mean_conservation(self):
        model = config()
        out = simulate(model)
        time = np.array(out["time_s"])
        mode = out["modal_coefficients"]["a"]
        np.testing.assert_allclose(mode["0"]["cos"], 1)
        for n, component, initial in ((2, "cos", 2), (3, "sin", 3)):
            nodes = model["ring"]["nodes"]
            rate = .2*4*math.sin(math.pi*n/nodes)**2/(2*math.pi/nodes)**2
            np.testing.assert_allclose(mode[str(n)][component], initial*np.exp(-rate*time), rtol=1e-8)

    def test_grid_refinement_and_no_spectral_mixing_by_diffusion(self):
        exact = 2*math.exp(-.2*4*2)
        coarse, refined = [simulate(config(n))["modal_coefficients"]["a"] for n in (16, 32)]
        error = lambda out: abs(out["2"]["cos"][-1]-exact)
        self.assertLess(error(refined), .3*error(coarse))
        self.assertTrue(all(x == 0 for x in refined["3"]["cos"]))

    def test_camera_exposure_and_delay_apply_after_mode_projection(self):
        model = config()
        model["observations"]["camera"].update(exposure_s=.25, delay_s=.5)
        out = simulate(model)["observations"]["camera"]
        t = np.asarray(out["time_s"])
        self.assertEqual(t[0], .75)
        rate = .2*4*math.sin(2*math.pi/32)**2/(2*math.pi/32)**2
        expected = 8/rate*(np.exp(-rate*(t-.75))-np.exp(-rate*(t-.5)))
        np.testing.assert_allclose(out["values"], expected, rtol=1e-8)
        self.assertEqual(compare(model, deepcopy(model))["channels"]["camera"]["rmse"], 0)

    def test_nyquist_and_inconsistent_detectors_rejected(self):
        model = config(16)
        model["ring"]["modes"].append(8)
        with self.assertRaisesRegex(ValueError, "mode must"):
            simulate(model)
        left, right = config(), config()
        right["reporters"]["shape"]["mode"] = 3
        with self.assertRaisesRegex(ValueError, "same reporter"):
            compare(left, right)
