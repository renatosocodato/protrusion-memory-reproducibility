"""Contract checks for the installable exact model and reporter."""

from dataclasses import replace
import json
from pathlib import Path
import unittest

import numpy as np

from protrusion_memory.model import Model, Pulse, Reporter, compare_reporter, probabilities, simulate_pair


class ExactPairTests(unittest.TestCase):
    def test_default_matches_published_pair(self):
        a, b, optical = simulate_pair()
        np.testing.assert_allclose(optical, b[:, 0], atol=1e-12, rtol=0)
        np.testing.assert_allclose(a[:, 1:], b[:, 1:], atol=1e-12, rtol=0)
        self.assertAlmostEqual(np.max(abs(a[:, 0] - b[:, 0])), 0.12590721500489369)

    def test_custom_near_equal_rates_and_overlapping_pulses(self):
        model = Model(fast_s=10, slow_s=10.000001, deformation_s=10, duration_s=12,
                      gate=Pulse(2, 7, 0.2), actin=Pulse(4, 8, 0.3))
        for order in ("uv", "vu"):
            a, b, optical = simulate_pair(model, order)
            np.testing.assert_allclose(optical, b[:, 0], atol=1e-12, rtol=0)
            np.testing.assert_allclose(a[:, 1:], b[:, 1:], atol=1e-12, rtol=0)
            np.testing.assert_array_equal(a[:3], 0)
            gate_start = 2 if order == "uv" else 4
            np.testing.assert_array_equal(a[:gate_start + 1, 0], 0)
            self.assertGreater(a[gate_start + 1, 0], 0)


class ReporterTests(unittest.TestCase):
    def test_frozen_figure_values(self):
        table = Path(__file__).resolve().parents[1] / "figures/data/independent_gate_readout.json"
        if not table.is_file():
            self.skipTest("frozen figure inputs are distributed in the complete paper bundle")
        for row in json.loads(table.read_text())["rows"]:
            with self.subTest(row=row):
                result = compare_reporter(reporter=Reporter(
                    response_s=row["tau"], cadence_s=row["cadence"], mixing_bound=row["bound"],
                ))
                if result["status"] == "touching_means":
                    self.assertAlmostEqual(result["either_declaration_at_touching_mean"],
                                           row["correct"] + row["wrong"])
                else:
                    self.assertAlmostEqual(result["closest_mean_probabilities"]["correct"], row["correct"])
                    self.assertAlmostEqual(result["closest_mean_probabilities"]["wrong"], row["wrong"])

    def test_overlap_and_no_signal_have_no_classification_probability(self):
        touch = compare_reporter(reporter=Reporter(mixing_bound=0.68))
        self.assertEqual(touch["gap_sigma"], 0)
        self.assertAlmostEqual(touch["either_declaration_at_touching_mean"], 0.05)
        overlap = compare_reporter(reporter=Reporter(mixing_bound=0.8))
        self.assertEqual(overlap["status"], "overlapping_means")
        self.assertEqual(overlap["gap_sigma"], 0)
        self.assertIsNone(overlap["closest_mean_probabilities"])
        model = Model(gate=Pulse(20, 80, 0), actin=Pulse(140, 200, 0))
        self.assertEqual(compare_reporter(model)["status"], "zero_template")

    def test_decimal_touching_boundary(self):
        result = compare_reporter(Model(fast_s=0.1, slow_s=0.4), Reporter(mixing_bound=0.75))
        self.assertEqual(result["status"], "touching_means")
        self.assertEqual(result["gap_sigma"], 0)
        self.assertAlmostEqual(result["either_declaration_at_touching_mean"], 0.05)

    def test_signal_design_changes_the_decision(self):
        base = Reporter()
        frequent = compare_reporter(reporter=base)
        rare = compare_reporter(reporter=replace(base, cadence_s=60))
        self.assertGreater(frequent["gap_sigma"], rare["gap_sigma"])
        for distance in (0, 3, 6):
            self.assertAlmostEqual(sum(probabilities(distance)), 1)

    def test_invalid_input_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "noise_sd"):
            Reporter(noise_sd=0)
        with self.assertRaisesRegex(ValueError, "cadence"):
            compare_reporter(reporter=Reporter(cadence_s=302))
        with self.assertRaisesRegex(ValueError, "pulse start"):
            Pulse(20.5, 80, 0.3)


if __name__ == "__main__":
    unittest.main()
