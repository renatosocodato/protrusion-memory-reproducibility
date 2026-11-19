"""Mechanical calibration uses a third, independent preparation split."""

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from benchmarks.delivery_calibration.ring_controls import CASES, MECHANICAL, VERSION, example, write_plan
from protrusion_memory import mode_test


class DeliveryCalibration(unittest.TestCase):
    def test_transferred_mechanics_and_unqualified_control(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch)
            response = np.array([.36, 0.])
            blocked = np.zeros(2)
            report = example(path, CASES["transferred_leakage"], response, blocked, np.random.SeedSequence(73))
            self.assertEqual(report["mechanical_calibration"]["status"], "qualified")
            self.assertGreater(report["mechanical_calibration"]["mechanical_allowance"], .24)
            interval = report["shape_interval_with_leakage"]
            self.assertLess(interval[0], 0)
            self.assertGreater(interval[1], 0)
            self.assertEqual(report["decision"], "deformation_unresolved")

            report = example(path, CASES["incomplete_block"], response, np.array([.22, .1]),
                             np.random.SeedSequence(74))
            self.assertEqual(report["decision"], "mechanical_control_unqualified")
            self.assertIsNone(report["shape_interval_with_leakage"])

    def test_frozen_control_is_bound_to_exact_gain_fit_and_preparations(self):
        with tempfile.TemporaryDirectory() as scratch:
            path = Path(scratch)
            example(path, CASES["mediated"], np.array([.36, .25]), np.zeros(2), np.random.SeedSequence(75))
            plan = mode_test.Plan.load(path / "plan.toml")
            gains, reserved, control = (path / name for name in ("angular.json", "reserved.csv", "mechanical.json"))
            with self.assertRaisesRegex(ValueError, "requires a separate mechanical calibration"):
                mode_test.score(plan, reserved, gains, VERSION)
            frozen = json.loads(control.read_text())
            frozen["angular_calibration_sha256"] = "other"
            control.write_text(json.dumps(frozen))
            with self.assertRaisesRegex(ValueError, "does not match"):
                mode_test.score(plan, reserved, gains, VERSION, control)
            with self.assertRaisesRegex(ValueError, "disjoint"):
                write_plan(path / "plan.toml", mechanical=[*MECHANICAL[:-1], plan.reserved[0]])
