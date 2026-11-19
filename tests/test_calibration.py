"""Calibration and independent-scoring contracts."""

import csv
import json
from pathlib import Path
import tempfile
import unittest

from protrusion_memory.calibration import Plan, fit, score


ENDPOINT = """schema = "ode-v1"
[time]
stop_s = 1
sample_s = 1
[states]
x = 0
[parameters]
shift = 0
treatment = 0
[equations]
x = "0"
[observations.endpoint]
unit = "assigned"
expression = "shift * treatment"
"""
PLAN = """schema = "calibration-v1"
[split]
calibration = ["a", "b"]
heldout = ["c", "d"]
[scenarios.treated]
treatment = 1
[channels.endpoint]
unit = "assigned"
scale = 1
[models.stationary]
path = "endpoint.toml"
[models.modulated]
path = "endpoint.toml"
[models.modulated.fit.shift]
initial = 0
lower = 0
upper = 3
"""


class CalibrationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.plan_file = self.root / "plan.toml"
        self.plan_file.write_text(PLAN)
        self.model_file = self.root / "endpoint.toml"
        self.model_file.write_text(ENDPOINT)
        self.plan = Plan(self.plan_file)
        self.train = self.csv("train.csv", [("a", "1", 1), ("a", "2", 1), ("b", "1", 1)])
        self.test = self.csv("test.csv", [("d", "1", 1), ("c", "1", 1)])

    def csv(self, filename, data):
        path = self.root / filename
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["preparation", "unit", "scenario", "channel", "time_s", "value", "measurement_unit"])
            for p, u, value in data:
                writer.writerow([p, u, "treated", "endpoint", 0, value, "assigned"])
        return path

    def frozen(self):
        path = self.root / "fit.json"
        path.write_text(json.dumps(fit(self.plan, self.train, "1.3.0")))
        return path

    def test_recovery_starting_at_bound_and_equal_preparations(self):
        result = fit(self.plan, self.train, "1.3.0")
        self.assertAlmostEqual(result["candidates"]["modulated"]["parameters"]["shift"], 1)
        unbalanced = self.csv("unbalanced.csv", [("a", str(i), 0) for i in range(1, 20)] + [("b", "1", 2)])
        result = fit(self.plan, unbalanced, "1.3.0")
        self.assertAlmostEqual(result["candidates"]["modulated"]["parameters"]["shift"], 1)

    def test_score_uses_frozen_fit_and_stationary_can_win(self):
        frozen = self.frozen()
        before = frozen.read_bytes()
        report = score(self.plan, self.test, frozen, "1.3.0")
        self.assertAlmostEqual(report["modulated_minus_stationary_scaled_rmse"], -1)
        zero = self.csv("zero.csv", [("c", "1", 0), ("d", "1", 0)])
        report = score(self.plan, zero, frozen, "1.3.0")
        self.assertAlmostEqual(report["modulated_minus_stationary_scaled_rmse"], 1)
        self.assertEqual(frozen.read_bytes(), before)

    def test_split_and_hash_guards(self):
        frozen = self.frozen()
        with self.assertRaisesRegex(ValueError, "unexpected preparation"):
            fit(self.plan, self.test, "1.3.0")
        with self.assertRaisesRegex(ValueError, "unexpected preparation"):
            score(self.plan, self.train, frozen, "1.3.0")
        self.model_file.write_text(ENDPOINT + "\n# changed\n")
        with self.assertRaisesRegex(ValueError, "frozen fit does not match"):
            score(Plan(self.plan_file), self.test, frozen, "1.3.0")

    def test_duplicate_records_and_overlapping_splits_are_rejected(self):
        duplicate = self.csv("duplicate.csv", [("a", "1", 1), ("a", "1", 2), ("b", "1", 1)])
        with self.assertRaisesRegex(ValueError, "duplicate observation"):
            fit(self.plan, duplicate, "1.3.0")
        self.plan_file.write_text(PLAN.replace('heldout = ["c", "d"]', 'heldout = ["a", "d"]'))
        with self.assertRaisesRegex(ValueError, "overlap"):
            Plan(self.plan_file)

    def test_scenario_cannot_override_a_fitted_parameter(self):
        self.plan_file.write_text(PLAN.replace("treatment = 1", "treatment = 1\nshift = 1"))
        with self.assertRaisesRegex(ValueError, "overrides a fitted parameter"):
            Plan(self.plan_file)

    def test_operator_mismatch_is_rejected(self):
        (self.root / "other.toml").write_text(ENDPOINT + 'operator = "interval_rate"\n')
        self.plan_file.write_text(PLAN.replace('[models.modulated]\npath = "endpoint.toml"',
                                               '[models.modulated]\npath = "other.toml"'))
        with self.assertRaisesRegex(ValueError, "same observation operator"):
            Plan(self.plan_file)


if __name__ == "__main__":
    unittest.main()
