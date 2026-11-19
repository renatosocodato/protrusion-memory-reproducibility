"""Preparation isolation, time availability and forecast fitting."""

import csv
import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

from protrusion_memory.forecast import Plan, eligible_origins, fit, score


PLAN = """schema = "forecast-v1"
sample_s = 300
horizon_frames = 1
prefix_frames = 2
motion_delay_frames = 1
motion_unit = "um2"
[split]
calibration = ["a"]
heldout = ["b"]
"""


class ForecastTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.plan_path = self.root / "plan.toml"
        self.plan_path.write_text(PLAN)
        self.plan = Plan(self.plan_path)

    def write_csv(self, name, p, duplicate=False):
        path = self.root / name
        exp = [0, 0, 0, 1, 2, -1, 0.5, -0.5, 1, 0, 2]
        rows = [[p, "cell", t, math.exp(x), 10 + (2 * exp[t - 1] if t else 0)]
                for t, x in enumerate(exp)]
        if duplicate:
            rows.append(rows[-1])
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["preparation", "unit", "frame", "actin", "protrusion"])
            writer.writerows(rows)
        return path

    def test_recovers_known_future_signal_and_scores_frozen_fit(self):
        fit_data = self.write_csv("train.csv", "a")
        test_data = self.write_csv("test.csv", "b")
        frozen = self.root / "fit.json"
        result = fit(self.plan, fit_data, "1.4.0")
        self.assertAlmostEqual(result["coefficients"]["actin_modulated"][0], 0)
        self.assertAlmostEqual(result["coefficients"]["actin_modulated"][1], 2)
        frozen.write_text(json.dumps(result))
        tested = score(self.plan, test_data, frozen, "1.4.0")
        self.assertLess(tested["candidates"]["actin_modulated"]["scaled_rmse"], 1e-12)
        self.assertGreater(tested["candidates"]["stationary"]["scaled_rmse"], 0)
        with self.assertRaisesRegex(ValueError, "unexpected preparation"):
            score(self.plan, fit_data, frozen, "1.4.0")

    def test_future_measurements_do_not_change_origin_features(self):
        series = {t: (t * 0.01, t + 10) for t in range(8)}
        key = ("prep", "cell")
        before = {row.frame: row for row in eligible_origins({key: series}, self.plan)}[3]
        series[4] = (500, 10000)
        series[5] = (400, 9999)
        after = {row.frame: row for row in eligible_origins({key: series}, self.plan)}[3]
        self.assertEqual((before.baseline, before.actin, before.motion),
                         (after.baseline, after.actin, after.motion))
        self.assertNotEqual(before.target, after.target)

    def test_prefix_availability_and_no_bridging_gaps(self):
        plan = SimpleNamespace(prefix=3, motion_delay=2, horizon=1)
        series = {t: (0.1, t + 10) for t in (0, 1, 2, 3, 4, 5, 7, 8, 9)}
        origins = eligible_origins({("prep", "cell"): series}, plan)
        self.assertEqual([row.frame for row in origins], [4])
        self.assertEqual(origins[0].baseline, 11)

    def test_duplicate_frames_and_overlapping_splits_fail(self):
        duplicate = self.write_csv("duplicate.csv", "a", duplicate=True)
        with self.assertRaisesRegex(ValueError, "duplicate frame"):
            fit(self.plan, duplicate, "1.4.0")
        self.plan_path.write_text(PLAN.replace('heldout = ["b"]', 'heldout = ["a"]'))
        with self.assertRaisesRegex(ValueError, "overlap"):
            Plan(self.plan_path)


if __name__ == "__main__":
    unittest.main()
