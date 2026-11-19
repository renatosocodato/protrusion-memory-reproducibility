"""Independent calibration and signed preparation-level intervention controls."""

import csv
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from protrusion_memory.mode_test import ARMS, CHANNELS, Plan, calibrate, score


PLAN = """schema = "mode-test-v1"
mode = 2
axis = [1, 0]
window_s = [10, 30]
[split]
calibration = ["a", "b", "c"]
reserved = ["d", "e", "f"]
[criteria]
actin_min = 0.15
shape_min = 0.08
drive_max = 0.075
drive_gain = 0.55
delivery_bound = 0.025
[resampling]
seed = 15
draws = 200
"""


class ModeTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.matrix = np.array([[1, -.2], [.3, 1.1]])
        self.plan_file = self.root/"plan.toml"
        self.plan_file.write_text(PLAN)
        self.plan = Plan.load(self.plan_file)
        self.gain_data = self.root/"gains.csv"
        self.gain_file = self.root/"frozen.json"
        self.responses = self.root/"response.csv"

    def standards(self, vectors):
        with self.gain_data.open("w", newline="") as stream:
            records = csv.writer(stream)
            records.writerow(["preparation","mode","channel","known_cos","known_sin","read_cos","read_sin"])
            for name in self.plan.calibration:
                for c in CHANNELS:
                    for v in vectors:
                        records.writerow([name,2,c,*v,*(self.matrix @ v)])
        measured = calibrate(self.plan, self.gain_data, "1.9.0")
        self.gain_file.write_text(json.dumps(measured))
        return measured

    def measurements(self, h, u, a=.4, missing=False):
        with self.responses.open("w", newline="") as stream:
            records = csv.writer(stream)
            records.writerow(["preparation","mode","arm",*(f"{c}_{d}" for c in CHANNELS for d in ("cos","sin"))])
            for name in self.plan.reserved:
                for arm in ARMS[1:] if missing else ARMS:
                    sign = 1 if arm.endswith("plus") else -1
                    change = np.array([[a,0], [h,.01], [u,0]]) if arm.startswith("active") else np.zeros((3,2))
                    same = np.array([[1,2], [.1,.2], [2,.4]])
                    detector = ((same+sign*change) @ self.matrix.T).ravel()
                    records.writerow([name,2,arm,*detector])
        return score(self.plan, self.responses, self.gain_file, "1.9.0")

    def test_rotated_gain_and_triangle_bound(self):
        gains = self.standards(([1,0], [0,1], [-1,0]))
        np.testing.assert_allclose(gains["gain_matrices"]["shape"], self.matrix)
        report = self.measurements(h=.25, u=.04)
        self.assertEqual(report["decision"], "deformation_beyond_leakage")
        np.testing.assert_allclose(report["shape_interval_with_leakage"], [.203,.297])
        self.assertAlmostEqual(report["by_preparation"]["d"]["shape_perpendicular_odd"], .01)
        self.assertEqual(self.measurements(h=.25, u=.15)["decision"], "drive_not_matched")
        self.assertEqual(self.measurements(h=.25, u=.04, a=.01)["decision"], "actin_unresolved")

    def test_direct_only_and_missing_arm(self):
        self.standards(([1,0], [0,1], [-1,0]))
        self.assertEqual(self.measurements(h=0, u=.001)["decision"], "deformation_bounded_small")
        with self.assertRaisesRegex(ValueError, "four intervention arms"):
            self.measurements(h=0, u=0, missing=True)

    def test_singular_gain_and_overlap_fail(self):
        with self.assertRaisesRegex(ValueError, "span both"):
            self.standards(([1,0], [2,0], [-1,0]))
        self.plan_file.write_text(PLAN.replace('["d", "e", "f"]', '["a", "e", "f"]'))
        with self.assertRaisesRegex(ValueError, "disjoint"):
            Plan.load(self.plan_file)

    def test_drive_error_is_propagated_to_matching_and_leakage(self):
        self.plan_file.write_text(PLAN.replace("delivery_bound = 0.025", "delivery_bound = 0.025\ndrive_error = 0.04"))
        self.plan = Plan.load(self.plan_file)
        self.standards(([1,0], [0,1], [-1,0]))
        result = self.measurements(h=.25, u=.04)
        self.assertEqual(result["decision"], "drive_not_matched")
        self.assertAlmostEqual(result["mean_leakage_allowance"], .025+.55*.08)
