"""Targeted regression checks for the signed, pulse-conditioned edge comparison."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from protrusion_memory.calibration import weights
from protrusion_memory.forecast import Plan, edge_origins


PLAN = '''schema = "edge-forecast-v1"
sample_s = 5
horizon_frames = 1
prefix_frames = 2
motion_delay_frames = 0
actin_delay_frames = 1
motion_unit = "um/min"
[split]
calibration = ["a", "b"]
heldout = ["c", "d"]
[perturbations.a]
onset_frame = 3
offset_frame = 5
[perturbations.b]
onset_frame = 3
offset_frame = 5
[perturbations.c]
onset_frame = 4
offset_frame = 6
[perturbations.d]
onset_frame = 4
offset_frame = 6
'''


class EdgeForecastTests(unittest.TestCase):
    def test_edge_features_use_only_available_measurements(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "plan.toml"
            path.write_text(PLAN)
            plan = Plan(path)
        records = {("a", "one", "left"): {i: (i * 1.0, -i * 2.0) for i in range(8)}}
        series = records["a", "one", "left"]
        by_frame = {row.frame: row for row in edge_origins(records, plan)}
        row = by_frame[3]
        self.assertEqual((row.actin, row.motion, row.context, row.target),
                         (1.5, -5.0, (-3.0, 1.0, 0.0), -8.0))
        series[3] = (999.0, -6.0)
        series[4] = (999.0, -100.0)
        after = edge_origins(records, plan)[0]
        self.assertEqual((after.actin, after.motion, after.context), (row.actin, row.motion, row.context))
        self.assertNotEqual(after.target, row.target)
        del series[4]
        self.assertEqual([r.frame for r in edge_origins(records, plan)], [6])

    def test_trailing_velocity_can_be_undefined_at_first_frame(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "plan.toml"
            path.write_text(PLAN)
            plan = Plan(path)
        series = {i: (i * 1.0, -i * 2.0) for i in range(8)}
        series[0] = (0.0, None)
        rows = edge_origins({("a", "one", "left"): series}, plan)
        self.assertEqual(rows[0].frame, 3)
        self.assertEqual(rows[0].baseline, -3.0)
        series[1] = (1.0, None)
        with self.assertRaisesRegex(ValueError, "no eligible edge forecasts"):
            edge_origins({("a", "one", "left"): series}, plan)

    def test_edge_regions_do_not_inflate_cells_or_preparations(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "plan.toml"
            path.write_text(PLAN)
            plan = Plan(path)
        records = {key: {i: (i * 1.0, -i * 2.0) for i in range(n)}
                   for key, n in [(("a", "one", "left"), 8), (("a", "one", "right"), 7),
                                  (("a", "two", "left"), 7), (("b", "one", "left"), 10)]}
        rows = edge_origins(records, plan)
        self.assertAlmostEqual(sum(w for r, w in zip(rows, weights(rows)) if r.preparation == "a"), 0.5)
        totals = [sum(w for r, w in zip(rows, weights(rows))
                      if (r.preparation, r.unit, r.channel) == key) for key in records]
        np.testing.assert_allclose(totals, [0.125, 0.125, 0.25, 0.5])
