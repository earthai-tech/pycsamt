# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Regression: metric questions must compute on the *selected* line (labelled by
line name), and auto-detected line groups must be persisted as the loaded
groups — not collapsed to the upload folder.

Guards the bug where "strike of line 2" → pick L22 still reported
"strike for am_edi_folder_… from 128 stations" (whole folder).
"""

from __future__ import annotations

import importlib.util
import os
import tempfile
import unittest

_HAS_DASH = importlib.util.find_spec("dash") is not None


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestMetricLineResolution(unittest.TestCase):
    def _store(self, **extra):
        s = {
            "path": "/up/am_edi_folder_x",
            "groups": {
                "L18": ["a", "b"],
                "L22": ["c", "d", "e"],
                "L34": ["f"],
            },
        }
        s.update(extra)
        return s

    def test_selected_line_wins(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _resolve_metric_targets,
        )

        t = _resolve_metric_targets(
            "give me the strike direction of line 2",
            self._store(selected_lines=["L22"]),
            {},
            False,
        )
        self.assertEqual([(l, len(s)) for l, s in t], [("L22", 3)])

    def test_all_lines_overrides_stale_selection(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _resolve_metric_targets,
        )

        t = _resolve_metric_targets(
            "strike for all lines",
            self._store(selected_lines=["L22"]),
            {},
            True,
        )
        self.assertEqual(sorted(l for l, _ in t), ["L18", "L22", "L34"])

    def test_named_line_resolves_without_picker(self):
        from pycsamt.app.agent_master.callbacks.chat import (
            _resolve_metric_targets,
        )

        t = _resolve_metric_targets("strike of L34", self._store(), {}, False)
        self.assertEqual([(l, len(s)) for l, s in t], [("L34", 1)])


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestDetectLinesToFiles(unittest.TestCase):
    def test_groups_by_station_id_prefix(self):
        from pycsamt.app.agent_master.callbacks.edi import (
            _detect_lines_to_files,
        )

        d = tempfile.mkdtemp()
        for fn in ["22-001A.edi", "22-002A.edi", "18-001.edi", "34-009.edi"]:
            with open(os.path.join(d, fn), "w") as fh:
                fh.write("x")
        g = _detect_lines_to_files(d)
        self.assertEqual(
            {k: len(v) for k, v in g.items()}, {"L22": 2, "L18": 1, "L34": 1}
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
