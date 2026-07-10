# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Param-modal line/station selectors: options, auto-injection, dependence."""
from __future__ import annotations

import importlib.util
import unittest

_HAS_DASH = importlib.util.find_spec("dash") is not None


@unittest.skipUnless(_HAS_DASH, "Dash not installed")
class TestParamLineStation(unittest.TestCase):
    def _mod(self):
        from pycsamt.app.agent_master.callbacks import params
        return params

    def test_line_station_options(self):
        p = self._mod()
        groups = {
            "L18": ["/d/22-001A.edi", "/d/22-002A.edi"],
            "L22": ["/d/x1.edi", "/d/x2.edi", "/d/x3.edi"],
        }
        line_opts, station_opts, line_to_st = p._line_station_options(
            groups, ""
        )
        self.assertEqual([o["value"] for o in line_opts], ["L18", "L22"])
        self.assertEqual(len(station_opts), 5)
        self.assertEqual(set(line_to_st), {"L18", "L22"})

    def test_line_field_injected_only_when_multiline(self):
        p = self._mod()
        fields = [dict(p._PLOT_FIELD_STATIONS)]
        # >1 line → a "lines" selector is prepended
        multi, _ = p._prepare_dynamic_fields(
            fields, {"A": ["/d/a.edi"], "B": ["/d/b.edi"]}, "", []
        )
        self.assertEqual([f["key"] for f in multi], ["lines", "stations"])
        # single line → no line selector
        one, _ = p._prepare_dynamic_fields(
            fields, {"A": ["/d/a.edi"]}, "", []
        )
        self.assertEqual([f["key"] for f in one], ["stations"])

    def test_stations_depend_on_selected_lines(self):
        p = self._mod()
        line_to_st = {"L18": ["s1", "s2"], "L22": ["s3", "s4", "s5"]}
        # one line → only its stations
        opts = p._station_options_for_lines(["L22"], line_to_st)
        self.assertEqual([o["value"] for o in opts], ["s3", "s4", "s5"])
        # no selection → all stations
        allo = p._station_options_for_lines([], line_to_st)
        self.assertEqual(len(allo), 5)

    def test_stations_field_is_multiselect(self):
        p = self._mod()
        self.assertEqual(p._PLOT_FIELD_STATIONS["type"], "multiselect")
        self.assertEqual(p._PLOT_FIELD_LINES["type"], "multiselect")


if __name__ == "__main__":
    unittest.main(verbosity=2)
