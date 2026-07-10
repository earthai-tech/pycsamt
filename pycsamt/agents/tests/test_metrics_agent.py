# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MetricsAgent — inline computed-value answers about a line."""
from __future__ import annotations

import os
import unittest

from pycsamt.agents.metrics import (  # noqa: E402
    METRIC_KINDS,
    MetricsAgent,
    looks_like_metric_query,
    parse_metric_request,
)

_DATA = os.path.join("data", "3edis")
_HAS_DATA = os.path.isdir(_DATA)


class TestMetricParsing(unittest.TestCase):
    def test_parse_single(self):
        kinds, all_lines = parse_metric_request("what is the strike of L22PLT?")
        self.assertEqual(kinds, ["strike"])
        self.assertFalse(all_lines)

    def test_parse_all_lines(self):
        kinds, all_lines = parse_metric_request("azimuth for all lines")
        self.assertEqual(kinds, ["azimuth"])
        self.assertTrue(all_lines)

    def test_parse_summary(self):
        kinds, _ = parse_metric_request("tell me about this line")
        self.assertEqual(kinds, ["summary"])

    def test_parse_summary_yields_to_specific(self):
        # a summary word next to a specific metric prefers the specific one
        kinds, _ = parse_metric_request("give me an overview of the strike")
        self.assertEqual(kinds, ["strike"])

    def test_metric_query_detection(self):
        self.assertTrue(looks_like_metric_query("what's the strike of L22PLT?"))
        self.assertTrue(looks_like_metric_query("how many stations in the survey"))
        self.assertTrue(looks_like_metric_query("azimuth of all lines"))

    def test_not_metric_query(self):
        # visual / workflow phrasings must not be treated as value questions
        self.assertFalse(looks_like_metric_query("plot the strike rose"))
        self.assertFalse(looks_like_metric_query("run phase tensor analysis"))
        self.assertFalse(looks_like_metric_query("strike analyzer"))
        self.assertFalse(looks_like_metric_query("hello"))

    def test_kinds_complete(self):
        self.assertEqual(len(METRIC_KINDS), 10)


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestMetricsAgent(unittest.TestCase):
    def _run(self, **params):
        r = MetricsAgent().execute({"path": _DATA, **params})
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertTrue(r.summary.strip())
        return r

    def test_each_kind(self):
        for k in METRIC_KINDS:
            r = self._run(kinds=[k], label="3edis")
            self.assertIn(k, (r.data or {}).get("values", {}))

    def test_multi_kind(self):
        r = self._run(kinds=["stations", "periods"], label="L")
        self.assertIn("stations", r.data["values"])
        self.assertIn("periods", r.data["values"])

    def test_summary(self):
        r = self._run(kinds=["summary"])
        self.assertTrue(len(r.summary) > 20)

    def test_no_data(self):
        self.assertEqual(MetricsAgent().execute({"kinds": ["strike"]}).status,
                         "failed")

    def test_station_filter(self):
        r = self._run(kinds=["stations"], stations="new_E1_1")
        self.assertIn("1 station", r.data["values"]["stations"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
