# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ToolAgent (strike / dimensionality / validator)."""
from __future__ import annotations

import os
import unittest

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pycsamt.agents.tooling import (  # noqa: E402
    ToolAgent,
    TOOL_KINDS,
    _df_to_text,
    _circular_strike_mean,
    _ll_to_utm,
    _get_latlon,
    _safe_filename,
    _EXPORT_BUNDLES,
)

_DATA = os.path.join("data", "3edis")
_HAS_DATA = os.path.isdir(_DATA)


class TestToolHelpers(unittest.TestCase):
    def test_kinds(self):
        self.assertEqual(
            set(TOOL_KINDS),
            {"strike", "dimensionality", "validator",
             "coords", "elevation", "converter", "batch_export",
             "freq_editor", "layered_model"},
        )

    def test_layered_model_dataless(self):
        # layered_model needs no EDI dataset.
        import matplotlib.pyplot as plt
        r = ToolAgent().execute({
            "kind": "layered_model", "preset": "custom",
            "resistivities": "100, 10, 500", "thicknesses": "300, 800",
        })
        self.assertEqual(r.status, "success", msg=r.summary)
        figs = r.data.get("figures") or {}
        self.assertTrue(figs)
        for f in figs.values():
            self.assertIsInstance(f, plt.Figure)
        plt.close("all")

    def test_layered_model_preset(self):
        import matplotlib.pyplot as plt
        r = ToolAgent().execute({
            "kind": "layered_model", "preset": "random", "n_layers": 5,
        })
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("5 layers", r.data["table_text"] + r.summary)
        plt.close("all")

    def test_circular_strike_mean(self):
        import numpy as np
        self.assertAlmostEqual(_circular_strike_mean([10, 10, 10]), 10.0, 1)
        self.assertTrue(np.isnan(_circular_strike_mean([])))

    def test_df_to_text(self):
        import pandas as pd
        df = pd.DataFrame({"a": [1.23456, 2.0], "b": ["x", "y"]})
        txt = _df_to_text(df, columns=["a", "b"], ndigits=2)
        self.assertIn("1.23", txt)
        self.assertIn("x", txt)

    def test_ll_to_utm(self):
        # London ~ UTM zone 30N; easting/northing must be finite & in range.
        import numpy as np
        e, n, z = _ll_to_utm(51.5074, -0.1278, None, "N", "WGS84")
        self.assertEqual(int(z), 30)
        self.assertTrue(np.isfinite(e) and np.isfinite(n))
        self.assertGreater(n, 5e6)  # northing for ~51°N

    def test_get_latlon(self):
        class _Obj:
            lat, lon = 12.5, -3.25
        la, lo = _get_latlon(_Obj())
        self.assertAlmostEqual(la, 12.5)
        self.assertAlmostEqual(lo, -3.25)
        # missing coords → (None, None)
        self.assertEqual(_get_latlon(object()), (None, None))

    def test_safe_filename(self):
        self.assertEqual(_safe_filename("Rho/Phi (xy)"), "Rho_Phi__xy_")
        self.assertEqual(_safe_filename("   "), "figure")

    def test_export_bundles(self):
        # every bundle resolves to known PlotAgent kinds
        from pycsamt.agents.plotting import PLOT_KINDS
        for kinds in _EXPORT_BUNDLES.values():
            for k in kinds:
                self.assertIn(k, PLOT_KINDS)

    def test_unknown_kind(self):
        self.assertEqual(
            ToolAgent().execute({"path": _DATA, "kind": "nope"}).status,
            "failed",
        )

    def test_missing_data(self):
        self.assertEqual(
            ToolAgent().execute({"kind": "strike"}).status, "failed"
        )


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestToolAgent(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def _run(self, kind, **params):
        r = ToolAgent().execute({"path": _DATA, "kind": kind, **params})
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertTrue((r.data or {}).get("table_text"))
        for fig in (r.data.get("figures") or {}).values():
            self.assertIsInstance(fig, plt.Figure)
        return r

    def test_strike(self):
        self._run("strike", method="consensus")

    def test_dimensionality(self):
        self._run("dimensionality")

    def test_validator(self):
        r = self._run("validator")
        # validator returns a table, no figure
        self.assertEqual((r.data.get("figures") or {}), {})

    # ── Wave C: data / IO tools ──────────────────────────────────────────────
    def test_coords(self):
        # No figure; succeeds even when stations lack coordinates.
        r = self._run("coords", datum="WGS84")
        self.assertEqual((r.data.get("figures") or {}), {})

    def test_converter_csv(self):
        import tempfile
        out = tempfile.mkdtemp(prefix="wc_conv_")
        r = self._run("converter", format="csv", output_dir=out)
        self.assertTrue(os.path.isfile(os.path.join(out, "survey_stations.csv")))
        self.assertEqual((r.data.get("figures") or {}), {})

    def test_converter_json(self):
        import tempfile
        out = tempfile.mkdtemp(prefix="wc_conv_")
        self._run("converter", format="json", output_dir=out)
        self.assertTrue(
            os.path.isfile(os.path.join(out, "survey_stations.json"))
        )

    def test_batch_export(self):
        import tempfile
        out = tempfile.mkdtemp(prefix="wc_fig_")
        r = self._run("batch_export", plots="overview",
                      format="png", dpi=80, output_dir=out)
        saved = [f for f in os.listdir(out) if f.endswith(".png")]
        self.assertTrue(saved, msg="no figures written")
        self.assertTrue((r.data.get("figures") or {}))

    def test_elevation_offline(self):
        # Stations without coordinates → graceful, no network call.
        r = ToolAgent().execute({"path": _DATA, "kind": "elevation"})
        self.assertEqual(r.status, "success")

    def test_elevation_fetch_mocked(self):
        # Inject synthetic coords + a stubbed API so no network is touched.
        import pycsamt.agents.tooling as T
        import pycsamt.gis.utils as G
        orig_coords = T._station_coords
        orig_api = G.get_elevation_from_api
        T._station_coords = lambda sites: [
            ("A", 12.0, -3.0), ("B", 12.1, -3.1), ("C", None, None),
        ]
        G.get_elevation_from_api = lambda lat, lon, api_name=None: (
            [100.0, 200.0]
        )
        try:
            r = ToolAgent().execute({"path": _DATA, "kind": "elevation"})
        finally:
            T._station_coords = orig_coords
            G.get_elevation_from_api = orig_api
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("100", r.data["table_text"])
        self.assertIn("200", r.data["table_text"])

    def test_freq_editor(self):
        # Mutating tool: succeeds and hands back an edited survey for the
        # post-processing modal.
        r = ToolAgent().execute({
            "path": _DATA, "kind": "freq_editor",
            "mode": "recover", "method": "composite", "threshold": 0.5,
        })
        self.assertEqual(r.status, "success", msg=r.summary)
        self.assertIn("corrected_sites", r.data)
        self.assertIsNotNone(r.data["corrected_sites"])
        for fig in (r.data.get("figures") or {}).values():
            self.assertIsInstance(fig, plt.Figure)


if __name__ == "__main__":
    unittest.main(verbosity=2)
