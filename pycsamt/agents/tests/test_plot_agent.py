# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for PlotAgent (rho/phi, phase & phase-tensor pseudo-sections)."""

from __future__ import annotations

import os
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from pycsamt.agents.plotting import (  # noqa: E402
    PLOT_KINDS,
    PlotAgent,
    _as_list,
    _has_tipper,
    _norm_components,
    _norm_parts,
    _period_range,
    _truthy,
)
from pycsamt.emtools._core import ensure_sites  # noqa: E402

_DATA = os.path.join("data", "3edis")
_HAS_DATA = os.path.isdir(_DATA)
_TIP = "edi_out"
_HAS_TIP = (
    os.path.isdir(_TIP)
    and _has_tipper(ensure_sites(_TIP, recursive=True, verbose=0))
    if os.path.isdir(_TIP)
    else False
)


class TestPlotHelpers(unittest.TestCase):
    def test_as_list(self):
        self.assertEqual(_as_list("a, b ,c"), ["a", "b", "c"])
        self.assertEqual(_as_list(["x", "y"]), ["x", "y"])
        self.assertEqual(_as_list(None), [])
        self.assertEqual(_as_list(""), [])

    def test_norm_components(self):
        self.assertEqual(_norm_components(""), ("xy", "yx"))
        self.assertEqual(_norm_components("both"), ("xy", "yx"))
        self.assertEqual(_norm_components("xy"), ("xy",))
        self.assertEqual(_norm_components("xy,yx"), ("xy", "yx"))
        self.assertEqual(_norm_components("det"), ("det",))
        # junk falls back
        self.assertEqual(_norm_components("zz"), ("xy", "yx"))

    def test_norm_parts(self):
        self.assertEqual(_norm_parts(""), ("real", "imag"))
        self.assertEqual(_norm_parts("real"), ("real",))
        self.assertEqual(_norm_parts("imag"), ("imag",))
        self.assertEqual(_norm_parts("real,imag"), ("real", "imag"))

    def test_period_range(self):
        self.assertIsNone(_period_range({}))
        self.assertIsNone(_period_range({"period_min": "", "period_max": ""}))
        self.assertEqual(
            _period_range({"period_min": "0.01", "period_max": "10"}),
            (0.01, 10.0),
        )
        # single-sided bound tolerated
        lo, hi = _period_range({"period_min": "0.1"})
        self.assertEqual(lo, 0.1)
        self.assertGreater(hi, lo)

    def test_truthy(self):
        for v in ("on", "yes", "true", "1", True):
            self.assertTrue(_truthy(v))
        for v in ("off", "no", "false", "0", "", False, None):
            self.assertFalse(_truthy(v))

    def test_unknown_kind_fails(self):
        r = PlotAgent().execute({"path": _DATA, "kind": "nope"})
        self.assertEqual(r.status, "failed")

    def test_missing_data_fails(self):
        r = PlotAgent().execute({"kind": "rhophi"})
        self.assertEqual(r.status, "failed")


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestPlotAgentFigures(unittest.TestCase):
    """Each plot kind produces a real matplotlib figure from sample data."""

    def tearDown(self):
        plt.close("all")

    def _run(self, kind, **params):
        r = PlotAgent().execute({"path": _DATA, "kind": kind, **params})
        self.assertEqual(r.status, "success", msg=r.summary)
        figs = (r.data or {}).get("figures", {})
        self.assertTrue(figs, "no figure produced")
        for fig in figs.values():
            self.assertIsInstance(fig, plt.Figure)
            self.assertTrue(fig.axes)
        return r

    def test_kinds_constant(self):
        self.assertEqual(
            set(PLOT_KINDS),
            {
                "rhophi",
                "phase_psection",
                "pt_psection",
                "tipper",
                "pt_map",
                "station_response",
                "strike_profile",
                "pt_strip",
                "pt_strip_grid",
            },
        )

    def test_pt_map(self):
        self._run("pt_map", period="1", color_by="skew")

    def test_station_response_default_and_multi(self):
        r = self._run("station_response", components="xy,yx")
        self.assertEqual(len((r.data or {}).get("figures", {})), 2)

    def test_strike_profile(self):
        for m in ("consensus", "sweep", "pt"):
            self._run("strike_profile", method=m)

    def test_tipper_unavailable_is_clear(self):
        # the AMT sample (3edis) carries no tipper → clean 'no_tipper' fail
        r = PlotAgent().execute({"path": _DATA, "kind": "tipper"})
        self.assertEqual(r.status, "failed")
        self.assertEqual(r.get("reason"), "no_tipper")

    def test_rhophi(self):
        self._run("rhophi", components="xy,yx", publication="on")

    def test_phase_psection(self):
        self._run("phase_psection", components="xy,yx")

    def test_pt_psection(self):
        self._run("pt_psection", color_by="skew")

    def test_pt_strip_default_station(self):
        r = self._run("pt_strip")
        self.assertTrue(r.warnings, "expected a 'no station given' warning")

    def test_pt_strip_explicit_station(self):
        self._run("pt_strip", stations="new_E1_1")

    def test_pt_strip_grid(self):
        self._run("pt_strip_grid", per_line="2")


@unittest.skipUnless(_HAS_TIP, "no EDI dataset with tipper available")
class TestTipperPlots(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def _run(self, **params):
        r = PlotAgent().execute({"path": _TIP, "kind": "tipper", **params})
        self.assertEqual(r.status, "success", msg=r.summary)
        figs = (r.data or {}).get("figures", {})
        self.assertTrue(figs)
        for fig in figs.values():
            self.assertIsInstance(fig, plt.Figure)
        return r

    def test_components(self):
        self._run(view="components", parts="real,imag")

    def test_arrows(self):
        self._run(view="arrows", convention="park", period="1")


if __name__ == "__main__":
    unittest.main(verbosity=2)
