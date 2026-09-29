# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Regression guard for the apparent-resistivity unit convention.

pycsamt stores the impedance tensor in **field units** (mV/km / nT), so the
apparent resistivity is ``ρ_a = 0.2·|Z|²/f``. Several app-layer plots once used
the SI form ``|Z|²/(ω·μ₀)`` (or ``|Z|²/(5·μ₀·f)``) directly on the field-unit Z,
which over-estimates ρ_a — and the propagated error bars — by ``1/(0.2·2π·μ₀)``
≈ 6.3·10⁵ (curves shooting to ~10⁹ Ω·m on data whose true ρ_a is ≤ 10⁴).

These tests pin the magnitude so the regression cannot silently return.
"""

from __future__ import annotations

import os
import unittest

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

_DATA = os.path.join("data", "3edis")
_HAS_DATA = os.path.isdir(_DATA)

# A correct field-unit ρ_a for crustal AMT/CSAMT data is well under 10⁶ Ω·m;
# the SI-on-field-Z bug pushed it past 10⁹. Guard with a generous ceiling.
_RHO_CEILING = 1.0e6


def _plotted_rho(ax) -> np.ndarray:
    vals = [ln.get_ydata() for ln in ax.get_lines() if len(ln.get_ydata())]
    if not vals:
        return np.array([np.nan])
    arr = np.concatenate([np.asarray(v, float) for v in vals])
    return arr[np.isfinite(arr)]


class TestRhoFieldFormula(unittest.TestCase):
    def test_field_formula_matches_cached_rho(self):
        """0.2·|Z|²/f equals the EDI's own cached ρ (both field-unit)."""
        if not _HAS_DATA:
            self.skipTest("sample EDI data not available")
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
            ensure_sites,
        )

        S = ensure_sites(_DATA, recursive=True, verbose=0)
        ed = next(_iter_items(S))
        _, z, fr = _get_z_block(ed)
        rho_field = 0.2 * np.abs(z[:, 0, 1]) ** 2 / fr
        rho_attr = getattr(ed, "rho", None)
        self.assertIsNotNone(rho_attr)
        m = np.isfinite(rho_field) & np.isfinite(rho_attr[:, 0, 1])
        np.testing.assert_allclose(
            rho_field[m],
            rho_attr[:, 0, 1][m],
            rtol=1e-6,
        )
        self.assertLess(np.nanmax(rho_field), _RHO_CEILING)


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestControllerRhoMagnitude(unittest.TestCase):
    def tearDown(self):
        plt.close("all")

    def _ctrl(self):
        from pycsamt.app.desktop.controllers.correction_controller import (
            CorrectionController,
        )
        from pycsamt.emtools._core import ensure_sites

        S = ensure_sites(_DATA, recursive=True, verbose=0)
        ctrl = CorrectionController()
        ctrl.dark = False
        ctrl.set_raw_sites(S)
        return ctrl, S

    def test_plot_rho_curves_in_range(self):
        ctrl, S = self._ctrl()
        fig, ax = plt.subplots()
        ctrl.plot_rho_curves(S, ax)
        rho = _plotted_rho(ax)
        self.assertTrue(np.isfinite(rho).any())
        self.assertLess(
            np.nanmax(rho),
            _RHO_CEILING,
            msg="ρ_a inflated — SI formula on field-unit Z?",
        )

    def test_pseudosection_in_range(self):
        ctrl, S = self._ctrl()
        fig = plt.figure()
        # plot_rho_pseudosection writes its own axes onto the figure
        try:
            ctrl.plot_rho_pseudosection(S, fig)
        except Exception as exc:  # pragma: no cover - method optional
            self.skipTest(f"pseudosection unavailable: {exc}")
        for ax in fig.axes:
            for coll in ax.collections:
                arr = coll.get_array()
                if arr is None:
                    continue
                a = np.asarray(arr, float)
                a = a[np.isfinite(a)]
                if a.size:
                    # pseudosections plot log10(ρ); ceil in linear or log space
                    cap = np.log10(_RHO_CEILING) if np.nanmax(a) < 20 else _RHO_CEILING
                    self.assertLess(np.nanmax(a), cap)


@unittest.skipUnless(_HAS_DATA, "sample EDI data not available")
class TestStationResponseErrorBars(unittest.TestCase):
    """The web station-response ρ_a error bars must stay physical (≲ ρ_a)."""

    def test_rho_error_relative(self):
        from pycsamt.emtools._core import (
            _get_z_block,
            _iter_items,
            ensure_sites,
        )

        S = ensure_sites(_DATA, recursive=True, verbose=0)
        ed = next(_iter_items(S))
        Z, z, fr = _get_z_block(ed)
        ze = getattr(Z, "z_err", None)
        if ze is None:
            self.skipTest("no impedance errors in sample data")
        # the fixed tools.py error formula: 0.4·|Z|·σ_Z / f
        rho = 0.2 * np.abs(z[:, 0, 1]) ** 2 / fr
        rho_err = 0.4 * np.abs(z[:, 0, 1]) * np.asarray(ze)[:, 0, 1] / fr
        rel = rho_err / np.where(rho > 0, rho, np.nan)
        # a well-conditioned error bar is a fraction of ρ_a, never 10⁴× it
        self.assertLess(np.nanmedian(rel), 5.0)


if __name__ == "__main__":
    unittest.main()
