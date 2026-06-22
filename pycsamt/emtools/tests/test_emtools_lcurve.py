# -*- coding: utf-8 -*-
"""Tests for pycsamt.emtools.lcurve"""
from __future__ import annotations

import numpy as np
import pytest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.lcurve import lcurve_table, plot_lcurve


# ─────────────────────────────────────────────────────────────────────────────
# Synthetic L-curve data
# ─────────────────────────────────────────────────────────────────────────────

def _ldata(n: int = 20):
    """Return (misfit, roughness, lambda) arrays describing a clean L-curve."""
    lam = np.logspace(3, -3, n)
    rough  = lam ** -1.2   # roughness decreases with λ
    misfit = lam **  0.8   # misfit decreases as λ decreases
    return misfit, rough, lam


# ─────────────────────────────────────────────────────────────────────────────
# lcurve_table
# ─────────────────────────────────────────────────────────────────────────────

class TestLcurveTable:

    def test_returns_dataframe(self):
        import pandas as pd
        m, r, l = _ldata()
        df = lcurve_table(m, r, l)
        assert isinstance(df, pd.DataFrame)

    def test_expected_columns(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l)
        for col in ("misfit", "rough", "lam", "curv"):
            assert col in df.columns

    def test_row_count(self):
        n = 20
        m, r, l = _ldata(n)
        df = lcurve_table(m, r, l)
        assert len(df) == n

    def test_corner_idx_attr(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l)
        assert "corner_idx" in df.attrs
        assert 0 <= df.attrs["corner_idx"] < len(df)

    def test_no_lambda(self):
        """Lambda is optional."""
        m, r, _ = _ldata()
        df = lcurve_table(m, r)
        assert isinstance(df, df.__class__)

    def test_return_dict(self):
        m, r, l = _ldata()
        d = lcurve_table(m, r, l, return_dict=True)
        assert isinstance(d, dict)
        assert "corner" in d
        assert "rough" in d

    def test_curvature_method(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l, method="curvature")
        assert not df.empty

    def test_maxdist_method(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l, method="maxdist")
        assert not df.empty

    def test_positive_values_only(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l)
        assert (df["misfit"].dropna() > 0).all()
        assert (df["rough"].dropna() > 0).all()

    def test_slope_column_present(self):
        m, r, l = _ldata()
        df = lcurve_table(m, r, l)
        assert "slope" in df.columns


# ─────────────────────────────────────────────────────────────────────────────
# plot_lcurve
# ─────────────────────────────────────────────────────────────────────────────

class TestPlotLcurve:

    def test_returns_axes(self):
        m, r, l = _ldata()
        ax = plot_lcurve(m, r, l)
        plt.close("all")
        assert ax is not None

    def test_external_ax(self):
        m, r, l = _ldata()
        fig, ax_in = plt.subplots()
        ax = plot_lcurve(m, r, l, ax=ax_in)
        plt.close("all")
        assert ax is not None

    def test_multiple_curves(self):
        m1, r1, l1 = _ldata(15)
        m2, r2, l2 = _ldata(20)
        ax = plot_lcurve([m1, m2], [r1, r2], [l1, l2],
                         labels=["run-A", "run-B"])
        plt.close("all")
        assert ax is not None

    def test_no_corner_marker(self):
        m, r, l = _ldata()
        ax = plot_lcurve(m, r, l, show_corner=False)
        plt.close("all")
        assert ax is not None

    def test_no_inset(self):
        m, r, l = _ldata()
        ax = plot_lcurve(m, r, l, show_inset=False)
        plt.close("all")
        assert ax is not None

    def test_no_lambda_no_crash(self):
        m, r, _ = _ldata()
        ax = plot_lcurve(m, r)
        plt.close("all")
        assert ax is not None
