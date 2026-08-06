"""Tests for pycsamt.emtools.lcurve"""

from __future__ import annotations

from pathlib import Path

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.emtools.lcurve import (
    LCurveData,
    lcurve_from_mare2dem,
    lcurve_from_modem,
    lcurve_from_occam2d,
    lcurve_table,
    plot_lcurve,
)

_REPO_ROOT = Path(__file__).resolve().parents[3]
_OCCAM2D_LOG = _REPO_ROOT / "data" / "occam2D" / "LogFile.logfile"
_MODEM_LOG = (
    _REPO_ROOT
    / "data"
    / "modem"
    / "willy_27freq_watex_line02_sample"
    / "Modular_NLCG.log"
)
_MARE2DEM_LOG = (
    _REPO_ROOT / "data" / "mare2dem" / "demo_mt_inversion" / "demo.logfile"
)

# ─────────────────────────────────────────────────────────────────────────────
# Synthetic L-curve data
# ─────────────────────────────────────────────────────────────────────────────


def _ldata(n: int = 20):
    """Return (misfit, roughness, lambda) arrays describing a clean L-curve."""
    lam = np.logspace(3, -3, n)
    rough = lam**-1.2  # roughness decreases with λ
    misfit = lam**0.8  # misfit decreases as λ decreases
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
        ax = plot_lcurve([m1, m2], [r1, r2], [l1, l2], labels=["run-A", "run-B"])
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

    def test_target_misfit_line(self):
        m, r, l = _ldata()
        ax = plot_lcurve(m, r, l, target_misfit=1.0, target_label="target")
        plt.close("all")
        assert any(
            np.allclose(line.get_ydata(), 1.0) for line in ax.get_lines()
        )

    def test_label_every(self):
        m, r, l = _ldata()
        ax = plot_lcurve(m, r, l, label_every=3, label_prefix="lambda=")
        plt.close("all")
        assert any(
            t.get_text().startswith("lambda=") for t in ax.texts
        )


# ─────────────────────────────────────────────────────────────────────────────
# Real-inversion adapters (Occam2D, ModEM, MARE2DEM)
# ─────────────────────────────────────────────────────────────────────────────


@pytest.mark.skipif(
    not _OCCAM2D_LOG.exists(), reason="bundled Occam2D log not available"
)
class TestLcurveFromOccam2d:
    def test_returns_lcurvedata(self):
        sweep = lcurve_from_occam2d(_OCCAM2D_LOG)
        assert isinstance(sweep, LCurveData)
        assert sweep.backend == "occam2d"

    def test_arrays_aligned_and_positive(self):
        sweep = lcurve_from_occam2d(_OCCAM2D_LOG)
        n = sweep.misfit.size
        assert n > 0
        assert sweep.rough.size == n
        assert sweep.lam.size == n
        assert sweep.iterations.size == n
        assert np.all(sweep.misfit > 0)
        assert np.all(sweep.rough > 0)

    def test_misfit_decreases_overall(self):
        sweep = lcurve_from_occam2d(_OCCAM2D_LOG)
        order = np.argsort(sweep.iterations)
        assert sweep.misfit[order][-1] < sweep.misfit[order][0]

    def test_table_and_plot(self):
        sweep = lcurve_from_occam2d(_OCCAM2D_LOG)
        df = sweep.table()
        assert not df.empty
        ax = sweep.plot()
        plt.close("all")
        assert ax is not None


@pytest.mark.skipif(
    not _MODEM_LOG.exists(), reason="bundled ModEM log not available"
)
class TestLcurveFromModem:
    def test_returns_lcurvedata(self):
        sweep = lcurve_from_modem(_MODEM_LOG)
        assert isinstance(sweep, LCurveData)
        assert sweep.backend == "modem"

    def test_arrays_aligned_and_positive(self):
        sweep = lcurve_from_modem(_MODEM_LOG)
        n = sweep.misfit.size
        assert n > 0
        assert sweep.rough.size == n
        assert sweep.lam.size == n
        assert np.all(sweep.misfit > 0)
        assert np.all(sweep.rough > 0)

    def test_lambda_decreases_with_damping_schedule(self):
        sweep = lcurve_from_modem(_MODEM_LOG)
        order = np.argsort(sweep.iterations)
        assert sweep.lam[order][-1] < sweep.lam[order][0]

    def test_table_and_plot(self):
        sweep = lcurve_from_modem(_MODEM_LOG)
        df = sweep.table()
        assert not df.empty
        ax = sweep.plot()
        plt.close("all")
        assert ax is not None


@pytest.mark.skipif(
    not _MARE2DEM_LOG.exists(), reason="bundled MARE2DEM log not available"
)
class TestLcurveFromMare2dem:
    def test_returns_lcurvedata(self):
        sweep = lcurve_from_mare2dem(_MARE2DEM_LOG)
        assert isinstance(sweep, LCurveData)
        assert sweep.backend == "mare2dem"

    def test_arrays_aligned_and_positive(self):
        sweep = lcurve_from_mare2dem(_MARE2DEM_LOG)
        n = sweep.misfit.size
        assert n > 0
        assert sweep.rough.size == n
        assert sweep.lam.size == n
        assert np.all(sweep.misfit > 0)
        assert np.all(sweep.rough > 0)
        # "Optimal Mu" is log10(mu) in the raw log; the adapter must
        # convert it back to linear scale.
        assert np.all(sweep.lam > 0)

    def test_table_and_plot(self):
        sweep = lcurve_from_mare2dem(_MARE2DEM_LOG)
        df = sweep.table()
        assert not df.empty
        ax = sweep.plot()
        plt.close("all")
        assert ax is not None


class TestLcurveFromMissingFile:
    def test_occam2d_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            lcurve_from_occam2d(tmp_path / "does_not_exist.logfile")

    def test_modem_missing_file_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            lcurve_from_modem(tmp_path / "does_not_exist.log")

    def test_mare2dem_missing_file_returns_empty(self, tmp_path):
        # Mare2DEMLog does not raise for a missing path; it silently
        # leaves `iterations` empty, so the adapter mirrors that instead
        # of pretending it validated the file.
        sweep = lcurve_from_mare2dem(tmp_path / "does_not_exist.logfile")
        assert isinstance(sweep, LCurveData)
        assert sweep.misfit.size == 0
