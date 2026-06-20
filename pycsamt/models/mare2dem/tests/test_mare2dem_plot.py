# -*- coding: utf-8 -*-
"""Real-data integration tests for models/mare2dem/plot.py.

Data sources
------------
  MT inversion  — data/mare2dem/demo_mt_inversion/
      demo_mt_synth.emdata  · demo.{0,6}.resistivity  · demo.6.resp
      13 MT receivers · 11 frequencies · 6 iterations · RMS ≈ 1.002
      PlotConvergence ✓  PlotSurveyLayout ✓

  CSEM forward  — data/mare2dem/demo_csem/
      demo_csem.emdata  · demo.0.resistivity  · demo.0.resp
      14 CSEM transmitters · 0 MT receivers · no logfile (forward only)
      PlotSurveyLayout ✓  (PlotConvergence not applicable — no iterations)

  Hill model    — data/mare2dem/hill/
      hill.emdata  · hill.0.resistivity  · hill.0.resp
      228 MT receivers · forward problem only
      PlotSurveyLayout ✓

PlotModel and PlotResponse are not yet implemented in the MARE2DEM port
(they raise NotImplementedError) and are therefore excluded from this suite.

Run:
    pytest -v pycsamt/models/mare2dem/tests/test_mare2dem_plot.py
"""
from __future__ import annotations

import pytest
from pathlib import Path

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.figure

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_DATA_ROOT = Path(__file__).parents[4] / "data" / "mare2dem"
_MT_DIR    = _DATA_ROOT / "demo_mt_inversion"
_CSEM_DIR  = _DATA_ROOT / "demo_csem"
_HILL_DIR  = _DATA_ROOT / "hill"

_SKIP_MT   = pytest.mark.skipif(not _MT_DIR.exists(),
                                  reason=f"MARE2DEM MT data not found: {_MT_DIR}")
_SKIP_CSEM = pytest.mark.skipif(not _CSEM_DIR.exists(),
                                  reason=f"MARE2DEM CSEM data not found: {_CSEM_DIR}")
_SKIP_HILL = pytest.mark.skipif(not _HILL_DIR.exists(),
                                  reason=f"MARE2DEM hill data not found: {_HILL_DIR}")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def result_mt():
    from pycsamt.models.mare2dem.results import InversionResult
    return InversionResult(workdir=_MT_DIR)


@pytest.fixture(scope="module")
def result_csem():
    from pycsamt.models.mare2dem.results import InversionResult
    return InversionResult(workdir=_CSEM_DIR)


@pytest.fixture(scope="module")
def result_hill():
    from pycsamt.models.mare2dem.results import InversionResult
    return InversionResult(workdir=_HILL_DIR)


@pytest.fixture(autouse=True)
def close_figs():
    yield
    plt.close("all")


def _is_figure(obj) -> bool:
    return isinstance(obj, matplotlib.figure.Figure)


# ===========================================================================
# GROUP 1 — InversionResult loading
# ===========================================================================

class TestMare2DEMResults:

    @_SKIP_MT
    def test_mt_result_loads(self, result_mt):
        assert result_mt is not None
        assert result_mt.data is not None
        assert result_mt.log is not None

    @_SKIP_MT
    def test_mt_n_iterations(self, result_mt):
        assert result_mt.n_iterations == 6

    @_SKIP_MT
    def test_mt_final_rms(self, result_mt):
        assert result_mt.final_rms == pytest.approx(1.002, abs=0.01)

    @_SKIP_MT
    def test_mt_n_receivers(self, result_mt):
        assert result_mt.data.n_mt_receivers == 13

    @_SKIP_MT
    def test_mt_n_frequencies(self, result_mt):
        assert result_mt.data.n_mt_frequencies == 11

    @_SKIP_MT
    def test_mt_converged(self, result_mt):
        # Demo converges to RMS ≈ 1 within 6 iterations
        assert result_mt.log.converged is True or result_mt.final_rms < 1.1

    @_SKIP_CSEM
    def test_csem_result_loads(self, result_csem):
        assert result_csem is not None
        assert result_csem.data is not None

    @_SKIP_CSEM
    def test_csem_no_iterations(self, result_csem):
        assert result_csem.n_iterations == 0

    @_SKIP_CSEM
    def test_csem_has_transmitters(self, result_csem):
        assert result_csem.data.n_csem_transmitters >= 1

    @_SKIP_HILL
    def test_hill_result_loads(self, result_hill):
        assert result_hill is not None
        assert result_hill.data is not None

    @_SKIP_HILL
    def test_hill_mt_receivers(self, result_hill):
        assert result_hill.data.n_mt_receivers == 228

    @_SKIP_HILL
    def test_hill_no_iterations(self, result_hill):
        assert result_hill.n_iterations == 0


# ===========================================================================
# GROUP 2 — PlotConvergence  (requires logfile with iteration records)
# ===========================================================================

class TestPlotConvergence:

    @_SKIP_MT
    def test_returns_figure(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        assert _is_figure(fig)

    @_SKIP_MT
    def test_figure_has_axes(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        assert len(fig.get_axes()) >= 1

    @_SKIP_MT
    def test_y_label_is_rms(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        ylabel = fig.get_axes()[0].get_ylabel().lower()
        assert "rms" in ylabel or "misfit" in ylabel

    @_SKIP_MT
    def test_curve_has_6_points(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        ax  = fig.get_axes()[0]
        lines = ax.get_lines()
        assert len(lines) >= 1
        assert len(lines[0].get_ydata()) == 6

    @_SKIP_MT
    def test_rms_decreases_overall(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        rms = fig.get_axes()[0].get_lines()[0].get_ydata()
        assert float(rms[0]) > float(rms[-1])

    @_SKIP_MT
    def test_final_rms_near_one(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        rms = fig.get_axes()[0].get_lines()[0].get_ydata()
        assert abs(float(rms[-1]) - 1.0) < 0.1

    @_SKIP_MT
    def test_figure_size_reasonable(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        w, h = fig.get_size_inches()
        assert w >= 4 and h >= 2

    @_SKIP_MT
    def test_x_axis_is_iteration_number(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence
        fig = PlotConvergence(result_mt).plot()
        ax  = fig.get_axes()[0]
        lines = ax.get_lines()
        x = lines[0].get_xdata()
        # x-values should be 1 … n_iterations
        assert int(x[0]) >= 1
        assert int(x[-1]) == result_mt.n_iterations

    def test_no_log_raises(self):
        """PlotConvergence rejects objects that are not a log or InversionResult."""
        from pycsamt.models.mare2dem.plot import PlotConvergence

        with pytest.raises((ValueError, AttributeError, TypeError)):
            PlotConvergence(object()).plot()

    @_SKIP_CSEM
    def test_forward_only_result_raises(self, result_csem):
        """CSEM-only forward dataset has no logfile → should raise."""
        from pycsamt.models.mare2dem.plot import PlotConvergence
        with pytest.raises((ValueError, AttributeError)):
            PlotConvergence(result_csem).plot()


# ===========================================================================
# GROUP 3 — PlotSurveyLayout  (requires EMDataFile)
# ===========================================================================

class TestPlotSurveyLayout:

    @_SKIP_MT
    def test_mt_survey_returns_axes(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_mt.data).plot(ax=ax)
        assert ax is not None

    @_SKIP_MT
    def test_mt_survey_has_xlabel(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_mt.data).plot(ax=ax)
        xlabel = ax.get_xlabel().lower()
        assert "easting" in xlabel or "distance" in xlabel or "x" in xlabel

    @_SKIP_MT
    def test_mt_survey_no_crash_on_replot(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        plotter = PlotSurveyLayout(result_mt.data)
        for _ in range(2):
            fig, ax = plt.subplots()
            plotter.plot(ax=ax)
            plt.close(fig)

    @_SKIP_CSEM
    def test_csem_survey_returns_axes(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_csem.data).plot(ax=ax)
        assert ax is not None

    @_SKIP_CSEM
    def test_csem_survey_xlabel(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_csem.data).plot(ax=ax)
        xlabel = ax.get_xlabel().lower()
        assert "easting" in xlabel or "distance" in xlabel or "x" in xlabel

    @_SKIP_HILL
    def test_hill_survey_returns_axes(self, result_hill):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_hill.data).plot(ax=ax)
        assert ax is not None

    @_SKIP_HILL
    def test_hill_survey_228_receivers_visible(self, result_hill):
        """Hill dataset has 228 MT receivers; the plot must not raise."""
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        fig, ax = plt.subplots(figsize=(12, 4))
        PlotSurveyLayout(result_hill.data).plot(ax=ax)
        assert ax is not None

    def test_survey_accepts_ax_none(self):
        """PlotSurveyLayout must work when ax is not provided (creates its own)."""
        if not _MT_DIR.exists():
            pytest.skip("MT data not available")
        from pycsamt.models.mare2dem.results import InversionResult
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        r = InversionResult(workdir=_MT_DIR)
        PlotSurveyLayout(r.data).plot()   # no ax argument

    def test_survey_invalid_data_raises(self):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout
        with pytest.raises((AttributeError, TypeError, ValueError)):
            PlotSurveyLayout(None).plot()


# ===========================================================================
# GROUP 4 — PlotModel / PlotResponse  (not yet implemented)
# ===========================================================================

class TestNotYetImplemented:

    @_SKIP_MT
    def test_plot_model_raises_not_implemented(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotModel
        with pytest.raises(NotImplementedError):
            fig, ax = plt.subplots()
            PlotModel(result_mt).plot(ax=ax)

    @_SKIP_MT
    def test_plot_response_raises_not_implemented(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotResponse
        with pytest.raises(NotImplementedError):
            fig, ax = plt.subplots()
            PlotResponse(result_mt).plot(ax=ax)


# ===========================================================================
# GROUP 5 — InversionResult API consistency checks
# ===========================================================================

class TestResultAPIs:

    @_SKIP_MT
    def test_log_iteration_records_length(self, result_mt):
        iters = result_mt.log.iterations
        assert len(iters) == result_mt.n_iterations

    @_SKIP_MT
    def test_log_rms_values_positive(self, result_mt):
        for rec in result_mt.log.iterations:
            assert rec.rms > 0

    @_SKIP_MT
    def test_log_rms_monotone_decreasing(self, result_mt):
        """RMS should generally decrease (not strictly, but start > end)."""
        rms_vals = [rec.rms for rec in result_mt.log.iterations]
        assert rms_vals[0] > rms_vals[-1]

    @_SKIP_MT
    def test_data_n_data_positive(self, result_mt):
        assert result_mt.data.n_data > 0

    @_SKIP_CSEM
    def test_csem_data_n_data(self, result_csem):
        assert result_csem.data.n_data > 0

    @_SKIP_HILL
    def test_hill_data_n_data(self, result_hill):
        assert result_hill.data.n_data > 0
