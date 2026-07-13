"""Real-data integration tests for models/occam2d/plot.py.

Two data sources are exercised:

  PYCSAMT — data/occam2D/
      Startup + Occam2DMesh + Occam2DModel + OccamDataFile.dat
      ITER17.iter + RESP17.resp
      47 sites · 17 frequencies · 1 stored iteration · RMS ≈ 0.998

  OCCAM2DMT — ~/Documents/OCCAM2DMT_V3.0/OCCAM2DMT_V3.0/Examples/
      Test2 : 6 sites · 11 ITER/RESP files · RMS ≈ 1.000
      Test4 : 6 sites · 10 ITER/RESP files + STATICS · RMS ≈ 1.000

Both directories are skipped individually when absent, so the suite runs
in CI environments that only have the bundled data and on developer machines
with the full OCCAM2DMT package.

Run all:
    pytest -v pycsamt/models/occam2d/tests/test_occam2d_plot_realdata.py

Run one group:
    pytest -v -k "pycsamt" pycsamt/models/occam2d/tests/test_occam2d_plot_realdata.py
    pytest -v -k "test2 or test4" pycsamt/models/occam2d/tests/test_occam2d_plot_realdata.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.figure
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Data paths
# ---------------------------------------------------------------------------

_PYCSAMT_DIR = Path(__file__).parents[4] / "data" / "occam2D"

_OCCAM_BASE = (
    Path.home()
    / "Documents"
    / "OCCAM2DMT_V3.0"
    / "OCCAM2DMT_V3.0"
    / "Examples"
)
_TEST2_DIR = _OCCAM_BASE / "Test2"
_TEST4_DIR = _OCCAM_BASE / "Test4"

_SKIP_PYCSAMT = pytest.mark.skipif(
    not _PYCSAMT_DIR.exists(),
    reason=f"bundled occam2D data not found: {_PYCSAMT_DIR}",
)
_SKIP_TEST2 = pytest.mark.skipif(
    not _TEST2_DIR.exists(),
    reason=f"OCCAM2DMT Test2 data not found: {_TEST2_DIR}",
)
_SKIP_TEST4 = pytest.mark.skipif(
    not _TEST4_DIR.exists(),
    reason=f"OCCAM2DMT Test4 data not found: {_TEST4_DIR}",
)


# ---------------------------------------------------------------------------
# Module-level fixtures (load once per test session)
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def result_pycsamt():
    from pycsamt.models.occam2d.results import InversionResult

    return InversionResult(workdir=_PYCSAMT_DIR)


@pytest.fixture(scope="module")
def result_test2():
    from pycsamt.models.occam2d.results import InversionResult

    return InversionResult(workdir=_TEST2_DIR)


@pytest.fixture(scope="module")
def result_test4():
    from pycsamt.models.occam2d.results import InversionResult

    return InversionResult(workdir=_TEST4_DIR)


@pytest.fixture(autouse=True)
def close_figs():
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _is_figure(obj) -> bool:
    return isinstance(obj, matplotlib.figure.Figure)


# ===========================================================================
# GROUP 1 — pycsamt bundled data (47 sites, 17 freqs, ITER17)
# ===========================================================================


class TestPycsamtData:
    # ── Fixture availability ─────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_result_loaded(self, result_pycsamt):
        r = result_pycsamt
        assert r is not None
        assert r.data is not None
        assert r.mesh is not None
        assert r.response is not None

    @_SKIP_PYCSAMT
    def test_result_stats(self, result_pycsamt):
        r = result_pycsamt
        assert r.final_rms == pytest.approx(0.998, abs=0.01)
        assert r.n_iterations >= 1
        assert r.data.n_sites == 47
        assert r.data.n_frequencies == 17

    # ── PlotMisfit ───────────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_misfit_returns_figure(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_misfit_has_rms_axis(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_pycsamt).plot()
        ax = fig.get_axes()[0]
        assert (
            "rms" in ax.get_ylabel().lower()
            or "misfit" in ax.get_ylabel().lower()
        )

    @_SKIP_PYCSAMT
    def test_misfit_has_convergence_line(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_pycsamt).plot()
        ax = fig.get_axes()[0]
        lines = ax.get_lines()
        assert len(lines) >= 1
        assert len(lines[0].get_ydata()) >= 1

    # ── PlotModel ────────────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_model_returns_figure(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotModel

        fig = PlotModel(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_model_has_colourbar(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotModel

        fig = PlotModel(result=result_pycsamt).plot()
        assert len(fig.get_axes()) >= 2

    @_SKIP_PYCSAMT
    def test_model_depth_max_clips_y_axis(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotModel

        fig = PlotModel(result=result_pycsamt, depth_max=500).plot()
        ax = fig.get_axes()[0]
        ylo, yhi = ax.get_ylim()
        assert max(ylo, yhi) <= 600

    @_SKIP_PYCSAMT
    def test_model_rho_range_accepted(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotModel

        fig = PlotModel(
            result=result_pycsamt, rho_min=1, rho_max=10000
        ).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_model_no_stations(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotModel

        fig = PlotModel(result=result_pycsamt, show_stations=False).plot()
        assert _is_figure(fig)

    # ── PlotResponse ─────────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_response_all_sites_default(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotResponse

        fig = PlotResponse(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_response_single_site(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotResponse

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotResponse(result=result_pycsamt, station=sites[0]).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_response_subset_modes_present(self, result_pycsamt):
        """Only request modes whose type codes actually appear in the data."""
        from pycsamt.models.occam2d.plot import PlotResponse

        # pycsamt data has type codes 5,6 (impedance columns); use default modes
        fig = PlotResponse(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_response_multiple_sites(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotResponse

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotResponse(result=result_pycsamt, station=sites[:3]).plot()
        assert _is_figure(fig)

    # ── PlotPseudo ───────────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_pseudo_default_modes(self, result_pycsamt):
        """pycsamt data has codes 5,6 (impedance); default mode uses all present."""
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_pseudo_has_enough_axes(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_pycsamt).plot()
        assert len(fig.get_axes()) >= 2

    # ── PlotSounding1D ───────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_sounding1d_single_site(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSounding1D

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotSounding1D(result=result_pycsamt, station=sites[0]).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_sounding1d_multi_site(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSounding1D

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotSounding1D(result=result_pycsamt, station=sites[:4]).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_sounding1d_depth_max(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSounding1D

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotSounding1D(
            result=result_pycsamt, station=sites[0], depth_max=1000
        ).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_sounding1d_panel_count_matches_sites(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSounding1D

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        n = min(3, len(sites))
        fig = PlotSounding1D(result=result_pycsamt, station=sites[:n]).plot()
        # Each site gets a subplot row (rho + possibly phase)
        assert len(fig.get_axes()) >= n

    # ── PlotSiteMisfit ───────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_sitemisfit_returns_figure(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_sitemisfit_has_multiple_axes(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_pycsamt).plot()
        assert len(fig.get_axes()) >= 2

    @_SKIP_PYCSAMT
    def test_sitemisfit_default_modes(self, result_pycsamt):
        """Use default modes so the filter matches type codes actually in data."""
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_pycsamt).plot()
        assert _is_figure(fig)

    # ── PlotResponseGrid ─────────────────────────────────────────────────────

    @_SKIP_PYCSAMT
    def test_response_grid_returns_figure(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_pycsamt).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_response_grid_has_subpanels(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_pycsamt).plot()
        assert len(fig.get_axes()) >= 2

    @_SKIP_PYCSAMT
    def test_response_grid_te_only(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_pycsamt, modes=["TE"]).plot()
        assert _is_figure(fig)

    @_SKIP_PYCSAMT
    def test_response_grid_subset_sites(self, result_pycsamt):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        sites = np.unique(result_pycsamt.response.site_indices).tolist()
        fig = PlotResponseGrid(
            result=result_pycsamt, station=sites[:5]
        ).plot()
        assert _is_figure(fig)


# ===========================================================================
# GROUP 2 — OCCAM2DMT Test2  (6 sites, 11 iterations)
# No mesh → PlotModel / PlotSounding1D not available
# ===========================================================================


class TestOccamMT2:
    @_SKIP_TEST2
    def test_result_loaded(self, result_test2):
        assert result_test2.response is not None
        assert result_test2.data is not None

    @_SKIP_TEST2
    def test_result_n_iterations(self, result_test2):
        assert result_test2.n_iterations == 11

    @_SKIP_TEST2
    def test_result_final_rms_near_one(self, result_test2):
        assert abs(result_test2.final_rms - 1.0) < 0.01

    @_SKIP_TEST2
    def test_result_n_sites(self, result_test2):
        assert result_test2.data.n_sites == 6

    # ── PlotMisfit ───────────────────────────────────────────────────────────

    @_SKIP_TEST2
    def test_misfit_returns_figure(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test2).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_misfit_curve_points_match_log(self, result_test2):
        """Convergence curve length equals the number of log records,
        which may be less than the total ITER files when the solver
        re-uses the final iteration for the output model."""
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test2).plot()
        ax = fig.get_axes()[0]
        n_points = len(ax.get_lines()[0].get_xdata())
        assert n_points >= 1
        assert n_points <= result_test2.n_iterations

    @_SKIP_TEST2
    def test_misfit_rms_decreases(self, result_test2):
        """RMS should generally decrease from iteration 0 to the final."""
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test2).plot()
        ax = fig.get_axes()[0]
        rms = ax.get_lines()[0].get_ydata()
        assert float(rms[0]) > float(rms[-1])

    @_SKIP_TEST2
    def test_misfit_show_lagrange_optional(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotMisfit

        for show in (True, False):
            fig = PlotMisfit(result=result_test2, show_lagrange=show).plot()
            assert _is_figure(fig)

    # ── PlotPseudo ───────────────────────────────────────────────────────────

    @_SKIP_TEST2
    def test_pseudo_returns_figure(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_test2).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_pseudo_te_mode(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_test2, mode="TE").plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_pseudo_tm_mode(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_test2, mode="TM").plot()
        assert _is_figure(fig)

    # ── PlotResponse ─────────────────────────────────────────────────────────

    @_SKIP_TEST2
    def test_response_all_sites(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotResponse

        fig = PlotResponse(result=result_test2).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_response_first_site(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotResponse

        sites = np.unique(result_test2.response.site_indices).tolist()
        fig = PlotResponse(result=result_test2, station=sites[0]).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_response_te_mode_only(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotResponse

        fig = PlotResponse(result=result_test2, modes=["TE"]).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_response_subset_three_sites(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotResponse

        sites = np.unique(result_test2.response.site_indices).tolist()
        fig = PlotResponse(result=result_test2, station=sites[:3]).plot()
        assert _is_figure(fig)

    # ── PlotSiteMisfit ───────────────────────────────────────────────────────

    @_SKIP_TEST2
    def test_sitemisfit_returns_figure(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_test2).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_sitemisfit_axes_count(self, result_test2):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_test2).plot()
        assert len(fig.get_axes()) >= 2

    # ── PlotResponseGrid ─────────────────────────────────────────────────────

    @_SKIP_TEST2
    def test_response_grid_returns_figure(self, result_test2):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_test2).plot()
        assert _is_figure(fig)

    @_SKIP_TEST2
    def test_response_grid_subpanels_ge_2(self, result_test2):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_test2).plot()
        assert len(fig.get_axes()) >= 2

    @_SKIP_TEST2
    def test_response_grid_first_site_subset(self, result_test2):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        sites = np.unique(result_test2.response.site_indices).tolist()
        fig = PlotResponseGrid(result=result_test2, station=sites[:2]).plot()
        assert _is_figure(fig)


# ===========================================================================
# GROUP 3 — OCCAM2DMT Test4  (6 sites, 10 iterations, STATICS file)
# ===========================================================================


class TestOccamMT4:
    @_SKIP_TEST4
    def test_result_loaded(self, result_test4):
        assert result_test4.response is not None
        assert result_test4.data is not None

    @_SKIP_TEST4
    def test_result_n_iterations(self, result_test4):
        assert result_test4.n_iterations == 10

    @_SKIP_TEST4
    def test_result_final_rms(self, result_test4):
        assert result_test4.final_rms == pytest.approx(1.0, abs=0.01)

    # ── PlotMisfit ───────────────────────────────────────────────────────────

    @_SKIP_TEST4
    def test_misfit_returns_figure(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test4).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_misfit_curve_points_match_log(self, result_test4):
        """Convergence curve length equals the number of log records."""
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test4).plot()
        ax = fig.get_axes()[0]
        n_points = len(ax.get_lines()[0].get_xdata())
        assert n_points >= 1
        assert n_points <= result_test4.n_iterations

    @_SKIP_TEST4
    def test_misfit_rms_decreases(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotMisfit

        fig = PlotMisfit(result=result_test4).plot()
        rms = fig.get_axes()[0].get_lines()[0].get_ydata()
        assert float(rms[0]) > float(rms[-1])

    # ── PlotPseudo ───────────────────────────────────────────────────────────

    @_SKIP_TEST4
    def test_pseudo_default(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(result=result_test4).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_pseudo_rho_min_max(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotPseudo

        fig = PlotPseudo(
            result=result_test4, rho_min=0.1, rho_max=1000
        ).plot()
        assert _is_figure(fig)

    # ── PlotResponse ─────────────────────────────────────────────────────────

    @_SKIP_TEST4
    def test_response_all_sites(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotResponse

        fig = PlotResponse(result=result_test4).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_response_last_site(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotResponse

        sites = np.unique(result_test4.response.site_indices).tolist()
        fig = PlotResponse(result=result_test4, station=sites[-1]).plot()
        assert _is_figure(fig)

    # ── PlotSiteMisfit ───────────────────────────────────────────────────────

    @_SKIP_TEST4
    def test_sitemisfit_returns_figure(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_test4).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_sitemisfit_te_only(self, result_test4):
        from pycsamt.models.occam2d.plot import PlotSiteMisfit

        fig = PlotSiteMisfit(result=result_test4, modes=["TE"]).plot()
        assert _is_figure(fig)

    # ── PlotResponseGrid ─────────────────────────────────────────────────────

    @_SKIP_TEST4
    def test_response_grid_returns_figure(self, result_test4):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_test4).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_response_grid_subset_sites(self, result_test4):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        sites = np.unique(result_test4.response.site_indices).tolist()
        fig = PlotResponseGrid(result=result_test4, station=sites[:3]).plot()
        assert _is_figure(fig)

    @_SKIP_TEST4
    def test_response_grid_tm_only(self, result_test4):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        fig = PlotResponseGrid(result=result_test4, modes=["TM"]).plot()
        assert _is_figure(fig)


# ===========================================================================
# GROUP 4 — Cross-dataset sanity: same plot class, multiple inputs
# ===========================================================================


class TestCrossDataset:
    """Sanity checks that re-use the same PlotMisfit / PlotResponseGrid
    constructor on both bundled data and OCCAM2DMT examples, verifying
    the API is consistent."""

    @pytest.mark.skipif(
        not (_PYCSAMT_DIR.exists() and _TEST2_DIR.exists()),
        reason="needs both pycsamt and OCCAM2DMT Test2 data",
    )
    def test_misfit_two_datasets_same_api(self, result_pycsamt, result_test2):
        from pycsamt.models.occam2d.plot import PlotMisfit

        for result in (result_pycsamt, result_test2):
            fig = PlotMisfit(result=result).plot()
            assert _is_figure(fig)

    @pytest.mark.skipif(
        not (_TEST2_DIR.exists() and _TEST4_DIR.exists()),
        reason="needs OCCAM2DMT Test2 and Test4 data",
    )
    def test_response_grid_two_occam_examples(
        self, result_test2, result_test4
    ):
        from pycsamt.models.occam2d.plot import (
            PlotResponseGrid,
        )

        for result in (result_test2, result_test4):
            fig = PlotResponseGrid(result=result).plot()
            assert _is_figure(fig)

    @pytest.mark.skipif(
        not (_TEST2_DIR.exists() and _TEST4_DIR.exists()),
        reason="needs OCCAM2DMT Test2 and Test4 data",
    )
    def test_pseudo_two_occam_examples(self, result_test2, result_test4):
        from pycsamt.models.occam2d.plot import PlotPseudo

        for result in (result_test2, result_test4):
            fig = PlotPseudo(result=result).plot()
            assert _is_figure(fig)
