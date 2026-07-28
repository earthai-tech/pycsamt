"""Tests for Phase-5 plot classes (PlotMisfit, PlotModel, PlotResponse, PlotPseudo)."""

import matplotlib
import pytest

matplotlib.use("Agg")  # non-interactive backend — safe for CI
from pathlib import Path

import matplotlib.pyplot as plt

DATA_DIR = Path(__file__).parents[4] / "data" / "occam2D"

pytestmark = pytest.mark.skipif(
    not DATA_DIR.exists(), reason="sample occam2D data not available"
)


@pytest.fixture(scope="module")
def result():
    from pycsamt.models.occam2d.results import InversionResult

    return InversionResult(workdir=DATA_DIR)


@pytest.fixture(autouse=True)
def close_figs():
    """Close every figure created during a test to avoid resource warnings."""
    yield
    plt.close("all")


# ---------------------------------------------------------------------------
# InversionResult.data loaded
# ---------------------------------------------------------------------------


def test_result_data_loaded(result):
    assert result.data is not None


def test_result_data_n_sites(result):
    assert result.data.n_sites == 47


def test_result_data_n_frequencies(result):
    assert result.data.n_frequencies == 17


# ---------------------------------------------------------------------------
# PlotMisfit
# ---------------------------------------------------------------------------


def test_misfit_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotMisfit

    fig = PlotMisfit(result=result).plot()
    assert hasattr(fig, "savefig")


def test_misfit_with_lagrange(result):
    from pycsamt.models.occam2d.plot import PlotMisfit

    fig = PlotMisfit(result=result, show_lagrange=True).plot()
    assert len(fig.axes) >= 2


def test_misfit_no_roughness(result):
    from pycsamt.models.occam2d.plot import PlotMisfit

    fig = PlotMisfit(result=result, show_roughness=False, target_line=False).plot()
    assert fig is not None


def test_misfit_no_log_raises():
    from pycsamt.models.occam2d.plot import PlotMisfit

    # Construct a minimal result with no log
    class _FakeResult:
        log = None

    with pytest.raises(RuntimeError, match="no log data"):
        PlotMisfit(result=_FakeResult()).plot()


# ---------------------------------------------------------------------------
# PlotModel
# ---------------------------------------------------------------------------


def test_model_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotModel

    fig = PlotModel(result=result).plot()
    assert hasattr(fig, "savefig")


def test_model_dynamic_section_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotModel

    fig = PlotModel(result=result, section="dynamic").plot()
    assert hasattr(fig, "savefig")
    assert len(fig.axes) >= 2


def test_model_no_stations(result):
    from pycsamt.models.occam2d.plot import PlotModel

    fig = PlotModel(result=result, show_stations=False).plot()
    assert fig is not None


def test_model_depth_max(result):
    from pycsamt.models.occam2d.plot import PlotModel

    fig = PlotModel(result=result, depth_max=50000.0).plot()
    assert fig is not None


def test_model_unit_m(result):
    from pycsamt.models.occam2d.plot import PlotModel

    fig = PlotModel(result=result, profile_distance_unit="m").plot()
    assert fig is not None


def test_model_no_rho2d_raises():
    from pycsamt.models.occam2d.plot import PlotModel

    class _FakeResult:
        rho_2d = None

    with pytest.raises(RuntimeError, match="no rho_2d"):
        PlotModel(result=_FakeResult()).plot()


def test_extract_profile(result):
    from pycsamt.models.occam2d.plot import PlotModel

    pm = PlotModel(result=result, profile_distance_unit="km")
    x, z, rho = pm.extract_profile(-50.0, 50.0)
    assert x.shape[0] > 0
    assert rho.shape == (result.mesh.n_zcells, x.shape[0])


# ---------------------------------------------------------------------------
# PlotResponse
# ---------------------------------------------------------------------------


def test_response_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotResponse

    fig = PlotResponse(result=result, max_stations=3).plot()
    assert hasattr(fig, "savefig")


def test_response_tm_only(result):
    from pycsamt.models.occam2d.plot import PlotResponse

    fig = PlotResponse(result=result, modes=["TM"], max_stations=3).plot()
    assert fig is not None


def test_response_selected_stations(result):
    from pycsamt.models.occam2d.plot import PlotResponse

    stations = result.data.sites[:3]
    fig = PlotResponse(result=result, stations=stations).plot()
    assert fig is not None


def test_response_no_data_raises():
    from pycsamt.models.occam2d.plot import PlotResponse

    class _FakeResult:
        response = None
        data = None

    with pytest.raises(RuntimeError, match="no response data"):
        PlotResponse(result=_FakeResult()).plot()


# ---------------------------------------------------------------------------
# PlotPseudo
# ---------------------------------------------------------------------------


def test_pseudo_tm_rho_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotPseudo

    fig = PlotPseudo(result=result, mode="TM", data_type="rho").plot()
    assert hasattr(fig, "savefig")


def test_pseudo_tm_phase(result):
    from pycsamt.models.occam2d.plot import PlotPseudo

    fig = PlotPseudo(result=result, mode="TM", data_type="phase").plot()
    assert fig is not None


def test_pseudo_bad_mode_raises(result):
    from pycsamt.models.occam2d.plot import PlotPseudo

    with pytest.raises((ValueError, RuntimeError)):
        PlotPseudo(result=result, mode="TE", data_type="rho").plot()


def test_pseudo_no_data_raises():
    from pycsamt.models.occam2d.plot import PlotPseudo

    class _FakeResult:
        data = None

    with pytest.raises(RuntimeError, match="no data blocks"):
        PlotPseudo(result=_FakeResult()).plot()


# ---------------------------------------------------------------------------
# PlotSounding1D
# ---------------------------------------------------------------------------


def test_sounding1d_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    fig = PlotSounding1D(result=result, max_stations=4).plot()
    assert hasattr(fig, "savefig")


def test_sounding1d_overlay(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    fig = PlotSounding1D(result=result, overlay=True, max_stations=5).plot()
    assert len(fig.axes) == 1


def test_sounding1d_selected_stations(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    stations = result.data.sites[:3]
    fig = PlotSounding1D(result=result, stations=list(stations)).plot()
    assert fig is not None


def test_sounding1d_depth_max(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    fig = PlotSounding1D(result=result, depth_max=30000.0, max_stations=4).plot()
    assert fig is not None


def test_sounding1d_no_rho2d_raises():
    from pycsamt.models.occam2d.plot import PlotSounding1D

    class _FakeResult:
        rho_2d = None

    with pytest.raises(RuntimeError, match="no rho_2d"):
        PlotSounding1D(result=_FakeResult()).plot()


def test_sounding1d_no_offsets_raises(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    class _FakeResult:
        rho_2d = result.rho_2d
        mesh = result.mesh
        best_iter = None
        data = None

    with pytest.raises(RuntimeError, match="no station offsets"):
        PlotSounding1D(result=_FakeResult()).plot()


def test_sounding1d_panel_count(result):
    from pycsamt.models.occam2d.plot import PlotSounding1D

    stations = result.data.sites[:6]
    fig = PlotSounding1D(result=result, stations=list(stations)).plot()
    # 6 stations → 2 rows × 3 cols (or similar) → visible axes = 6
    visible = [ax for ax in fig.axes if ax.get_visible()]
    assert len(visible) == 6


# ---------------------------------------------------------------------------
# PlotSiteMisfit
# ---------------------------------------------------------------------------


def test_site_misfit_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotSiteMisfit

    fig = PlotSiteMisfit(result=result).plot()
    assert hasattr(fig, "savefig")


def test_site_misfit_two_panels(result):
    from pycsamt.models.occam2d.plot import PlotSiteMisfit

    fig = PlotSiteMisfit(result=result, show_residual_map=True).plot()
    assert len(fig.axes) >= 2


def test_site_misfit_bar_only(result):
    from pycsamt.models.occam2d.plot import PlotSiteMisfit

    fig = PlotSiteMisfit(result=result, show_residual_map=False).plot()
    assert len(fig.axes) == 1


def test_site_misfit_tm_only(result):
    from pycsamt.models.occam2d.plot import PlotSiteMisfit

    fig = PlotSiteMisfit(result=result, modes=["TM"]).plot()
    assert fig is not None


def test_site_misfit_no_response_raises():
    from pycsamt.models.occam2d.plot import PlotSiteMisfit

    class _FakeResult:
        response = None
        data = None

    with pytest.raises(RuntimeError, match="no response data"):
        PlotSiteMisfit(result=_FakeResult()).plot()


# ---------------------------------------------------------------------------
# PlotResponseGrid
# ---------------------------------------------------------------------------


def test_response_grid_returns_figure(result):
    from pycsamt.models.occam2d.plot import PlotResponseGrid

    fig = PlotResponseGrid(result=result, max_stations=6, n_cols=3).plot()
    assert hasattr(fig, "savefig")


def test_response_grid_rms_in_title(result):
    from pycsamt.models.occam2d.plot import PlotResponseGrid

    fig = PlotResponseGrid(result=result, max_stations=2, n_cols=2).plot()
    titles = [ax.get_title() for ax in fig.axes if ax.get_title()]
    assert any("[" in t for t in titles), "expected RMS annotation in panel titles"


def test_response_grid_te_only(result):
    from pycsamt.models.occam2d.plot import PlotResponseGrid

    fig = PlotResponseGrid(result=result, modes=["TE"], max_stations=4).plot()
    assert fig is not None


def test_response_grid_selected_stations(result):
    from pycsamt.models.occam2d.plot import PlotResponseGrid

    stations = result.data.sites[:4]
    fig = PlotResponseGrid(result=result, stations=list(stations), n_cols=4).plot()
    assert fig is not None


def test_response_grid_no_response_raises():
    from pycsamt.models.occam2d.plot import PlotResponseGrid

    class _FakeResult:
        response = None
        data = None

    with pytest.raises(RuntimeError, match="no response data"):
        PlotResponseGrid(result=_FakeResult()).plot()
