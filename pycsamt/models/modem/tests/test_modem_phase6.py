"""Phase 6 tests — PlotMisfit, PlotModel2D, PlotModel3D, PlotResponse, PlotPseudo."""

from pathlib import Path

import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

_EX2D = (
    Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "2D_MT" / "BLOCK2"
)
_EX3D = (
    Path(__file__).parents[4] / "ModEMv626" / "ModEM" / "examples" / "3D_MT" / "BLOCK2"
)
_BUNDLED_3D = (
    Path(__file__).parents[4] / "data" / "modem" / "willy_27freq_watex_line02_sample"
)

_RESULT_3D = _BUNDLED_3D if _BUNDLED_3D.exists() else _EX3D


@pytest.fixture(scope="module")
def result_2d():
    from pycsamt.models.modem.results import InversionResult

    if not _EX2D.exists():
        pytest.skip("2D ModEM example data not available")
    return InversionResult(_EX2D)


@pytest.fixture(scope="module")
def result_3d():
    from pycsamt.models.modem.results import InversionResult

    if not _RESULT_3D.exists():
        pytest.skip("3D ModEM example data not available")
    return InversionResult(_RESULT_3D)


# ---------------------------------------------------------------------------
# PlotMisfit
# ---------------------------------------------------------------------------


def test_plot_misfit_returns_figure(result_3d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotMisfit

    fig = PlotMisfit(result=result_3d).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    matplotlib.pyplot.close(fig)


def test_plot_misfit_no_result_raises():
    from pycsamt.models.modem.plot import PlotMisfit

    with pytest.raises(ValueError, match="No InversionResult"):
        PlotMisfit().plot()


def test_plot_misfit_no_log_raises(tmp_path):
    from pycsamt.models.modem.plot import PlotMisfit
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    with pytest.raises(ValueError, match="no log"):
        PlotMisfit(result=r).plot()


# ---------------------------------------------------------------------------
# PlotModel2D
# ---------------------------------------------------------------------------


def test_plot_model2d_returns_figure(result_2d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotModel2D

    fig = PlotModel2D(result=result_2d, which="final").plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    matplotlib.pyplot.close(fig)


def test_plot_model2d_dynamic_section(result_2d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotModel2D

    fig = PlotModel2D(result=result_2d, section="dynamic").plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    assert len(fig.axes) >= 2
    matplotlib.pyplot.close(fig)


def test_plot_model2d_initial(result_2d):
    from pycsamt.models.modem.plot import PlotModel2D

    fig = PlotModel2D(result=result_2d, which="initial", depth_max=5000.0).plot()
    assert fig is not None
    matplotlib.pyplot.close(fig)


def test_plot_model2d_no_model_raises(tmp_path):
    from pycsamt.models.modem.plot import PlotModel2D
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    with pytest.raises(ValueError, match="No final model"):
        PlotModel2D(result=r).plot()


# ---------------------------------------------------------------------------
# PlotModel3D
# ---------------------------------------------------------------------------


def test_plot_model3d_returns_figure(result_3d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotModel3D

    fig = PlotModel3D(result=result_3d, n_cols=2).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    matplotlib.pyplot.close(fig)


def test_plot_model3d_custom_depths(result_3d):
    from pycsamt.models.modem.plot import PlotModel3D

    depths = [1000.0, 5000.0, 20000.0]
    fig = PlotModel3D(result=result_3d, depths=depths, n_cols=3).plot()
    assert fig is not None
    matplotlib.pyplot.close(fig)


# ---------------------------------------------------------------------------
# PlotResponse
# ---------------------------------------------------------------------------


def test_plot_response_returns_figure(result_3d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotResponse

    fig = PlotResponse(result=result_3d, max_stations=4).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    matplotlib.pyplot.close(fig)


def test_plot_response_no_data_raises(tmp_path):
    from pycsamt.models.modem.plot import PlotResponse
    from pycsamt.models.modem.results import InversionResult

    r = InversionResult(tmp_path)
    with pytest.raises(ValueError, match="no data_obs"):
        PlotResponse(result=r).plot()


# ---------------------------------------------------------------------------
# PlotPseudo
# ---------------------------------------------------------------------------


def test_plot_pseudo_returns_figure(result_3d):
    import matplotlib.figure

    from pycsamt.models.modem.plot import PlotPseudo

    # discover available component from first block
    comp = result_3d.data_obs.blocks[0]["rows"][0][5] if result_3d.data_obs else "TE"
    fig = PlotPseudo(result=result_3d, component=comp).plot()
    assert isinstance(fig, matplotlib.figure.Figure)
    matplotlib.pyplot.close(fig)
