"""Coverage-gap tests for models/mare2dem/plot.py.

Complements test_mare2dem_plot.py by exercising PlotRxParams,
PlotTxParams, plot_poly, and PlotModel's real triangulated-mesh
rendering path -- none of which had any coverage before. Uses the
same real bundled MARE2DEM datasets (data/mare2dem/*) as the sibling
file; a small synthetic Triangle .node/.ele pair is generated for the
mesh-rendering path since none of the bundled datasets ship one.
"""

from __future__ import annotations

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")
import matplotlib.figure
import matplotlib.pyplot as plt

from pycsamt.models.mare2dem.tests.conftest import (
    CSEM_DIR,
    HILL_DIR,
    INVERSION_DIR,
)

_MT_DIR = INVERSION_DIR
_CSEM_DIR = CSEM_DIR
_HILL_DIR = HILL_DIR

_SKIP_MT = pytest.mark.skipif(
    not _MT_DIR.exists(), reason=f"MARE2DEM MT data not found: {_MT_DIR}"
)
_SKIP_CSEM = pytest.mark.skipif(
    not _CSEM_DIR.exists(),
    reason=f"MARE2DEM CSEM data not found: {_CSEM_DIR}",
)
_SKIP_HILL = pytest.mark.skipif(
    not _HILL_DIR.exists(),
    reason=f"MARE2DEM hill data not found: {_HILL_DIR}",
)


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
# PlotConvergence extra branches
# ===========================================================================


class TestPlotConvergenceExtra:
    @_SKIP_MT
    def test_ax_provided_reuses_figure(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence

        fig0, ax = plt.subplots()
        fig = PlotConvergence(result_mt).plot(ax=ax)
        assert fig is fig0

    @_SKIP_MT
    def test_target_rms_line_drawn(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotConvergence

        fig = PlotConvergence(result_mt).plot(target_rms=1.0)
        ax = fig.get_axes()[0]
        assert len(ax.get_lines()) >= 2

    @_SKIP_MT
    def test_savefig_writes_file(self, result_mt, tmp_path):
        from pycsamt.models.mare2dem.plot import PlotConvergence

        out = tmp_path / "conv.png"
        PlotConvergence(result_mt).plot(savefig=out)
        assert out.exists()


# ===========================================================================
# PlotSurveyLayout extra branches
# ===========================================================================


class TestPlotSurveyLayoutExtra:
    @_SKIP_HILL
    def test_no_csem_branch(self, result_hill):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout

        assert result_hill.data.csem is None or not len(
            result_hill.data.csem.receivers or []
        )
        fig, ax = plt.subplots()
        PlotSurveyLayout(result_hill.data).plot(ax=ax)

    @_SKIP_CSEM
    def test_csem_only_no_mt_branch(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout

        fig, ax = plt.subplots()
        PlotSurveyLayout(result_csem.data).plot(ax=ax)

    @_SKIP_MT
    def test_savefig_writes_file(self, result_mt, tmp_path):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout

        out = tmp_path / "layout.png"
        PlotSurveyLayout(result_mt.data).plot(savefig=out)
        assert out.exists()

    @_SKIP_MT
    def test_units_meters(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotSurveyLayout

        fig, ax = plt.subplots()
        PlotSurveyLayout(result_mt.data).plot(ax=ax, units="m")
        assert "(m)" in ax.get_xlabel()


# ===========================================================================
# PlotRxParams — previously fully uncovered
# ===========================================================================


class TestPlotRxParams:
    @_SKIP_MT
    def test_mt_receivers_returns_figure(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        fig = PlotRxParams(result_mt.data).plot()
        assert _is_figure(fig)
        assert len(fig.axes) == 6

    @_SKIP_CSEM
    def test_csem_receivers_returns_figure(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        fig = PlotRxParams(result_csem.data).plot()
        assert _is_figure(fig)

    @_SKIP_MT
    def test_savefig_writes_file(self, result_mt, tmp_path):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        out = tmp_path / "rx.png"
        PlotRxParams(result_mt.data).plot(savefig=out)
        assert out.exists()

    @_SKIP_MT
    def test_units_meters(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        fig = PlotRxParams(result_mt.data).plot(units="m")
        assert "(m)" in fig.axes[-1].get_xlabel()

    @_SKIP_MT
    def test_given_fig_reuses_axes(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        fig0, axes = plt.subplots(6, 1)
        fig = PlotRxParams(result_mt.data).plot(fig=fig0)
        assert fig is fig0

    def test_no_receivers_raises(self):
        from pycsamt.models.mare2dem.plot import PlotRxParams

        class _Empty:
            csem = None
            mt = None

        with pytest.raises(ValueError):
            PlotRxParams(_Empty()).plot()


# ===========================================================================
# PlotTxParams — previously fully uncovered
# ===========================================================================


class TestPlotTxParams:
    @_SKIP_CSEM
    def test_csem_transmitters_returns_figure(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotTxParams

        fig = PlotTxParams(result_csem.data).plot()
        assert _is_figure(fig)
        assert len(fig.axes) >= 1

    @_SKIP_CSEM
    def test_savefig_writes_file(self, result_csem, tmp_path):
        from pycsamt.models.mare2dem.plot import PlotTxParams

        out = tmp_path / "tx.png"
        PlotTxParams(result_csem.data).plot(savefig=out)
        assert out.exists()

    @_SKIP_CSEM
    def test_units_meters(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotTxParams

        fig = PlotTxParams(result_csem.data).plot(units="m")
        assert "(m)" in fig.axes[-1].get_xlabel()

    @_SKIP_CSEM
    def test_given_fig_reuses_axes(self, result_csem):
        from pycsamt.models.mare2dem.plot import PlotTxParams

        fig0, axes = plt.subplots(5, 1)
        fig = PlotTxParams(result_csem.data).plot(fig=fig0)
        assert fig is fig0

    @_SKIP_MT
    def test_no_transmitters_raises(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotTxParams

        with pytest.raises(ValueError):
            PlotTxParams(result_mt.data).plot()


# ===========================================================================
# plot_poly — previously fully uncovered
# ===========================================================================


class TestPlotPoly:
    @_SKIP_MT
    def test_returns_axes(self):
        from pycsamt.models.mare2dem.plot import plot_poly

        ax = plot_poly(_MT_DIR / "demo.poly")
        assert ax is not None
        assert len(ax.get_lines()) >= 1

    @_SKIP_HILL
    def test_given_ax_reused(self):
        from pycsamt.models.mare2dem.plot import plot_poly

        fig, ax0 = plt.subplots()
        ax = plot_poly(_HILL_DIR / "hill.poly", ax=ax0)
        assert ax is ax0

    @_SKIP_MT
    def test_savefig_writes_file(self, tmp_path):
        from pycsamt.models.mare2dem.plot import plot_poly

        out = tmp_path / "poly.png"
        plot_poly(_MT_DIR / "demo.poly", savefig=out)
        assert out.exists()

    def test_empty_pslg_returns_ax_unchanged(self, tmp_path):
        from pycsamt.models.mare2dem.plot import plot_poly

        empty_poly = tmp_path / "empty.poly"
        empty_poly.write_text("0 2 0 0\n0\n0\n")
        fig, ax0 = plt.subplots()
        ax = plot_poly(empty_poly, ax=ax0)
        assert ax is ax0
        assert len(ax.get_lines()) == 0

    @_SKIP_MT
    def test_custom_style(self):
        from pycsamt.models.mare2dem.plot import plot_poly

        ax = plot_poly(_MT_DIR / "demo.poly", linewidth=2.0, color="red")
        assert ax.get_lines()[0].get_color() == "red"


# ===========================================================================
# PlotModel — __init__ branches + real triangulated-mesh rendering
# ===========================================================================


class TestPlotModelExtra:
    @_SKIP_MT
    def test_init_from_resistivity_file_directly(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotModel
        from pycsamt.models.mare2dem.iotools.resistivity import (
            ResistivityFile,
        )

        # InversionResult.model is populated with a raw ResistivityFile
        # in practice (not the ResistivityModel wrapper) -- this is the
        # normal path, exercising the `isinstance(rf_like, ResistivityFile)`
        # branch of PlotModel.__init__.
        assert isinstance(result_mt.model, ResistivityFile)
        fig = PlotModel(result_mt).plot()
        assert _is_figure(fig)

    @_SKIP_MT
    def test_init_from_resistivity_model_wrapper(self, result_mt, monkeypatch):
        from pycsamt.models.mare2dem.mesh import ResistivityModel
        from pycsamt.models.mare2dem.plot import PlotModel

        wrapped = ResistivityModel(_MT_DIR / "demo.0.resistivity")
        monkeypatch.setattr(result_mt, "model", wrapped)
        fig = PlotModel(result_mt).plot()
        assert _is_figure(fig)

    def test_init_invalid_type_raises(self):
        from pycsamt.models.mare2dem.plot import PlotModel

        with pytest.raises(TypeError):
            PlotModel(12345)

    @_SKIP_MT
    def test_histogram_ax_provided_reuses_figure(self, result_mt):
        from pycsamt.models.mare2dem.mesh import ResistivityModel
        from pycsamt.models.mare2dem.plot import PlotModel

        fig0, ax = plt.subplots()
        fig = PlotModel(ResistivityModel()).plot(ax=ax)
        assert fig is fig0

    def test_mesh_rendering_with_synthetic_triangulation(self, tmp_path):
        """Copy a real .resistivity next to a hand-built .node/.ele pair
        so PlotModel._load_mesh finds a mesh and takes the tripcolor
        path instead of the histogram fallback (never exercised by any
        bundled dataset, which all lack the mesh output files)."""
        import shutil

        from pycsamt.models.mare2dem.iotools.poly import write_triangulation
        from pycsamt.models.mare2dem.mesh import ResistivityModel
        from pycsamt.models.mare2dem.plot import PlotModel

        if not _MT_DIR.exists():
            pytest.skip("MARE2DEM MT data not found")

        src = _MT_DIR / "demo.0.resistivity"
        dst = tmp_path / "demo.0.resistivity"
        shutil.copy(src, dst)

        model = ResistivityModel(dst)
        stem = __import__("pathlib").Path(model._rf.poly_file).stem
        n_reg = len(model._rf.resistivity)
        assert n_reg > 0

        nodes = np.array(
            [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
        )
        triangles = np.array([[0, 1, 2], [1, 3, 2]])
        region_attrs = np.array([1, min(2, n_reg)])
        write_triangulation(
            nodes, triangles, region_attrs, tmp_path / f"{stem}.node"
        )

        fig = PlotModel(model).plot()
        assert _is_figure(fig)
        ax = fig.get_axes()[0]
        assert ax.collections  # tripcolor added a PolyCollection

    def test_mesh_rendering_with_explicit_vmin_vmax(self, tmp_path):
        import shutil

        from pycsamt.models.mare2dem.iotools.poly import write_triangulation
        from pycsamt.models.mare2dem.mesh import ResistivityModel
        from pycsamt.models.mare2dem.plot import PlotModel

        if not _MT_DIR.exists():
            pytest.skip("MARE2DEM MT data not found")

        src = _MT_DIR / "demo.0.resistivity"
        dst = tmp_path / "demo.0.resistivity"
        shutil.copy(src, dst)
        model = ResistivityModel(dst)
        stem = __import__("pathlib").Path(model._rf.poly_file).stem

        nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        triangles = np.array([[0, 1, 2]])
        write_triangulation(
            nodes, triangles, np.array([1]), tmp_path / f"{stem}.node"
        )

        fig = PlotModel(model).plot(vmin=0.5, vmax=3.5)
        assert _is_figure(fig)

    def test_no_mesh_files_falls_back_to_histogram(self, tmp_path):
        import shutil

        from pycsamt.models.mare2dem.mesh import ResistivityModel
        from pycsamt.models.mare2dem.plot import PlotModel

        if not _MT_DIR.exists():
            pytest.skip("MARE2DEM MT data not found")

        src = _MT_DIR / "demo.0.resistivity"
        dst = tmp_path / "demo.0.resistivity"
        shutil.copy(src, dst)
        model = ResistivityModel(dst)

        fig = PlotModel(model).plot()
        assert _is_figure(fig)
        ax = fig.get_axes()[0]
        assert "distribution" in ax.get_title().lower()


# ===========================================================================
# PlotResponse extra branches
# ===========================================================================


class TestPlotResponseExtra:
    @_SKIP_MT
    def test_station_filter_matches_real_receiver(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotResponse

        if result_mt.response is None:
            pytest.skip("No response file in MT fixture.")
        name = result_mt.data.mt.receiver_name[0]
        fig = PlotResponse(result_mt).plot(station=name)
        assert _is_figure(fig)
        assert len(fig.axes) == 2

    @_SKIP_MT
    def test_station_filter_unknown_raises(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotResponse

        with pytest.raises(ValueError, match="not found"):
            PlotResponse(result_mt).plot(station="__no_such_station__")

    @_SKIP_MT
    def test_no_predicted_response_still_plots_observed(
        self, result_mt, monkeypatch
    ):
        from pycsamt.models.mare2dem.plot import PlotResponse

        if result_mt.response is None:
            pytest.skip("No response file in MT fixture.")
        monkeypatch.setattr(result_mt, "response", None)
        fig = PlotResponse(result_mt).plot(max_rx=1)
        assert _is_figure(fig)

    @_SKIP_MT
    def test_no_frequencies_raises(self, result_mt):
        from pycsamt.models.mare2dem.plot import PlotResponse

        real_em = result_mt.data
        original = real_em.mt.frequencies
        real_em.mt.frequencies = np.array([])
        try:
            with pytest.raises(ValueError, match="frequencies"):
                PlotResponse(result_mt).plot()
        finally:
            # restore for other tests sharing the module-scoped fixture
            real_em.mt.frequencies = original
