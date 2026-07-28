# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for TDEMController — loads TDEM survey data and drives all TDEM
canvas renders.

Real data
---------
data/TEMAVG/JIANGSU/ — 55 real TEMAVG stations (AVG/Z/LOG triplets).
A small 3-station subset (TEM100, TEM1020, TEM1060) is copied into a
session tmp dir so loading is fast (~0.3s) while still exercising the
real ``read_temavg_survey`` / ``read_temavg_soundings`` parsers and the
real ``pycsamt.tdem`` plot classes end to end — no synthetic data needed
since the real fixture set is small and readily available.

Strategy
--------
* TDEMController.draw() must never *raise*; failures (missing coordinate
  data, unknown class names, no data loaded) are annotated into the
  figure instead.
* Every catalogue entry in TDEM_GROUPS is exercised against the real
  small survey.
* load_folder() success/failure paths and progress_cb invocation are
  checked directly against real data and a real empty/missing folder.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from pycsamt.app.desktop.controllers.tdem_controller import (
    DASHBOARD_PLOTS,
    DECAY_PLOTS,
    MAP_PLOTS,
    SECTION_PLOTS,
    TDEM_GROUPS,
    TDEMController,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_JIANGSU = _ROOT / "data" / "TEMAVG" / "JIANGSU"
_HAS_JIANGSU = _JIANGSU.exists() and any(_JIANGSU.glob("*.AVG"))
_STEMS = ["TEM100", "TEM1020", "TEM1060"]


# ── Fixtures ──────────────────────────────────────────────────────────────────


@pytest.fixture(scope="session")
def tdem_dir(tmp_path_factory):
    """
    A small real TDEM survey folder (3 stations, AVG/Z/LOG) copied out of
    data/TEMAVG/JIANGSU/ so loading is fast but still 100% real data.
    """
    if not _HAS_JIANGSU:
        pytest.skip("data/TEMAVG/JIANGSU not available")
    dest = tmp_path_factory.mktemp("tdem_jiangsu_subset")
    for stem in _STEMS:
        for f in _JIANGSU.glob(stem + ".*"):
            shutil.copy(f, dest / f.name)
    return dest


@pytest.fixture(scope="session")
def loaded_ctrl(tdem_dir):
    """A TDEMController with the small real survey loaded once per session."""
    ctrl = TDEMController()
    ok = ctrl.load_folder(str(tdem_dir))
    assert ok is True
    return ctrl


def _fig():
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


def _all_texts(fig):
    return [t.get_text() for ax in fig.axes for t in ax.texts]


# ── Construction / defaults ────────────────────────────────────────────────────


class TestConstruction:
    def test_creates(self):
        ctrl = TDEMController()
        assert ctrl is not None

    def test_default_dark(self):
        assert TDEMController().dark is True

    def test_default_not_loaded(self):
        assert TDEMController().is_loaded() is False

    def test_default_summary_empty(self):
        assert TDEMController().summary == ""

    def test_default_folder_empty(self):
        assert TDEMController()._folder == ""

    def test_default_soundings_empty_list(self):
        assert TDEMController()._soundings == []


# ── load_folder — success path ─────────────────────────────────────────────────


class TestLoadFolderSuccess:
    def test_returns_true(self, tdem_dir):
        ctrl = TDEMController()
        assert ctrl.load_folder(str(tdem_dir)) is True

    def test_is_loaded_after_success(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert ctrl.is_loaded() is True

    def test_survey_populated(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert ctrl._survey is not None

    def test_soundings_populated(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert len(ctrl._soundings) > 0

    def test_folder_recorded(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert ctrl._folder == str(tdem_dir)

    def test_summary_reflects_counts(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert str(tdem_dir) in ctrl.summary
        assert f"AVG files: {ctrl._survey.n_avg_files}" in ctrl.summary
        assert f"Z files: {ctrl._survey.n_z_files}" in ctrl.summary
        assert f"Soundings: {len(ctrl._soundings)}" in ctrl.summary

    def test_progress_cb_called_10_50_100(self, tdem_dir):
        calls = []
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir), progress_cb=calls.append)
        assert calls == [10, 50, 100]

    def test_no_progress_cb_does_not_raise(self, tdem_dir):
        ctrl = TDEMController()
        assert ctrl.load_folder(str(tdem_dir), progress_cb=None) is True


# ── load_folder — failure paths ────────────────────────────────────────────────


class TestLoadFolderFailure:
    def test_nonexistent_folder_returns_false(self):
        ctrl = TDEMController()
        assert ctrl.load_folder("Z:/no/such/tdem/folder/xyz") is False

    def test_nonexistent_folder_not_loaded(self):
        ctrl = TDEMController()
        ctrl.load_folder("Z:/no/such/tdem/folder/xyz")
        assert ctrl.is_loaded() is False

    def test_nonexistent_folder_summary_has_error(self):
        ctrl = TDEMController()
        ctrl.load_folder("Z:/no/such/tdem/folder/xyz")
        assert ctrl.summary.startswith("Load error:")

    def test_empty_folder_returns_false(self, tmp_path):
        ctrl = TDEMController()
        assert ctrl.load_folder(str(tmp_path)) is False

    def test_empty_folder_not_loaded(self, tmp_path):
        ctrl = TDEMController()
        ctrl.load_folder(str(tmp_path))
        assert ctrl.is_loaded() is False

    def test_empty_folder_summary_has_error(self, tmp_path):
        ctrl = TDEMController()
        ctrl.load_folder(str(tmp_path))
        assert "Load error" in ctrl.summary

    def test_failure_progress_cb_only_called_with_10(self, tmp_path):
        calls = []
        ctrl = TDEMController()
        ctrl.load_folder(str(tmp_path), progress_cb=calls.append)
        assert calls == [10]

    def test_soundings_failure_does_not_block_survey_success(
        self, tdem_dir, monkeypatch
    ):
        """
        If read_temavg_soundings() raises but read_temavg_survey() succeeds,
        load_folder() still returns True with an empty soundings list.
        """
        import pycsamt.tdem as tdem_mod

        def _boom(*_a, **_kw):
            raise RuntimeError("boom")

        monkeypatch.setattr(tdem_mod, "read_temavg_soundings", _boom)

        calls = []
        ctrl = TDEMController()
        ok = ctrl.load_folder(str(tdem_dir), progress_cb=calls.append)

        assert ok is True
        assert ctrl._soundings == []
        assert calls == [10, 50, 100]
        assert "Soundings: 0" in ctrl.summary


# ── is_loaded() ────────────────────────────────────────────────────────────────


class TestIsLoaded:
    def test_false_before_load(self):
        assert TDEMController().is_loaded() is False

    def test_true_after_successful_load(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        assert ctrl.is_loaded() is True

    def test_false_after_failed_load(self, tmp_path):
        ctrl = TDEMController()
        ctrl.load_folder(str(tmp_path))
        assert ctrl.is_loaded() is False


# ── _build_summary() ────────────────────────────────────────────────────────────


class TestBuildSummary:
    def test_no_survey_message(self):
        ctrl = TDEMController()
        ctrl._build_summary()
        assert ctrl.summary == "No data loaded"

    def test_populated_after_load(self, tdem_dir):
        ctrl = TDEMController()
        ctrl.load_folder(str(tdem_dir))
        # Re-invoking directly must reproduce the same summary deterministically.
        prev = ctrl.summary
        ctrl._build_summary()
        assert ctrl.summary == prev


# ── draw() — no-data / error-annotation paths ───────────────────────────────────


class TestDrawNoData:
    def test_no_data_does_not_raise(self):
        ctrl = TDEMController()
        fig = _fig()
        ctrl.draw("PlotDecayCurve", True, "soundings", fig)
        _close()

    def test_no_data_shows_load_message(self):
        ctrl = TDEMController()
        fig = _fig()
        ctrl.draw("PlotDecayCurve", True, "soundings", fig)
        texts = " ".join(_all_texts(fig)).lower()
        assert "load" in texts and "browse" in texts
        _close()

    def test_no_data_returns_none(self):
        ctrl = TDEMController()
        fig = _fig()
        ret = ctrl.draw("PlotDecayCurve", True, "soundings", fig)
        assert ret is None
        _close()

    def test_unknown_class_does_not_raise(self, loaded_ctrl):
        fig = _fig()
        loaded_ctrl.draw("NotARealPlotClass", False, "survey", fig)
        _close()

    def test_unknown_class_annotates_figure(self, loaded_ctrl):
        fig = _fig()
        loaded_ctrl.draw("NotARealPlotClass", False, "survey", fig)
        texts = " ".join(_all_texts(fig)).lower()
        assert "not found" in texts and "notarealplotclass" in texts
        _close()

    def test_unknown_class_returns_none(self, loaded_ctrl):
        fig = _fig()
        ret = loaded_ctrl.draw("NotARealPlotClass", False, "survey", fig)
        assert ret is None
        _close()


# ── draw() — real catalogue, has_ax=True in-place plots ─────────────────────────


def _flatten(*groups):
    out = []
    for g in groups:
        out.extend(g)
    return out


_HAS_AX_TRUE_PLOTS = [
    (label, cls, has_ax, key)
    for (label, cls, has_ax, key) in _flatten(DECAY_PLOTS, SECTION_PLOTS, MAP_PLOTS)
    if has_ax
]
_HAS_AX_FALSE_PLOTS = [
    (label, cls, has_ax, key)
    for (label, cls, has_ax, key) in _flatten(
        DECAY_PLOTS, SECTION_PLOTS, MAP_PLOTS, DASHBOARD_PLOTS
    )
    if not has_ax
]


class TestDrawHasAxTrue:
    @pytest.mark.parametrize(
        "label,class_name,has_ax,data_key",
        _HAS_AX_TRUE_PLOTS,
        ids=[p[1] for p in _HAS_AX_TRUE_PLOTS],
    )
    def test_does_not_raise(self, loaded_ctrl, label, class_name, has_ax, data_key):
        fig = _fig()
        loaded_ctrl.draw(class_name, has_ax, data_key, fig)  # must never raise
        _close()

    @pytest.mark.parametrize(
        "label,class_name,has_ax,data_key",
        _HAS_AX_TRUE_PLOTS,
        ids=[p[1] for p in _HAS_AX_TRUE_PLOTS],
    )
    def test_returns_none_and_draws_into_fig(
        self, loaded_ctrl, label, class_name, has_ax, data_key
    ):
        fig = _fig()
        ret = loaded_ctrl.draw(class_name, has_ax, data_key, fig)
        assert ret is None
        assert len(fig.axes) >= 1
        _close()


# ── draw() — real catalogue, has_ax=False own-figure plots ──────────────────────


class TestDrawHasAxFalse:
    @pytest.mark.parametrize(
        "label,class_name,has_ax,data_key",
        _HAS_AX_FALSE_PLOTS,
        ids=[p[1] for p in _HAS_AX_FALSE_PLOTS],
    )
    def test_does_not_raise(self, loaded_ctrl, label, class_name, has_ax, data_key):
        fig = _fig()
        loaded_ctrl.draw(class_name, has_ax, data_key, fig)  # must never raise
        _close()

    def test_apparent_resistivity_returns_populated_figure(self, loaded_ctrl):
        """PlotTransformedRho has real soundings and coordinate-free data,
        so it should produce a genuine multi-axes figure (not an error)."""
        fig = _fig()
        ret = loaded_ctrl.draw("PlotTransformedRho", False, "soundings", fig)
        assert ret is not None
        assert len(ret.axes) >= 1
        texts = " ".join(t.get_text() for ax in ret.axes for t in ax.texts)
        assert "error" not in texts.lower()
        _close()

    def test_full_dashboard_returns_populated_figure(self, loaded_ctrl):
        fig = _fig()
        ret = loaded_ctrl.draw("PlotTEMDashboard", False, "dashboard", fig)
        assert ret is not None
        assert len(ret.axes) >= 1
        _close()

    def test_survey_overview_missing_coords_shows_no_figure_message(self, loaded_ctrl):
        """The small JIANGSU subset carries no coordinate rows, so
        PlotSurveyOverview (which needs PlotSurveyMap internally) fails
        inside _call_figure_cls and draw() falls back to the
        'No figure produced' annotation instead of raising."""
        fig = _fig()
        ret = loaded_ctrl.draw("PlotSurveyOverview", False, "survey", fig)
        assert ret is None
        texts = " ".join(_all_texts(fig)).lower()
        assert "no figure produced" in texts
        _close()


# ── draw() — error-annotation path (real data, real failure) ────────────────────


class TestDrawErrorAnnotation:
    def test_survey_map_missing_coords_annotates_error(self, loaded_ctrl):
        """PlotSurveyMap (has_ax=True) requires coordinate rows; the small
        JIANGSU subset has none, so draw() must catch the real ValueError
        and annotate it rather than raising."""
        fig = _fig()
        ret = loaded_ctrl.draw("PlotSurveyMap", True, "survey", fig)
        assert ret is None
        texts = " ".join(_all_texts(fig)).lower()
        assert "plotsurveymap error" in texts
        _close()

    def test_elevation_profile_missing_coords_annotates_error(self, loaded_ctrl):
        fig = _fig()
        ret = loaded_ctrl.draw("PlotElevationProfile", True, "survey", fig)
        assert ret is None
        texts = " ".join(_all_texts(fig)).lower()
        assert "plotelevationprofile error" in texts
        _close()


# ── draw() — signature-dispatch edge cases ──────────────────────────────────────


class TestDrawSignatureDispatch:
    def test_axes_kwarg_branch_is_used_when_has_ax_true(self, loaded_ctrl):
        """PlotTransformedRho.plot(axes=...) uses 'axes', not 'ax'; forcing
        has_ax=True exercises the 'axes in sig.parameters' branch inside
        draw() instead of the normal has_ax=False catalogue path."""
        fig = _fig()
        loaded_ctrl.draw("PlotTransformedRho", True, "soundings", fig)
        # Whatever happens (success or annotated error), it must not raise
        # and must leave at least one axes behind.
        assert len(fig.axes) >= 1
        _close()

    def test_dashboard_has_ax_true_defensive_branch(self, loaded_ctrl):
        """The catalogue always uses has_ax=False for the dashboard, but
        draw() has a defensive has_ax=True branch for it. PlotTEMDashboard
        .plot() takes no ax/axes argument, so this must fail gracefully
        into the error annotation rather than raising."""
        fig = _fig()
        ret = loaded_ctrl.draw("PlotTEMDashboard", True, "dashboard", fig)
        assert ret is None
        texts = " ".join(_all_texts(fig)).lower()
        assert "plottemdashboard error" in texts
        _close()


# ── draw() — tight_layout failure is swallowed, never raises ────────────────────


class TestTightLayoutGuard:
    """
    draw() wraps both its ``fig.tight_layout(...)`` calls (has_ax=True path
    and has_ax=False path) in ``except Exception: pass``. Force
    tight_layout() itself to raise to exercise those guards directly.
    """

    def test_has_ax_true_tight_layout_failure_swallowed(self, loaded_ctrl, monkeypatch):
        import matplotlib.figure

        def _boom(self, *_a, **_kw):
            raise RuntimeError("tight_layout boom")

        monkeypatch.setattr(matplotlib.figure.Figure, "tight_layout", _boom)
        fig = _fig()
        ret = loaded_ctrl.draw("PlotDecayCurve", True, "soundings", fig)
        assert ret is None
        _close()

    def test_has_ax_false_tight_layout_failure_swallowed(
        self, loaded_ctrl, monkeypatch
    ):
        import matplotlib.figure

        def _boom(self, *_a, **_kw):
            raise RuntimeError("tight_layout boom")

        monkeypatch.setattr(matplotlib.figure.Figure, "tight_layout", _boom)
        fig = _fig()
        ret = loaded_ctrl.draw("PlotTransformedRho", False, "soundings", fig)
        assert ret is not None  # own-figure path still returns the figure
        _close()


# ── draw() — dark / light theming ────────────────────────────────────────────────


class TestDrawThemes:
    def test_dark_and_light_do_not_raise(self, loaded_ctrl):
        for dark in (True, False):
            loaded_ctrl.dark = dark
            fig = _fig()
            loaded_ctrl.draw("PlotDecayCurve", True, "soundings", fig)
            _close()
        loaded_ctrl.dark = True  # restore default for other tests

    def test_no_data_light_mode_does_not_raise(self):
        ctrl = TDEMController()
        ctrl.dark = False
        fig = _fig()
        ctrl.draw("PlotDecayCurve", True, "soundings", fig)
        _close()


# ── _call_figure_cls() — direct unit tests ───────────────────────────────────────


class TestCallFigureCls:
    def test_dashboard_data_none_returns_figure(self, loaded_ctrl):
        import pycsamt.tdem as tdem

        result = loaded_ctrl._call_figure_cls(tdem.PlotTEMDashboard, None)
        assert result is not None
        assert hasattr(result, "get_axes")
        _close()

    def test_non_dashboard_data_returns_figure(self, loaded_ctrl):
        import pycsamt.tdem as tdem

        result = loaded_ctrl._call_figure_cls(
            tdem.PlotTransformedRho, loaded_ctrl._soundings
        )
        assert result is not None
        assert hasattr(result, "get_axes")
        _close()

    def test_raising_cls_returns_none(self, loaded_ctrl):
        class _Boom:
            def __init__(self, data):
                pass

            def plot(self):
                raise RuntimeError("kaboom")

        result = loaded_ctrl._call_figure_cls(_Boom, "data")
        assert result is None

    def test_fallback_new_fignum_returned_when_no_get_axes(self, loaded_ctrl):
        """If plot() creates a new pyplot figure but returns something
        without get_axes, _call_figure_cls falls back to the newest
        matplotlib fignum created during the call."""

        class _CreatesFigNoReturn:
            def __init__(self, data):
                pass

            def plot(self):
                fig = plt.figure()
                fig.add_subplot(111).plot([1, 2, 3])
                return None  # no get_axes attribute

        result = loaded_ctrl._call_figure_cls(_CreatesFigNoReturn, "data")
        assert result is not None
        assert hasattr(result, "get_axes")
        _close()

    def test_no_figure_and_no_get_axes_returns_none(self, loaded_ctrl):
        class _NoFigure:
            def __init__(self, data):
                pass

            def plot(self):
                return "not a figure"

        result = loaded_ctrl._call_figure_cls(_NoFigure, "data")
        assert result is None


# ── Catalogue sanity ──────────────────────────────────────────────────────────


class TestCatalogue:
    def test_tdem_groups_structure(self):
        assert len(TDEM_GROUPS) == 4
        for group_name, plots in TDEM_GROUPS:
            assert isinstance(group_name, str)
            for entry in plots:
                assert len(entry) == 4
                label, class_name, has_ax, data_key = entry
                assert isinstance(label, str)
                assert isinstance(class_name, str)
                assert isinstance(has_ax, bool)
                assert data_key in ("survey", "soundings", "dashboard")

    def test_all_catalogue_classes_exist_in_tdem_module(self):
        import pycsamt.tdem as tdem

        for _group_name, plots in TDEM_GROUPS:
            for _label, class_name, _has_ax, _data_key in plots:
                assert hasattr(
                    tdem, class_name
                ), f"{class_name} missing from pycsamt.tdem"
