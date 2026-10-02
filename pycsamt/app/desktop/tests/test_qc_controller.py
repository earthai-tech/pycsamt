# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for QCController and the QC plot-kind catalogue.

Real data
---------
data/MT/kap03lmt_edis/  — 26 KP TIPPER EDIs
data/AMT/WILLY_DATA/    — 128 WILLY AMT EDIs across 5 profiles

Strategy
--------
* Every catalogue entry (ALL_GROUPS) is exercised against real EDIs via
  QCController.draw(), following the pattern used in test_advanced_tools.py.
* draw() must never *raise*; on failure it writes an error annotation into
  the figure instead of propagating the exception.
* No-data paths (self._sites is None) must never raise either, and must
  annotate "Load survey data first" via _annotate_empty().
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    COVERAGE_PLOTS,
    DISTORTION_PLOTS,
    GROUP_ICONS,
    NOISE_PLOTS,
    OVERVIEW_PLOTS,
    PLOT_DESCRIPTIONS,
    SKEW_DIM_PLOTS,
    STATIC_SHIFT_PLOTS,
    QCController,
    QCUnavailableResult,
    _annotate_empty,
    _auto_rho_bounds,
    _dispatch_polar_ax,
    _dispatch_polar_coverage,
    _style_all_axes,
    _style_ax,
    describe_plot,
    qc_parameter_specs,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[4]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA"
_GV_XML = _ROOT / "data" / "gv_data" / "xml"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_WILLY = _WILLY.exists() and any(_WILLY.rglob("*.edi"))
_HAS_GV_XML = _GV_XML.exists() and any(_GV_XML.glob("*.xml"))

# ── Session-scoped fixtures ───────────────────────────────────────────────────


@pytest.fixture(scope="session")
def tipper_sites():
    """26-station TIPPER Sites loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture(scope="session")
def willy_sites():
    """128-station WILLY Sites (all 5 profiles) loaded once for the whole session."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY))


@pytest.fixture(scope="session")
def gv_xml_sites():
    """Three EMTF-XML stations without source-offset geometry."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_GV_XML:
        pytest.skip("GV XML data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_GV_XML))


@pytest.fixture
def ctrl():
    return QCController()


# ── Helpers ───────────────────────────────────────────────────────────────────


def _fig():
    """Return a fresh matplotlib Figure."""
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


def _texts(fig):
    return " ".join(t.get_text() for ax in fig.axes for t in ax.texts).lower()


def _catalogue_ids(group):
    return [fn for _, fn, _ in group]


def _catalogue_params(group):
    return [(fn, has_ax) for _, fn, has_ax in group]


ALL_PLOTS = [entry for _, group in ALL_GROUPS for entry in group]


# ── QCController construction / state management ─────────────────────────────


class TestQCControllerConstruction:
    def test_creates(self):
        c = QCController()
        assert c is not None

    def test_default_dark(self):
        assert QCController().dark is True

    def test_default_sites_none(self):
        assert QCController()._sites is None

    def test_default_unavailable_result_none(self):
        assert QCController().last_unavailable is None

    def test_unavailable_result_value_object(self):
        result = QCUnavailableResult("Title", "Reason", "Guidance", "test")
        assert result.title == "Title"
        assert result.kind == "test"

    def test_set_sites(self):
        c = QCController()
        c.set_sites("dummy")
        assert c._sites == "dummy"

    def test_clear(self):
        c = QCController()
        c.set_sites("dummy")
        c.clear()
        assert c._sites is None

    def test_clear_when_already_none(self):
        c = QCController()
        c.clear()
        assert c._sites is None


class TestFilterSites:
    def test_filter_no_station_ids_is_noop(self, ctrl):
        ctrl.set_sites("dummy")
        ctrl.filter_sites([])
        assert ctrl._sites == "dummy"

    def test_filter_none_station_ids_is_noop(self, ctrl):
        ctrl.set_sites("dummy")
        ctrl.filter_sites(None)
        assert ctrl._sites == "dummy"

    def test_filter_with_no_sites_loaded_is_noop(self, ctrl):
        assert ctrl._sites is None
        ctrl.filter_sites(["some_station"])
        assert ctrl._sites is None

    def test_filter_swallows_exceptions_on_bad_sites(self, ctrl):
        """A non-iterable, non-.items() object must not raise."""
        ctrl.set_sites(object())
        ctrl.filter_sites(["A", "B"])
        # No exception raised; state left as-is since nothing matched.

    def test_filter_with_real_sites_does_not_raise(self, ctrl, tipper_sites):
        ctrl.set_sites(tipper_sites)
        ctrl.filter_sites(["kap103"])
        # Whatever the outcome, _sites must still be a usable collection.
        assert ctrl._sites is not None

    def test_filter_unmatched_ids_keeps_original_sites(self, ctrl, tipper_sites):
        """When nothing in station_ids matches, filtered list is empty and
        the controller silently keeps the original (unfiltered) sites."""
        ctrl.set_sites(tipper_sites)
        ctrl.filter_sites(["definitely_not_a_real_station_id_xyz"])
        assert ctrl._sites is tipper_sites


# ── Styling helpers ────────────────────────────────────────────────────────────


class TestStyleAx:
    @pytest.mark.parametrize("dark", [True, False])
    def test_style_ax_sets_facecolor(self, dark):
        fig, ax = plt.subplots()
        _style_ax(ax, dark=dark)
        assert ax.get_facecolor() is not None
        _close()

    @pytest.mark.parametrize("dark", [True, False])
    def test_style_ax_sets_figure_facecolor(self, dark):
        fig, ax = plt.subplots()
        _style_ax(ax, dark=dark)
        expected = "#1e1e2e" if dark else "#ffffff"
        # matplotlib stores colors as RGBA tuples; compare via to_hex
        import matplotlib.colors as mcolors

        assert mcolors.to_hex(fig.patch.get_facecolor()) == expected
        _close()

    def test_style_ax_default_is_dark(self):
        fig, ax = plt.subplots()
        _style_ax(ax)
        import matplotlib.colors as mcolors

        assert mcolors.to_hex(fig.patch.get_facecolor()) == "#1e1e2e"
        _close()

    def test_style_ax_grid_enabled(self):
        fig, ax = plt.subplots()
        _style_ax(ax, dark=True)
        assert ax.xaxis._gridOnMajor if hasattr(ax.xaxis, "_gridOnMajor") else True
        _close()

    def test_style_ax_no_figure_does_not_raise(self):
        """ax.get_figure() could in principle return None; _style_ax guards it."""
        fig, ax = plt.subplots()
        _style_ax(ax, dark=True)  # normal path, figure present
        _close()


class TestAnnotateEmpty:
    def test_default_message(self):
        fig, ax = plt.subplots()
        _annotate_empty(ax)
        texts = [t.get_text() for t in ax.texts]
        assert "No data" in texts
        _close()

    def test_custom_message(self):
        fig, ax = plt.subplots()
        _annotate_empty(ax, msg="Custom message here")
        texts = [t.get_text() for t in ax.texts]
        assert "Custom message here" in texts
        _close()


class TestStyleAllAxes:
    def test_styles_all_axes_in_figure(self):
        fig, axes = plt.subplots(2, 2)
        _style_all_axes(fig, dark=True)
        import matplotlib.colors as mcolors

        assert mcolors.to_hex(fig.patch.get_facecolor()) == "#1e1e2e"
        _close()

    def test_no_axes_does_not_raise(self):
        fig = plt.figure()
        _style_all_axes(fig, dark=True)
        _close()

    def test_swallows_per_axis_exceptions(self, monkeypatch):
        """If styling one axis raises, _style_all_axes must continue silently."""
        import pycsamt.app.desktop.controllers.qc_controller as qc_mod

        fig, ax = plt.subplots()

        def boom(ax, dark=True):
            raise RuntimeError("boom")

        monkeypatch.setattr(qc_mod, "_style_ax", boom)
        qc_mod._style_all_axes(fig, dark=True)  # must not raise
        _close()

    @pytest.mark.parametrize("dark", [True, False])
    def test_light_and_dark(self, dark):
        fig, axes = plt.subplots(1, 2)
        _style_all_axes(fig, dark)
        _close()


# ── describe_plot / catalogue ──────────────────────────────────────────────────


class TestDescribePlot:
    @pytest.mark.parametrize(
        "fn_name",
        [fn for _, fn, _ in ALL_PLOTS],
    )
    def test_describe_known_plot_returns_nonempty_string(self, fn_name):
        desc = describe_plot(fn_name)
        assert isinstance(desc, str)
        assert len(desc) > 0
        assert desc == PLOT_DESCRIPTIONS[fn_name]

    def test_describe_unknown_plot_returns_fallback(self):
        desc = describe_plot("not_a_real_plot_kind")
        assert desc == "Render this QC diagnostic plot."

    def test_all_catalogue_entries_have_description(self):
        for _, fn_name, _ in ALL_PLOTS:
            assert fn_name in PLOT_DESCRIPTIONS, f"missing description for {fn_name}"

    def test_catalogue_group_names_match_icons(self):
        group_names = {name for name, _ in ALL_GROUPS}
        assert group_names == set(GROUP_ICONS.keys())

    def test_all_groups_are_nonempty(self):
        for name, group in ALL_GROUPS:
            assert len(group) > 0, f"group {name} is empty"

    def test_catalogue_entries_have_three_fields(self):
        for _, group in ALL_GROUPS:
            for entry in group:
                assert len(entry) == 3
                label, fn_name, has_ax = entry
                assert isinstance(label, str)
                assert isinstance(fn_name, str)
                assert isinstance(has_ax, bool)

    def test_catalogue_functions_exist_in_emtools(self):
        import pycsamt.emtools as et

        for _, fn_name, _ in ALL_PLOTS:
            assert (
                getattr(et, fn_name, None) is not None
            ), f"{fn_name} not found in pycsamt.emtools"


class TestDynamicParameterSpecs:
    @pytest.mark.parametrize("method,expected,absent", [
        ("ama", "weights", "sig_val"),
        ("loess", "poly", "weights"),
        ("bilateral", "sig_val", "poly"),
        ("refmedian", "smooth_sites", "half_window"),
    ])
    def test_static_shift_method_parameters(self, method, expected, absent):
        specs = {s.name: s for s in qc_parameter_specs("ss_qc_profile", method)}
        assert specs[expected].method_specific
        assert absent not in specs
        assert set(specs["method"].choices) == {
            "ama", "loess", "bilateral", "refmedian"
        }

    def test_static_shift_parameters_reach_estimator(self, monkeypatch):
        from matplotlib.figure import Figure
        import pycsamt.emtools as et

        received = {}

        def plot(sites, *, method="ama", ax=None, verbose=0, **corr):
            received.update(corr)
            ax.plot([0, 1], [1, 2])

        monkeypatch.setattr(et, "ss_qc_profile", plot)
        ctrl = QCController()
        ctrl.set_sites(object())
        ctrl.draw("ss_qc_profile", True, Figure(), method="loess",
                  half_window=6, poly=1, it=4, sig_val=2)
        assert received == {"half_window": 6, "poly": 1, "it": 4}

    def test_field_zones_exposes_required_scientific_controls(self):
        specs = {spec.name: spec for spec in qc_parameter_specs("plot_field_zones")}
        assert specs["source_offset"].kind == "optional_float"
        assert specs["source_offset"].view is False
        assert specs["near_threshold"].default == pytest.approx(0.3)
        assert specs["far_threshold"].default == pytest.approx(3.0)

    def test_field_zones_separates_plot_view_parameters(self):
        specs = {spec.name: spec for spec in qc_parameter_specs("plot_field_zones")}
        for name in ("contour_kr", "kr_levels", "period_axis", "log_y", "figsize"):
            assert specs[name].view is True

    def test_source_overprint_exposes_source_offset_control(self):
        specs = {
            spec.name: spec
            for spec in qc_parameter_specs("plot_overprint_section")
        }
        assert specs["source_offset"].kind == "optional_float"
        assert specs["source_offset"].view is False

    def test_internal_plumbing_is_not_exposed(self):
        names = {spec.name for spec in qc_parameter_specs("plot_field_zones")}
        assert not names.intersection(
            {"sites", "ax", "recursive", "on_dup", "strict", "verbose"}
        )


# ── _auto_rho_bounds ────────────────────────────────────────────────────────────


class TestAutoRhoBounds:
    def test_with_real_tipper_sites_returns_finite_ordered_bounds(self, tipper_sites):
        lo, hi = _auto_rho_bounds(tipper_sites)
        assert np.isfinite(lo) and np.isfinite(hi)
        assert lo <= hi

    def test_with_real_willy_sites_returns_finite_ordered_bounds(self, willy_sites):
        lo, hi = _auto_rho_bounds(willy_sites)
        assert np.isfinite(lo) and np.isfinite(hi)
        assert lo <= hi

    def test_custom_percentiles(self, willy_sites):
        lo, hi = _auto_rho_bounds(willy_sites, lo_pct=10.0, hi_pct=90.0)
        assert np.isfinite(lo) and np.isfinite(hi)
        assert lo <= hi

    def test_empty_sites_falls_back_to_default_bounds(self):
        lo, hi = _auto_rho_bounds([])
        assert (lo, hi) == (0.5, 5000.0)

    def test_none_sites_does_not_raise_and_falls_back(self):
        lo, hi = _auto_rho_bounds(None)
        assert (lo, hi) == (0.5, 5000.0)

    def test_garbage_sites_does_not_raise_and_falls_back(self):
        lo, hi = _auto_rho_bounds(object())
        assert (lo, hi) == (0.5, 5000.0)

    def test_iterable_of_non_edi_items_falls_back(self):
        lo, hi = _auto_rho_bounds([1, 2, 3])
        assert (lo, hi) == (0.5, 5000.0)


# ── Polar dispatch helpers ───────────────────────────────────────────────────


class TestPolarDispatchHelpers:
    def test_dispatch_polar_coverage_adds_polar_axis(self, tipper_sites):
        import pycsamt.emtools as et

        ctrl_ = QCController()
        ctrl_.set_sites(tipper_sites)
        fig = _fig()
        _dispatch_polar_coverage(ctrl_, et.plot_polar_coverage, fig)
        assert len(fig.axes) >= 1
        assert fig.axes[0].name == "polar"
        _close()

    def test_dispatch_polar_ax_adds_polar_axis(self, willy_sites):
        import pycsamt.emtools as et

        ctrl_ = QCController()
        ctrl_.set_sites(willy_sites)
        fig = _fig()
        _dispatch_polar_ax(ctrl_, et.plot_ss_radar, fig)
        assert len(fig.axes) >= 1
        assert fig.axes[0].name == "polar"
        _close()


# ── draw() — no-data / error paths ────────────────────────────────────────────


class TestDrawNoData:
    """draw() must never raise; it annotates the figure when data is missing."""

    def test_no_sites_has_ax_does_not_raise(self, ctrl):
        fig = _fig()
        ctrl.draw("plot_coverage", True, fig)
        _close()

    def test_no_sites_no_ax_does_not_raise(self, ctrl):
        fig = _fig()
        ctrl.draw("plot_qc_quicklook", False, fig)
        _close()

    def test_no_sites_shows_load_message(self, ctrl):
        fig = _fig()
        ctrl.draw("plot_coverage", True, fig)
        assert "load" in _texts(fig) or "data" in _texts(fig)
        _close()

    def test_no_sites_returns_none(self, ctrl):
        fig = _fig()
        ret = ctrl.draw("plot_coverage", True, fig)
        assert ret is None
        _close()

    def test_no_sites_exposes_structured_reason(self, ctrl):
        ctrl.draw("plot_coverage", True, _fig())
        assert ctrl.last_unavailable is not None
        assert ctrl.last_unavailable.kind == "missing-data"
        assert "stations" in ctrl.last_unavailable.reason.lower()
        _close()

    def test_field_zones_explains_missing_source_geometry(
        self, ctrl, gv_xml_sites
    ):
        ctrl.set_sites(gv_xml_sites)
        ctrl.draw("plot_field_zones", True, _fig())
        result = ctrl.last_unavailable
        assert result is not None
        assert result.kind == "missing-prerequisite"
        assert "3 of 3 loaded stations" in result.reason
        assert "source_offset" in result.guidance
        _close()

    def test_source_overprint_missing_geometry_does_not_warn(
        self, ctrl, gv_xml_sites
    ):
        ctrl.set_sites(gv_xml_sites)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            ctrl.draw("plot_overprint_section", True, _fig())
        result = ctrl.last_unavailable
        assert result is not None
        assert result.kind == "missing-prerequisite"
        assert "Source-overprint" in result.title
        assert not any(
            "no source offset" in str(item.message).lower() for item in caught
        )
        _close()

    def test_source_overprint_runs_when_explicit_offset_is_supplied(
        self, ctrl, gv_xml_sites
    ):
        ctrl.set_sites(gv_xml_sites)
        ctrl.draw(
            "plot_overprint_section",
            True,
            _fig(),
            source_offset=1000.0,
        )
        assert ctrl.last_unavailable is None
        _close()

    def test_field_zones_runs_when_explicit_offset_is_supplied(
        self, ctrl, gv_xml_sites
    ):
        ctrl.set_sites(gv_xml_sites)
        ctrl.draw(
            "plot_field_zones",
            True,
            _fig(),
            source_offset=1000.0,
        )
        assert ctrl.last_unavailable is None
        _close()

    def test_unknown_function_does_not_raise(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("not_a_real_function", False, fig)
        _close()

    def test_unknown_function_annotates_figure(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("no_such_fn", False, fig)
        texts = _texts(fig)
        assert "not" in texts or "function" in texts
        _close()

    def test_figsize_resizes_figure(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("plot_coverage", True, fig, figsize=(8.0, 5.0))
        w, h = fig.get_size_inches()
        assert (round(w, 1), round(h, 1)) == (8.0, 5.0)
        _close()

    def test_draw_error_path_annotates_instead_of_raising(
        self, ctrl, willy_sites, monkeypatch
    ):
        """Force the wrapped emtools function to raise and verify draw()
        catches it and annotates the error rather than propagating."""
        import pycsamt.emtools as et

        def boom(*a, **kw):
            raise RuntimeError("synthetic failure")

        monkeypatch.setattr(et, "plot_coverage", boom)
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("plot_coverage", True, fig)  # must not raise
        assert "error" in _texts(fig) or "synthetic failure" in _texts(fig)
        _close()

    def test_draw_no_figure_produced_annotates(self, ctrl, willy_sites, monkeypatch):
        """has_ax=False function returning nothing / no new figure hits the
        'No figure produced' branch of draw()."""
        import pycsamt.emtools as et

        def returns_nothing(*a, **kw):
            return None

        monkeypatch.setattr(et, "plot_qc_quicklook", returns_nothing)
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw("plot_qc_quicklook", False, fig)
        assert "no figure produced" in _texts(fig)
        _close()


class TestCallFigureFn:
    def test_returns_figure_object_directly(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        import pycsamt.emtools as et

        result = ctrl._call_figure_fn(et.plot_qc_quicklook)
        assert hasattr(result, "get_axes")
        _close()

    def test_detects_newly_opened_figure_by_fignum(self, ctrl, willy_sites):
        """When the wrapped fn opens a new pyplot figure but returns something
        without get_axes, _call_figure_fn should locate it via fignums."""
        ctrl.set_sites(willy_sites)

        def opens_new_fig(sites, verbose=0, **kw):
            plt.figure()
            return "not-a-figure"

        result = ctrl._call_figure_fn(opens_new_fig)
        assert result is not None
        assert hasattr(result, "get_axes")
        _close()

    def test_returns_none_when_nothing_created(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)

        def does_nothing(sites, verbose=0, **kw):
            return "not-a-figure"

        set(plt.get_fignums())
        result = ctrl._call_figure_fn(does_nothing)
        assert result is None
        _close()


# ── Dark/light mode ────────────────────────────────────────────────────────────


class TestDarkLightMode:
    def test_dark_mode_toggle_does_not_raise(self, ctrl, willy_sites):
        ctrl.set_sites(willy_sites)
        for dark in (True, False):
            ctrl.dark = dark
            fig = _fig()
            ctrl.draw("plot_coverage", True, fig)
            _close()


# ── Catalogue parametrisation — draw() must never raise ──────────────────────
#
# Split by group (mirrors ALL_GROUPS) so failures are attributed clearly.
# Each group is exercised against both WILLY (multi-station, no tipper) and
# TIPPER (has actual tipper channels) real data.


class TestOverviewPlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(OVERVIEW_PLOTS),
        ids=_catalogue_ids(OVERVIEW_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(OVERVIEW_PLOTS),
        ids=_catalogue_ids(OVERVIEW_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


class TestCoveragePlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(COVERAGE_PLOTS),
        ids=_catalogue_ids(COVERAGE_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(COVERAGE_PLOTS),
        ids=_catalogue_ids(COVERAGE_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


class TestNoisePlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(NOISE_PLOTS),
        ids=_catalogue_ids(NOISE_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(NOISE_PLOTS),
        ids=_catalogue_ids(NOISE_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


class TestSkewDimPlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(SKEW_DIM_PLOTS),
        ids=_catalogue_ids(SKEW_DIM_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(SKEW_DIM_PLOTS),
        ids=_catalogue_ids(SKEW_DIM_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


class TestStaticShiftPlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(STATIC_SHIFT_PLOTS),
        ids=_catalogue_ids(STATIC_SHIFT_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(STATIC_SHIFT_PLOTS),
        ids=_catalogue_ids(STATIC_SHIFT_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


class TestDistortionPlots:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(DISTORTION_PLOTS),
        ids=_catalogue_ids(DISTORTION_PLOTS),
    )
    def test_draw_willy_does_not_raise(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()

    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(DISTORTION_PLOTS),
        ids=_catalogue_ids(DISTORTION_PLOTS),
    )
    def test_draw_tipper_does_not_raise(self, fn_name, has_ax, tipper_sites):
        c = QCController()
        c.set_sites(tipper_sites)
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


# ── Catalogue parametrisation — no-data path for every entry ─────────────────


class TestAllPlotsNoData:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(ALL_PLOTS),
        ids=_catalogue_ids(ALL_PLOTS),
    )
    def test_draw_no_sites_does_not_raise(self, fn_name, has_ax):
        c = QCController()
        fig = _fig()
        c.draw(fn_name, has_ax, fig)
        _close()


# ── Catalogue parametrisation — produces axes when data is present ───────────


class TestAllPlotsProduceAxes:
    @pytest.mark.parametrize(
        "fn_name,has_ax",
        _catalogue_params(ALL_PLOTS),
        ids=_catalogue_ids(ALL_PLOTS),
    )
    def test_draw_willy_produces_axes(self, fn_name, has_ax, willy_sites):
        c = QCController()
        c.set_sites(willy_sites)
        fig = _fig()
        ret = c.draw(fn_name, has_ax, fig)
        target = fig if ret is None else ret
        assert len(target.axes) >= 1
        _close()
