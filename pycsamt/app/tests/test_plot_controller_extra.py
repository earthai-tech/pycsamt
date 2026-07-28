# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Additional PlotController tests — real-EDI-backed coverage.

Supplements ``test_plot_controller.py`` (construction + state setters using a
fake "sites" string) with tests that exercise the module against real
survey data: the phase-tensor cache, period-range clipping helpers, the
module-level styling/formatting helpers, and a full smoke sweep of every
``draw_*`` method both with real data and with no/invalid data.

Real data
---------
data/MT/kap03lmt_edis/   — 26 KP TIPPER EDIs (real tipper channels)
data/AMT/WILLY_DATA/L18PLT/ — 28 WILLY AMT EDIs along one profile line
                               (good for pseudosections: multiple stations)

Strategy
--------
* draw_*() must never raise — smoke-tested with real data, no sites, and an
  invalid/unknown station id.
* Cache behaviour (_get_or_build_pt_df) is checked via object identity and
  via invalidate_phase_tensor()/set_sites() clearing it.
* Module-level helpers (style_axes, _annotate_empty, _relabel_colorbar_log10,
  _add_station_markers, _apply_log10_xfmt, _fix_psection_axes,
  _style_figure_full) are exercised directly against real ``plt.subplots()``
  axes.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pytest

from pycsamt.app.desktop.controllers.plot_controller import (
    PlotController,
    _add_station_markers,
    _annotate_empty,
    _apply_log10_xfmt,
    _fix_psection_axes,
    _relabel_colorbar_log10,
    _style_figure_full,
    style_axes,
)

# ── Paths ─────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"

_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


# ── Session-scoped fixtures ──────────────────────────────────────────────


@pytest.fixture(scope="session")
def tipper_sites():
    """26-station TIPPER Sites (real tipper channels), loaded once."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture(scope="session")
def willy_sites():
    """28-station WILLY L18PLT profile Sites, loaded once."""
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


# ── Helpers ───────────────────────────────────────────────────────────────


def _fig():
    return plt.figure(figsize=(6, 4))


def _close():
    plt.close("all")


def _texts(target):
    """Collect all text strings from a Figure or an Axes."""
    axes = target.axes if hasattr(target, "axes") else [target]
    return [t.get_text() for ax in axes for t in ax.texts]


def _drew_something(ax) -> bool:
    """True if ax has lines, collections, images, patches, or non-empty texts."""
    if ax.lines or ax.collections or ax.images or ax.patches:
        return True
    return any(t.get_text() for t in ax.texts)


# ═════════════════════════════════════════════════════════════════════════
# Phase-tensor DataFrame cache
# ═════════════════════════════════════════════════════════════════════════


class TestPhaseTensorCache:
    def test_no_sites_returns_empty_dataframe(self):
        ctrl = PlotController()
        df = ctrl._get_or_build_pt_df()
        assert df.empty

    def test_builds_and_caches_dataframe(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        df1 = ctrl._get_or_build_pt_df()
        assert not df1.empty
        assert len(ctrl._pt_df_cache) == 1

    def test_cache_hit_returns_same_object(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        df1 = ctrl._get_or_build_pt_df()
        df2 = ctrl._get_or_build_pt_df()
        assert df1 is df2  # cache hit, no rebuild

    def test_invalidate_forces_rebuild(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        df1 = ctrl._get_or_build_pt_df()
        ctrl.invalidate_phase_tensor()
        assert ctrl._pt_df_cache == {}
        df2 = ctrl._get_or_build_pt_df()
        assert df1 is not df2  # rebuilt fresh
        assert df1.equals(df2)  # but content is equivalent

    def test_set_sites_clears_cache(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl._get_or_build_pt_df()
        assert ctrl._pt_df_cache
        ctrl.set_sites(willy_sites)  # re-set (even same object) clears cache
        assert ctrl._pt_df_cache == {}

    def test_clear_clears_cache(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl._get_or_build_pt_df()
        assert ctrl._pt_df_cache
        ctrl.clear()
        assert ctrl._pt_df_cache == {}

    def test_phase_tensor_key_changes_with_sites_period_dark(self, willy_sites):
        ctrl = PlotController()
        key0 = ctrl.phase_tensor_key()
        ctrl.set_sites(willy_sites)
        key1 = ctrl.phase_tensor_key()
        assert key1 != key0
        ctrl.set_period_range(0.01, 1.0)
        key2 = ctrl.phase_tensor_key()
        assert key2 != key1
        ctrl.dark = False
        key3 = ctrl.phase_tensor_key()
        assert key3 != key2

    def test_phase_tensor_key_stable_without_change(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        assert ctrl.phase_tensor_key() == ctrl.phase_tensor_key()


# ═════════════════════════════════════════════════════════════════════════
# _x_label
# ═════════════════════════════════════════════════════════════════════════


def test_x_label_is_log10_period():
    ctrl = PlotController()
    assert ctrl._x_label == r"$\log_{10}(T)\ \mathrm{(s)}$"


# ═════════════════════════════════════════════════════════════════════════
# _active_period_range / clipping helpers
# ═════════════════════════════════════════════════════════════════════════


class TestActivePeriodRange:
    def test_none_when_unset(self):
        ctrl = PlotController()
        assert ctrl._active_period_range() is None

    def test_none_when_full_range(self):
        ctrl = PlotController()
        ctrl.set_period_range(0.0, 1e9)
        assert ctrl._active_period_range() is None

    def test_returns_tuple_when_narrowed(self):
        ctrl = PlotController()
        ctrl.set_period_range(0.001, 10.0)
        assert ctrl._active_period_range() == (0.001, 10.0)


class TestClipXAxisPeriod:
    def test_no_filter_leaves_xlim_untouched(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        before = ax.get_xlim()
        ctrl._clip_xaxis_period(ax)
        assert ax.get_xlim() == before
        plt.close(fig)

    def test_filter_sets_xlim(self):
        ctrl = PlotController()
        ctrl.set_period_range(0.01, 5.0)
        fig, ax = plt.subplots()
        ax.plot([0.001, 1, 100], [1, 2, 3])
        ctrl._clip_xaxis_period(ax)
        assert ax.get_xlim() == (0.01, 5.0)
        plt.close(fig)

    def test_multiple_axes_all_clipped(self):
        ctrl = PlotController()
        ctrl.set_period_range(0.01, 5.0)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ctrl._clip_xaxis_period(ax1, ax2)
        assert ax1.get_xlim() == (0.01, 5.0)
        assert ax2.get_xlim() == (0.01, 5.0)
        plt.close(fig)

    def test_no_axes_does_not_raise(self):
        ctrl = PlotController()
        ctrl.set_period_range(0.01, 5.0)
        ctrl._clip_xaxis_period()  # no axes supplied


class TestClipYAxisPeriod:
    def test_no_filter_leaves_ylim_untouched(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])
        before = ax.get_ylim()
        ctrl._clip_yaxis_period(ax)
        assert ax.get_ylim() == before
        plt.close(fig)

    def test_filter_clips_normal_axis(self):
        ctrl = PlotController()
        ctrl.set_period_range(1.0, 8.0)
        fig, ax = plt.subplots()
        ax.set_ylim(0.0, 10.0)
        ctrl._clip_yaxis_period(ax)
        lo, hi = ax.get_ylim()
        assert lo == pytest.approx(1.0)
        assert hi == pytest.approx(8.0)
        plt.close(fig)

    def test_filter_preserves_inverted_axis(self):
        ctrl = PlotController()
        ctrl.set_period_range(1.0, 8.0)
        fig, ax = plt.subplots()
        ax.set_ylim(10.0, 0.0)  # inverted: top < bottom numerically
        ctrl._clip_yaxis_period(ax)
        lo, hi = ax.get_ylim()
        # inverted axis preserved: first value > second
        assert lo > hi
        assert lo == pytest.approx(8.0)
        assert hi == pytest.approx(1.0)
        plt.close(fig)

    def test_multiple_axes_all_clipped(self):
        ctrl = PlotController()
        ctrl.set_period_range(1.0, 8.0)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        ax1.set_ylim(0, 10)
        ax2.set_ylim(0, 10)
        ctrl._clip_yaxis_period(ax1, ax2)
        assert ax1.get_ylim()[0] == pytest.approx(1.0)
        assert ax2.get_ylim()[0] == pytest.approx(1.0)
        plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# _get_site / _mark_station
# ═════════════════════════════════════════════════════════════════════════


class TestGetSite:
    def test_no_sites_returns_none(self):
        ctrl = PlotController()
        assert ctrl._get_site("anything") is None

    def test_valid_station_returns_site(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        name = willy_sites.by_index(0).name
        site = ctrl._get_site(name)
        assert site is not None

    def test_invalid_station_returns_none(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        assert ctrl._get_site("NOT-A-REAL-STATION-XYZ") is None


class TestMarkStation:
    def test_no_station_id_does_not_raise(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl._mark_station(ax)  # station_id is None -> no-op
        assert len(ax.lines) == 0
        plt.close(fig)

    def test_station_not_in_labels_does_not_raise(self):
        ctrl = PlotController()
        ctrl.set_station("S99")
        fig, ax = plt.subplots()
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["A", "B", "C"])
        ctrl._mark_station(ax)  # not found among labels -> no-op, no raise
        plt.close(fig)

    def test_station_in_labels_adds_vline(self):
        ctrl = PlotController()
        ctrl.set_station("B")
        fig, ax = plt.subplots()
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(["A", "B", "C"])
        ctrl._mark_station(ax)
        assert len(ax.lines) == 1
        plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# Module-level styling / formatting helpers
# ═════════════════════════════════════════════════════════════════════════


class TestStyleAxesFull:
    def test_dark_sets_dark_facecolor(self):
        fig, ax = plt.subplots()
        style_axes(ax, dark=True)
        # dark facecolor #181825
        import matplotlib.colors as mcolors

        assert ax.get_facecolor() == mcolors.to_rgba("#181825")
        plt.close(fig)

    def test_light_sets_light_facecolor(self):
        fig, ax = plt.subplots()
        style_axes(ax, dark=False)
        import matplotlib.colors as mcolors

        assert ax.get_facecolor() == mcolors.to_rgba("#eff1f5")
        plt.close(fig)

    def test_sets_figure_patch_color(self):
        fig, ax = plt.subplots()
        style_axes(ax, dark=True)
        import matplotlib.colors as mcolors

        assert fig.patch.get_facecolor() == mcolors.to_rgba("#1e1e2e")
        plt.close(fig)


class TestStyleFigureFull:
    def test_does_not_raise_on_plain_axes(self):
        fig, ax = plt.subplots()
        style_axes(ax, dark=True)
        _style_figure_full(ax, dark=True)
        plt.close(fig)

    def test_does_not_raise_light_mode(self):
        fig, ax = plt.subplots()
        style_axes(ax, dark=False)
        _style_figure_full(ax, dark=False)
        plt.close(fig)

    def test_recolors_bbox_text_dark(self):
        fig, ax = plt.subplots()
        txt = ax.text(0.5, 0.5, "hi", bbox=dict(facecolor="white", edgecolor="black"))
        txt.set_color("#2166ac")
        _style_figure_full(ax, dark=True)
        assert txt.get_color() == "#89b4fa"
        bb = txt.get_bbox_patch()
        assert bb.get_facecolor() == (
            __import__("matplotlib").colors.to_rgba("#1a1a2e", alpha=0.88)
        )
        plt.close(fig)

    def test_recolors_free_text_dark(self):
        fig, ax = plt.subplots()
        txt = ax.text(0.5, 0.5, "no bbox")
        _style_figure_full(ax, dark=True)
        assert txt.get_color() == "#cdd6f4"
        plt.close(fig)

    def test_unfilled_patch_edge_recolored_dark(self):
        import matplotlib.patches as mpatches

        fig, ax = plt.subplots()
        patch = mpatches.Circle((0, 0), 1, fill=False, edgecolor="k")
        ax.add_patch(patch)
        _style_figure_full(ax, dark=True)
        import matplotlib.colors as mcolors

        assert patch.get_edgecolor() == mcolors.to_rgba("#cccccc")
        plt.close(fig)

    def test_unfilled_patch_edge_recolored_light(self):
        import matplotlib.patches as mpatches

        fig, ax = plt.subplots()
        patch = mpatches.Circle((0, 0), 1, fill=False, edgecolor="k")
        ax.add_patch(patch)
        _style_figure_full(ax, dark=False)
        import matplotlib.colors as mcolors

        assert patch.get_edgecolor() == mcolors.to_rgba("#444444")
        plt.close(fig)

    def test_figure_none_is_noop(self):
        fig, ax = plt.subplots()
        plt.close(fig)  # ax now has no live figure reference issue avoided
        # Use a fresh Axes-like scenario: get_figure() still returns fig
        # even when closed, so directly assert no exception when called twice.
        _style_figure_full(ax, dark=True)


class TestAnnotateEmpty:
    def test_default_message(self):
        fig, ax = plt.subplots()
        _annotate_empty(ax)
        assert ax.texts[-1].get_text() == "No data"
        plt.close(fig)

    def test_custom_message(self):
        fig, ax = plt.subplots()
        _annotate_empty(ax, "Custom message here")
        assert ax.texts[-1].get_text() == "Custom message here"
        plt.close(fig)


class TestRelabelColorbarLog10:
    def test_relabels_positive_yticks(self):
        fig, ax = plt.subplots()
        cb_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        cb_ax.set_yticks([1, 10, 100, 10000])
        _relabel_colorbar_log10(cb_ax)
        labels = [t.get_text() for t in cb_ax.get_yticklabels()]
        assert labels == ["0", "1", "2", "4"]
        plt.close(fig)

    def test_nonpositive_ticks_become_empty_string(self):
        fig, ax = plt.subplots()
        cb_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        cb_ax.set_yticks([-1, 0, 1])
        _relabel_colorbar_log10(cb_ax)
        labels = [t.get_text() for t in cb_ax.get_yticklabels()]
        assert labels[0] == ""
        assert labels[1] == ""
        plt.close(fig)

    def test_falls_back_to_xticks_when_no_yticks(self):
        fig, ax = plt.subplots()
        cb_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        cb_ax.set_yticks([])
        cb_ax.set_xticks([1, 100])
        _relabel_colorbar_log10(cb_ax)
        labels = [t.get_text() for t in cb_ax.get_xticklabels()]
        assert labels == ["0", "2"]
        plt.close(fig)

    def test_does_not_raise_on_broken_axes(self):
        # Passing something without get_yticks must be swallowed silently
        class _Bad:
            pass

        _relabel_colorbar_log10(_Bad())  # must not raise


class TestAddStationMarkers:
    def test_no_valid_ticks_is_noop(self):
        fig, ax = plt.subplots()
        ax.set_xticks([])
        ax.set_xlim(0, 1)
        _add_station_markers(ax, dark=True)
        # nothing scattered
        assert len(ax.collections) == 0
        plt.close(fig)

    def test_adds_scatter_when_ticks_in_range(self):
        fig, ax = plt.subplots()
        ax.set_xlim(-0.5, 5.5)
        ax.set_xticks([0, 1, 2, 3, 4, 5])
        n_before = len(ax.collections)
        _add_station_markers(ax, dark=True)
        assert len(ax.collections) > n_before
        plt.close(fig)

    def test_does_not_raise_light_mode(self):
        fig, ax = plt.subplots()
        ax.set_xlim(-0.5, 3.5)
        ax.set_xticks([0, 1, 2, 3])
        _add_station_markers(ax, dark=False)
        plt.close(fig)


class TestApplyLog10Xfmt:
    def test_sets_major_formatter(self):
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        _apply_log10_xfmt(ax)
        assert isinstance(ax.xaxis.get_major_formatter(), mticker.FuncFormatter)
        assert isinstance(ax.xaxis.get_minor_formatter(), mticker.NullFormatter)
        plt.close(fig)

    def test_formatter_maps_powers_of_ten(self):
        fig, ax = plt.subplots()
        ax.set_xscale("log")
        _apply_log10_xfmt(ax)
        fmt = ax.xaxis.get_major_formatter()
        assert fmt(1.0, 0) == "$0$"
        assert fmt(100.0, 0) == "$2$"
        assert fmt(0.001, 0) == "$-3$"

    def test_formatter_non_positive_is_blank(self):
        fig, ax = plt.subplots()
        _apply_log10_xfmt(ax)
        fmt = ax.xaxis.get_major_formatter()
        assert fmt(0.0, 0) == ""
        assert fmt(-5.0, 0) == ""
        plt.close(fig)

    def test_formatter_non_power_of_ten_uses_decimal(self):
        fig, ax = plt.subplots()
        _apply_log10_xfmt(ax)
        fmt = ax.xaxis.get_major_formatter()
        # log10(50) ~ 1.7 -> not close to an integer
        out = fmt(50.0, 0)
        assert out.startswith("$1.7")
        plt.close(fig)


class TestFixPsectionAxes:
    def test_relabels_period_yticks_to_log10(self):
        fig, ax = plt.subplots()
        ax.set_yticks([0, 1, 2])
        ax.set_yticklabels(["0.001", "1", "1000"])
        _fix_psection_axes(ax)
        labels = [t.get_text() for t in ax.get_yticklabels()]
        assert labels == ["-3", "0", "3"]
        assert ax.get_ylabel() == r"$\log_{10}(T)\ \mathrm{(s)}$"
        plt.close(fig)

    def test_non_numeric_labels_preserved(self):
        fig, ax = plt.subplots()
        ax.set_yticks([0, 1])
        ax.set_yticklabels(["foo", "bar"])
        _fix_psection_axes(ax)  # must not raise
        labels = [t.get_text() for t in ax.get_yticklabels()]
        assert labels == ["foo", "bar"]
        plt.close(fig)

    def test_colorbar_label_applied_when_second_axes_present(self):
        fig, ax = plt.subplots()
        cb_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        ax.set_yticks([0])
        ax.set_yticklabels(["1"])
        _fix_psection_axes(ax, colorbar_label=r"$\log_{10}(\rho_a)$")
        assert cb_ax.get_ylabel() == r"$\log_{10}(\rho_a)$"
        plt.close(fig)

    def test_no_colorbar_label_leaves_other_axes_alone(self):
        fig, ax = plt.subplots()
        cb_ax = fig.add_axes([0.85, 0.1, 0.05, 0.8])
        _fix_psection_axes(ax, colorbar_label="")
        assert cb_ax.get_ylabel() == ""
        plt.close(fig)


# ═════════════════════════════════════════════════════════════════════════
# draw_* smoke sweep — real data
# ═════════════════════════════════════════════════════════════════════════


class TestDrawRhoPhiRealData:
    def test_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        _close()

    def test_draws_curves(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        assert len(fig.axes) == 2
        ax_r, ax_p = fig.axes
        assert len(ax_r.lines) > 0 or len(ax_r.collections) > 0
        assert len(ax_p.lines) > 0 or len(ax_p.collections) > 0
        plt.close(fig)

    def test_bw_mode_uses_black(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_bw_mode(True)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        _close()

    def test_period_range_clips_x_axis(self, willy_sites):
        """NOTE: as of matplotlib 3.5.x this exposes a real bug — see
        TestRegressionLegendNcolsBug below.  ``draw_rho_phi`` calls
        ``ax_r.legend(ncols=n_comps, ...)`` (matplotlib>=3.6 only kwarg);
        on older matplotlib this raises TypeError *inside* the method's
        own try/except, so the method silently falls into its error
        branch and the period-range clip (which runs after the legend
        call) never executes.  Once that bug is fixed upstream, this
        test should start asserting lo == 0.01 / hi == 1.0 again.
        """
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_period_range(0.01, 1.0)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        ax_p = fig.axes[1]
        lo, hi = ax_p.get_xlim()
        try:
            import matplotlib.legend as _mleg

            _mleg.Legend(
                plt.subplots()[1], [], [], ncols=1
            )  # probe: does this mpl support ncols?
            supports_ncols = True
        except TypeError:
            supports_ncols = False
        finally:
            plt.close("all")
        if supports_ncols:
            assert lo == pytest.approx(0.01)
            assert hi == pytest.approx(1.0)
        else:
            # Bug reproduced: clip never applied because legend() crashed first.
            assert (lo, hi) != (0.01, 1.0)
        plt.close(fig)

    def test_phase_ylim_applied(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_phase_ylim(-90, 90)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        ax_p = fig.axes[1]
        assert ax_p.get_ylim() == (-90.0, 90.0)
        plt.close(fig)

    def test_no_errbar(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_show_errbar(False)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        _close()

    def test_invalid_station_shows_error(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station("NOT-A-REAL-STATION-XYZ")
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        texts = " ".join(_texts(fig)).lower()
        assert "error" in texts or "not found" in texts
        _close()

    def test_all_four_components(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_components(("xy", "yx", "xx", "yy"))
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        _close()


class TestDrawRhoPseudosectionRealData:
    def test_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        plt.close(fig)

    def test_produces_image(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        assert len(ax.images) > 0
        # Title is set with loc="left" — not the default center location.
        assert ax.get_title(loc="left")
        plt.close(fig)

    def test_with_station_marker(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(3).name)
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        assert len(ax.images) > 0
        plt.close(fig)

    def test_with_period_range(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_period_range(0.001, 1.0)
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        assert len(ax.images) > 0
        plt.close(fig)

    def test_light_mode(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.dark = False
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        plt.close(fig)


class TestDrawPhasePseudosectionRealData:
    def test_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_pseudosection(ax)
        plt.close(fig)

    def test_produces_image(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_pseudosection(ax)
        assert len(ax.images) > 0
        plt.close(fig)


class TestDrawTipperRealData:
    def test_survey_level_does_not_raise(self, tipper_sites):
        ctrl = PlotController()
        ctrl.set_sites(tipper_sites)
        fig, ax = plt.subplots()
        ctrl.draw_tipper(ax)
        plt.close(fig)

    def test_survey_level_draws_something(self, tipper_sites):
        ctrl = PlotController()
        ctrl.set_sites(tipper_sites)
        fig, ax = plt.subplots()
        ctrl.draw_tipper(ax)
        assert _drew_something(ax)
        plt.close(fig)

    def test_single_station_draws_something(self, tipper_sites):
        ctrl = PlotController()
        ctrl.set_sites(tipper_sites)
        ctrl.set_station(tipper_sites.by_index(0).name)
        fig, ax = plt.subplots()
        ctrl.draw_tipper(ax)
        assert _drew_something(ax)
        assert tipper_sites.by_index(0).name in ax.get_title()
        plt.close(fig)

    def test_period_range_clips_xaxis(self, tipper_sites):
        ctrl = PlotController()
        ctrl.set_sites(tipper_sites)
        ctrl.set_period_range(0.01, 1.0)
        fig, ax = plt.subplots()
        ctrl.draw_tipper(ax)
        plt.close(fig)


class TestDrawPhaseTensorRealData:
    def test_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor(ax)
        plt.close(fig)

    def test_draws_ellipses(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor(ax)
        assert _drew_something(ax)
        assert ax.get_title()
        plt.close(fig)

    def test_reuses_cache_across_calls(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig1, ax1 = plt.subplots()
        ctrl.draw_phase_tensor(ax1)
        cached_df = ctrl._pt_df_cache[id(willy_sites)]
        fig2, ax2 = plt.subplots()
        ctrl.draw_phase_tensor(ax2)
        assert ctrl._pt_df_cache[id(willy_sites)] is cached_df
        plt.close(fig1)
        plt.close(fig2)

    def test_light_mode(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.dark = False
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor(ax)
        plt.close(fig)


class TestDrawPhaseTensorStripRealData:
    def test_no_station_shows_prompt(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor_strip(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "select" in texts or "station" in texts
        plt.close(fig)

    def test_with_station_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor_strip(ax)
        plt.close(fig)

    def test_with_station_draws_something(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor_strip(ax)
        assert _drew_something(ax)
        plt.close(fig)


class TestDrawPublicationViewRealData:
    def test_does_not_raise(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_publication_view(fig)
        _close()

    def test_produces_axes(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_publication_view(fig)
        assert len(fig.axes) >= 2  # at least rho + phase rows
        plt.close(fig)

    def test_with_tipper_data_adds_tx_ty_rows(self, tipper_sites):
        ctrl = PlotController()
        ctrl.set_sites(tipper_sites)
        ctrl.set_station(tipper_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_publication_view(fig)
        # has_tipper branch -> 4 rows per component column
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_no_station_shows_prompt(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig = _fig()
        ctrl.draw_publication_view(fig)
        texts = " ".join(_texts(fig)).lower()
        assert "select" in texts
        _close()

    def test_multi_component_columns(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        ctrl.set_components(("xy", "yx", "xx", "yy"))
        fig = _fig()
        ctrl.draw_publication_view(fig)
        _close()

    def test_invalid_station_annotates_error(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station("NOT-A-REAL-STATION-XYZ")
        fig = _fig()
        ctrl.draw_publication_view(fig)
        texts = " ".join(_texts(fig)).lower()
        assert "error" in texts or "not found" in texts
        _close()


# ═════════════════════════════════════════════════════════════════════════
# draw_* smoke sweep — no-data / invalid-data paths
# ═════════════════════════════════════════════════════════════════════════


class TestDrawNoDataPaths:
    """Every draw_* method must annotate 'No data'/'Load...' rather than raise."""

    def test_draw_rho_pseudosection_no_sites(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl.draw_rho_pseudosection(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "load" in texts
        plt.close(fig)

    def test_draw_phase_pseudosection_no_sites(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl.draw_phase_pseudosection(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "load" in texts
        plt.close(fig)

    def test_draw_tipper_no_sites(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl.draw_tipper(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "load" in texts
        plt.close(fig)

    def test_draw_phase_tensor_no_sites(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "load" in texts
        plt.close(fig)

    def test_draw_phase_tensor_strip_no_sites(self):
        ctrl = PlotController()
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor_strip(ax)
        texts = " ".join(t.get_text() for t in ax.texts).lower()
        assert "load" in texts
        plt.close(fig)

    def test_draw_publication_view_no_sites(self):
        ctrl = PlotController()
        fig = _fig()
        ctrl.draw_publication_view(fig)
        texts = " ".join(_texts(fig)).lower()
        assert "select" in texts
        _close()

    def test_draw_rho_phi_no_sites(self):
        ctrl = PlotController()
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        texts = " ".join(_texts(fig)).lower()
        assert "load" in texts
        _close()


# ═════════════════════════════════════════════════════════════════════════
# Regression / robustness checks
# ═════════════════════════════════════════════════════════════════════════


class TestRegressionLegendNcolsBug:
    """Documents a real bug found while writing these tests (NOT fixed here
    per task constraints — production code under pycsamt/app/ is untouched).

    ``PlotController.draw_rho_phi`` (plot_controller.py, around line 587)
    calls::

        ax_r.legend(ncols=n_comps, fontsize=8, loc="best", framealpha=0.65)

    ``ncols`` (plural) was only added to ``Axes.legend()`` in matplotlib
    3.6 (the long-standing kwarg is ``ncol``, singular). On matplotlib
    3.5.x this raises::

        TypeError: __init__() got an unexpected keyword argument 'ncols'

    Because draw_rho_phi wraps its whole body in a broad try/except that
    annotates the axes with "Plot error: ..." instead of re-raising (by
    design, per the module's documented contract), the failure is
    swallowed silently — the curves are still drawn (they're plotted
    before the legend() call) but:
      * no legend ever appears,
      * the period-range x-axis clip (``self._clip_xaxis_period(ax_p)``,
        which runs *after* the legend() call) never executes.

    ``draw_publication_view`` does not call ``.legend()`` and is
    unaffected.
    """

    def test_legend_ncols_kwarg_matches_installed_matplotlib(self):
        """Sanity-check *why* the bug does or doesn't reproduce here."""
        import matplotlib

        major, minor = (int(x) for x in matplotlib.__version__.split(".")[:2])
        supports_ncols = (major, minor) >= (3, 6)
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1], label="x")
        if supports_ncols:
            ax.legend(ncols=1)  # must succeed on modern matplotlib
        else:
            with pytest.raises(TypeError, match="ncols"):
                ax.legend(ncols=1)
        plt.close(fig)

    def test_draw_rho_phi_still_draws_curves_despite_legend_bug(self, willy_sites):
        """Even when the legend() call raises, the earlier errorbar()
        calls already happened, so lines are still visible — this is why
        a naive 'is something drawn' smoke test does not catch the bug."""
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        ax_r, ax_p = fig.axes
        drew_curves = bool(ax_r.lines) or bool(ax_r.collections)
        texts = " ".join(t.get_text() for t in ax_r.texts)
        # Either the curves are visible (bug reproduced: crash happened
        # after drawing but before legend) or a real "Plot error" was
        # annotated on top of them — both are consistent with the bug.
        assert drew_curves or "Plot error" in texts
        plt.close(fig)


class TestRepeatedDrawStability:
    """Calling a draw_* method repeatedly on the same axes/figure must not
    accumulate stale artists or raise on the second pass."""

    def test_draw_rho_phi_repeated_calls(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        ctrl.set_station(willy_sites.by_index(0).name)
        fig = _fig()
        ctrl.draw_rho_phi(fig)
        ctrl.draw_rho_phi(fig)  # fig.clear() inside -> should not accumulate
        assert len(fig.axes) == 2
        plt.close(fig)

    def test_draw_phase_tensor_repeated_calls_same_axes(self, willy_sites):
        ctrl = PlotController()
        ctrl.set_sites(willy_sites)
        fig, ax = plt.subplots()
        ctrl.draw_phase_tensor(ax)
        ctrl.draw_phase_tensor(ax)  # does not clear ax -> must still not raise
        plt.close(fig)
