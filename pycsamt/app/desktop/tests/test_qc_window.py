# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for the QC studio (controllers/qc_studio + windows/qc_window).

Real data: two profiles (L18, L22) of the bundled Baohuashan CSAMT survey
(``data/AMT/WILLY_DATA``, Kouabena 2025).
"""

from __future__ import annotations

import glob
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers import qc_studio as st
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.controllers.inversion_engines import defaults

_ROOT = Path(__file__).resolve().parents[4]
_WILLY = _ROOT / "data" / "AMT" / "WILLY_DATA"
pytestmark = pytest.mark.skipif(not (_WILLY / "L22PLT").exists(),
                                reason="Baohuashan data missing")


@pytest.fixture(scope="module")
def sites():
    from pycsamt.site.base import to_sites

    files = sorted(glob.glob(str(_WILLY / "L18PLT" / "*.edi"))
                   + glob.glob(str(_WILLY / "L22PLT" / "*.edi")))
    return to_sites(files)


@pytest.fixture(scope="module")
def lines(sites):
    return {n: "L" + n.split("-")[0] for n in st.station_names(sites)}


@pytest.fixture
def win(qapp, sites, lines):
    from pycsamt.app.desktop.windows.qc_window import QCDashboardWindow

    w = QCDashboardWindow(parent=None)
    w.show()
    w.set_lines(lines)
    w.set_sites(sites)
    yield w
    w.close()
    plt.close("all")


def _render(v, sites, values=None, **kw):
    vals = dict(defaults(st.options_for(v)), **(values or {}))
    kwargs, err = st.to_kwargs(v, vals)
    assert err is None, err
    return st.render(v, sites, kwargs, **kw)


# ── catalogue & options ───────────────────────────────────────────────────
class TestCatalogue:
    def test_every_view_exists_in_emtools(self):
        import pycsamt.emtools as et

        assert all(callable(getattr(et, v.fn, None)) for v in st.VIEWS)
        assert {v.group for v in st.VIEWS} == set(st.GROUPS)
        assert "Overview" not in st.GROUPS

    def test_colour_maps_offer_jet_and_jet_r(self):
        v = st.view("plot_snr_section")
        cmap = next(f for f in st.options_for(v) if f.key == "cmap")
        values = [c for c, _l in cmap.choices]
        assert {"jet", "jet_r", "RdYlBu_r"} <= set(values)
        assert cmap.default == values[0] == "RdYlGn"  # its own default first

    def test_confidence_method_offers_only_valid_criteria(self):
        for fn in ("plot_confidence_profile", "plot_confidence_band_summary",
                   "plot_frequency_confidence_psection"):
            f = next(f for f in st.options_for(st.view(fn))
                     if f.key == "method")
            assert {c for c, _l in f.choices} == {"presence", "composite"}

    def test_fixed_string_options_are_dropdowns(self):
        fields = {f.key: f for f in
                  st.options_for(st.view("plot_confidence_grid_map"))}
        assert fields["interpolation"].kind == "choice"
        assert {c for c, _ in fields["interpolation"].choices} == {
            "linear", "cubic"}
        assert fields["map_aspect"].kind == "choice"

    def test_cosmetic_options_fold_under_more_settings(self):
        fields = st.options_for(st.view("plot_confidence_distribution"))
        colour = [f for f in fields if f.key.endswith("color")]
        assert colour and all(f.advanced and f.section == "View"
                              for f in colour)

    def test_to_kwargs_parses_and_validates(self):
        v = st.view("plot_confidence_profile")
        kw, err = st.to_kwargs(v, {"ylim": "0.2, 1.1", "annotate_low_step":
                                   "", "ci_hi": 0.9, "ci_lo": 0.8})
        assert err is None and kw["ylim"] == (0.2, 1.1)
        assert "annotate_low_step" not in kw  # blank = function default
        _kw, err = st.to_kwargs(v, {"ci_hi": 0.8, "ci_lo": 0.9})
        assert "below" in err
        _kw, err = st.to_kwargs(v, {"annotate_low_step": "x"})
        assert "number" in err

    def test_bool_or_auto_options(self):
        v = st.view("plot_confidence_heatmap")
        kw, _ = st.to_kwargs(v, {"annotate": "false"})
        assert kw["annotate"] is False
        kw, _ = st.to_kwargs(v, {"annotate": "auto"})
        assert kw["annotate"] == "auto"


# ── rendering ─────────────────────────────────────────────────────────────
class TestRender:
    @pytest.mark.parametrize("v", [v for v in st.VIEWS if not v.fn in (
        "plot_overprint_section", "plot_field_zones",
        "nr_qc_harmonic_waterfall")], ids=lambda v: v.fn)
    def test_every_view_draws_one_line(self, sites, lines, v):
        fig, un = _render(v, sites, lines=lines, line="L22",
                          station="22-5U" if v.station else "")
        assert un is None, un and un.reason
        assert figure_blank_reason(fig) is None
        plt.close("all")

    def test_source_geometry_plots_explain_or_draw(self, sites, lines):
        v = st.view("plot_field_zones")
        _fig, un = _render(v, sites, lines=lines, line="L22")
        assert un is not None and "source" in un.reason
        fig, un = _render(v, sites, {"source_offset": "5000"}, lines=lines,
                          line="L22")
        assert un is None and figure_blank_reason(fig) is None
        plt.close("all")

    def test_specific_reason_passed_through(self, sites):
        _fig, un = _render(st.view("nr_qc_harmonic_waterfall"), sites)
        assert "harmonic" in un.reason  # the plot's own sentence
        plt.close("all")

    def test_panels_one_per_line(self, sites, lines):
        fig, un = _render(st.view("plot_confidence_profile"), sites,
                          lines=lines, layout="panels")
        assert un is None
        tags = [t.get_text() for ax in fig.axes for t in ax.texts
                if t.get_text().startswith("Line ")]
        assert sorted(tags) == ["Line L18", "Line L22"]
        assert fig._suptitle is not None
        plt.close("all")

    def test_together_colours_each_line(self, sites, lines):
        fig, _un = _render(st.view("plot_confidence_profile"), sites,
                           lines=lines, layout="together")
        legend = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
        assert "L18" in legend and "L22" in legend
        plt.close("all")

    def test_strike_rose_by_line_gets_groups(self, sites, lines):
        fig, un = _render(st.view("plot_strike_rose_by_line"), sites,
                          lines=lines)
        assert un is None and figure_blank_reason(fig) is None
        plt.close("all")


# ── summary ───────────────────────────────────────────────────────────────
class TestScorecard:
    def test_columns_and_status(self, sites, lines):
        df = st.scorecard(sites, lines)
        assert len(df) == len(st.station_names(sites))
        assert {"Line", "Station", "Status", "Confidence ratio",
                "Composite score", "Coverage", "SNR (median)",
                "Flags"} <= set(df.columns)
        assert set(df["Status"]) <= {"pass", "warn", "fail"}
        assert set(df["Line"]) == {"L18", "L22"}

    def test_thresholds_drive_status(self, sites):
        strict = st.scorecard(sites, ci_lo=1.01, ci_hi=1.02)
        assert set(strict["Status"]) == {"fail"}
        noisy = st.scorecard(sites, min_snr=1e6)
        assert "pass" not in set(noisy["Status"])
        assert "Line" not in noisy.columns  # no line labels given


# ── window ────────────────────────────────────────────────────────────────
class TestWindow:
    def test_empty_window_explains(self, qapp):
        from pycsamt.app.desktop.windows.qc_window import QCDashboardWindow

        w = QCDashboardWindow()
        try:
            assert not w.has_data
            assert "No survey" in w._data_status.text()
        finally:
            w.close()

    def test_scope_lists_lines_and_stations(self, win, sites):
        assert "2 lines" in win._data_status.text()
        items = [win._line_combo.itemData(i)
                 for i in range(win._line_combo.count())]
        assert items == ["", "L18", "L22"]
        win._line_combo.setCurrentIndex(1)
        assert all(win._station_combo.itemData(i).startswith("18-")
                   for i in range(win._station_combo.count()))

    def test_scope_rows_follow_the_view(self, win):
        win.select_view("plot_confidence_profile")
        assert win._layout_combo.isVisibleTo(win)
        assert not win._station_combo.isVisibleTo(win)
        win.select_view("plot_station_confidence_dashboard")
        assert win._station_combo.isVisibleTo(win)
        assert not win._layout_combo.isVisibleTo(win)
        win._line_combo.setCurrentIndex(1)
        win.select_view("plot_confidence_profile")
        assert not win._layout_combo.isVisibleTo(win)  # one line picked

    def test_views_draw(self, win):
        for fn in ("plot_confidence_profile", "plot_snr_section",
                   "ss_qc_psection", "plot_station_confidence_dashboard"):
            win.select_view(fn)
            win._on_draw()
            assert win._status.text().startswith("✓"), (fn,
                                                         win._status.text())

    def test_invalid_option_shows_reason(self, win):
        win.select_view("plot_confidence_profile")
        fa, _fv = win._current_forms()
        fa.set_values({"ci_lo": 0.99, "ci_hi": 0.5})
        win._on_draw()
        assert "Waiting" in win._status.text()

    def test_static_shift_method_rebuilds_its_parameters(self, win):
        win.select_view("ss_qc_profile")
        fa, _fv = win._current_forms()
        keys_ama = set(fa.widgets)
        fa.set_values({"method": "loess"})
        fa2, _ = win._current_forms()
        assert fa2 is not fa and fa2.values()["method"] == "loess"
        assert set(fa2.widgets) != keys_ama

    def test_filter_hides_diagnostics(self, win):
        win._filter.setText("strike")
        visible = [win._view_list.item(r).text().strip()
                   for r in range(win._view_list.count())
                   if not win._view_list.item(r).isHidden()]
        assert "STRIKE" in visible and "Strike profile" in visible
        assert "Confidence profile" not in visible

    def test_summary_tab(self, win, tmp_path):
        win._tabs.setCurrentIndex(1)
        df = win._summary
        assert df is not None and win._summary_table.rowCount() == len(df)
        assert "pass" in win._summary_counts.text()
        win._status_filter.setCurrentIndex(
            win._status_filter.findData("warn"))
        assert win._summary_table.rowCount() == (df["Status"] == "warn").sum()
        path = win.export_summary(str(tmp_path / "qc.csv"))
        assert Path(path).read_text().startswith("Line,Station,Status")
        win._thr["ci_lo"].setValue(0.999)
        win._thr["ci_hi"].setValue(1.0)
        assert win._summary is not df  # thresholds recompute it

    def test_open_station_from_summary(self, win):
        win.open_station("22-15U")
        assert win._view().fn == "plot_station_confidence_dashboard"
        assert win._station_combo.currentData() == "22-15U"
        assert win._line_combo.currentData() == "L22"

    def test_quicklook_tab(self, win):
        win._tabs.setCurrentIndex(2)
        assert not win._quicklook_stale
        assert win._quick.figure.axes

    def test_figures_stay_white_on_dark_theme(self, win):
        win.set_dark_mode(True)
        assert win._ctrl.dark is False
