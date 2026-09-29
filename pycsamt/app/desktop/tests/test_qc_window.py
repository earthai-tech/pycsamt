# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for QCDashboardWindow (pycsamt.app.desktop.windows.qc_window).

Strategy mirrors test_tdem_window.py / test_advanced_window.py:
* QCController.draw() is monkeypatched directly for run/export tests so
  no real QC data or matplotlib rendering pipeline is required.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.qc_controller import (
    ALL_GROUPS,
    QCUnavailableResult,
)
from pycsamt.app.desktop.windows.qc_window import QCDashboardWindow


@pytest.fixture
def win(qapp):
    w = QCDashboardWindow(parent=None)
    w.show()
    yield w
    w.close()


def _select_category(win, label):
    for row, (grp_label, _plots) in enumerate(ALL_GROUPS):
        if grp_label == label:
            win._combo_category.setCurrentIndex(row)
            return row
    raise AssertionError(f"category {label!r} not found")


def _select_plot(win, label):
    for row in range(win._combo_plot.count()):
        if win._combo_plot.itemText(row) == label:
            win._combo_plot.setCurrentIndex(row)
            return row
    raise AssertionError(f"plot {label!r} not found")


def test_static_shift_methods_rebuild_controls_and_remember_values(win):
    from PySide6.QtWidgets import QComboBox

    _select_category(win, "Static Shift")
    _select_plot(win, "SS QC profile")
    method = win._parameter_widgets["method"][0]
    assert isinstance(method, QComboBox)
    assert [method.itemData(i) for i in range(method.count())] == [
        "ama", "loess", "bilateral", "refmedian"
    ]
    win._parameter_widgets["half_window"][0].setValue(5)
    method.setCurrentIndex(method.findData("loess"))
    assert not method.isVisible()  # old form disappears before deferred delete
    assert "poly" in win._parameter_widgets
    assert "weights" not in win._parameter_widgets
    assert win._parameter_widgets["half_window"][0].value() == 3
    win._parameter_widgets["it"][0].setValue(4)
    method = win._parameter_widgets["method"][0]
    method.setCurrentIndex(method.findData("bilateral"))
    assert "sig_dist" in win._parameter_widgets
    assert "poly" not in win._parameter_widgets
    win._parameter_widgets["sig_val"][0].setText("-1")
    assert win._collect_plot_kwargs()[1] is not None
    method = win._parameter_widgets["method"][0]
    method.setCurrentIndex(method.findData("refmedian"))
    assert "smooth_sites" in win._parameter_widgets
    assert "half_window" not in win._parameter_widgets
    assert win._collect_plot_kwargs()[1] is None
    method = win._parameter_widgets["method"][0]
    method.setCurrentIndex(method.findData("ama"))
    assert win._parameter_widgets["half_window"][0].value() == 5
    kwargs, error = win._collect_plot_kwargs()
    assert error is None
    assert "poly" not in kwargs and "sig_val" not in kwargs
    method = win._parameter_widgets["method"][0]
    method.setCurrentIndex(method.findData("loess"))
    assert win._parameter_widgets["it"][0].value() == 4


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "QC Dashboard" in win.windowTitle()

    def test_category_combo_populated(self, win):
        assert win._combo_category.count() == len(ALL_GROUPS)

    def test_plot_combo_populated_for_first_category(self, win):
        _label, plots = ALL_GROUPS[0]
        assert win._combo_plot.count() == len(plots)

    def test_desc_label_populated_on_init(self, win):
        assert win._desc_lbl.text() != ""

    def test_status_label_empty_by_default(self, win):
        assert win._status_lbl.text() == ""

    def test_initial_content_is_explanatory_view(self, win):
        assert win._result_stack.currentWidget() is win._unavailable_view
        assert "Load survey data" in win._unavailable_view._title.text()

    def test_category_icons_include_null_and_non_null(self, win, monkeypatch):
        # Force _icon() to always return a null QIcon so the "no icon"
        # branch of _populate_category_combo is exercised too.
        from PySide6.QtGui import QIcon

        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.qc_window._icon",
            lambda name: QIcon(),
        )
        win._combo_category.clear()
        win._populate_category_combo()
        assert win._combo_category.count() == len(ALL_GROUPS)

    def test_canvas_fit_keeps_title_and_labels_inside_view(self, win, qapp):
        win._result_stack.setCurrentWidget(win._plot_page)
        win.resize(900, 600)
        qapp.processEvents()
        figure = win._canvas.figure
        figure.clear()
        axes = figure.add_subplot(111)
        axes.set_title("QC field-zone classification")
        axes.set_xlabel("Station")
        axes.set_ylabel("Period (s)")
        win._canvas.fit_to_view()
        win._canvas._canvas.draw()
        renderer = win._canvas._canvas.get_renderer()
        figure_box = figure.bbox
        assert axes.title.get_window_extent(renderer).y1 <= figure_box.y1 + 1
        assert axes.xaxis.label.get_window_extent(renderer).y0 >= -1
        assert axes.yaxis.label.get_window_extent(renderer).x0 >= -1


# ── Category / plot switching ────────────────────────────────────────────────


class TestCategorySwitching:
    def test_switch_category_populates_plot_combo(self, win):
        row = _select_category(win, "Coverage")
        _label, plots = ALL_GROUPS[row]
        assert win._combo_plot.count() == len(plots)

    def test_on_category_changed_out_of_range_noop(self, win):
        win._on_category_changed(-1)
        win._on_category_changed(999)

    def test_switch_category_resets_plot_index(self, win):
        _select_category(win, "Noise / SNR")
        assert win._combo_plot.currentIndex() == 0

    def test_switch_category_updates_description(self, win):
        _select_category(win, "Skew / Dim")
        assert win._desc_lbl.text() != ""

    def test_plot_changed_updates_description(self, win):
        row = _select_category(win, "Static Shift")
        _label, plots = ALL_GROUPS[row]
        if len(plots) > 1:
            win._combo_plot.setCurrentIndex(1)
        assert win._desc_lbl.text() != ""

    def test_update_desc_bad_index_swallowed(self, win):
        win._update_desc(999, 999)
        assert win._desc_lbl.text() == ""

    def test_field_zone_builds_dynamic_controls(self, win):
        _select_category(win, "Distortion")
        _select_plot(win, "Field zones")
        assert "source_offset" in win._parameter_widgets
        assert "near_threshold" in win._parameter_widgets
        assert "far_threshold" in win._parameter_widgets
        source_widget, kind = win._parameter_widgets["source_offset"]
        assert kind == "optional_float"
        assert win._controls_group.isAncestorOf(source_widget)

    def test_source_overprint_builds_source_offset_control(self, win):
        _select_category(win, "Distortion")
        _select_plot(win, "Near-surface overprint section")
        source_widget, kind = win._parameter_widgets["source_offset"]
        assert kind == "optional_float"
        assert win._controls_group.isAncestorOf(source_widget)

    def test_plot_view_parameters_use_standard_group(self, win):
        _select_category(win, "Distortion")
        _select_plot(win, "Field zones")
        contour_widget, _kind = win._parameter_widgets["contour_kr"]
        assert win._view_group.isAncestorOf(contour_widget)
        assert win._view_group.objectName() == "ParamsGroup"
        assert win._view_group.isVisible()

    def test_actions_group_is_last_parameter_section(self, win):
        layout = win._params_layout
        assert layout.indexOf(win._actions_group) > layout.indexOf(win._view_group)

    def test_plot_has_hover_refresh_overlay(self, win):
        assert win._plot_page.isAncestorOf(win._btn_overlay_refresh)
        assert win._btn_overlay_refresh.objectName() == "PlotRefreshOverlay"
        assert "Hard refresh" in win._btn_overlay_refresh.toolTip()

    def test_field_zone_source_offset_is_collected(self, win):
        _select_category(win, "Distortion")
        _select_plot(win, "Field zones")
        source_widget, _kind = win._parameter_widgets["source_offset"]
        source_widget.setText("1250")
        kwargs, error = win._collect_plot_kwargs()
        assert error is None
        assert kwargs["source_offset"] == pytest.approx(1250.0)

    def test_selection_requests_automatic_render(self, win, monkeypatch):
        requests = []
        monkeypatch.setattr(
            win, "_request_render", lambda delay_ms=0: requests.append(delay_ms)
        )
        _select_category(win, "Coverage")
        _select_plot(win, "Coverage pseudosection")
        assert requests

    def test_parameter_edit_requests_debounced_render(self, win, monkeypatch):
        _select_category(win, "Distortion")
        _select_plot(win, "Field zones")
        requests = []
        monkeypatch.setattr(
            win, "_request_render", lambda delay_ms=0: requests.append(delay_ms)
        )
        source_widget, _kind = win._parameter_widgets["source_offset"]
        source_widget.setText("900")
        assert requests == [350]


# ── Run / Export ──────────────────────────────────────────────────────────────


class TestRunExport:
    def test_run_no_category_or_plot_noop(self, win):
        win._combo_category.setCurrentIndex(-1)
        win._on_run()  # guarded early return, must not raise

    def test_run_without_sites_shows_message(self, win):
        _select_category(win, "Overview")
        assert win._ctrl._sites is None
        win._on_run()
        assert win._status_lbl.text() == "Load survey data first."

    def test_run_success_draws_in_place(self, win, monkeypatch):
        _select_category(win, "Overview")
        win._ctrl._sites = object()
        monkeypatch.setattr(
            win._ctrl, "draw", lambda fn_name, has_ax, fig, **kw: None
        )
        win._on_run()
        assert win._status_lbl.text() == "Done."
        assert win._btn_run.isEnabled()
        assert win._btn_overlay_refresh.isEnabled()
        assert win._result_stack.currentWidget() is win._plot_page

    def test_run_success_replaces_figure(self, win, monkeypatch):
        import matplotlib.figure

        _select_category(win, "Overview")
        win._ctrl._sites = object()
        new_fig = matplotlib.figure.Figure()
        monkeypatch.setattr(
            win._ctrl, "draw", lambda fn_name, has_ax, fig, **kw: new_fig
        )
        win._on_run()
        assert win._status_lbl.text() == "Done."

    def test_unavailable_result_replaces_canvas(self, win, monkeypatch):
        _select_category(win, "Overview")
        win._ctrl._sites = object()

        def _unavailable(fn_name, has_ax, fig, **kw):
            win._ctrl.last_unavailable = QCUnavailableResult(
                "Cannot compute this result",
                "Required input is missing.",
                "Load compatible data.",
            )
            return None

        monkeypatch.setattr(win._ctrl, "draw", _unavailable)
        win._on_run()
        assert win._result_stack.currentWidget() is win._unavailable_view
        assert "Cannot compute" in win._unavailable_view._title.text()
        assert win._btn_export.isEnabled() is False
        assert win._status_lbl.text() == "Result unavailable."

    def test_run_exception_reported(self, win, monkeypatch):
        _select_category(win, "Overview")
        win._ctrl._sites = object()

        def _boom(fn_name, has_ax, fig, **kw):
            raise RuntimeError("draw boom")

        monkeypatch.setattr(win._ctrl, "draw", _boom)
        win._on_run()
        assert "Error: draw boom" in win._status_lbl.text()
        assert win._btn_run.isEnabled()

    def test_run_plot_row_out_of_range_noop(self, win):
        row = _select_category(win, "Overview")
        _label, plots = ALL_GROUPS[row]
        win._ctrl._sites = object()
        win._combo_plot.blockSignals(True)
        win._combo_plot.addItem("extra")
        win._combo_plot.setCurrentIndex(len(plots))  # stale, out of range
        win._combo_plot.blockSignals(False)
        win._on_run()  # must not raise despite plot_row >= len(plots)

    def test_run_status_shows_running_message(self, win, monkeypatch):
        _select_category(win, "Overview")
        win._ctrl._sites = object()
        seen = []

        def _record(fn_name, has_ax, fig, **kw):
            seen.append(fn_name)
            return None

        monkeypatch.setattr(win._ctrl, "draw", _record)
        _label, plots = ALL_GROUPS[0]
        _plot_label, fn_name, _has_ax = plots[0]
        win._on_run()
        assert seen == [fn_name]

    def test_on_export_opens_dialog(self, win, monkeypatch):
        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append(figure)

            def exec(self):
                calls.append("exec")

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        win._on_export()
        assert calls[-1] == "exec"


# ── Auto-render ───────────────────────────────────────────────────────────────


class TestAutoRender:
    def test_auto_render_noop_when_already_rendered(self, win):
        win._auto_rendered = True
        win._ctrl._sites = object()
        win._auto_render_if_ready()  # no crash, still marked rendered

    def test_auto_render_noop_without_sites(self, win):
        win._auto_rendered = False
        win._ctrl._sites = None
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_auto_render_noop_when_not_visible(self, win, monkeypatch):
        win._auto_rendered = False
        win._ctrl._sites = object()
        monkeypatch.setattr(win, "isVisible", lambda: False)
        win._auto_render_if_ready()
        assert win._auto_rendered is False

    def test_auto_render_triggers_timer(self, win, monkeypatch):
        win._auto_rendered = False
        win._ctrl._sites = object()
        calls = []
        monkeypatch.setattr(
            "pycsamt.app.desktop.windows.qc_window.QTimer.singleShot",
            staticmethod(lambda ms, fn: calls.append(fn)),
        )
        win._auto_render_if_ready()
        assert win._auto_rendered is True
        assert len(calls) == 1

    def test_show_event_triggers_auto_render(self, win, monkeypatch):
        from PySide6.QtGui import QShowEvent

        win._ctrl._sites = object()
        win._auto_rendered = False
        calls = []
        monkeypatch.setattr(
            win, "_auto_render_if_ready", lambda: calls.append(1)
        )
        win.showEvent(QShowEvent())
        assert calls == [1]

    def test_set_sites_resets_and_triggers_auto_render(self, win, monkeypatch):
        win._auto_rendered = True
        calls = []
        monkeypatch.setattr(
            win, "_auto_render_if_ready", lambda: calls.append(1)
        )
        win.set_sites(object())
        assert win._auto_rendered is False
        assert calls == [1]


# ── Public API ────────────────────────────────────────────────────────────────


class TestPublicApi:
    def test_set_dark_mode_delegates(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False
        win.set_dark_mode(True)
        assert win._ctrl.dark is True

    def test_set_sites_stores_on_base_and_controller(self, win):
        sites = object()
        win.set_sites(sites)
        assert win._sites is sites
        assert win._ctrl._sites is sites
