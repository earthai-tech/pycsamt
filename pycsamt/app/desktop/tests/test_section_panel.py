# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for SectionPanel (Phase 5)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.panels.section_panel import (
    SectionPanel,
)
from pycsamt.app.desktop.widgets.colorbar_widget import (
    ColorbarWidget,
)
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas


@pytest.fixture
def panel(qapp):
    p = SectionPanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_section_panel_creates(qapp):
    p = SectionPanel()
    assert p is not None
    p.close()


def test_has_load_button(panel):
    assert panel._btn_load is not None


def test_has_compare_button(panel):
    assert panel._btn_compare is not None
    assert not panel._btn_compare.isEnabled()  # disabled until result loaded


def test_has_export_button(panel):
    assert panel._btn_export is not None
    assert not panel._btn_export.isEnabled()


def test_has_canvas(panel):
    assert isinstance(panel._canvas, MplCanvas)


def test_has_colorbar(panel):
    assert isinstance(panel._cbar, ColorbarWidget)


def test_has_summary_text(panel):
    assert panel._summary is not None


# ── Controls ──────────────────────────────────────────────────────────────


def test_rho_min_default(panel):
    assert panel._rho_min.value() == pytest.approx(1.0)


def test_rho_max_default(panel):
    assert panel._rho_max.value() == pytest.approx(1000.0)


def test_depth_max_special_value(panel):
    assert panel._depth_max.value() == 0.0
    assert panel._depth_max.specialValueText() == "Auto"


def test_show_stations_default_checked(panel):
    assert panel._show_stations.isChecked()


def test_cmap_combo_has_options(panel):
    assert panel._cmap_combo.count() > 0


# ── No-data state ─────────────────────────────────────────────────────────


def test_clear_disables_buttons(panel):
    panel.clear()
    assert not panel._btn_compare.isEnabled()
    assert not panel._btn_export.isEnabled()


def test_redraw_without_result_does_not_raise(panel):
    panel._redraw()


def test_dark_mode_toggle_does_not_raise(panel):
    panel.set_dark_mode(False)
    panel.set_dark_mode(True)


# ── Mock result ───────────────────────────────────────────────────────────


class _MockMesh:
    x_nodes = [0, 500, 1000, 1500, 2000]
    z_nodes = [0, 100, 300, 600]
    station_x = [500, 1000, 1500]


class _MockResult:
    workdir = "/tmp/mock_inversion"
    n_iterations = 42
    final_rms = 0.987
    mesh = _MockMesh()
    summary = "MockResult: 42 iterations, RMS=0.987"

    import numpy as np

    rho_2d = np.log10(np.full((3, 4), 100.0))  # 3 z-cells × 4 x-cells


def test_set_result_enables_buttons(panel):
    panel.set_result(_MockResult())
    assert panel._btn_compare.isEnabled()
    assert panel._btn_export.isEnabled()


def test_set_result_stores_result(panel):
    result = _MockResult()
    panel.set_result(result)
    assert panel._result is result


def test_set_result_does_not_raise(panel):
    panel.set_result(_MockResult())  # redraw with mock data


def test_clear_after_result_removes_result(panel):
    panel.set_result(_MockResult())
    panel.clear()
    assert panel._result is None
    assert not panel._btn_compare.isEnabled()
