# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for QCPanel.

QCController's own draw() dispatch is exhaustively covered in
test_qc_controller.py; here we only need panel-level integration:
category/item wiring, refresh/set_sites/set_dark_mode orchestration, and
the export/empty-state paths.
"""

from __future__ import annotations

from pathlib import Path
from unittest import mock

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.qc_controller import ALL_GROUPS
from pycsamt.app.desktop.panels.qc_panel import QCPanel

_ROOT = Path(__file__).parents[3]
_TIPPER = _ROOT / "data" / "MT" / "kap03lmt_edis"
_HAS_TIPPER = _TIPPER.exists() and any(_TIPPER.glob("*.edi"))


@pytest.fixture(scope="session")
def tipper_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_TIPPER:
        pytest.skip("TIPPER data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_TIPPER))


@pytest.fixture
def panel(qapp):
    p = QCPanel()
    yield p
    p.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_creates(qapp):
    p = QCPanel()
    assert p is not None
    p.close()


def test_categories_populated(panel):
    labels = [group_label for group_label, _ in ALL_GROUPS]
    assert panel._bar.current_category() == 0
    assert panel._bar.current_item() >= 0
    assert len(labels) == len(ALL_GROUPS)


def test_canvas_created(panel):
    assert panel._canvas is not None


def test_initial_draw_is_empty_annotation(panel):
    ax = panel._canvas.figure.axes[0]
    assert len(ax.texts) >= 1


# ── set_sites / set_dark_mode ─────────────────────────────────────────────


def test_set_sites_updates_controller(panel, tipper_sites):
    panel.set_sites(tipper_sites)
    assert panel._ctrl._sites is tipper_sites


def test_set_sites_does_not_raise_with_real_data(panel, tipper_sites):
    panel.set_sites(tipper_sites)  # must not raise (exercises refresh() too)


def test_set_dark_mode_toggles_and_refreshes(panel):
    panel.set_dark_mode(False)
    assert panel._ctrl.dark is False
    panel.set_dark_mode(True)
    assert panel._ctrl.dark is True


# ── refresh() guards ───────────────────────────────────────────────────────


def test_refresh_out_of_range_category_noop(panel, monkeypatch):
    monkeypatch.setattr(panel._bar, "current_category", lambda: 999)
    panel.refresh()  # must not raise


def test_refresh_negative_category_noop(panel, monkeypatch):
    monkeypatch.setattr(panel._bar, "current_category", lambda: -1)
    panel.refresh()  # must not raise


def test_refresh_out_of_range_item_noop(panel, monkeypatch):
    monkeypatch.setattr(panel._bar, "current_item", lambda: 999)
    panel.refresh()  # must not raise


def test_refresh_negative_item_noop(panel, monkeypatch):
    monkeypatch.setattr(panel._bar, "current_item", lambda: -1)
    panel.refresh()  # must not raise


# ── category / item change wiring ─────────────────────────────────────────


@pytest.mark.parametrize("cat_idx", range(len(ALL_GROUPS)))
def test_on_category_changed_populates_items_and_renders(panel, cat_idx, tipper_sites):
    panel.set_sites(tipper_sites)
    panel._on_category_changed(cat_idx)
    _, plot_list = ALL_GROUPS[cat_idx]
    assert panel._bar._combo_item.count() == len(plot_list)
    assert panel._bar.current_item() == 0


def test_on_item_changed_calls_refresh(panel):
    with mock.patch.object(panel, "refresh") as m:
        panel._on_item_changed(1)
        m.assert_called_once()


def test_populate_items_invalid_index_noop(panel):
    panel._populate_items(-1)
    panel._populate_items(999)  # must not raise


# ── export ────────────────────────────────────────────────────────────────


def test_on_export_opens_export_dialog(panel, monkeypatch):
    captured = []

    class _FakeDialog:
        def __init__(self, *a, **kw):
            captured.append(kw)

        def exec(self):
            return None

    monkeypatch.setattr(
        "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog", _FakeDialog
    )
    panel._on_export()
    assert captured
    assert captured[0]["figure"] is panel._canvas.figure


def test_export_clicked_signal_triggers_on_export(panel, monkeypatch):
    with mock.patch.object(panel, "_on_export") as m:
        panel._bar.export_clicked.emit()
        m.assert_called_once()


# ── refresh_clicked wiring ────────────────────────────────────────────────


def test_refresh_clicked_signal_triggers_refresh(panel):
    with mock.patch.object(panel, "refresh") as m:
        panel._bar.refresh_clicked.emit()
        m.assert_called_once()


# ── _draw_empty ───────────────────────────────────────────────────────────


def test_draw_empty_directly(panel):
    panel._draw_empty()
    ax = panel._canvas.figure.axes[0]
    assert any("QC" in t.get_text() or "data" in t.get_text() for t in ax.texts)
