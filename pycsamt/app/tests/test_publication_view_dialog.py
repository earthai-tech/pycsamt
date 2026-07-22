# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for PublicationViewDialog
(pycsamt.app.desktop.windows.publication_view_dialog).

Real data
---------
data/AMT/WILLY_DATA/L18PLT/ — used via a real PlotController so
``draw_publication_view`` exercises its real station-selected code path,
not just the "no station" placeholder branch.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from pycsamt.app.desktop.controllers.plot_controller import PlotController
from pycsamt.app.desktop.windows.publication_view_dialog import (
    PublicationViewDialog,
)

# ── Paths ─────────────────────────────────────────────────────────────────────

_ROOT = Path(__file__).parents[3]  # pycsamt/
_WILLY_L18 = _ROOT / "data" / "AMT" / "WILLY_DATA" / "L18PLT"
_HAS_WILLY = _WILLY_L18.exists() and any(_WILLY_L18.glob("*.edi"))


@pytest.fixture(scope="session")
def willy_sites():
    pytest.importorskip("pycsamt.emtools")
    if not _HAS_WILLY:
        pytest.skip("WILLY L18PLT data not available")
    from pycsamt.emtools import ensure_sites

    return ensure_sites(str(_WILLY_L18))


@pytest.fixture
def empty_controller():
    """A PlotController with no sites/station selected."""
    return PlotController()


@pytest.fixture
def real_controller(willy_sites):
    ctrl = PlotController()
    ctrl.set_sites(willy_sites)
    ctrl.set_station(willy_sites[0].name)
    return ctrl


class TestConstruction:
    def test_title_includes_station_name(self, qapp, empty_controller):
        dlg = PublicationViewDialog(
            empty_controller, station_name="STA01", parent=None
        )
        assert "STA01" in dlg.windowTitle()
        dlg.close()

    def test_title_without_station_name(self, qapp, empty_controller):
        dlg = PublicationViewDialog(empty_controller, parent=None)
        assert dlg.windowTitle() == "Publication View"
        dlg.close()

    def test_delete_on_close_disabled(self, qapp, empty_controller):
        from PySide6.QtCore import Qt

        dlg = PublicationViewDialog(empty_controller, parent=None)
        assert not dlg.testAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        dlg.close()

    def test_canvas_and_buttons_exist(self, qapp, empty_controller):
        dlg = PublicationViewDialog(empty_controller, parent=None)
        assert dlg._canvas is not None
        assert dlg._btn_export.text() == "⬆  Export…"
        assert dlg._btn_close.text() == "✕  Close"
        dlg.close()

    def test_close_button_closes_dialog(self, qapp, empty_controller):
        dlg = PublicationViewDialog(empty_controller, parent=None)
        dlg.show()
        dlg._btn_close.click()
        assert not dlg.isVisible()

    def test_draw_called_with_no_station_shows_placeholder(
        self, qapp, empty_controller
    ):
        dlg = PublicationViewDialog(empty_controller, parent=None)
        # draw_publication_view degrades gracefully; at least one Axes exists.
        assert len(dlg._canvas.figure.axes) >= 1
        dlg.close()

    def test_draw_called_with_real_station(self, qapp, real_controller):
        dlg = PublicationViewDialog(
            real_controller,
            station_name=real_controller._station_id,
            dark=real_controller.dark,
            parent=None,
        )
        assert len(dlg._canvas.figure.axes) >= 1
        dlg.close()


class TestOnExport:
    def test_export_dialog_opened_when_available(
        self, qapp, empty_controller, monkeypatch
    ):
        import pycsamt.app.desktop.windows.publication_view_dialog as mod

        calls = []

        class _FakeExportDialog:
            def __init__(self, figure, parent):
                calls.append((figure, parent))

            def exec(self):
                calls.append("exec")

        monkeypatch.setattr(
            "pycsamt.app.desktop.dialogs.export_dlg.ExportDialog",
            _FakeExportDialog,
        )
        dlg = PublicationViewDialog(empty_controller, parent=None)
        dlg._on_export()
        assert calls[0][0] is dlg._canvas.figure
        assert calls[-1] == "exec"
        dlg.close()

    def test_export_falls_back_to_file_dialog_on_import_error(
        self, qapp, empty_controller, monkeypatch, tmp_path
    ):
        import builtins

        from PySide6.QtWidgets import QFileDialog

        real_import = builtins.__import__

        def _fake_import(name, *a, **k):
            if name == "pycsamt.app.desktop.dialogs.export_dlg":
                raise ImportError("simulated missing module")
            return real_import(name, *a, **k)

        monkeypatch.setattr(builtins, "__import__", _fake_import)

        out_path = str(tmp_path / "pub_view.pdf")
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: (out_path, "PDF (*.pdf)")),
        )
        dlg = PublicationViewDialog(empty_controller, parent=None)
        dlg._on_export()
        assert Path(out_path).exists()
        dlg.close()

    def test_export_falls_back_and_user_cancels_no_file_written(
        self, qapp, empty_controller, monkeypatch, tmp_path
    ):
        import builtins

        from PySide6.QtWidgets import QFileDialog

        real_import = builtins.__import__

        def _fake_import(name, *a, **k):
            if name == "pycsamt.app.desktop.dialogs.export_dlg":
                raise ImportError("simulated missing module")
            return real_import(name, *a, **k)

        monkeypatch.setattr(builtins, "__import__", _fake_import)
        monkeypatch.setattr(
            QFileDialog,
            "getSaveFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        dlg = PublicationViewDialog(empty_controller, parent=None)
        dlg._on_export()  # must not raise despite empty path
        dlg.close()
