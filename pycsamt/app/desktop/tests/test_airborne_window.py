# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Tests for AirborneWindow (pycsamt.app.desktop.windows.airborne_window).

Strategy
--------
* Category/plot combo wiring and the dynamic parameter form are exercised
  against the real ``CATALOGUE`` with no data loaded — ``AirborneController
  .generate()`` degrades to a placeholder figure rather than raising, so
  this covers navigation safely and cheaply.
* ``QFileDialog`` is monkeypatched for Load-button tests, matching the
  pattern used throughout the other window test files (offscreen Qt
  cannot show real native dialogs).
* One real end-to-end test loads real bundled ZTEM EMTF-XML data through
  the actual Load button handler and renders a real diagnostic.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

from PySide6.QtWidgets import QFileDialog

from pycsamt.app.desktop.controllers.airborne_controller import CATALOGUE, CATEGORIES
from pycsamt.app.desktop.windows.airborne_window import AirborneWindow

_ROOT = Path(__file__).resolve().parents[4]
_ZTEM_DIR = _ROOT / "data" / "ZTEM" / "gold_springs_nv"
_HAS_ZTEM = _ZTEM_DIR.exists() and any(_ZTEM_DIR.glob("*.xml"))


@pytest.fixture
def win(qapp):
    w = AirborneWindow(parent=None)
    w.show()  # isVisible() assertions need the ancestor chain actually shown
    yield w
    w.close()


# ── Construction ────────────────────────────────────────────────────────────


class TestConstruction:
    def test_window_title(self, win):
        assert "Airborne" in win.windowTitle()

    def test_category_combo_populated(self, win):
        assert win._combo_category.count() == len(CATEGORIES)

    def test_default_category_is_ztem(self, win):
        assert win._combo_category.currentText() == "ZTEM"

    def test_plot_combo_populated_for_default_category(self, win):
        assert win._combo_plot.count() == len(CATALOGUE["ZTEM"])

    def test_initial_status_label(self, win):
        assert "no airborne data" in win._data_status.text().lower()

    def test_mobilemt_note_mentions_generic_adapter(self, win):
        idx = CATEGORIES.index("MobileMT")
        win._combo_category.setCurrentIndex(idx)
        assert "generic-adapter" in win._category_note.text().lower()
        assert "raw vendor mobilemt" in win._category_note.text().lower()


# ── Category / plot navigation ────────────────────────────────────────────────


class TestNavigation:
    def test_switching_category_repopulates_plot_combo(self, win):
        idx = CATEGORIES.index("AFMAG")
        win._combo_category.setCurrentIndex(idx)
        assert win._combo_plot.count() == len(CATALOGUE["AFMAG"])

    def test_switching_plot_rebuilds_param_form(self, win):
        idx = CATEGORIES.index("AFMAG")
        win._combo_category.setCurrentIndex(idx)
        motion_idx = win._combo_plot.findText("Motion susceptibility map")
        assert motion_idx >= 0
        win._combo_plot.setCurrentIndex(motion_idx)
        assert "inclination" in win._param_widgets
        assert "declination" in win._param_widgets
        assert "roll_amplitude_deg" in win._param_widgets
        assert "pitch_amplitude_deg" in win._param_widgets

    def test_zero_param_entry_shows_no_params_label(self, win):
        idx = win._combo_plot.findText("Flight lines map")
        assert idx >= 0
        win._combo_plot.setCurrentIndex(idx)
        assert win._no_params_lbl.isVisible()
        assert win._param_widgets == {}

    def test_every_category_and_plot_is_reachable(self, win):
        for cat_idx, cat in enumerate(CATEGORIES):
            win._combo_category.setCurrentIndex(cat_idx)
            assert win._combo_plot.count() == len(CATALOGUE[cat])
            for plot_idx in range(win._combo_plot.count()):
                win._combo_plot.setCurrentIndex(plot_idx)  # must not raise


# ── Dynamic parameter form widget kinds ───────────────────────────────────────


class TestParamWidgets:
    def test_get_param_values_reads_combo_spin_dspin_check(self, win):
        idx = win._combo_plot.findText("Divergence profile")
        assert idx >= 0
        win._combo_plot.setCurrentIndex(idx)
        vals = win._get_param_values()
        assert vals["component"] == "tzx"
        assert vals["part"] == "real"
        assert vals["spacing_m"] == pytest.approx(200.0)

        win._param_widgets["spacing_m"].setValue(250.0)
        win._param_widgets["part"].setCurrentText("imag")
        vals = win._get_param_values()
        assert vals["spacing_m"] == pytest.approx(250.0)
        assert vals["part"] == "imag"


# ── Draw (no data loaded) ─────────────────────────────────────────────────────


class TestDrawNoData:
    def test_draw_without_data_shows_placeholder(self, win):
        win._on_draw()
        texts = " ".join(
            t.get_text() for t in win._canvas.figure.axes[0].texts
        )
        assert "no airborne data" in texts.lower()


# ── Load / Clear ───────────────────────────────────────────────────────────────


class TestLoadClear:
    def test_load_cancelled_directory_and_file_dialogs(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog, "getExistingDirectory",
            staticmethod(lambda *a, **k: ""),
        )
        monkeypatch.setattr(
            QFileDialog, "getOpenFileName",
            staticmethod(lambda *a, **k: ("", "")),
        )
        win._on_load()
        assert "no airborne data" in win._data_status.text().lower()

    def test_load_success_updates_status(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog, "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/dir"),
        )
        monkeypatch.setattr(win._ctrl, "load", lambda path: 42)
        win._on_load()
        assert "42" in win._data_status.text()

    def test_load_failure_reported(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog, "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/dir"),
        )

        def _boom(path):
            raise ValueError("bad xml")

        monkeypatch.setattr(win._ctrl, "load", _boom)
        win._on_load()
        assert "Load failed" in win._data_status.text()

    def test_clear_resets_status(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog, "getExistingDirectory",
            staticmethod(lambda *a, **k: "/fake/dir"),
        )
        monkeypatch.setattr(win._ctrl, "load", lambda path: 5)
        win._on_load()
        win._on_clear()
        assert "no airborne data" in win._data_status.text().lower()


# ── Real end-to-end ────────────────────────────────────────────────────────────


@pytest.mark.skipif(not _HAS_ZTEM, reason="ZTEM sample data not available")
class TestRealEndToEnd:
    def test_load_real_ztem_dir_and_draw(self, win, monkeypatch):
        monkeypatch.setattr(
            QFileDialog, "getExistingDirectory",
            staticmethod(lambda *a, **k: str(_ZTEM_DIR)),
        )
        win._on_load()
        assert "station" in win._data_status.text().lower()
        assert win._ctrl.has_data

        idx = win._combo_plot.findText("Flight lines map")
        win._combo_plot.setCurrentIndex(idx)
        win._on_draw()
        ax = win._canvas.figure.axes[0]
        assert len(win._canvas.figure.axes) == 1
        # A real render has a real title (station/line count) and at
        # least one navigation trace per detected flight line, not an
        # error/placeholder message.
        assert "stations" in ax.get_title().lower()
        assert len(ax.lines) > 0


# ── Theme ────────────────────────────────────────────────────────────────────


class TestTheme:
    def test_set_dark_mode_forwards_to_controller(self, win):
        win.set_dark_mode(False)
        assert win._ctrl.dark is False
        win.set_dark_mode(True)
        assert win._ctrl.dark is True
