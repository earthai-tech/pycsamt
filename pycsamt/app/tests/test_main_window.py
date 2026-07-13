# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for MainWindow (Phase 1) — runs headless via Qt offscreen platform."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")


@pytest.fixture
def window(qapp, monkeypatch):
    # Isolate from real on-disk session so theme/state defaults are predictable.
    from pycsamt.app.desktop.models.session import (
        SessionState,
    )

    monkeypatch.setattr(SessionState, "load", classmethod(lambda cls: cls()))
    from pycsamt.app.desktop.main_window import MainWindow

    win = MainWindow()
    yield win
    win.close()


# ── Construction ──────────────────────────────────────────────────────────


def test_main_window_creates_without_error(qapp):
    from pycsamt.app.desktop.main_window import MainWindow

    win = MainWindow()
    assert win is not None
    win.close()


def test_window_title(window):
    assert window.windowTitle() == "pycsamt"


def test_minimum_size(window):
    assert window.minimumWidth() >= 1000
    assert window.minimumHeight() >= 600


# ── Dock widgets ──────────────────────────────────────────────────────────
# MainWindow uses a floating-window architecture: map/profile/agent are
# independent QWidget windows, not dock widgets.  Only the Log panel remains
# as a proper QDockWidget.


def test_station_panel_exists(window):
    from pycsamt.app.desktop.panels.station_panel import (
        StationPanel,
    )

    assert isinstance(window._station_panel, StationPanel)


def test_map_window_exists(window):
    from pycsamt.app.desktop.windows.map_window import (
        MapViewerWindow,
    )

    assert isinstance(window._map_win, MapViewerWindow)


def test_profile_window_exists(window):
    from pycsamt.app.desktop.windows.profile_window import (
        ProfileViewerWindow,
    )

    assert isinstance(window._profile_win, ProfileViewerWindow)


def test_agent_master_launcher_exists(window):
    # the embedded agent window was replaced by the Agent Master
    # web app, launched on demand through a bridge handler
    assert callable(window._open_agent_master)


def test_log_dock_exists(window):
    from PySide6.QtWidgets import QDockWidget

    assert isinstance(window._log_dock, QDockWidget)


def test_log_panel_is_attached(window):
    from pycsamt.app.desktop.panels.log_panel import LogPanel

    assert isinstance(window._log_panel, LogPanel)


# ── All panel windows exist ───────────────────────────────────────────────


def test_all_panel_windows_exist(window):
    """All scientific panel windows must be created at startup."""
    for attr in (
        "_map_win",
        "_profile_win",
        "_qc_win",
        "_correction_win",
        "_advanced_win",
        "_tdem_win",
        "_pipeline_win",
        "_forward_win",
        "_inversion_win",
        "_interp_win",
    ):
        assert hasattr(window, attr), f"MainWindow missing {attr}"
        assert getattr(window, attr) is not None


def test_map_window_is_widget(window):
    from PySide6.QtWidgets import QWidget

    assert isinstance(window._map_win, QWidget)


def test_profile_window_is_widget(window):
    from PySide6.QtWidgets import QWidget

    assert isinstance(window._profile_win, QWidget)


def test_map_window_has_parent(window):
    assert window._map_win.parent() is window


def test_profile_window_has_parent(window):
    assert window._profile_win.parent() is window


# ── Central widget (station list + overview) ─────────────────────────────


def test_central_widget_is_not_splitter(window):
    """Central widget is a QWidget container, not a bare QSplitter."""
    from PySide6.QtWidgets import QSplitter

    assert not isinstance(window.centralWidget(), QSplitter)


def test_central_widget_is_widget(window):
    """Central widget must exist and be a QWidget."""
    from PySide6.QtWidgets import QWidget

    cw = window.centralWidget()
    assert isinstance(cw, QWidget)


# ── Status bar ────────────────────────────────────────────────────────────


def test_status_bar_exists(window):
    assert window.statusBar() is not None


def test_set_file_label(window):
    window.set_file_label("EDI: /data/WILLY")
    assert window._status_file_lbl.text() == "EDI: /data/WILLY"


def test_set_freq_label(window):
    window.set_freq_label("Freq: 3.5 Hz")
    assert window._status_freq_lbl.text() == "Freq: 3.5 Hz"


def test_set_status_message(window):
    window.set_status("Loading…", timeout_ms=0)
    # showMessage is fire-and-forget; just confirm no exception


# ── Logging ───────────────────────────────────────────────────────────────


def test_log_delegates_to_log_panel(window):
    window.log("test message")
    text = window._log_panel._text.toPlainText()
    assert "test message" in text


def test_log_contains_ready_message(window):
    text = window._log_panel._text.toPlainText()
    assert "ready" in text.lower()


# ── Theme ─────────────────────────────────────────────────────────────────


def test_default_theme_is_dark(window):
    assert window._session.theme == "dark"


def test_apply_light_theme(window):
    window._apply_theme("light")
    assert window._session.theme == "light"


def test_toggle_theme_switches_dark_to_light(window):
    window._apply_theme("dark")
    window._toggle_theme()
    assert window._session.theme == "light"


def test_toggle_theme_switches_light_to_dark(window):
    window._apply_theme("light")
    window._toggle_theme()
    assert window._session.theme == "dark"


# ── Session persistence ───────────────────────────────────────────────────


def test_save_layout_populates_geometry(window):
    window._save_layout()
    assert window._session.dock_geometry is not None
    assert len(window._session.dock_geometry) > 0


def test_save_layout_populates_state(window):
    window._save_layout()
    assert window._session.dock_state is not None
    assert len(window._session.dock_state) > 0


def test_close_saves_session(window, tmp_path, monkeypatch):
    target = tmp_path / "session.json"
    monkeypatch.setattr(
        "pycsamt.app.desktop.models.session._SESSION_PATH", target
    )
    window.close()
    # Session file should have been written
    assert target.exists()


# ── Public API ────────────────────────────────────────────────────────────


def test_public_api_methods_exist(window):
    assert callable(window.log)
    assert callable(window.set_status)
    assert callable(window.set_file_label)
    assert callable(window.set_freq_label)
