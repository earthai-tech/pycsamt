# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AppController — central application state and callback bus.

No Qt imports: this class is usable from the Dash web app as-is.
The desktop app registers Qt-connected callables; the web app registers
Dash callback targets.  Both share identical business logic.
"""

from __future__ import annotations

from typing import Callable

from pycsamt.app.desktop.models.session import SessionState


class AppController:
    """
    Central hub holding the current survey state and notifying listeners.

    Attributes
    ----------
    session : SessionState
        Persistent session (theme, recent files, last selection…)
    sites : object or None
        The currently loaded ``Sites`` collection (set by DataController).
    selected_station : str or None
        Station ID currently highlighted across all panels.
    """

    def __init__(self, session: SessionState | None = None) -> None:
        self.session: SessionState = session or SessionState.load()
        self.sites = None
        self.all_sites = None
        self.active_lines: set[str] = set(self.session.active_lines)
        self.primary_line: str | None = self.session.primary_line
        self.selected_station: str | None = self.session.selected_station

        # ── Callback lists (listeners register here) ──────────────────
        self._on_data_loaded: list[Callable] = []
        self._on_station_selected: list[Callable] = []
        self._on_status_message: list[Callable] = []
        self._on_active_lines_changed: list[Callable] = []

    # ── Registration ──────────────────────────────────────────────────

    def on_data_loaded(self, callback: Callable) -> None:
        """Register *callback(sites)* to fire after successful data load."""
        self._on_data_loaded.append(callback)

    def on_station_selected(self, callback: Callable) -> None:
        """Register *callback(station_id: str)* to fire on selection change."""
        self._on_station_selected.append(callback)

    def on_status_message(self, callback: Callable) -> None:
        """Register *callback(message: str)* for status-bar updates."""
        self._on_status_message.append(callback)

    def on_active_lines_changed(self, callback: Callable) -> None:
        self._on_active_lines_changed.append(callback)

    # ── State mutations (called by workers / panels) ───────────────────

    def set_sites(self, sites) -> None:
        """Store the loaded Sites and notify all data-loaded listeners."""
        self.sites = sites
        self.all_sites = sites
        for cb in self._on_data_loaded:
            try:
                cb(sites)
            except Exception:
                pass

    def set_active_scope(
        self, sites, active_lines: set[str], primary_line: str | None
    ) -> None:
        """Update the application-wide filtered survey without reloading data."""
        self.sites = sites
        self.active_lines = set(active_lines)
        self.primary_line = primary_line
        self.session.active_lines = sorted(active_lines)
        self.session.primary_line = primary_line
        for callback in self._on_active_lines_changed:
            try:
                callback(sites, set(active_lines), primary_line)
            except Exception:
                pass

    def select_station(self, station_id: str) -> None:
        """Select a station and notify all listeners."""
        self.selected_station = station_id
        self.session.selected_station = station_id
        for cb in self._on_station_selected:
            try:
                cb(station_id)
            except Exception:
                pass

    def post_status(self, message: str) -> None:
        """Broadcast a status message (e.g. to the main-window status bar)."""
        for cb in self._on_status_message:
            try:
                cb(message)
            except Exception:
                pass

    # ── Convenience ───────────────────────────────────────────────────

    def add_recent_file(self, path: str) -> None:
        recent = self.session.recent_files
        if path in recent:
            recent.remove(path)
        recent.insert(0, path)
        self.session.recent_files = recent[:20]

    @property
    def n_stations(self) -> int:
        if self.sites is None:
            return 0
        try:
            return len(self.sites)
        except TypeError:
            return 0
