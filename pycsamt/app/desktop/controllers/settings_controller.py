# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SettingsController — pure-Python controller for all PYCSAMT_* API singletons.

Responsibilities
----------------
* Per-tab apply / reset methods called by each SettingsPage.
* Snapshot and restore all singletons (Cancel support in APIConfigDialog).
* Serialise / deserialise settings to JSON (Save Profile / Load Profile).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


class SettingsController:
    """
    Mediates between the settings UI and the PYCSAMT_* API singletons.

    Each ``apply_<tab>(**kw)`` method writes only the keys present in *kw*,
    so pages can pass a partial dict.  ``snapshot()`` / ``restore()`` use the
    same key-space, making JSON persistence trivial.
    """

    SETTINGS_PATH = Path.home() / ".pycsamt" / "settings.json"

    # ── Snapshot / restore ────────────────────────────────────────────────────

    def snapshot(self) -> dict[str, Any]:
        """Return a JSON-serialisable snapshot of all configurable singletons."""
        snap: dict[str, Any] = {}

        # View controls
        try:
            from pycsamt.api.control import (
                PYCSAMT_CONTROL as C,
            )
            snap["view_controls"] = {
                "rho_view":    C.rho.view,
                "phase_range": list(C.phase.range),
                "phase_unit":  C.phase.unit,
                "phase_wrap":  C.phase.wrap,
                "x_view":      C.x.view,
            }
        except Exception:
            pass

        # Station rendering (pseudosection preset)
        try:
            from pycsamt.api.station import (
                PYCSAMT_STATION_RENDERING as SR,
            )
            ps = SR.pseudosection
            snap["station"] = {
                "side":          ps.side,
                "show_markers":  ps.show_markers,
                "marker_symbol": ps.marker.marker,
                "marker_size":   ps.marker.size,
                "marker_offset": ps.marker.offset,
                "max_labels":    ps.max_labels,
            }
        except Exception:
            pass

        # Section style (pseudosection preset axis)
        try:
            from pycsamt.api.section import (
                PYCSAMT_SECTION as SEC,
            )
            ax = SEC.pseudosection.axis
            snap["section"] = {
                "y_direction":  ax.y_direction,
                "station_side": ax.station_side,
            }
        except Exception:
            pass

        # Topography
        try:
            from pycsamt.topo import PYCSAMT_TOPO as T
            snap["topography"] = {
                "enabled":     T.enabled,
                "exaggeration": T.exaggeration,
                "marker_pad":  T.marker_pad_fraction,
            }
        except Exception:
            pass

        return snap

    def restore(self, snap: dict[str, Any]) -> None:
        """Restore all singletons from a snapshot produced by :meth:`snapshot`."""
        if "view_controls" in snap:
            self.apply_view_controls(**snap["view_controls"])
        if "station" in snap:
            self.apply_station(**snap["station"])
        if "section" in snap:
            self.apply_section(**snap["section"])
        if "topography" in snap:
            self.apply_topography(**snap["topography"])

    # ── Per-tab apply ─────────────────────────────────────────────────────────

    def apply_view_controls(self, **kw: Any) -> None:
        """Write view-control fields to :data:`PYCSAMT_CONTROL`."""
        try:
            from pycsamt.api.control import (
                PYCSAMT_CONTROL as C,
            )
            if "rho_view" in kw:
                C.rho.view = kw["rho_view"]
            if "phase_range" in kw:
                C.phase.range = tuple(kw["phase_range"])
            if "phase_unit" in kw:
                C.phase.unit = kw["phase_unit"]
            if "phase_wrap" in kw:
                C.phase.wrap = bool(kw["phase_wrap"])
            if "x_view" in kw:
                C.x.view = kw["x_view"]
        except Exception:
            pass

    def apply_station(self, **kw: Any) -> None:
        """Write station-marker fields to the pseudosection preset of
        :data:`PYCSAMT_STATION_RENDERING`."""
        try:
            from pycsamt.api.station import (
                PYCSAMT_STATION_RENDERING as SR,
            )
            ps = SR.pseudosection
            if "side" in kw:
                ps.side = kw["side"]
            if "show_markers" in kw:
                ps.show_markers = bool(kw["show_markers"])
            if "marker_symbol" in kw:
                ps.marker.marker = kw["marker_symbol"]
            if "marker_size" in kw:
                ps.marker.size = float(kw["marker_size"])
            if "marker_offset" in kw:
                ps.marker.offset = float(kw["marker_offset"])
            if "max_labels" in kw:
                ps.max_labels = int(kw["max_labels"])
        except Exception:
            pass

    def apply_section(self, **kw: Any) -> None:
        """Write section-axis fields to all pseudosection-family presets of
        :data:`PYCSAMT_SECTION`."""
        try:
            from pycsamt.api.section import (
                PYCSAMT_SECTION as SEC,
            )
            for preset in ("pseudosection", "dashboard", "compact", "publication", "dynamic"):
                try:
                    ax = SEC.style_for(preset).axis
                    if "y_direction" in kw:
                        ax.y_direction = kw["y_direction"]
                    if "station_side" in kw:
                        ax.station_side = kw["station_side"]
                except Exception:
                    pass
        except Exception:
            pass

    def apply_topography(self, **kw: Any) -> None:
        """Write fields to :data:`PYCSAMT_TOPO`."""
        try:
            from pycsamt.topo import PYCSAMT_TOPO as T
            if "enabled" in kw:
                T.enabled = bool(kw["enabled"])
            if "exaggeration" in kw:
                T.exaggeration = float(kw["exaggeration"])
            if "marker_pad" in kw:
                T.marker_pad_fraction = float(kw["marker_pad"])
        except Exception:
            pass

    # ── Reset ─────────────────────────────────────────────────────────────────

    def reset_tab(self, tab: str) -> None:
        """Reset the singleton(s) used by *tab* to package defaults."""
        try:
            if tab == "view_controls":
                from pycsamt.api.control import (
                    PYCSAMT_CONTROL,
                )
                PYCSAMT_CONTROL.reset()
            elif tab == "pseudosections":
                from pycsamt.api.section import (
                    PYCSAMT_SECTION,
                )
                from pycsamt.api.station import (
                    PYCSAMT_STATION_RENDERING,
                )
                PYCSAMT_STATION_RENDERING.reset()
                PYCSAMT_SECTION.reset()
            elif tab == "topography":
                from pycsamt.topo import PYCSAMT_TOPO
                PYCSAMT_TOPO.reset()
            elif tab == "display":
                from pycsamt.api.style import PYCSAMT_STYLE
                PYCSAMT_STYLE.reset()
            elif tab == "interpretation":
                from pycsamt.api.interp import PYCSAMT_INTERP
                PYCSAMT_INTERP.reset()
        except Exception:
            pass

    def reset_all(self) -> None:
        """Reset all API singletons to package defaults."""
        for tab in ("view_controls", "pseudosections", "topography", "display", "interpretation"):
            self.reset_tab(tab)

    # ── Persistence ───────────────────────────────────────────────────────────

    def save(self, path: Path | str | None = None) -> None:
        """Save current settings to *path* (default: ``~/.pycsamt/settings.json``)."""
        path = Path(path or self.SETTINGS_PATH)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.snapshot(), indent=2))

    def load(self, path: Path | str | None = None) -> bool:
        """Load settings from *path*.  Returns True if successful."""
        path = Path(path or self.SETTINGS_PATH)
        if not path.exists():
            return False
        try:
            snap = json.loads(path.read_text())
            self.restore(snap)
            return True
        except Exception:
            return False
