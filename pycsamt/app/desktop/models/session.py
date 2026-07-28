# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SessionState — JSON-serialisable application session.

Phase 1: Persists last-opened files, dock layout bytes, theme name,
frequency range, and selected station across restarts.  Written to
~/.pycsamt/session.json via save() / load() class methods.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

_SESSION_PATH = Path.home() / ".pycsamt" / "session.json"


@dataclass
class SessionState:
    """Persistent session state for the pycsamt desktop application."""

    theme: str = "dark"
    recent_files: list[str] = field(default_factory=list)
    last_data_dir: str = ""
    selected_station: str | None = None
    freq_min_hz: float | None = None
    freq_max_hz: float | None = None
    overlay: str = "Apparent Resistivity"
    dock_geometry: str | None = None  # base64-encoded QMainWindow geometry
    dock_state: str | None = None  # base64-encoded QMainWindow state (docks/toolbars)

    # ── Per-window geometries (profile / map / qc / agent) ────────────
    # Each value is {"geometry": "<base64>", "visible": bool}
    window_geometries: dict = field(default_factory=dict)

    # ── Phase 5: solver paths + LLM key + advanced prefs ──────────────
    api_key: str = ""  # Anthropic / OpenAI API key
    occam2d_binary: str = ""  # path to Occam2D executable
    modem_binary: str = ""  # path to ModEM executable
    mare2dem_binary: str = ""  # path to MARE2DEM executable
    inversion_workdir: str = ""  # default working directory for inversions
    log_level: str = "WARNING"  # Python logging level name
    tile_provider: str = "OpenStreetMap.Mapnik"  # contextily tile source
    max_recent_files: int = 20  # cap on Recent Files list

    # ------------------------------------------------------------------
    # Persistence
    # ------------------------------------------------------------------

    def save(self, path: Path | None = None) -> None:
        if path is None:
            path = _SESSION_PATH
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(asdict(self), fh, indent=2)

    @classmethod
    def load(cls, path: Path = _SESSION_PATH) -> SessionState:
        if not path.exists():
            return cls()
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
            return cls(
                **{k: v for k, v in data.items() if k in cls.__dataclass_fields__}
            )
        except Exception:
            return cls()
