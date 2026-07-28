# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationDetailCard — shows per-station metadata when a row is selected.

Layout (outer):
    ┌──────────────────────────────────────────┐
    │  ScrollArea (expands)                    │
    │  ┌────────────────────────────────────┐  │
    │  │  ● 18-023A                         │  │
    │  │  ─────────────────────────────     │  │
    │  │  Latitude    32.140 °N             │  │
    │  │  Longitude  119.129 °E             │  │
    │  │  Elevation      69  m              │  │
    │  │  Frequencies    53                 │  │
    │  │  Period range   0.001–1024 s       │  │
    │  │  Tipper        No                  │  │
    │  │  Components  Zxx  Zxy  Zyx  Zyy   │  │
    │  │  ─────────────────────────────     │  │
    │  │  Quality   ● ● ● ○ ○              │  │
    │  └────────────────────────────────────┘  │
    ├──────────────────────────────────────────┤
    │  [Open Profile]        [Show on Map]     │  ← always visible
    └──────────────────────────────────────────┘
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

_ICONS = Path(__file__).parent.parent / "resources" / "icons"


def _lbl(text: str, obj_name: str = "") -> QLabel:
    lbl = QLabel(text)
    if obj_name:
        lbl.setObjectName(obj_name)
    return lbl


class StationDetailCard(QWidget):
    """
    Shows metadata for the currently-selected station.

    The metadata block is inside a QScrollArea so it scrolls when the
    panel is too small to display all rows.  The action buttons are
    pinned below the scroll area so they are always reachable.

    Signals
    -------
    open_profile_requested(station_id)
        Emitted when the user clicks [Open Profile].
    show_on_map_requested(station_id)
        Emitted when the user clicks [Show on Map].
    """

    open_profile_requested = Signal(str)
    show_on_map_requested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._station_id: str = ""
        self._build_ui()
        self.clear()

    # ── Construction ──────────────────────────────────────────────

    def _build_ui(self) -> None:
        self.setObjectName("StationDetailCard")
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # ── Scroll area (all metadata) ─────────────────────────
        scroll = QScrollArea(self)
        scroll.setObjectName("DetailScrollArea")
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setFrameShape(QFrame.Shape.NoFrame)

        body = QWidget()
        body.setObjectName("DetailBody")
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(12, 12, 12, 8)
        body_layout.setSpacing(6)

        # Station name headline
        self._lbl_name = QLabel("No station selected")
        self._lbl_name.setObjectName("StationName")
        self._lbl_name.setWordWrap(True)
        body_layout.addWidget(self._lbl_name)

        sep1 = QFrame()
        sep1.setFrameShape(QFrame.Shape.HLine)
        sep1.setObjectName("Separator")
        body_layout.addWidget(sep1)

        # Key / value grid
        grid = QGridLayout()
        grid.setSpacing(4)
        grid.setColumnMinimumWidth(0, 100)

        def row(r: int, key: str, obj: str = "DetailValue") -> QLabel:
            k = _lbl(key, "DetailKey")
            v = _lbl("—", obj)
            v.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            v.setWordWrap(True)
            grid.addWidget(
                k,
                r,
                0,
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
            )
            grid.addWidget(
                v,
                r,
                1,
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
            )
            return v

        self._lbl_lat = row(0, "Latitude")
        self._lbl_lon = row(1, "Longitude")
        self._lbl_elev = row(2, "Elevation")
        self._lbl_nfreq = row(3, "Frequencies")
        self._lbl_frange = row(4, "Period range")
        self._lbl_tipper = row(5, "Tipper")
        self._lbl_comps = row(6, "Components")

        body_layout.addLayout(grid)

        sep2 = QFrame()
        sep2.setFrameShape(QFrame.Shape.HLine)
        sep2.setObjectName("Separator")
        body_layout.addWidget(sep2)

        # Quality dots
        q_row = QHBoxLayout()
        q_row.setSpacing(6)
        q_row.addWidget(_lbl("Quality:", "DetailKey"))
        self._lbl_quality = QLabel("○ ○ ○ ○ ○")
        self._lbl_quality.setObjectName("QualityDots")
        q_row.addWidget(self._lbl_quality)
        q_row.addStretch()
        body_layout.addLayout(q_row)

        body_layout.addStretch(1)
        scroll.setWidget(body)
        outer.addWidget(scroll, stretch=1)

        # ── Pinned footer: action buttons ──────────────────────
        footer_sep = QFrame()
        footer_sep.setFrameShape(QFrame.Shape.HLine)
        footer_sep.setObjectName("Separator")
        outer.addWidget(footer_sep)

        btn_row = QHBoxLayout()
        btn_row.setContentsMargins(12, 6, 12, 8)
        btn_row.setSpacing(6)
        self._btn_profile = QPushButton("Open Profile")
        self._btn_profile.setObjectName("CardButton")
        self._btn_profile.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._btn_profile.clicked.connect(
            lambda: self.open_profile_requested.emit(self._station_id)
        )
        self._btn_map = QPushButton("Show on Map")
        self._btn_map.setObjectName("CardButton")
        self._btn_map.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._btn_map.clicked.connect(
            lambda: self.show_on_map_requested.emit(self._station_id)
        )
        btn_row.addWidget(self._btn_profile)
        btn_row.addWidget(self._btn_map)
        outer.addLayout(btn_row)

    # ── Public API ────────────────────────────────────────────────

    def update_station(self, station_id: str, sites) -> None:
        """Populate the card from a Sites collection for *station_id*."""
        self._station_id = station_id
        try:
            site = sites.get(station_id)
            s = site.summary()
        except Exception:
            self._lbl_name.setText(station_id)
            return

        import numpy as np

        name = s.get("name", station_id)
        lat = s.get("lat", None)
        lon = s.get("lon", None)
        elev = s.get("elev", None)
        nfreq = s.get("nfreq", 0)
        has_tip = bool(s.get("tipper", False))
        comps = s.get("components", [])

        self._lbl_name.setText(name)

        self._lbl_lat.setText(f"{lat:+.4f} °" if lat is not None else "—")
        self._lbl_lon.setText(f"{lon:+.4f} °" if lon is not None else "—")
        self._lbl_elev.setText(f"{elev:.0f} m" if elev is not None else "—")
        self._lbl_nfreq.setText(str(nfreq))

        # Period range
        try:
            freq = site.freq
            if freq is not None and len(freq) > 0:
                f_sorted = np.sort(freq)
                t_min = 1.0 / f_sorted[-1]
                t_max = 1.0 / f_sorted[0]
                self._lbl_frange.setText(f"{t_min:.2e} – {t_max:.2e} s")
            else:
                self._lbl_frange.setText("—")
        except Exception:
            self._lbl_frange.setText("—")

        self._lbl_tipper.setText("Yes" if has_tip else "No")
        self._lbl_comps.setText("  ".join(comps) if comps else "—")

        # Quality dots (rough heuristic based on nfreq coverage)
        quality = min(5, max(1, nfreq // 10))
        dots = "● " * quality + "○ " * (5 - quality)
        self._lbl_quality.setText(dots.strip())

        self._btn_profile.setEnabled(True)
        self._btn_map.setEnabled(True)

    def clear(self) -> None:
        """Reset to empty / no-selection state."""
        self._station_id = ""
        self._lbl_name.setText("No station selected")
        for lbl in (
            self._lbl_lat,
            self._lbl_lon,
            self._lbl_elev,
            self._lbl_nfreq,
            self._lbl_frange,
            self._lbl_tipper,
            self._lbl_comps,
        ):
            lbl.setText("—")
        self._lbl_quality.setText("○ ○ ○ ○ ○")
        self._btn_profile.setEnabled(False)
        self._btn_map.setEnabled(False)
