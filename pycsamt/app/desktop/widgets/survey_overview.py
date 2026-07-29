# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SurveyOverviewWidget — read-only aggregate statistics for the loaded dataset.

Shown in the right column of the main window above the StationDetailCard.
Populated by ``update_survey(df, paths)`` after data is loaded; reset to
the placeholder state by ``clear()``.
"""

from __future__ import annotations

import os
from pathlib import Path

import pandas as pd
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)


class SurveyOverviewWidget(QWidget):
    """
    Compact, read-only survey statistics card.

    Call ``update_survey(df, paths)`` to populate it and ``clear()`` to
    return to the empty placeholder state.

    Parameters
    ----------
    parent : QWidget, optional
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("SurveyOverview")
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self._build_ui()
        self.clear()

    # ── Construction ──────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(12, 10, 12, 8)
        root.setSpacing(6)

        # Title row
        title_row = QHBoxLayout()
        title_row.setSpacing(6)
        icon = QLabel("◈")
        icon.setObjectName("SurveyIcon")
        title = QLabel("Survey Overview")
        title.setObjectName("SurveyOverviewTitle")
        title_row.addWidget(icon)
        title_row.addWidget(title)
        title_row.addStretch()
        root.addLayout(title_row)

        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        root.addWidget(sep)

        # Empty-state placeholder
        self._placeholder = QLabel("No survey loaded\nOpen EDI files to begin")
        self._placeholder.setObjectName("SurveyPlaceholder")
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._placeholder.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        root.addWidget(self._placeholder)

        # Stats grid (hidden until data is loaded)
        self._stats_frame = QWidget()
        self._stats_frame.setObjectName("SurveyStatsFrame")
        grid = QGridLayout(self._stats_frame)
        grid.setContentsMargins(0, 2, 0, 2)
        grid.setVerticalSpacing(3)
        grid.setHorizontalSpacing(10)
        grid.setColumnMinimumWidth(0, 100)

        def _row(r: int, key: str) -> QLabel:
            k = QLabel(key)
            k.setObjectName("SurveyStatKey")
            v = QLabel("—")
            v.setObjectName("SurveyStatValue")
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

        self._v_stations = _row(0, "Stations")
        self._v_lat = _row(1, "Latitude")
        self._v_lon = _row(2, "Longitude")
        self._v_elev = _row(3, "Elevation")
        self._v_nfreq = _row(4, "Freq. count")
        self._v_tipper = _row(5, "With tipper")
        self._v_format = _row(6, "Format")
        self._v_source = _row(7, "Source")

        root.addWidget(self._stats_frame)
        root.addStretch(1)

    # ── Public API ────────────────────────────────────────────────

    def update_survey(
        self,
        df: pd.DataFrame,
        paths: list[str] | None = None,
    ) -> None:
        """
        Populate the panel from the station summary DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            Columns: ID, Latitude, Longitude, Elevation, N_freq, Tipper
        paths : list of str, optional
            Original file paths — used to derive format and source dir.
        """
        n = len(df)
        if n == 0:
            self.clear()
            return

        self._placeholder.setVisible(False)
        self._stats_frame.setVisible(True)

        # Stations
        self._v_stations.setText(str(n))

        # Latitude range
        lat = df["Latitude"].dropna()
        if len(lat) >= 2:
            self._v_lat.setText(f"{lat.min():.3f}° – {lat.max():.3f}°")
        elif len(lat) == 1:
            self._v_lat.setText(f"{lat.iloc[0]:.3f}°")
        else:
            self._v_lat.setText("—")

        # Longitude range
        lon = df["Longitude"].dropna()
        if len(lon) >= 2:
            self._v_lon.setText(f"{lon.min():.3f}° – {lon.max():.3f}°")
        elif len(lon) == 1:
            self._v_lon.setText(f"{lon.iloc[0]:.3f}°")
        else:
            self._v_lon.setText("—")

        # Elevation range
        elev = df["Elevation"].dropna()
        if len(elev) >= 2:
            self._v_elev.setText(f"{elev.min():.0f} – {elev.max():.0f} m")
        elif len(elev) == 1:
            self._v_elev.setText(f"{elev.iloc[0]:.0f} m")
        else:
            self._v_elev.setText("—")

        # Frequency count (min/max/median)
        nf = df["N_freq"].dropna()
        if len(nf) > 0 and nf.max() > 0:
            lo, hi, med = int(nf.min()), int(nf.max()), int(nf.median())
            if lo == hi:
                self._v_nfreq.setText(str(lo))
            else:
                self._v_nfreq.setText(f"{lo} – {hi}  (med {med})")
        else:
            self._v_nfreq.setText("—")

        # Tipper count
        n_tip = int(df["Tipper"].sum()) if "Tipper" in df.columns else 0
        self._v_tipper.setText(f"{n_tip} / {n}")

        # Format (inferred from file extensions)
        if paths:
            exts = {Path(p).suffix.lower().lstrip(".").upper() for p in paths}
            self._v_format.setText("  ".join(sorted(exts)) or "—")
        else:
            self._v_format.setText("EDI")

        # Source directory
        if paths:
            dirs = {str(Path(p).parent) for p in paths}
            if len(dirs) == 1:
                src = list(dirs)[0]
                parts = Path(src).parts
                if len(parts) > 4:
                    src = os.path.join("…", *parts[-3:])
                self._v_source.setText(src)
            else:
                self._v_source.setText(f"{len(dirs)} directories")
        else:
            self._v_source.setText("—")

    def clear(self) -> None:
        """Return to the empty placeholder state."""
        self._placeholder.setVisible(True)
        self._stats_frame.setVisible(False)
