# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
StationResponseDialog — per-station impedance tensor + tipper viewer.

Shows the full 4-component (Zxx/Zxy/Zyx/Zyy) apparent-resistivity, phase,
and tipper response for any station in the loaded survey.  Calls
``pycsamt.emtools.inspect.plot_station_response`` which produces a 3-row ×
N-column figure that follows the package colour scheme.

Usage
-----
    from pycsamt.app.desktop.tools.station_response_tool import StationResponseDialog
    dlg = StationResponseDialog(sites, parent=self)
    dlg.exec()
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QGroupBox,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

# ── Background worker ─────────────────────────────────────────────────────────


class _ResponseWorker(QThread):
    done = Signal(object)  # Figure
    error = Signal(str)

    def __init__(
        self,
        sites,
        station: str,
        components: tuple,
        show_tipper: bool,
        show_errors: bool,
    ):
        super().__init__()
        self._sites = sites
        self._station = station
        self._components = components
        self._tipper = show_tipper
        self._errors = show_errors

    def run(self):
        try:
            from pycsamt.emtools.inspect import (
                plot_station_response,
            )

            fig = plot_station_response(
                self._sites,
                station=self._station,
                components=self._components,
                show_tipper=self._tipper,
                show_error_bars=self._errors,
                verbose=0,
            )
            self.done.emit(fig)
        except Exception as exc:
            self.error.emit(str(exc))


# ── Dialog ────────────────────────────────────────────────────────────────────


class StationResponseDialog(QDialog):
    """
    Per-station impedance tensor + tipper viewer.

    Parameters
    ----------
    sites : any
        Loaded survey (path, SitesCollection, list of EDI objects …).
    parent : QWidget, optional
    """

    def __init__(self, sites=None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Station Response Inspector")
        self.setMinimumSize(900, 620)
        self._sites = sites
        self._worker = None
        self._station_names: list[str] = []
        self._build_ui()
        self._populate_stations()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: controls ────────────────────────────────────────────────
        ctrl = QWidget()
        ctrl.setFixedWidth(240)
        ctrl_lay = QVBoxLayout(ctrl)
        ctrl_lay.setContentsMargins(6, 6, 6, 6)
        ctrl_lay.setSpacing(10)

        grp_station = QGroupBox("Station")
        form_s = QFormLayout(grp_station)
        self._station_combo = QComboBox()
        self._station_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        form_s.addRow("Select:", self._station_combo)
        ctrl_lay.addWidget(grp_station)

        grp_comp = QGroupBox("Components")
        form_c = QFormLayout(grp_comp)
        self._cb_xx = QCheckBox("Zxx")
        self._cb_xy = QCheckBox("Zxy")
        self._cb_xy.setChecked(True)
        self._cb_yx = QCheckBox("Zyx")
        self._cb_yx.setChecked(True)
        self._cb_yy = QCheckBox("Zyy")
        for cb in (self._cb_xx, self._cb_xy, self._cb_yx, self._cb_yy):
            form_c.addRow(cb)
        ctrl_lay.addWidget(grp_comp)

        grp_opt = QGroupBox("Options")
        form_o = QFormLayout(grp_opt)
        self._cb_tipper = QCheckBox("Show tipper")
        self._cb_tipper.setChecked(True)
        self._cb_errors = QCheckBox("Show error bars")
        self._cb_errors.setChecked(True)
        form_o.addRow(self._cb_tipper)
        form_o.addRow(self._cb_errors)
        ctrl_lay.addWidget(grp_opt)

        self._run_btn = QPushButton("Plot Response")
        self._run_btn.clicked.connect(self._on_plot)
        ctrl_lay.addWidget(self._run_btn)

        self._status_lbl = QLabel("")
        self._status_lbl.setWordWrap(True)
        ctrl_lay.addWidget(self._status_lbl)
        ctrl_lay.addStretch()

        splitter.addWidget(ctrl)

        # ── Right: canvas ─────────────────────────────────────────────────
        self._canvas = MplCanvas(parent=self)
        splitter.addWidget(self._canvas)
        splitter.setStretchFactor(1, 1)

        root.addWidget(splitter, stretch=1)

        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Populate ──────────────────────────────────────────────────────────────

    def _populate_stations(self) -> None:
        if self._sites is None:
            self._status_lbl.setText("No survey loaded.")
            self._run_btn.setEnabled(False)
            return
        try:
            from pycsamt.emtools._core import (
                _iter_items,
                _unwrap,
            )

            for ed in _iter_items(self._sites):
                try:
                    ed = _unwrap(ed)
                except Exception:
                    pass
                name = (
                    getattr(ed, "station", None)
                    or getattr(ed, "id", None)
                    or "?"
                )
                self._station_names.append(str(name))
        except Exception as exc:
            self._status_lbl.setText(f"Cannot read stations: {exc}")
            self._run_btn.setEnabled(False)
            return

        self._station_combo.addItems(self._station_names)
        if self._station_names:
            self._status_lbl.setText(
                f"{len(self._station_names)} stations loaded."
            )

    # ── Plot ──────────────────────────────────────────────────────────────────

    def _on_plot(self) -> None:
        station = self._station_combo.currentText()
        if not station:
            self._status_lbl.setText("Select a station first.")
            return

        components = tuple(
            c
            for c, cb in (
                ("xx", self._cb_xx),
                ("xy", self._cb_xy),
                ("yx", self._cb_yx),
                ("yy", self._cb_yy),
            )
            if cb.isChecked()
        )
        if not components:
            self._status_lbl.setText("Select at least one component.")
            return

        self._run_btn.setEnabled(False)
        self._status_lbl.setText(f"Computing {station}…")

        self._worker = _ResponseWorker(
            self._sites,
            station,
            components,
            self._cb_tipper.isChecked(),
            self._cb_errors.isChecked(),
        )
        self._worker.done.connect(self._on_done)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_done(self, fig) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText("Done.")
        self._canvas.show_figure(fig)

    def _on_error(self, msg: str) -> None:
        self._run_btn.setEnabled(True)
        self._status_lbl.setText(f"Error: {msg}")
