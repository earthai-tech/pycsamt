# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Evidence panel of the Interpretation Studio.

Everything interpretation is built from, in one column: the resistivity
model (any PCSF/PCSM file or inversion run folder, with its line for
multi-line models), boreholes (CSV, LAS, PCBH), structural measurements,
the rock database, and the monitoring inputs (repeat surveys, a second
model to fuse).

Model loading can be slow (a ModEM 3-D run is converted first, a MARE2DEM
mesh is resampled), so the panel only *asks* for it
(:attr:`EvidencePanel.model_requested`); the window loads in the
background.  Everything else acts on the controller directly and reports
through :attr:`EvidencePanel.changed`.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidget,
    QMenu,
    QPushButton,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.windows._base import make_group

_MODEL_FILTER = "PCSF / PCSM model (*.pcsf *.pcsm *.pcsm.gz);;All files (*)"


def _small(text: str, slot, tip: str = "") -> QPushButton:
    b = QPushButton(text)
    b.setObjectName("FileListBtn")
    if tip:
        b.setToolTip(tip)
    b.clicked.connect(slot)
    return b


def _info(text: str = "") -> QLabel:
    lbl = QLabel(text)
    lbl.setObjectName("InfoLabel")
    lbl.setWordWrap(True)
    return lbl


class EvidencePanel(QWidget):
    """Model + field evidence for the interpretation.

    Signals
    -------
    model_requested(str, str)
        Path of a model file / run folder and the line ("" = first).
    inversion_requested()
        "From Inversion Studio" was chosen.
    changed(str)
        Evidence changed; the message says what happened.
    """

    model_requested = Signal(str, str)
    inversion_requested = Signal()
    changed = Signal(str)

    def __init__(self, ctrl, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ctrl = ctrl
        self._model_path = ""
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 4, 0)
        v.setSpacing(6)

        # ── model ────────────────────────────────────────────────────────
        grp, lay = make_group("Model")
        self.model_card = QLabel()
        self.model_card.setObjectName("ModelStatusCard")
        self.model_card.setWordWrap(True)
        self.model_card.setTextFormat(Qt.TextFormat.RichText)
        self.model_card.setAlignment(Qt.AlignmentFlag.AlignTop)
        lay.addWidget(self.model_card)
        self.btn_load = QToolButton()
        self.btn_load.setText("Load model  ▾")
        self.btn_load.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.btn_load.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        menu = QMenu(self.btn_load)
        menu.addAction("PCSF / PCSM file…", self._pick_model_file)
        menu.addAction("Inversion run folder…", self._pick_model_folder)
        menu.addSeparator()
        menu.addAction("From Inversion Studio", self.inversion_requested.emit)
        self.btn_load.setMenu(menu)
        lay.addWidget(self.btn_load)
        row = QHBoxLayout()
        self.line_lbl = QLabel("Line:")
        self.line_combo = QComboBox()
        self.line_combo.setToolTip("Profile line of a multi-line / 3-D model")
        self.line_combo.activated.connect(self._on_line)
        row.addWidget(self.line_lbl)
        row.addWidget(self.line_combo, 1)
        lay.addLayout(row)
        self.model_notes = _info()
        lay.addWidget(self.model_notes)
        v.addWidget(grp)

        # ── boreholes ────────────────────────────────────────────────────
        grp, lay = make_group("Boreholes")
        self.bh_list = QListWidget()
        self.bh_list.setObjectName("StackList")
        self.bh_list.setMaximumHeight(90)
        lay.addWidget(self.bh_list)
        row = QHBoxLayout()
        row.addWidget(_small("+ CSV", self.add_borehole_csv,
                             "Borehole log CSV (top, bottom, lithology…)"))
        row.addWidget(_small("+ LAS", self.add_borehole_las, "LAS 2.0 log"))
        row.addWidget(_small("+ PCBH", self.add_borehole_pcbh,
                             "PCBH borehole file (.pcbh.json)"))
        lay.addLayout(row)
        lay.addWidget(_small("Remove selected", self.remove_borehole))
        v.addWidget(grp)

        # ── structure ────────────────────────────────────────────────────
        grp, lay = make_group("Structure")
        self.struct_status = _info("No structural data")
        lay.addWidget(self.struct_status)
        row = QHBoxLayout()
        row.addWidget(_small(
            "+ Planar", lambda: self._add_struct("planar"),
            "CSV: x, kind, strike_deg, dip_deg, dip_direction_deg[, z…]"))
        row.addWidget(_small(
            "+ Linear", lambda: self._add_struct("linear"),
            "CSV: x, kind, trend_deg, plunge_deg[, z…]"))
        lay.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(_small(
            "+ Faults", lambda: self._add_struct("faults"),
            "CSV: x, dip_deg, downthrown_side[, sense, throw_m…]"))
        row.addWidget(_small("Clear", self.clear_structure))
        lay.addLayout(row)
        v.addWidget(grp)

        # ── rock database ────────────────────────────────────────────────
        grp, lay = make_group("Rock database")
        self.db_lbl = _info("Default rock library")
        lay.addWidget(self.db_lbl)
        row = QHBoxLayout()
        row.addWidget(_small("Default", self.use_default_db))
        row.addWidget(_small("Load CSV…", self.load_db_csv,
                             "Resistivity ranges per rock type"))
        lay.addLayout(row)
        v.addWidget(grp)

        # ── monitoring ───────────────────────────────────────────────────
        grp, lay = make_group("Monitoring")
        lay.addWidget(_info("Repeat surveys (same grid as the model):"))
        self.survey_list = QListWidget()
        self.survey_list.setObjectName("StackList")
        self.survey_list.setMaximumHeight(70)
        lay.addWidget(self.survey_list)
        row = QHBoxLayout()
        row.addWidget(_small("+ Survey…", self.add_survey,
                             "A later model of the same survey"))
        row.addWidget(_small("Clear", self.clear_surveys))
        lay.addLayout(row)
        self.second_lbl = _info("Second model (to fuse): none")
        lay.addWidget(self.second_lbl)
        lay.addWidget(_small("Load second model…", self.load_second,
                             "A deeper-reaching model, e.g. MT under AMT"))
        v.addWidget(grp)
        v.addStretch(1)
        self.refresh()

    # ── model ────────────────────────────────────────────────────────────
    def _pick_model_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Open model",
                                              self._start_dir(),
                                              _MODEL_FILTER)
        if path:
            self.model_requested.emit(path, "")

    def _pick_model_folder(self) -> None:
        path = QFileDialog.getExistingDirectory(
            self, "Open an Occam1D / Occam2D / ModEM / MARE2DEM run folder",
            self._start_dir())
        if path:
            self.model_requested.emit(path, "")

    def _on_line(self, _i: int) -> None:
        info = self._ctrl.state.model_info
        line = self.line_combo.currentText()
        if info is not None and line and line != info.line:
            self.model_requested.emit(info.path, line)

    def _start_dir(self) -> str:
        return str(Path(self._model_path).parent) if self._model_path else ""

    # ── boreholes ────────────────────────────────────────────────────────
    def _open(self, title: str, filt: str) -> str:
        path, _ = QFileDialog.getOpenFileName(self, title, self._start_dir(),
                                              filt)
        return path

    def add_borehole_csv(self, path: str = "") -> None:
        path = path or self._open("Borehole CSV",
                                  "CSV files (*.csv);;All files (*)")
        if path:
            self._report(lambda: self._ctrl.add_borehole_csv(path),
                         "Borehole {} added.")

    def add_borehole_las(self, path: str = "") -> None:
        path = path or self._open("Borehole LAS",
                                  "LAS files (*.las *.LAS);;All files (*)")
        if path:
            self._report(lambda: self._ctrl.add_borehole_las(path),
                         "Borehole {} added.")

    def add_borehole_pcbh(self, path: str = "", positions=None) -> None:
        """PCBH collars are map coordinates: boreholes named like a model
        station sit at that station, the others are asked for."""
        path = path or self._open(
            "PCBH borehole file",
            "PCBH (*.pcbh.json *.json);;All files (*)")
        if not path:
            return
        try:
            ids, pos = self._ctrl.pcbh_positions(path)
        except Exception as exc:
            self.changed.emit(f"Could not read {Path(path).name}: {exc}")
            return
        pos.update(positions or {})
        for bid in ids:
            if bid in pos:
                continue
            x, ok = QInputDialog.getDouble(
                self, "Borehole position",
                f"Distance of borehole {bid} along the model profile (m):",
                0.0, -1e7, 1e7, 1)
            if not ok:
                self.changed.emit("PCBH import cancelled.")
                return
            pos[bid] = x
        self._report(lambda: ", ".join(self._ctrl.add_borehole_pcbh(path,
                                                                    pos)),
                     "Borehole(s) {} added.")

    def remove_borehole(self) -> None:
        item = self.bh_list.currentItem()
        if item is not None:
            name = item.text()  # refresh() deletes the item
            self._ctrl.remove_borehole(name)
            self.refresh()
            self.changed.emit(f"Borehole {name} removed.")

    # ── structure ────────────────────────────────────────────────────────
    def _add_struct(self, kind: str, path: str = "") -> None:
        path = path or self._open(f"Structural {kind} CSV",
                                  "CSV files (*.csv);;All files (*)")
        if path:
            fn = getattr(self._ctrl, f"add_structural_{kind}_csv")
            self._report(lambda: fn(path), f"{{}} {kind} record(s) added.")

    def clear_structure(self) -> None:
        self._ctrl.clear_structural_model()
        self.refresh()
        self.changed.emit("Structural data cleared.")

    # ── rock database ────────────────────────────────────────────────────
    def use_default_db(self) -> None:
        self._ctrl.set_rock_db_default()
        self._db_name = ""
        self.refresh()
        self.changed.emit("Using the default rock library.")

    def load_db_csv(self, path: str = "") -> None:
        path = path or self._open("Rock database CSV",
                                  "CSV files (*.csv);;All files (*)")
        if path:
            try:
                self._ctrl.set_rock_db_csv(path)
            except Exception as exc:
                self.changed.emit(f"Rock database not loaded: {exc}")
                return
            self._db_name = Path(path).name
            self.refresh()
            self.changed.emit(f"Rock database {self._db_name} loaded.")

    # ── monitoring ───────────────────────────────────────────────────────
    def _pick_any_model(self, title: str) -> str:
        """A PCSF/PCSM file, or (cancel) an inversion run folder."""
        path, _ = QFileDialog.getOpenFileName(
            self, f"{title} — PCSF/PCSM file (cancel for a run folder)",
            self._start_dir(), _MODEL_FILTER)
        if not path:
            path = QFileDialog.getExistingDirectory(
                self, f"{title} — inversion run folder", self._start_dir())
        return path

    def add_survey(self, path: str = "", label: str = "") -> None:
        if self._ctrl.state.model is None:
            self.changed.emit("Load the baseline model first.")
            return
        path = path or self._pick_any_model("Repeat survey model")
        if not path:
            return
        if not label:
            label, ok = QInputDialog.getText(
                self, "Survey label", "Label (e.g. the survey year):",
                text=Path(path).stem)
            if not ok:
                return
        self._report(lambda: self._ctrl.add_timelapse_survey(path, label),
                     "Repeat survey added ({} surveys).")

    def clear_surveys(self) -> None:
        self._ctrl.clear_timelapse()
        self.refresh()
        self.changed.emit("Repeat surveys cleared.")

    def load_second(self, path: str = "") -> None:
        path = path or self._pick_any_model("Second model")
        if path:
            self._report(lambda: self._ctrl.set_secondary_model(path),
                         "Second model: {}.")

    # ── state ────────────────────────────────────────────────────────────
    def _report(self, fn, ok_msg: str) -> None:
        try:
            out = fn()
        except Exception as exc:
            self.changed.emit(f"Failed: {exc}")
            return
        self.refresh()
        self.changed.emit(ok_msg.format(out))

    def set_model_path(self, path: str) -> None:
        self._model_path = path

    def refresh(self) -> None:
        st = self._ctrl.state
        self._refresh_model_card()
        self.bh_list.clear()
        for bh in st.boreholes:
            self.bh_list.addItem(bh.name)
        sm = st.structural_model
        if sm is None or not (sm.faults or sm.planar or sm.linear):
            self.struct_status.setText("No structural data")
        else:
            self.struct_status.setText(
                f"{len(sm.planar)} planar · {len(sm.linear)} linear · "
                f"{len(sm.faults)} fault(s)")
        name = getattr(self, "_db_name", "")
        self.db_lbl.setText(f"Custom: {name}" if name
                            else "Default rock library")
        self.survey_list.clear()
        labels = st.timelapse_labels or []
        for i, _m in enumerate(st.timelapse_surveys):
            self.survey_list.addItem(labels[i] if i < len(labels)
                                     else f"survey {i + 1}")
        sec = st.secondary_model
        self.second_lbl.setText(
            "Second model (to fuse): " + (
                f"{getattr(sec, 'method', 'model')} · "
                f"{sec.n_x} × {sec.n_z}" if sec is not None else "none"))

    def _refresh_model_card(self) -> None:
        st = self._ctrl.model_status
        info = self._ctrl.state.model_info
        lines = list(info.lines) if info is not None else []
        many = len(lines) > 1
        self.line_lbl.setVisible(many)
        self.line_combo.setVisible(many)
        if many:
            self.line_combo.blockSignals(True)
            self.line_combo.clear()
            self.line_combo.addItems(lines)
            if info.line in lines:
                self.line_combo.setCurrentText(info.line)
            self.line_combo.blockSignals(False)
        self.model_notes.setText("\n".join(info.notes) if info else "")
        self.model_notes.setVisible(bool(info and info.notes))
        if not st["loaded"]:
            self.model_card.setText(
                "<b>No model</b><br><span style='color:#6b7280'>Load a "
                "PCSF/PCSM file or an inversion run folder.</span>")
            return
        src = info.source if info is not None else st["method"]
        depth = st["depth_max"]
        prof = st["profile_m"]
        rms = st["rms"]
        n_st = len(getattr(self._ctrl.state.model, "station_names", None)
                   or [])
        bits = [f"{st['n_x']} × {st['n_z']} cells"]
        if n_st:
            bits.append(f"{n_st} stations")
        if depth:
            bits.append(f"to {depth:,.0f} m")
        if prof:
            bits.append(f"{prof / 1e3:.1f} km")
        if rms is not None:
            try:
                bits.append(f"RMS {float(rms):.2f}")
            except (TypeError, ValueError):
                pass
        self.model_card.setText(
            f"<b>✓ {src}</b><br><span style='color:#6b7280'>"
            f"{' · '.join(bits)}</span>")


__all__ = ["EvidencePanel"]
