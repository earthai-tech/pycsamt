# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
RecomputeDialog — modal dialog for recomputing EDI files via EDIRecomputer.

Layout
──────
  ┌──────────────────────────────────────────────────────────────────┐
  │  ⚙  Recompute EDI Files                               [× Close]  │
  ├──────────────────────────────────────────────────────────────────┤
  │  Source                                                           │
  │    ● Use loaded stations   ○ Select folder / files               │
  │    [Folder path ____________________________] [Browse…]           │
  ├──────────────────────────────────────────────────────────────────┤
  │  Processing options                                               │
  │    Rotation angle (°): [___]   Components: [Z + Tip ▾]           │
  │    Frequency range:    Min [_____] Hz   Max [_____] Hz           │
  │    Fill missing:       [none ▾]                                   │
  │    ☑ Recompute apparent-resistivity / phase from impedance        │
  ├──────────────────────────────────────────────────────────────────┤
  │  Output                                                           │
  │    Output folder: [recomputed_edis___________] [Browse…]         │
  │    Filename template: [{source_stem}.edi______]                   │
  │    ☑ Preserve line subdirectories   ☑ Overwrite   ☑ Manifest     │
  │    ☐ Write to disk (uncheck → in-memory only)                    │
  ├──────────────────────────────────────────────────────────────────┤
  │  Progress                                                         │
  │    [██████████████████░░░░░░░░░░░░░░░]   34 / 128                │
  │    Currently:  L18PLT / S034                                     │
  │    ✔ 33 ok     ✕ 1 failed     ○ 94 pending                      │
  │  ┌─────────────────────────────────────────────────────────────┐ │
  │  │ ✔ S001   ✔ S002   ✕ S003 (file exists)  ✔ S004  …         │ │
  │  └─────────────────────────────────────────────────────────────┘ │
  ├──────────────────────────────────────────────────────────────────┤
  │                 [ Close ]        [ ▶  Recompute ]                 │
  └──────────────────────────────────────────────────────────────────┘
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont, QIcon
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QSpacerItem,
    QVBoxLayout,
    QWidget,
)

_ICONS = Path(__file__).resolve().parents[1] / "resources" / "icons"


def _icon(name: str) -> QIcon:
    for candidate in (name, f"{name}.svg", f"{name}.png"):
        path = _ICONS / candidate
        if path.exists():
            return QIcon(str(path))
    return QIcon()


class RecomputeDialog(QDialog):
    """
    Dialog for configuring and running an EDI recomputation workflow.

    Parameters
    ----------
    sites : Sites or None
        Currently loaded Sites collection.  When provided the dialog
        defaults to "Use loaded stations".
    parent : QWidget or None
        Parent widget.

    Signals
    -------
    recompute_committed(object)
        Emitted with the ``EDIRecomputeResult`` once recomputation finishes
        successfully and the user elects to load the result back.
    """

    recompute_committed = Signal(object)  # EDIRecomputeResult

    # ── Construction ──────────────────────────────────────────────────

    def __init__(
        self,
        sites: Any = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._sites = sites
        self._worker = None
        self._result = None
        self._done = 0
        self._total = 0
        self._ok = 0
        self._failed = 0

        self.setWindowTitle("Recompute EDI Files")
        self.setWindowIcon(_icon("recompute"))
        self.setMinimumSize(640, 680)
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Expanding
        )
        self._build_ui()
        self._refresh_source_ui()

    # ── UI ────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(16, 14, 16, 12)
        root.setSpacing(12)

        root.addWidget(self._make_source_group())
        root.addWidget(self._make_options_group())
        root.addWidget(self._make_output_group())
        root.addWidget(self._make_progress_group())

        root.addSpacerItem(
            QSpacerItem(
                0, 4, QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed
            )
        )
        root.addLayout(self._make_buttons())

    # ── Source group ──────────────────────────────────────────────────

    def _make_source_group(self) -> QGroupBox:
        box = QGroupBox("Source")
        v = QVBoxLayout(box)
        v.setSpacing(6)

        self._rb_loaded = QRadioButton("Use loaded stations")
        self._rb_folder = QRadioButton("Select folder / files")
        self._rb_loaded.setChecked(self._sites is not None)
        self._rb_folder.setChecked(self._sites is None)
        self._rb_loaded.setEnabled(self._sites is not None)
        self._rb_loaded.toggled.connect(self._refresh_source_ui)

        v.addWidget(self._rb_loaded)
        v.addWidget(self._rb_folder)

        folder_row = QHBoxLayout()
        self._folder_edit = QLineEdit()
        self._folder_edit.setPlaceholderText("Path to EDI folder or file…")
        self._folder_edit.setObjectName("RecomputeFolderEdit")
        btn_browse = QPushButton("Browse…")
        btn_browse.setObjectName("BrowseButton")
        btn_browse.setFixedWidth(84)
        btn_browse.clicked.connect(self._browse_source)
        folder_row.addWidget(self._folder_edit)
        folder_row.addWidget(btn_browse)
        v.addLayout(folder_row)

        return box

    # ── Options group ─────────────────────────────────────────────────

    def _make_options_group(self) -> QGroupBox:
        box = QGroupBox("Processing Options")
        v = QVBoxLayout(box)
        v.setSpacing(8)

        # Rotation
        rot_row = QHBoxLayout()
        rot_row.addWidget(QLabel("Rotation angle (°):"))
        self._rot_spin = QDoubleSpinBox()
        self._rot_spin.setRange(-360.0, 360.0)
        self._rot_spin.setDecimals(1)
        self._rot_spin.setSpecialValueText("none")
        self._rot_spin.setValue(self._rot_spin.minimum())
        self._rot_spin.setFixedWidth(90)
        rot_row.addWidget(self._rot_spin)
        rot_row.addSpacing(16)
        rot_row.addWidget(QLabel("Components:"))
        self._comp_combo = QComboBox()
        self._comp_combo.addItems(["Z + Tip (both)", "Z only", "Tip only"])
        self._comp_combo.setFixedWidth(140)
        rot_row.addWidget(self._comp_combo)
        rot_row.addStretch()
        v.addLayout(rot_row)

        # Frequency range
        freq_row = QHBoxLayout()
        freq_row.addWidget(QLabel("Frequency range:"))
        freq_row.addWidget(QLabel("Min"))
        self._fmin_edit = QLineEdit()
        self._fmin_edit.setPlaceholderText("Hz")
        self._fmin_edit.setFixedWidth(80)
        freq_row.addWidget(self._fmin_edit)
        freq_row.addWidget(QLabel("Max"))
        self._fmax_edit = QLineEdit()
        self._fmax_edit.setPlaceholderText("Hz")
        self._fmax_edit.setFixedWidth(80)
        freq_row.addWidget(self._fmax_edit)
        freq_row.addStretch()
        v.addLayout(freq_row)

        # Fill missing
        fill_row = QHBoxLayout()
        fill_row.addWidget(QLabel("Fill missing values:"))
        self._fill_combo = QComboBox()
        self._fill_combo.addItems(["none", "nan", "zero"])
        self._fill_combo.setFixedWidth(100)
        fill_row.addWidget(self._fill_combo)
        fill_row.addStretch()
        v.addLayout(fill_row)

        # Rho/phase checkbox
        self._chk_resphase = QCheckBox(
            "Recompute apparent-resistivity / phase from impedance"
        )
        self._chk_resphase.setChecked(True)
        v.addWidget(self._chk_resphase)

        return box

    # ── Output group ──────────────────────────────────────────────────

    def _make_output_group(self) -> QGroupBox:
        box = QGroupBox("Output")
        v = QVBoxLayout(box)
        v.setSpacing(8)

        # Output folder
        out_row = QHBoxLayout()
        out_row.addWidget(QLabel("Output folder:"))
        self._out_edit = QLineEdit()
        self._out_edit.setPlaceholderText(
            "recomputed_edis  (auto next to source)"
        )
        out_row.addWidget(self._out_edit)
        btn_out = QPushButton("Browse…")
        btn_out.setObjectName("BrowseButton")
        btn_out.setFixedWidth(84)
        btn_out.clicked.connect(self._browse_output)
        out_row.addWidget(btn_out)
        v.addLayout(out_row)

        # Filename template
        tpl_row = QHBoxLayout()
        tpl_row.addWidget(QLabel("Filename template:"))
        self._tpl_edit = QLineEdit("{source_stem}.edi")
        tpl_row.addWidget(self._tpl_edit)
        v.addLayout(tpl_row)

        # Checkboxes row
        chk_row = QHBoxLayout()
        chk_row.setSpacing(20)
        self._chk_preserve = QCheckBox("Preserve line subdirs")
        self._chk_preserve.setChecked(True)
        self._chk_overwrite = QCheckBox("Overwrite")
        self._chk_overwrite.setChecked(False)
        self._chk_manifest = QCheckBox("Write manifest CSV")
        self._chk_manifest.setChecked(True)
        self._chk_write = QCheckBox("Write to disk")
        self._chk_write.setChecked(True)
        chk_row.addWidget(self._chk_preserve)
        chk_row.addWidget(self._chk_overwrite)
        chk_row.addWidget(self._chk_manifest)
        chk_row.addWidget(self._chk_write)
        chk_row.addStretch()
        v.addLayout(chk_row)

        return box

    # ── Progress group ────────────────────────────────────────────────

    def _make_progress_group(self) -> QGroupBox:
        box = QGroupBox("Progress")
        v = QVBoxLayout(box)
        v.setSpacing(6)

        # Main progress bar
        bar_row = QHBoxLayout()
        self._pbar = QProgressBar()
        self._pbar.setRange(0, 100)
        self._pbar.setValue(0)
        self._pbar.setObjectName("RecomputeProgressBar")
        self._pbar.setTextVisible(False)
        bar_row.addWidget(self._pbar)
        self._pbar_lbl = QLabel("—")
        self._pbar_lbl.setFixedWidth(80)
        self._pbar_lbl.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        bar_row.addWidget(self._pbar_lbl)
        v.addLayout(bar_row)

        # Current station label
        self._cur_lbl = QLabel("Ready")
        self._cur_lbl.setObjectName("RecomputeCurrentLabel")
        v.addWidget(self._cur_lbl)

        # Stats row
        stats_row = QHBoxLayout()
        self._ok_lbl = QLabel("✔  0 ok")
        self._failed_lbl = QLabel("✕  0 failed")
        self._pending_lbl = QLabel("○  — pending")
        self._ok_lbl.setObjectName("RecomputeOkLabel")
        self._failed_lbl.setObjectName("RecomputeFailedLabel")
        self._pending_lbl.setObjectName("RecomputePendingLabel")
        stats_row.addWidget(self._ok_lbl)
        stats_row.addSpacing(20)
        stats_row.addWidget(self._failed_lbl)
        stats_row.addSpacing(20)
        stats_row.addWidget(self._pending_lbl)
        stats_row.addStretch()
        v.addLayout(stats_row)

        # Per-station result log
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setObjectName("Separator")
        v.addWidget(sep)

        self._log_list = QListWidget()
        self._log_list.setObjectName("RecomputeLogList")
        self._log_list.setMinimumHeight(90)
        self._log_list.setMaximumHeight(140)
        self._log_list.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        self._log_list.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        v.addWidget(self._log_list)

        return box

    # ── Buttons ───────────────────────────────────────────────────────

    def _make_buttons(self) -> QHBoxLayout:
        h = QHBoxLayout()
        h.setSpacing(8)
        h.addStretch()

        self._btn_close = QPushButton("Close")
        self._btn_close.setObjectName("RecomputeCloseBtn")
        self._btn_close.clicked.connect(self._on_close_clicked)

        self._btn_run = QPushButton("Recompute")
        self._btn_run.setIcon(_icon("recompute"))
        self._btn_run.setObjectName("RecomputeRunBtn")
        self._btn_run.setDefault(True)
        self._btn_run.clicked.connect(self._on_run_clicked)

        # "Load Result" button — appears after a successful run
        self._btn_load = QPushButton("Load Recomputed")
        self._btn_load.setIcon(_icon("recompute"))
        self._btn_load.setObjectName("RecomputeLoadBtn")
        self._btn_load.setVisible(False)
        self._btn_load.clicked.connect(self._on_load_result)

        h.addWidget(self._btn_close)
        h.addWidget(self._btn_load)
        h.addWidget(self._btn_run)
        return h

    # ── Source-mode refresh ───────────────────────────────────────────

    def _refresh_source_ui(self) -> None:
        use_folder = self._rb_folder.isChecked()
        self._folder_edit.setEnabled(use_folder)
        # Find the Browse button (sibling in same layout)
        # It is always enabled when folder mode is selected
        for w in self.findChildren(QPushButton, "BrowseButton"):
            if w.text() == "Browse…" and w.parent() is not None:
                # First Browse button is for source
                break

    # ── Browse slots ──────────────────────────────────────────────────

    def _browse_source(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, "Select EDI survey folder", str(Path.home())
        )
        if folder:
            self._folder_edit.setText(folder)

    def _browse_output(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, "Select output folder", str(Path.home())
        )
        if folder:
            self._out_edit.setText(folder)

    # ── Run / cancel ─────────────────────────────────────────────────

    def _on_run_clicked(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.requestInterruption()
            self._worker.wait(2000)
            self._set_running(False)
            self._cur_lbl.setText("Cancelled.")
            return

        source = self._resolve_source()
        if source is None:
            return

        kw = self._collect_kwargs()
        self._reset_progress()
        self._set_running(True)

        from pycsamt.app.desktop.workers.recompute_worker import (
            RecomputeWorker,
        )

        self._worker = RecomputeWorker(source, parent=self, **kw)
        self._worker.station_done.connect(self._on_station_done)
        self._worker.finished.connect(self._on_finished)
        self._worker.error.connect(self._on_error)
        self._worker.start()

    def _on_close_clicked(self) -> None:
        if self._worker is not None and self._worker.isRunning():
            self._worker.requestInterruption()
            self._worker.wait(2000)
        self.reject()

    # ── Source resolution ─────────────────────────────────────────────

    def _resolve_source(self) -> Any:
        if self._rb_loaded.isChecked() and self._sites is not None:
            return self._sites
        folder = self._folder_edit.text().strip()
        if not folder:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(
                self,
                "No source",
                "Please select a folder/file or use loaded stations.",
            )
            return None
        p = Path(folder)
        if not p.exists():
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(
                self, "Path not found", f"Path does not exist:\n{folder}"
            )
            return None
        return p

    # ── Kwargs collector ──────────────────────────────────────────────

    def _collect_kwargs(self) -> dict:
        # Rotation
        rot_val = self._rot_spin.value()
        rotate_angle = (
            None if rot_val == self._rot_spin.minimum() else rot_val
        )
        comp_idx = self._comp_combo.currentIndex()
        rotate_components = {
            0: ("Z", "Tip"),
            1: ("Z",),
            2: ("Tip",),
        }[comp_idx]

        # Frequency
        fmin = self._parse_float(self._fmin_edit.text())
        fmax = self._parse_float(self._fmax_edit.text())

        # Fill
        fill_text = self._fill_combo.currentText()
        fill_missing = None if fill_text == "none" else fill_text

        # Output
        out_text = self._out_edit.text().strip()
        output_root = out_text if out_text else None

        tpl = self._tpl_edit.text().strip() or "{source_stem}.edi"

        return dict(
            output_root=output_root,
            preserve_line_dirs=self._chk_preserve.isChecked(),
            template=tpl,
            overwrite=self._chk_overwrite.isChecked(),
            write=self._chk_write.isChecked(),
            manifest_csv=self._chk_manifest.isChecked(),
            rotate_angle=rotate_angle,
            rotate_components=rotate_components,
            fmin=fmin,
            fmax=fmax,
            fill_missing_values=fill_missing,
            recompute_resphase=self._chk_resphase.isChecked(),
        )

    # ── Progress callbacks ────────────────────────────────────────────

    def _on_station_done(
        self,
        done: int,
        total: int,
        station: str,
        status: str,
        message: str = "",
    ) -> None:
        self._done = done
        self._total = total

        if status == "ok":
            self._ok += 1
            icon = "✔"
            color = "#a6e3a1"  # green
        elif status == "failed":
            self._failed += 1
            icon = "✕"
            color = "#f38ba8"  # red
        else:
            icon = "–"
            color = "#a6adc8"  # muted

        pct = int(done / total * 100) if total > 0 else 0
        self._pbar.setValue(pct)
        self._pbar_lbl.setText(f"{done} / {total}")
        self._cur_lbl.setText(f"Processing:  {station}")

        pending = max(0, total - done)
        self._ok_lbl.setText(f"✔  {self._ok} ok")
        self._failed_lbl.setText(f"✕  {self._failed} failed")
        self._pending_lbl.setText(f"○  {pending} pending")

        # Append station row to log
        item = QListWidgetItem(f"  {icon}  {station}")
        item.setForeground(QColor(color))
        item.setFont(QFont("Courier New", 11))
        self._log_list.addItem(item)

        # Append error detail on a second line when failed
        if status == "failed" and message:
            short = message if len(message) <= 90 else message[:87] + "…"
            err_item = QListWidgetItem(f"      ↳ {short}")
            err_item.setForeground(QColor("#f38ba8"))
            err_item.setFont(QFont("Courier New", 9))
            self._log_list.addItem(err_item)

        self._log_list.scrollToBottom()

    def _on_finished(self, result) -> None:
        self._result = result
        self._set_running(False)
        self._pbar.setValue(100)
        n_ok = len([r for r in result.records if r.status == "ok"])
        n_failed = len(result.failed)
        self._cur_lbl.setText(f"Done — {n_ok} recomputed, {n_failed} failed.")

        # Diagnose common all-fail patterns and surface a helpful hint
        if n_failed > 0 and n_ok == 0:
            failed_msgs = [r.message or "" for r in result.failed]
            if all(
                "FileExistsError" in m or "already exists" in m.lower()
                for m in failed_msgs
                if m
            ):
                out = result.output_root or "the selected output folder"
                hint = QListWidgetItem(
                    f"  ℹ  Output files already exist in {out} — enable Overwrite "
                    "or choose a new output folder."
                )
                hint.setForeground(QColor("#fab387"))  # peach / warning
                hint.setFont(QFont("Courier New", 10))
                self._log_list.addItem(hint)
                self._log_list.scrollToBottom()
            elif failed_msgs:
                first = failed_msgs[0]
                first = first if len(first) <= 120 else first[:117] + "…"
                hint = QListWidgetItem(f"  ℹ  First failure: {first}")
                hint.setForeground(QColor("#fab387"))
                hint.setFont(QFont("Courier New", 10))
                self._log_list.addItem(hint)
                self._log_list.scrollToBottom()

        # Show "Load Recomputed" button so user can pull result back into the app
        self._btn_load.setVisible(n_ok > 0)

    def _on_error(self, message: str) -> None:
        self._set_running(False)
        self._cur_lbl.setText(f"Error: {message}")
        item = QListWidgetItem(f"  ✕  {message}")
        item.setForeground(QColor("#f38ba8"))
        self._log_list.addItem(item)

    # ── Load result back ──────────────────────────────────────────────

    def _on_load_result(self) -> None:
        if self._result is not None:
            self.recompute_committed.emit(self._result)
            self.accept()

    # ── Helpers ───────────────────────────────────────────────────────

    def _reset_progress(self) -> None:
        self._done = 0
        self._total = 0
        self._ok = 0
        self._failed = 0
        self._pbar.setValue(0)
        self._pbar_lbl.setText("—")
        self._cur_lbl.setText("Starting…")
        self._ok_lbl.setText("✔  0 ok")
        self._failed_lbl.setText("✕  0 failed")
        self._pending_lbl.setText("○  — pending")
        self._log_list.clear()
        self._btn_load.setVisible(False)
        self._result = None

    def _set_running(self, running: bool) -> None:
        self._btn_run.setText("Stop" if running else "Recompute")
        # Disable option groups while running
        for gb_name in ("Source", "Processing Options", "Output"):
            for gb in self.findChildren(QGroupBox):
                if gb.title() == gb_name:
                    gb.setEnabled(not running)
                    break

    @staticmethod
    def _parse_float(text: str) -> float | None:
        t = text.strip()
        if not t:
            return None
        try:
            return float(t)
        except ValueError:
            return None
