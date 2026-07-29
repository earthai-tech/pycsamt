# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InversionWizardDialog — step-by-step inversion launcher (Phase 5).

Four pages driven by a QStackedWidget + Prev / Next navigation:

  Page 0 — Data      : select stations / profile, set output workdir
  Page 1 — Config    : OccamConfig parameters (modes, mesh, inversion)
  Page 2 — Run       : launch button, live log, progress, stop
  Page 3 — Results   : summary + "Open in Viewer" button

The dialog builds Occam2D input files via ``InputBuilder``, runs the
inversion via ``InversionWorker``, and returns an ``InversionResult``
via ``get_result()``.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Slot
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)


class InversionWizardDialog(QDialog):
    """
    Inversion wizard (Occam2D).

    Parameters
    ----------
    sites : Sites or None
        Currently loaded survey data.
    session : SessionState
        Application session (read: solver path, workdir; write: nothing).
    parent : QWidget, optional
    """

    def __init__(self, sites=None, session=None, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Inversion Wizard — Occam2D")
        self.setMinimumSize(640, 520)

        self._sites = sites
        self._session = session
        self._result = None
        self._worker = None

        self._build_ui()

    # ── Result accessor ───────────────────────────────────────────────

    def get_result(self):
        """Return the ``InversionResult`` after a successful run, or None."""
        return self._result

    # ── Build UI ──────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)

        # Page header
        self._header_lbl = QLabel()
        self._header_lbl.setObjectName("WizardHeader")
        root.addWidget(self._header_lbl)

        # Stacked pages
        self._stack = QStackedWidget()
        root.addWidget(self._stack)

        self._page_data = self._build_page_data()
        self._page_config = self._build_page_config()
        self._page_run = self._build_page_run()
        self._page_result = self._build_page_result()
        for page in (
            self._page_data,
            self._page_config,
            self._page_run,
            self._page_result,
        ):
            self._stack.addWidget(page)

        # Navigation buttons
        nav = QHBoxLayout()
        self._btn_prev = QPushButton("← Back")
        self._btn_prev.setEnabled(False)
        self._btn_prev.clicked.connect(self._go_prev)
        nav.addWidget(self._btn_prev)
        nav.addStretch()

        self._btn_next = QPushButton("Next →")
        self._btn_next.clicked.connect(self._go_next)
        nav.addWidget(self._btn_next)

        self._btn_close = QPushButton("Close")
        self._btn_close.clicked.connect(self.accept)
        self._btn_close.setVisible(False)
        nav.addWidget(self._btn_close)

        root.addLayout(nav)

        self._update_header()

    # ── Page 0: Data ──────────────────────────────────────────────────

    def _build_page_data(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)

        layout.addWidget(
            QLabel("Select the stations / profile to invert and the output directory.")
        )

        # Station list
        layout.addWidget(QLabel("Available stations:"))
        self._station_list = QListWidget()
        self._station_list.setSelectionMode(QListWidget.SelectionMode.ExtendedSelection)
        if self._sites is not None:
            for site in self._sites.as_list():
                self._station_list.addItem(QListWidgetItem(site.name))
        else:
            item = QListWidgetItem("(no stations loaded)")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self._station_list.addItem(item)
        layout.addWidget(self._station_list)

        # Modes
        modes_box = QGroupBox("Impedance modes")
        modes_row = QHBoxLayout(modes_box)
        self._mode_te = QCheckBox("TE (xy)")
        self._mode_te.setChecked(True)
        self._mode_tm = QCheckBox("TM (yx)")
        self._mode_tm.setChecked(True)
        modes_row.addWidget(self._mode_te)
        modes_row.addWidget(self._mode_tm)
        layout.addWidget(modes_box)

        # Output working directory
        wd_default = ""
        if self._session and self._session.inversion_workdir:
            wd_default = self._session.inversion_workdir
        wd_box = QGroupBox("Output working directory")
        wd_row = QHBoxLayout(wd_box)
        self._workdir_edit = QLineEdit(wd_default)
        self._workdir_edit.setPlaceholderText("Select output directory…")
        btn_wd = QPushButton("Browse…")
        btn_wd.clicked.connect(self._browse_workdir)
        wd_row.addWidget(self._workdir_edit)
        wd_row.addWidget(btn_wd)
        layout.addWidget(wd_box)

        return w

    def _browse_workdir(self) -> None:
        d = QFileDialog.getExistingDirectory(
            self,
            "Select output directory",
            self._workdir_edit.text() or str(Path.home()),
        )
        if d:
            self._workdir_edit.setText(d)

    # ── Page 1: Config ────────────────────────────────────────────────

    def _build_page_config(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)
        layout.addWidget(QLabel("Configure Occam2D inversion parameters."))

        # Mesh
        mesh_box = QGroupBox("Mesh")
        mesh_form = QFormLayout(mesh_box)

        self._n_layers = QSpinBox()
        self._n_layers.setRange(5, 100)
        self._n_layers.setValue(30)
        mesh_form.addRow("Layers:", self._n_layers)

        self._n_airlayers = QSpinBox()
        self._n_airlayers.setRange(1, 20)
        self._n_airlayers.setValue(5)
        mesh_form.addRow("Air layers:", self._n_airlayers)

        self._cell_size_h = QDoubleSpinBox()
        self._cell_size_h.setRange(10, 10_000)
        self._cell_size_h.setValue(100.0)
        self._cell_size_h.setSuffix(" m")
        mesh_form.addRow("Cell size (horizontal):", self._cell_size_h)

        self._depth_scale = QDoubleSpinBox()
        self._depth_scale.setRange(1.01, 2.0)
        self._depth_scale.setSingleStep(0.01)
        self._depth_scale.setValue(1.2)
        mesh_form.addRow("Depth scale factor:", self._depth_scale)

        layout.addWidget(mesh_box)

        # Inversion
        inv_box = QGroupBox("Inversion")
        inv_form = QFormLayout(inv_box)

        self._max_iter = QSpinBox()
        self._max_iter.setRange(1, 500)
        self._max_iter.setValue(100)
        inv_form.addRow("Max iterations:", self._max_iter)

        self._target_rms = QDoubleSpinBox()
        self._target_rms.setRange(0.1, 50.0)
        self._target_rms.setValue(1.0)
        inv_form.addRow("Target RMS:", self._target_rms)

        self._initial_rho = QDoubleSpinBox()
        self._initial_rho.setRange(1.0, 100_000.0)
        self._initial_rho.setValue(100.0)
        self._initial_rho.setSuffix(" Ω·m")
        inv_form.addRow("Initial resistivity:", self._initial_rho)

        self._error_floor_rho = QDoubleSpinBox()
        self._error_floor_rho.setRange(0.001, 1.0)
        self._error_floor_rho.setValue(0.05)
        inv_form.addRow("Error floor (ρ):", self._error_floor_rho)

        self._error_floor_phase = QDoubleSpinBox()
        self._error_floor_phase.setRange(0.001, 10.0)
        self._error_floor_phase.setValue(0.5)
        self._error_floor_phase.setSuffix(" °")
        inv_form.addRow("Error floor (φ):", self._error_floor_phase)

        layout.addWidget(inv_box)
        return w

    # ── Page 2: Run ───────────────────────────────────────────────────

    def _build_page_run(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)

        self._run_log = QPlainTextEdit()
        self._run_log.setReadOnly(True)
        self._run_log.setMaximumBlockCount(2000)
        layout.addWidget(self._run_log)

        self._run_progress = QProgressBar()
        self._run_progress.setRange(0, 100)
        self._run_progress.setVisible(False)
        layout.addWidget(self._run_progress)

        btn_row = QHBoxLayout()
        self._btn_launch = QPushButton("▶  Launch Occam2D")
        self._btn_launch.clicked.connect(self._on_launch)
        btn_row.addWidget(self._btn_launch)

        self._btn_stop_inv = QPushButton("■  Stop")
        self._btn_stop_inv.setEnabled(False)
        self._btn_stop_inv.clicked.connect(self._on_stop_inversion)
        btn_row.addWidget(self._btn_stop_inv)
        layout.addLayout(btn_row)

        return w

    # ── Page 3: Results ───────────────────────────────────────────────

    def _build_page_result(self) -> QWidget:
        w = QWidget()
        layout = QVBoxLayout(w)

        layout.addWidget(QLabel("Inversion complete."))

        self._result_text = QPlainTextEdit()
        self._result_text.setReadOnly(True)
        layout.addWidget(self._result_text)

        self._btn_open_viewer = QPushButton("Open in 2D Section Viewer")
        self._btn_open_viewer.setEnabled(False)
        layout.addWidget(self._btn_open_viewer)

        return w

    # ── Navigation ────────────────────────────────────────────────────

    _HEADERS = [
        "Step 1 / 4 — Data Selection",
        "Step 2 / 4 — Inversion Configuration",
        "Step 3 / 4 — Run Inversion",
        "Step 4 / 4 — Results",
    ]

    def _update_header(self) -> None:
        idx = self._stack.currentIndex()
        self._header_lbl.setText(f"<b>{self._HEADERS[idx]}</b>")
        self._btn_prev.setEnabled(idx > 0)
        is_last = idx == self._stack.count() - 1
        self._btn_next.setVisible(not is_last)
        self._btn_close.setVisible(is_last)

    def _go_prev(self) -> None:
        i = self._stack.currentIndex()
        if i > 0:
            self._stack.setCurrentIndex(i - 1)
            self._update_header()

    def _go_next(self) -> None:
        i = self._stack.currentIndex()
        if i < self._stack.count() - 1:
            if not self._validate_page(i):
                return
            self._stack.setCurrentIndex(i + 1)
            self._update_header()

    def _validate_page(self, page_idx: int) -> bool:
        if page_idx == 0:
            if not self._workdir_edit.text().strip():
                self._workdir_edit.setFocus()
                return False
        return True

    # ── Launch / Stop ─────────────────────────────────────────────────

    def _on_launch(self) -> None:
        from pycsamt.app.desktop.workers.inversion_worker import (
            InversionWorker,
        )
        from pycsamt.models.occam2d import (
            InputBuilder,
            OccamConfig,
        )

        workdir = self._workdir_edit.text().strip()
        if not workdir:
            self._run_log.appendPlainText("ERROR: no working directory specified.")
            return

        Path(workdir).mkdir(parents=True, exist_ok=True)

        # Build config from form values
        modes = []
        if self._mode_te.isChecked():
            modes.append("te")
        if self._mode_tm.isChecked():
            modes.append("tm")

        config = OccamConfig(
            modes=modes,
            n_layers=self._n_layers.value(),
            n_airlayers=self._n_airlayers.value(),
            cell_size_horizontal=self._cell_size_h.value(),
            depth_scale=self._depth_scale.value(),
            max_iterations=self._max_iter.value(),
            target_misfit=self._target_rms.value(),
            initial_rho=self._initial_rho.value(),
            error_floor_rho=self._error_floor_rho.value(),
            error_floor_phase=self._error_floor_phase.value(),
        )

        # Build input files
        selected = self._selected_sites()
        try:
            self._run_log.appendPlainText("Building Occam2D input files…")
            builder = InputBuilder(source=selected, workdir=workdir, config=config)
            builder.build()
            self._run_log.appendPlainText("Input files written.")
        except Exception as exc:
            self._run_log.appendPlainText(f"InputBuilder error: {exc}")
            return

        # Start worker
        binary = (
            self._session.occam2d_binary
            if self._session and self._session.occam2d_binary
            else None
        )
        self._worker = InversionWorker(
            workdir=workdir,
            binary_path=binary,
            max_iter=self._max_iter.value(),
            target_misfit=self._target_rms.value(),
            parent=self,
        )
        self._worker.stdout_line.connect(self._run_log.appendPlainText)
        self._worker.progress.connect(self._run_progress.setValue)
        self._worker.finished.connect(self._on_inversion_finished)
        self._worker.error.connect(self._on_inversion_error)

        self._run_progress.setValue(0)
        self._run_progress.setVisible(True)
        self._btn_launch.setEnabled(False)
        self._btn_stop_inv.setEnabled(True)
        self._worker.start()

    def _on_stop_inversion(self) -> None:
        if self._worker:
            self._worker.cancel()
            self._worker.quit()
        self._btn_launch.setEnabled(True)
        self._btn_stop_inv.setEnabled(False)
        self._run_progress.setVisible(False)
        self._run_log.appendPlainText("Stopped.")

    @Slot(object)
    def _on_inversion_finished(self, result) -> None:
        self._result = result
        self._run_progress.setVisible(False)
        self._btn_launch.setEnabled(True)
        self._btn_stop_inv.setEnabled(False)
        self._btn_next.click()  # advance to Results page

        try:
            self._result_text.setPlainText(result.summary)
        except Exception:
            self._result_text.setPlainText("Inversion completed successfully.")
        self._btn_open_viewer.setEnabled(True)

    @Slot(str)
    def _on_inversion_error(self, message: str) -> None:
        self._run_log.appendPlainText(f"ERROR: {message}")
        self._run_progress.setVisible(False)
        self._btn_launch.setEnabled(True)
        self._btn_stop_inv.setEnabled(False)

    # ── Helpers ───────────────────────────────────────────────────────

    def _selected_sites(self):
        """Return a Sites subset from the list selection, or all sites."""
        if self._sites is None:
            return None
        selected_names = {item.text() for item in self._station_list.selectedItems()}
        if not selected_names:
            return self._sites
        try:
            return self._sites.select(list(selected_names))
        except Exception:
            return self._sites

    def _build_occam_config(self):
        """Return an ``OccamConfig`` from current form values."""
        from pycsamt.models.occam2d import OccamConfig

        modes = []
        if self._mode_te.isChecked():
            modes.append("te")
        if self._mode_tm.isChecked():
            modes.append("tm")
        return OccamConfig(
            modes=modes,
            n_layers=self._n_layers.value(),
            n_airlayers=self._n_airlayers.value(),
            cell_size_horizontal=self._cell_size_h.value(),
            depth_scale=self._depth_scale.value(),
            max_iterations=self._max_iter.value(),
            target_misfit=self._target_rms.value(),
            initial_rho=self._initial_rho.value(),
            error_floor_rho=self._error_floor_rho.value(),
            error_floor_phase=self._error_floor_phase.value(),
        )
