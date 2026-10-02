# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
SolverBuilderWindow — compile Occam2D, ModEM 2-D/3-D and MARE2DEM.

┌──────────────┬──────────────────────────────────────────────────────────┐
│ SOLVERS      │  ModEM 3-D                                    ● Built    │
│ ▣ Occam2D    │  3-D MT NLCG inversion (Mod3DMT).                        │
│ ▣ ModEM 2-D  │  Source   (•) Bundled source   ( ) My own folder […]     │
│ ▣ ModEM 3-D  │  Dependencies  ✓ gfortran  ✓ make  ✓ BLAS  ✓ source      │
│ ▣ MARE2DEM   │  Setup   (•) Install missing automatically ( ) Myself    │
│              │  [Install dependencies]              [ Build ModEM 3-D ] │
│              │  Compiling ▓▓▓▓▓▓▓░░░ 68 %   elapsed 00:42      [Stop]   │
│              │  ✓ Built: …\\Mod3DMT.exe   [Copy path] [Open folder]      │
│              │  ▸ Show details (log)                                    │
└──────────────┴──────────────────────────────────────────────────────────┘

The engine is :mod:`pycsamt.models.solver_build` (root-free toolchains via
micromamba, Python-driven ``make`` for the Fortran solvers, a WSL2 build for
MARE2DEM on Windows).  A successful build is registered so the Inversion
window's *Binary* field picks it up automatically.
"""

from __future__ import annotations

import time
from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QFont, QGuiApplication, QIcon
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QFileDialog,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPlainTextEdit,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.compact_button import compact_button
from pycsamt.models import solver_build as sb

_ICONS = Path(__file__).parent.parent / "resources" / "icons"

# Filled pills, >= 4.5:1 contrast with white text (same palette as Pipeline
# Studio's status pills).
_PILL = {
    "built": ("Built", "#2a7f3f"),
    "missing": ("Not built", "#5f6b7a"),
    "checking": ("Checking…", "#3b5bdb"),
    "busy": ("Working…", "#1864ab"),
    "error": ("Error", "#c92a2a"),
}
_FOUND, _MISSING = "#2a7f3f", "#c92a2a"


def _pill(label: QLabel, state: str) -> None:
    text, colour = _PILL[state]
    label.setText(text)
    label.setStyleSheet(
        f"QLabel {{ background: {colour}; color: white; border-radius: 8px;"
        " padding: 1px 8px; font-weight: 600; font-size: 11px; }")


class _SolverRow(QWidget):
    def __init__(self, spec: sb.SolverSpec) -> None:
        super().__init__()
        # transparent so the list's selection colour shows through
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet("_SolverRow, QLabel { background: transparent; }")
        h = QHBoxLayout(self)
        h.setContentsMargins(6, 5, 6, 5)
        v = QVBoxLayout()
        v.setSpacing(0)
        title = QLabel(spec.label)
        title.setStyleSheet("font-weight: 600;")
        sub = QLabel(spec.summary)
        sub.setObjectName("InfoLabel")
        sub.setWordWrap(True)
        sub.setSizePolicy(QSizePolicy.Policy.Ignored,
                          QSizePolicy.Policy.Preferred)
        v.addWidget(title)
        v.addWidget(sub)
        h.addLayout(v, 1)
        self.pill = QLabel()
        self.pill.setSizePolicy(QSizePolicy.Policy.Fixed,
                                QSizePolicy.Policy.Fixed)
        h.addWidget(self.pill, 0, Qt.AlignmentFlag.AlignTop)
        _pill(self.pill, "missing")


class SolverBuilderWindow(QWidget):
    """Tools ▸ Solver Builder.

    Signals
    -------
    binary_built(str, str)
        ``(solver key, binary)`` after a successful build.
    """

    binary_built = Signal(str, str)

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (Qt.WindowType.Window | Qt.WindowType.WindowCloseButtonHint
                 | Qt.WindowType.WindowMinimizeButtonHint
                 | Qt.WindowType.WindowMaximizeButtonHint)
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Solver Builder")
        for name in ("tools.svg", "advanced-tools.svg"):
            if (_ICONS / name).exists():
                self.setWindowIcon(QIcon(str(_ICONS / name)))
                break
        self._key = "occam2d"
        self._readiness: dict[str, sb.Readiness] = {}
        self._worker = None
        self._t0 = 0.0
        self._rows: dict[str, _SolverRow] = {}
        self._build_ui()
        self._fit_to_screen(1040, 720)
        for key in sb.SOLVERS:  # every card shows Built / Not built at once
            self._update_head_pill(key)
        self._list.setCurrentRow(0)

    # ══ UI ════════════════════════════════════════════════════════════════

    def _build_ui(self) -> None:
        root = QHBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        root.addWidget(split)

        left = QWidget()
        left.setMinimumWidth(220)
        left.setMaximumWidth(300)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        t = QLabel("SOLVERS")
        t.setObjectName("PanelTitle")
        lv.addWidget(t)
        self._list = QListWidget()
        self._list.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        for key, spec in sb.SOLVERS.items():
            row = _SolverRow(spec)
            item = QListWidgetItem()
            item.setData(Qt.ItemDataRole.UserRole, key)
            from PySide6.QtCore import QSize

            item.setSizeHint(QSize(1, row.sizeHint().height()))
            self._list.addItem(item)
            self._list.setItemWidget(item, row)
            self._rows[key] = row
        self._list.currentRowChanged.connect(self._on_solver_selected)
        lv.addWidget(self._list, 1)
        where = QLabel(f"Toolchains are installed per user in:\n"
                       f"{sb.toolchain_root()}")
        where.setObjectName("InfoLabel")
        where.setWordWrap(True)
        lv.addWidget(where)
        split.addWidget(left)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        page = QWidget()
        self._page = QVBoxLayout(page)
        self._page.setContentsMargins(12, 8, 12, 8)
        self._page.setSpacing(10)
        self._build_header()
        self._build_source_group()
        self._build_deps_group()
        self._build_setup_group()
        self._build_actions()
        self._build_result()
        self._build_log()
        self._page.addStretch(1)
        scroll.setWidget(page)
        split.addWidget(scroll)
        split.setStretchFactor(1, 1)
        split.setSizes([250, 780])

        self._elapsed_timer = QTimer(self)
        self._elapsed_timer.setInterval(1000)
        self._elapsed_timer.timeout.connect(self._tick)

    def _build_header(self) -> None:
        h = QHBoxLayout()
        v = QVBoxLayout()
        self._title = QLabel()
        self._title.setObjectName("SolverTitle")
        self._title.setStyleSheet(
            "QLabel#SolverTitle { font-size: 17px; font-weight: 600; }")
        self._summary = QLabel()
        self._summary.setWordWrap(True)
        self._license = QLabel()
        self._license.setObjectName("InfoLabel")
        self._license.setWordWrap(True)
        v.addWidget(self._title)
        v.addWidget(self._summary)
        v.addWidget(self._license)
        h.addLayout(v, 1)
        self._head_pill = QLabel()
        h.addWidget(self._head_pill, 0, Qt.AlignmentFlag.AlignTop)
        self._page.addLayout(h)

    def _build_source_group(self) -> None:
        g = QGroupBox("Source code")
        v = QVBoxLayout(g)
        self._rb_bundled = QRadioButton()
        self._rb_custom = QRadioButton("Use my own source folder")
        grp = QButtonGroup(self)
        grp.addButton(self._rb_bundled)
        grp.addButton(self._rb_custom)
        self._rb_bundled.setChecked(True)
        v.addWidget(self._rb_bundled)
        self._bundled_path = QLabel()
        self._bundled_path.setObjectName("InfoLabel")
        self._bundled_path.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)
        self._bundled_path.setContentsMargins(22, 0, 0, 0)
        v.addWidget(self._bundled_path)
        v.addWidget(self._rb_custom)
        row = QHBoxLayout()
        row.setContentsMargins(22, 0, 0, 0)
        self._src_edit = QLineEdit()
        self._src_edit.setPlaceholderText(
            "Folder containing the solver's Makefile and sources")
        browse = compact_button(QPushButton("…"), 28)
        browse.setToolTip("Browse for the source folder")
        browse.clicked.connect(self._browse_source)
        row.addWidget(self._src_edit, 1)
        row.addWidget(browse)
        v.addLayout(row)
        self._rb_custom.toggled.connect(self._on_source_mode)
        self._src_edit.editingFinished.connect(self._recheck)
        self._on_source_mode(False)
        self._page.addWidget(g)

    def _build_deps_group(self) -> None:
        g = QGroupBox("Dependencies")
        v = QVBoxLayout(g)
        self._deps = QTreeWidget()
        self._deps.setColumnCount(3)
        self._deps.setHeaderLabels(["Status", "Requirement", "Details"])
        self._deps.setRootIsDecorated(False)
        self._deps.setUniformRowHeights(True)
        hdr = self._deps.header()
        hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hdr.setStretchLastSection(True)
        self._deps.setMinimumHeight(150)
        v.addWidget(self._deps)
        row = QHBoxLayout()
        self._deps_summary = QLabel("Checking…")
        self._deps_summary.setWordWrap(True)
        row.addWidget(self._deps_summary, 1)
        self._btn_recheck = QPushButton("Check again")
        self._btn_recheck.clicked.connect(self._recheck)
        row.addWidget(self._btn_recheck)
        v.addLayout(row)
        self._page.addWidget(g)

    def _build_setup_group(self) -> None:
        g = QGroupBox("Setup")
        v = QVBoxLayout(g)
        self._rb_auto = QRadioButton(
            "Install missing dependencies automatically (recommended) — "
            "private toolchain, no administrator rights")
        self._rb_manual = QRadioButton("I will install them myself")
        grp = QButtonGroup(self)
        grp.addButton(self._rb_auto)
        grp.addButton(self._rb_manual)
        self._rb_auto.setChecked(True)
        v.addWidget(self._rb_auto)
        v.addWidget(self._rb_manual)
        self._manual = QLabel()
        self._manual.setWordWrap(True)
        self._manual.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)
        self._manual.setObjectName("InfoLabel")
        v.addWidget(self._manual)
        self._rb_auto.toggled.connect(lambda _on: self._update_state())
        self._page.addWidget(g)

    def _build_actions(self) -> None:
        box = QFrame()
        box.setObjectName("BuildCard")
        v = QVBoxLayout(box)
        v.setContentsMargins(0, 0, 0, 0)
        row = QHBoxLayout()
        self._chk_clean = QCheckBox("Clean build (recompile everything)")
        row.addWidget(self._chk_clean)
        row.addStretch(1)
        self._btn_install = QPushButton("Install dependencies")
        self._btn_install.clicked.connect(lambda: self._start("install"))
        self._btn_build = QPushButton("Build")
        self._btn_build.setObjectName("SolverBuildButton")
        self._btn_build.setMinimumWidth(170)
        # the one primary action on the page: accent-filled
        self._btn_build.setStyleSheet(
            "QPushButton#SolverBuildButton { background: #1864ab; "
            "color: white; font-weight: 600; border: none; "
            "border-radius: 6px; padding: 6px 16px; }"
            "QPushButton#SolverBuildButton:hover { background: #1a73c2; }"
            "QPushButton#SolverBuildButton:disabled { background: #9aa7b5; }")
        self._btn_build.clicked.connect(lambda: self._start("build"))
        self._btn_stop = QPushButton("Stop")
        self._btn_stop.clicked.connect(self._stop)
        for b in (self._btn_install, self._btn_build, self._btn_stop):
            row.addWidget(b)
        v.addLayout(row)
        prow = QHBoxLayout()
        self._stage = QLabel("")
        f = QFont(self._stage.font())
        f.setBold(True)
        self._stage.setFont(f)
        prow.addWidget(self._stage, 1)
        self._elapsed = QLabel("")
        self._elapsed.setObjectName("InfoLabel")
        prow.addWidget(self._elapsed)
        v.addLayout(prow)
        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setTextVisible(True)
        self._progress.setVisible(False)
        v.addWidget(self._progress)
        self._page.addWidget(box)

    def _build_result(self) -> None:
        self._result = QGroupBox("Result")
        v = QVBoxLayout(self._result)
        self._result_msg = QLabel()
        self._result_msg.setWordWrap(True)
        self._result_msg.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse)
        v.addWidget(self._result_msg)
        row = QHBoxLayout()
        self._btn_copy = QPushButton("Copy path")
        self._btn_copy.clicked.connect(self._copy_path)
        self._btn_open = QPushButton("Open folder")
        self._btn_open.clicked.connect(self._open_folder)
        row.addWidget(self._btn_copy)
        row.addWidget(self._btn_open)
        row.addStretch(1)
        v.addLayout(row)
        self._result.setVisible(False)
        self._page.addWidget(self._result)

    def _build_log(self) -> None:
        self._btn_log = QPushButton("▸  Show details")
        self._btn_log.setCheckable(True)
        self._btn_log.setFlat(True)
        self._btn_log.toggled.connect(self._toggle_log)
        self._page.addWidget(self._btn_log, 0, Qt.AlignmentFlag.AlignLeft)
        self._log = QPlainTextEdit()
        self._log.setReadOnly(True)
        self._log.setMaximumBlockCount(20000)
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.StyleHint.Monospace)
        self._log.setFont(mono)
        self._log.setMinimumHeight(220)
        self._log.setVisible(False)
        self._page.addWidget(self._log)

    def _fit_to_screen(self, w: int, h: int) -> None:
        screen = self.screen() or QGuiApplication.primaryScreen()
        if screen is not None:
            a = screen.availableGeometry()
            w, h = min(w, int(a.width() * 0.9)), min(h, int(a.height() * 0.9))
        self.resize(w, h)

    # ══ Public API ════════════════════════════════════════════════════════

    def open_solver(self, key: str) -> None:
        keys = list(sb.SOLVERS)
        if key in keys:
            self._list.setCurrentRow(keys.index(key))

    # ══ Selection & checks ════════════════════════════════════════════════

    def _on_solver_selected(self, row: int) -> None:
        keys = list(sb.SOLVERS)
        if not 0 <= row < len(keys):
            return
        self._key = keys[row]
        spec = sb.SOLVERS[self._key]
        self._title.setText(spec.label)
        self._summary.setText(spec.summary)
        self._license.setText(spec.license_note)
        self._license.setVisible(bool(spec.license_note))
        self._btn_build.setText(f"Build {spec.label}")
        src = spec.default_source_dir()
        if spec.vendored:
            self._rb_bundled.setText("Use pyCSAMT's bundled source "
                                     "(recommended)")
            self._bundled_path.setText(str(src))
        else:
            bundled = (src / "Makefile").is_file()
            self._rb_bundled.setText(
                "Use the local copy" if bundled else
                "Download automatically from bitbucket.org (recommended)")
            self._bundled_path.setText(
                str(src) if bundled else "bitbucket.org/mare2dem — the "
                "source is not bundled with pyCSAMT (separate license).")
        self._result.setVisible(False)
        self._update_head_pill()
        cached = self._readiness.get(self._key)
        if cached is not None:
            self._show_readiness(cached)
        else:
            self._recheck()

    def _source_dir(self):
        if self._rb_custom.isChecked() and self._src_edit.text().strip():
            return Path(self._src_edit.text().strip())
        return None

    def _recheck(self) -> None:
        if self._busy():
            return
        self._deps.clear()
        self._deps_summary.setText("Checking dependencies…")
        _pill(self._rows[self._key].pill, "checking")
        from pycsamt.app.desktop.workers.solver_build_worker import (
            SolverBuildWorker,
        )

        w = SolverBuildWorker(self._key, mode="check",
                              source_dir=self._source_dir(), parent=self)
        key = self._key
        w.checks_ready.connect(lambda r, k=key: self._on_checks(k, r))
        w.finished.connect(lambda: self._update_state())
        self._worker = w
        self._update_state()
        w.start()

    def _on_checks(self, key: str, readiness) -> None:
        self._readiness[key] = readiness
        self._update_head_pill(key)
        if key == self._key:
            self._show_readiness(readiness)
        self._worker = None
        self._update_state()

    def _show_readiness(self, r: sb.Readiness) -> None:
        self._deps.clear()
        for c in r.checks:
            item = QTreeWidgetItem(["Found" if c.ok else "Missing", c.label,
                                    c.detail])
            item.setForeground(0, Qt.GlobalColor.white)
            item.setBackground(0, Qt.GlobalColor.transparent)
            from PySide6.QtGui import QBrush, QColor

            item.setBackground(0, QBrush(QColor(_FOUND if c.ok
                                               else _MISSING)))
            item.setToolTip(2, c.detail)
            if not c.ok and c.hint:
                item.setToolTip(1, c.hint)
            self._deps.addTopLevelItem(item)
        miss = r.missing
        if not miss:
            self._deps_summary.setText("✓ All dependencies found — ready "
                                       "to build.")
        elif r.can_auto_install:
            self._deps_summary.setText(
                f"{len(miss)} missing — they can be installed automatically "
                "(Build does it for you), or install them yourself.")
        else:
            self._deps_summary.setText(
                f"{len(miss)} missing, including some that must be "
                "installed manually (see Setup).")
        hints = []
        for c in miss:
            if c.hint and c.hint not in hints:
                hints.append(c.hint)
        self._manual.setText("<br>".join(f"• {h}" for h in hints))
        self._update_state()

    def _update_head_pill(self, key: str | None = None) -> None:
        key = key or self._key
        built = sb.find_binary(key, None)
        state = "built" if built else "missing"
        _pill(self._rows[key].pill, state)
        if key == self._key:
            _pill(self._head_pill, state)
            self._head_pill.setToolTip(built or "No binary built yet")

    # ══ Actions ═══════════════════════════════════════════════════════════

    def _busy(self) -> bool:
        # The worker reference is the busy flag (cleared in _on_checks /
        # _on_done). isRunning() is still False right after creation, which
        # left Build clickable and Stop disabled during a build.
        return self._worker is not None

    def _update_state(self) -> None:
        busy = self._busy()
        r = self._readiness.get(self._key)
        auto = self._rb_auto.isChecked()
        self._manual.setVisible(bool(r and r.missing) and (
            not auto or not r.can_auto_install))
        can_build = bool(r) and (r.ready or (auto and r.can_auto_install))
        self._btn_build.setEnabled(not busy and can_build)
        self._btn_install.setEnabled(
            not busy and bool(r) and bool(r.missing) and auto
            and r.can_auto_install)
        self._btn_stop.setEnabled(busy and getattr(self._worker, "mode", "")
                                  != "check")
        self._btn_recheck.setEnabled(not busy)
        self._list.setEnabled(not busy)
        for w in (self._rb_bundled, self._rb_custom, self._src_edit,
                  self._rb_auto, self._rb_manual, self._chk_clean):
            w.setEnabled(not busy)

    def _start(self, mode: str) -> None:
        if self._busy():
            return
        from pycsamt.app.desktop.workers.solver_build_worker import (
            SolverBuildWorker,
        )

        self._log.clear()
        self._result.setVisible(False)
        self._progress.setVisible(True)
        self._progress.setRange(0, 0)
        self._stage.setText("Starting…")
        self._t0 = time.monotonic()
        self._elapsed_timer.start()
        self._tick()
        _pill(self._rows[self._key].pill, "busy")
        _pill(self._head_pill, "busy")
        w = SolverBuildWorker(self._key, mode=mode,
                              source_dir=self._source_dir(),
                              clean=self._chk_clean.isChecked(),
                              auto_install=self._rb_auto.isChecked(),
                              parent=self)
        w.log_line.connect(self._log.appendPlainText)
        w.stage.connect(self._stage.setText)
        w.progress.connect(self._on_progress)
        w.done.connect(self._on_done)
        self._worker = w
        self._update_state()
        w.start()

    def _on_progress(self, value: int) -> None:
        if value < 0:
            self._progress.setRange(0, 0)  # busy indicator
        else:
            self._progress.setRange(0, 100)
            self._progress.setValue(value)

    def _stop(self) -> None:
        if self._busy():
            self._worker.stop()
            self._stage.setText("Stopping…")

    def _on_done(self, ok: bool, message: str, binary: str) -> None:
        self._elapsed_timer.stop()
        self._tick()
        self._progress.setRange(0, 100)
        self._progress.setValue(100 if ok else self._progress.value())
        mode = self._worker.mode if self._worker else "build"
        self._worker = None
        if ok and binary:
            self._stage.setText("✓ Build complete")
            self._result_msg.setText(
                f"<b>{Path(binary.replace('wsl:', '')).name}</b> is ready:"
                f"<br><code>{binary}</code><br><br>The Inversion window will "
                "use it automatically (you can still browse to another "
                "binary there).")
            self._result.setVisible(True)
            self._last_binary = binary
            self._btn_open.setEnabled(not sb.is_wsl_binary(binary))
            self.binary_built.emit(self._key, binary)
        elif ok:
            self._stage.setText(f"✓ {message}")
        else:
            self._stage.setText(f"✕ {message}")
            _pill(self._rows[self._key].pill, "error")
            _pill(self._head_pill, "error")
            if not self._btn_log.isChecked():
                self._btn_log.setChecked(True)  # show what went wrong
        if ok:
            self._update_head_pill()
        if ok and mode != "check":
            # an install/build may have added dependencies: refresh the
            # checklist.  Not after a failure -- that would reset the card's
            # "Error" state to "Not built" and hide what just happened.
            self._readiness.pop(self._key, None)
            self._recheck()
        self._update_state()

    def _tick(self) -> None:
        s = int(time.monotonic() - self._t0) if self._t0 else 0
        self._elapsed.setText(f"elapsed {s // 60:02d}:{s % 60:02d}")

    # ══ Small helpers ═════════════════════════════════════════════════════

    def _on_source_mode(self, custom: bool) -> None:
        self._src_edit.setEnabled(custom)
        # back to the bundled source: its dependencies/readiness may differ
        if not custom and hasattr(self, "_deps") and self.isVisible():
            self._recheck()

    def _browse_source(self) -> None:
        path = QFileDialog.getExistingDirectory(self, "Solver source folder")
        if path:
            self._src_edit.setText(path)
            self._rb_custom.setChecked(True)
            self._recheck()

    def _copy_path(self) -> None:
        QGuiApplication.clipboard().setText(getattr(self, "_last_binary", ""))

    def _open_folder(self) -> None:
        b = getattr(self, "_last_binary", "")
        if b and not sb.is_wsl_binary(b):
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(Path(b).parent)))

    def _toggle_log(self, on: bool) -> None:
        self._log.setVisible(on)
        self._btn_log.setText(("▾" if on else "▸") + "  Show details")

    # ══ Session ═══════════════════════════════════════════════════════════

    def save_geometry_to(self, store: dict) -> None:
        store["solver_builder"] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
        }

    def restore_geometry_from(self, store: dict) -> None:
        entry = store.get("solver_builder") or {}
        geo = entry.get("geometry")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass

    def closeEvent(self, event) -> None:
        if self._busy():
            event.ignore()  # a build keeps running; just hide the window
            self.hide()
            return
        self.hide()
        event.ignore()


__all__ = ["SolverBuilderWindow"]
