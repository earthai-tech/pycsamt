# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
InterpretationWindow — the Interpretation Studio (desktop v2.6).

┌ Interpretation Studio ────────────────────────────────────────────────────┐
│ [Geology|Structure|Hydrology|Monitoring|Uncertainty|Diagnostics]          │
│                                   Station [S12 ▾]  Export ▾  Library ▸   │
├──────────────┬──────────────────────────────────────────────┬─────────────┤
│ EVIDENCE     │ Classify & calibrate  [▶ Run]   ✓ 47 stations │ PINNED      │
│ Model  ✓     ├──────────────┬───────────────────────────────┤ [thumb]     │
│ Boreholes    │ VIEWS        │                               │ [thumb]     │
│ Structure    │ ● Model      │   figure (or "why not" card)  │             │
│ Rock DB      │ ○ Strat log  │                               │             │
│ Monitoring   │ SETTINGS     │                               │             │
├──────────────┴──────────────┴───────────────────────────────┴─────────────┤
│ Console                                                                   │
└───────────────────────────────────────────────────────────────────────────┘

The tabs, their steps, settings and views are declared in
:mod:`pycsamt.app.desktop.controllers.interp_studio`; every computation
and plot is :class:`~pycsamt.app.desktop.controllers.interp_controller.
InterpController`'s.  Each view shows whether it can draw now (●) or what
it still needs (○, with the reason on hover and on its card).  Figures
are publication-white; pinned ones go to the Library for side-by-side
review and batch export.
"""

from __future__ import annotations

import io
from pathlib import Path

from PySide6.QtCore import QByteArray, QSize, Qt, QThread, Signal
from PySide6.QtGui import QColor, QIcon, QKeySequence, QPixmap, QShortcut
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QPushButton,
    QScrollArea,
    QSplitter,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.controllers.interp_controller import InterpController
from pycsamt.app.desktop.controllers.interp_studio import (
    TABS,
    StudioView,
    view_status,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import _icon, make_group
from pycsamt.app.desktop.windows.interp.evidence import EvidencePanel
from pycsamt.app.desktop.windows.inversion.console import (
    ConsolePanel,
    ConsoleWindow,
)
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm
from pycsamt.app.desktop.windows.inversion.window import _segmented, _vline

_KEY_ROLE = Qt.ItemDataRole.UserRole
_READY = QColor("#2a7f3f")
_MISSING = QColor("#8a94a3")


class _Task(QThread):
    """Runs one controller computation off the GUI thread."""

    done = Signal(object)  # the return value, or the exception
    def __init__(self, fn, parent=None) -> None:
        super().__init__(parent)
        self._fn = fn

    def run(self) -> None:
        try:
            out = self._fn()
        except Exception as exc:  # reported in the console
            out = exc
        self.done.emit(out)


def _dot(colour: QColor) -> QIcon:
    pix = QPixmap(10, 10)
    pix.fill(Qt.GlobalColor.transparent)
    from PySide6.QtGui import QPainter

    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setBrush(colour)
    p.setPen(Qt.PenStyle.NoPen)
    p.drawEllipse(1, 1, 8, 8)
    p.end()
    return QIcon(pix)


def _thumbnail(fig) -> QIcon:
    try:
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=28, facecolor="white")
        pix = QPixmap()
        pix.loadFromData(buf.getvalue())
        return QIcon(pix.scaled(120, 80, Qt.AspectRatioMode.KeepAspectRatio,
                                Qt.TransformationMode.SmoothTransformation))
    except Exception:
        return QIcon()


class InterpretationWindow(QWidget):
    """Interpretation Studio: evidence, interpretation steps, views."""

    panel_closed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        flags = (Qt.WindowType.Window | Qt.WindowType.WindowCloseButtonHint
                 | Qt.WindowType.WindowMinimizeButtonHint
                 | Qt.WindowType.WindowMaximizeButtonHint)
        super().__init__(parent, flags)
        self.setWindowTitle("pycsamt — Interpretation Studio")
        ic = _icon("interpret")
        if not ic.isNull():
            self.setWindowIcon(ic)
        self.resize(1380, 860)
        self._session_key = "interpretation"
        self._sites = None
        self._dark = False
        self._ctrl = InterpController()
        self._ctrl.dark = False  # figures are for publication: white
        self._task: _Task | None = None
        self._task_name = ""
        self._tab_key = TABS[0].key
        self._forms: dict[str, SettingsForm] = {}
        self._view_row: dict[str, int] = {}  # last view per tab
        self._figure = None  # figure on the canvas
        self._pinned: list = []  # (label, Figure)
        self._build_ui()
        self._select_tab(0)

    # ══ UI ═══════════════════════════════════════════════════════════════
    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)
        root.addWidget(self._build_header())

        self._vsplit = QSplitter(Qt.Orientation.Vertical)
        self._vsplit.setChildrenCollapsible(False)
        body = QSplitter(Qt.Orientation.Horizontal)
        body.setChildrenCollapsible(False)

        self._evidence = EvidencePanel(self._ctrl)
        self._evidence.model_requested.connect(self.load_model)
        self._evidence.inversion_requested.connect(self._from_inversion)
        self._evidence.changed.connect(self._on_evidence_changed)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(self._evidence)
        scroll.setMinimumWidth(250)
        scroll.setMaximumWidth(340)
        body.addWidget(scroll)
        body.addWidget(self._build_center())
        self._lib_panel = self._build_library()
        body.addWidget(self._lib_panel)
        body.setStretchFactor(1, 1)
        body.setSizes([280, 900, 200])
        self._vsplit.addWidget(body)

        self._console = ConsolePanel()
        self._console.pop_out_requested.connect(self._pop_out_console)
        self._console_win = ConsoleWindow(self)
        self._console_win.setWindowTitle("Interpretation console")
        self._console_win.dock_requested.connect(self._dock_console)
        self._vsplit.addWidget(self._console)
        self._vsplit.setStretchFactor(0, 4)
        self._vsplit.setStretchFactor(1, 1)
        self._vsplit.setSizes([680, 130])
        root.addWidget(self._vsplit, 1)
        QShortcut(QKeySequence("Ctrl+L"), self,
                  activated=lambda: self._btn_library.toggle())
        QShortcut(QKeySequence("Ctrl+R"), self, activated=self._run_step)

    def _build_header(self) -> QWidget:
        bar = QWidget()
        h = QHBoxLayout(bar)
        h.setContentsMargins(2, 0, 2, 0)
        h.setSpacing(8)
        box, self._tab_group, self._tab_buttons = _segmented(
            tuple(t.label for t in TABS), bar)
        self._tab_group.idClicked.connect(self._select_tab)
        h.addWidget(box)
        h.addStretch(1)
        self._data_lbl = QLabel("")
        self._data_lbl.setObjectName("InfoLabel")
        h.addWidget(self._data_lbl)
        h.addWidget(_vline())
        h.addWidget(QLabel("Station:"))
        self._station_combo = QComboBox()
        self._station_combo.setMinimumWidth(110)
        self._station_combo.setToolTip("Station used by per-station views "
                                       "(logs, depth profiles, LAS export)")
        self._station_combo.activated.connect(self._on_station)
        h.addWidget(self._station_combo)
        self._btn_export = QToolButton()
        self._btn_export.setText("Export  ▾")
        self._btn_export.setPopupMode(
            QToolButton.ToolButtonPopupMode.InstantPopup)
        menu = QMenu(self._btn_export)
        menu.addAction("Current figure…", self._export_current_figure)
        menu.addAction("Pinned figures to a folder…", self._export_pinned)
        menu.addSeparator()
        menu.addAction("Logs — Oasis Montaj XYZ…",
                       lambda: self._export("xyz", "XYZ (*.xyz)"))
        menu.addAction("Log — LAS 2.0 (station)…",
                       lambda: self._export("las", "LAS (*.las)"))
        menu.addAction("Logs — CSV table…",
                       lambda: self._export("csv", "CSV (*.csv)"))
        menu.addAction("Model — VTK…",
                       lambda: self._export("vtk", "VTK (*.vtk)"))
        self._btn_export.setMenu(menu)
        h.addWidget(self._btn_export)
        self._btn_library = QToolButton()
        self._btn_library.setCheckable(True)
        self._btn_library.setChecked(True)
        self._btn_library.setText("Library ▸")
        self._btn_library.setToolTip("Show / hide pinned figures (Ctrl+L)")
        self._btn_library.toggled.connect(
            lambda on: self._lib_panel.setVisible(on))
        h.addWidget(self._btn_library)
        return bar

    def _build_center(self) -> QWidget:
        center = QWidget()
        cv = QVBoxLayout(center)
        cv.setContentsMargins(0, 0, 0, 0)
        cv.setSpacing(4)

        step = QWidget()
        step.setObjectName("StepBar")
        sh = QHBoxLayout(step)
        sh.setContentsMargins(0, 0, 0, 0)
        self._step_desc = QLabel("")
        self._step_desc.setObjectName("InfoLabel")
        self._step_desc.setWordWrap(True)
        sh.addWidget(self._step_desc, 1)
        self._step_status = QLabel("")
        self._step_status.setObjectName("InfoLabel")
        sh.addWidget(self._step_status)
        self._btn_step = QPushButton("▶  Run")
        self._btn_step.setToolTip("Run this tab's step (Ctrl+R)")
        self._btn_step.clicked.connect(self._run_step)
        sh.addWidget(self._btn_step)
        cv.addWidget(step)

        split = QSplitter(Qt.Orientation.Horizontal)
        split.setChildrenCollapsible(False)
        left = QWidget()
        left.setMinimumWidth(200)
        left.setMaximumWidth(300)
        lv = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.setSpacing(4)
        cap = QLabel("VIEWS")
        cap.setObjectName("InfoLabel")
        lv.addWidget(cap)
        self._view_list = QListWidget()
        self._view_list.setObjectName("EngineList")
        self._view_list.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._view_list.setTextElideMode(Qt.TextElideMode.ElideRight)
        self._view_list.currentItemChanged.connect(
            lambda cur, _prev: cur is not None and self._render(
                cur.data(_KEY_ROLE)))
        lv.addWidget(self._view_list, 2)
        self._grp_settings, gl = make_group("Settings")
        self._form_stack = QStackedWidget()
        self._no_settings = QLabel("No settings for this tab.")
        self._no_settings.setObjectName("InfoLabel")
        self._form_stack.addWidget(self._no_settings)
        for t in TABS:
            if t.fields:
                form = SettingsForm(list(t.fields))
                self._forms[t.key] = form
                self._form_stack.addWidget(form)
        gl.addWidget(self._form_stack)
        lv.addWidget(self._grp_settings)
        split.addWidget(left)

        right = QWidget()
        rv = QVBoxLayout(right)
        rv.setContentsMargins(0, 0, 0, 0)
        rv.setSpacing(2)
        self._view = CanvasResultView(
            self, empty_title="Load a model to begin",
            empty_reason="Open a PCSF/PCSM file or an inversion run folder "
                         "(Evidence ▸ Model ▸ Load model).",
            empty_guidance="EDI/XML data loaded in the main window feed "
                           "the Structure and Diagnostics views.")
        rv.addWidget(self._view, 1)
        tools = QHBoxLayout()
        tools.addStretch(1)
        self._btn_refresh = QPushButton("⟳  Redraw")
        self._btn_refresh.setObjectName("FileListBtn")
        self._btn_refresh.clicked.connect(self._rerender)
        self._btn_pin = QPushButton("📌  Pin to library")
        self._btn_pin.setObjectName("FileListBtn")
        self._btn_pin.clicked.connect(self._pin_current)
        tools.addWidget(self._btn_refresh)
        tools.addWidget(self._btn_pin)
        rv.addLayout(tools)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        split.setSizes([230, 700])
        cv.addWidget(split, 1)
        return center

    def _build_library(self) -> QWidget:
        w = QWidget()
        w.setMinimumWidth(170)
        w.setMaximumWidth(260)
        v = QVBoxLayout(w)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        cap = QLabel("PINNED FIGURES")
        cap.setObjectName("InfoLabel")
        v.addWidget(cap)
        self._gallery = QListWidget()
        self._gallery.setObjectName("GalleryList")
        self._gallery.setViewMode(QListWidget.ViewMode.IconMode)
        self._gallery.setIconSize(QSize(120, 80))
        self._gallery.setResizeMode(QListWidget.ResizeMode.Adjust)
        self._gallery.setMovement(QListWidget.Movement.Static)
        self._gallery.setWordWrap(True)
        self._gallery.setSpacing(4)
        self._gallery.itemClicked.connect(self._show_pinned)
        v.addWidget(self._gallery, 1)
        row = QHBoxLayout()
        b = QPushButton("Remove")
        b.setObjectName("FileListBtn")
        b.clicked.connect(self._unpin)
        row.addWidget(b)
        b = QPushButton("Export all…")
        b.setObjectName("FileListBtn")
        b.clicked.connect(self._export_pinned)
        row.addWidget(b)
        v.addLayout(row)
        return w

    # ══ tabs & views ═════════════════════════════════════════════════════
    @property
    def _tab(self):
        return next(t for t in TABS if t.key == self._tab_key)

    def _select_tab(self, i: int) -> None:
        t = TABS[i]
        if t.key != self._tab_key and self._view_list.currentRow() >= 0:
            self._view_row[self._tab_key] = self._view_list.currentRow()
        self._tab_key = t.key
        self._tab_buttons[i].setChecked(True)
        self._step_desc.setText(t.step_help or
                                "Views of the data and the model; nothing "
                                "to compute on this tab.")
        self._btn_step.setVisible(bool(t.step))
        # (button text keeps "&&": a single "&" is a Qt mnemonic)
        self._btn_step.setText(f"▶  {t.step_label}" if t.step else "▶  Run")
        form = self._forms.get(t.key)
        self._form_stack.setCurrentWidget(form or self._no_settings)
        self._grp_settings.setVisible(form is not None)
        self._update_step_status()
        self._fill_views(select=self._view_row.get(t.key, 0))

    def _fill_views(self, select: int | None = None) -> None:
        row = self._view_list.currentRow() if select is None else select
        self._view_list.blockSignals(True)
        self._view_list.clear()
        for v in self._tab.views:
            ok, why = view_status(self._ctrl, v)
            item = QListWidgetItem(_dot(_READY if ok else _MISSING), v.label)
            item.setData(_KEY_ROLE, v.method)
            item.setToolTip(v.help + ("" if ok else f"\n\nNeeds: {why}"))
            if not ok:
                item.setForeground(_MISSING)
            self._view_list.addItem(item)
        self._view_list.blockSignals(False)
        if self._view_list.count():
            row = min(max(row, 0), self._view_list.count() - 1)
            self._view_list.setCurrentRow(row)
            self._render(self._view_list.item(row).data(_KEY_ROLE))

    def _studio_view(self, method: str) -> StudioView | None:
        for t in TABS:
            for v in t.views:
                if v.method == method:
                    return v
        return None

    def _render(self, method: str) -> None:
        v = self._studio_view(method)
        if v is None:
            return
        ok, why = view_status(self._ctrl, v)
        if not ok:
            self._show_card(v.label, why)
            return
        kw = {}
        if v.per_station and self._station_combo.currentText():
            kw["station"] = self._station_combo.currentText()
        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            fig = self._ctrl.generate(method, **kw)
        finally:
            self.unsetCursor()
        msg = figure_blank_reason(fig)
        if msg is not None:
            self._close_figure(fig)
            self._show_card(v.label, msg.replace("⚠", "").replace(
                "✕", "").strip())
            return
        self._close_figure(self._figure)
        self._figure = fig
        self._figure_label = v.label
        self._view.canvas.show_figure(fig)
        self._view.show_canvas()
        self._btn_pin.setEnabled(True)

    def _show_card(self, title: str, reason: str) -> None:
        self._close_figure(self._figure)
        self._figure = None
        self._btn_pin.setEnabled(False)
        self._view.show_unavailable(title, reason)

    def _close_figure(self, fig) -> None:
        if fig is None or any(f is fig for _l, f in self._pinned):
            return
        import matplotlib.pyplot as plt

        plt.close(fig)

    def _rerender(self) -> None:
        item = self._view_list.currentItem()
        if item is not None:
            self._render(item.data(_KEY_ROLE))

    def current_view(self) -> str:
        item = self._view_list.currentItem()
        return item.data(_KEY_ROLE) if item is not None else ""

    def select_view(self, method: str) -> None:
        """Switch to the tab holding *method* and show it."""
        for i, t in enumerate(TABS):
            keys = [v.method for v in t.views]
            if method in keys:
                self._view_row[t.key] = keys.index(method)
                self._select_tab(i)
                return

    def _on_station(self, _i: int) -> None:
        v = self._studio_view(self.current_view())
        if v is not None and v.per_station:
            self._rerender()

    # ══ steps ════════════════════════════════════════════════════════════
    def _step_call(self, key: str):
        """The controller call for tab *key*'s step, with its settings."""
        c = self._ctrl
        vals = self._forms[key].values() if key in self._forms else {}
        if key == "geology":
            return c.run_geological
        if key == "hydrology":
            def hydro():
                c.set_petro_config(**vals)
                return c.run_hydro()
            return hydro
        if key == "uncertainty":
            petro = self._forms["hydrology"].values()

            def mc():
                c.set_petro_config(**petro)
                return c.run_monte_carlo(
                    n_samples=int(vals["n_samples"]),
                    rho_w_range=(vals["rho_w_lo"], vals["rho_w_hi"]),
                    m_range=(vals["m_lo"], vals["m_hi"]),
                    n_range=(vals["n_lo"], vals["n_hi"]),
                    phi_range=(vals["phi_lo"], vals["phi_hi"]))
            return mc
        if key == "monitoring":
            return lambda: c.run_fusion(vals["primary_max_depth"],
                                        vals["secondary_min_depth"],
                                        vals["blend"])
        return None

    def _run_step(self) -> None:
        t = self._tab
        if not t.step or self._busy:
            return
        if self._ctrl.state.model is None:
            self._console.append("Load a model first (Evidence ▸ Model).")
            self._step_status.setText("Load a model first")
            return
        label = t.step_label.replace("&&", "&")
        self._start_task(label, self._step_call(t.key), self._on_step_done)

    def _on_step_done(self, out) -> None:
        text = f"Error: {out}" if isinstance(out, Exception) else str(out)
        self._console.append(f"{self._task_name}: {text}")
        self._step_status.setText(text if len(text) < 70
                                  else text[:67] + "…")
        self._step_status.setToolTip(text)
        self._evidence.refresh()
        self._refresh_stations()
        self._fill_views()

    def _update_step_status(self) -> None:
        st = self._ctrl.state
        done = {"geology": f"✓ {len(st.strat_logs)} logs"
                if st.strat_logs else "",
                "hydrology": "✓ estimated" if st.hydro_result is not None
                else "",
                "uncertainty": "✓ ensemble ready"
                if st.mc_result is not None else "",
                "monitoring": "✓ fused" if st.fusion_model is not None
                else ""}.get(self._tab_key, "")
        self._step_status.setText(done)
        self._step_status.setToolTip("")

    # ══ background tasks ═════════════════════════════════════════════════
    @property
    def _busy(self) -> bool:
        return self._task is not None and self._task.isRunning()

    def _start_task(self, name: str, fn, on_done) -> None:
        self._task_name = name
        self._btn_step.setEnabled(False)
        self._evidence.btn_load.setEnabled(False)
        self._step_status.setText(f"{name}…")
        self._console.append(f"{name}…")
        self._task = _Task(fn, self)
        self._task.done.connect(on_done)
        self._task.finished.connect(self._task_finished)
        self._task.start()

    def _task_finished(self) -> None:
        self._btn_step.setEnabled(True)
        self._evidence.btn_load.setEnabled(True)

    def wait_for_task(self, ms: int = 600_000) -> bool:
        """Block until the running task ends (scripts and tests)."""
        if self._task is None:
            return True
        ok = self._task.wait(ms)
        from PySide6.QtWidgets import QApplication

        QApplication.processEvents()
        return ok

    # ══ model ════════════════════════════════════════════════════════════
    def load_model(self, path: str, line: str = "") -> None:
        """Load a PCSF/PCSM file or run folder (in the background)."""
        if self._busy:
            return
        self._evidence.set_model_path(path)
        self._start_task(f"Loading {Path(path).name}",
                         lambda: self._ctrl.load_model_source(
                             path, line=line or None),
                         self._on_model_loaded)

    def _on_model_loaded(self, out) -> None:
        if isinstance(out, Exception):
            self._console.append(f"Model not loaded: {out}")
            self._step_status.setText("Model not loaded")
            self._step_status.setToolTip(str(out))
            return
        self._console.append(f"Model: {out.source}")
        for note in out.notes:
            self._console.append(f"  {note}")
        self._update_status_card()
        self._fill_views()

    def set_model(self, model, source: str = "") -> None:
        """Use an in-memory ResistivityModel."""
        self._ctrl.set_model(model)
        self._ctrl.state.model_info = None
        if source:
            self._console.append(f"Model: {source}")
        self._update_status_card()
        self._fill_views()

    def receive_inversion(self, payload: dict) -> None:
        """A model sent by the Inversion Studio (``{"engine", "result",
        "path"}``)."""
        path = payload.get("path")
        if path:
            self.load_model(str(path))
            return
        result = payload.get("result")
        if result is None:
            return
        if hasattr(result, "rho_2d"):
            self.set_model(result, "from the Inversion Studio")
            return
        result_dir = getattr(result, "result_dir", None)
        if result_dir:
            self.load_model(str(result_dir))

    def _from_inversion(self) -> None:
        inv = getattr(self.parent(), "_inversion_win", None)
        run = getattr(inv, "_loaded", None) if inv is not None else None
        if run is not None and getattr(run, "path", None):
            self.load_model(str(run.path))
            return
        model = getattr(inv, "_result_model", None) if inv else None
        if model is not None:
            self.set_model(model, "from the Inversion Studio")
            return
        self._console.append("The Inversion Studio has no result open: run "
                             "an inversion or open a run in View Results.")

    def _on_evidence_changed(self, msg: str) -> None:
        self._console.append(msg)
        self._fill_views()

    # ══ library ══════════════════════════════════════════════════════════
    def _pin_current(self) -> None:
        if self._figure is None:
            return
        label = f"{self._tab.label} · {self._figure_label}"
        st = self._station_combo.currentText()
        v = self._studio_view(self.current_view())
        if v is not None and v.per_station and st:
            label += f" · {st}"
        self._pinned.append((label, self._figure))
        item = QListWidgetItem(_thumbnail(self._figure), label)
        item.setToolTip(label)
        self._gallery.addItem(item)
        self._btn_library.setChecked(True)
        self._console.append(f"Pinned: {label}")

    def _show_pinned(self, item: QListWidgetItem) -> None:
        i = self._gallery.row(item)
        if 0 <= i < len(self._pinned):
            label, fig = self._pinned[i]
            self._close_figure(self._figure)
            self._figure = fig
            self._figure_label = label
            self._view.canvas.show_figure(fig)
            self._view.show_canvas()

    def _unpin(self) -> None:
        i = self._gallery.currentRow()
        if 0 <= i < len(self._pinned):
            _label, fig = self._pinned.pop(i)
            self._gallery.takeItem(i)
            if fig is not self._figure:
                self._close_figure(fig)

    def pinned_labels(self) -> list[str]:
        return [label for label, _f in self._pinned]

    # ══ export ═══════════════════════════════════════════════════════════
    def _export_current_figure(self) -> None:
        if self._figure is None:
            self._console.append("No figure to export.")
            return
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog

        ExportDialog(figure=self._figure, parent=self).exec()

    def export_pinned(self, folder, fmt: str = "png", dpi: int = 300) -> list:
        """Save every pinned figure into *folder*; returns the paths."""
        out = []
        folder = Path(folder)
        folder.mkdir(parents=True, exist_ok=True)
        for i, (label, fig) in enumerate(self._pinned, 1):
            stem = "".join(ch if ch.isalnum() else "_" for ch in label)
            path = folder / f"{i:02d}_{stem.strip('_')}.{fmt}"
            fig.savefig(path, dpi=dpi, facecolor="white",
                        bbox_inches="tight")
            out.append(path)
        return out

    def _export_pinned(self) -> None:
        if not self._pinned:
            self._console.append("Pin figures first (📌 under the figure).")
            return
        d = QFileDialog.getExistingDirectory(self, "Export pinned figures")
        if d:
            paths = self.export_pinned(d)
            self._console.append(f"Exported {len(paths)} figure(s) to {d}")

    def _export(self, kind: str, filt: str) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, f"Export {kind.upper()}", "", f"{filt};;All files (*)")
        if not path:
            return
        fn = getattr(self._ctrl, f"export_{kind}")
        msg = (fn(path, station=self._station_combo.currentText())
               if kind == "las" else fn(path))
        self._console.append(msg)

    # ══ console ══════════════════════════════════════════════════════════
    def _pop_out_console(self) -> None:
        if self._console.parent() is self._console_win:
            self._console_win.raise_()
            return
        self._console_win.hold(self._console)
        self._console.btn_pop.setText("Dock")
        self._console.pop_out_requested.disconnect()
        self._console.pop_out_requested.connect(self._dock_console)
        self._console_win.show()

    def _dock_console(self) -> None:
        if self._console.parent() is not self._console_win:
            return
        self._vsplit.addWidget(self._console)
        self._console.btn_pop.setText("Pop out")
        self._console.pop_out_requested.disconnect()
        self._console.pop_out_requested.connect(self._pop_out_console)
        self._console_win.hide()
        self._vsplit.setSizes([680, 130])

    # ══ public API (main window) ═════════════════════════════════════════
    def set_sites(self, sites) -> None:
        self._sites = sites
        self._ctrl.set_sites(sites)
        self._update_status_card()
        self._fill_views()

    def set_dark_mode(self, dark: bool) -> None:
        # the UI follows the app theme; figures stay publication-white
        self._dark = dark
        self._ctrl.dark = False

    def _update_status_card(self) -> None:
        """Refresh everything that depends on the model / data."""
        self._evidence.refresh()
        self._refresh_stations()
        n = self._n_sites()
        self._data_lbl.setText(f"{n} EDI/XML stations" if n else
                               "No EDI/XML data")
        self._update_step_status()

    def _n_sites(self) -> int:
        sites = self._ctrl.state.sites
        if sites is None:
            return 0
        try:
            from pycsamt.emtools._core import _iter_items

            return len(list(_iter_items(sites)))
        except Exception:
            return 0

    def _station_names(self) -> list[str]:
        model = self._ctrl.state.model
        names = list(getattr(model, "station_names", None) or []) \
            if model is not None else []
        names += [b.name for b in self._ctrl.state.boreholes
                  if b.name not in names]
        if names:
            return [str(n) for n in names]
        sites = self._ctrl.state.sites
        if sites is None:
            return []
        try:
            from pycsamt.emtools._core import _iter_items, _name

            return [_name(ed, i) for i, ed in enumerate(_iter_items(sites))]
        except Exception:
            return []

    def _refresh_stations(self) -> None:
        cur = self._station_combo.currentText()
        names = self._station_names()
        self._station_combo.blockSignals(True)
        self._station_combo.clear()
        self._station_combo.addItems(names)
        if cur in names:
            self._station_combo.setCurrentText(cur)
        self._station_combo.blockSignals(False)

    # ── session ───────────────────────────────────────────────────────
    def save_geometry_to(self, store: dict) -> None:
        store[self._session_key] = {
            "geometry": self.saveGeometry().toBase64().data().decode(),
            "visible": self.isVisible(),
            "library_visible": self._btn_library.isChecked(),
            "tab": self._tab_key,
            "settings": {k: f.values() for k, f in self._forms.items()},
            "console_theme": self._console.view.theme,
        }

    def restore_geometry_from(self, store: dict) -> None:
        entry = store.get(self._session_key)
        if not entry:
            return
        geo = entry.get("geometry")
        if geo:
            try:
                self.restoreGeometry(QByteArray.fromBase64(geo.encode()))
            except Exception:
                pass
        if "library_visible" in entry:
            self._btn_library.setChecked(bool(entry["library_visible"]))
        for key, vals in (entry.get("settings") or {}).items():
            if key in self._forms and isinstance(vals, dict):
                self._forms[key].set_values(vals)
        if entry.get("console_theme"):
            self._console.set_theme(entry["console_theme"])
        keys = [t.key for t in TABS]
        if entry.get("tab") in keys:
            self._select_tab(keys.index(entry["tab"]))

    def closeEvent(self, event) -> None:  # noqa: N802
        """Hide instead of destroying: a step keeps running."""
        self.hide()
        event.ignore()
        self.panel_closed.emit()


__all__ = ["InterpretationWindow"]
