# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
TDEMWindow — time-domain EM studio.

Left
    **Data** — pick the source (TEMAVG survey folder, Geosoft, AMIRA,
    Zonge, WalkTEM, XYZ) and its acquisition settings, then load (in the
    background: a large TEMAVG survey takes seconds).
    **Selection** — a profile and its soundings.  Curve views draw only the
    checked soundings, section views only the chosen profile ("Whole
    survey" draws every profile).
    **Views** — soundings, profiles and survey views, each with its own
    options (:mod:`pycsamt.app.desktop.controllers.tdem_studio`).
    **Convert to EDI** — TEM -> impedance (method, frequency convention,
    phase, loop correction, transmitter waveform), then save the EDIs or
    send the sites to the main survey (append or replace, undoable).
Right
    **Plot** · **Soundings** (what was loaded) · **Converted** (the sites
    of the last conversion).
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QThread, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QPushButton,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import tdem_studio as st
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import PanelWindow, make_group
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm

_KEY = Qt.ItemDataRole.UserRole
_WHOLE = ""  # profile combo value: the whole survey


def _compact(combo: QComboBox) -> None:
    """Let a combo with long item texts shrink to the panel width."""
    combo.setSizeAdjustPolicy(
        QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    combo.setMinimumContentsLength(8)


def _show_page(stack: QStackedWidget, page: QWidget) -> None:
    """Show *page* and size the stack to it, not to its tallest page."""
    for i in range(stack.count()):
        w = stack.widget(i)
        pol = (QSizePolicy.Policy.Preferred if w is page
               else QSizePolicy.Policy.Ignored)
        w.setSizePolicy(pol, pol)
    stack.setCurrentWidget(page)
    stack.adjustSize()


class _LoadWorker(QThread):
    done = Signal(object)
    failed = Signal(str)

    def __init__(self, source, path, values, parent=None):
        super().__init__(parent)
        self._args = (source, path, values)

    def run(self) -> None:
        try:
            self.done.emit(st.load(*self._args))
        except Exception as exc:  # reported in the Data card
            self.failed.emit(str(exc))


class TDEMWindow(PanelWindow):
    """TDEM studio: load, select, view, convert to EDI."""

    #: converted sites and mode ("append" | "replace")
    send_to_survey = Signal(object, str)

    def __init__(self, parent: QWidget | None = None) -> None:
        self._data: st.TDEMData | None = None
        self._converted = None
        self._forms: dict[str, SettingsForm] = {}
        self._worker: _LoadWorker | None = None
        super().__init__(title="TDEM Studio", session_key="tdem_window",
                         params_width=300, icon_name="tdem", parent=parent)
        self.resize(1240, 820)
        self._on_source(0)
        self._refresh_data()

    # ══ left panel ═══════════════════════════════════════════════════════
    def _build_params(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Data")
        self._source_combo = QComboBox()
        for s in st.SOURCES:
            self._source_combo.addItem(s.label, s.key)
        self._source_combo.currentIndexChanged.connect(self._on_source)
        _compact(self._source_combo)
        lay.addWidget(self._source_combo)
        self._acq_stack = QStackedWidget()
        self._acq_forms: dict[str, SettingsForm] = {}
        for s in st.SOURCES:
            form = SettingsForm(list(s.fields)) if s.fields else QLabel(
                "Columns: site, time, data (and error).")
            self._acq_forms[s.key] = form
            self._acq_stack.addWidget(form)
        lay.addWidget(self._acq_stack)
        row = QHBoxLayout()
        self._btn_load = QPushButton("Load…")
        self._btn_load.clicked.connect(self._on_load)
        row.addWidget(self._btn_load, 1)
        clear = QPushButton("Clear")
        clear.clicked.connect(self._on_clear)
        row.addWidget(clear)
        lay.addLayout(row)
        self._data_status = QLabel("")
        self._data_status.setObjectName("InfoLabel")
        self._data_status.setWordWrap(True)
        self._data_status.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(self._data_status)
        layout.addWidget(grp)

        grp, lay = make_group("Selection")
        row = QHBoxLayout()
        row.addWidget(QLabel("Profile:"))
        self._profile_combo = QComboBox()
        self._profile_combo.setToolTip(
            "Sections draw this profile; the Soundings list shows its "
            "points")
        self._profile_combo.currentIndexChanged.connect(self._on_profile)
        _compact(self._profile_combo)
        row.addWidget(self._profile_combo, 1)
        lay.addLayout(row)
        self._snd_list = QListWidget()
        self._snd_list.setMinimumHeight(120)
        self._snd_list.itemChanged.connect(self._on_checks)
        lay.addWidget(self._snd_list)
        row = QHBoxLayout()
        for text, slot in (("All", lambda: self._check_all(True)),
                           ("None", lambda: self._check_all(False))):
            b = QPushButton(text)
            b.clicked.connect(slot)
            row.addWidget(b)
        lay.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel("Check every"))
        self._every = QSpinBox()
        self._every.setRange(1, 500)
        self._every.setValue(5)
        self._every.setToolTip("Check every n-th sounding of the profile")
        row.addWidget(self._every)
        b = QPushButton("Pick")
        b.clicked.connect(lambda: self.check_every(self._every.value()))
        row.addWidget(b)
        row.addStretch(1)
        lay.addLayout(row)
        self._sel_status = QLabel("")
        self._sel_status.setObjectName("InfoLabel")
        lay.addWidget(self._sel_status)
        layout.addWidget(grp)

        grp, lay = make_group("Views")
        self._view_list = QListWidget()
        self._view_list.setMinimumHeight(200)
        group = None
        for v in st.VIEWS:
            if v.group != group:
                group = v.group
                head = QListWidgetItem(group.upper())
                head.setFlags(Qt.ItemFlag.NoItemFlags)
                self._view_list.addItem(head)
            it = QListWidgetItem("   " + v.label)
            it.setData(_KEY, v.key)
            it.setToolTip(v.help)
            self._view_list.addItem(it)
        self._view_list.currentItemChanged.connect(self._on_view)
        lay.addWidget(self._view_list)
        self._view_desc = QLabel("")
        self._view_desc.setObjectName("InfoLabel")
        self._view_desc.setWordWrap(True)
        lay.addWidget(self._view_desc)
        layout.addWidget(grp)

        grp, lay = make_group("Options")
        self._form_stack = QStackedWidget()
        self._no_opts = QLabel("This view has no options.")
        self._no_opts.setObjectName("InfoLabel")
        self._form_stack.addWidget(self._no_opts)
        for v in st.VIEWS:
            if v.fields:
                form = SettingsForm(list(v.fields))
                form.changed.connect(self._schedule)
                self._forms[v.key] = form
                self._form_stack.addWidget(form)
        lay.addWidget(self._form_stack)
        layout.addWidget(grp)

        row = QHBoxLayout()
        self._btn_draw = QPushButton("▶  Draw")
        self._btn_draw.setToolTip("Render the selected view")
        self._btn_draw.clicked.connect(self._on_draw)
        row.addWidget(self._btn_draw, 1)
        self._chk_auto = QCheckBox("Auto")
        self._chk_auto.setChecked(True)
        self._chk_auto.setToolTip("Redraw when the view, selection or an "
                                  "option changes")
        row.addWidget(self._chk_auto)
        layout.addLayout(row)

        grp, lay = make_group("Convert to EDI")
        self._conv_form = SettingsForm(list(st.CONVERT_FIELDS))
        lay.addWidget(self._conv_form)
        self._btn_convert = QPushButton("Convert checked soundings")
        self._btn_convert.clicked.connect(self._on_convert)
        lay.addWidget(self._btn_convert)
        row = QHBoxLayout()
        self._btn_save = QPushButton("Save EDIs…")
        self._btn_save.clicked.connect(self._on_save)
        row.addWidget(self._btn_save)
        send = QToolButton()
        send.setText("Send to survey ▾")
        send.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        send.setToolTip("Hand the converted sites to the main window "
                        "(Edit ▸ Undo reverts it)")
        menu = QMenu(send)
        menu.addAction("Append to the survey",
                       lambda: self.send("append"))
        menu.addAction("Replace the survey",
                       lambda: self.send("replace"))
        send.setMenu(menu)
        self._btn_send = send
        row.addWidget(send, 1)
        lay.addLayout(row)
        layout.addWidget(grp)

        self._status = QLabel("")
        self._status.setObjectName("InfoLabel")
        self._status.setWordWrap(True)
        layout.addWidget(self._status)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(350)
        self._timer.timeout.connect(self._on_draw)

    # ══ right panel ══════════════════════════════════════════════════════
    def _build_content(self, layout: QVBoxLayout) -> None:
        self._tabs = QTabWidget()
        self._canvas_view = CanvasResultView(
            toolbar=True, empty_title="No TDEM data",
            empty_reason="Pick a source, then Load.",
            empty_guidance="A TEMAVG survey folder holds .AVG / .Z / .LOG "
                           "files, one per profile.")
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(self._on_draw,
                                          tooltip="Render the selected view")
        self._tabs.addTab(self._canvas_view, "Plot")
        self._snd_table = self._table()
        self._tabs.addTab(self._snd_table, "Soundings")
        self._conv_table = self._table()
        self._tabs.addTab(self._conv_table, "Converted")
        layout.addWidget(self._tabs)

    @staticmethod
    def _table() -> QTableWidget:
        t = QTableWidget()
        t.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        t.setAlternatingRowColors(True)
        return t

    # ══ data ═════════════════════════════════════════════════════════════
    def _source(self) -> st.Source:
        return st.source(self._source_combo.currentData())

    def _on_source(self, _i: int) -> None:
        _show_page(self._acq_stack, self._acq_forms[self._source().key])
        self._btn_load.setText("Load folder…" if self._source().folder
                               else "Load file…")

    def _acq_values(self) -> dict:
        form = self._acq_forms[self._source().key]
        return form.values() if isinstance(form, SettingsForm) else {}

    def load(self, path: str, source: str | None = None) -> int:
        """Load *path* synchronously; returns the sounding count."""
        if source:
            self.set_source(source)
        self._set_data(st.load(self._source().key, path,
                               self._acq_values()))
        return self._data.n

    def set_source(self, key: str) -> None:
        i = self._source_combo.findData(key)
        if i >= 0:
            self._source_combo.setCurrentIndex(i)

    def _on_load(self) -> None:
        src = self._source()
        if src.folder:
            path = QFileDialog.getExistingDirectory(
                self, "Select TEMAVG survey folder")
        else:
            path, _ = QFileDialog.getOpenFileName(
                self, f"Select {src.label}", "", src.filter)
        if not path:
            return
        self._data_status.setText(f"Loading {path} …")
        self._btn_load.setEnabled(False)
        self._worker = _LoadWorker(src.key, path, self._acq_values(), self)
        self._worker.done.connect(self._set_data)
        self._worker.failed.connect(self._on_load_failed)
        self._worker.finished.connect(lambda: self._btn_load.setEnabled(True))
        self._worker.start()

    def _on_load_failed(self, msg: str) -> None:
        self._data_status.setText(f"Load failed: {msg}")

    def _set_data(self, data) -> None:
        self._data = data
        self._converted = None
        self._refresh_data()

    def _on_clear(self) -> None:
        self._data = None
        self._converted = None
        self._refresh_data()

    @property
    def has_data(self) -> bool:
        return self._data is not None and self._data.n > 0

    def _refresh_data(self) -> None:
        d = self._data
        if not self.has_data:
            self._data_status.setText("No TDEM data loaded.")
        else:
            kind = next(s.label for s in st.SOURCES if s.key == d.source)
            self._data_status.setText(
                f"<b>{d.n} soundings</b> · {len(d.profiles)} profile"
                f"{'s' if len(d.profiles) != 1 else ''}<br>"
                f"<small>{kind}</small>")
        self._profile_combo.blockSignals(True)
        self._profile_combo.clear()
        if self.has_data:
            if d.survey is not None:
                self._profile_combo.addItem("Whole survey", _WHOLE)
            for name, idx in d.profiles.items():
                self._profile_combo.addItem(f"{name}  ({len(idx)})", name)
            if d.survey is not None and len(d.profiles) > 1:
                self._profile_combo.setCurrentIndex(1)
        self._profile_combo.blockSignals(False)
        self._fill_soundings_table()
        self._fill_converted_table()
        self._on_profile(self._profile_combo.currentIndex())
        self._update_actions()
        if not self.has_data:
            self._canvas_view.show_unavailable(
                "No TDEM data", "Pick a source, then Load.")
        elif self._view() is None:
            self.select_view("decay")
        else:
            self._schedule(immediate=True)

    # ══ selection ════════════════════════════════════════════════════════
    def profile(self) -> str:
        return self._profile_combo.currentData() or _WHOLE

    def select_profile(self, name: str) -> None:
        i = self._profile_combo.findData(name)
        if i >= 0:
            self._profile_combo.setCurrentIndex(i)

    def _profile_indices(self) -> list[int]:
        if not self.has_data:
            return []
        p = self.profile()
        if p:
            return list(self._data.profiles.get(p, []))
        return list(range(self._data.n))

    def _on_profile(self, _i: int) -> None:
        self._snd_list.blockSignals(True)
        self._snd_list.clear()
        idx = self._profile_indices()
        default = set(st.selection_every(idx, max(1, len(idx) // 4)))
        for i in idx:
            s = self._data.soundings[i]
            it = QListWidgetItem(getattr(s, "station_name", "") or f"S{i}")
            it.setData(_KEY, i)
            it.setFlags(it.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(Qt.CheckState.Checked if i in default and
                             len(default) <= 6 else Qt.CheckState.Unchecked)
            self._snd_list.addItem(it)
        if idx and not self.selected():
            self._snd_list.item(0).setCheckState(Qt.CheckState.Checked)
        self._snd_list.blockSignals(False)
        self._on_checks()

    def selected(self) -> list[int]:
        """Indices of the checked soundings."""
        out = []
        for r in range(self._snd_list.count()):
            it = self._snd_list.item(r)
            if it.checkState() == Qt.CheckState.Checked:
                out.append(it.data(_KEY))
        return out

    def _set_checked(self, keep) -> None:
        self._snd_list.blockSignals(True)
        for r in range(self._snd_list.count()):
            it = self._snd_list.item(r)
            it.setCheckState(Qt.CheckState.Checked if keep(it.data(_KEY))
                             else Qt.CheckState.Unchecked)
        self._snd_list.blockSignals(False)
        self._on_checks()

    def _check_all(self, on: bool) -> None:
        self._set_checked(lambda _i: on)

    def check_every(self, step: int) -> None:
        chosen = set(st.selection_every(self._profile_indices(), step))
        self._set_checked(lambda i: i in chosen)

    def _on_checks(self, *_a) -> None:
        n = len(self.selected())
        self._sel_status.setText(f"{n} of {self._snd_list.count()} "
                                 "soundings checked")
        self._update_actions()
        v = self._view()
        if v is not None and (v.needs == "soundings" or v.key == "dashboard"):
            self._schedule()

    # ══ views ════════════════════════════════════════════════════════════
    def _view(self) -> st.TDEMView | None:
        it = self._view_list.currentItem()
        key = it.data(_KEY) if it is not None else None
        return next((v for v in st.VIEWS if v.key == key), None)

    def select_view(self, key: str) -> None:
        for r in range(self._view_list.count()):
            if self._view_list.item(r).data(_KEY) == key:
                self._view_list.setCurrentRow(r)
                return

    def _on_view(self, _cur, _prev) -> None:
        v = self._view()
        if v is None:
            return
        self._view_desc.setText(v.help)
        _show_page(self._form_stack, self._forms.get(v.key, self._no_opts))
        self._schedule(immediate=True)

    def _schedule(self, *_a, immediate: bool = False) -> None:
        if not self._chk_auto.isChecked() or not self.has_data:
            return
        if immediate:
            self._on_draw()
        else:
            self._timer.start()

    def _on_draw(self) -> None:
        v = self._view()
        if v is None or not self.has_data:
            return
        form = self._forms.get(v.key)
        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            fig = st.render(v.key, self._data, selected=self.selected(),
                            profile=self.profile(),
                            values=form.values() if form else {})
        except Exception as exc:
            self._canvas_view.show_unavailable(v.label, f"{exc}")
            self._status.setText(f"✕ {v.label}: {exc}")
            return
        finally:
            self.unsetCursor()
        why = figure_blank_reason(fig)
        if why is not None:
            import matplotlib.pyplot as plt

            plt.close(fig)
            self._canvas_view.show_unavailable(v.label, why)
            self._status.setText(f"{v.label}: nothing to draw.")
            return
        self._canvas.show_figure(fig)
        self._canvas_view.show_canvas()
        self._tabs.setCurrentIndex(0)
        self._status.setText(f"✓ {v.label}")

    # ══ conversion ═══════════════════════════════════════════════════════
    def _update_actions(self) -> None:
        self._btn_convert.setEnabled(self.has_data and bool(self.selected()))
        self._btn_convert.setText(
            f"Convert {len(self.selected())} checked sounding"
            f"{'s' if len(self.selected()) != 1 else ''}")
        ok = self._converted is not None
        self._btn_save.setEnabled(ok)
        self._btn_send.setEnabled(ok)

    def convert(self):
        """Convert the checked soundings; returns the Sites."""
        snds = [self._data.soundings[i] for i in self.selected()]
        self._converted = st.convert(snds, self._conv_form.values())
        self._fill_converted_table()
        self._update_actions()
        return self._converted

    def _on_convert(self) -> None:
        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            sites = self.convert()
        except Exception as exc:
            self._status.setText(f"✕ Conversion failed: {exc}")
            return
        finally:
            self.unsetCursor()
        self._tabs.setCurrentIndex(2)
        self._status.setText(f"✓ Converted {len(sites)} sounding(s) to "
                             "impedance — Save or Send to survey.")

    def save_edis(self, folder: str) -> list:
        paths = self._converted.write(folder, exist_ok=True)
        self._status.setText(f"Saved {len(paths)} EDI file(s) to {folder}")
        return paths

    def _on_save(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "Save EDIs to")
        if folder:
            try:
                self.save_edis(folder)
            except Exception as exc:
                self._status.setText(f"✕ Save failed: {exc}")

    def send(self, mode: str = "append") -> None:
        if self._converted is None:
            return
        self.send_to_survey.emit(self._converted, mode)
        self._status.setText(
            f"Sent {len(self._converted)} site(s) to the survey ({mode}).")

    # ══ tables ═══════════════════════════════════════════════════════════
    @staticmethod
    def _fill(table: QTableWidget, cols, rows) -> None:
        table.setColumnCount(len(cols))
        table.setHorizontalHeaderLabels(cols)
        table.setRowCount(len(rows))
        for r, vals in enumerate(rows):
            for c, v in enumerate(vals):
                text = f"{v:.6g}" if isinstance(v, float) else str(v)
                table.setItem(r, c, QTableWidgetItem(text))
        table.resizeColumnsToContents()

    def _fill_soundings_table(self) -> None:
        rows = []
        if self.has_data:
            for i, s in enumerate(self._data.soundings[:5000]):
                rows.append([
                    getattr(s, "station_name", "") or f"S{i}",
                    st._profile_of(getattr(s, "station_name", "")),
                    float(getattr(s, "x", float("nan"))),
                    float(getattr(s, "y", float("nan"))),
                    float(getattr(s, "elevation", float("nan"))),
                    int(getattr(s, "n_gates", 0))])
        self._fill(self._snd_table, ["Sounding", "Profile", "X", "Y",
                                     "Elevation", "Gates"], rows)

    def _fill_converted_table(self) -> None:
        rows = []
        for site in (self._converted or []):
            f = getattr(site, "freq", None)
            rows.append([site.name, "—" if f is None else len(f),
                         "—" if f is None or not len(f) else float(min(f)),
                         "—" if f is None or not len(f) else float(max(f))])
        self._fill(self._conv_table, ["Site", "Frequencies", "f min (Hz)",
                                      "f max (Hz)"], rows)

    # ══ theme ════════════════════════════════════════════════════════════
    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)


__all__ = ["TDEMWindow"]
