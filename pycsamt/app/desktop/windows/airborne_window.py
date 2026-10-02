# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
AirborneWindow — Airborne EM studio (ZTEM, AFMAG, AirMt, MobileMT).

Left
    **Data** — load EMTF-XML (a folder or one file); the card shows the
    technology detected, stations, flight lines and frequency band.
    **Flight line** — restrict every view to one line.
    **Views** — only the views that can draw the loaded technology
    (profiles, sections, maps, motion, tables; see
    :mod:`pycsamt.app.desktop.controllers.airborne_studio`).
    **Options** — the view's own settings, with the data's real
    frequencies / stations as choices.
    **Survey geometry** — geomagnetic field and aircraft attitude, shared
    by the motion views (shown when the view needs it).
Right
    **Plot** · **Table** (diagnostic tables, CSV export) · **Stations**
    (what was loaded).

MobileMT is generic-adapter-only: only already-decoded EMTF-XML is read;
raw vendor MobileMT files are a permanent restriction, not a gap.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
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
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import airborne_studio as st
from pycsamt.app.desktop.controllers.airborne_controller import (
    AirborneController,
)
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.windows._base import PanelWindow, make_group
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm

_KEY = Qt.ItemDataRole.UserRole
_MOBILEMT_NOTE = ("MobileMT: only already-decoded EMTF-XML is read — raw "
                  "vendor MobileMT files are never read (a permanent "
                  "restriction).")


class AirborneWindow(PanelWindow):
    """Airborne EM studio: technology-aware views, options and tables."""

    def __init__(self, parent: QWidget | None = None) -> None:
        self._ctrl = AirborneController()
        self._ctrl.dark = False  # figures are for publication: white
        self._ctx: dict = {}
        self._forms: dict[str, SettingsForm] = {}
        self._tables: list = []
        super().__init__(title="Airborne EM", session_key="airborne_window",
                         params_width=300, icon_name="induction",
                         parent=parent)
        self.resize(1200, 780)
        self._refresh_data()

    # ══ left panel ═══════════════════════════════════════════════════════
    def _build_params(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Data")
        self._data_status = QLabel("No airborne data loaded")
        self._data_status.setObjectName("InfoLabel")
        self._data_status.setWordWrap(True)
        self._data_status.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(self._data_status)
        row = QHBoxLayout()
        btn = QToolButton()
        btn.setText("Load EMTF-XML  ▾")
        btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        btn.setToolTip("EMTF-XML files, one per station. " + _MOBILEMT_NOTE)
        menu = QMenu(btn)
        menu.addAction("Folder…", self._on_load)
        menu.addAction("Single file…", self._on_load_file)
        btn.setMenu(menu)
        self._btn_load = btn
        clear = QPushButton("Clear")
        clear.clicked.connect(self._on_clear)
        row.addWidget(btn, 1)
        row.addWidget(clear)
        lay.addLayout(row)
        row = QHBoxLayout()
        row.addWidget(QLabel("Flight line:"))
        self._line_combo = QComboBox()
        self._line_combo.setToolTip("Restrict every view to one line")
        self._line_combo.currentIndexChanged.connect(
            lambda _i: self._schedule())
        row.addWidget(self._line_combo, 1)
        lay.addLayout(row)
        layout.addWidget(grp)

        grp, lay = make_group("Views")
        self._view_list = QListWidget()
        self._view_list.setMinimumHeight(230)
        self._view_list.currentItemChanged.connect(self._on_view)
        lay.addWidget(self._view_list)
        self._view_desc = QLabel("")
        self._view_desc.setObjectName("InfoLabel")
        self._view_desc.setWordWrap(True)
        lay.addWidget(self._view_desc)
        layout.addWidget(grp)

        self._grp_opts, lay = make_group("Options")
        self._form_stack = QStackedWidget()
        self._no_opts = QLabel("This view has no options.")
        self._no_opts.setObjectName("InfoLabel")
        self._form_stack.addWidget(self._no_opts)
        lay.addWidget(self._form_stack)
        layout.addWidget(self._grp_opts)

        self._grp_geo, lay = make_group("Survey geometry")
        self._geo_form = SettingsForm(list(st.GEOMETRY_FIELDS))
        self._geo_form.changed.connect(self._schedule)
        lay.addWidget(self._geo_form)
        layout.addWidget(self._grp_geo)

        row = QHBoxLayout()
        self._btn_draw = QPushButton("▶  Draw")
        self._btn_draw.setToolTip("Render the selected view (F5)")
        self._btn_draw.clicked.connect(self._on_draw)
        row.addWidget(self._btn_draw, 1)
        self._chk_auto = QCheckBox("Auto")
        self._chk_auto.setChecked(True)
        self._chk_auto.setToolTip("Redraw when the view or an option "
                                  "changes")
        row.addWidget(self._chk_auto)
        layout.addLayout(row)
        self._status = QLabel("")
        self._status.setObjectName("InfoLabel")
        self._status.setWordWrap(True)
        layout.addWidget(self._status)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self._on_draw)

    # ══ right panel ══════════════════════════════════════════════════════
    def _build_content(self, layout: QVBoxLayout) -> None:
        self._tabs = QTabWidget()
        self._canvas_view = CanvasResultView(
            toolbar=True, empty_title="No airborne data",
            empty_reason="Load EMTF-XML (a folder or one file).",
            empty_guidance=_MOBILEMT_NOTE)
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(self._on_draw,
                                          tooltip="Render the selected view")
        self._tabs.addTab(self._canvas_view, "Plot")

        tab = QWidget()
        v = QVBoxLayout(tab)
        row = QHBoxLayout()
        self._table_pick = QComboBox()
        self._table_pick.currentIndexChanged.connect(self._show_table)
        row.addWidget(self._table_pick, 1)
        b = QPushButton("Export CSV…")
        b.clicked.connect(self._export_table)
        row.addWidget(b)
        v.addLayout(row)
        self._table = QTableWidget()
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setAlternatingRowColors(True)
        v.addWidget(self._table, 1)
        self._tabs.addTab(tab, "Table")

        self._stations = QTableWidget()
        self._stations.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers)
        self._stations.setAlternatingRowColors(True)
        self._tabs.addTab(self._stations, "Stations")
        layout.addWidget(self._tabs)

    # ══ data ═════════════════════════════════════════════════════════════
    def load(self, path: str) -> int:
        """Load an EMTF-XML folder or file; returns the station count."""
        n = self._ctrl.load(path)
        self._refresh_data()
        return n

    def _on_load(self) -> None:
        path = QFileDialog.getExistingDirectory(self,
                                                "Select EMTF-XML folder")
        if path:
            self._load_reporting(path)

    def _on_load_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select EMTF-XML file", "",
            "EMTF-XML (*.xml);;All files (*)")
        if path:
            self._load_reporting(path)

    def _load_reporting(self, path: str) -> None:
        try:
            self.load(path)
        except Exception as exc:
            self._data_status.setText(f"Load failed: {exc}")

    def _on_clear(self) -> None:
        self._ctrl.clear()
        self._refresh_data()

    def _asites(self):
        s = self._ctrl.state.asites
        line = self._line_combo.currentData()
        if s is not None and line:
            s = s.select(predicate=lambda x: x.line_id == line)
        return s

    def _refresh_data(self) -> None:
        s = self._ctrl.state.asites
        has = self._ctrl.has_data
        self._ctx = st.data_context(s) if has else {}
        # data card
        if not has:
            self._data_status.setText("No airborne data loaded")
        else:
            techs = ", ".join(st.TECH_LABELS.get(t, t)
                              for t in self._ctx["techs"]) or "unknown"
            f = self._ctx["freqs"]
            band = (f"{min(f):.4g} – {max(f):.4g} Hz ({len(f)})"
                    if f else "no frequencies")
            lines = len(self._ctx["lines"])
            self._data_status.setText(
                f"<b>{techs}</b><br>{len(s)} stations · "
                f"{lines or 'no'} flight line{'s' if lines != 1 else ''}"
                f"<br>{band}")
        # flight lines
        self._line_combo.blockSignals(True)
        self._line_combo.clear()
        self._line_combo.addItem("All lines", "")
        for ln in self._ctx.get("lines", []):
            self._line_combo.addItem(ln, ln)
        self._line_combo.blockSignals(False)
        self._line_combo.setEnabled(len(self._ctx.get("lines", [])) > 1)
        # views + their forms (choices depend on the data)
        for form in self._forms.values():
            self._form_stack.removeWidget(form)
            form.deleteLater()
        self._forms.clear()
        self._view_list.blockSignals(True)
        self._view_list.clear()
        views = st.views_for(self._ctx.get("techs", []))
        for group in st.GROUPS:
            members = [v for v in views if v.group == group]
            if not members:
                continue
            head = QListWidgetItem(group.upper())
            head.setFlags(Qt.ItemFlag.NoItemFlags)
            self._view_list.addItem(head)
            for v in members:
                it = QListWidgetItem("   " + v.label)
                it.setData(_KEY, v.key)
                it.setToolTip(v.help)
                self._view_list.addItem(it)
        self._view_list.blockSignals(False)
        self._fill_stations()
        if views:
            self.select_view(views[0].key)
        else:
            self._grp_geo.setVisible(False)
            self._canvas_view.show_unavailable(
                "No airborne data" if not has else "No view for this data",
                "Load EMTF-XML (a folder or one file)." if not has else
                "The technology of these files is not recognised.",
                _MOBILEMT_NOTE)

    def _fill_stations(self) -> None:
        s = self._ctrl.state.asites
        cols = ["Station", "Line", "Technology", "Latitude", "Longitude",
                "Elevation", "Frequencies"]
        self._stations.setColumnCount(len(cols))
        self._stations.setHorizontalHeaderLabels(cols)
        rows = list(s) if self._ctrl.has_data else []
        self._stations.setRowCount(len(rows))
        for r, site in enumerate(rows):
            try:
                lat, lon, elev = site.coords
            except Exception:
                lat = lon = elev = float("nan")
            f = getattr(site, "freq", None)
            vals = [site.name, site.line_id or "—",
                    st.TECH_LABELS.get(site.technology, site.technology
                                       or "—"),
                    f"{lat:.6f}", f"{lon:.6f}", f"{elev:.1f}",
                    "—" if f is None else str(len(f))]
            for c, text in enumerate(vals):
                self._stations.setItem(r, c, QTableWidgetItem(str(text)))
        self._stations.resizeColumnsToContents()

    # ══ views ════════════════════════════════════════════════════════════
    def _view(self) -> st.AirView | None:
        it = self._view_list.currentItem()
        key = it.data(_KEY) if it is not None else None
        return next((v for v in st.VIEWS if v.key == key), None)

    def select_view(self, key: str) -> None:
        for r in range(self._view_list.count()):
            if self._view_list.item(r).data(_KEY) == key:
                self._view_list.setCurrentRow(r)
                return

    def _on_view(self, cur, _prev) -> None:
        v = self._view()
        if v is None:
            return
        self._view_desc.setText(v.help)
        if v.key not in self._forms:
            fields = st.options_for(v, self._ctx)
            if fields:
                form = SettingsForm(fields)
                form.changed.connect(self._schedule)
                self._forms[v.key] = form
                self._form_stack.addWidget(form)
        form = self._forms.get(v.key)
        self._form_stack.setCurrentWidget(form or self._no_opts)
        self._grp_geo.setVisible(st.needs_geometry(v))
        self._schedule(immediate=True)

    def _schedule(self, *_a, immediate: bool = False) -> None:
        if not self._chk_auto.isChecked() or not self._ctrl.has_data:
            return
        if immediate:
            self._on_draw()
        else:
            self._timer.start()

    def _on_draw(self) -> None:
        v = self._view()
        if v is None or not self._ctrl.has_data:
            return
        form = self._forms.get(v.key)
        values = form.values() if form else {}
        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            out = st.run_view(v, self._asites(), values,
                              self._geo_form.values())
        except Exception as exc:
            self._canvas_view.show_unavailable(v.label, f"{exc}")
            self._status.setText(f"✕ {v.label}: {exc}")
            return
        finally:
            self.unsetCursor()
        if v.kind == "table":
            self._set_tables(out)
            self._tabs.setCurrentIndex(1)
            self._status.setText(f"{v.label}: {len(self._tables)} "
                                 "table(s).")
            return
        why = figure_blank_reason(out)
        if why is not None:
            import matplotlib.pyplot as plt

            plt.close(out)
            self._canvas_view.show_unavailable(v.label, why)
            self._status.setText(f"{v.label}: nothing to draw.")
            return
        self._canvas.show_figure(out)
        self._canvas_view.show_canvas()
        self._tabs.setCurrentIndex(0)
        self._status.setText(f"✓ {v.label}")

    # ══ tables ═══════════════════════════════════════════════════════════
    def _set_tables(self, tables) -> None:
        self._tables = list(tables)
        self._table_pick.blockSignals(True)
        self._table_pick.clear()
        for title, df in self._tables:
            self._table_pick.addItem(f"{title}  ({len(df)} rows)")
        self._table_pick.blockSignals(False)
        self._show_table(0)

    def _show_table(self, i: int) -> None:
        if not 0 <= i < len(self._tables):
            self._table.setRowCount(0)
            return
        _title, df = self._tables[i]
        self._table.setColumnCount(len(df.columns))
        self._table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        self._table.setRowCount(len(df))
        for r in range(len(df)):
            for c, col in enumerate(df.columns):
                v = df.iloc[r, c]
                text = f"{v:.6g}" if isinstance(v, float) else str(v)
                self._table.setItem(r, c, QTableWidgetItem(text))
        self._table.resizeColumnsToContents()

    def export_table(self, path: str, index: int | None = None) -> str:
        i = self._table_pick.currentIndex() if index is None else index
        _title, df = self._tables[i]
        df.to_csv(path, index=False)
        return path

    def _export_table(self) -> None:
        if not self._tables:
            self._status.setText("Draw a table view first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export table", "",
                                              "CSV (*.csv)")
        if path:
            self.export_table(path)
            self._status.setText(f"Exported {path}")

    # ══ theme ════════════════════════════════════════════════════════════
    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = False


__all__ = ["AirborneWindow"]
