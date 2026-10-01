# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
QCDashboardWindow — transfer-function QC studio.

Left
    **Survey** — stations, lines and the pass / warn / fail count.
    **Scope** — one line or all of them (profiles and pseudo-sections then
    draw one panel per line, or all lines together), and the station for
    single-station views.
    **Diagnostics** — Confidence · Coverage · Noise / SNR · Dimensionality
    & skew · Static shift · Distortion & source · Strike (filterable).
    **Analysis** and **Plot view** — the diagnostic's own parameters as
    drop-downs, spin boxes and fields (colour maps from the shared
    catalogue, stations from the survey).
Right
    **Plot** · **Summary** (one row per station, status from the
    confidence ratio and the QC flags; sortable, CSV export; double-click a
    station to open its dashboard) · **Quick-look** (the multi-panel
    survey overview).

See :mod:`pycsamt.app.desktop.controllers.qc_studio`.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers import qc_studio as st
from pycsamt.app.desktop.controllers.correction_views import (
    figure_blank_reason,
)
from pycsamt.app.desktop.controllers.qc_controller import (
    QC_STATIC_SHIFT_PLOTS,
    QCController,
)
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas
from pycsamt.app.desktop.windows._base import PanelWindow, make_group
from pycsamt.app.desktop.windows.inversion.forms import SettingsForm

_KEY = Qt.ItemDataRole.UserRole
_STATUS_COLORS = {"pass": "#2e9d57", "warn": "#d9922e", "fail": "#d62728"}


class _NumItem(QTableWidgetItem):
    """Sorts by the number it shows."""

    def __lt__(self, other):  # noqa: D105
        try:
            return float(self.data(_KEY)) < float(other.data(_KEY))
        except (TypeError, ValueError):
            return super().__lt__(other)


class QCDashboardWindow(PanelWindow):
    """QC studio: scoped diagnostics, full parameters, station summary."""

    def __init__(self, parent: QWidget | None = None) -> None:
        self._ctrl = QCController()
        self._ctrl.dark = False  # figures are for publication: white
        self._sites = None
        self._lines: dict[str, str] = {}
        self._forms: dict[tuple, tuple[SettingsForm, SettingsForm]] = {}
        self._method: dict[str, str] = {}
        self._summary = None
        self._summary_stale = True
        self._quicklook_stale = True
        super().__init__(title="QC Studio", session_key="qc_dashboard",
                         params_width=310, icon_name="qc", parent=parent)
        self.resize(1240, 820)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(350)
        self._timer.timeout.connect(self._on_draw)
        self._refresh_scope()
        self.select_view("plot_confidence_profile")

    # ══ left panel ═══════════════════════════════════════════════════════
    def _build_params(self, layout: QVBoxLayout) -> None:
        grp, lay = make_group("Survey")
        self._data_status = QLabel("")
        self._data_status.setObjectName("InfoLabel")
        self._data_status.setWordWrap(True)
        self._data_status.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(self._data_status)
        layout.addWidget(grp)

        grp, lay = make_group("Scope")
        form = QFormLayout()
        form.setSpacing(5)
        self._line_combo = QComboBox()
        self._line_combo.setToolTip("Draw one line, or every line in scope")
        self._line_combo.currentIndexChanged.connect(self._on_line)
        form.addRow("Line:", self._line_combo)
        self._layout_combo = QComboBox()
        self._layout_combo.addItem("One panel per line", "panels")
        self._layout_combo.addItem("All lines together", "together")
        self._layout_combo.setToolTip(
            "Profiles and pseudo-sections assume one line: several lines "
            "are drawn one panel each, or together on one axes")
        self._layout_combo.currentIndexChanged.connect(
            lambda _i: self._schedule())
        self._layout_row = QLabel("Several lines:")
        form.addRow(self._layout_row, self._layout_combo)
        self._station_combo = QComboBox()
        self._station_combo.currentIndexChanged.connect(
            lambda _i: self._schedule())
        self._station_row = QLabel("Station:")
        form.addRow(self._station_row, self._station_combo)
        lay.addLayout(form)
        layout.addWidget(grp)

        grp, lay = make_group("Diagnostics")
        self._filter = QLineEdit()
        self._filter.setPlaceholderText("Filter diagnostics…")
        self._filter.setClearButtonEnabled(True)
        self._filter.textChanged.connect(self._apply_filter)
        lay.addWidget(self._filter)
        self._view_list = QListWidget()
        self._view_list.setMinimumHeight(260)
        for group in st.GROUPS:
            head = QListWidgetItem(group.upper())
            head.setFlags(Qt.ItemFlag.NoItemFlags)
            head.setData(_KEY, None)
            self._view_list.addItem(head)
            for v in (v for v in st.VIEWS if v.group == group):
                it = QListWidgetItem("   " + v.label)
                it.setData(_KEY, v.fn)
                it.setToolTip(v.description)
                self._view_list.addItem(it)
        self._view_list.currentItemChanged.connect(self._on_view)
        lay.addWidget(self._view_list)
        self._view_desc = QLabel("")
        self._view_desc.setObjectName("InfoLabel")
        self._view_desc.setWordWrap(True)
        lay.addWidget(self._view_desc)
        layout.addWidget(grp)

        self._grp_analysis, lay = make_group("Analysis")
        self._analysis_stack = QStackedWidget()
        lay.addWidget(self._analysis_stack)
        layout.addWidget(self._grp_analysis)
        self._grp_view, lay = make_group("Plot view")
        self._view_stack = QStackedWidget()
        lay.addWidget(self._view_stack)
        layout.addWidget(self._grp_view)
        self._empty_a = QLabel("Standard analytical defaults.")
        self._empty_a.setObjectName("InfoLabel")
        self._empty_v = QLabel("")
        self._analysis_stack.addWidget(self._empty_a)
        self._view_stack.addWidget(self._empty_v)

        row = QHBoxLayout()
        self._btn_draw = QPushButton("▶  Draw")
        self._btn_draw.setToolTip("Render the selected diagnostic")
        self._btn_draw.clicked.connect(self._on_draw)
        row.addWidget(self._btn_draw, 1)
        self._chk_auto = QCheckBox("Auto")
        self._chk_auto.setChecked(True)
        self._chk_auto.setToolTip("Redraw when the diagnostic, scope or an "
                                  "option changes")
        row.addWidget(self._chk_auto)
        layout.addLayout(row)
        self._status = QLabel("")
        self._status.setObjectName("InfoLabel")
        self._status.setWordWrap(True)
        layout.addWidget(self._status)

    # ══ right panel ══════════════════════════════════════════════════════
    def _build_content(self, layout: QVBoxLayout) -> None:
        self._tabs = QTabWidget()
        self._canvas_view = CanvasResultView(
            toolbar=True, empty_title="No survey loaded",
            empty_reason="Load EDI or EMTF-XML data in the main window.")
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(self._on_draw,
                                          tooltip="Redraw this diagnostic")
        self._tabs.addTab(self._canvas_view, "Plot")

        tab = QWidget()
        v = QVBoxLayout(tab)
        row = QHBoxLayout()
        self._thr: dict[str, QDoubleSpinBox] = {}
        for key, label, val, lo, hi, step, tip in (
                ("ci_hi", "Safe ≥", 0.95, 0.0, 1.0, 0.01,
                 "Confidence ratio for a safe station"),
                ("ci_lo", "Recoverable ≥", 0.85, 0.0, 1.0, 0.01,
                 "Below this a station fails"),
                ("min_snr", "Min SNR", 2.0, 0.0, 1e4, 0.5,
                 "Median SNR below this warns"),
                ("max_skew", "Max |skew|", 6.0, 0.0, 90.0, 0.5,
                 "Median |β| (°) above this warns")):
            row.addWidget(QLabel(label))
            spin = QDoubleSpinBox()
            spin.setRange(lo, hi)
            spin.setSingleStep(step)
            spin.setDecimals(2)
            spin.setValue(val)
            spin.setToolTip(tip)
            spin.valueChanged.connect(self._invalidate_summary)
            self._thr[key] = spin
            row.addWidget(spin)
        row.addStretch(1)
        v.addLayout(row)
        row = QHBoxLayout()
        self._status_filter = QComboBox()
        for text, key in (("All stations", ""), ("Fail", "fail"),
                          ("Warn", "warn"), ("Pass", "pass")):
            self._status_filter.addItem(text, key)
        self._status_filter.currentIndexChanged.connect(self._fill_summary)
        row.addWidget(self._status_filter)
        self._summary_counts = QLabel("")
        self._summary_counts.setTextFormat(Qt.TextFormat.RichText)
        row.addWidget(self._summary_counts, 1)
        b = QPushButton("Export CSV…")
        b.clicked.connect(self._export_summary)
        row.addWidget(b)
        v.addLayout(row)
        self._summary_table = QTableWidget()
        self._summary_table.setEditTriggers(
            QTableWidget.EditTrigger.NoEditTriggers)
        self._summary_table.setSelectionBehavior(
            QTableWidget.SelectionBehavior.SelectRows)
        self._summary_table.setAlternatingRowColors(True)
        self._summary_table.setSortingEnabled(True)
        self._summary_table.cellDoubleClicked.connect(self._on_summary_pick)
        v.addWidget(self._summary_table, 1)
        hint = QLabel("Double-click a station to open its confidence "
                      "dashboard.")
        hint.setObjectName("InfoLabel")
        v.addWidget(hint)
        self._tabs.addTab(tab, "Summary")

        self._quick = MplCanvas(toolbar=True)
        self._tabs.addTab(self._quick, "Quick-look")
        self._tabs.currentChanged.connect(self._on_tab)
        layout.addWidget(self._tabs)

    # ══ data ═════════════════════════════════════════════════════════════
    def set_lines(self, lines: dict | None) -> None:
        """``{station: line}`` of the survey (from the station table)."""
        self._lines = {str(k): str(v) for k, v in (lines or {}).items()
                       if v is not None and str(v) not in ("", "nan")}
        if self._sites is not None:
            self._refresh_scope()

    def set_sites(self, sites) -> None:
        super().set_sites(sites)
        self._sites = sites
        self._ctrl.set_sites(sites)
        self._summary_stale = self._quicklook_stale = True
        self._refresh_scope()
        if self.isVisible():
            self._schedule(immediate=True)

    @property
    def has_data(self) -> bool:
        return self._sites is not None and bool(st.station_names(self._sites))

    def _effective_lines(self) -> dict:
        """Line labels restricted to the stations in scope (none when
        every station sits on one line)."""
        names = st.station_names(self._sites)
        lines = {n: self._lines[n] for n in names if n in self._lines}
        return lines if len(set(lines.values())) > 1 else {}

    def _refresh_scope(self) -> None:
        names = st.station_names(self._sites) if self._sites is not None \
            else []
        groups = st.line_groups_of(self._sites, self._effective_lines()) \
            if names else {}
        # survey card
        if not names:
            self._data_status.setText("No survey loaded — load data in the "
                                      "main window.")
        else:
            nl = len([g for g in groups if g])
            self._data_status.setText(
                f"<b>{len(names)} stations</b>"
                + (f" · {nl} lines" if nl > 1 else "")
                + "<br><small>Summary tab: pass / warn / fail per "
                  "station</small>")
        # line combo
        cur = self._line_combo.currentData()
        self._line_combo.blockSignals(True)
        self._line_combo.clear()
        many = len([g for g in groups if g]) > 1
        self._line_combo.addItem("All lines" if many else "Survey", "")
        if many:
            for g, members in groups.items():
                self._line_combo.addItem(f"{g}  ({len(members)})", g)
        i = self._line_combo.findData(cur)
        self._line_combo.setCurrentIndex(max(i, 0))
        self._line_combo.setEnabled(many)
        self._line_combo.blockSignals(False)
        self._fill_stations()
        self._update_scope_rows()
        if not names:
            self._canvas_view.show_unavailable(
                "No survey loaded",
                "The QC studio works on the survey loaded in the main "
                "window.", "Load EDI or EMTF-XML data, then pick a "
                           "diagnostic.")

    def _scope_line(self) -> str:
        return self._line_combo.currentData() or ""

    def _fill_stations(self) -> None:
        groups = st.line_groups_of(self._sites, self._effective_lines()) \
            if self._sites is not None else {}
        line = self._scope_line()
        names = groups.get(line, []) if line else [
            n for members in groups.values() for n in members]
        cur = self._station_combo.currentData()
        self._station_combo.blockSignals(True)
        self._station_combo.clear()
        for n in names:
            self._station_combo.addItem(n, n)
        i = self._station_combo.findData(cur)
        self._station_combo.setCurrentIndex(max(i, 0))
        self._station_combo.blockSignals(False)

    def _on_line(self, _i: int) -> None:
        self._fill_stations()
        self._update_scope_rows()
        self._schedule()

    def _update_scope_rows(self) -> None:
        v = self._view()
        several = (not self._scope_line()
                   and self._line_combo.count() > 1)
        show_layout = bool(v and v.per_line and several)
        self._layout_row.setVisible(show_layout)
        self._layout_combo.setVisible(show_layout)
        show_station = bool(v and v.station)
        self._station_row.setVisible(show_station)
        self._station_combo.setVisible(show_station)

    # ══ views ════════════════════════════════════════════════════════════
    def _view(self) -> st.QCView | None:
        it = self._view_list.currentItem()
        key = it.data(_KEY) if it is not None else None
        return next((v for v in st.VIEWS if v.fn == key), None)

    def select_view(self, fn: str) -> None:
        for r in range(self._view_list.count()):
            if self._view_list.item(r).data(_KEY) == fn:
                self._view_list.setCurrentRow(r)
                return

    def _apply_filter(self, text: str) -> None:
        text = text.strip().lower()
        head = None
        shown = False
        for r in range(self._view_list.count()):
            it = self._view_list.item(r)
            if it.data(_KEY) is None:
                if head is not None:
                    head.setHidden(not shown)
                head, shown = it, False
                continue
            hit = (not text or text in it.text().lower()
                   or text in (it.toolTip() or "").lower())
            it.setHidden(not hit)
            shown = shown or hit
        if head is not None:
            head.setHidden(not shown)

    def _forms_for(self, v: st.QCView) -> tuple[SettingsForm, SettingsForm]:
        method = self._method.get(v.fn) if v.fn in QC_STATIC_SHIFT_PLOTS \
            else None
        key = (v.fn, method)
        if key not in self._forms:
            fields = st.options_for(v, method=method)
            fa = SettingsForm([f for f in fields if f.section == "Analysis"])
            fv = SettingsForm([f for f in fields if f.section == "View"])
            for form in (fa, fv):
                form.changed.connect(self._on_option)
            self._analysis_stack.addWidget(fa)
            self._view_stack.addWidget(fv)
            self._forms[key] = (fa, fv)
        return self._forms[key]

    def _current_forms(self):
        v = self._view()
        return self._forms_for(v) if v is not None else (None, None)

    def _on_view(self, _cur, _prev) -> None:
        v = self._view()
        if v is None:
            return
        self._view_desc.setText(v.description)
        self._show_forms(v)
        self._update_scope_rows()
        self._schedule(immediate=True)

    def _show_forms(self, v: st.QCView) -> None:
        fa, fv = self._forms_for(v)
        self._analysis_stack.setCurrentWidget(fa if fa.fields
                                              else self._empty_a)
        self._view_stack.setCurrentWidget(fv)
        self._grp_view.setVisible(bool(fv.fields))

    def _on_option(self) -> None:
        v = self._view()
        if v is not None and v.fn in QC_STATIC_SHIFT_PLOTS:
            fa, _fv = self._forms_for(v)
            method = fa.values().get("method")
            if method and method != self._method.get(v.fn, "ama"):
                self._method[v.fn] = method  # method-specific parameters
                self._forms_for(v)[0].set_values({"method": method})
                self._show_forms(v)
        self._schedule()

    def _schedule(self, *_a, immediate: bool = False) -> None:
        if not self._chk_auto.isChecked() or not self.has_data:
            return
        if not self.isVisible():
            return
        if immediate:
            self._timer.stop()
            self._on_draw()
        else:
            self._timer.start()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if self.has_data and self._chk_auto.isChecked():
            QTimer.singleShot(0, self._on_draw)

    def _on_run(self) -> None:  # the main window's refresh hook
        self._on_draw()

    def _on_draw(self) -> None:
        v = self._view()
        if v is None or not self.has_data:
            return
        fa, fv = self._forms_for(v)
        values = {**fa.values(), **fv.values()}
        kwargs, err = st.to_kwargs(v, values, self._method.get(v.fn))
        if err:
            self._canvas_view.show_unavailable(
                "Adjust the options", err,
                "Correct the value in Analysis or Plot view; the plot "
                "updates on its own.")
            self._status.setText("Waiting for valid options.")
            return
        self._status.setText(f"Drawing {v.label}…")
        self.setCursor(Qt.CursorShape.WaitCursor)
        self.repaint()
        try:
            fig, unavailable = st.render(
                v, self._sites, kwargs, lines=self._effective_lines(),
                line=self._scope_line(),
                layout=self._layout_combo.currentData(),
                station=self._station_combo.currentData() or "",
                ctrl=self._ctrl)
        except Exception as exc:
            self._canvas_view.show_unavailable(v.label, str(exc))
            self._status.setText(f"✕ {v.label}: {exc}")
            return
        finally:
            self.unsetCursor()
        import matplotlib.pyplot as plt

        if unavailable is not None:
            plt.close(fig)
            self._canvas_view.show_unavailable(
                unavailable.title, unavailable.reason, unavailable.guidance)
            self._status.setText(f"{v.label}: unavailable.")
            return
        why = figure_blank_reason(fig)
        if why is not None:
            plt.close(fig)
            self._canvas_view.show_unavailable(v.label, why)
            self._status.setText(f"{v.label}: nothing to draw.")
            return
        self._canvas.show_figure(fig)
        self._canvas_view.show_canvas()
        self._tabs.setCurrentIndex(0)
        self._status.setText(f"✓ {v.label}")

    # ══ summary ══════════════════════════════════════════════════════════
    def _invalidate_summary(self, *_a) -> None:
        self._summary_stale = True
        if self._tabs.currentIndex() == 1:
            self.refresh_summary()

    def _on_tab(self, i: int) -> None:
        if i == 1 and self._summary_stale:
            self.refresh_summary()
        elif i == 2 and self._quicklook_stale:
            self.refresh_quicklook()

    def refresh_summary(self):
        """Recompute the station summary; returns the DataFrame."""
        if not self.has_data:
            self._summary = None
            self._fill_summary()
            return None
        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            self._summary = st.scorecard(
                self._sites, self._effective_lines(),
                ci_hi=self._thr["ci_hi"].value(),
                ci_lo=self._thr["ci_lo"].value(),
                min_snr=self._thr["min_snr"].value(),
                max_skew=self._thr["max_skew"].value())
        except Exception as exc:
            self._summary = None
            self._summary_counts.setText(f"Summary failed: {exc}")
            return None
        finally:
            self.unsetCursor()
        self._summary_stale = False
        self._fill_summary()
        return self._summary

    def _fill_summary(self, *_a) -> None:
        t = self._summary_table
        t.setSortingEnabled(False)
        df = self._summary
        if df is None or df.empty:
            t.setRowCount(0)
            t.setColumnCount(0)
            self._summary_counts.setText("")
            t.setSortingEnabled(True)
            return
        counts = df["Status"].value_counts().to_dict()
        self._summary_counts.setText("  ".join(
            f"<span style='color:{_STATUS_COLORS[k]}'><b>{counts.get(k, 0)}"
            f"</b> {k}</span>" for k in ("pass", "warn", "fail")))
        want = self._status_filter.currentData()
        rows = df[df["Status"] == want] if want else df
        t.setColumnCount(len(df.columns))
        t.setHorizontalHeaderLabels(list(df.columns))
        t.setRowCount(len(rows))
        for r, (_i, rec) in enumerate(rows.iterrows()):
            for c, col in enumerate(df.columns):
                val = rec[col]
                if isinstance(val, float):
                    item = _NumItem("—" if val != val else f"{val:.3g}")
                    item.setData(_KEY, val)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                else:
                    item = QTableWidgetItem(str(val))
                if col == "Status":
                    item.setForeground(QColor(_STATUS_COLORS.get(val,
                                                                 "#888")))
                t.setItem(r, c, item)
        t.resizeColumnsToContents()
        t.setSortingEnabled(True)

    def _on_summary_pick(self, row: int, _col: int) -> None:
        cols = [self._summary_table.horizontalHeaderItem(c).text()
                for c in range(self._summary_table.columnCount())]
        station = self._summary_table.item(row, cols.index("Station")).text()
        self.open_station(station)

    def open_station(self, station: str) -> None:
        """Scope to *station*'s line and open its confidence dashboard."""
        line = self._effective_lines().get(station, "")
        i = self._line_combo.findData(line)
        if i >= 0:
            self._line_combo.setCurrentIndex(i)
        self.select_view("plot_station_confidence_dashboard")
        j = self._station_combo.findData(station)
        if j >= 0:
            self._station_combo.setCurrentIndex(j)
        self._on_draw()

    def export_summary(self, path: str) -> str:
        if self._summary is None:
            self.refresh_summary()
        self._summary.to_csv(path, index=False)
        return path

    def _export_summary(self) -> None:
        if not self.has_data:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export QC summary",
                                              "qc_summary.csv",
                                              "CSV (*.csv)")
        if path:
            self.export_summary(path)
            self._status.setText(f"Exported {path}")

    # ══ quick-look ═══════════════════════════════════════════════════════
    def refresh_quicklook(self) -> None:
        if not self.has_data:
            return
        import pycsamt.emtools as et

        self.setCursor(Qt.CursorShape.WaitCursor)
        try:
            self._quick.show_figure(et.plot_qc_quicklook(self._sites,
                                                         verbose=0))
            self._quicklook_stale = False
        except Exception as exc:
            self._status.setText(f"Quick-look failed: {exc}")
        finally:
            self.unsetCursor()

    # ══ theme ════════════════════════════════════════════════════════════
    def set_dark_mode(self, dark: bool) -> None:
        super().set_dark_mode(dark)
        self._ctrl.dark = False


__all__ = ["QCDashboardWindow"]
