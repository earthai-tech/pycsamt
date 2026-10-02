# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
ResPhaseSectionPanel — MTPy-style apparent-resistivity/phase pseudo-section.

Thin desktop wrapper around
:class:`pycsamt.emtools.resphase_psection.PlotResPhasePseudoSection`: one
column per impedance component, resistivity stacked above phase, a shared
period axis and shared colour scales, and the pyCSAMT station-marker
convention on the station axis. This is a richer, MTPy-style alternative to
the ``Pseudosection ρₐ`` / ``Pseudosection φ`` tabs (which draw one quantity
at a time via ``pycsamt.emtools.pseudosection``); both stay available side
by side rather than replacing one with the other.

Wired as the "Res/Phase Section" tab inside ``ProfilePanel``.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.plot_controller import style_axes
from pycsamt.app.desktop.widgets.canvas_stack import CanvasResultView

_STATION_SIDES = ["top", "bottom", "none"]
_PHASE_RANGES = ["0-90", "-45-45", "-90-90", "-180-180", "auto"]


class ResPhaseSectionPanel(QWidget):
    """Self-contained control strip + canvas for the res/phase pseudo-section."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._sites = None
        self._dark: bool = True
        # Redraw is deferred while this tab is not the one currently shown
        # by the owning ProfilePanel (see set_active()) -- the pseudo-
        # section rebuild is not free, and eagerly redrawing on every
        # set_sites() call would repeat that cost even for a survey the
        # user never opens this tab for.
        self._active = False
        self._dirty = False
        self._build_ui()

    # ── Construction ──────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 4)
        bar.setSpacing(8)

        bar.addWidget(QLabel("Stations:"))
        self._side_combo = QComboBox()
        self._side_combo.addItems(_STATION_SIDES)
        bar.addWidget(self._side_combo)

        bar.addWidget(QLabel("Phase range:"))
        self._phase_range_combo = QComboBox()
        self._phase_range_combo.addItems(_PHASE_RANGES)
        bar.addWidget(self._phase_range_combo)

        self._share_period_cb = QCheckBox("Share period")
        self._share_period_cb.setChecked(True)
        self._share_period_cb.setToolTip(
            "Resample every group onto one shared period grid. Turn off"
            " when stacking lines from different instruments so each"
            " keeps its own period window."
        )
        bar.addWidget(self._share_period_cb)

        self._panel_labels_cb = QCheckBox("Panel labels")
        bar.addWidget(self._panel_labels_cb)

        self._grid_cb = QCheckBox("Grid")
        bar.addWidget(self._grid_cb)

        bar.addStretch()

        self._btn_redraw = QPushButton("↻ Redraw")
        self._btn_redraw.clicked.connect(self._redraw)
        bar.addWidget(self._btn_redraw)

        self._btn_export = QPushButton("⬆ Export…")
        self._btn_export.setEnabled(False)
        self._btn_export.clicked.connect(self._on_export)
        bar.addWidget(self._btn_export)

        root.addLayout(bar)

        self._canvas_view = CanvasResultView(
            self, toolbar=True,
            empty_title="No res/phase section yet",
            empty_reason="Load survey data to view the res/phase pseudo-section.",
        )
        self._canvas = self._canvas_view.canvas
        self._canvas.set_refresh_callback(
            self._redraw, tooltip="Redraw the res/phase pseudo-section"
        )
        root.addWidget(self._canvas_view)

        for combo in (self._side_combo, self._phase_range_combo):
            combo.currentIndexChanged.connect(self._redraw)
        for cb in (self._share_period_cb, self._panel_labels_cb, self._grid_cb):
            cb.toggled.connect(self._redraw)

        self._draw_empty()

    # ── Public API ─────────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        self._sites = sites
        self._btn_export.setEnabled(sites is not None)
        self._request_redraw()

    def set_dark_mode(self, dark: bool) -> None:
        self._dark = dark
        self._request_redraw()

    # ── Lazy redraw ──────────────────────────────────────────────────────
    #
    # Driven explicitly by the owning ProfilePanel (set_active(), called
    # from its tab-switch handler) rather than Qt widget-visibility/
    # showEvent: a QTabWidget page's real on-screen visibility only
    # resolves once every ancestor up to the top-level window has been
    # shown, which makes showEvent unreliable before the floating window's
    # first .show() and in tests.

    def set_active(self, active: bool) -> None:
        """Called by ProfilePanel when this tab becomes/stops being current."""
        self._active = active
        if active and self._dirty:
            self._dirty = False
            self._redraw()

    def _request_redraw(self) -> None:
        """Redraw now if this is the active tab, else defer until it is."""
        if self._active:
            self._redraw()
        else:
            self._dirty = True

    # ── Drawing ────────────────────────────────────────────────────────

    def _draw_empty(self) -> None:
        self._canvas_view.show_unavailable(
            "No res/phase section yet",
            "Load survey data to view the res/phase pseudo-section.",
        )

    def _redraw(self, *_args) -> None:
        if self._sites is None:
            self._draw_empty()
            return
        try:
            from pycsamt.emtools.resphase_psection import (
                PlotResPhasePseudoSection,
            )

            phase_range = self._phase_range_combo.currentText()
            fig = PlotResPhasePseudoSection(
                self._sites,
                station_side=self._side_combo.currentText(),
                share_period=self._share_period_cb.isChecked(),
                panel_labels=self._panel_labels_cb.isChecked(),
                grid=self._grid_cb.isChecked(),
                phase_range=phase_range,
            ).plot()
            self._style_full_figure(fig)
            self._canvas.show_figure(fig)
            self._canvas_view.show_canvas()
        except Exception as exc:
            self._canvas_view.show_unavailable("Res/phase section error", str(exc))

    def _style_full_figure(self, fig) -> None:
        """Apply dark/light styling to every axes the library figure owns."""
        for ax in fig.axes:
            style_axes(ax, self._dark)
        fig.patch.set_facecolor("#1e1e2e" if self._dark else "#ffffff")

    # ── Export ────────────────────────────────────────────────────────

    def _on_export(self) -> None:
        from pycsamt.app.desktop.dialogs.export_dlg import ExportDialog

        ExportDialog(figure=self._canvas.figure, parent=self).exec()
