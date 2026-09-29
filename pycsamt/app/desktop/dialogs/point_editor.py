# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PointEditorDialog — Edit ▸ Frequencies ▸ Point Editor.

One station at a time: apparent resistivity (top) and phase (bottom) of
the XY and YX impedances against frequency, with error bars.

Selecting
    Click a point to select it (Ctrl+click adds); drag a box on either
    panel to select every visible point inside; Esc clears.
Editing (buttons or keys)
    **M** mask · **D** delete the frequency row · **I** interpolate from
    neighbours · **R** restore the original value.  "All stations" applies
    the same action at the same frequencies on every station.
Static shift
    Shift + drag a resistivity curve up or down (or type a factor) scales
    that whole curve -- the standard static-shift correction; phase is
    unchanged.
Navigation
    ← / → previous / next station; Ctrl+Z undoes the last edit here.

Masked points are shown as grey crosses, interpolated ones as hollow
diamonds.  Nothing reaches the survey until **Apply**, which hands the
edited stations back as one undoable edit; every change is also written to
the station's EDI ``INFO`` block (see :mod:`~pycsamt.app.desktop.
controllers.point_edits`).
"""

from __future__ import annotations

import numpy as np
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QToolButton,
    QVBoxLayout,
)

from pycsamt.app.desktop.controllers import point_edits as pe
from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_COLOURS = {"xy": "#1f77b4", "yx": "#d62728"}


class PointEditorDialog(QDialog):
    """Mask, delete, interpolate, restore or static-shift data points."""

    def __init__(self, sites, *, station: str | None = None,
                 parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Point Editor")
        self.resize(1000, 700)
        self._orig = {str(s.name): s for s in sites}  # as opened
        self._order = list(self._orig)
        self._cur = {k: v for k, v in self._orig.items()}  # edited copies
        self._undo: list[tuple[str, dict]] = []  # (label, snapshot of _cur)
        self._interp: dict[str, set] = {}  # station -> {(row, comp)}
        self._sel: set[tuple[int, str]] = set()  # (row, comp)
        self._drag = None
        self._build_ui()
        start = station if station in self._orig else self._order[0]
        self._station.setCurrentText(start)
        self._draw()

    # ── UI ────────────────────────────────────────────────────────────
    def _build_ui(self) -> None:
        v = QVBoxLayout(self)
        top = QHBoxLayout()
        self._prev = QToolButton()
        self._prev.setText("◀")
        self._prev.setToolTip("Previous station (←)")
        self._prev.clicked.connect(lambda: self.step_station(-1))
        self._station = QComboBox()
        self._station.addItems(self._order)
        self._station.currentTextChanged.connect(self._on_station)
        self._next = QToolButton()
        self._next.setText("▶")
        self._next.setToolTip("Next station (→)")
        self._next.clicked.connect(lambda: self.step_station(+1))
        top.addWidget(QLabel("Station:"))
        top.addWidget(self._prev)
        top.addWidget(self._station, 1)
        top.addWidget(self._next)
        top.addSpacing(12)
        self._show = {}
        for comp in ("xy", "yx"):
            cb = QCheckBox(comp.upper())
            cb.setChecked(True)
            cb.setStyleSheet(f"QCheckBox {{ color: {_COLOURS[comp]}; "
                             "font-weight: 600; }")
            cb.toggled.connect(lambda _on: self._draw())
            self._show[comp] = cb
            top.addWidget(cb)
        self._errbars = QCheckBox("Error bars")
        self._errbars.setChecked(True)
        self._errbars.toggled.connect(lambda _on: self._draw())
        top.addWidget(self._errbars)
        v.addLayout(top)

        self._canvas = MplCanvas(self, toolbar=True)
        v.addWidget(self._canvas, 1)
        c = self._canvas._canvas
        c.mpl_connect("button_press_event", self._on_press)
        c.mpl_connect("motion_notify_event", self._on_motion)
        c.mpl_connect("button_release_event", self._on_release)

        act = QHBoxLayout()
        self._sel_lbl = QLabel("Nothing selected")
        self._sel_lbl.setObjectName("InfoLabel")
        act.addWidget(self._sel_lbl, 1)
        self._all_stations = QCheckBox("All stations")
        self._all_stations.setToolTip("Apply mask / delete / interpolate / "
                                      "restore at the same frequencies on "
                                      "every station")
        act.addWidget(self._all_stations)
        for text, key, tip, slot in (
                ("Mask", "M", "Hide the selected values (kept as NaN)",
                 lambda: self.edit("mask")),
                ("Delete frequency", "D", "Remove the whole frequency rows",
                 lambda: self.edit("delete")),
                ("Interpolate", "I", "Replace with values interpolated "
                 "from the valid neighbours", lambda: self.edit("interp")),
                ("Restore", "R", "Put back the original values",
                 lambda: self.edit("restore"))):
            b = QPushButton(f"{text}  ({key})")
            b.setToolTip(tip)
            b.clicked.connect(slot)
            act.addWidget(b)
            QShortcut(QKeySequence(key), self, activated=slot)
        v.addLayout(act)

        shift = QHBoxLayout()
        shift.addWidget(QLabel("Static shift:"))
        self._shift_comp = QComboBox()
        self._shift_comp.addItems(["XY", "YX"])
        self._shift_factor = QDoubleSpinBox()
        self._shift_factor.setRange(0.001, 1000.0)
        self._shift_factor.setDecimals(3)
        self._shift_factor.setValue(1.0)
        self._shift_factor.setPrefix("ρa × ")
        b = QPushButton("Apply shift")
        b.clicked.connect(lambda: self.shift(
            self._shift_comp.currentText().lower(),
            self._shift_factor.value()))
        for w in (self._shift_comp, self._shift_factor, b):
            shift.addWidget(w)
        hint = QLabel("or Shift + drag a resistivity curve")
        hint.setObjectName("InfoLabel")
        shift.addWidget(hint)
        shift.addStretch(1)
        self._log_lbl = QLabel("")
        self._log_lbl.setObjectName("InfoLabel")
        self._log_lbl.setWordWrap(True)
        v.addLayout(shift)
        v.addWidget(self._log_lbl)

        bb = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self._btn_apply = bb.addButton("Apply to survey",
                                       QDialogButtonBox.ButtonRole.AcceptRole)
        self._btn_apply.clicked.connect(self.accept)
        bb.rejected.connect(self.reject)
        v.addWidget(bb)
        QShortcut(QKeySequence("Left"), self, activated=lambda:
                  self.step_station(-1))
        QShortcut(QKeySequence("Right"), self, activated=lambda:
                  self.step_station(+1))
        QShortcut(QKeySequence("Escape"), self, activated=self.clear_selection)
        QShortcut(QKeySequence.StandardKey.Undo, self, activated=self.undo)
        self._update_status()

    # ── state ─────────────────────────────────────────────────────────
    @property
    def station(self) -> str:
        return self._station.currentText()

    def edited_stations(self) -> list[str]:
        return [k for k in self._order if self._cur[k] is not self._orig[k]]

    def edited_sites(self):
        """The whole survey with the edited stations swapped in."""
        from pycsamt.site.base import Sites

        return Sites([self._cur[k] for k in self._order])

    def _visible(self) -> list[str]:
        return [c for c in ("xy", "yx") if self._show[c].isChecked()]

    def step_station(self, step: int) -> None:
        i = self._order.index(self.station) + step
        if 0 <= i < len(self._order):
            self._station.setCurrentIndex(i)

    def _on_station(self, _name: str) -> None:
        self._sel.clear()
        self._draw()

    def clear_selection(self) -> None:
        self._sel.clear()
        self._draw()

    def select(self, rows, comps=None) -> None:
        """Select *rows* of *comps* (default: the visible components)."""
        for r in rows:
            for c in comps or self._visible():
                self._sel.add((int(r), c))
        self._draw()

    # ── editing ───────────────────────────────────────────────────────
    def _push_undo(self, label: str) -> None:
        self._undo.append((label, dict(self._cur),
                           {k: set(v) for k, v in self._interp.items()}))

    def undo(self) -> None:
        if not self._undo:
            return
        label, cur, interp = self._undo.pop()
        self._cur, self._interp = cur, interp
        self._sel.clear()
        self._log_lbl.setText(f"Undone: {label}")
        self._draw()

    def edit(self, op: str) -> None:
        """Apply *op* (mask, delete, interp, restore) to the selection."""
        if not self._sel:
            self._log_lbl.setText("Select points first (click or drag a "
                                  "box).")
            return
        rows = sorted({r for r, _c in self._sel})
        comps = sorted({c for _r, c in self._sel})
        freqs = pe.curves(self._cur[self.station])["freq"][rows]
        targets = (self._order if self._all_stations.isChecked()
                   else [self.station])
        label = {"mask": "mask", "delete": "delete", "interp":
                 "interpolate", "restore": "restore"}[op]
        self._push_undo(f"{label} {len(rows)} frequency(ies)")
        errors, skipped, done = [], 0, 0
        for name in targets:
            site = self._cur[name]
            r = (rows if name == self.station else pe.nearest_rows(
                pe.curves(site)["freq"], freqs))
            if not r:
                skipped += 1  # the station has none of these frequencies
                continue
            try:
                if op == "mask":
                    new = pe.mask(site, r, comps)
                elif op == "delete":
                    new = pe.delete(site, r)
                    self._interp.pop(name, None)
                elif op == "interp":
                    new = pe.interpolate(site, r, comps)
                    self._interp.setdefault(name, set()).update(
                        (x, c) for x in r for c in comps)
                else:
                    new = pe.restore(site, self._orig[name], r, comps)
                    self._interp.get(name, set()).difference_update(
                        (x, c) for x in r for c in comps)
            except Exception as exc:
                errors.append(f"{name}: {exc}")
                continue
            self._cur[name] = new
            done += 1
        self._sel.clear()
        note = (f" ({skipped} station(s) do not have these frequencies)"
                if skipped else "")
        self._log_lbl.setText(
            ("⚠ " + "; ".join(errors[:3])) if errors else
            f"{label.capitalize()}: {len(rows)} frequency(ies) on "
            f"{done} station(s){note}.")
        self._draw()

    def shift(self, comp: str, factor: float) -> None:
        if abs(factor - 1.0) < 1e-9:
            return
        self._push_undo(f"static shift {comp.upper()}")
        try:
            self._cur[self.station] = pe.static_shift(
                self._cur[self.station], comp, factor)
        except Exception as exc:
            self._undo.pop()
            QMessageBox.warning(self, "Static shift", str(exc))
            return
        self._log_lbl.setText(f"Static shift {comp.upper()}: ρa × "
                              f"{factor:.3g} on {self.station}.")
        self._draw()

    # ── drawing ───────────────────────────────────────────────────────
    def _draw(self) -> None:
        if not hasattr(self, "_canvas"):
            return
        fig = self._canvas.figure
        fig.clear()
        self._ax_rho, self._ax_phi = fig.subplots(2, 1, sharex=True,
                                                  gridspec_kw={"height_ratios"
                                                               : [3, 2]})
        site = self._cur.get(self.station)
        if site is None:
            return
        cv = pe.curves(site)
        f = cv["freq"]
        interp = self._interp.get(self.station, set())
        self._points = {}  # (row, comp) -> (f, rho, phi)
        for comp in self._visible():
            d = cv[comp]
            col = _COLOURS[comp]
            ok = d["valid"]
            for r in np.nonzero(ok)[0]:
                self._points[(int(r), comp)] = (f[r], d["rho"][r], d["phi"][r])
            normal = ok.copy()
            for r, c in interp:
                if c == comp and r < normal.size:
                    normal[r] = False
            for ax, key, err in ((self._ax_rho, "rho", "rho_err"),
                                 (self._ax_phi, "phi", "phi_err")):
                if self._errbars.isChecked():
                    ax.errorbar(f[normal], d[key][normal],
                                yerr=d[err][normal], fmt="none",
                                ecolor=col, alpha=0.35, lw=0.8)
                ax.plot(f[normal], d[key][normal], "o", ms=4, color=col,
                        label=comp.upper() if key == "rho" else None)
                iv = ok & ~normal
                if iv.any():
                    ax.plot(f[iv], d[key][iv], "D", ms=6, mfc="none",
                            mec=col, label=f"{comp.upper()} interpolated"
                            if key == "rho" else None)
            # masked rows: grey crosses at the original value
            orig = pe.curves(self._orig[self.station])
            mrows = [k for k in range(f.size) if not ok[k]]
            if mrows:
                of = orig["freq"]
                for r in mrows:
                    k = pe.nearest_rows(of, [f[r]], rtol=1e-6)
                    if k and orig[comp]["valid"][k[0]]:
                        self._ax_rho.plot(f[r], orig[comp]["rho"][k[0]], "x",
                                          color="#9ca3af", ms=6)
                        self._ax_phi.plot(f[r], orig[comp]["phi"][k[0]], "x",
                                          color="#9ca3af", ms=6)
        for (r, comp) in self._sel:
            if (r, comp) in self._points:
                fx, rh, ph = self._points[(r, comp)]
                self._ax_rho.plot(fx, rh, "o", ms=10, mfc="none", mec="#111",
                                  mew=1.6)
                self._ax_phi.plot(fx, ph, "o", ms=10, mfc="none", mec="#111",
                                  mew=1.6)
        self._ax_rho.set_xscale("log")
        self._ax_rho.set_yscale("log")
        self._ax_rho.set_ylabel("ρa (Ω·m)")
        self._ax_phi.set_ylabel("Phase (°)")
        self._ax_phi.set_xlabel("Frequency (Hz)")
        self._ax_phi.set_ylim(0, 90)
        self._ax_rho.invert_xaxis()  # high frequency (shallow) on the left
        self._ax_rho.set_title(f"{self.station} — {int(np.isfinite(f).sum())}"
                               " frequencies", fontsize=10)
        handles, labels = self._ax_rho.get_legend_handles_labels()
        if handles:
            self._ax_rho.legend(fontsize=8, loc="best")
        for ax in (self._ax_rho, self._ax_phi):
            ax.grid(True, which="both", alpha=0.25)
        self._canvas.draw()
        self._update_status()

    def _update_status(self) -> None:
        n = len(self._sel)
        self._sel_lbl.setText(f"{n} point(s) selected" if n
                              else "Nothing selected — click or drag a box")
        edited = self.edited_stations() if hasattr(self, "_cur") else []
        self._btn_apply.setEnabled(bool(edited))
        self._btn_apply.setText(f"Apply to survey ({len(edited)} station"
                                f"{'s' if len(edited) != 1 else ''})"
                                if edited else "Apply to survey")

    # ── mouse ─────────────────────────────────────────────────────────
    def _nearest(self, event):
        """(row, comp) of the point nearest the click, within 8 px."""
        best, dist = None, 8.0
        for key, (fx, rh, ph) in self._points.items():
            y = rh if event.inaxes is self._ax_rho else ph
            px = event.inaxes.transData.transform((fx, y))
            d = float(np.hypot(px[0] - event.x, px[1] - event.y))
            if d < dist:
                best, dist = key, d
        return best

    def _on_press(self, event) -> None:
        if event.inaxes not in (getattr(self, "_ax_rho", None),
                                getattr(self, "_ax_phi", None)):
            return
        if self._canvas._toolbar is not None and \
                getattr(self._canvas._toolbar, "mode", ""):
            return  # zoom / pan tool active
        shift = "shift" in (event.key or "")
        if shift and event.inaxes is self._ax_rho:
            hit = self._nearest(event)
            if hit is not None:
                self._drag = ("shift", hit[1], event.ydata)
            return
        self._drag = ("box", event.inaxes, event.xdata, event.ydata,
                      "control" in (event.key or ""))

    def _on_motion(self, event) -> None:
        if not self._drag or event.inaxes is None or event.ydata is None:
            return
        if self._drag[0] == "shift" and event.inaxes is self._ax_rho:
            factor = event.ydata / self._drag[2]
            self._log_lbl.setText(f"Static shift {self._drag[1].upper()}: "
                                  f"ρa × {factor:.3g} (release to apply)")

    def _on_release(self, event) -> None:
        drag, self._drag = self._drag, None
        if not drag or event.ydata is None:
            return
        if drag[0] == "shift":
            if event.inaxes is self._ax_rho and drag[2]:
                self.shift(drag[1], float(event.ydata / drag[2]))
            return
        _, ax, x0, y0, add = drag
        if event.inaxes is not ax or x0 is None:
            return
        small = (abs(np.log10(max(event.xdata, 1e-12) / max(x0, 1e-12)))
                 < 0.02)
        if small:  # a click: toggle the nearest point
            hit = self._nearest(event)
            if not add:
                keep = {hit} if hit in self._sel else set()
                self._sel = {s for s in self._sel if s in keep}
            if hit is not None:
                self._sel ^= {hit}
            self._draw()
            return
        lo_x, hi_x = sorted((x0, event.xdata))
        lo_y, hi_y = sorted((y0, event.ydata))
        if not add:
            self._sel.clear()
        for key, (fx, rh, ph) in self._points.items():
            y = rh if ax is self._ax_rho else ph
            if lo_x <= fx <= hi_x and lo_y <= y <= hi_y:
                self._sel.add(key)
        self._draw()


__all__ = ["PointEditorDialog"]
