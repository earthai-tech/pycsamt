# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
LayeredModelDialog — interactive 1-D layered earth model builder.

Lets the user define N earth layers (resistivity + thickness) via an editable
table, with live preview of the resistivity–depth profile.  Factory presets
(Random, Blocky, Smooth) populate the table automatically.

Calls ``pycsamt.forward.synthetic.LayeredModel``.

Usage
-----
    from pycsamt.app.desktop.tools.layered_model_tool import LayeredModelDialog
    dlg = LayeredModelDialog(parent=self)
    dlg.exec()
    model = dlg.model   # LayeredModel or None
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

_PRESETS = ["— choose preset —", "Random", "Blocky", "Smooth"]
_DEFAULT_LAYERS = [
    (100.0, 300.0),
    (10.0, 800.0),
    (500.0, None),  # halfspace
]


# ── Dialog ────────────────────────────────────────────────────────────────────


class LayeredModelDialog(QDialog):
    """
    Interactive 1-D earth model builder.

    Attributes
    ----------
    model : LayeredModel or None
        The last built model; available after the dialog is closed.

    Parameters
    ----------
    parent : QWidget, optional
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Layered Model Builder")
        self.setMinimumSize(820, 520)
        self.model = None
        self._build_ui()
        self._load_default()

    # ── Build ─────────────────────────────────────────────────────────────────

    def _build_ui(self) -> None:
        root = QVBoxLayout(self)
        root.setSpacing(8)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ── Left: layer table + controls ──────────────────────────────────
        left = QWidget()
        left.setFixedWidth(380)
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(6, 6, 6, 6)
        left_lay.setSpacing(8)

        # Preset row
        preset_row = QHBoxLayout()
        preset_row.addWidget(QLabel("Preset:"))
        self._preset_combo = QComboBox()
        self._preset_combo.addItems(_PRESETS)
        self._preset_combo.currentIndexChanged.connect(self._on_preset)
        preset_row.addWidget(self._preset_combo)

        self._n_spin = QSpinBox()
        self._n_spin.setRange(2, 20)
        self._n_spin.setValue(3)
        self._n_spin.setPrefix("N = ")
        preset_row.addWidget(self._n_spin)
        left_lay.addLayout(preset_row)

        # Layer table: columns = [Layer, ρ (Ω·m), Thickness (m)]
        self._table = QTableWidget(0, 3)
        self._table.setHorizontalHeaderLabels(
            ["Layer", "ρ (Ω·m)", "Thickness (m)"]
        )
        self._table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.ResizeToContents
        )
        self._table.horizontalHeader().setSectionResizeMode(
            1, QHeaderView.ResizeMode.Stretch
        )
        self._table.horizontalHeader().setSectionResizeMode(
            2, QHeaderView.ResizeMode.Stretch
        )
        self._table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection
        )
        left_lay.addWidget(self._table)

        # Add / Remove row buttons
        btn_row = QHBoxLayout()
        add_btn = QPushButton("+ Add layer")
        add_btn.clicked.connect(self._add_row)
        rem_btn = QPushButton("− Remove last")
        rem_btn.clicked.connect(self._remove_row)
        btn_row.addWidget(add_btn)
        btn_row.addWidget(rem_btn)
        left_lay.addLayout(btn_row)

        # Plot button
        self._plot_btn = QPushButton("Preview model")
        self._plot_btn.clicked.connect(self._on_preview)
        left_lay.addWidget(self._plot_btn)

        self._status_lbl = QLabel("")
        self._status_lbl.setWordWrap(True)
        left_lay.addWidget(self._status_lbl)
        left_lay.addStretch()

        splitter.addWidget(left)

        # ── Right: canvas ─────────────────────────────────────────────────
        self._canvas = MplCanvas(parent=self)
        splitter.addWidget(self._canvas)
        splitter.setStretchFactor(1, 1)

        root.addWidget(splitter, stretch=1)

        box = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Close
        )
        box.accepted.connect(self._on_ok)
        box.rejected.connect(self.reject)
        root.addWidget(box)

    # ── Table helpers ─────────────────────────────────────────────────────────

    def _load_default(self) -> None:
        self._table.setRowCount(0)
        for i, (rho, thick) in enumerate(_DEFAULT_LAYERS):
            self._append_row(i + 1, rho, thick)
        self._on_preview()

    def _append_row(self, idx: int, rho: float, thick) -> None:
        r = self._table.rowCount()
        self._table.insertRow(r)
        label = "Halfspace" if thick is None else f"Layer {idx}"
        self._table.setItem(r, 0, QTableWidgetItem(label))
        self._table.setItem(r, 1, QTableWidgetItem(f"{rho:.1f}"))
        self._table.setItem(
            r, 2, QTableWidgetItem("" if thick is None else f"{thick:.1f}")
        )

    def _add_row(self) -> None:
        n = self._table.rowCount()
        self._append_row(n, 100.0, 500.0)
        self._on_preview()

    def _remove_row(self) -> None:
        n = self._table.rowCount()
        if n > 2:
            self._table.removeRow(n - 1)
            self._on_preview()

    # ── Presets ───────────────────────────────────────────────────────────────

    def _on_preset(self, idx: int) -> None:
        preset = self._preset_combo.currentText()
        n = self._n_spin.value()
        if preset == "Random":
            self._apply_model_preset("random", n)
        elif preset == "Blocky":
            self._apply_model_preset("blocky", n)
        elif preset == "Smooth":
            self._apply_model_preset("smooth", n)

    def _apply_model_preset(self, kind: str, n: int) -> None:
        try:
            from pycsamt.forward.synthetic import LayeredModel

            if kind == "random":
                m = LayeredModel.random(n)
            elif kind == "blocky":
                m = LayeredModel.blocky(n)
            else:
                m = LayeredModel.smooth(n)
        except Exception as exc:
            self._status_lbl.setText(f"Preset error: {exc}")
            return

        self._table.setRowCount(0)
        for i, rho in enumerate(m.resistivity):
            thick = m.thickness[i] if i < len(m.thickness) else None
            self._append_row(i + 1, rho, thick)
        self._on_preview()

    # ── Build model from table ─────────────────────────────────────────────────

    def _read_model(self):
        """Parse table → LayeredModel, or raise ValueError."""
        from pycsamt.forward.synthetic import LayeredModel

        rhos, thicks = [], []
        n = self._table.rowCount()
        for r in range(n):
            rho_item = self._table.item(r, 1)
            thick_item = self._table.item(r, 2)
            try:
                rhos.append(float(rho_item.text()))
            except (ValueError, AttributeError):
                raise ValueError(f"Invalid resistivity in row {r + 1}")
            thick_text = thick_item.text().strip() if thick_item else ""
            if r < n - 1:
                try:
                    thicks.append(float(thick_text))
                except ValueError:
                    raise ValueError(f"Invalid thickness in row {r + 1}")
        return LayeredModel(resistivity=rhos, thickness=thicks)

    # ── Preview ───────────────────────────────────────────────────────────────

    def _on_preview(self) -> None:
        try:
            model = self._read_model()
        except Exception as exc:
            self._status_lbl.setText(str(exc))
            return

        fig, ax = plt.subplots(figsize=(4, 5))
        model.plot(ax=ax)
        fig.tight_layout()
        self._canvas.show_figure(fig)
        self._status_lbl.setText(
            f"{model.n_layers} layers  |  "
            f"ρ: {model.resistivity.min():.1f}–{model.resistivity.max():.1f} Ω·m  |  "
            f"depth: {model.depth[-1]:.0f} m"
        )

    # ── OK ────────────────────────────────────────────────────────────────────

    def _on_ok(self) -> None:
        try:
            self.model = self._read_model()
            self.accept()
        except Exception as exc:
            self._status_lbl.setText(f"Cannot build model: {exc}")
