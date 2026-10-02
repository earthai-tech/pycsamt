# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Settings forms generated from engine :class:`Field` declarations.

Each engine declares its settings once (``Engine.fields()``); the studio
builds one ``SettingsForm`` per step section (Data / Mesh / Settings /
Run), so a new engine or a new option needs no window code.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QLineEdit,
    QSizePolicy,
    QSpinBox,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from pycsamt.app.desktop.controllers.inversion_engines import Field


def _make_widget(f: Field) -> QWidget:
    if f.kind == "bool":
        w = QCheckBox()
        w.setChecked(bool(f.default))
    elif f.kind == "choice":
        w = QComboBox()
        for value, label in f.choices:
            w.addItem(label, value)
        w.setCurrentIndex(max(w.findData(f.default), 0))
    elif f.kind == "text":  # lists, pairs, optional numbers (blank = auto)
        w = QLineEdit(str(f.default if f.default is not None else ""))
        w.setPlaceholderText("Automatic")
    elif f.kind == "int":
        w = QSpinBox()
        w.setRange(int(f.lo if f.lo is not None else -10**9),
                   int(f.hi if f.hi is not None else 10**9))
        w.setValue(int(f.default))
        if f.unit:
            w.setSuffix(f" {f.unit}")
    else:
        w = QDoubleSpinBox()
        w.setDecimals(f.decimals)
        w.setRange(float(f.lo if f.lo is not None else -1e12),
                   float(f.hi if f.hi is not None else 1e12))
        if f.step:
            w.setSingleStep(float(f.step))
        w.setValue(float(f.default))
        if f.unit:
            w.setSuffix(f" {f.unit}")
        if f.auto_zero:
            w.setSpecialValueText("Auto")
    if isinstance(w, (QSpinBox, QDoubleSpinBox, QComboBox, QLineEdit)):
        w.setSizePolicy(QSizePolicy.Policy.Expanding,
                        QSizePolicy.Policy.Fixed)
    if f.help:
        w.setToolTip(f.help)
    return w


def _value(w: QWidget):
    if isinstance(w, QCheckBox):
        return w.isChecked()
    if isinstance(w, QComboBox):
        return w.currentData()
    if isinstance(w, QLineEdit):
        return w.text()
    return w.value()


def _set(w: QWidget, value) -> None:
    if isinstance(w, QCheckBox):
        w.setChecked(bool(value))
    elif isinstance(w, QComboBox):
        i = w.findData(value)
        if i >= 0:
            w.setCurrentIndex(i)
    elif isinstance(w, QLineEdit):
        w.setText(str(value))
    elif isinstance(w, QSpinBox):
        w.setValue(int(value))
    else:
        w.setValue(float(value))


def _changed_signal(w: QWidget):
    if isinstance(w, QCheckBox):
        return w.toggled
    if isinstance(w, QComboBox):
        return w.currentIndexChanged
    if isinstance(w, QLineEdit):
        return w.editingFinished
    return w.valueChanged


class SettingsForm(QWidget):
    """Form for the fields of one section; advanced ones fold away."""

    changed = Signal()

    def __init__(self, fields: list[Field],
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.fields = list(fields)
        self.widgets: dict[str, QWidget] = {}
        v = QVBoxLayout(self)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        basic = QFormLayout()
        basic.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        basic.setSpacing(5)
        v.addLayout(basic)
        adv_fields = [f for f in self.fields if f.advanced]
        self._adv_box = QWidget()
        adv = QFormLayout(self._adv_box)
        adv.setContentsMargins(0, 0, 0, 0)
        adv.setSpacing(5)
        for f in self.fields:
            w = _make_widget(f)
            self.widgets[f.key] = w
            _changed_signal(w).connect(lambda *_: self.changed.emit())
            (adv if f.advanced else basic).addRow(f"{f.label}:", w)
        self.btn_more = QToolButton()
        self.btn_more.setCheckable(True)
        self.btn_more.setAutoRaise(True)
        self.btn_more.setToolButtonStyle(
            Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.btn_more.setArrowType(Qt.ArrowType.RightArrow)
        self.btn_more.setText(f"More settings ({len(adv_fields)})")
        self.btn_more.toggled.connect(self._toggle_adv)
        self.btn_more.setVisible(bool(adv_fields))
        v.addWidget(self.btn_more)
        v.addWidget(self._adv_box)
        self._adv_box.setVisible(False)

    def _toggle_adv(self, on: bool) -> None:
        self.btn_more.setArrowType(Qt.ArrowType.DownArrow if on
                                   else Qt.ArrowType.RightArrow)
        self._adv_box.setVisible(on)

    def values(self) -> dict:
        return {k: _value(w) for k, w in self.widgets.items()}

    def set_values(self, values: dict) -> None:
        for k, v in values.items():
            w = self.widgets.get(k)
            if w is not None and v is not None:
                w.blockSignals(True)
                try:
                    _set(w, v)
                finally:
                    w.blockSignals(False)
        self.changed.emit()

    def reset(self) -> None:
        self.set_values({f.key: f.default for f in self.fields})

    def set_enabled(self, on: bool) -> None:
        for w in self.widgets.values():
            w.setEnabled(on)


__all__ = ["SettingsForm"]
