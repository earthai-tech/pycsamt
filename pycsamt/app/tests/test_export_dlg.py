# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for ExportDialog (Phase 4)."""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pycsamt.app.desktop.dialogs.export_dlg import (
    _FORMATS,
    ExportDialog,
)


@pytest.fixture
def simple_fig():
    fig, ax = plt.subplots()
    ax.plot([1, 2, 3], [4, 5, 6])
    yield fig
    plt.close(fig)


@pytest.fixture
def dlg(qapp, simple_fig):
    d = ExportDialog(figure=simple_fig)
    yield d
    d.close()


# ── Construction ──────────────────────────────────────────────────────────

def test_export_dlg_creates(qapp, simple_fig):
    d = ExportDialog(figure=simple_fig)
    assert d is not None
    d.close()


def test_format_combo_has_all_formats(dlg):
    combo_items = [dlg._fmt_combo.itemText(i) for i in range(dlg._fmt_combo.count())]
    for fmt_key in _FORMATS:
        assert fmt_key in combo_items


def test_dpi_spinbox_default_300(dlg):
    assert dlg._dpi_spin.value() == 300


def test_dpi_spinbox_range(dlg):
    assert dlg._dpi_spin.minimum() == 36
    assert dlg._dpi_spin.maximum() == 1200


def test_path_field_has_default(dlg):
    assert dlg._path_edit.text() != ""


def test_path_extension_updates_on_format_change(dlg):
    dlg._fmt_combo.setCurrentText("PDF  (vector)")
    dlg._update_path_extension("PDF  (vector)")
    assert dlg._path_edit.text().endswith(".pdf")


# ── Export to disk ────────────────────────────────────────────────────────

def test_export_png_creates_file(qapp, simple_fig, tmp_path):
    out = tmp_path / "test.png"
    d = ExportDialog(figure=simple_fig, default_path=str(out))
    d._fmt_combo.setCurrentText("PNG  (raster, lossless)")
    d._dpi_spin.setValue(72)
    d._path_edit.setText(str(out))
    d._on_export()
    assert out.exists()
    d.close()


def test_export_svg_creates_file(qapp, simple_fig, tmp_path):
    out = tmp_path / "test.svg"
    d = ExportDialog(figure=simple_fig, default_path=str(out))
    d._fmt_combo.setCurrentText("SVG  (vector)")
    d._path_edit.setText(str(out))
    d._on_export()
    assert out.exists()
    d.close()


def test_export_pdf_creates_file(qapp, simple_fig, tmp_path):
    out = tmp_path / "test.pdf"
    d = ExportDialog(figure=simple_fig, default_path=str(out))
    d._fmt_combo.setCurrentText("PDF  (vector)")
    d._path_edit.setText(str(out))
    d._on_export()
    assert out.exists()
    d.close()
