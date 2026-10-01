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


def test_dpi_spinbox_default_publication(dlg):
    assert dlg._dpi_spin.value() == 600


def test_dpi_spinbox_range(dlg):
    assert dlg._dpi_spin.minimum() == 36
    assert dlg._dpi_spin.maximum() == 2400


def test_path_field_has_default(dlg):
    assert dlg._path_edit.text() != ""


def test_path_extension_updates_on_format_change(dlg):
    dlg._fmt_combo.setCurrentText("PDF  (vector)")
    dlg._update_path_extension("PDF  (vector)")
    assert dlg._path_edit.text().endswith(".pdf")


def test_export_reads_api_defaults(qapp, simple_fig, monkeypatch):
    from pycsamt.api.plot import PLOT_CONFIG

    with PLOT_CONFIG.context(fmt="svg", dpi=1200, transparent=True):
        d = ExportDialog(figure=simple_fig)
        assert d._dpi_spin.value() == 1200
        assert d._path_edit.text().endswith(".svg")
        saved = {}
        monkeypatch.setattr(simple_fig, "savefig", lambda path, **kw: saved.update(kw))
        d._on_export()
        assert saved["format"] == "svg"
        assert saved["transparent"] is True
        assert saved["dpi"] == 1200
        d.close()


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


# ── Figure sources (Export toolbar button) ────────────────────────────────


def _window_with_figure(title, draw=True):
    from PySide6.QtWidgets import QTabWidget, QVBoxLayout, QWidget

    from pycsamt.app.desktop.widgets.mpl_canvas import MplCanvas

    win = QWidget()
    win.setWindowTitle(title)
    tabs = QTabWidget(win)
    canvas = MplCanvas()
    tabs.addTab(canvas, "Plot")
    QVBoxLayout(win).addWidget(tabs)
    if draw:
        canvas.figure.clf()
        ax = canvas.figure.add_subplot(111)
        ax.plot([1, 2], [2, 1])
        ax.set_title("Decay curves")
    win.show()
    return win, canvas.figure


def test_collect_sources_finds_drawn_canvases_only(qapp):
    from pycsamt.app.desktop.dialogs.export_dlg import collect_sources

    a, fig_a = _window_with_figure("TDEM Studio")
    b, _fig_b = _window_with_figure("Empty", draw=False)
    hidden, _fig_h = _window_with_figure("Hidden")
    hidden.hide()
    try:
        sources = collect_sources([a, b, hidden])
        assert [s.figure for s in sources] == [fig_a]
        assert sources[0].label == "TDEM Studio  ▸  Plot  ▸  Decay curves"
        assert sources[0].kind == "mpl"
    finally:
        for w in (a, b, hidden):
            w.close()


def test_collect_sources_puts_first_window_first(qapp):
    from pycsamt.app.desktop.dialogs.export_dlg import collect_sources

    a, _fa = _window_with_figure("A")
    b, fig_b = _window_with_figure("B")
    try:
        assert collect_sources([a, b], first=b)[0].figure is fig_b
    finally:
        a.close()
        b.close()


def test_dialog_with_sources_switches_figure(qapp, tmp_path):
    from pycsamt.app.desktop.dialogs.export_dlg import (
        FigureSource,
        ExportDialog as Dlg,
    )

    f1, ax1 = plt.subplots(figsize=(4, 3))
    ax1.plot([1, 2])
    f2, ax2 = plt.subplots(figsize=(8, 5))
    ax2.plot([2, 1])
    srcs = [FigureSource("Win ▸ Plot ▸ One", figure=f1),
            FigureSource("Win ▸ Plot ▸ Two", figure=f2)]
    d = Dlg(sources=srcs)
    try:
        assert d._source_combo.count() == 2
        assert d._figure is f1 and "4.0 × 3.0 in" in d._size_lbl.text()
        d._source_combo.setCurrentIndex(1)
        assert d._figure is f2 and "8.0 × 5.0 in" in d._size_lbl.text()
        assert "pycsamt_win_two" in d._path_edit.text()
        d._dpi_spin.setValue(100)
        assert "800 × 500 px" in d._size_lbl.text()
        out = tmp_path / "two.png"
        d._path_edit.setText(str(out))
        d._on_export()
        assert out.exists()
    finally:
        d.close()
        plt.close(f1)
        plt.close(f2)


def test_dialog_plotly_source_exports_html(qapp, tmp_path):
    from pycsamt.app.desktop.dialogs.export_dlg import (
        _PLOTLY_FORMATS,
        FigureSource,
        ExportDialog as Dlg,
    )

    class _View:
        def export_html(self, path):
            open(path, "w").write("<html></html>")

    d = Dlg(sources=[FigureSource("PCSF 3-D ▸ 3-D scene", plotly=_View())])
    try:
        items = [d._fmt_combo.itemText(i) for i in range(d._fmt_combo.count())]
        assert items == list(_PLOTLY_FORMATS)
        assert not d._dpi_spin.isEnabled()
        d._fmt_combo.setCurrentText("HTML (interactive 3-D)")
        assert d._path_edit.text().endswith(".html")
        out = tmp_path / "scene.html"
        d._path_edit.setText(str(out))
        d._on_export()
        assert out.exists()
    finally:
        d.close()
