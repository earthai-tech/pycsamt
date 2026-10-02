# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for PlotlyView (pycsamt.app.desktop.widgets.plotly_view).

Renders a real Plotly figure to a temp HTML file and loads it through
QWebEngineView -- these tests check the file-based render/export/cleanup
plumbing, not the browser's own rendering (out of reach under pytest).
"""

from __future__ import annotations

import pytest

pytest.importorskip("PySide6", reason="PySide6 required")
pytest.importorskip(
    "PySide6.QtWebEngineWidgets", reason="QtWebEngine required"
)
pytest.importorskip("plotly", reason="plotly required")

from pycsamt.app.desktop.widgets.plotly_view import PlotlyView


def _sample_figure():
    import plotly.graph_objects as go

    return go.Figure(
        data=[go.Scatter3d(x=[0, 1, 2], y=[0, 1, 2], z=[0, 1, 2], mode="markers")]
    )


def test_constructs_with_a_placeholder(qapp):
    view = PlotlyView()
    assert view.figure is None
    assert view._html_path.exists()
    assert "Nothing to display" in view._html_path.read_text(encoding="utf-8")
    view.close()


def test_set_figure_renders_real_html(qapp):
    view = PlotlyView()
    fig = _sample_figure()
    view.set_figure(fig)
    assert view.figure is fig
    html = view._html_path.read_text(encoding="utf-8")
    assert "Scatter3d" in html or "scatter3d" in html
    assert len(html) > 1_000_000  # inline plotly.js bundle is a few MB
    view.close()


def test_clear_resets_to_placeholder(qapp):
    view = PlotlyView()
    view.set_figure(_sample_figure())
    view.clear()
    assert view.figure is None
    assert "Nothing to display" in view._html_path.read_text(encoding="utf-8")
    view.close()


def test_export_html_writes_a_standalone_file(qapp, tmp_path):
    view = PlotlyView()
    view.set_figure(_sample_figure())
    out_path = tmp_path / "exported.html"
    view.export_html(str(out_path))
    assert out_path.exists()
    assert out_path.stat().st_size > 1_000_000
    view.close()


def test_export_html_raises_without_a_figure(qapp, tmp_path):
    view = PlotlyView()
    with pytest.raises(ValueError):
        view.export_html(str(tmp_path / "out.html"))
    view.close()


def test_close_cleans_up_temp_directory(qapp):
    view = PlotlyView()
    tmp_dir = view._tmp_dir
    assert tmp_dir.exists()
    view.close()
    assert not tmp_dir.exists()


def test_each_instance_gets_its_own_temp_file(qapp):
    view1 = PlotlyView()
    view2 = PlotlyView()
    assert view1._html_path != view2._html_path
    assert view1._tmp_dir != view2._tmp_dir
    view1.close()
    view2.close()
