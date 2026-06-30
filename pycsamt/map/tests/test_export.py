# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for map figure export helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.map import ExportOptions
from pycsamt.map.export import (
    export_figure,
    figure_to_dict,
    save_png,
    write_dict,
    write_html,
    write_image,
    write_json,
)


class _HtmlFigure:
    def __init__(self) -> None:
        self.calls = []

    def write_html(self, *args, **kwargs) -> None:
        self.calls.append((args, kwargs))


class _ImageFigure:
    def __init__(self) -> None:
        self.calls = []

    def write_image(self, *args, **kwargs) -> None:
        self.calls.append((args, kwargs))


class _KaleidoFigure:
    def write_image(self, *_args, **_kwargs) -> None:
        raise ValueError("kaleido is required")


class _MatplotlibFigure:
    def __init__(self) -> None:
        self.calls = []

    def savefig(self, *args, **kwargs) -> None:
        if "scale" in kwargs:
            raise AssertionError("scale leaked to savefig")
        self.calls.append((args, kwargs))


class _JsonFigure:
    def __init__(self) -> None:
        self.path = ""

    def write_json(self, path: str) -> None:
        self.path = path


class _DictFigure:
    def to_dict(self):
        return {"data": [], "layout": {"title": "map"}}


def test_write_html_routes_to_plotly_method() -> None:
    fig = _HtmlFigure()
    out = write_html(fig, "map.html", include_plotlyjs="cdn")
    assert out == Path("map.html")
    assert fig.calls[0][1]["include_plotlyjs"] == "cdn"


def test_write_image_explains_kaleido() -> None:
    with pytest.raises(ImportError, match="kaleido"):
        write_image(_KaleidoFigure(), "map.png")


def test_save_png_filters_plotly_kwargs() -> None:
    fig = _MatplotlibFigure()
    out = save_png(
        fig,
        "map.png",
        scale=3,
        width=800,
        dpi=150,
    )
    assert out == Path("map.png")
    assert fig.calls[0][1] == {"format": "png", "dpi": 150}


def test_export_figure_uses_format() -> None:
    fig = _HtmlFigure()
    out = export_figure(
        fig,
        ExportOptions(path="map", format="html"),
    )
    assert out == Path("map")
    assert fig.calls


def test_export_figure_rejects_missing_format() -> None:
    with pytest.raises(ValueError, match="extension"):
        export_figure(
            _HtmlFigure(),
            ExportOptions(path="map"),
        )


def test_write_json_prefers_write_json_method() -> None:
    fig = _JsonFigure()
    out = write_json(fig, "map.json")
    assert out == Path("map.json")
    assert fig.path == "map.json"


def test_figure_to_dict_and_write_dict(monkeypatch) -> None:
    writes = {}

    def fake_write_text(self, text, encoding=None):
        writes[str(self)] = (text, encoding)

    monkeypatch.setattr(Path, "write_text", fake_write_text)
    fig = _DictFigure()
    assert figure_to_dict(fig)["layout"]["title"] == "map"
    out = write_dict(fig, "map.dict")
    assert out == Path("map.dict")
    assert '"title": "map"' in writes["map.dict"][0]
