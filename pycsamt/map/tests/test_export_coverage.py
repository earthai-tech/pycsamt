# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional coverage for :mod:`pycsamt.map.export`.

``test_export.py`` covers the primary happy/error paths for each public
function individually. This file rounds out the branches it misses:
the non-kaleido re-raise in ``write_image``, its normal return, the
Plotly-fallback branch of ``save_png``, the json/dict/unknown-suffix
branches of ``export_figure``, the ``to_json``/dict-fallback branches
of ``write_json``, and the ``to_plotly_json`` branches of
``figure_to_dict``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from pycsamt.map import ExportOptions
from pycsamt.map.export import (
    export_figure,
    figure_to_dict,
    save_png,
    write_image,
    write_json,
)


class _ImageFigure:
    def __init__(self, exc=None):
        self._exc = exc
        self.calls = []

    def write_image(self, path, **kwargs):
        if self._exc is not None:
            raise self._exc
        self.calls.append((path, kwargs))


class _PlotlyLikeFigure:
    """No savefig — save_png must fall back to write_image."""

    def __init__(self):
        self.calls = []

    def write_image(self, path, **kwargs):
        self.calls.append((path, kwargs))


class _ToJsonFigure:
    def __init__(self):
        self.text = None

    def to_json(self):
        return '{"a": 1}'


class _ToDictOnlyFigure:
    def to_dict(self):
        return {"layout": {"title": "m"}}


class _ToPlotlyJsonFigure:
    def to_plotly_json(self):
        return {"data": [1, 2, 3]}


class _NonDictToDictFigure:
    """to_dict returns something that isn't a dict -> falls through."""

    def to_dict(self):
        return "not-a-dict"

    def to_plotly_json(self):
        return {"ok": True}


class _NeitherFigure:
    pass


class _MatplotlibFigure:
    def __init__(self) -> None:
        self.calls = []

    def savefig(self, *args, **kwargs) -> None:
        self.calls.append((args, kwargs))


def test_write_image_returns_path_on_success(tmp_path) -> None:
    fig = _ImageFigure()
    out = write_image(fig, tmp_path / "a.png")
    assert out == tmp_path / "a.png"
    assert fig.calls


def test_write_image_reraises_non_kaleido_error(tmp_path) -> None:
    fig = _ImageFigure(exc=RuntimeError("some other failure"))
    with pytest.raises(RuntimeError, match="some other failure"):
        write_image(fig, tmp_path / "b.png")


def test_save_png_delegates_to_write_image_for_plotly_like(tmp_path) -> None:
    fig = _PlotlyLikeFigure()
    out = save_png(fig, tmp_path / "c.png", scale=2)
    assert out == tmp_path / "c.png"
    assert fig.calls
    assert fig.calls[0][1]["format"] == "png"


def test_export_figure_png_with_width_and_height(tmp_path) -> None:
    fig = _MatplotlibFigure()
    out = export_figure(
        fig,
        ExportOptions(path=tmp_path / "d.png", width=640, height=480),
    )
    assert out == tmp_path / "d.png"
    assert fig.calls


def test_export_figure_json_suffix(tmp_path) -> None:
    fig = _ToDictOnlyFigure()
    out = export_figure(fig, ExportOptions(path=tmp_path / "e.json"))
    assert out == tmp_path / "e.json"
    assert out.read_text(encoding="utf-8")


def test_export_figure_dict_suffix(tmp_path) -> None:
    fig = _ToDictOnlyFigure()
    out = export_figure(fig, ExportOptions(path=tmp_path / "f.dict"))
    assert out == tmp_path / "f.dict"


def test_export_figure_unknown_suffix_falls_back_to_write_image(
    tmp_path,
) -> None:
    fig = _ImageFigure()
    out = export_figure(fig, ExportOptions(path=tmp_path / "g.svg"))
    assert out == tmp_path / "g.svg"
    assert fig.calls


def test_write_json_uses_to_json_when_no_write_json(tmp_path) -> None:
    fig = _ToJsonFigure()
    out = write_json(fig, tmp_path / "h.json")
    assert out.read_text(encoding="utf-8") == '{"a": 1}'


def test_write_json_falls_back_to_figure_to_dict(tmp_path) -> None:
    fig = _ToDictOnlyFigure()
    out = write_json(fig, tmp_path / "i.json")
    assert '"title": "m"' in out.read_text(encoding="utf-8")


def test_figure_to_dict_uses_to_plotly_json_when_no_to_dict() -> None:
    fig = _ToPlotlyJsonFigure()
    assert figure_to_dict(fig) == {"data": [1, 2, 3]}


def test_figure_to_dict_falls_through_non_dict_to_dict() -> None:
    fig = _NonDictToDictFigure()
    assert figure_to_dict(fig) == {"ok": True}


def test_figure_to_dict_raises_when_neither_available() -> None:
    with pytest.raises(TypeError, match="to_dict or to_plotly_json"):
        figure_to_dict(_NeitherFigure())
