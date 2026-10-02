# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Pattern packs (Geology Studio "Patterns" tab) callback tests."""

from __future__ import annotations

from types import SimpleNamespace

import pytest

fitz = pytest.importorskip("fitz")

from pycsamt.app.mapview.callbacks import patterns as patterns_cb  # noqa: E402


def _capture():
    captured: dict = {}

    class _App:
        def callback(self, *a, **k):
            def deco(fn):
                captured[fn.__name__] = fn
                return fn

            return deco

    patterns_cb.register_patterns(_App())
    return captured


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_PATTERN_CACHE", str(tmp_path / "cache"))
    yield


def _make_pdf(path) -> None:
    doc = fitz.open()
    page = doc.new_page(width=200, height=200)
    page.draw_rect(fitz.Rect(20, 20, 180, 180))
    doc.save(str(path))
    doc.close()


def test_registers_expected_callbacks():
    captured = _capture()
    assert {
        "pack_options", "do_import", "sheet_options", "figure",
        "save_crop", "grid", "select", "assign",
    } <= set(captured)


def test_pack_options_lists_builtin_by_default():
    captured = _capture()
    options, value = captured["pack_options"](None, None, None)
    ids = {o["value"] for o in options}
    assert "builtin" in ids
    assert value == "builtin"


def test_import_from_path_populates_pack_select(tmp_path, monkeypatch):
    captured = _capture()
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")

    from pycsamt.app.mapview._ids import IDs

    monkeypatch.setattr(
        patterns_cb, "ctx",
        SimpleNamespace(triggered_id=IDs.GEO_BTN_IMPORT_PACK_PATH),
    )
    refresh, pack_id, status = captured["do_import"](
        None, 1, None, str(src), "My Pack", 0,
    )
    assert refresh == 1
    assert pack_id == "my-pack"
    assert "Imported" in status.children


def test_import_reports_error_for_missing_path(monkeypatch):
    captured = _capture()
    from pycsamt.app.mapview._ids import IDs

    monkeypatch.setattr(
        patterns_cb, "ctx",
        SimpleNamespace(triggered_id=IDs.GEO_BTN_IMPORT_PACK_PATH),
    )
    refresh, pack_id, status = captured["do_import"](
        None, 1, None, "", None, 0,
    )
    from dash import no_update

    assert refresh is no_update
    assert "path" in status.children.lower() or status.children


def test_sheet_options_lists_pages_of_the_selected_pack(tmp_path, monkeypatch):
    from pycsamt.app._patterns import import_pack_path

    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    summary = import_pack_path(str(src), name="Sheet Pack")

    captured = _capture()
    options, value = captured["sheet_options"](summary["id"], 0)
    assert len(options) == 1
    assert value == options[0]["value"]


def test_sheet_figure_embeds_image(tmp_path):
    from pycsamt.app._patterns import import_pack_path

    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    summary = import_pack_path(str(src), name="Fig Pack")

    captured = _capture()
    fig = captured["figure"](summary["sheets"][0]["id"], summary["id"])
    assert len(fig.layout.images) == 1


def test_save_crop_requires_a_drawn_shape(tmp_path):
    from pycsamt.app._patterns import import_pack_path

    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    summary = import_pack_path(str(src), name="Crop Pack")

    captured = _capture()
    refresh, status = captured["save_crop"](
        1, {}, summary["id"], summary["sheets"][0]["id"], "Swatch", 0,
    )
    from dash import no_update

    assert refresh is no_update
    assert "draw" in status.children.lower()


def test_save_crop_with_a_drawn_shape_creates_a_tile(tmp_path):
    from pycsamt.app._patterns import import_pack_path, load_pack

    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    summary = import_pack_path(str(src), name="Crop Pack 2")

    captured = _capture()
    relayout = {"shapes": [{"x0": 10, "y0": 10, "x1": 100, "y1": 100}]}
    refresh, status = captured["save_crop"](
        1, relayout, summary["id"], summary["sheets"][0]["id"],
        "My Swatch", 0,
    )
    assert refresh == 1
    assert "Saved" in status.children
    assert len(load_pack(summary["id"]).tiles) == 1


def test_grid_lists_tiles_for_the_selected_pack():
    captured = _capture()
    cells = captured["grid"]("builtin", 0, None)
    assert len(cells) >= 6


def test_swatch_select_updates_store_from_triggered_id(monkeypatch):
    captured = _capture()
    monkeypatch.setattr(
        patterns_cb, "ctx",
        SimpleNamespace(triggered_id={"type": "geo-swatch-btn", "pack": "builtin", "tile": "builtin-01"}),
    )
    result = captured["select"]([1])
    assert result == {"pack": "builtin", "tile": "builtin-01"}


def test_assign_requires_a_selected_row_and_swatch():
    captured = _capture()
    from dash import no_update

    rows = [{"name": "Clay", "pattern_id": ""}]
    result, status = captured["assign"](1, rows, [], {"pack": "builtin", "tile": "builtin-01"})
    assert result is no_update
    assert "select a row" in status.children.lower()

    result, status = captured["assign"](1, rows, [0], None)
    assert result is no_update
    assert "swatch" in status.children.lower()


def test_assign_patches_the_selected_legend_row():
    captured = _capture()
    rows = [
        {"name": "Clay", "pattern_id": ""},
        {"name": "Granite", "pattern_id": ""},
    ]
    result, status = captured["assign"](
        1, rows, [1], {"pack": "builtin", "tile": "builtin-01"}
    )
    assert result[0]["pattern_id"] == ""
    assert result[1]["pattern_id"] == "builtin-01"
    assert result[1]["pattern_source"] == "builtin"
    assert "Assigned" in status.children
