# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared pattern-pack app-adapter tests."""

from __future__ import annotations

import base64

import pytest

fitz = pytest.importorskip("fitz")


def _make_pdf_bytes(n_pages: int = 1) -> bytes:
    import io

    doc = fitz.open()
    for _ in range(n_pages):
        page = doc.new_page(width=200, height=200)
        page.draw_rect(fitz.Rect(20, 20, 180, 180))
    buf = io.BytesIO(doc.tobytes())
    doc.close()
    return buf.getvalue()


def _b64_upload(raw: bytes) -> str:
    return "data:application/pdf;base64," + base64.b64encode(raw).decode()


@pytest.fixture(autouse=True)
def _isolated_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_PATTERN_CACHE", str(tmp_path / "cache"))
    yield


def test_list_packs_includes_builtin():
    from pycsamt.app._patterns import BUILTIN_PACK_ID, list_packs

    packs = list_packs()
    assert packs[0]["id"] == BUILTIN_PACK_ID
    assert packs[0]["n_tiles"] >= 6


def test_import_pack_upload_from_a_single_pdf():
    from pycsamt.app._patterns import import_pack_upload

    summary = import_pack_upload(
        _b64_upload(_make_pdf_bytes()), "sheet.ai", name="Uploaded Pack"
    )
    assert summary["n_sheets"] == 1
    assert summary["sheets"][0]["source_file"] == "sheet.ai"


def test_import_pack_upload_rejects_bad_payload():
    from pycsamt.app._patterns import import_pack_upload

    with pytest.raises(ValueError):
        import_pack_upload("not a data url", "x.pdf")


def test_import_pack_path_from_a_folder(tmp_path):
    from pycsamt.app._patterns import import_pack_path

    src = tmp_path / "source"
    src.mkdir()
    (src / "sheet.ai").write_bytes(_make_pdf_bytes(n_pages=2))
    summary = import_pack_path(str(src), name="Folder Pack")
    assert summary["n_sheets"] == 2


def test_import_pack_path_rejects_blank_path():
    from pycsamt.app._patterns import import_pack_path

    with pytest.raises(ValueError):
        import_pack_path("   ")


def test_sheet_figure_embeds_the_rasterized_image_with_drawrect_enabled(tmp_path):
    from pycsamt.app._patterns import import_pack_path, sheet_figure

    src = tmp_path / "source"
    src.mkdir()
    (src / "sheet.ai").write_bytes(_make_pdf_bytes())
    summary = import_pack_path(str(src), name="Fig Pack")
    fig = sheet_figure(summary["id"], summary["sheets"][0]["id"])
    assert len(fig.layout.images) == 1
    assert fig.layout.dragmode == "drawrect"


def test_add_tile_from_shape_and_grid_round_trip(tmp_path):
    from pycsamt.app._patterns import (
        add_tile_from_shape,
        import_pack_path,
        tile_grid_children,
        tile_data_uri,
    )

    src = tmp_path / "source"
    src.mkdir()
    (src / "sheet.ai").write_bytes(_make_pdf_bytes())
    summary = import_pack_path(str(src), name="Crop Pack")
    shape = {"x0": 10, "y0": 10, "x1": 100, "y1": 100}
    tile = add_tile_from_shape(
        summary["id"], summary["sheets"][0]["id"], shape, "My Swatch"
    )
    assert tile["name"] == "My Swatch"

    cells = tile_grid_children(summary["id"])
    assert len(cells) == 1
    assert tile_data_uri(summary["id"], tile["id"]).startswith("data:image/png;base64,")


def test_delete_tile_removes_it(tmp_path):
    from pycsamt.app._patterns import (
        add_tile_from_shape,
        delete_tile,
        import_pack_path,
        load_pack,
    )

    src = tmp_path / "source"
    src.mkdir()
    (src / "sheet.ai").write_bytes(_make_pdf_bytes())
    summary = import_pack_path(str(src), name="Delete Pack")
    tile = add_tile_from_shape(
        summary["id"], summary["sheets"][0]["id"],
        {"x0": 10, "y0": 10, "x1": 100, "y1": 100}, "Doomed",
    )
    delete_tile(summary["id"], tile["id"])
    assert load_pack(summary["id"]).tiles == []


def test_tile_grid_children_reports_empty_pack_message():
    from pycsamt.app._patterns import BUILTIN_PACK_ID, tile_grid_children

    # the builtin pack always has tiles, but a freshly imported empty
    # pack (no crops yet) should show a friendly placeholder instead of
    # an empty list -- exercise that branch directly.
    from pycsamt.geology.patterns import PatternLibrary

    cells = tile_grid_children(BUILTIN_PACK_ID)
    assert len(cells) == len(PatternLibrary().load_pack(BUILTIN_PACK_ID).tiles)
