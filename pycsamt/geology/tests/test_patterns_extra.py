# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Additional pycsamt.geology.patterns coverage: cache-dir resolution,
corrupt/missing manifests, PNG/zip import branches, tile deletion
edge cases, and the small standalone helpers not exercised by
test_patterns.py's end-to-end PDF round trip."""

from __future__ import annotations

import io
import json
import sys
import zipfile
from pathlib import Path

import pytest

from pycsamt.geology.patterns import (
    BUILTIN_PACK_ID,
    PatternLibrary,
    PatternTile,
    _builtin_tile_png_bytes,
    _now_slug,
    _png_size,
    _rasterize_pdf,
    tile_stencil_array,
)

fitz = pytest.importorskip("fitz")


def _library(tmp_path) -> PatternLibrary:
    return PatternLibrary(cache_dir=tmp_path / "cache")


def _make_pdf(path, n_pages: int = 1) -> None:
    doc = fitz.open()
    for _ in range(n_pages):
        page = doc.new_page(width=200, height=200)
        page.draw_rect(fitz.Rect(20, 20, 180, 180))
    doc.save(str(path))
    doc.close()


def _make_png_bytes(size=(16, 16)) -> bytes:
    from PIL import Image

    buf = io.BytesIO()
    Image.new("RGBA", size, (10, 20, 30, 255)).save(buf, format="PNG")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# cache-dir resolution
# ---------------------------------------------------------------------------


def test_default_cache_dir_falls_back_to_home_when_no_arg_or_env(monkeypatch):
    monkeypatch.delenv("PYCSAMT_PATTERN_CACHE", raising=False)
    lib = PatternLibrary()
    assert lib.cache_dir == Path.home() / ".pycsamt" / "geology_patterns"


# ---------------------------------------------------------------------------
# list_packs: entries without a manifest, or with a corrupt one
# ---------------------------------------------------------------------------


def test_list_packs_skips_entry_without_manifest(tmp_path):
    lib = _library(tmp_path)
    (lib.cache_dir / "no-manifest-here").mkdir(parents=True)
    packs = lib.list_packs()
    assert [p["id"] for p in packs] == [BUILTIN_PACK_ID]


def test_list_packs_skips_entry_with_corrupt_manifest_json(tmp_path):
    lib = _library(tmp_path)
    bad = lib.cache_dir / "corrupt-pack"
    bad.mkdir(parents=True)
    (bad / "library.json").write_text("{not valid json", encoding="utf-8")
    packs = lib.list_packs()
    assert [p["id"] for p in packs] == [BUILTIN_PACK_ID]


# ---------------------------------------------------------------------------
# import_pack: PNG sheet, path-based zip, raw-bytes non-zip
# ---------------------------------------------------------------------------


def test_import_pack_from_single_png_file(tmp_path):
    png_path = tmp_path / "swatch.png"
    png_path.write_bytes(_make_png_bytes())
    lib = _library(tmp_path)
    pack = lib.import_pack(png_path, name="PNG Pack")
    assert len(pack.sheets) == 1
    assert pack.sheets[0].source_file == "swatch.png"
    assert pack.sheets[0].width == 16


def test_import_pack_from_zip_path_not_bytes(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    zip_path = tmp_path / "pack.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(src / "sheet.ai", "sheet.ai")

    lib = _library(tmp_path)
    pack = lib.import_pack(zip_path, name="Zip Path Pack")
    assert len(pack.sheets) == 1


def test_import_pack_from_raw_png_bytes_not_zip(tmp_path):
    lib = _library(tmp_path)
    pack = lib.import_pack(
        _make_png_bytes(), name="Raw Bytes Pack", filename="raw.png"
    )
    assert len(pack.sheets) == 1


def test_import_pack_zip_skips_directory_and_unsupported_members(tmp_path):
    zip_path = tmp_path / "mixed.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.writestr("a_dir/", "")  # directory entry
        zf.writestr("readme.txt", "not a pattern sheet")  # unsupported
        zf.writestr("good.png", _make_png_bytes())  # supported

    lib = _library(tmp_path)
    pack = lib.import_pack(zip_path, name="Mixed Zip")
    assert len(pack.sheets) == 1
    assert pack.sheets[0].source_file == "good.png"


def test_import_pack_slug_falls_back_to_timestamp_when_name_has_no_letters(
    tmp_path,
):
    png_path = tmp_path / "swatch.png"
    png_path.write_bytes(_make_png_bytes())
    lib = _library(tmp_path)
    pack = lib.import_pack(png_path, name="!!!")
    assert pack.id.startswith("pack-")


# ---------------------------------------------------------------------------
# delete_tile edge cases
# ---------------------------------------------------------------------------


def test_delete_tile_unknown_id_is_a_noop(tmp_path):
    png_path = tmp_path / "swatch.png"
    png_path.write_bytes(_make_png_bytes())
    lib = _library(tmp_path)
    pack = lib.import_pack(png_path, name="Deletable")
    lib.delete_tile(pack.id, "tile-does-not-exist")  # must not raise


def test_delete_tile_with_falsy_file_skips_unlink(tmp_path):
    png_path = tmp_path / "swatch.png"
    png_path.write_bytes(_make_png_bytes())
    lib = _library(tmp_path)
    pack = lib.import_pack(png_path, name="Ghost Tile")
    pack.tiles.append(PatternTile(id="tile-0001", name="ghost", file=""))
    lib._write_manifest(pack)

    lib.delete_tile(pack.id, "tile-0001")  # must not attempt to unlink ""
    reloaded = lib.load_pack(pack.id)
    assert reloaded.tile("tile-0001") is None


# ---------------------------------------------------------------------------
# built-in materialization cache hit
# ---------------------------------------------------------------------------


def test_materialize_builtin_reuses_existing_manifest_on_second_call(tmp_path):
    lib = _library(tmp_path)
    first = lib.load_pack(BUILTIN_PACK_ID)
    second = lib.load_pack(BUILTIN_PACK_ID)
    assert [t.id for t in first.tiles] == [t.id for t in second.tiles]


def test_materialize_builtin_remakes_when_stale_manifest_has_fewer_tiles(
    tmp_path,
):
    lib = _library(tmp_path)
    directory = lib.cache_dir / BUILTIN_PACK_ID
    directory.mkdir(parents=True)
    (directory / "library.json").write_text(
        json.dumps(
            {
                "pack_id": BUILTIN_PACK_ID,
                "name": "Built-in patterns",
                "source": "",
                "imported_at": "2020-01-01T00:00:00Z",
                "sheets": [],
                "tiles": [],  # stale: fewer than the real builtin tile count
            }
        ),
        encoding="utf-8",
    )
    pack = lib.load_pack(BUILTIN_PACK_ID)
    assert len(pack.tiles) > 0


# ---------------------------------------------------------------------------
# small standalone helpers
# ---------------------------------------------------------------------------


def test_builtin_tile_png_bytes_unknown_spec_still_returns_blank_png():
    data = _builtin_tile_png_bytes("not-a-known-spec")
    assert data[:8] == b"\x89PNG\r\n\x1a\n"


def test_now_slug_is_a_compact_digit_timestamp():
    slug = _now_slug()
    assert slug.isdigit()
    assert len(slug) == 14


def test_png_size_reads_dimensions(tmp_path):
    path = tmp_path / "sized.png"
    path.write_bytes(_make_png_bytes(size=(32, 24)))
    assert _png_size(path) == (32, 24)


def test_tile_stencil_array_shape_and_range(tmp_path):
    path = tmp_path / "stencil.png"
    path.write_bytes(_make_png_bytes(size=(8, 8)))
    stencil = tile_stencil_array(path)
    assert stencil.shape == (8, 8)
    assert stencil.min() >= 0.0
    assert stencil.max() <= 1.0


def test_rasterize_pdf_raises_import_error_when_pymupdf_unavailable(
    monkeypatch, tmp_path
):
    monkeypatch.setitem(sys.modules, "pymupdf", None)
    monkeypatch.setitem(sys.modules, "fitz", None)
    pdf_path = tmp_path / "whatever.pdf"
    pdf_path.write_bytes(b"%PDF-fake")
    with pytest.raises(ImportError, match="requires PyMuPDF"):
        _rasterize_pdf(pdf_path)
