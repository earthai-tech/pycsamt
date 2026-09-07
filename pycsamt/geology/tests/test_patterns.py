# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Tests for pycsamt.geology.patterns — pattern/swatch pack import,
persistence, and cropping.

Uses a genuine end-to-end import/crop/reload round trip against a
synthetic PDF built with PyMuPDF itself, not mocks -- the same
"real file, real filesystem" philosophy as test_rock_providers.py.
"""

from __future__ import annotations

import shutil

import pytest

from pycsamt.geology.patterns import (
    BUILTIN_PACK_ID,
    PatternImportError,
    PatternLibrary,
    builtin_pack,
)

fitz = pytest.importorskip("fitz")


def _make_pdf(path, n_pages: int = 1) -> None:
    doc = fitz.open()
    for _ in range(n_pages):
        page = doc.new_page(width=200, height=200)
        page.draw_rect(fitz.Rect(20, 20, 180, 180))
    doc.save(str(path))
    doc.close()


def _library(tmp_path) -> PatternLibrary:
    return PatternLibrary(cache_dir=tmp_path / "cache")


# ---------------------------------------------------------------------------
# built-in pack
# ---------------------------------------------------------------------------


def test_builtin_pack_has_named_tiles_and_no_sheets():
    pack = builtin_pack()
    assert pack.id == BUILTIN_PACK_ID
    assert len(pack.tiles) >= 6
    assert pack.sheets == []
    assert all(t.name for t in pack.tiles)


def test_list_packs_always_includes_builtin(tmp_path):
    lib = _library(tmp_path)
    packs = lib.list_packs()
    assert packs[0]["id"] == BUILTIN_PACK_ID
    assert packs[0]["n_tiles"] == len(builtin_pack().tiles)


def test_load_builtin_materializes_real_tile_files(tmp_path):
    lib = _library(tmp_path)
    pack = lib.load_pack(BUILTIN_PACK_ID)
    for tile in pack.tiles:
        assert (pack.directory / tile.file).is_file()


# ---------------------------------------------------------------------------
# import from a folder / single file / zip
# ---------------------------------------------------------------------------


def test_import_pack_from_folder_rasterizes_every_pdf(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet_a.ai", n_pages=1)
    _make_pdf(src / "sheet_b.ai", n_pages=2)

    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Test Pack")
    assert pack.id == "test-pack"
    # 1 page from sheet_a + 2 pages from sheet_b
    assert len(pack.sheets) == 3
    for sheet in pack.sheets:
        assert (pack.directory / sheet.file).is_file()
        assert sheet.width > 0 and sheet.height > 0


def test_import_pack_from_single_pdf_file(tmp_path):
    pdf_path = tmp_path / "single.ai"
    _make_pdf(pdf_path, n_pages=1)

    lib = _library(tmp_path)
    pack = lib.import_pack(pdf_path, name="Single Sheet")
    assert len(pack.sheets) == 1
    assert pack.sheets[0].source_file == "single.ai"


def test_import_pack_from_zip_bytes(tmp_path):
    import zipfile

    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    zip_path = tmp_path / "pack.zip"
    with zipfile.ZipFile(zip_path, "w") as zf:
        zf.write(src / "sheet.ai", "sheet.ai")

    lib = _library(tmp_path)
    pack = lib.import_pack(
        zip_path.read_bytes(), name="Zipped Pack", filename="pack.zip"
    )
    assert len(pack.sheets) == 1


def test_import_pack_rejects_missing_source(tmp_path):
    lib = _library(tmp_path)
    with pytest.raises(PatternImportError):
        lib.import_pack(tmp_path / "does-not-exist")


def test_import_pack_rejects_source_with_no_usable_sheets(tmp_path):
    src = tmp_path / "empty_source"
    src.mkdir()
    (src / "readme.txt").write_text("nothing usable here")
    lib = _library(tmp_path)
    with pytest.raises(PatternImportError):
        lib.import_pack(src)


def test_import_pack_rejects_bad_zip(tmp_path):
    lib = _library(tmp_path)
    with pytest.raises(PatternImportError):
        lib.import_pack(b"not a zip", filename="bad.zip")


# ---------------------------------------------------------------------------
# persistence across the original source's deletion
# ---------------------------------------------------------------------------


def test_pack_survives_deletion_of_the_original_source(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")

    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Persistent Pack")
    shutil.rmtree(src)

    reloaded = lib.load_pack(pack.id)
    assert len(reloaded.sheets) == 1
    assert (reloaded.directory / reloaded.sheets[0].file).is_file()


def test_pack_survives_across_a_new_library_instance(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")

    lib1 = _library(tmp_path)
    pack = lib1.import_pack(src, name="Cross Instance")

    lib2 = _library(tmp_path)
    reloaded = lib2.load_pack(pack.id)
    assert reloaded.name == "Cross Instance"
    assert len(reloaded.sheets) == 1


# ---------------------------------------------------------------------------
# named crops (tiles)
# ---------------------------------------------------------------------------


def test_add_tile_from_crop_creates_a_named_swatch(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Crop Pack")

    tile = lib.add_tile_from_crop(
        pack.id, pack.sheets[0].id, (10, 10, 100, 100), "My Swatch"
    )
    assert tile.name == "My Swatch"
    assert (pack.directory / tile.file).is_file()

    reloaded = lib.load_pack(pack.id)
    assert len(reloaded.tiles) == 1
    assert reloaded.tiles[0].bbox == (10, 10, 100, 100)


def test_add_tile_from_crop_rejects_tiny_crop(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Tiny Crop Pack")
    with pytest.raises(PatternImportError):
        lib.add_tile_from_crop(pack.id, pack.sheets[0].id, (0, 0, 1, 1), "x")


def test_add_tile_from_crop_rejects_blank_name(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Blank Name Pack")
    with pytest.raises(PatternImportError):
        lib.add_tile_from_crop(
            pack.id, pack.sheets[0].id, (10, 10, 100, 100), "  "
        )


def test_add_tile_from_crop_rejects_builtin_pack(tmp_path):
    lib = _library(tmp_path)
    with pytest.raises(PatternImportError):
        lib.add_tile_from_crop(
            BUILTIN_PACK_ID, "sheet-0001", (0, 0, 50, 50), "x"
        )


def test_delete_tile_removes_it_and_its_file(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    lib = _library(tmp_path)
    pack = lib.import_pack(src, name="Delete Pack")
    tile = lib.add_tile_from_crop(
        pack.id, pack.sheets[0].id, (10, 10, 100, 100), "Doomed"
    )
    file_path = pack.directory / tile.file
    assert file_path.is_file()

    lib.delete_tile(pack.id, tile.id)
    assert not file_path.is_file()
    reloaded = lib.load_pack(pack.id)
    assert reloaded.tiles == []


def test_list_packs_reports_imported_packs_alongside_builtin(tmp_path):
    src = tmp_path / "source"
    src.mkdir()
    _make_pdf(src / "sheet.ai")
    lib = _library(tmp_path)
    lib.import_pack(src, name="Listed Pack")

    ids = {p["id"] for p in lib.list_packs()}
    assert BUILTIN_PACK_ID in ids
    assert "listed-pack" in ids


def test_env_var_overrides_default_cache_dir(tmp_path, monkeypatch):
    monkeypatch.setenv("PYCSAMT_PATTERN_CACHE", str(tmp_path / "env-cache"))
    lib = PatternLibrary()
    assert lib.cache_dir == tmp_path / "env-cache"
