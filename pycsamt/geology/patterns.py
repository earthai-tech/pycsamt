# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Pattern/swatch packs for the Interpretation legend's ``pattern_id`` /
``pattern_source`` fields (see :class:`pycsamt.geology.lithology.RockEntry`).

pycsamt never vendors a third-party geologic pattern set (e.g. the USGS
FGDC "Digital Cartographic Standard for Geologic Map Symbolization" AI
pack) -- only the machinery to **import one the user already downloaded**
and use it. :class:`PatternLibrary` copies whatever is handed to it
(a folder, a ``.zip``, or a single PDF/AI sheet) into a small persistent
cache so the app keeps working after the original download is deleted,
mirroring the fetch/cache convention
:class:`pycsamt.geology.rock_providers.RemoteRockPropertyProvider`
already uses for the rock database (``~/.pycsamt/geology_patterns``,
override with ``$PYCSAMT_PATTERN_CACHE``).

A named, third-party pattern *sheet* (a page of many drawn symbols, as
FGDC-style packs ship them) carries no machine-readable per-symbol crop
box or label -- that information exists only as a human reading the
page. So this module deliberately does not attempt to auto-slice one:
:meth:`PatternLibrary.import_pack` rasterizes each PDF/AI page to a
reference image (via `PyMuPDF <https://pymupdf.readthedocs.io/>`_, an
optional dependency -- install the ``patterns`` extra), and
:meth:`PatternLibrary.add_tile_from_crop` lets the application crop a
named rectangle out of that image by hand, the same way a person would
pick a swatch out of a printed legend sheet. :func:`builtin_pack` ships
a small procedurally generated pattern set (diagonal lines, cross-hatch,
dots, brick, ...) that needs no import at all.
"""

from __future__ import annotations

import io
import json
import os
import re
import zipfile
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Union

__all__ = [
    "PatternPack",
    "PatternSheet",
    "PatternTile",
    "PatternLibrary",
    "PatternImportError",
    "builtin_pack",
    "tile_stencil_array",
]

PathLike = Union[str, Path]

BUILTIN_PACK_ID = "builtin"
_RASTER_SUFFIXES = {".pdf", ".ai"}
_MANIFEST_NAME = "library.json"
_DPI = 200


class PatternImportError(ValueError):
    """Raised when a pattern-pack source cannot be imported."""


@dataclass
class PatternSheet:
    """One rasterized reference page (a source PDF/AI page rendered to
    PNG), the surface :meth:`PatternLibrary.add_tile_from_crop` crops
    named swatches out of."""

    id: str
    file: str  # path relative to the pack directory
    source_file: str
    page: int
    width: int
    height: int


@dataclass
class PatternTile:
    """One named, usable swatch -- a small PNG, optionally a plain solid
    colour instead of an image."""

    id: str
    name: str
    file: str  # path relative to the pack directory ("" for a solid colour)
    color: str | None = None
    sheet_id: str | None = None
    bbox: tuple[int, int, int, int] | None = None


@dataclass
class PatternPack:
    """A named collection of :class:`PatternSheet` / :class:`PatternTile`,
    persisted under one directory in the pattern cache."""

    id: str
    name: str
    source: str
    imported_at: str
    directory: Path
    sheets: list[PatternSheet] = field(default_factory=list)
    tiles: list[PatternTile] = field(default_factory=list)

    def tile(self, tile_id: str) -> PatternTile | None:
        return next((t for t in self.tiles if t.id == tile_id), None)

    def sheet(self, sheet_id: str) -> PatternSheet | None:
        return next((s for s in self.sheets if s.id == sheet_id), None)


# ---------------------------------------------------------------------------
# PatternLibrary
# ---------------------------------------------------------------------------


class PatternLibrary:
    """Manages persisted pattern packs under one cache directory."""

    def __init__(self, cache_dir: PathLike | None = None) -> None:
        self.cache_dir = self._resolve_cache_dir(cache_dir)

    @staticmethod
    def _resolve_cache_dir(cache_dir: PathLike | None) -> Path:
        if cache_dir is not None:
            return Path(cache_dir)
        env = os.environ.get("PYCSAMT_PATTERN_CACHE")
        if env:
            return Path(env)
        return Path.home() / ".pycsamt" / "geology_patterns"

    # ------------------------------------------------------------------
    # discovery
    # ------------------------------------------------------------------

    def list_packs(self) -> list[dict[str, Any]]:
        """Return ``{id, name, n_sheets, n_tiles}`` for every persisted
        pack, plus the always-available built-in pack first."""
        out = [
            {
                "id": BUILTIN_PACK_ID,
                "name": "Built-in patterns",
                "n_sheets": 0,
                "n_tiles": len(builtin_pack().tiles),
            }
        ]
        if not self.cache_dir.is_dir():
            return out
        for entry in sorted(self.cache_dir.iterdir()):
            manifest = entry / _MANIFEST_NAME
            if entry.name == BUILTIN_PACK_ID or not manifest.is_file():
                continue
            try:
                data = json.loads(manifest.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            out.append(
                {
                    "id": data.get("pack_id", entry.name),
                    "name": data.get("name", entry.name),
                    "n_sheets": len(data.get("sheets", [])),
                    "n_tiles": len(data.get("tiles", [])),
                }
            )
        return out

    def load_pack(self, pack_id: str) -> PatternPack:
        """Load a persisted pack by id, materializing :data:`builtin_pack`
        into the cache on first request so every pack (built-in or
        imported) is read through the same on-disk path."""
        if pack_id == BUILTIN_PACK_ID:
            return self._materialize_builtin()
        directory = self.cache_dir / _slug(pack_id)
        manifest_path = directory / _MANIFEST_NAME
        if not manifest_path.is_file():
            raise PatternImportError(f"unknown pattern pack {pack_id!r}")
        data = json.loads(manifest_path.read_text(encoding="utf-8"))
        return _pack_from_manifest(data, directory)

    # ------------------------------------------------------------------
    # import
    # ------------------------------------------------------------------

    def import_pack(
        self,
        source: PathLike | bytes,
        *,
        name: str | None = None,
        filename: str | None = None,
    ) -> PatternPack:
        """Copy *source* into the persistent cache and rasterize any
        PDF/AI sheet it contains.

        *source* is a folder, a ``.zip`` archive, a single PDF/AI file
        path, or raw ``bytes`` (with *filename* giving its name/kind, for
        an uploaded file). Returns the resulting :class:`PatternPack`;
        raises :class:`PatternImportError` on anything unusable, and
        raises :class:`ImportError` (never silently degrades) when a
        PDF/AI sheet needs rasterizing but ``pymupdf`` is not installed.
        """
        pack_name = name or (
            Path(filename).stem if filename
            else Path(str(source)).stem if isinstance(source, (str, Path))
            else "pattern-pack"
        )
        pack_id = _slug(pack_name) or f"pack-{_now_slug()}"
        directory = self.cache_dir / pack_id
        source_dir = directory / "source"
        sheets_dir = directory / "sheets"
        source_dir.mkdir(parents=True, exist_ok=True)
        sheets_dir.mkdir(parents=True, exist_ok=True)

        raster_files = list(
            self._stage_source(source, filename, source_dir)
        )
        if not raster_files:
            raise PatternImportError(
                "no .pdf/.ai/.png/.jpg pattern sheet found in this source"
            )

        sheets: list[PatternSheet] = []
        n = 0
        for staged in raster_files:
            if staged.suffix.lower() in _RASTER_SUFFIXES:
                for page_png, page_no, w, h in _rasterize_pdf(staged):
                    n += 1
                    sheet_id = f"sheet-{n:04d}"
                    out_path = sheets_dir / f"{sheet_id}.png"
                    out_path.write_bytes(page_png)
                    sheets.append(
                        PatternSheet(
                            id=sheet_id, file=f"sheets/{sheet_id}.png",
                            source_file=staged.name, page=page_no,
                            width=w, height=h,
                        )
                    )
            elif staged.suffix.lower() in {".png", ".jpg", ".jpeg"}:
                n += 1
                sheet_id = f"sheet-{n:04d}"
                out_path = sheets_dir / f"{sheet_id}.png"
                out_path.write_bytes(staged.read_bytes())
                w, h = _png_size(out_path)
                sheets.append(
                    PatternSheet(
                        id=sheet_id, file=f"sheets/{sheet_id}.png",
                        source_file=staged.name, page=1, width=w, height=h,
                    )
                )

        pack = PatternPack(
            id=pack_id, name=pack_name, source=filename or str(source),
            imported_at=_now(), directory=directory, sheets=sheets, tiles=[],
        )
        self._write_manifest(pack)
        return pack

    def _stage_source(
        self,
        source: PathLike | bytes,
        filename: str | None,
        source_dir: Path,
    ) -> list[Path]:
        """Copy *source* under *source_dir* and return the staged
        PDF/AI/PNG/JPG files worth rasterizing."""
        if isinstance(source, bytes):
            name = filename or "upload"
            if name.lower().endswith(".zip") or source[:2] == b"PK":
                return self._extract_zip_bytes(source, source_dir)
            target = source_dir / Path(name).name
            target.write_bytes(source)
            return [target] if target.suffix.lower() in (
                _RASTER_SUFFIXES | {".png", ".jpg", ".jpeg"}
            ) else []

        path = Path(source)
        if not path.exists():
            raise PatternImportError(f"source not found: {path}")
        if path.is_dir():
            staged = []
            for p in sorted(path.rglob("*")):
                if p.is_file() and p.suffix.lower() in (
                    _RASTER_SUFFIXES | {".png", ".jpg", ".jpeg"}
                ):
                    target = source_dir / p.relative_to(path)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(p.read_bytes())
                    staged.append(target)
            return staged
        if path.suffix.lower() == ".zip":
            return self._extract_zip_bytes(path.read_bytes(), source_dir)
        target = source_dir / path.name
        target.write_bytes(path.read_bytes())
        return [target] if target.suffix.lower() in (
            _RASTER_SUFFIXES | {".png", ".jpg", ".jpeg"}
        ) else []

    @staticmethod
    def _extract_zip_bytes(raw: bytes, source_dir: Path) -> list[Path]:
        staged: list[Path] = []
        try:
            with zipfile.ZipFile(io.BytesIO(raw)) as zf:
                for info in zf.infolist():
                    if info.is_dir():
                        continue
                    name = Path(info.filename)
                    if name.suffix.lower() not in (
                        _RASTER_SUFFIXES | {".png", ".jpg", ".jpeg"}
                    ):
                        continue
                    # Flatten into source_dir, guarding against a
                    # zip-slip path escaping it.
                    safe_name = _safe_member_name(name)
                    target = source_dir / safe_name
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(zf.read(info))
                    staged.append(target)
        except zipfile.BadZipFile as error:
            raise PatternImportError("not a valid .zip archive") from error
        return staged

    # ------------------------------------------------------------------
    # tiles (named crops)
    # ------------------------------------------------------------------

    def add_tile_from_crop(
        self,
        pack_id: str,
        sheet_id: str,
        bbox: tuple[float, float, float, float],
        name: str,
    ) -> PatternTile:
        """Crop *bbox* (pixel ``x0, y0, x1, y1``) out of *sheet_id* and
        save it as a new named tile in the pack."""
        try:
            from PIL import Image
        except ImportError as error:  # pragma: no cover - Pillow is core
            raise ImportError("Pillow is required to crop a pattern sheet") from error

        if pack_id == BUILTIN_PACK_ID:
            raise PatternImportError("the built-in pack cannot be edited")
        pack = self.load_pack(pack_id)
        sheet = pack.sheet(sheet_id)
        if sheet is None:
            raise PatternImportError(f"unknown sheet {sheet_id!r}")
        if not str(name).strip():
            raise PatternImportError("a swatch needs a name")

        x0, y0, x1, y1 = bbox
        x0, x1 = sorted((max(0, int(x0)), max(0, int(x1))))
        y0, y1 = sorted((max(0, int(y0)), max(0, int(y1))))
        if x1 - x0 < 4 or y1 - y0 < 4:
            raise PatternImportError("crop is too small")

        image = Image.open(pack.directory / sheet.file)
        crop = image.crop(
            (x0, y0, min(x1, image.width), min(y1, image.height))
        )
        tile_id = f"tile-{len(pack.tiles) + 1:04d}"
        tiles_dir = pack.directory / "tiles"
        tiles_dir.mkdir(exist_ok=True)
        out_path = tiles_dir / f"{tile_id}.png"
        crop.convert("RGBA").save(out_path)

        tile = PatternTile(
            id=tile_id, name=str(name).strip(),
            file=f"tiles/{tile_id}.png", sheet_id=sheet_id,
            bbox=(x0, y0, x1, y1),
        )
        pack.tiles.append(tile)
        self._write_manifest(pack)
        return tile

    def delete_tile(self, pack_id: str, tile_id: str) -> None:
        if pack_id == BUILTIN_PACK_ID:
            raise PatternImportError("the built-in pack cannot be edited")
        pack = self.load_pack(pack_id)
        tile = pack.tile(tile_id)
        if tile is None:
            return
        pack.tiles = [t for t in pack.tiles if t.id != tile_id]
        if tile.file:
            (pack.directory / tile.file).unlink(missing_ok=True)
        self._write_manifest(pack)

    # ------------------------------------------------------------------
    # helpers
    # ------------------------------------------------------------------

    def _write_manifest(self, pack: PatternPack) -> None:
        pack.directory.mkdir(parents=True, exist_ok=True)
        payload = {
            "pack_id": pack.id,
            "name": pack.name,
            "source": pack.source,
            "imported_at": pack.imported_at,
            "sheets": [vars(s) for s in pack.sheets],
            "tiles": [
                {**vars(t), "bbox": list(t.bbox) if t.bbox else None}
                for t in pack.tiles
            ],
        }
        (pack.directory / _MANIFEST_NAME).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8"
        )

    def _materialize_builtin(self) -> PatternPack:
        pack = builtin_pack()
        directory = self.cache_dir / BUILTIN_PACK_ID
        tiles_dir = directory / "tiles"
        tiles_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = directory / _MANIFEST_NAME
        if manifest_path.is_file():
            data = json.loads(manifest_path.read_text(encoding="utf-8"))
            if len(data.get("tiles", [])) == len(pack.tiles):
                return _pack_from_manifest(data, directory)
        for tile in pack.tiles:
            (directory / tile.file).write_bytes(
                _builtin_tile_png_bytes(tile.name)
            )
        pack.directory = directory
        self._write_manifest(pack)
        return pack


# ---------------------------------------------------------------------------
# built-in procedurally generated patterns
# ---------------------------------------------------------------------------

_BUILTIN_SPECS = (
    "diagonal-lines", "cross-hatch", "dots", "brick", "dashed-lines",
    "grid", "triangles", "stipple",
)


def builtin_pack() -> PatternPack:
    """The always-available, zero-setup pattern pack -- a small set of
    neutral procedurally generated hatch/stipple/brick/dot tiles, not
    derived from any third-party source."""
    tiles = [
        PatternTile(
            id=f"builtin-{i + 1:02d}", name=name,
            file=f"tiles/builtin-{i + 1:02d}.png",
        )
        for i, name in enumerate(_BUILTIN_SPECS)
    ]
    return PatternPack(
        id=BUILTIN_PACK_ID, name="Built-in patterns", source="",
        imported_at=_now(), directory=Path(), sheets=[], tiles=tiles,
    )


def _builtin_tile_png_bytes(spec_name: str, size: int = 64) -> bytes:
    from PIL import Image, ImageDraw

    img = Image.new("RGBA", (size, size), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)
    fg = (30, 30, 30, 255)
    if spec_name == "diagonal-lines":
        for x in range(-size, size * 2, 10):
            draw.line([(x, 0), (x + size, size)], fill=fg, width=2)
    elif spec_name == "cross-hatch":
        for x in range(-size, size * 2, 10):
            draw.line([(x, 0), (x + size, size)], fill=fg, width=1)
            draw.line([(x + size, 0), (x, size)], fill=fg, width=1)
    elif spec_name == "dots":
        for y in range(4, size, 12):
            for x in range(4, size, 12):
                draw.ellipse([x - 2, y - 2, x + 2, y + 2], fill=fg)
    elif spec_name == "brick":
        for y in range(0, size, 16):
            draw.line([(0, y), (size, y)], fill=fg, width=2)
        for i, y in enumerate(range(0, size, 16)):
            off = 0 if (i % 2 == 0) else 16
            for x in range(off, size, 32):
                draw.line([(x, y), (x, y + 16)], fill=fg, width=2)
    elif spec_name == "dashed-lines":
        for y in range(4, size, 10):
            for x in range(0, size, 12):
                draw.line([(x, y), (x + 6, y)], fill=fg, width=2)
    elif spec_name == "grid":
        for x in range(0, size, 10):
            draw.line([(x, 0), (x, size)], fill=fg, width=1)
        for y in range(0, size, 10):
            draw.line([(0, y), (size, y)], fill=fg, width=1)
    elif spec_name == "triangles":
        for y in range(0, size, 16):
            for x in range(0, size, 16):
                draw.polygon(
                    [(x, y + 14), (x + 7, y + 2), (x + 14, y + 14)],
                    outline=fg,
                )
    elif spec_name == "stipple":
        import random

        rng = random.Random(hash(spec_name) & 0xFFFF)
        for _ in range(140):
            x, y = rng.randrange(size), rng.randrange(size)
            draw.point((x, y), fill=fg)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


# ---------------------------------------------------------------------------
# PDF/AI rasterization
# ---------------------------------------------------------------------------


def _rasterize_pdf(path: Path) -> list[tuple[bytes, int, int, int]]:
    """Render every page of a PDF/AI(PDF-compatible) file to PNG bytes
    via PyMuPDF. Raises :class:`ImportError` (not a silent skip) if
    ``pymupdf`` is not installed -- rasterizing is the whole point of
    importing a sheet-based pack."""
    try:
        import pymupdf as fitz  # the modern import name; `fitz` itself
        # is a deprecated compatibility alias PyMuPDF >=1.23 still ships
        # but warns on at import time.
    except ImportError:
        try:
            import fitz  # pragma: no cover - pre-1.23 PyMuPDF only
        except ImportError as error:
            raise ImportError(
                "Importing a PDF/AI pattern sheet requires PyMuPDF -- "
                "install it with `pip install pycsamt[patterns]` "
                "(or `pip install pymupdf`)."
            ) from error

    out = []
    zoom = _DPI / 72.0
    matrix = fitz.Matrix(zoom, zoom)
    with fitz.open(path) as doc:
        for i, page in enumerate(doc, start=1):
            pix = page.get_pixmap(matrix=matrix)
            out.append((pix.tobytes("png"), i, pix.width, pix.height))
    return out


def _png_size(path: Path) -> tuple[int, int]:
    from PIL import Image

    with Image.open(path) as img:
        return img.width, img.height


def tile_stencil_array(path: PathLike):
    """Read a tile PNG as an ink-density *stencil*: an ``H x W`` array
    of alpha-weighted darkness in ``[0, 1]`` (``0`` = no ink / paper,
    ``1`` = full ink), read from :attr:`PatternTile.file`.

    The tile's own colour is discarded on purpose -- a pattern is
    applied tinted by the *legend's* assigned colour (see
    :func:`pycsamt.map.styles.pattern_band_stops`), not by whatever
    colour the source pack happened to draw its swatch in, so the same
    hatch/stipple/brick shape reads consistently whatever unit it is
    assigned to.
    """
    import numpy as np
    from PIL import Image

    with Image.open(path) as img:
        rgba = np.asarray(img.convert("RGBA"), dtype=np.float32)
    r, g, b, a = rgba[..., 0], rgba[..., 1], rgba[..., 2], rgba[..., 3]
    luminance = (0.299 * r + 0.587 * g + 0.114 * b) / 255.0
    alpha = a / 255.0
    return alpha * (1.0 - luminance)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------


def _pack_from_manifest(data: dict, directory: Path) -> PatternPack:
    sheets = [PatternSheet(**s) for s in data.get("sheets", [])]
    tiles = [
        PatternTile(
            **{**t, "bbox": tuple(t["bbox"]) if t.get("bbox") else None}
        )
        for t in data.get("tiles", [])
    ]
    return PatternPack(
        id=data.get("pack_id", directory.name),
        name=data.get("name", directory.name),
        source=data.get("source", ""),
        imported_at=data.get("imported_at", ""),
        directory=directory,
        sheets=sheets,
        tiles=tiles,
    )


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9_-]+", "-", str(text).strip().casefold()).strip("-")[:60]


def _safe_member_name(name: Path) -> Path:
    parts = [p for p in name.parts if p not in ("..", "", ".")]
    return Path(*parts) if parts else Path(name.name)


def _now() -> str:
    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def _now_slug() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%d%H%M%S")
