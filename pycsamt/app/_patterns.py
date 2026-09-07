# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared pattern-pack upload, browsing, cropping, and Plotly adapters
for pyCSAMT applications.

Bridges :mod:`pycsamt.geology.patterns` (the persistent
:class:`~pycsamt.geology.patterns.PatternLibrary`) to Map View's Dash
layer: base64 upload decoding, a rasterized-sheet figure with Plotly's
draw-rectangle tool enabled (so a swatch can be cropped by hand out of a
third-party pattern sheet the same way a person would pick one off a
printed legend), a browsable swatch grid, and data-URI helpers for
displaying tile/sheet PNGs without writing them anywhere web-servable.
"""

from __future__ import annotations

import base64
import binascii
from typing import Any

from pycsamt.geology.patterns import (
    BUILTIN_PACK_ID,
    PatternLibrary,
    PatternPack,
)

__all__ = [
    "BUILTIN_PACK_ID",
    "list_packs",
    "load_pack",
    "import_pack_upload",
    "import_pack_path",
    "add_tile_from_shape",
    "delete_tile",
    "sheet_data_uri",
    "tile_data_uri",
    "sheet_figure",
    "tile_grid_children",
    "tile_stencil_array",
]

_MAX_UPLOAD_BYTES = 60 * 1024 * 1024  # a full FGDC-style AI pack is tens of MB


def _library() -> PatternLibrary:
    return PatternLibrary()


def list_packs() -> list[dict[str, Any]]:
    """``[{id, name, n_sheets, n_tiles}, ...]``, built-in pack first."""
    return _library().list_packs()


def load_pack(pack_id: str) -> PatternPack:
    return _library().load_pack(pack_id)


def _decode_b64(contents: str) -> tuple[bytes, str]:
    if not isinstance(contents, str) or "," not in contents:
        raise ValueError("invalid upload payload")
    header, encoded = contents.split(",", 1)
    if ";base64" not in header.lower():
        raise ValueError("upload must be base64 encoded")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("invalid base64 upload") from error
    if len(raw) > _MAX_UPLOAD_BYTES:
        raise ValueError(f"upload exceeds {_MAX_UPLOAD_BYTES} bytes")
    return raw, header


def import_pack_upload(
    contents: str, filename: str | None = None, *, name: str | None = None
) -> dict[str, Any]:
    """Decode a browser upload (a ``.zip`` pack or a single ``.pdf``/
    ``.ai`` sheet) and import it into the persistent pattern library."""
    raw, _header = _decode_b64(contents)
    pack = _library().import_pack(raw, name=name, filename=filename)
    return _pack_summary(pack)


def import_pack_path(path: str, *, name: str | None = None) -> dict[str, Any]:
    """Import a pack from a local folder, ``.zip``, or single sheet path
    -- the mapview server runs locally, so a path the user can browse to
    on disk is a legitimate alternative to a (possibly very large)
    browser upload."""
    if not str(path or "").strip():
        raise ValueError("enter a folder, .zip, or .pdf/.ai file path")
    pack = _library().import_pack(path, name=name)
    return _pack_summary(pack)


def add_tile_from_shape(
    pack_id: str, sheet_id: str, shape: dict[str, Any], name: str
) -> dict[str, Any]:
    """Crop a swatch out of *sheet_id* using a Plotly draw-rectangle
    shape (``{"x0":..., "y0":..., "x1":..., "y1":...}`` from
    ``relayoutData``/``figure.layout.shapes`` -- image traces map data
    coordinates 1:1 to pixel rows/columns, so no axis flip is needed)."""
    x0, y0, x1, y1 = shape["x0"], shape["y0"], shape["x1"], shape["y1"]
    lib = _library()
    tile = lib.add_tile_from_crop(
        pack_id, sheet_id, (x0, y0, x1, y1), name
    )
    return {
        "id": tile.id, "name": tile.name, "file": tile.file,
        "sheet_id": tile.sheet_id, "bbox": tile.bbox,
    }


def delete_tile(pack_id: str, tile_id: str) -> None:
    _library().delete_tile(pack_id, tile_id)


# ---------------------------------------------------------------------------
# image data URIs
# ---------------------------------------------------------------------------


def sheet_data_uri(pack_id: str, sheet_id: str) -> str:
    pack = load_pack(pack_id)
    sheet = pack.sheet(sheet_id)
    if sheet is None:
        raise ValueError(f"unknown sheet {sheet_id!r}")
    raw = (pack.directory / sheet.file).read_bytes()
    return "data:image/png;base64," + base64.b64encode(raw).decode()


def tile_data_uri(pack_id: str, tile_id: str) -> str:
    pack = load_pack(pack_id)
    tile = pack.tile(tile_id)
    if tile is None or not tile.file:
        raise ValueError(f"unknown tile {tile_id!r}")
    raw = (pack.directory / tile.file).read_bytes()
    return "data:image/png;base64," + base64.b64encode(raw).decode()


def tile_stencil_array(pack_id: str, tile_id: str):
    """``H x W`` ink-density array (``[0, 1]``) for *tile_id* -- what
    :func:`pycsamt.map.volume` actually samples to texture-fill a
    geology band with this swatch. See
    :func:`pycsamt.geology.patterns.tile_stencil_array`."""
    from pycsamt.geology.patterns import tile_stencil_array as _stencil

    pack = load_pack(pack_id)
    tile = pack.tile(tile_id)
    if tile is None or not tile.file:
        raise ValueError(f"unknown tile {tile_id!r}")
    return _stencil(pack.directory / tile.file)


# ---------------------------------------------------------------------------
# figures / grid
# ---------------------------------------------------------------------------


def sheet_figure(pack_id: str, sheet_id: str):
    """A crop-ready figure for *sheet_id*: the rasterized sheet as a
    background image, with Plotly's rectangle-draw tool the caller wires
    up via ``dcc.Graph(config={"modeBarButtonsToAdd": ["drawrect"]})``.
    """
    import plotly.graph_objects as go

    pack = load_pack(pack_id)
    sheet = pack.sheet(sheet_id)
    if sheet is None:
        raise ValueError(f"unknown sheet {sheet_id!r}")
    uri = sheet_data_uri(pack_id, sheet_id)
    fig = go.Figure()
    fig.add_layout_image(
        dict(
            source=uri, xref="x", yref="y", x=0, y=0,
            sizex=sheet.width, sizey=sheet.height,
            sizing="stretch", layer="below",
        )
    )
    fig.update_xaxes(
        range=[0, sheet.width], visible=False, constrain="domain",
    )
    fig.update_yaxes(
        range=[sheet.height, 0], visible=False,
        scaleanchor="x", scaleratio=1,
    )
    fig.update_layout(
        dragmode="drawrect",
        newshape=dict(line_color="#ef4444", fillcolor="rgba(239,68,68,0.15)"),
        margin=dict(l=0, r=0, t=0, b=0),
        height=420,
    )
    return fig


def tile_grid_children(pack_id: str, *, selected_tile_id: str | None = None):
    """Return the Illustrator-swatch-panel-style grid: one tile per
    clickable cell (name, thumbnail, resistivity-agnostic id) using
    Dash pattern-matching ids ``{"type": "geo-swatch-btn", "pack":
    pack_id, "tile": tile.id}`` so the caller wires a single
    ``Input(..., ALL)`` callback regardless of how many tiles exist."""
    from dash import html

    pack = load_pack(pack_id)
    if not pack.tiles:
        return [html.Div("No swatches in this pack yet", className="text-muted")]
    cells = []
    for tile in pack.tiles:
        is_selected = tile.id == selected_tile_id
        cells.append(
            html.Button(
                [
                    html.Img(
                        src=tile_data_uri(pack_id, tile.id),
                        style={
                            "width": "48px", "height": "48px",
                            "objectFit": "cover", "borderRadius": "3px",
                            "display": "block", "margin": "0 auto",
                        },
                    ),
                    html.Div(
                        tile.name,
                        style={
                            "fontSize": "9.5px", "textAlign": "center",
                            "marginTop": "2px", "whiteSpace": "nowrap",
                            "overflow": "hidden", "textOverflow": "ellipsis",
                            "maxWidth": "60px",
                        },
                    ),
                ],
                id={"type": "geo-swatch-btn", "pack": pack_id, "tile": tile.id},
                n_clicks=0,
                title=tile.name,
                className="mv-swatch-cell" + (
                    " mv-swatch-cell-selected" if is_selected else ""
                ),
                style={
                    "border": "2px solid #3b82f6" if is_selected else
                    "1px solid var(--mv-border, #ccc)",
                    "borderRadius": "5px", "padding": "4px",
                    "background": "none", "cursor": "pointer",
                },
            )
        )
    return cells


def _pack_summary(pack: PatternPack) -> dict[str, Any]:
    return {
        "id": pack.id,
        "name": pack.name,
        "n_sheets": len(pack.sheets),
        "n_tiles": len(pack.tiles),
        "sheets": [
            {"id": s.id, "source_file": s.source_file, "page": s.page}
            for s in pack.sheets
        ],
    }
