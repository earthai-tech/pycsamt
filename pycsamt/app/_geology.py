# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCGL (Interpretation legend) upload, auto-suggest, and Plotly
render-band adapters for pyCSAMT applications.

Mirrors :mod:`pycsamt.app._borehole` / :mod:`pycsamt.app._points`: this
module owns the bridge between the versioned PCGL document
(:mod:`pycsamt.format.geology`) and Map View's Dash layer -- a JSON-safe
application store, a flat-table Studio draft (PCGL is already one flat
table, so no multi-table bridge like PCBH's is needed), and the plain
``(rho_min, rho_max, hex_color)`` band list
:func:`pycsamt.map.volume.build_3d_map` consumes via
``VolumeMapOptions.geology``.
"""

from __future__ import annotations

import base64
import binascii
import json
from typing import Any

import numpy as np

from pycsamt.format.geology import (
    GeologyLegend,
    GeologyLegendValidationError,
    legend_from_dict,
    legend_to_dict,
)
from pycsamt.geology.lithology import RockDatabase, RockEntry

__all__ = [
    "TABLE_COLUMNS",
    "decode_geology_upload",
    "legend_from_store",
    "store_from_legend",
    "geology_bands_from_store",
    "table_rows_from_legend",
    "legend_from_table_rows",
    "auto_suggest_legend",
    "legend_preview_figure",
    "borehole_lithology_colors",
    "sync_legend_colors_with_boreholes",
    "geology_pattern_stencils_from_store",
]

_MAX_BYTES = 2 * 1024 * 1024

# Columns shown/edited in the Interpretation Studio's Legend table, in
# display order. ``code`` is deliberately omitted -- it is an internal
# LAS-export id, assigned automatically, not something a user curates.
TABLE_COLUMNS = (
    "name", "rho_min", "rho_max", "color",
    "pattern_id", "pattern_source", "description", "source",
)


def decode_geology_upload(
    contents: str, filename: str | None = None
) -> dict[str, Any]:
    """Validate a Dash data URL and return a JSON-safe PCGL store.

    Accepts a canonical ``.pcgl.json`` document or a plain CSV with the
    columns :meth:`pycsamt.geology.lithology.RockDatabase.from_csv`
    understands (only ``name, rho_min, rho_max`` required).
    """
    name = filename or "legend.pcgl.json"
    lower = name.lower()
    if not isinstance(contents, str) or "," not in contents:
        raise ValueError("invalid PCGL upload payload")
    header, encoded = contents.split(",", 1)
    if ";base64" not in header.lower():
        raise ValueError("PCGL upload must be base64 encoded")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("invalid base64 PCGL upload") from error
    if len(raw) > _MAX_BYTES:
        raise ValueError(f"PCGL upload exceeds {_MAX_BYTES} bytes")

    if lower.endswith((".csv", ".tsv", ".txt")):
        import tempfile
        from pathlib import Path

        folder = tempfile.mkdtemp(prefix="pcgl-")
        path = Path(folder) / name
        path.write_bytes(raw)
        legend = GeologyLegend.from_csv(path, title=Path(name).stem)
    else:
        try:
            payload = json.loads(raw.decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as error:
            raise ValueError("PCGL upload is not valid UTF-8 JSON") from error
        legend = legend_from_dict(payload)
    return store_from_legend(legend, filename=name)


def store_from_legend(
    legend: GeologyLegend, *, filename: str = "legend.pcgl.json"
) -> dict[str, Any]:
    """Build the JSON-safe application store for *legend*."""
    return {
        "filename": filename,
        "document": legend_to_dict(legend),
        "n_entries": len(legend.entries),
    }


def legend_from_store(store: dict[str, Any] | None) -> GeologyLegend | None:
    """Restore a validated :class:`GeologyLegend` from an application
    store (``{"document": <pcgl dict>}`` or the raw pcgl dict itself)."""
    if not store:
        return None
    payload = store.get("document", store)
    if not isinstance(payload, dict):
        raise ValueError("invalid PCGL application store")
    return legend_from_dict(payload)


def geology_bands_from_store(
    store: dict[str, Any] | GeologyLegend | None,
) -> tuple[tuple[float, float, str, str], ...] | None:
    """Return the ``(rho_min, rho_max, color, name)`` band list
    ``VolumeMapOptions.geology`` expects, or ``None`` when *store* carries
    no usable legend -- the caller (Map View's ``render()``) leaves the
    continuous colour ramp untouched in that case.

    The 4th element (the entry's name) lets
    :func:`pycsamt.map.styles.geology_colorbar_ticks` replace the
    colourbar's numeric Ω·m ticks with rock names once applied.
    """
    legend = (
        store if isinstance(store, GeologyLegend) else legend_from_store(store)
    )
    if legend is None or not legend.entries:
        return None
    return tuple(
        (float(e.rho_min), float(e.rho_max), e.color or "#AAAAAA", e.name)
        for e in legend.entries
    )


# ---------------------------------------------------------------------------
# Studio table <-> GeologyLegend bridge
# ---------------------------------------------------------------------------


def table_rows_from_legend(
    legend: GeologyLegend | None,
) -> list[dict[str, Any]]:
    """Flatten *legend* into the row dicts the Studio's editable
    ``dash_table.DataTable`` (:data:`TABLE_COLUMNS`) expects."""
    if legend is None:
        return []
    return [
        {col: getattr(entry, col) for col in TABLE_COLUMNS}
        for entry in legend.entries
    ]


def legend_from_table_rows(
    rows: list[dict[str, Any]] | None,
    *,
    document_id: str = "pcgl:studio",
    created_by: str = "pyCSAMT Interpretation Studio",
    title: str = "",
    previous: GeologyLegend | None = None,
) -> GeologyLegend:
    """Build (and validate) a :class:`GeologyLegend` from Studio table
    rows. Raises :class:`GeologyLegendValidationError` on invalid rows
    (via :meth:`GeologyLegend.validate`) so the caller can surface it."""
    entries: list[RockEntry] = []
    for i, row in enumerate(rows or []):
        name = str((row or {}).get("name") or "").strip()
        if not name:
            continue
        try:
            rho_min = float(row.get("rho_min"))
            rho_max = float(row.get("rho_max"))
        except (TypeError, ValueError):
            continue
        entries.append(
            RockEntry(
                name=name,
                rho_min=rho_min,
                rho_max=rho_max,
                color=str(row.get("color") or "#AAAAAA") or "#AAAAAA",
                description=str(row.get("description") or ""),
                code=i + 1,
                source=str(row.get("source") or ""),
                pattern_id=str(row.get("pattern_id") or ""),
                pattern_source=str(row.get("pattern_source") or ""),
            )
        )
    legend = GeologyLegend(
        document_id=(previous.document_id if previous else document_id),
        created_at=(previous.created_at if previous else _now()),
        created_by=created_by,
        db=RockDatabase(entries),
        title=title or (previous.title if previous else ""),
        description=previous.description if previous else "",
    )
    legend.validate()
    return legend


# ---------------------------------------------------------------------------
# Auto-suggest
# ---------------------------------------------------------------------------


def auto_suggest_legend(
    rho_values: Any,
    *,
    db: RockDatabase | None = None,
    n_bins: int = 8,
    title: str = "Auto-suggested",
) -> GeologyLegend:
    """Propose a legend covering the *actual* resistivity range of the
    currently rendered model, instead of a rock database's full
    (much wider) literature range.

    Builds ``n_bins`` log-spaced bins across the finite, positive values
    in *rho_values*, classifies each bin's geometric-mean resistivity
    against *db* (default :meth:`RockDatabase.default`), and merges
    adjacent bins that land on the same rock entry -- the same merge
    logic :class:`~pycsamt.geology.lithology.StratigraphicLog` uses for a
    depth column, applied here to a resistivity histogram instead.

    Raises :class:`GeologyLegendValidationError` if *rho_values* has
    fewer than two distinct finite, positive samples (nothing to bin).
    """
    source_db = db if db is not None else RockDatabase.default()
    arr = np.asarray(rho_values, dtype=float).ravel()
    arr = arr[np.isfinite(arr) & (arr > 0)]
    if arr.size < 2 or float(arr.min()) == float(arr.max()):
        raise GeologyLegendValidationError(
            "need at least two distinct positive resistivity samples to "
            "auto-suggest a legend"
        )
    lo, hi = float(arr.min()), float(arr.max())
    edges = np.logspace(np.log10(lo), np.log10(hi), int(max(2, n_bins)) + 1)

    entries: list[RockEntry] = []
    i = 0
    n = len(edges) - 1
    while i < n:
        e0 = source_db.classify(float(np.sqrt(edges[i] * edges[i + 1])))
        j = i + 1
        while j < n:
            mid = float(np.sqrt(edges[j] * edges[j + 1]))
            if source_db.classify(mid).name != e0.name:
                break
            j += 1
        entries.append(
            RockEntry(
                name=e0.name,
                rho_min=float(edges[i]),
                rho_max=float(edges[j]),
                color=e0.color,
                description=e0.description,
                code=len(entries) + 1,
                source=e0.source,
            )
        )
        i = j
    legend = GeologyLegend(
        document_id="pcgl:auto-suggest",
        created_at=_now(),
        created_by="auto-suggest",
        db=RockDatabase(entries),
        title=title,
        description=(
            f"Auto-suggested from {arr.size} resistivity sample(s), "
            f"{lo:.3g}-{hi:.3g} Ω·m."
        ),
    )
    legend.validate()
    return legend


# ---------------------------------------------------------------------------
# preview
# ---------------------------------------------------------------------------


def legend_preview_figure(legend: GeologyLegend | None):
    """A horizontal stacked-bar reading of *legend* -- each entry drawn
    as its own resistivity-range segment, coloured and labelled by name.
    Used by both the Interpretation Studio's Legend-tab preview and the
    Geology section's own canvas (``GEO_VIEW_MODE == "legend"``)."""
    import plotly.graph_objects as go

    fig = go.Figure()
    if legend is None or not legend.entries:
        fig.add_annotation(
            text="No legend yet — import, load the default rock "
            "database, or auto-suggest one",
            showarrow=False,
            xref="paper", yref="paper", x=0.5, y=0.5,
        )
        fig.update_layout(
            xaxis=dict(visible=False), yaxis=dict(visible=False),
            margin=dict(l=10, r=10, t=10, b=10),
        )
        return fig
    for entry in legend.entries:
        fig.add_trace(
            go.Bar(
                x=[entry.rho_max - entry.rho_min],
                base=[entry.rho_min],
                y=["legend"],
                orientation="h",
                marker=dict(color=entry.color),
                name=entry.name,
                hovertemplate=(
                    f"{entry.name}<br>{entry.rho_min:g}"
                    f"–{entry.rho_max:g} Ω·m<extra></extra>"
                ),
            )
        )
    fig.update_layout(
        barmode="stack",
        xaxis=dict(type="log", title="Resistivity (Ω·m)"),
        yaxis=dict(visible=False),
        showlegend=True,
        margin=dict(l=10, r=10, t=10, b=40),
    )
    return fig


# ---------------------------------------------------------------------------
# colour consistency with loaded boreholes
# ---------------------------------------------------------------------------


def borehole_lithology_colors(document: Any) -> dict[str, str]:
    """Map each PCBH lithology-vocabulary entry's name (casefolded,
    trimmed) to the exact colour :mod:`pycsamt.format.borehole` actually
    draws it with.

    A borehole is *observed* (drilled), so where the same unit name
    appears in both an interpretation legend and a loaded borehole's own
    log, the borehole's colour is the one worth keeping -- letting the
    two overlays be read as agreeing or disagreeing on that unit at a
    glance, rather than differing only because each was coloured by an
    independent scheme. Matching is exact (case/whitespace-insensitive)
    only, deliberately not fuzzy: a real digitized log can carry two
    genuinely different labels that merely look similar (e.g.
    "Pyrite-mineralized granodiorite porphyry" vs. "Pyrite-ized
    granodiorite porphyry"), and silently merging those would defeat the
    very comparison this is for.
    """
    from pycsamt.format.borehole import deterministic_color

    out: dict[str, str] = {}
    for entry in getattr(document, "lithologies", None) or []:
        name = str(getattr(entry, "name", "") or "").strip().casefold()
        if not name:
            continue
        color = getattr(entry, "color", None) or deterministic_color(
            getattr(entry, "code", name)
        )
        out[name] = color
    return out


def sync_legend_colors_with_boreholes(
    rows: list[dict[str, Any]] | None,
    borehole_colors: dict[str, str],
) -> tuple[list[dict[str, Any]], int]:
    """Overwrite each Studio Legend row's ``color`` with the matching
    borehole lithology's colour (see :func:`borehole_lithology_colors`),
    leaving rows with no matching borehole unit untouched -- a legend
    resistivity-classifies the *whole* model, so it routinely has units
    no drilled hole happened to intersect, and those keep whatever
    colour they already had.

    Returns ``(updated_rows, n_matched)``.
    """
    updated: list[dict[str, Any]] = []
    n_matched = 0
    for row in rows or []:
        row = dict(row)
        key = str(row.get("name") or "").strip().casefold()
        color = borehole_colors.get(key) if key else None
        if color:
            row["color"] = color
            n_matched += 1
        updated.append(row)
    return updated, n_matched


# ---------------------------------------------------------------------------
# pattern-texture fill (fence / depth-slice only, see pycsamt.map.volume)
# ---------------------------------------------------------------------------


def geology_pattern_stencils_from_store(
    store: dict[str, Any] | GeologyLegend | None,
) -> dict[str, Any]:
    """``{legend entry name: ink-density stencil array}`` for every
    entry with an assigned pattern (``pattern_id``/``pattern_source`` --
    see the Interpretation Studio's Patterns tab) -- the
    ``VolumeMapOptions.geology_patterns`` dict Map View's ``geology_fill
    ="pattern"`` texture rendering consumes.

    An entry whose pattern can no longer be loaded (its pack was
    deleted, the tile removed, ...) is skipped rather than raised -- a
    stale pattern reference should degrade that one band back to its
    solid legend colour, not break the whole render.
    """
    legend = (
        store if isinstance(store, GeologyLegend) else legend_from_store(store)
    )
    if legend is None:
        return {}
    from pycsamt.app._patterns import tile_stencil_array

    out: dict[str, Any] = {}
    for entry in legend.entries:
        if not entry.pattern_id or not entry.pattern_source:
            continue
        try:
            out[entry.name] = tile_stencil_array(
                entry.pattern_source, entry.pattern_id
            )
        except Exception:  # noqa: BLE001
            continue
    return out


def _now() -> str:
    from datetime import datetime, timezone

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
