# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCGS (structural geology) upload, preview, and Plotly
scene-trace adapters for pyCSAMT applications.

Mirrors :mod:`pycsamt.app._borehole` / :mod:`pycsamt.app._geology`: this
module bridges the versioned PCGS document
(:mod:`pycsamt.format.structure`) and Map View's Dash layer -- a
JSON-safe application store, three flat-table Studio drafts (planar /
linear / faults, matching :class:`pycsamt.geology.structural.
StructuralModel`'s own three lists), a 2-D preview figure, and the
scene traces that place a fault trace / measurement into the main 3-D
view the same way a borehole is placed
(:func:`pycsamt.app._borehole.scene_borehole_traces`).

Structural ``x`` is already profile-relative -- exactly the along-line
axis a fence panel uses for its own x-axis -- so placing it in the 3-D
scene needs no lat/lon reprojection, only the same per-line
``offset_shift``/``line_spacing`` panel math the borehole/point overlays
already use.
"""

from __future__ import annotations

import base64
import binascii
import json
import math
from typing import Any, Callable

from pycsamt.format.structure import (
    StructModel,
    read_structure,
    structure_from_dict,
    structure_to_dict,
)
from pycsamt.geology.structural import (
    FaultTrace,
    LinearMeasurement,
    StructuralMeasurement,
    StructuralModel,
)

__all__ = [
    "PLANAR_COLUMNS",
    "LINEAR_COLUMNS",
    "FAULT_COLUMNS",
    "decode_structure_json_upload",
    "decode_structure_csv_upload",
    "structure_from_store",
    "store_from_structure",
    "table_rows_from_model",
    "model_from_table_rows",
    "structure_section_figure",
    "structure_scene_traces",
]

_MAX_BYTES = 2 * 1024 * 1024

PLANAR_COLUMNS = (
    "x", "line", "kind", "strike_deg", "dip_deg", "dip_direction_deg",
    "z", "station", "confidence", "notes",
)
LINEAR_COLUMNS = (
    "x", "line", "kind", "trend_deg", "plunge_deg", "z", "station",
    "confidence", "notes",
)
FAULT_COLUMNS = (
    "x", "line", "dip_deg", "downthrown_side", "sense", "throw_m",
    "strike_deg", "z_top", "confidence", "evidence", "notes",
)

_SENSE_COLOR = {
    "normal": "#3b82f6",
    "reverse": "#ef4444",
    "strike_slip": "#a855f7",
    "unknown": "#64748b",
}


# ---------------------------------------------------------------------------
# upload decoding
# ---------------------------------------------------------------------------


def _decode_b64(contents: str, max_bytes: int = _MAX_BYTES) -> bytes:
    if not isinstance(contents, str) or "," not in contents:
        raise ValueError("invalid upload payload")
    header, encoded = contents.split(",", 1)
    if ";base64" not in header.lower():
        raise ValueError("upload must be base64 encoded")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("invalid base64 upload") from error
    if len(raw) > max_bytes:
        raise ValueError(f"upload exceeds {max_bytes} bytes")
    return raw


def decode_structure_json_upload(
    contents: str, filename: str | None = None
) -> dict[str, Any]:
    """Validate a Dash data URL and return a JSON-safe PCGS store."""
    name = filename or "structure.pcgs.json"
    raw = _decode_b64(contents)
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("PCGS upload is not valid UTF-8 JSON") from error
    doc = structure_from_dict(payload)
    return store_from_structure(doc, filename=name)


def decode_structure_csv_upload(
    contents: str, filename: str, *, kind: str
) -> list[dict[str, Any]]:
    """Decode one planar/linear/faults CSV upload into Studio table rows.

    *kind* selects which of the three tables the CSV feeds
    (``"planar"``, ``"linear"``, or ``"faults"``).
    """
    import tempfile
    from pathlib import Path

    raw = _decode_b64(contents)
    folder = tempfile.mkdtemp(prefix="pcgs-")
    path = Path(folder) / (filename or f"{kind}.csv")
    path.write_bytes(raw)
    if kind == "planar":
        model = StructuralModel.from_csv(planar_path=path)
        return _rows(model.planar, PLANAR_COLUMNS)
    if kind == "linear":
        model = StructuralModel.from_csv(linear_path=path)
        return _rows(model.linear, LINEAR_COLUMNS)
    if kind == "faults":
        model = StructuralModel.from_csv(faults_path=path)
        return _rows(model.faults, FAULT_COLUMNS)
    raise ValueError(f"unknown kind {kind!r}")


def structure_from_store(store: dict[str, Any] | None) -> StructModel | None:
    """Restore a validated :class:`StructModel` from an application
    store (``{"document": <pcgs dict>}`` or the raw pcgs dict itself)."""
    if not store:
        return None
    payload = store.get("document", store)
    if not isinstance(payload, dict):
        raise ValueError("invalid PCGS application store")
    return structure_from_dict(payload)


def store_from_structure(
    doc: StructModel, *, filename: str = "structure.pcgs.json"
) -> dict[str, Any]:
    """Build the JSON-safe application store for *doc*."""
    return {
        "filename": filename,
        "document": structure_to_dict(doc),
        "n_planar": len(doc.planar),
        "n_linear": len(doc.linear),
        "n_faults": len(doc.faults),
    }


# ---------------------------------------------------------------------------
# Studio table <-> StructuralModel bridge
# ---------------------------------------------------------------------------


def _rows(items, columns: tuple[str, ...]) -> list[dict[str, Any]]:
    return [{col: getattr(item, col, None) for col in columns} for item in items]


def table_rows_from_model(
    model: StructuralModel | None,
) -> dict[str, list[dict[str, Any]]]:
    """Flatten *model* into the three Studio table row-lists."""
    if model is None:
        return {"planar": [], "linear": [], "faults": []}
    return {
        "planar": _rows(model.planar, PLANAR_COLUMNS),
        "linear": _rows(model.linear, LINEAR_COLUMNS),
        "faults": _rows(model.faults, FAULT_COLUMNS),
    }


def _clean_str(value: Any) -> str | None:
    text = str(value).strip() if value is not None else ""
    return text or None


def _clean_float(value: Any) -> float | None:
    if value is None or (isinstance(value, str) and not value.strip()):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def model_from_table_rows(
    planar_rows: list[dict[str, Any]] | None,
    linear_rows: list[dict[str, Any]] | None,
    fault_rows: list[dict[str, Any]] | None,
) -> StructuralModel:
    """Build a :class:`StructuralModel` from the three Studio tables.

    Rows missing a required numeric field are silently skipped (an
    in-progress "Add row" blank); rows with a required field that fails
    validation (e.g. an inconsistent strike/dip-direction pair) raise
    ``ValueError``, surfaced by the caller.
    """
    planar: list[StructuralMeasurement] = []
    for row in planar_rows or []:
        x, strike, dip, dipdir = (
            _clean_float(row.get(k))
            for k in ("x", "strike_deg", "dip_deg", "dip_direction_deg")
        )
        kind = _clean_str(row.get("kind"))
        if x is None or strike is None or dip is None or dipdir is None or not kind:
            continue
        planar.append(
            StructuralMeasurement(
                x=x, kind=kind, strike_deg=strike, dip_deg=dip,
                dip_direction_deg=dipdir,
                z=_clean_float(row.get("z")),
                station=_clean_str(row.get("station")),
                confidence=_clean_float(row.get("confidence")) or 1.0,
                notes=str(row.get("notes") or ""),
                line=_clean_str(row.get("line")),
            )
        )
    linear: list[LinearMeasurement] = []
    for row in linear_rows or []:
        x, trend, plunge = (
            _clean_float(row.get(k)) for k in ("x", "trend_deg", "plunge_deg")
        )
        kind = _clean_str(row.get("kind"))
        if x is None or trend is None or plunge is None or not kind:
            continue
        linear.append(
            LinearMeasurement(
                x=x, kind=kind, trend_deg=trend, plunge_deg=plunge,
                z=_clean_float(row.get("z")),
                station=_clean_str(row.get("station")),
                confidence=_clean_float(row.get("confidence")) or 1.0,
                notes=str(row.get("notes") or ""),
                line=_clean_str(row.get("line")),
            )
        )
    faults: list[FaultTrace] = []
    for row in fault_rows or []:
        x, dip = (_clean_float(row.get(k)) for k in ("x", "dip_deg"))
        side = _clean_str(row.get("downthrown_side"))
        if x is None or dip is None or not side:
            continue
        faults.append(
            FaultTrace(
                x=x, dip_deg=dip, downthrown_side=side,
                sense=_clean_str(row.get("sense")) or "unknown",
                throw_m=_clean_float(row.get("throw_m")),
                strike_deg=_clean_float(row.get("strike_deg")),
                z_top=_clean_float(row.get("z_top")),
                confidence=_clean_float(row.get("confidence")) or 1.0,
                evidence=str(row.get("evidence") or ""),
                notes=str(row.get("notes") or ""),
                line=_clean_str(row.get("line")),
            )
        )
    return StructuralModel(planar=planar, linear=linear, faults=faults)


# ---------------------------------------------------------------------------
# 2-D preview (Geology section canvas + Studio Structure-tab preview)
# ---------------------------------------------------------------------------


def structure_section_figure(
    model: StructuralModel | None,
    *,
    theme: str = "light",
    depth_extent: float = 300.0,
):
    """A profile-position / depth preview of *model* -- the same view
    (before it goes into the 3-D scene) whether or not the survey has
    more than one line.

    A fault trace is drawn as a straight line through ``(x, z_top)``,
    tilted by its apparent dip (steeper = closer to vertical); a planar
    measurement is a short tilted tick through ``(x, z)`` labelled
    ``strike/dip``; a linear measurement is a diamond marker labelled
    ``trend/plunge``. This is a deliberately simplified 2-D reading of
    a 3-D attitude -- exactly the same simplification
    :class:`~pycsamt.geology.structural.FaultTrace` itself documents for
    ``dip_deg`` (the section's *apparent* dip, not necessarily the true
    3-D dip).
    """
    import plotly.graph_objects as go

    fig = go.Figure()
    dark = theme == "dark"
    grid = "#313244" if dark else "#e5e7eb"
    text_color = "#cdd6f4" if dark else "#374151"

    model = model or StructuralModel()
    half = max(float(depth_extent) * 0.06, 5.0)

    for fault in model.faults:
        z_top = float(fault.z_top or 0.0)
        z_bot = z_top + float(depth_extent)
        dip = max(1.0, min(89.0, float(fault.dip_deg)))
        direction = 1.0 if fault.downthrown_side == "right" else -1.0
        dx = direction * float(depth_extent) / math.tan(math.radians(dip))
        color = _SENSE_COLOR.get(fault.sense, _SENSE_COLOR["unknown"])
        label = (
            f"fault {fault.x:.0f} m — {fault.sense}, dip {fault.dip_deg:.0f}°"
            + (f", throw {fault.throw_m:.0f} m" if fault.throw_m else "")
        )
        fig.add_trace(
            go.Scatter(
                x=[fault.x, fault.x + dx],
                y=[z_top, z_bot],
                mode="lines",
                line=dict(color=color, width=3),
                name=label,
                hovertemplate=label + "<extra></extra>",
                showlegend=False,
            )
        )

    for m in model.planar:
        z = float(m.z or 0.0)
        dip = max(1.0, min(89.0, float(m.dip_deg)))
        dx = half / math.tan(math.radians(dip))
        label = (
            f"{m.kind} @ {m.x:.0f} m — {m.strike_deg:.0f}/"
            f"{m.dip_deg:.0f}→{m.dip_direction_deg:.0f}"
        )
        fig.add_trace(
            go.Scatter(
                x=[m.x - dx, m.x + dx],
                y=[z - half, z + half],
                mode="lines+markers",
                line=dict(color="#22c55e", width=2),
                marker=dict(size=4, color="#22c55e"),
                name=label,
                hovertemplate=label + "<extra></extra>",
                showlegend=False,
            )
        )

    for m in model.linear:
        z = float(m.z or 0.0)
        label = f"{m.kind} @ {m.x:.0f} m — {m.trend_deg:.0f}/{m.plunge_deg:.0f}"
        fig.add_trace(
            go.Scatter(
                x=[m.x],
                y=[z],
                mode="markers",
                marker=dict(size=10, color="#f59e0b", symbol="diamond"),
                name=label,
                hovertemplate=label + "<extra></extra>",
                showlegend=False,
            )
        )

    fig.update_layout(
        template="plotly_dark" if dark else "plotly_white",
        xaxis=dict(title="Profile position (m)", gridcolor=grid),
        yaxis=dict(
            title="Depth (m)", autorange="reversed", gridcolor=grid,
        ),
        font=dict(color=text_color),
        margin=dict(l=50, r=10, t=10, b=40),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    if not (model.faults or model.planar or model.linear):
        fig.add_annotation(
            text="No structural data yet — import or add a row",
            showarrow=False,
            xref="paper", yref="paper", x=0.5, y=0.5,
            font=dict(color=text_color, size=13),
        )
    return fig


# ---------------------------------------------------------------------------
# 3-D scene traces
# ---------------------------------------------------------------------------


def structure_scene_traces(
    model: StructuralModel | None,
    *,
    line_offsets: dict[str, float],
    default_line: str | None = None,
    azimuth_deg: float = 0.0,
    surface: Callable[[float], float] | None = None,
    depth_extent: float = 300.0,
    opacity: float = 0.85,
):
    """Return 3-D scene traces placing *model* into the fence/block/iso
    scene.

    ``line_offsets`` maps a line name to its already-normalised
    cross-strike scene offset (the same ``v`` a fence panel places
    itself at -- see :func:`pycsamt.map.volume._line_offset` /
    :func:`pycsamt.app._borehole.scene_borehole_traces`'s caller for how
    that offset is computed). Items whose ``line`` is unset use
    *default_line*'s offset (typically the first/only line). Items on an
    unknown line are skipped rather than guessed at.

    ``surface`` is a ``u -> elevation`` interpolator (see
    :func:`pycsamt.map.borehole_align.surface_from_sections`); ``None``
    uses a flat ``elevation = 0`` datum, matching the ``"zero"`` datum
    boreholes fall back to without topography.
    """
    import numpy as np
    import plotly.graph_objects as go

    az = math.radians(float(azimuth_deg))
    sin_az, cos_az = math.sin(az), math.cos(az)

    def _scene_xy(u: float, line: str | None) -> tuple[float, float] | None:
        key = line if line in line_offsets else default_line
        if key is None or key not in line_offsets:
            return None
        v = float(line_offsets[key])
        return u + v * sin_az, v * cos_az

    def _elev(u: float) -> float:
        return float(surface(u)) if surface is not None else 0.0

    model = model or StructuralModel()
    traces: list[Any] = []

    for fault in model.faults:
        xy = _scene_xy(float(fault.x), fault.line)
        if xy is None:
            continue
        x0, y0 = xy
        z_top = _elev(float(fault.x)) - float(fault.z_top or 0.0)
        z_bot = z_top - float(depth_extent)
        dip = max(1.0, min(89.0, float(fault.dip_deg)))
        direction = 1.0 if fault.downthrown_side == "right" else -1.0
        dx = direction * float(depth_extent) / math.tan(math.radians(dip))
        x1 = x0 + dx * cos_az
        y1 = y0 - dx * sin_az
        color = _SENSE_COLOR.get(fault.sense, _SENSE_COLOR["unknown"])
        label = (
            f"Fault @ {fault.x:.0f} m ({fault.sense}, "
            f"dip {fault.dip_deg:.0f}°)"
        )
        traces.append(
            go.Mesh3d(
                x=[x0, x1, x1, x0],
                y=[y0, y1, y1, y0],
                z=[z_top, z_top, z_bot, z_bot],
                i=[0, 0], j=[1, 2], k=[2, 3],
                color=color,
                opacity=float(opacity),
                flatshading=True,
                name=label,
                hovertemplate=label + "<extra></extra>",
                showscale=False,
            )
        )

    tick_x, tick_y, tick_z, tick_text = [], [], [], []
    for m in model.planar:
        xy = _scene_xy(float(m.x), m.line)
        if xy is None:
            continue
        x0, y0 = xy
        z0 = _elev(float(m.x)) - float(m.z or 0.0)
        tick_x.append(x0)
        tick_y.append(y0)
        tick_z.append(z0)
        tick_text.append(
            f"{m.kind}: {m.strike_deg:.0f}/{m.dip_deg:.0f}"
            f"→{m.dip_direction_deg:.0f}"
        )
    if tick_x:
        traces.append(
            go.Scatter3d(
                x=tick_x, y=tick_y, z=tick_z,
                mode="markers",
                marker=dict(size=5, color="#22c55e", symbol="square"),
                text=tick_text,
                hovertemplate="%{text}<extra></extra>",
                name="Planar measurements",
                showlegend=False,
            )
        )

    lin_x, lin_y, lin_z, lin_text = [], [], [], []
    for m in model.linear:
        xy = _scene_xy(float(m.x), m.line)
        if xy is None:
            continue
        x0, y0 = xy
        z0 = _elev(float(m.x)) - float(m.z or 0.0)
        lin_x.append(x0)
        lin_y.append(y0)
        lin_z.append(z0)
        lin_text.append(f"{m.kind}: {m.trend_deg:.0f}/{m.plunge_deg:.0f}")
    if lin_x:
        traces.append(
            go.Scatter3d(
                x=lin_x, y=lin_y, z=lin_z,
                mode="markers",
                marker=dict(size=5, color="#f59e0b", symbol="diamond"),
                text=lin_text,
                hovertemplate="%{text}<extra></extra>",
                name="Linear measurements",
                showlegend=False,
            )
        )

    return traces
