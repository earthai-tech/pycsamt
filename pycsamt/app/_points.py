# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCPT (points / targets) upload and Plotly adapters."""

from __future__ import annotations

import base64
import binascii
import json
from typing import Any

from pycsamt.format import (
    PointSet,
    point_set_from_dict,
    point_set_to_dict,
)

__all__ = [
    "decode_points_upload",
    "pointset_from_store",
    "point_map_markers",
    "point_scene_traces",
]

_MAX_BYTES = 4 * 1024 * 1024
_KIND_COLOR = {
    "target": "#ef4444",
    "planned_borehole": "#f59e0b",
    "sample": "#22c55e",
    "anomaly": "#a855f7",
    "poi": "#3b82f6",
    "other": "#64748b",
}


def decode_points_upload(
    contents: str, filename: str | None = None
) -> dict[str, Any]:
    """Validate a Dash data URL and return a JSON-safe PCPT store."""
    name = filename or "points.pcpt.json"
    lower = name.lower()
    if not isinstance(contents, str) or "," not in contents:
        raise ValueError("invalid PCPT upload payload")
    header, encoded = contents.split(",", 1)
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("invalid base64 PCPT upload") from error
    if len(raw) > _MAX_BYTES:
        raise ValueError(f"PCPT upload exceeds {_MAX_BYTES} bytes")

    if lower.endswith((".csv", ".tsv", ".txt")):
        import tempfile
        from pathlib import Path

        from pycsamt.format import points_from_csv

        folder = tempfile.mkdtemp(prefix="pcpt-")
        path = Path(folder) / name
        path.write_bytes(raw)
        point_set = points_from_csv(path)
    elif lower.endswith(".xlsx"):
        from pycsamt.format import points_from_xlsx

        point_set = points_from_xlsx(raw)
    else:
        point_set = point_set_from_dict(json.loads(raw.decode("utf-8")))
    return {
        "filename": name,
        "document": point_set_to_dict(point_set),
        "n_points": len(point_set.points),
    }


def pointset_from_store(store: dict[str, Any] | None) -> PointSet | None:
    """Restore a validated :class:`PointSet` from an application store."""
    if not store:
        return None
    payload = store.get("document", store)
    if not isinstance(payload, dict):
        raise ValueError("invalid PCPT application store")
    return point_set_from_dict(payload)


def point_map_markers(
    store_or_set: dict[str, Any] | PointSet | None,
    *,
    show_labels: bool = True,
):
    """Return a basemap scatter trace for PCPT points, coloured by kind."""
    import plotly.graph_objects as go

    point_set = _resolve(store_or_set)
    if point_set is None:
        return []
    lonlat = point_set.lonlat()
    if not lonlat:
        return []
    by_id = {p.id: p for p in point_set.points}
    lats, lons, texts, hovers, colors = [], [], [], [], []
    for pid, (lat, lon) in lonlat.items():
        point = by_id.get(pid)
        lats.append(lat)
        lons.append(lon)
        texts.append(point.name or pid if point else pid)
        colors.append(
            (point.color if point and point.color else None)
            or _KIND_COLOR.get(point.kind if point else "poi", "#3b82f6")
        )
        hovers.append(_hover(point, pid))
    scatter = getattr(go, "Scattermap", None) or go.Scattermapbox
    return [
        scatter(
            lat=lats,
            lon=lons,
            mode="markers+text" if show_labels else "markers",
            text=texts if show_labels else None,
            textposition="top center",
            marker={"size": 10, "color": colors, "symbol": "circle"},
            name="Targets",
            hovertext=hovers,
            hoverinfo="text",
            showlegend=False,
        )
    ]


def point_scene_traces(
    store_or_set: dict[str, Any] | PointSet | None,
    frame,
    *,
    surface=None,
    datum: str = "surface",
    azimuth_deg: float = 0.0,
    offset_shift: float = 0.0,
    offset_scale: float = 1.0,
    show_labels: bool = True,
):
    """Return 3-D scene traces (marker + depth stem) for PCPT points.

    ``frame`` is a :class:`pycsamt.map.geometry.SurveyFrame`. Points with
    no locatable position are skipped. ``offset_shift`` / ``offset_scale``
    match the cross-strike normalisation the 3-D volume applies to its
    line panels (see
    :func:`pycsamt.map.borehole_align.align_boreholes_to_scene`).
    """
    import math

    import plotly.graph_objects as go

    point_set = _resolve(store_or_set)
    if point_set is None or frame is None:
        return []
    lonlat = point_set.lonlat()
    by_id = {p.id: p for p in point_set.points}
    az = math.radians(float(azimuth_deg))
    sin_az, cos_az = math.sin(az), math.cos(az)

    xs, ys, zs, texts, hovers, colors = [], [], [], [], [], []
    stems: list[Any] = []
    for pid, (lat, lon) in lonlat.items():
        point = by_id.get(pid)
        u, v = frame.project(lat, lon)
        v = (v - float(offset_shift)) * float(offset_scale)
        x = u + v * sin_az
        y = v * cos_az
        if datum == "surface" and surface is not None:
            z_top = float(surface(x))
        elif datum == "collar_z" and point and point.z is not None:
            z_top = float(point.z)
        else:
            z_top = 0.0
        top = z_top - float(point.depth_top or 0.0) if point else z_top
        colour = (
            (point.color if point and point.color else None)
            or _KIND_COLOR.get(point.kind if point else "poi", "#3b82f6")
        )
        xs.append(x)
        ys.append(y)
        zs.append(top)
        texts.append((point.name or pid) if point else pid)
        hovers.append(_hover(point, pid))
        colors.append(colour)
        if point and point.depth_bottom is not None:
            bottom = z_top - float(point.depth_bottom)
            stems.append(
                go.Scatter3d(
                    x=[x, x],
                    y=[y, y],
                    z=[top, bottom],
                    mode="lines",
                    line={"color": colour, "width": 5},
                    hoverinfo="skip",
                    showlegend=False,
                )
            )
    if not xs:
        return []
    marker_trace = go.Scatter3d(
        x=xs,
        y=ys,
        z=zs,
        mode="markers+text" if show_labels else "markers",
        text=texts if show_labels else None,
        textposition="top center",
        marker={"size": 5, "color": colors, "symbol": "diamond"},
        name="Targets",
        hovertext=hovers,
        hoverinfo="text",
        showlegend=False,
    )
    return [marker_trace, *stems]


def _resolve(value: dict[str, Any] | PointSet | None) -> PointSet | None:
    if isinstance(value, PointSet):
        return value
    return pointset_from_store(value)


def _hover(point, pid: str) -> str:
    if point is None:
        return pid
    bits = [point.name or pid, f"kind: {point.kind}"]
    if point.depth_top is not None:
        window = f"{point.depth_top:g}"
        if point.depth_bottom is not None:
            window += f"–{point.depth_bottom:g}"
        bits.append(f"depth {window} m")
    if point.note:
        bits.append(point.note)
    return "<br>".join(bits)
