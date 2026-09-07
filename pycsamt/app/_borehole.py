# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Shared PCBH upload and Plotly adapters for pyCSAMT applications."""

from __future__ import annotations

import base64
import binascii
import json
from typing import Any

from pycsamt.format.borehole import (
    DEFAULT_MAX_BYTES,
    DisplayRadiusPolicy,
    PCBHDocument,
    build_render_model,
    document_from_builder,
    document_to_builder,
    pcbh_from_dict,
    pcbh_to_dict,
)

__all__ = [
    "decode_pcbh_upload",
    "document_from_store",
    "pcbh_plotly_traces",
    "add_pcbh_to_figure",
    "add_embedded_pcbh_to_pcsf_figure",
    "builder_draft_from_document",
    "document_from_builder_draft",
    "strip_log_figure",
    "scene_borehole_traces",
    "collar_map_markers",
    "document_collar_lonlat",
]


def decode_pcbh_upload(
    contents: str,
    filename: str | None = None,
    *,
    max_bytes: int = DEFAULT_MAX_BYTES,
) -> dict[str, Any]:
    """Validate a Dash data URL and return a JSON-safe PCBH store."""
    name = filename or "boreholes.pcbh.json"
    if not name.lower().endswith((".pcbh.json", ".json")):
        raise ValueError("PCBH uploads must use .pcbh.json or .json")
    if not isinstance(contents, str) or "," not in contents:
        raise ValueError("invalid PCBH upload payload")
    header, encoded = contents.split(",", 1)
    if ";base64" not in header.lower():
        raise ValueError("PCBH upload must be base64 encoded")
    if len(encoded) > ((max_bytes + 2) // 3) * 4 + 4:
        raise ValueError(f"PCBH upload exceeds {max_bytes} bytes")
    try:
        raw = base64.b64decode(encoded, validate=True)
    except (binascii.Error, ValueError) as error:
        raise ValueError("invalid base64 PCBH upload") from error
    if len(raw) > max_bytes:
        raise ValueError(f"PCBH upload exceeds {max_bytes} bytes")
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ValueError("PCBH upload is not valid UTF-8 JSON") from error
    document = pcbh_from_dict(payload)
    return {
        "filename": name,
        "document": pcbh_to_dict(document),
        "n_boreholes": len(document.boreholes),
    }


def document_from_store(store: dict[str, Any] | None) -> PCBHDocument | None:
    """Restore a validated document from an application store."""
    if not store:
        return None
    payload = store.get("document", store)
    if not isinstance(payload, dict):
        raise ValueError("invalid PCBH application store")
    return pcbh_from_dict(payload)


def pcbh_plotly_traces(
    document: PCBHDocument,
    *,
    family: str = "lithology",
    visible: bool = True,
    show_labels: bool = True,
    opacity: float = 0.9,
    radius_mode: str = "auto",
    radius: float | None = None,
    selected_ids: set[str] | None = None,
) -> list[Any]:
    """Translate the viewer-neutral PCBH model to Plotly 3-D traces."""
    if not visible:
        return []
    import plotly.graph_objects as go

    policy = DisplayRadiusPolicy(
        mode=radius_mode,
        fixed_radius=radius if radius_mode == "fixed" else None,
        exaggeration=(radius or 1.0)
        if radius_mode == "exaggeration"
        else 1.0,
    )
    model = build_render_model(
        document,
        family=family,
        selected_ids=selected_ids,
        radius_policy=policy,
    )
    alpha = min(1.0, max(0.0, float(opacity)))
    traces: list[Any] = []
    for hole in model.boreholes:
        for segment in hole.interval_segments:
            points = segment.points
            traces.append(
                go.Scatter3d(
                    x=[point.x for point in points],
                    y=[point.y for point in points],
                    z=[point.z for point in points],
                    mode="lines",
                    name=segment.metadata.get("label")
                    or segment.metadata.get("code")
                    or family,
                    legendgroup=f"pcbh-{family}",
                    line={
                        "color": segment.color,
                        "width": max(
                            3.0,
                            min(18.0, segment.display_radius * 4),
                        ),
                    },
                    opacity=alpha,
                    customdata=[
                        [segment.borehole_id, point.md] for point in points
                    ],
                    hovertemplate=(
                        "Borehole %{customdata[0]}<br>"
                        "MD %{customdata[1]:.2f}<br>"
                        + f"{family}: {_segment_label(segment.metadata)}"
                        + "<extra></extra>"
                    ),
                    showlegend=False,
                )
            )
        if not hole.interval_segments:
            points = hole.centerline.points
            traces.append(
                go.Scatter3d(
                    x=[point.x for point in points],
                    y=[point.y for point in points],
                    z=[point.z for point in points],
                    mode="lines",
                    name=hole.borehole_id,
                    line={"color": hole.centerline.color, "width": 4},
                    opacity=alpha,
                    hovertemplate=f"{hole.borehole_id}<extra></extra>",
                    showlegend=False,
                )
            )
        for contact in hole.contacts:
            traces.append(
                go.Scatter3d(
                    x=[contact.position[0]],
                    y=[contact.position[1]],
                    z=[contact.position[2]],
                    mode="markers",
                    marker={
                        "size": 5,
                        "color": "#ffffff",
                        "line": {"width": 2, "color": "#111827"},
                    },
                    name="Contact",
                    customdata=[[contact.borehole_id, contact.md]],
                    hovertemplate=(
                        "Borehole %{customdata[0]}<br>"
                        "Contact MD %{customdata[1]:.2f}<extra></extra>"
                    ),
                    showlegend=False,
                )
            )
        for glyph in hole.structure_glyphs:
            traces.append(
                go.Scatter3d(
                    x=[glyph.position[0]],
                    y=[glyph.position[1]],
                    z=[glyph.position[2]],
                    mode="markers",
                    marker={"size": 6, "color": "#ec4899", "symbol": "x"},
                    name=glyph.kind,
                    customdata=[[glyph.borehole_id, glyph.md]],
                    hovertemplate=(
                        "%{customdata[0]}<br>Structure "
                        "%{customdata[1]:.2f} MD<extra></extra>"
                    ),
                    showlegend=False,
                )
            )
        collar = hole.collar
        traces.append(
            go.Scatter3d(
                x=[collar.position[0]],
                y=[collar.position[1]],
                z=[collar.position[2]],
                mode="markers+text" if show_labels else "markers",
                text=[collar.label] if show_labels else None,
                textposition="top center",
                marker={
                    "size": 7 if collar.selected else 5,
                    "color": "#f59e0b",
                    "symbol": "diamond",
                },
                name=collar.label,
                customdata=[[collar.borehole_id]],
                hovertemplate="%{text}<br>%{customdata[0]}<extra></extra>",
                showlegend=False,
            )
        )
    return traces


def add_pcbh_to_figure(figure: Any, store: dict[str, Any] | None, **options):
    """Add PCBH traces to a Plotly figure and preserve its existing model."""
    document = document_from_store(store)
    if document is None:
        return figure
    for trace in pcbh_plotly_traces(document, **options):
        figure.add_trace(trace)
    return figure


def add_embedded_pcbh_to_pcsf_figure(
    figure: Any,
    model: Any,
    *,
    vertical_offset: float | None = None,
    opacity: float = 0.9,
):
    """Render an embedded PCBH in native/synthesized PCSF block space."""
    import plotly.graph_objects as go

    from pycsamt.format.borehole import align_pcbh_to_pcsf, extract_pcbh

    document = extract_pcbh(model)
    if document is None:
        return figure
    report = align_pcbh_to_pcsf(
        document,
        model,
        vertical_offset=vertical_offset,
    )
    alpha = min(1.0, max(0.0, float(opacity)))
    for result in report.bounds_results:
        points = report.trajectories[result.borehole_id]
        figure.add_trace(
            go.Scatter3d(
                x=[point.x for point in points],
                y=[point.y for point in points],
                z=[-point.z for point in points],
                mode="lines",
                line={"color": "#f59e0b", "width": 7},
                name=result.borehole_id,
                customdata=[
                    [result.borehole_id, point.md, result.relation]
                    for point in points
                ],
                hovertemplate=(
                    "Borehole %{customdata[0]}<br>"
                    "MD %{customdata[1]:.2f}<br>"
                    "Bounds: %{customdata[2]}<extra></extra>"
                ),
                opacity=alpha,
                showlegend=False,
            )
        )
    return figure


def _segment_label(metadata: dict[str, Any]) -> str:
    return metadata.get("label") or metadata.get("code") or "unknown"


# ---------------------------------------------------------------------------
# builder-draft bridge (shared by the web page and the Map View studio)
# ---------------------------------------------------------------------------


def builder_draft_from_document(
    store_or_document: dict[str, Any] | PCBHDocument | None,
) -> dict[str, Any] | None:
    """Return a flat editable builder draft for a document or store."""
    if store_or_document is None:
        return None
    document = (
        store_or_document
        if isinstance(store_or_document, PCBHDocument)
        else document_from_store(store_or_document)
    )
    if document is None:
        return None
    return document_to_builder(document)


def document_from_builder_draft(
    draft: dict[str, Any] | None,
) -> PCBHDocument | None:
    """Construct and validate a document from a builder draft."""
    if not draft:
        return None
    return document_from_builder(draft)


# ---------------------------------------------------------------------------
# strip-log ("linear") figure
# ---------------------------------------------------------------------------


def strip_log_figure(
    document: PCBHDocument | dict[str, Any] | None,
    *,
    hole_ids: list[str] | None = None,
    family: str = "lithology",
    theme: str = "light",
    model_tracks: dict[str, list[tuple[float, float]]] | None = None,
    show_legend: bool = True,
):
    """Return a side-by-side strip-log (depth column) figure.

    One column per borehole: stacked, vocabulary-coloured lithology
    intervals with measured depth increasing downward. ``model_tracks``
    optionally maps a borehole id to ``[(depth, resistivity), ...]`` from
    an inversion model sampled along the hole, drawn as a thin overlay
    line so the log can be compared with the imaged resistivity.
    """
    import plotly.graph_objects as go

    if not isinstance(document, PCBHDocument):
        document = document_from_store(document)
    figure = go.Figure()
    if document is None:
        figure.update_layout(
            template=_template(theme),
            annotations=[
                {
                    "text": "No boreholes — import or build one",
                    "showarrow": False,
                    "font": {"size": 13},
                }
            ],
        )
        return figure

    model = build_render_model(document, family=family)
    colours = {
        entry.code: entry.color
        for entry in document.lithologies
        if entry.color
    }
    holes = [
        hole
        for hole in model.boreholes
        if hole_ids is None or hole.borehole_id in hole_ids
    ]
    order = [hole.borehole_id for hole in holes]
    seen_labels: set[str] = set()
    max_depth = 0.0
    for hole in holes:
        for segment in hole.interval_segments:
            raw_label = _segment_label(segment.metadata)
            label = _short(raw_label, 34)
            colour = (
                colours.get(segment.metadata.get("code"))
                or segment.color
            )
            thickness = float(segment.to_md) - float(segment.from_md)
            max_depth = max(max_depth, float(segment.to_md))
            figure.add_trace(
                go.Bar(
                    x=[hole.borehole_id],
                    y=[thickness],
                    base=[float(segment.from_md)],
                    width=0.55,
                    marker={"color": colour, "line": {"width": 0}},
                    name=label,
                    legendgroup=label,
                    showlegend=show_legend and label not in seen_labels,
                    hovertemplate=(
                        f"{hole.borehole_id}<br>{raw_label}<br>"
                        f"%{{base:.1f}}–{segment.to_md:.1f} m<extra></extra>"
                    ),
                )
            )
            seen_labels.add(label)
        track = (model_tracks or {}).get(hole.borehole_id)
        if track:
            depths = [row[0] for row in track]
            rho = [row[1] for row in track]
            lo, span = min(rho), (max(rho) - min(rho)) or 1.0
            base = order.index(hole.borehole_id)
            scaled = [base - 0.275 + 0.55 * (v - lo) / span for v in rho]
            figure.add_trace(
                go.Scatter(
                    x=scaled,
                    y=depths,
                    mode="lines",
                    line={"color": "#111827", "width": 1.4},
                    name=f"{hole.borehole_id} model ρ",
                    customdata=rho,
                    hovertemplate=(
                        "model ρ %{customdata:.1f} Ω·m @ "
                        "%{y:.0f} m<extra></extra>"
                    ),
                    showlegend=False,
                )
            )

    figure.update_layout(
        template=_template(theme),
        barmode="overlay",
        bargap=0.4,
        showlegend=show_legend,
        legend={
            "title": {"text": family.title()},
            "font": {"size": 9},
            "y": 1,
            "yanchor": "top",
        },
        margin={"l": 56, "r": 8, "t": 30, "b": 30},
        xaxis={
            "type": "category",
            "categoryorder": "array",
            "categoryarray": order,
        },
        yaxis={
            "title": "Measured depth (m)",
            "autorange": "reversed",
            "range": [max_depth * 1.02, 0.0],
        },
    )
    return figure


# ---------------------------------------------------------------------------
# 3-D scene traces (from pycsamt.map.borehole_align.SceneAlignment)
# ---------------------------------------------------------------------------


def scene_borehole_traces(
    alignment,
    *,
    as_tubes: bool = True,
    opacity: float = 0.9,
    show_labels: bool = True,
    collar_size: float = 5.0,
    tube_sides: int = 8,
):
    """Translate an aligned scene into Plotly 3-D traces.

    ``alignment`` is a
    :class:`pycsamt.map.borehole_align.SceneAlignment`. Interval segments
    become coloured tubes (``Mesh3d``) or fat poly-lines; each hole also
    gets a collar diamond and, optionally, a label.
    """
    import numpy as np
    import plotly.graph_objects as go

    alpha = min(1.0, max(0.0, float(opacity)))
    traces: list[Any] = []
    for hole in alignment.placed:
        for segment in hole.segments:
            points = np.asarray(segment.points, dtype=float)
            if points.shape[0] < 2:
                continue
            label = segment.metadata.get("label") or segment.metadata.get(
                "code"
            )
            hover = (
                f"{hole.borehole_id}<br>{label}<br>"
                f"{segment.from_md:.1f}–{segment.to_md:.1f} m MD"
                "<extra></extra>"
            )
            if as_tubes:
                mesh = _tube_mesh(
                    points, hole.display_radius, tube_sides
                )
                if mesh is not None:
                    vertices, faces = mesh
                    traces.append(
                        go.Mesh3d(
                            x=vertices[:, 0],
                            y=vertices[:, 1],
                            z=vertices[:, 2],
                            i=faces[:, 0],
                            j=faces[:, 1],
                            k=faces[:, 2],
                            color=segment.color,
                            opacity=alpha,
                            name=label or hole.borehole_id,
                            hovertext=hover.replace("<extra></extra>", ""),
                            hoverinfo="text",
                            showlegend=False,
                        )
                    )
                    continue
            traces.append(
                go.Scatter3d(
                    x=points[:, 0],
                    y=points[:, 1],
                    z=points[:, 2],
                    mode="lines",
                    line={"color": segment.color, "width": 10},
                    opacity=alpha,
                    name=label or hole.borehole_id,
                    hovertemplate=hover,
                    showlegend=False,
                )
            )
        if not hole.segments and hole.centerline:
            line = np.asarray(hole.centerline, dtype=float)
            traces.append(
                go.Scatter3d(
                    x=line[:, 0],
                    y=line[:, 1],
                    z=line[:, 2],
                    mode="lines",
                    line={"color": "#111827", "width": 5},
                    name=hole.borehole_id,
                    hovertemplate=f"{hole.borehole_id}<extra></extra>",
                    showlegend=False,
                )
            )
        cx, cy, cz = hole.collar_scene
        traces.append(
            go.Scatter3d(
                x=[cx],
                y=[cy],
                z=[cz],
                mode="markers+text" if show_labels else "markers",
                text=[hole.name] if show_labels else None,
                textposition="top center",
                marker={
                    "size": float(collar_size),
                    "color": "#f59e0b",
                    "symbol": "diamond",
                },
                name=hole.name,
                hovertemplate=(
                    f"{hole.name}<br>{hole.borehole_id}"
                    f"<br>relation: {hole.relation}<extra></extra>"
                ),
                showlegend=False,
            )
        )
    return traces


def borehole_patch_traces(
    alignment,
    *,
    half_width: float = 10.0,
    opacity: float = 0.98,
):
    """Patch each hole's own logged interval onto the interpreted panel.

    ``alignment`` is a :class:`pycsamt.map.borehole_align.SceneAlignment`
    (the same one :func:`scene_borehole_traces` renders as tubes). Where
    a resistivity-band "Interpretation legend" cannot separate two
    facies (see :mod:`pycsamt.format.geology`), the borehole's actual
    logged lithology is ground truth at that exact position — this
    stamps it as an opaque coloured decal directly onto the fence/
    depth-slice panel, at the hole's true along-strike position and
    depth, rather than relying on the resistivity colour alone.

    Each decal is a flat rectangle (``go.Mesh3d``, two triangles)
    spanning ``2 * half_width`` metres along-strike, at the hole's own
    cross-strike offset, over the segment's ``[to_md, from_md]`` depth
    range — i.e. exactly the corridor :func:`pycsamt.map.borehole_align.
    align_boreholes_to_scene` already placed the hole's tube in, just
    wider and opaque so it reads as a patch rather than a thin line.
    """
    import numpy as np
    import plotly.graph_objects as go

    traces: list[Any] = []
    for hole in alignment.placed:
        if not hole.segments:
            continue
        cx, cy, _ = hole.collar_scene
        x0, x1 = cx - float(half_width), cx + float(half_width)
        for segment in hole.segments:
            if len(segment.points) < 2:
                continue
            zs = np.asarray(segment.points, dtype=float)[:, 2]
            z_top, z_bot = float(zs[0]), float(zs[-1])
            if z_bot > z_top:
                z_top, z_bot = z_bot, z_top
            label = segment.metadata.get("label") or segment.metadata.get(
                "code"
            ) or ""
            hover = (
                f"{hole.borehole_id}<br>{label}<br>"
                f"{segment.from_md:.1f}-{segment.to_md:.1f} m MD"
                "<br>(patched from the borehole log)<extra></extra>"
            )
            traces.append(
                go.Mesh3d(
                    x=[x0, x1, x1, x0],
                    y=[cy, cy, cy, cy],
                    z=[z_top, z_top, z_bot, z_bot],
                    i=[0, 0],
                    j=[1, 2],
                    k=[2, 3],
                    color=segment.color,
                    opacity=min(1.0, max(0.0, float(opacity))),
                    flatshading=True,
                    name=label or hole.borehole_id,
                    hovertext=hover.replace("<extra></extra>", ""),
                    hoverinfo="text",
                    showlegend=False,
                )
            )
    return traces


def _tube_mesh(points, radius: float, sides: int):
    """Return ``(vertices, faces)`` for a tube along *points*, or None."""
    import numpy as np

    points = np.asarray(points, dtype=float)
    if points.shape[0] < 2:
        return None
    radius = max(float(radius), 1e-6)
    verts: list[Any] = []
    faces: list[tuple[int, int, int]] = []
    for index in range(points.shape[0]):
        if index == 0:
            tangent = points[1] - points[0]
        elif index == points.shape[0] - 1:
            tangent = points[-1] - points[-2]
        else:
            tangent = points[index + 1] - points[index - 1]
        norm = np.linalg.norm(tangent)
        tangent = tangent / norm if norm else np.array([0.0, 0.0, 1.0])
        helper = (
            np.array([1.0, 0.0, 0.0])
            if abs(tangent[0]) < 0.9
            else np.array([0.0, 1.0, 0.0])
        )
        u = np.cross(tangent, helper)
        u /= np.linalg.norm(u) or 1.0
        v = np.cross(tangent, u)
        start = len(verts)
        for side in range(sides):
            angle = 2.0 * np.pi * side / sides
            verts.append(
                points[index]
                + radius * (np.cos(angle) * u + np.sin(angle) * v)
            )
        if index > 0:
            prev = start - sides
            for side in range(sides):
                nxt = (side + 1) % sides
                faces.append((prev + side, prev + nxt, start + side))
                faces.append((prev + nxt, start + nxt, start + side))
    return np.asarray(verts, dtype=float), np.asarray(faces, dtype=int)


# ---------------------------------------------------------------------------
# 2-D basemap collar / target markers
# ---------------------------------------------------------------------------


def document_collar_lonlat(
    document: PCBHDocument | dict[str, Any] | None,
) -> dict[str, tuple[float, float]]:
    """Return ``{borehole_id: (lat, lon)}`` for every locatable collar."""
    if not isinstance(document, PCBHDocument):
        document = document_from_store(document)
    if document is None:
        return {}
    from pycsamt.map.borehole_align import _collar_lonlat

    return _collar_lonlat(document)


def collar_map_markers(
    document: PCBHDocument | dict[str, Any] | None,
    *,
    show_labels: bool = True,
    color: str = "#f59e0b",
    as_target: bool = False,
):
    """Return a basemap scatter trace for borehole collars.

    Placed at the collar WGS84 position. ``as_target`` switches the
    diamond marker for an open down-triangle "drill target" glyph.
    """
    import plotly.graph_objects as go

    if not isinstance(document, PCBHDocument):
        document = document_from_store(document)
    if document is None:
        return []
    lonlat = document_collar_lonlat(document)
    if not lonlat:
        return []
    by_id = {hole.id: hole for hole in document.boreholes}
    lats, lons, texts, hovers = [], [], [], []
    for hid, (lat, lon) in lonlat.items():
        hole = by_id.get(hid)
        lats.append(lat)
        lons.append(lon)
        texts.append(hole.name if hole else hid)
        hovers.append(
            f"{hole.name if hole else hid}<br>{hid}"
            f"<br>{hole.kind if hole else ''}"
            f"<br>TD {hole.total_depth_md:.0f} m"
            if hole
            else hid
        )
    marker = {
        "size": 11,
        "color": color,
        "symbol": "triangle" if as_target else "circle",
    }
    scatter = getattr(go, "Scattermap", None) or go.Scattermapbox
    return [
        scatter(
            lat=lats,
            lon=lons,
            mode="markers+text" if show_labels else "markers",
            text=texts if show_labels else None,
            textposition="top center",
            marker=marker,
            name="Boreholes",
            hovertext=hovers,
            hoverinfo="text",
            showlegend=False,
        )
    ]


def _template(theme: str) -> str:
    return "plotly_dark" if str(theme).lower() == "dark" else "plotly_white"


def _short(text: str, limit: int) -> str:
    text = str(text)
    return text if len(text) <= limit else text[: limit - 1] + "…"
