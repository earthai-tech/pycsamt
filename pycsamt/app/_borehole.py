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
    pcbh_from_dict,
    pcbh_to_dict,
)

__all__ = [
    "decode_pcbh_upload",
    "document_from_store",
    "pcbh_plotly_traces",
    "add_pcbh_to_figure",
    "add_embedded_pcbh_to_pcsf_figure",
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
