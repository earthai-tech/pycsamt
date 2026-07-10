# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Bridge between :class:`pycsamt.map.MapView` and the Dash figures.

All plotting goes through ``MapView`` so the platform owns no plotting
logic of its own; these helpers only translate GUI control state into
``MapView`` calls and assemble the serialisable session store.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import replace
from typing import Any

import numpy as np

from pycsamt.map import MapView
from pycsamt.map._core import (
    MapData,
    _station_id_from_edi,
    frequency_axis,
)

VIEW_TITLES = {
    "map": "Map view",
    "pseudosection": "Pseudosection",
    "map3d": "3-D map",
}


# ── store assembly ─────────────────────────────────────


def store_from_view(view: MapView, *, data_dir: str = "[browsed]") -> dict:
    """Return a JSON-serialisable session store from *view*."""
    df = view.table()
    records = df.to_dict("records")
    line_counts = dict(
        Counter(s.line or "line" for s in view.data.stations)
    )
    freqs = sorted(
        {float(f) for f in frequency_axis(view.data) if f > 0},
        reverse=True,
    )
    return {
        "station_records": records,
        "n_stations": view.n_stations,
        "n_lines": len(line_counts),
        "line_counts": line_counts,
        "lines": list(view.lines),
        "frequencies": freqs,
        "has_geo": view.has_geo,
        "data_dir": data_dir,
    }


# ── view algebra ───────────────────────────────────────


def _carried_metadata(*sources: dict | None) -> dict:
    """Metadata to carry onto a derived MapData, later sources winning.

    Drops ``n_stations``/``n_profiles`` so :class:`MapData.__post_init__`
    recomputes them for the new station set, but preserves everything
    else — notably ``sections`` (precomputed inversion curtains, see
    :mod:`pycsamt.map.inversion`), which would otherwise silently
    disappear whenever a view is filtered, masked, or merged.
    """
    merged: dict[str, Any] = {}
    sections: dict[str, Any] = {}
    for meta in sources:
        for key, value in (meta or {}).items():
            if key in ("n_stations", "n_profiles"):
                continue
            if key == "sections" and isinstance(value, dict):
                sections.update(value)
                continue
            merged[key] = value
    if sections:
        merged["sections"] = sections
    return merged


def merge_views(old: MapView, new: MapView) -> MapView:
    """Append *new* into *old*; new stations win on ID collision."""
    by_id: dict[str, Any] = {}
    order: list[str] = []
    for station in (*old.data.stations, *new.data.stations):
        if station.id not in by_id:
            order.append(station.id)
        by_id[station.id] = station
    stations = tuple(
        replace(by_id[sid], index=i) for i, sid in enumerate(order)
    )
    edi_by_id: dict[str, Any] = {}
    for edi in (*old.data.iter_edis(), *new.data.iter_edis()):
        edi_by_id[_station_id_from_edi(edi)] = edi
    edis = tuple(
        edi_by_id[sid] for sid in order if sid in edi_by_id
    )
    metadata = _carried_metadata(old.data.metadata, new.data.metadata)
    data = MapData(sites=edis, stations=stations, profiles=(), metadata=metadata)
    return MapView(data, theme=old.theme, backend=old.backend)


def apply_settings(view, active_lines, masked):
    """Restrict a view to active lines and drop masked station ids."""
    view = restrict_to_lines(view, active_lines)
    if masked:
        view = exclude_stations(view, masked)
    return view


def exclude_stations(view: MapView, masked) -> MapView:
    """Return a view with the given station ids removed."""
    mset = {str(m) for m in (masked or [])}
    if not mset:
        return view
    keep = [s for s in view.data.stations if s.id not in mset]
    if len(keep) == len(view.data.stations):
        return view
    ids = {s.id for s in keep}
    edis = tuple(
        e for e in view.data.iter_edis()
        if _station_id_from_edi(e) in ids
    )
    metadata = _carried_metadata(view.data.metadata)
    data = MapData(sites=edis, stations=tuple(keep), profiles=(), metadata=metadata)
    return MapView(data, theme=view.theme, backend=view.backend)


def restrict_to_lines(
    view: MapView,
    active: list[str] | None,
) -> MapView:
    """Return a view limited to *active* lines (or *view* if all on)."""
    if not active:
        return view
    active_set = set(active)
    keep = [
        s
        for s in view.data.stations
        if (s.line or "line") in active_set
    ]
    if len(keep) == len(view.data.stations):
        return view
    ids = {s.id for s in keep}
    edis = tuple(
        e for e in view.data.iter_edis()
        if _station_id_from_edi(e) in ids
    )
    data = MapData(
        sites=edis,
        stations=tuple(keep),
        profiles=(),
        metadata=_carried_metadata(view.data.metadata),
    )
    return MapView(data, theme=view.theme, backend=view.backend)


# ── figure dispatch ────────────────────────────────────


def reproject_view(view, mode, zone, hem, epsg):
    """Interpret station coords in *mode*'s CRS and reproject to lon/lat.

    ``geo`` (or EPSG:4326) returns the view unchanged. ``utm`` maps the
    zone/hemisphere to EPSG 326xx/327xx; ``custom`` uses *epsg*.
    """
    if not mode or mode == "geo":
        return view
    code = _source_epsg(mode, zone, epsg, hem)
    if code in (4326, None):
        return view
    try:
        from pycsamt.map.overlays import (
            CRSConfig,
            transform_xy,
        )
    except Exception:
        return view

    stations = view.data.stations
    lon = np.array(
        [s.longitude if s.longitude is not None else np.nan
         for s in stations],
        dtype=float,
    )
    lat = np.array(
        [s.latitude if s.latitude is not None else np.nan
         for s in stations],
        dtype=float,
    )
    good = np.isfinite(lon) & np.isfinite(lat)
    if not good.any():
        return view
    try:
        xt, yt = transform_xy(
            lon, lat, crs=CRSConfig(source=code, target=4326)
        )
    except Exception:
        return view
    new_stations = tuple(
        replace(s, longitude=float(xt[i]), latitude=float(yt[i]))
        if good[i] else s
        for i, s in enumerate(stations)
    )
    data = MapData(
        sites=view.data.sites,
        stations=new_stations,
        profiles=(),
        metadata=dict(view.data.metadata),
    )
    return MapView(data, theme=view.theme, backend=view.backend)


def project_to_crs(lons, lats, mode, zone, hem, epsg):
    """Convert geographic lon/lat → the display CRS (easting, northing).

    Returns ``(east, north, epsg_code)`` as arrays, or ``(None, None,
    code)`` for geographic mode / on failure.
    """
    code = _source_epsg(mode, zone, epsg, hem)
    if not mode or mode == "geo" or code in (4326, None):
        return None, None, code
    try:
        from pycsamt.map.overlays import (
            CRSConfig,
            transform_xy,
        )

        east, north = transform_xy(
            np.asarray(lons, dtype=float),
            np.asarray(lats, dtype=float),
            crs=CRSConfig(source=4326, target=code),
        )
        return east, north, code
    except Exception:
        return None, None, code


def _source_epsg(mode, zone, epsg, hem):
    if mode == "utm":
        try:
            z = int(zone or 50)
        except (TypeError, ValueError):
            z = 50
        return (32600 + z) if str(hem or "N").upper() == "N" else (32700 + z)
    try:
        return int(epsg or 4326)
    except (TypeError, ValueError):
        return 4326


def figure_for(
    view_name: str,
    view: MapView | None,
    controls: dict | None,
    *,
    theme: str = "light",
    active_lines: list[str] | None = None,
    masked: list[str] | None = None,
    fit: int = 0,
) -> Any:
    """Build the figure for *view_name* from GUI *controls*."""
    if view is None or view.n_stations == 0:
        return empty_figure(theme)
    c = controls or {}
    # The basemap always uses the survey's geographic lon/lat. The CRS
    # panel is a *display* conversion (station inspector + table), not a
    # map reprojection — so we don't reproject the view here.
    view = apply_settings(view, active_lines, masked)
    if view.n_stations == 0:
        return empty_figure(theme, "All stations are hidden or masked.")
    name = (view_name or "map").lower()

    if name == "map":
        fig = view.station(
            theme=theme,
            overlay=c.get("overlay", "index"),
            component=c.get("component", "xy"),
            frequency=c.get("frequency"),
            cmap=c.get("cmap", "plasma"),
            log_color=bool(c.get("log", False)),
            basemap=c.get("basemap", "esri-satellite"),
            marker_size=int(c.get("marker_size", 10)),
            opacity=float(c.get("map_opacity", 92)) / 100.0,
            show_labels=bool(c.get("labels", True)),
            show_profiles=bool(c.get("profiles", True)),
            contour_image=bool(c.get("contour_enable", False)),
            contour_levels=int(c.get("contour_levels", 12)),
            contour_mode=c.get("contour_mode", "filled+lines"),
            contour_interp=c.get("contour_interp", "cubic"),
            contour_smooth=float(c.get("contour_smooth", 1.0)),
            contour_grid_res=int(c.get("contour_res", 150)),
        )
        # uirevision keeps the user's pan/zoom across control changes;
        # the Fit button bumps the token to re-fit to the data.
        fig.update_layout(uirevision=f"fit-{fit}")
    elif name == "pseudosection":
        fig = view.pseudosection(
            theme=theme,
            component=c.get("component", "xy"),
            quantity=c.get("quantity", "rho"),
            log_rho=bool(c.get("log", True)),
            cmap=c.get("cmap"),
        )
    elif name == "map3d":
        fig = view.map3d(
            theme=theme,
            mode=c.get("mode3d", "fence"),
            component=c.get("component", "xy"),
            quantity=c.get("quantity", "resistivity"),
            cmap=c.get("cmap", "RdYlBu_r"),
            opacity=float(c.get("opacity", 0.85)),
            azimuth=float(c.get("azimuth", 0.0)),
            line_spacing=float(c.get("line_spacing", 1.0)),
            n_slices=int(c.get("n_slices", 8)),
            surface_count=int(c.get("surface_count", 12)),
            show_contours=bool(c.get("contours", False)),
            depth_range=_pair(c.get("depth_lo"), c.get("depth_hi")),
            value_range=_pair(c.get("vmin"), c.get("vmax")),
            rho_range=_pair(c.get("rho_lo"), c.get("rho_hi")),
            log_color=c.get("scale", "log") == "log",
            topography=bool(c.get("topography", True)),
            show_terrain=bool(c.get("terrain", True)),
            aspectmode=c.get("aspect", "data"),
            x_unit=c.get("x_unit", "m"),
            depth_unit=c.get("depth_unit", "m"),
            smooth_sections=bool(c.get("smooth_sections", True)),
            section_res=int(c.get("section_res", 100)),
            show_stations=bool(c.get("show_stations", False)),
            station_labels=bool(c.get("station_labels", False)),
            station_symbol=c.get("station_symbol", "diamond"),
            station_size=int(c.get("station_size", 4)),
            station_color=c.get("station_color", "#1f2937"),
            show_labels=bool(c.get("labels", True)),
        )
        # uirevision keeps the user's camera orbit/zoom across control
        # changes; the 3-D toolbar's "Reset view" bumps the token to
        # snap the camera back to fit the data (mirrors the 2-D "Fit").
        fig.update_layout(uirevision=f"fit-{fit}")
    else:
        return empty_figure(theme, f"Unknown view: {view_name}")
    return _transparent(fig)


def _pair(lo, hi):
    """Return ``(lo, hi)`` if both finite and ordered, else ``None``."""
    try:
        lo_f, hi_f = float(lo), float(hi)
    except (TypeError, ValueError):
        return None
    if hi_f > lo_f >= 0:
        return (lo_f, hi_f)
    return None


def _transparent(fig):
    """Drop the solid plot/paper background so the panel shows through."""
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
    )
    try:  # 3-D scenes: clear the grey wall panes too
        clear = dict(
            showbackground=False,
            backgroundcolor="rgba(0,0,0,0)",
        )
        fig.update_scenes(
            bgcolor="rgba(0,0,0,0)",
            xaxis=clear,
            yaxis=clear,
            zaxis=clear,
        )
    except (ValueError, AttributeError):
        pass
    return fig


def empty_figure(theme: str = "light", msg: str = "Load EDI lines to begin"):
    """Return a themed placeholder figure."""
    import plotly.graph_objects as go

    from pycsamt.map.styles import theme_colors

    colors = theme_colors(theme)
    fig = go.Figure()
    fig.add_annotation(
        text=msg,
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color=colors["text"], size=15),
    )
    fig.update_layout(
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=0, r=0, t=10, b=0),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )
    return fig
