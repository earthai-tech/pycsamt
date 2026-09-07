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
    "bh": "Boreholes",
    "geology": "Geology",
}


# ── store assembly ─────────────────────────────────────


def store_from_view(view: MapView, *, data_dir: str = "[browsed]") -> dict:
    """Return a JSON-serialisable session store from *view*."""
    df = view.table()
    records = df.to_dict("records")
    line_counts = dict(Counter(s.line or "line" for s in view.data.stations))
    freqs = sorted(
        {float(f) for f in frequency_axis(view.data) if f > 0},
        reverse=True,
    )
    from pycsamt.map import inversion_depth_range

    depth_range = inversion_depth_range(view.data)
    return {
        "station_records": records,
        "n_stations": view.n_stations,
        "n_lines": len(line_counts),
        "line_counts": line_counts,
        "lines": list(view.lines),
        "frequencies": freqs,
        "has_geo": view.has_geo,
        "is_inversion": depth_range is not None,
        "depth_range": list(depth_range) if depth_range else None,
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
    edis = tuple(edi_by_id[sid] for sid in order if sid in edi_by_id)
    metadata = _carried_metadata(old.data.metadata, new.data.metadata)
    data = MapData(
        sites=edis, stations=stations, profiles=(), metadata=metadata
    )
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
        e for e in view.data.iter_edis() if _station_id_from_edi(e) in ids
    )
    metadata = _carried_metadata(view.data.metadata)
    data = MapData(
        sites=edis, stations=tuple(keep), profiles=(), metadata=metadata
    )
    return MapView(data, theme=view.theme, backend=view.backend)


def restrict_to_lines(
    view: MapView,
    active: list[str] | None,
) -> MapView:
    """Return a view limited to *active* lines (or *view* if all on)."""
    if not active:
        return view
    active_set = set(active)
    keep = [s for s in view.data.stations if (s.line or "line") in active_set]
    if len(keep) == len(view.data.stations):
        return view
    ids = {s.id for s in keep}
    edis = tuple(
        e for e in view.data.iter_edis() if _station_id_from_edi(e) in ids
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
        [s.longitude if s.longitude is not None else np.nan for s in stations],
        dtype=float,
    )
    lat = np.array(
        [s.latitude if s.latitude is not None else np.nan for s in stations],
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
        if good[i]
        else s
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
    boreholes: dict | None = None,
    geology: tuple[tuple[float, float, str], ...] | None = None,
    geology_patterns: dict | None = None,
    structure: dict | None = None,
    viewport: dict | None = None,
) -> Any:
    """Build the figure for *view_name* from GUI *controls*.

    ``viewport`` optionally carries the user's last camera / pan-zoom
    keyed by view name (``{"map": {...}, "map3d": {"camera": {...}}}`` —
    see ``pycsamt.app.mapview.callbacks.toolbar._register_viewport``). It
    is replayed onto the freshly-built figure so a control change never
    moves the scene; ``uirevision`` alone is unreliable here because the
    2-D ``map`` (MapLibre) subplot re-applies its ``center``/``zoom`` on
    every rebuild and ``dcc.Loading`` can blank the 3-D camera.

    ``boreholes`` optionally carries ``{"store": <pcbh store>, ...display
    options}``; collars are drawn on the 2-D map and scene-aligned tubes
    on the 3-D views.

    ``geology`` is the applied Interpretation legend as a plain
    ``(rho_min, rho_max, hex_color)`` band list (see
    ``pycsamt.app._geology.geology_bands_from_store``) — ``None`` (the
    default, no legend applied) leaves the continuous colour ramp exactly
    as before. Only the 3-D ``map3d`` views (block/fence/depth/iso) honour
    it; see ``pycsamt.map.config.VolumeMapOptions.geology``.

    ``geology_patterns`` optionally carries ``{legend entry name: ink-
    density stencil array}`` (see ``pycsamt.app._geology.
    geology_pattern_stencils_from_store``); combined with
    ``controls["geology_fill"] == "pattern"`` it texture-fills the
    fence / depth-slice modes instead of a flat colour per band. Block
    / iso-surface always stay solid regardless -- Plotly's Volume /
    Isosurface traces have no per-cell colour-axis override.

    ``structure`` optionally carries ``{"store": <pcgs store>}`` — the
    applied structural model (fault traces, planar/linear measurements),
    inserted into the 3-D scene the same way ``boreholes`` is; see
    ``pycsamt.app._structure.structure_scene_traces``.
    """
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
            depth=c.get("map_depth"),
            cmap=c.get("cmap", "plasma"),
            log_color=bool(c.get("log", False)),
            basemap=c.get("basemap", "esri-satellite"),
            marker_size=int(c.get("marker_size", 10)),
            opacity=float(c.get("map_opacity", 92)) / 100.0,
            show_labels=bool(c.get("labels", True)),
            show_markers=bool(c.get("map_stations", True)),
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
        _add_borehole_collars(fig, boreholes)
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
            crange_percentile=_pct_pair(
                c.get("crange_plo"), c.get("crange_phi")
            ),
            rho_display_max=_pos_or_none(c.get("rho_cutoff")),
            log_color=c.get("scale", "log") == "log",
            topography=bool(c.get("topography", True)),
            show_terrain=bool(c.get("terrain", True)),
            aspectmode=c.get("aspect", "data"),
            x_unit=c.get("x_unit", "m"),
            depth_unit=c.get("depth_unit", "m"),
            smooth_sections=bool(c.get("smooth_sections", True)),
            section_res=int(c.get("section_res", 100)),
            volume_smoothing=float(c.get("volume_smoothing", 0.0) or 0.0),
            show_stations=bool(c.get("show_stations", False)),
            station_labels=bool(c.get("station_labels", False)),
            station_label_angle=float(c.get("station_label_angle", 0.0)),
            station_label_fraction=float(
                c.get("station_label_fraction", 1.0) or 1.0
            ),
            station_label_names=(
                tuple(c["station_label_names"])
                if c.get("station_label_names")
                else None
            ),
            max_stations=c.get("station_max"),
            station_symbol=c.get("station_symbol", "diamond"),
            station_size=int(c.get("station_size", 4)),
            station_color=c.get("station_color", "#1f2937"),
            show_labels=bool(c.get("labels", True)),
            geology=tuple(geology) if geology else None,
            geology_legend=bool(c.get("geology_legend", True)),
            geology_legend_style=c.get("geology_legend_style", "swatch"),
            geology_fill=c.get("geology_fill", "solid"),
            geology_patterns=geology_patterns or None,
        )
        # uirevision keeps the user's camera orbit/zoom across control
        # changes; the 3-D toolbar's "Reset view" bumps the token to
        # snap the camera back to fit the data (mirrors the 2-D "Fit").
        fig.update_layout(uirevision=f"fit-{fit}")
        _add_borehole_scene(fig, view, c, boreholes)
        _add_structure_scene(fig, view, c, structure)
    else:
        return empty_figure(theme, f"Unknown view: {view_name}")
    _apply_viewport(fig, name, viewport)
    return _transparent(fig)


def _apply_viewport(fig: Any, name: str, viewport: dict | None) -> Any:
    """Replay the user's last camera / pan-zoom for *name* onto *fig*.

    A no-op when nothing is stored (fresh load, or just after Fit / Reset
    view, which clear the store). Values are the same ones the user is
    already looking at, so this never produces a visible jump — it only
    stops the rebuilt figure's defaults from winning.
    """
    vp = (viewport or {}).get(name) or {}
    if not vp:
        return fig
    try:
        if name == "map3d":
            camera = vp.get("camera")
            if camera:
                fig.update_layout(scene_camera=camera)
        elif name == "map":
            key = (
                "map"
                if any(
                    getattr(t, "type", "") in ("scattermap", "densitymap")
                    for t in fig.data
                )
                else "mapbox"
            )
            sub = {
                k: vp[k]
                for k in ("center", "zoom", "bearing", "pitch")
                if vp.get(k) is not None
            }
            if sub:
                fig.update_layout(**{key: sub})
        elif name == "pseudosection":
            axes = {}
            if vp.get("xrange") is not None:
                axes["xaxis"] = {"range": vp["xrange"], "autorange": False}
            if vp.get("yrange") is not None:
                axes["yaxis"] = {"range": vp["yrange"], "autorange": False}
            if axes:
                fig.update_layout(**axes)
    except (ValueError, AttributeError, TypeError, KeyError):
        pass
    return fig


# ── borehole overlays (PCBH) ───────────────────────────

_COMPASS_BEARING = {
    "N": 0.0, "NE": 45.0, "E": 90.0, "SE": 135.0,
    "S": 180.0, "SW": 225.0, "W": 270.0, "NW": 315.0,
}


def _hole_point_at_md(hole: Any, md: float, td: float):
    """Scene ``(x, y, z)`` a fraction ``md/td`` along a hole's centerline."""
    pts = np.asarray(hole.centerline, dtype=float)
    if pts.shape[0] < 2 or td <= 0:
        return None
    seg = np.linalg.norm(np.diff(pts, axis=0), axis=1)
    cum = np.concatenate([[0.0], np.cumsum(seg)])
    target = cum[-1] * min(max(md / td, 0.0), 1.0)
    return tuple(
        float(np.interp(target, cum, pts[:, k])) for k in range(3)
    )


def _add_borehole_annotations(
    fig: Any,
    alignment: Any,
    *,
    labels: bool,
    label_angle: float,
    label_size: float,
    depth_ticks: float,
) -> None:
    """Collar labels + measured-depth ticks as rotatable scene annotations.

    ``go.Scatter3d`` text has no rotation, so — exactly like the 3-D
    station labels in :mod:`pycsamt.map.volume` — these ride on
    ``layout.scene.annotations`` instead.
    """
    if not labels and depth_ticks <= 0:
        return
    notes: list[dict] = []
    for hole in alignment.placed:
        cx, cy, cz = hole.collar_scene
        if labels:
            notes.append(
                {
                    "x": cx, "y": cy, "z": cz,
                    "text": hole.name,
                    "showarrow": False,
                    "textangle": float(label_angle),
                    "yanchor": "bottom",
                    "yshift": 10,
                    "font": {"size": float(label_size), "color": "#b45309"},
                }
            )
        if depth_ticks > 0 and hole.centerline:
            td = max((s.to_md for s in hole.segments), default=0.0)
            if td <= 0:
                continue
            mark = depth_ticks
            while mark < td - 1e-6:
                point = _hole_point_at_md(hole, mark, td)
                if point is not None:
                    notes.append(
                        {
                            "x": point[0], "y": point[1], "z": point[2],
                            "text": f"{mark:.0f}",
                            "showarrow": False,
                            "xanchor": "left",
                            "xshift": 5,
                            "font": {
                                "size": max(7.0, float(label_size) - 3.0),
                                "color": "#6b7280",
                            },
                        }
                    )
                mark += depth_ticks
    if not notes:
        return
    existing = list(getattr(fig.layout.scene, "annotations", None) or [])
    fig.update_scenes(annotations=tuple(existing) + tuple(notes))


def _add_borehole_collars(fig: Any, boreholes: dict | None) -> Any:
    """Add basemap collar / target markers to the 2-D map figure."""
    if not boreholes:
        return fig
    labels = bool(boreholes.get("labels", True))
    if boreholes.get("store") and boreholes.get("show_on_map", True):
        try:
            from pycsamt.app._borehole import collar_map_markers

            for trace in collar_map_markers(
                boreholes["store"],
                show_labels=labels,
                as_target=bool(boreholes.get("as_target", False)),
            ):
                fig.add_trace(trace)
        except Exception:  # noqa: BLE001 - never break the map
            pass
    if boreholes.get("points_store") and boreholes.get(
        "points_visible", True
    ):
        try:
            from pycsamt.app._points import point_map_markers

            for trace in point_map_markers(
                boreholes["points_store"], show_labels=labels
            ):
                fig.add_trace(trace)
        except Exception:  # noqa: BLE001
            pass
    return fig


def _add_borehole_scene(
    fig: Any, view: MapView, c: dict, boreholes: dict | None
) -> Any:
    """Add scene-aligned borehole tubes and target points to the 3-D figure."""
    if not boreholes:
        return fig
    want_holes = bool(boreholes.get("store")) and bool(
        boreholes.get("visible", True)
    )
    want_points = bool(boreholes.get("points_store")) and bool(
        boreholes.get("points_visible", True)
    )
    if not want_holes and not want_points:
        return fig
    try:
        from pycsamt.app._borehole import (
            borehole_patch_traces,
            document_from_store,
            pcbh_plotly_traces,
            scene_borehole_traces,
        )
        from pycsamt.map.borehole_align import (
            align_boreholes_to_scene,
            surface_from_sections,
        )
        from pycsamt.map.geometry import survey_frame, survey_uv

        document = (
            document_from_store(boreholes["store"]) if want_holes else None
        )
        family = boreholes.get("family") or "lithology"
        opacity = float(boreholes.get("opacity", 0.9) or 0.9)
        as_tubes = bool(boreholes.get("as_tubes", True))
        labels = bool(boreholes.get("labels", True))
        lean_deg = float(boreholes.get("lean_deg", 0.0) or 0.0)
        lean_azimuth = _COMPASS_BEARING.get(
            str(boreholes.get("lean_dir", "N")).upper(), 0.0
        )
        label_angle = float(boreholes.get("label_angle", 0.0) or 0.0)
        label_size = float(boreholes.get("label_size", 11.0) or 11.0)
        collar_size = float(boreholes.get("collar_size", 5.0) or 5.0)
        depth_ticks = float(boreholes.get("depth_ticks", 0.0) or 0.0)
        patch_geology = bool(boreholes.get("patch_geology", False))
        patch_width = float(boreholes.get("patch_width", 20.0) or 20.0)
        radius = boreholes.get("radius")
        radius_policy = None
        if radius:
            from pycsamt.format.borehole import DisplayRadiusPolicy

            radius_policy = DisplayRadiusPolicy(
                mode="fixed", fixed_radius=float(radius)
            )

        ids, lats, lons, lines, elevs = [], [], [], [], []
        for station in view.data.stations:
            if station.latitude is None or station.longitude is None:
                continue
            ids.append(str(station.id))
            lats.append(float(station.latitude))
            lons.append(float(station.longitude))
            lines.append(station.line or "line")
            elevs.append(
                float(station.elevation)
                if station.elevation is not None
                else np.nan
            )
        frame = survey_frame(lats, lons, lines) if len(ids) >= 2 else None
        if frame is None:
            # No real geometry — fall back to the raw trajectory overlay.
            if document is not None:
                for trace in pcbh_plotly_traces(
                    document,
                    family=family,
                    opacity=opacity,
                    show_labels=labels,
                ):
                    fig.add_trace(trace)
            return fig

        uv = survey_uv(ids, lats, lons, lines)
        us = np.array([uv[i][0] for i in ids if i in uv], dtype=float)
        surface = None
        if bool(c.get("topography", True)) and np.isfinite(elevs).any():
            paired = [
                (uv[i][0], e)
                for i, e in zip(ids, elevs)
                if i in uv and np.isfinite(e)
            ]
            if len(paired) >= 2:
                surface = surface_from_sections(
                    [([p[0] for p in paired], [p[1] for p in paired])]
                )
        datum = (
            "surface"
            if (surface is not None and bool(c.get("topography", True)))
            else "zero"
        )
        depth_lo, depth_hi = c.get("depth_lo"), c.get("depth_hi")
        depth_range = (
            (float(depth_lo), float(depth_hi))
            if depth_lo is not None and depth_hi is not None
            else None
        )
        # The volume builder normalises the line panels: each line sits at
        # ``(median cross-strike v of its stations - front-most line's
        # median) * line_spacing`` (geometry.normalize_offsets +
        # resolve_offset). A borehole collar's raw projected v must go
        # through the *same* shift + stretch or the hole floats in front of
        # / behind its line.  See pycsamt.map.volume._line_offset.
        line_spacing = float(c.get("line_spacing", 1.0) or 1.0)
        line_v: dict[str, list[float]] = {}
        for sid, ln in zip(ids, lines):
            if sid in uv:
                line_v.setdefault(str(ln), []).append(uv[sid][1])
        medians = [float(np.nanmedian(v)) for v in line_v.values() if v]
        offset_shift = float(min(medians)) if medians else 0.0

        scene_bounds = None
        if us.size:
            vs = np.array([uv[i][1] for i in ids if i in uv], dtype=float)
            pad = 0.15 * max(float(np.ptp(us)), 1.0)
            v_lo = (float(vs.min()) - offset_shift) * line_spacing
            v_hi = (float(vs.max()) - offset_shift) * line_spacing
            scene_bounds = (
                float(us.min() - pad),
                float(us.max() + pad),
                min(v_lo, v_hi) - pad,
                max(v_lo, v_hi) + pad,
                -1e6,
                1e6,
            )
        azimuth = float(c.get("azimuth", 0.0) or 0.0)
        if document is not None:
            alignment = align_boreholes_to_scene(
                document,
                frame,
                family=family,
                azimuth_deg=azimuth,
                surface=surface,
                datum=datum,
                depth_range=depth_range,
                scene_bounds=scene_bounds,
                radius_policy=radius_policy,
                offset_shift=offset_shift,
                offset_scale=line_spacing,
                lean_deg=lean_deg,
                lean_azimuth_deg=lean_azimuth,
            )
            # Labels/depth ticks ride on scene.annotations (rotatable,
            # screen-facing) — the same trick volume.py uses for stations;
            # go.Scatter3d text cannot rotate. So suppress the trace text.
            for trace in scene_borehole_traces(
                alignment,
                as_tubes=as_tubes,
                opacity=opacity,
                show_labels=False,
                collar_size=collar_size,
            ):
                fig.add_trace(trace)
            if patch_geology:
                for trace in borehole_patch_traces(
                    alignment, half_width=patch_width / 2.0
                ):
                    fig.add_trace(trace)
            _add_borehole_annotations(
                fig,
                alignment,
                labels=labels,
                label_angle=label_angle,
                label_size=label_size,
                depth_ticks=depth_ticks,
            )
        if want_points:
            from pycsamt.app._points import point_scene_traces

            for trace in point_scene_traces(
                boreholes["points_store"],
                frame,
                surface=surface,
                datum=datum,
                azimuth_deg=azimuth,
                offset_shift=offset_shift,
                offset_scale=line_spacing,
                show_labels=labels,
            ):
                fig.add_trace(trace)
    except Exception:  # noqa: BLE001 - never break the scene for an overlay
        pass
    return fig


def _add_structure_scene(fig: Any, view: MapView, c: dict, structure: dict | None) -> Any:
    """Add applied structural-geology traces (fault planes, planar/linear
    measurement glyphs) to the 3-D figure — the same per-line offset math
    :func:`_add_borehole_scene` uses, but simpler: a structural item's
    ``x`` is already a profile position on its own line (like a fence
    panel's own along-strike axis), not a real-world lat/lon needing
    reprojection.
    """
    if not structure or not structure.get("store"):
        return fig
    try:
        from pycsamt.app._structure import (
            structure_from_store,
            structure_scene_traces,
        )
        from pycsamt.map.borehole_align import surface_from_sections
        from pycsamt.map.geometry import survey_frame, survey_uv

        doc = structure_from_store(structure["store"])
        if doc is None or len(doc) == 0:
            return fig

        ids, lats, lons, lines, elevs = [], [], [], [], []
        for station in view.data.stations:
            if station.latitude is None or station.longitude is None:
                continue
            ids.append(str(station.id))
            lats.append(float(station.latitude))
            lons.append(float(station.longitude))
            lines.append(station.line or "line")
            elevs.append(
                float(station.elevation)
                if station.elevation is not None
                else np.nan
            )
        if len(ids) < 2:
            return fig
        frame = survey_frame(lats, lons, lines)
        uv = survey_uv(ids, lats, lons, lines)
        surface = None
        if bool(c.get("topography", True)) and np.isfinite(elevs).any():
            paired = [
                (uv[i][0], e)
                for i, e in zip(ids, elevs)
                if i in uv and np.isfinite(e)
            ]
            if len(paired) >= 2:
                surface = surface_from_sections(
                    [([p[0] for p in paired], [p[1] for p in paired])]
                )
        # Same per-line offset normalisation the volume builder and
        # _add_borehole_scene use, see pycsamt.map.volume._line_offset.
        line_spacing = float(c.get("line_spacing", 1.0) or 1.0)
        line_v: dict[str, list[float]] = {}
        for sid, ln in zip(ids, lines):
            if sid in uv:
                line_v.setdefault(str(ln), []).append(uv[sid][1])
        medians = {
            ln: float(np.nanmedian(v)) for ln, v in line_v.items() if v
        }
        if not medians:
            return fig
        offset_shift = min(medians.values())
        line_offsets = {
            ln: (v - offset_shift) * line_spacing
            for ln, v in medians.items()
        }
        default_line = next(iter(line_offsets))
        azimuth = float(c.get("azimuth", 0.0) or 0.0)
        depth_lo, depth_hi = c.get("depth_lo"), c.get("depth_hi")
        depth_extent = (
            float(depth_hi) - float(depth_lo)
            if depth_lo is not None and depth_hi is not None
            else 300.0
        )
        for trace in structure_scene_traces(
            doc.model,
            line_offsets=line_offsets,
            default_line=default_line,
            azimuth_deg=azimuth,
            surface=surface,
            depth_extent=max(depth_extent, 10.0),
        ):
            fig.add_trace(trace)
    except Exception:  # noqa: BLE001 - never break the scene for an overlay
        pass
    return fig


def _pair(lo, hi):
    """Return ``(lo, hi)`` if both finite and ordered, else ``None``."""
    try:
        lo_f, hi_f = float(lo), float(hi)
    except (TypeError, ValueError):
        return None
    if hi_f > lo_f >= 0:
        return (lo_f, hi_f)
    return None


def _pct_pair(lo, hi):
    """Percentile ``(lo, hi)`` in ``[0, 100]`` and ordered, else ``None``.

    ``None`` restores the raw min/max auto colour range (see
    :attr:`pycsamt.map.config.VolumeMapOptions.crange_percentile`).
    """
    try:
        lo_f, hi_f = float(lo), float(hi)
    except (TypeError, ValueError):
        return None
    if 0.0 <= lo_f < hi_f <= 100.0:
        return (lo_f, hi_f)
    return None


def _pos_or_none(value):
    """Return ``float(value)`` when it is a positive number, else ``None``."""
    try:
        v = float(value)
    except (TypeError, ValueError):
        return None
    return v if v > 0 else None


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
