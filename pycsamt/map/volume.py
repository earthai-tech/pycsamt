# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""3-D map API for fence, block, and depth-slice views."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from ._backends import require_plotly
from ._core import (
    MapData,
    ensure_map_data,
    pseudosection_table,
)
from .config import VolumeMapOptions
from .geometry import (
    normalize_offsets,
    resolve_offset,
    survey_uv,
)
from .styles import theme_colors, to_plotly_cmap


class Map3D:
    """Builder object for 3-D survey maps."""

    def __init__(
        self,
        data: Any,
        *,
        options: VolumeMapOptions | None = None,
        **ensure_kwargs: Any,
    ) -> None:
        self.data: MapData = ensure_map_data(
            data,
            **ensure_kwargs,
        )
        self.options = options or VolumeMapOptions()

    def with_mode(self, mode: str) -> Map3D:
        """Return a copy with another 3-D mode."""
        opts = replace(self.options, mode=mode)  # type: ignore[arg-type]
        new = object.__new__(Map3D)
        new.data = self.data
        new.options = opts
        return new

    def with_options(self, **kwargs: Any) -> Map3D:
        """Return a copy with updated options."""
        new = object.__new__(Map3D)
        new.data = self.data
        new.options = replace(self.options, **kwargs)
        return new

    def with_quantity(self, quantity: str) -> Map3D:
        """Return a copy using another mapped quantity."""
        return self.with_options(quantity=quantity)

    def with_component(self, component: str) -> Map3D:
        """Return a copy using another component."""
        return self.with_options(component=component)

    def figure(self) -> Any:
        """Build and return the 3-D map figure."""
        return build_3d_map(self.data, self.options)


VolumeMap = Map3D


def plot_3d_map(
    data: Any,
    *,
    options: VolumeMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a 3-D map."""
    return Map3D(data, options=options, **kwargs).figure()


def plot_volume_map(
    data: Any,
    *,
    options: VolumeMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a 3-D volume map."""
    return plot_3d_map(data, options=options, **kwargs)


def build_3d_map(
    data: MapData,
    options: VolumeMapOptions,
) -> Any:
    """Build the concrete 3-D figure."""
    profiles = _profile_grids(data, options)
    colors = theme_colors(options.theme)
    if not profiles:
        return _empty_3d_figure(colors)
    builders = {
        "block": _block_figure,
        "depth": _depth_figure,
        "surface": _surface_figure,
        "fence": _fence_figure,
    }
    builder = builders.get(options.mode)
    if builder is None:
        msg = f"Unknown 3-D map mode: {options.mode}"
        raise ValueError(msg)
    fig = builder(profiles, options, colors)
    if options.show_stations:
        _add_station_markers(fig, profiles, options)
    _apply_axis_units(fig, options)
    return fig


def _unit_scale(unit) -> float:
    """Metres -> display factor ('m' -> 1, 'km' -> 1e-3)."""
    return 0.001 if str(unit).lower() == "km" else 1.0


def _apply_axis_units(fig, options) -> None:
    """Rescale trace coordinates from metres to the chosen display units.

    Profile (x) and line-offset (y) use ``x_unit``; depth (z) uses
    ``depth_unit``.  The axes are independent, so x can stay in metres
    while depth is shown in km without squashing the scene.
    """
    sx = _unit_scale(getattr(options, "x_unit", "m"))
    sz = _unit_scale(getattr(options, "depth_unit", "m"))
    if sx == 1.0 and sz == 1.0:
        return
    for trace in fig.data:
        for attr, scale in (("x", sx), ("y", sx), ("z", sz)):
            values = getattr(trace, attr, None)
            if values is None:
                continue
            arr = np.asarray(values, dtype=float)
            if arr.size:
                setattr(trace, attr, arr * scale)


def _profile_grids(
    data: MapData,
    options: VolumeMapOptions,
) -> dict[str, dict[str, np.ndarray]]:
    """Per-line (x, z, rho, ...) grids for the 3-D builders.

    Merges two sources: lines with a precomputed section (e.g. a
    ModEM 3-D volume already sliced per line — see
    :mod:`pycsamt.map.inversion`) are used directly; any other line
    falls back to the EDI pseudosection/skin-depth path. Both share
    one real-geometry projection (:func:`_dataset_uv`) so mixed
    EDI + inversion multi-line views stay consistently placed.
    """
    uv = _dataset_uv(data)
    out = _edi_profile_grids(data, options, uv)
    out.update(_inversion_profile_grids(data, options, uv))
    return out


def _edi_profile_grids(
    data: MapData,
    options: VolumeMapOptions,
    uv: dict[str, tuple[float, float]],
) -> dict[str, dict[str, np.ndarray]]:
    quantity = _volume_quantity(options.quantity)
    value_table = pseudosection_table(
        data,
        quantity=quantity,
        component=options.component,
    )
    rho_table = pseudosection_table(
        data,
        quantity="rho",
        component=options.component,
    )
    if value_table.empty or rho_table.empty:
        return {}
    station_line = {
        station.id: station.line or "line"
        for station in data.stations
    }
    station_elev = {
        station.id: (
            float(station.elevation)
            if station.elevation is not None
            else np.nan
        )
        for station in data.stations
    }
    value_table = _prepare_volume_table(
        value_table,
        station_line,
        options,
    )
    rho_table = _prepare_volume_table(
        rho_table,
        station_line,
        options,
    )
    out: dict[str, dict[str, np.ndarray]] = {}
    for line, group in value_table.groupby("line"):
        rho_group = rho_table[rho_table["line"] == line]
        if rho_group.empty:
            continue
        value_piv = group.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()
        rho_piv = rho_group.pivot_table(
            index="period",
            columns="station",
            values="value",
            aggfunc="median",
        ).sort_index()
        rho_piv = rho_piv.reindex_like(value_piv)
        if value_piv.empty or rho_piv.empty:
            continue
        periods = value_piv.index.to_numpy(dtype=float)
        values = value_piv.to_numpy(dtype=float)
        rho = rho_piv.to_numpy(dtype=float)
        rho = np.where(rho > 0, rho, np.nan)
        median_rho = np.nanmedian(rho, axis=1)
        depth = 503.0 * np.sqrt(median_rho * periods)
        depth = np.where(np.isfinite(depth), depth, periods)
        if options.depth_range:
            lo, hi = options.depth_range
            keep = (depth >= lo) & (depth <= hi)
            depth = depth[keep]
            rho = rho[keep, :]
            values = values[keep, :]
            periods = periods[keep]
        if values.size == 0:
            continue
        stations = np.array(
            list(value_piv.columns),
            dtype=object,
        )
        x, y_offset = _station_uv(group, stations, uv)
        out[str(line)] = {
            "x": x,
            "y_offset": y_offset,
            "z": depth,
            "period": periods,
            "rho": rho,
            "value": values,
            "elev": np.array(
                [station_elev.get(str(s), np.nan) for s in stations],
                dtype=float,
            ),
            "quantity": np.array([quantity], dtype=object),
            "stations": stations,
        }
    return out


def _inversion_profile_grids(
    data: MapData,
    options: VolumeMapOptions,
    uv: dict[str, tuple[float, float]],
) -> dict[str, dict[str, np.ndarray]]:
    """Per-line grids for lines with a precomputed section.

    Consumes ``data.metadata["sections"]`` (see
    :func:`pycsamt.map.inversion.load_modem_lines`) directly instead
    of the EDI pseudosection/skin-depth path — an inversion result is
    already a real ``(x, z, rho)`` section, not raw impedance spectra
    to be gridded.
    """
    sections = (data.metadata or {}).get("sections")
    if not sections:
        return {}
    # A line survives (via merge/filter/mask) in ``metadata["sections"]``
    # even after its stations are removed from ``data.stations`` — skip
    # it here so a filtered-out or masked-away line doesn't keep
    # rendering its precomputed section.
    active_lines = {s.line for s in data.stations if s.line}
    # Elevation is read live from ``data.stations`` (falling back to the
    # section's own snapshot for stations no longer present there) so
    # MapView.with_elevations()/the topo-apply flow — which only update
    # station records, not the precomputed sections — still take effect.
    station_elev = {
        s.id: s.elevation for s in data.stations if s.elevation is not None
    }
    quantity = _volume_quantity(options.quantity)
    out: dict[str, dict[str, np.ndarray]] = {}
    for line, section in sections.items():
        if str(line) not in active_lines:
            continue
        stations = np.asarray(section["stations"], dtype=object)
        rho = np.asarray(section["rho"], dtype=float)
        z = np.asarray(section["z"], dtype=float)
        if stations.size == 0 or rho.size == 0:
            continue
        if options.depth_range:
            lo, hi = options.depth_range
            keep = (z >= lo) & (z <= hi)
            z = z[keep]
            rho = rho[keep, :]
        if rho.size == 0:
            continue
        x, y_offset = _section_uv(stations, uv)
        section_elev = section.get("elev", np.full(stations.size, np.nan))
        elev = np.array(
            [
                station_elev.get(str(sid), section_elev[i])
                for i, sid in enumerate(stations)
            ],
            dtype=float,
        )
        out[str(line)] = {
            "x": x,
            "y_offset": y_offset,
            "z": z,
            "period": np.full_like(z, np.nan),
            "rho": rho,
            "value": rho,
            "elev": elev,
            "quantity": np.array([quantity], dtype=object),
            "stations": stations,
        }
    return out


def _fence_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    az = np.deg2rad(float(options.azimuth))
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    cmin, cmax = _crange(options)
    for idx, (name, grid) in enumerate(profiles.items()):
        x, z_pos, values, elev = _prepare_section(grid, options)
        z = -z_pos
        xx, zz = np.meshgrid(x, z)
        offset = _line_offset(name, idx, real_offsets, unit, options)
        yy = np.full_like(
            xx,
            offset * np.cos(az),
        )
        xx = xx + offset * np.sin(az)
        if options.topography:
            zz = zz + elev[np.newaxis, :]
        fig.add_trace(
            go.Surface(
                x=xx,
                y=yy,
                z=zz,
                surfacecolor=_color_values(values, options),
                colorscale=to_plotly_cmap(options.cmap),
                cmin=cmin,
                cmax=cmax,
                opacity=float(options.opacity),
                name=name,
                showscale=idx == 0,
                colorbar=dict(title=dict(text=_colorbar_title(options), side="right")),
                contours=_surface_contours(options),
            )
        )
        if options.show_terrain and options.topography:
            fig.add_trace(
                _terrain_trace(xx[0], yy[0], elev, name)
            )
        if options.show_labels:
            fig.add_trace(
                go.Scatter3d(
                    x=[float(np.nanmean(xx))],
                    y=[float(np.nanmean(yy))],
                    z=[0.0],
                    text=[name],
                    mode="text",
                    showlegend=False,
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _prepare_section(grid, options):
    """Return ``(x, depth, values, elev)`` for one line's fence panel.

    Coordinates are sorted/de-duplicated so the section is a proper
    monotonic grid, then optionally resampled onto a denser regular
    grid (cubic/linear spline) so the panel reads as a smooth section
    instead of a handful of raw per-station stripes — matching the
    interpolated volume the web app's 3-D map view renders.
    """
    x = np.asarray(grid["x"], dtype=float)
    z = np.asarray(grid["z"], dtype=float)
    values = _filtered_values(grid, options)
    elev = _elev_for(grid, options)

    order = np.argsort(x)
    x, values, elev = x[order], values[:, order], elev[order]
    x, uniq = np.unique(x, return_index=True)
    values, elev = values[:, uniq], elev[uniq]

    zorder = np.argsort(z)
    z, values = z[zorder], values[zorder, :]
    z, zuniq = np.unique(z, return_index=True)
    values = values[zuniq, :]

    if not getattr(options, "smooth_sections", True):
        return x, z, values, elev
    if x.size < 3 or z.size < 3:
        return x, z, values, elev

    filled = _fill_nan_2d(z, x, values)
    if filled is None:
        return x, z, values, elev
    try:
        from scipy.interpolate import RectBivariateSpline

        res = max(int(getattr(options, "section_res", 100)), 2)
        kx = min(3, z.size - 1)
        ky = min(3, x.size - 1)
        spline = RectBivariateSpline(z, x, filled, kx=kx, ky=ky)
        xi = np.linspace(x.min(), x.max(), max(res, x.size))
        zi = np.linspace(z.min(), z.max(), max(res, z.size))
        vi = spline(zi, xi)
        # cubic splines can overshoot past the source range near sharp
        # gradients; clip back so the colorscale isn't stretched by
        # spline artefacts rather than real data.
        vi = np.clip(vi, float(np.nanmin(filled)), float(np.nanmax(filled)))
        elev_i = np.interp(xi, x, elev)
        return xi, zi, vi, elev_i
    except Exception:  # noqa: BLE001 - fall back to the raw grid
        return x, z, values, elev


def _fill_nan_2d(z, x, values):
    """Fill NaN gaps with nearest-neighbour values so a spline can fit."""
    values = np.asarray(values, dtype=float)
    good = np.isfinite(values)
    if good.sum() < 3:
        return None
    if good.all():
        return values
    zz, xx = np.meshgrid(z, x, indexing="ij")
    try:
        from scipy.interpolate import griddata

        filled = griddata(
            (zz[good], xx[good]), values[good],
            (zz, xx), method="nearest",
        )
        return np.where(good, values, filled)
    except Exception:  # noqa: BLE001
        return np.where(good, values, float(np.nanmedian(values[good])))


def _block_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    grid = _dense_volume_grid(profiles, options)
    if grid is not None:
        x_arr, y_arr, z_arr, rho_vol = grid
        X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing="ij")
        finite = rho_vol[np.isfinite(rho_vol)]
        if finite.size:
            iso_lo, iso_hi = _iso_range(finite, options)
            cmin, cmax = _crange(options)
            if cmin is None or cmax is None:
                cmin, cmax = float(finite.min()), float(finite.max())
            fig.add_trace(
                go.Volume(
                    x=X.ravel(),
                    y=Y.ravel(),
                    z=Z.ravel(),
                    value=rho_vol.ravel(),
                    isomin=iso_lo,
                    isomax=iso_hi,
                    # Fade rather than hard-cut so the block still reads
                    # as one solid shape instead of a jagged threshold.
                    opacityscale=[
                        [0.0, 0.0], [0.2, 0.3], [0.5, 0.7], [1.0, 1.0],
                    ],
                    opacity=float(options.opacity),
                    surface_count=max(2, int(options.surface_count)),
                    colorscale=to_plotly_cmap(options.cmap),
                    cmin=cmin,
                    cmax=cmax,
                    showscale=True,
                    colorbar=dict(title=dict(text=_colorbar_title(options), side="right")),
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _depth_figure(profiles, options, colors):
    go = require_plotly()

    depths = _slice_depths(profiles, options)
    az = np.deg2rad(float(options.azimuth))
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    fig = go.Figure()
    # Lines may hold different station counts; pad every row
    # to the widest so the per-depth surface grid is
    # rectangular. NaN gaps render as holes in Plotly.
    width = max(
        (np.asarray(g["x"]).size for g in profiles.values()),
        default=0,
    )
    if width == 0:
        return _empty_3d_figure(colors)
    cmin, cmax = _crange(options)
    for depth in depths:
        x_rows = []
        y_rows = []
        z_rows = []
        val_rows = []
        for line_idx, (name, grid) in enumerate(profiles.items()):
            values = _values_at_depth(
                grid,
                depth,
                options,
            )
            x = np.asarray(grid["x"], dtype=float)
            offset = _line_offset(name, line_idx, real_offsets, unit, options)
            y = np.full_like(
                x,
                offset * np.cos(az),
            )
            x_rows.append(
                _pad_row(x + offset * np.sin(az), width)
            )
            y_rows.append(_pad_row(y, width))
            z = _elev_for(grid, options) - float(depth)
            z_rows.append(
                _pad_row(z, width)
            )
            color = _color_values(values, options)
            val_rows.append(
                _pad_row(color, width)
            )
        fig.add_trace(
            go.Surface(
                x=np.vstack(x_rows),
                y=np.vstack(y_rows),
                z=np.vstack(z_rows),
                surfacecolor=np.vstack(val_rows),
                colorscale=to_plotly_cmap(options.cmap),
                cmin=cmin,
                cmax=cmax,
                opacity=float(options.opacity),
                showscale=True,
                colorbar=dict(title=dict(text=_colorbar_title(options), side="right")),
                contours=_surface_contours(options),
            )
        )
    if options.show_terrain and options.topography:
        for line_idx, (name, grid) in enumerate(profiles.items()):
            x = np.asarray(grid["x"], dtype=float)
            offset = _line_offset(name, line_idx, real_offsets, unit, options)
            fig.add_trace(
                _terrain_trace(
                    x + offset * np.sin(az),
                    np.full_like(x, offset * np.cos(az)),
                    _elev_for(grid, options),
                    name,
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _surface_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    grid = _dense_volume_grid(profiles, options)
    if grid is not None:
        x_arr, y_arr, z_arr, rho_vol = grid
        X, Y, Z = np.meshgrid(x_arr, y_arr, z_arr, indexing="ij")
        finite = rho_vol[np.isfinite(rho_vol)]
        if finite.size:
            lo, hi = _iso_range(finite, options)
            fig.add_trace(
                go.Isosurface(
                    x=X.ravel(),
                    y=Y.ravel(),
                    z=Z.ravel(),
                    value=rho_vol.ravel(),
                    isomin=lo,
                    isomax=hi,
                    surface_count=max(2, int(options.surface_count)),
                    opacity=float(options.opacity),
                    colorscale=to_plotly_cmap(options.cmap),
                    caps=dict(x_show=False, y_show=False),
                    colorbar=dict(title=dict(text=_colorbar_title(options), side="right")),
                    cmin=lo,
                    cmax=hi,
                )
            )
    _style_3d(fig, options, colors)
    return fig


def _prepare_volume_table(table, station_line, options):
    out = table.copy()
    out["line"] = (
        out["station"].map(station_line).fillna("line")
    )
    if options.period_range:
        lo, hi = options.period_range
        keep = (
            (out["period"] >= lo)
            & (out["period"] <= hi)
        )
        out = out[keep]
    return out


def _station_x(group, stations):
    x_map = (
        group.groupby("station")["distance"]
        .median()
        .to_dict()
    )
    x = np.array(
        [
            x_map.get(str(station), np.nan)
            for station in stations
        ],
        dtype=float,
    )
    if not np.isfinite(x).all():
        # No real geometry: fall back to an index spaced ~100 m apart so
        # the profile axis stays comparable to depth (metres), not 0,1,2…
        return np.arange(len(stations), dtype=float) * 100.0
    # ``distance`` is in km; the 3-D depth axis is in metres — convert so
    # the profile (x) and depth (z) axes share a unit and render to scale.
    return x * 1000.0


def _line_offset_unit(profiles) -> float:
    """Base Y-offset unit for fence/depth lines, in the same (metre)
    unit as the profile axis, so lines separate proportionally.

    Only used as a fallback when stations carry no real coordinates —
    see :func:`_line_real_offsets` for the geometry-based placement.
    """
    spans = []
    for grid in profiles.values():
        x = np.asarray(grid["x"], dtype=float)
        x = x[np.isfinite(x)]
        if x.size > 1:
            spans.append(float(x.max() - x.min()))
    return max(spans) if spans else 1000.0


def _dataset_uv(data: MapData) -> dict[str, tuple[float, float]]:
    """Adapt :class:`MapData` stations to :func:`geometry.survey_uv`.

    See :mod:`pycsamt.map.geometry` for the shared projection math
    (also used by ``pycsamt.app.web``'s 3-D map) — this just extracts
    the plain id/lat/lon/line arrays it needs from ``data.stations``.
    """
    ids: list[str] = []
    lats: list[float] = []
    lons: list[float] = []
    lines: list[str] = []
    for station in data.stations:
        if station.latitude is None or station.longitude is None:
            continue
        ids.append(station.id)
        lats.append(float(station.latitude))
        lons.append(float(station.longitude))
        lines.append(station.line or "line")
    return survey_uv(ids, lats, lons, lines)


def _station_uv(group, stations, uv):
    """Return ``(x, y_offset)`` for one line.

    ``x`` is the per-station along-strike distance (metres) and
    ``y_offset`` is the line's cross-strike position (metres), both
    derived from :func:`_dataset_uv` when every station in the line
    has real coordinates. Falls back to the path-length/index scheme
    in :func:`_station_x` (with ``y_offset=nan``, signalling callers
    to use the synthetic per-panel spacing) otherwise.
    """
    ids = [str(s) for s in stations]
    u = np.array(
        [uv[i][0] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    v = np.array(
        [uv[i][1] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    if u.size and np.isfinite(u).all():
        return u, float(np.nanmedian(v))
    return _station_x(group, stations), float("nan")


def _section_uv(stations, uv):
    """Same as :func:`_station_uv`, for a precomputed section with no
    pseudosection ``group`` to fall back on — synthetic index spacing
    is used instead when coordinates are unavailable."""
    ids = [str(s) for s in stations]
    u = np.array(
        [uv[i][0] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    v = np.array(
        [uv[i][1] if i in uv else np.nan for i in ids],
        dtype=float,
    )
    if u.size and np.isfinite(u).all():
        return u, float(np.nanmedian(v))
    return np.arange(len(ids), dtype=float) * 100.0, float("nan")


def _line_real_offsets(profiles) -> dict[str, float] | None:
    """Return ``{line: cross-strike offset (m)}`` when every line has
    real geometry, else ``None`` so callers fall back to synthetic,
    index-based spacing."""
    return normalize_offsets(
        {name: grid.get("y_offset") for name, grid in profiles.items()}
    )


def _line_offset(name, idx, real_offsets, unit, options) -> float:
    """Cross-strike placement for one line, in metres — see
    :func:`geometry.resolve_offset`."""
    return resolve_offset(
        name, idx, real_offsets, unit, float(options.line_spacing)
    )


def _filtered_values(grid, options):
    out = np.asarray(grid["value"], dtype=float).copy()
    if options.rho_range:
        lo, hi = options.rho_range
        rho = np.asarray(grid["rho"], dtype=float)
        out[(rho < lo) | (rho > hi)] = np.nan
    return out


def _color_values(values, options):
    arr = np.asarray(values, dtype=float)
    if _volume_quantity(options.quantity) == "rho" and options.log_color:
        return np.where(arr > 0, np.log10(arr), np.nan)
    return arr


def _crange(options) -> tuple[float | None, float | None]:
    """Return (cmin, cmax) in color space for the colorbar range."""
    if not options.value_range:
        return None, None
    lo, hi = options.value_range
    if _volume_quantity(options.quantity) == "rho" and options.log_color:
        lo = float(np.log10(lo)) if lo and lo > 0 else None
        hi = float(np.log10(hi)) if hi and hi > 0 else None
        return lo, hi
    return float(lo), float(hi)


def _elev_for(grid, options):
    """Return per-station surface elevation, or zeros when topo is off."""
    x = np.asarray(grid["x"], dtype=float)
    elev = np.asarray(grid.get("elev", []), dtype=float)
    if not options.topography or elev.size != x.size:
        return np.zeros_like(x)
    return np.where(np.isfinite(elev), elev, 0.0)


def _values_at_depth(grid, depth, options):
    z = np.asarray(grid["z"], dtype=float)
    values = _filtered_values(grid, options)
    order = np.argsort(z)
    z = z[order]
    values = values[order, :]
    out = np.empty(values.shape[1], dtype=float)
    for col in range(values.shape[1]):
        good = np.isfinite(z) & np.isfinite(values[:, col])
        if good.sum() == 0:
            out[col] = np.nan
            continue
        out[col] = np.interp(
            float(depth),
            z[good],
            values[good, col],
            left=np.nan,
            right=np.nan,
        )
    return out


def _marker_z_offset(profiles) -> float:
    """Small upward nudge (1% of the survey's own depth span) so
    station markers clear the panel surface instead of z-fighting it."""
    spans = []
    for grid in profiles.values():
        z = np.asarray(grid.get("z", []), dtype=float)
        z = z[np.isfinite(z)]
        if z.size:
            spans.append(float(z.max() - z.min()))
    span = max(spans) if spans else 0.0
    return max(span * 0.01, 1.0) if span > 0 else 1.0


def _add_station_markers(fig, profiles, options) -> None:
    """Overlay per-station markers at the survey surface across lines."""
    go = require_plotly()

    az = np.deg2rad(float(options.azimuth))
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    # A marker at exactly the same (x, y, z) as the panel's own top
    # edge z-fights with that surface in WebGL — same coordinates, so
    # the two triangles flicker/occlude depending on view angle. Lift
    # markers a hair above the surface, scaled to the survey's own
    # depth extent so it reads the same regardless of scale.
    z_offset = _marker_z_offset(profiles)
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    labels: list[str] = []
    hover: list[str] = []
    for line_idx, (name, grid) in enumerate(profiles.items()):
        x = np.asarray(grid["x"], dtype=float)
        elev = _elev_for(grid, options)
        offset = _line_offset(name, line_idx, real_offsets, unit, options)
        stations = grid.get("stations", np.arange(x.size))
        for j in range(x.size):
            xs.append(float(x[j] + offset * np.sin(az)))
            ys.append(float(offset * np.cos(az)))
            zs.append(float(elev[j]) + z_offset)
            labels.append(str(stations[j]))
            hover.append(f"{name} · {stations[j]}")
    if not xs:
        return
    show_labels = bool(options.station_labels)
    fig.add_trace(
        go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="markers+text" if show_labels else "markers",
            marker=dict(
                symbol=str(options.station_symbol),
                size=int(options.station_size),
                color=str(options.station_color),
                line=dict(width=0),
            ),
            text=labels if show_labels else None,
            textposition="top center",
            textfont=dict(size=9, color=str(options.station_color)),
            hovertext=hover,
            hovertemplate="%{hovertext}<extra></extra>",
            name="stations",
            showlegend=False,
        )
    )


def _terrain_trace(x, y, elev, name):
    """Return a Scatter3d line tracing the terrain top for one line."""
    go = require_plotly()

    return go.Scatter3d(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        z=np.asarray(elev, dtype=float),
        mode="lines",
        line=dict(color="#8d6e63", width=5),
        name=f"{name} terrain",
        showlegend=False,
        hoverinfo="skip",
    )


def _pad_row(arr, width):
    """Pad a 1-D array with NaN for ragged stacking."""
    arr = np.asarray(arr, dtype=float)
    if arr.size >= width:
        return arr[:width]
    out = np.full(width, np.nan, dtype=float)
    out[: arr.size] = arr
    return out


def _dense_volume_grid(profiles, options):
    """Build a regular ``(x, y, z)`` grid + resistivity volume for
    ``go.Volume``/``go.Isosurface``.

    Those trace types need genuine 3-D grid density to reconstruct
    smooth iso-surfaces. Each line's own ``(x, z)`` axes are too
    sparse and mutually inconsistent for that directly — this
    regrids every line onto one shared ``(x, z)`` reference, then
    densifies the cross-line (``y``) axis by interpolation, the same
    technique used by the ``pycsamt.app.web`` 3-D map's block mode.

    Returns ``None`` when there isn't enough data (or SciPy isn't
    available) to build a grid.
    """
    if not profiles:
        return None
    try:
        from scipy.interpolate import RegularGridInterpolator
    except Exception:  # noqa: BLE001 - scipy is an optional extra
        return None

    names = list(profiles.keys())
    n_lines = len(names)
    unit = _line_offset_unit(profiles)
    real_offsets = _line_real_offsets(profiles)
    y_vals = np.array(
        [
            _line_offset(name, idx, real_offsets, unit, options)
            for idx, name in enumerate(names)
        ],
        dtype=float,
    )

    ref_name = max(names, key=lambda n: np.asarray(profiles[n]["x"]).size)
    x_ref = np.unique(np.asarray(profiles[ref_name]["x"], dtype=float))
    z_ref = np.unique(np.asarray(profiles[ref_name]["z"], dtype=float))
    if x_ref.size < 2 or z_ref.size < 2:
        return None

    rho_vol = np.full((x_ref.size, n_lines, z_ref.size), np.nan, dtype=float)
    for idx, name in enumerate(names):
        grid = profiles[name]
        x = np.asarray(grid["x"], dtype=float)
        z = np.asarray(grid["z"], dtype=float)
        if x.size < 2 or z.size < 2:
            continue
        values = _color_values(
            np.asarray(grid["value"], dtype=float), options
        )
        xo, zo = np.argsort(x), np.argsort(z)
        x_sorted, z_sorted = x[xo], z[zo]
        if np.unique(x_sorted).size < 2 or np.unique(z_sorted).size < 2:
            continue
        vals_sorted = values[np.ix_(zo, xo)]
        try:
            interp = RegularGridInterpolator(
                (z_sorted, x_sorted), vals_sorted,
                bounds_error=False, fill_value=np.nan,
            )
        except ValueError:
            continue
        zz, xx = np.meshgrid(z_ref, x_ref, indexing="ij")
        sampled = interp(np.column_stack([zz.ravel(), xx.ravel()]))
        rho_vol[:, idx, :] = sampled.reshape(z_ref.size, x_ref.size).T

    if n_lines >= 2:
        order = np.argsort(y_vals)
        y_sorted = y_vals[order]
        n_y_dense = max(20, n_lines * 6)
        y_dense = np.linspace(y_sorted[0], y_sorted[-1], n_y_dense)
        dense = np.empty((x_ref.size, n_y_dense, z_ref.size), dtype=float)
        for ix in range(x_ref.size):
            for iz in range(z_ref.size):
                dense[ix, :, iz] = _interp_finite_1d(
                    y_dense, y_sorted, rho_vol[ix, order, iz]
                )
        rho_vol, y_vals = dense, y_dense

    return x_ref, y_vals, -np.abs(z_ref), rho_vol


def _interp_finite_1d(x_new, x, values):
    """1-D interpolation that skips NaNs instead of propagating them."""
    x = np.asarray(x, dtype=float)
    values = np.asarray(values, dtype=float)
    good = np.isfinite(x) & np.isfinite(values)
    if good.sum() >= 2:
        return np.interp(x_new, x[good], values[good])
    if good.sum() == 1:
        return np.full_like(x_new, values[good][0], dtype=float)
    return np.full_like(x_new, np.nan, dtype=float)


def _iso_range(vals, options):
    lo, hi = _iso_color_range(options)
    if lo is not None and hi is not None:
        return lo, hi
    arr = np.asarray(vals, dtype=float)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return 0.0, 1.0
    return (
        float(np.nanpercentile(finite, 20)),
        float(np.nanpercentile(finite, 80)),
    )


def _iso_color_range(options) -> tuple[float | None, float | None]:
    """``(isomin, isomax)`` in colour space, from the user-facing
    resistivity-range filter (``options.rho_range``) so it masks the
    block/iso-surface the same way it already masks fence/depth-slice
    cell values — without pre-filtering the dense grid to NaN, which
    would starve ``go.Volume``/``go.Isosurface`` of the density they
    need to interpolate a smooth shape.
    """
    if options.rho_range:
        lo, hi = options.rho_range
        if _volume_quantity(options.quantity) == "rho" and options.log_color:
            lo = float(np.log10(lo)) if lo and lo > 0 else None
            hi = float(np.log10(hi)) if hi and hi > 0 else None
        return lo, hi
    if options.iso_range:
        return options.iso_range
    return None, None


def _volume_quantity(quantity: str) -> str:
    if quantity.lower() in {"phase", "phi"}:
        return "phase"
    return "rho"


def _colorbar_title(options) -> str:
    if _volume_quantity(options.quantity) == "phase":
        return "φ (°)"
    return "log₁₀ ρ (Ω·m)" if options.log_color else "ρ (Ω·m)"


def _slice_depths(profiles, options):
    all_z = np.concatenate([
        np.asarray(grid["z"], dtype=float)
        for grid in profiles.values()
    ])
    all_z = all_z[np.isfinite(all_z)]
    if all_z.size == 0:
        return np.array([0.0])
    if options.depth_range:
        lo, hi = options.depth_range
    else:
        lo = float(np.nanmin(all_z))
        hi = float(np.nanmax(all_z))
    return np.linspace(lo, hi, max(1, int(options.n_slices)))


def _surface_contours(options):
    go = require_plotly()

    return go.surface.Contours(
        z=go.surface.contours.Z(
            show=bool(options.show_contours),
            usecolormap=True,
            project_z=False,
        )
    )


def _style_3d(fig, options, colors) -> None:
    title = options.title or f"pyCSAMT 3-D {options.mode} map"
    xu = getattr(options, "x_unit", "m")
    zu = getattr(options, "depth_unit", "m")
    z_label = "Elevation − depth" if options.topography else "Depth"
    z_title = f"{z_label} ({zu})"
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title=f"Profile distance ({xu})",
            yaxis_title=f"Line offset ({xu})",
            zaxis_title=z_title,
            bgcolor=colors["plot"],
            aspectmode=getattr(options, "aspectmode", "data") or "data",
        ),
        paper_bgcolor=colors["paper"],
        font=dict(color=colors["text"]),
        margin=dict(l=0, r=0, t=40, b=0),
    )


def _empty_3d_figure(colors):
    go = require_plotly()

    fig = go.Figure()
    fig.add_annotation(
        text="No 3-D map data available",
        x=0.5,
        y=0.5,
        xref="paper",
        yref="paper",
        showarrow=False,
        font=dict(color=colors["text"]),
    )
    fig.update_layout(
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["plot"],
    )
    return fig
