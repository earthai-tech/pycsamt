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
    return fig


def _profile_grids(
    data: MapData,
    options: VolumeMapOptions,
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
        out[str(line)] = {
            "x": _station_x(group, stations),
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


def _fence_figure(profiles, options, colors):
    go = require_plotly()

    fig = go.Figure()
    az = np.deg2rad(float(options.azimuth))
    cmin, cmax = _crange(options)
    for idx, (name, grid) in enumerate(profiles.items()):
        values = _filtered_values(grid, options)
        z = -np.asarray(grid["z"], dtype=float)
        x = np.asarray(grid["x"], dtype=float)
        xx, zz = np.meshgrid(x, z)
        offset = idx * float(options.line_spacing)
        yy = np.full_like(
            xx,
            offset * np.cos(az),
        )
        xx = xx + offset * np.sin(az)
        elev = _elev_for(grid, options)
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
                colorbar=dict(title=_colorbar_title(options)),
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


def _block_figure(profiles, options, colors):
    go = require_plotly()

    xs, ys, zs, vals = _volume_points(profiles, options)
    fig = go.Figure()
    if vals:
        cmin, cmax = _crange(options)
        fig.add_trace(
            go.Volume(
                x=xs,
                y=ys,
                z=zs,
                value=vals,
                cmin=cmin,
                cmax=cmax,
                opacity=float(options.opacity),
                surface_count=int(options.surface_count),
                colorscale=to_plotly_cmap(options.cmap),
                colorbar=dict(title=_colorbar_title(options)),
            )
        )
    _style_3d(fig, options, colors)
    return fig


def _depth_figure(profiles, options, colors):
    go = require_plotly()

    depths = _slice_depths(profiles, options)
    az = np.deg2rad(float(options.azimuth))
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
        for line_idx, grid in enumerate(profiles.values()):
            values = _values_at_depth(
                grid,
                depth,
                options,
            )
            x = np.asarray(grid["x"], dtype=float)
            offset = line_idx * float(options.line_spacing)
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
                colorbar=dict(title=_colorbar_title(options)),
                contours=_surface_contours(options),
            )
        )
    if options.show_terrain and options.topography:
        for line_idx, (name, grid) in enumerate(profiles.items()):
            x = np.asarray(grid["x"], dtype=float)
            offset = line_idx * float(options.line_spacing)
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

    xs, ys, zs, vals = _volume_points(profiles, options)
    fig = go.Figure()
    if vals:
        lo, hi = _iso_range(vals, options)
        fig.add_trace(
            go.Isosurface(
                x=xs,
                y=ys,
                z=zs,
                value=vals,
                isomin=lo,
                isomax=hi,
                surface_count=int(options.surface_count),
                opacity=float(options.opacity),
                colorscale=to_plotly_cmap(options.cmap),
                caps=dict(x_show=False, y_show=False),
                colorbar=dict(title=_colorbar_title(options)),
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
        return np.arange(len(stations), dtype=float)
    return x


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


def _add_station_markers(fig, profiles, options) -> None:
    """Overlay per-station markers at the survey surface across lines."""
    go = require_plotly()

    az = np.deg2rad(float(options.azimuth))
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    text: list[str] = []
    for line_idx, (name, grid) in enumerate(profiles.items()):
        x = np.asarray(grid["x"], dtype=float)
        elev = _elev_for(grid, options)
        offset = line_idx * float(options.line_spacing)
        stations = grid.get("stations", np.arange(x.size))
        for j in range(x.size):
            xs.append(float(x[j] + offset * np.sin(az)))
            ys.append(float(offset * np.cos(az)))
            zs.append(float(elev[j]))
            text.append(f"{name} · {stations[j]}")
    if not xs:
        return
    fig.add_trace(
        go.Scatter3d(
            x=xs,
            y=ys,
            z=zs,
            mode="markers",
            marker=dict(
                symbol=str(options.station_symbol),
                size=int(options.station_size),
                color=str(options.station_color),
                line=dict(width=0),
            ),
            text=text,
            hovertemplate="%{text}<extra></extra>",
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


def _volume_points(profiles, options):
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    vals: list[float] = []
    az = np.deg2rad(float(options.azimuth))
    for line_idx, grid in enumerate(profiles.values()):
        values = _filtered_values(grid, options)
        color_values = _color_values(values, options)
        elev = _elev_for(grid, options)
        offset = line_idx * float(options.line_spacing)
        x_shift = offset * np.sin(az)
        y_value = offset * np.cos(az)
        for iz, depth in enumerate(grid["z"]):
            for ix, xval in enumerate(grid["x"]):
                val = color_values[iz, ix]
                if not np.isfinite(val):
                    continue
                xs.append(float(xval + x_shift))
                ys.append(float(y_value))
                zs.append(float(elev[ix] - depth))
                vals.append(float(val))
    return xs, ys, zs, vals


def _iso_range(vals, options):
    if options.iso_range:
        return options.iso_range
    arr = np.asarray(vals, dtype=float)
    return (
        float(np.nanpercentile(arr, 20)),
        float(np.nanpercentile(arr, 80)),
    )


def _volume_quantity(quantity: str) -> str:
    if quantity.lower() in {"phase", "phi"}:
        return "phase"
    return "rho"


def _colorbar_title(options) -> str:
    if _volume_quantity(options.quantity) == "phase":
        return "Phase (deg)"
    return "log10 rho" if options.log_color else "rho (Ohm.m)"


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
    z_title = (
        "Elevation − depth (m)"
        if options.topography
        else "Depth (m)"
    )
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="Profile distance",
            yaxis_title="Line offset",
            zaxis_title=z_title,
            bgcolor=colors["plot"],
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
