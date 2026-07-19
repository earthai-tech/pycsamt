# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""2-D station-map API for pyCSAMT surveys."""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import numpy as np

from ._backends import (
    require_matplotlib_pyplot,
    require_plotly,
)
from ._core import (
    MapData,
    ensure_map_data,
    skin_depth_at_frequency,
    station_dataframe,
    value_at_frequency,
)
from .config import StationMapOptions
from .overlays import (
    build_basemap_layout,
    build_profile_line_overlay,
)
from .styles import (
    is_dark_theme,
    theme_colors,
    to_plotly_cmap,
)


class StationMap:
    """Builder object for station-map figures."""

    def __init__(
        self,
        data: Any,
        *,
        options: StationMapOptions | None = None,
        **ensure_kwargs: Any,
    ) -> None:
        self.data: MapData = ensure_map_data(
            data,
            **ensure_kwargs,
        )
        self.options = options or StationMapOptions()

    def with_overlay(
        self,
        overlay: str,
        **kwargs: Any,
    ) -> StationMap:
        """Return a copy with a different overlay."""
        opts = replace(
            self.options,
            overlay=overlay,
            **kwargs,
        )
        new = object.__new__(StationMap)
        new.data = self.data
        new.options = opts
        return new

    def with_options(self, **kwargs: Any) -> StationMap:
        """Return a copy with updated options."""
        new = object.__new__(StationMap)
        new.data = self.data
        new.options = replace(self.options, **kwargs)
        return new

    def figure(self) -> Any:
        """Build and return the station-map figure."""
        return build_station_map(self.data, self.options)


def plot_station_map(
    data: Any,
    *,
    options: StationMapOptions | None = None,
    **kwargs: Any,
) -> Any:
    """Build a 2-D station map."""
    return StationMap(
        data,
        options=options,
        **kwargs,
    ).figure()


def build_station_map(
    data: MapData,
    options: StationMapOptions | None = None,
) -> Any:
    """Build a Plotly 2-D station map."""
    opts = options or StationMapOptions()
    df = _station_values(data, opts)
    colors = theme_colors(opts.theme)
    if opts.backend == "matplotlib":
        return _matplotlib_station_map(df, opts, colors)
    if opts.backend != "plotly":
        msg = f"Unknown station map backend: {opts.backend}"
        raise ValueError(msg)
    dark = is_dark_theme(opts.theme)
    if df.empty:
        return _empty_station_figure(colors)

    if _has_geo(df):
        return _geo_station_map(df, opts, colors, dark)
    return _profile_station_map(df, opts, colors)


def _station_values(
    data: MapData,
    opts: StationMapOptions,
):
    df = station_dataframe(data)
    if df.empty:
        return df
    overlay = opts.overlay.lower()
    if overlay in {"index", "station"}:
        df["_value"] = df["Index"]
        df["_label"] = "Station index"
    elif overlay == "elevation":
        df["_value"] = df["Elevation"]
        df["_label"] = "Elevation (m)"
    elif overlay in {"rho", "resistivity"}:
        freq = getattr(opts, "frequency", None) or 1.0
        vals = value_at_frequency(
            data,
            frequency=float(freq),
            quantity="rho",
            component=opts.component,
            tolerance=opts.frequency_tolerance,
        )
        df["_value"] = df["ID"].map(vals)
        df["_label"] = "ρ (Ω·m)"
    elif overlay == "phase":
        freq = getattr(opts, "frequency", None) or 1.0
        vals = value_at_frequency(
            data,
            frequency=float(freq),
            quantity="phase",
            component=opts.component,
            tolerance=opts.frequency_tolerance,
        )
        df["_value"] = df["ID"].map(vals)
        df["_label"] = "φ (°)"
    elif overlay in {"depth", "skin_depth"}:
        freq = getattr(opts, "frequency", None) or 1.0
        vals = skin_depth_at_frequency(
            data,
            frequency=float(freq),
            component=opts.component,
            tolerance=opts.frequency_tolerance,
        )
        df["_value"] = df["ID"].map(vals)
        df["_label"] = "δ skin depth (m)"
    else:
        if opts.overlay in df:
            df["_value"] = df[opts.overlay]
        else:
            df["_value"] = df["Index"]
        df["_label"] = opts.overlay
    df["_plot_value"] = _plot_values(df["_value"], opts)
    return df


def _has_geo(df) -> bool:
    lat = np.asarray(df["Latitude"], dtype=float)
    lon = np.asarray(df["Longitude"], dtype=float)
    return bool(np.isfinite(lat).any() and np.isfinite(lon).any())


def _geo_station_map(df, opts, colors, dark: bool):
    go = require_plotly()

    cmin, cmax = _shared_color_range(df, opts)
    scale_shown = False
    fig = go.Figure()
    for line, group in df.groupby("Line", dropna=False):
        line_name = str(line or "line")
        if opts.line_filter and line_name not in opts.line_filter:
            continue
        lat = group["Latitude"].astype(float)
        lon = group["Longitude"].astype(float)
        marker_size = _marker_sizes(group["ID"], opts)
        if opts.show_contours:
            _add_density_layer(fig, group, opts)
        fig.add_trace(
            _scatter_map_trace(
                go,
                lat=lat,
                lon=lon,
                mode=_marker_mode(opts.show_labels),
                text=group["ID"],
                customdata=group["ID"],
                name=line_name,
                marker=dict(
                    size=marker_size,
                    color=group["_plot_value"],
                    colorscale=to_plotly_cmap(opts.cmap),
                    opacity=opts.opacity,
                    showscale=not scale_shown,
                    cmin=cmin,
                    cmax=cmax,
                    colorbar=dict(
                        title=dict(
                            text=_color_title(df, opts),
                            side="right",
                        ),
                    ),
                ),
                textposition="top right",
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "lat=%{lat:.5f}<br>lon=%{lon:.5f}"
                    "<extra></extra>"
                ),
            )
        )
        scale_shown = True
        if opts.show_profiles and len(group) > 1:
            fig.add_trace(
                build_profile_line_overlay(
                    lon,
                    lat,
                    geo=True,
                    name=f"{line_name} profile",
                    color=colors["accent"],
                    width=2,
                )
            )

    basemap = build_basemap_layout(
        df["Longitude"],
        df["Latitude"],
        dark=dark,
        style=opts.basemap,
        bearing=opts.bearing,
    )
    layers = list(basemap.layers)
    contour_layer = _contour_image_layer(df, opts)
    if contour_layer is not None:
        layers.append(contour_layer)
    map_layout = dict(
        style=basemap.style,
        center=basemap.center,
        zoom=basemap.zoom,
        bearing=basemap.bearing,
    )
    if layers:
        map_layout["layers"] = layers
    fig.update_layout(
        **{_map_layout_key(fig): map_layout},
        paper_bgcolor=colors["paper"],
        font=dict(color=colors["text"]),
        margin=dict(l=0, r=0, t=30, b=0),
        legend=dict(orientation="h"),
        title=opts.title or None,
    )
    return fig


def _contour_image_layer(df, opts):
    """Build a Surfer-style filled-contour image map layer, or None."""
    if not opts.contour_image:
        return None
    from .overlays import build_geo_contour_image

    overlay = build_geo_contour_image(
        df["Longitude"],
        df["Latitude"],
        df["_value"],
        cmap=opts.cmap,
        n_levels=opts.contour_levels,
        opacity=float(opts.contour_opacity),
        mode=opts.contour_mode,
        log_scale=bool(opts.log_color),
        grid_res=int(opts.contour_grid_res),
        interp=opts.contour_interp,
        smooth_sigma=float(opts.contour_smooth),
    )
    if overlay is None:
        return None
    return {
        "below": "traces",
        "sourcetype": "image",
        "source": overlay["image"],
        "coordinates": overlay["coordinates"],
    }


def _scatter_map_trace(go, **kwargs):
    if hasattr(go, "Scattermap"):
        return go.Scattermap(**kwargs)
    return go.Scattermapbox(**kwargs)


def _density_map_trace(go, **kwargs):
    if hasattr(go, "Densitymap"):
        return go.Densitymap(**kwargs)
    return go.Densitymapbox(**kwargs)


def _map_layout_key(fig) -> str:
    modern = {"scattermap", "densitymap"}
    if any(getattr(trace, "type", "") in modern for trace in fig.data):
        return "map"
    return "mapbox"


def _matplotlib_station_map(df, opts, colors):
    plt = require_matplotlib_pyplot()

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor(colors["paper"])
    ax.set_facecolor(colors["plot"])
    ax.tick_params(colors=colors["text"])
    if df.empty:
        ax.text(
            0.5,
            0.5,
            "No stations available",
            ha="center",
            va="center",
            color=colors["text"],
        )
        return fig
    x_col = "Longitude" if _has_geo(df) else "Index"
    y_col = "Latitude" if _has_geo(df) else "_value"
    sizes = (
        np.asarray(
            _marker_sizes(df["ID"], opts),
            dtype=float,
        )
        ** 2
    )
    sc = ax.scatter(
        df[x_col],
        df[y_col],
        c=df["_plot_value"],
        s=sizes,
        cmap=opts.cmap,
        alpha=float(opts.opacity),
        vmin=_cmin(opts),
        vmax=_cmax(opts),
    )
    if opts.show_labels:
        for _, row in df.iterrows():
            ax.text(
                row[x_col],
                row[y_col],
                str(row["ID"]),
                color=colors["text"],
                fontsize=8,
            )
    if opts.show_profiles and _has_geo(df):
        for _, group in df.groupby("Line", dropna=False):
            ax.plot(
                group["Longitude"],
                group["Latitude"],
                color=colors["accent"],
                linewidth=1.2,
            )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label(
        _color_title(df, opts),
        color=colors["text"],
    )
    ax.set_xlabel(x_col, color=colors["text"])
    ax.set_ylabel(y_col, color=colors["text"])
    if opts.title:
        ax.set_title(opts.title, color=colors["text"])
    return fig


def _profile_station_map(df, opts, colors):
    go = require_plotly()

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=df["Index"],
            y=df["_value"],
            mode=_marker_mode(opts.show_labels),
            text=df["ID"],
            marker=dict(
                size=_marker_sizes(df["ID"], opts),
                color=df["_plot_value"],
                colorscale=to_plotly_cmap(opts.cmap),
                showscale=True,
                cmin=_cmin(opts),
                cmax=_cmax(opts),
                colorbar=dict(
                    title=dict(text=_color_title(df, opts), side="right")
                ),
            ),
        )
    )
    fig.update_layout(
        xaxis_title="Station index",
        yaxis_title=str(df["_label"].iloc[0]),
        paper_bgcolor=colors["paper"],
        plot_bgcolor=colors["plot"],
        font=dict(color=colors["text"]),
        title=opts.title or None,
    )
    return fig


def _plot_values(values, opts):
    arr = np.asarray(values, dtype=float)
    if opts.log_color:
        arr = np.where(arr > 0, np.log10(arr), np.nan)
    return arr


def _color_title(df, opts) -> str:
    label = str(df["_label"].iloc[0])
    return f"log₁₀ {label}" if opts.log_color else label


def _cmin(opts) -> float | None:
    if not opts.value_range:
        return None
    lo, _ = opts.value_range
    if opts.log_color and lo > 0:
        return float(np.log10(lo))
    return float(lo)


def _cmax(opts) -> float | None:
    if not opts.value_range:
        return None
    _, hi = opts.value_range
    if opts.log_color and hi > 0:
        return float(np.log10(hi))
    return float(hi)


def _shared_color_range(df, opts) -> tuple[float | None, float | None]:
    """Resolve one cmin/cmax pair shared by every line trace.

    Falls back to the finite range of ``df["_plot_value"]`` so that
    colors stay comparable across lines when ``opts.value_range`` is
    not set, instead of each line auto-scaling to its own values.
    """
    cmin, cmax = _cmin(opts), _cmax(opts)
    if cmin is not None and cmax is not None:
        return cmin, cmax
    values = np.asarray(df["_plot_value"], dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return cmin, cmax
    if cmin is None:
        cmin = float(finite.min())
    if cmax is None:
        cmax = float(finite.max())
    return cmin, cmax


def _add_density_layer(fig, group, opts) -> None:
    go = require_plotly()

    values = np.asarray(group["_plot_value"], dtype=float)
    lat = np.asarray(group["Latitude"], dtype=float)
    lon = np.asarray(group["Longitude"], dtype=float)
    good = np.isfinite(values) & np.isfinite(lat) & np.isfinite(lon)
    if good.sum() < 3:
        return
    fig.add_trace(
        _density_map_trace(
            go,
            lat=lat[good],
            lon=lon[good],
            z=values[good],
            radius=max(10, int(opts.marker_size * 2)),
            colorscale=to_plotly_cmap(opts.cmap),
            opacity=float(opts.contour_opacity),
            showscale=False,
            hoverinfo="skip",
            name="overlay density",
        )
    )


def _marker_sizes(ids, opts) -> list[int]:
    selected = str(opts.selected_id) if opts.selected_id else None
    return [
        int(opts.marker_size * 1.8)
        if selected == str(sid)
        else int(opts.marker_size)
        for sid in ids
    ]


def _marker_mode(show_labels: bool) -> str:
    return "markers+text" if show_labels else "markers"


def _empty_station_figure(colors):
    go = require_plotly()

    fig = go.Figure()
    fig.add_annotation(
        text="No stations available",
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
