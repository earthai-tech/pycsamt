# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Overlay helpers for station, profile, and 3-D maps."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from ._backends import (
    require_plotly,
    require_pyproj_transformer,
    scipy_griddata,
)


@dataclass
class ContourOverlay:
    """Description of an interpolated contour overlay."""

    values: np.ndarray
    x: np.ndarray
    y: np.ndarray
    levels: int = 12
    cmap: str = "viridis"
    opacity: float = 0.65


@dataclass
class TopographyOverlay:
    """Topography aligned to map coordinates."""

    elevation: np.ndarray
    source: str = "stations"
    opacity: float = 0.70


@dataclass(frozen=True)
class CRSConfig:
    """Coordinate transform settings."""

    source: int | str
    target: int | str = 4326
    always_xy: bool = True


@dataclass(frozen=True)
class BasemapConfig:
    """Plotly geographic map layout settings."""

    style: str
    center: dict[str, float]
    zoom: int
    bearing: float = 0.0
    layers: tuple = ()


# ESRI ArcGIS raster tile services — free, no token required.
ESRI_TILES: dict[str, str] = {
    "esri-satellite": (
        "https://server.arcgisonline.com/ArcGIS/rest/services/"
        "World_Imagery/MapServer/tile/{z}/{y}/{x}"
    ),
    "esri-topo": (
        "https://server.arcgisonline.com/ArcGIS/rest/services/"
        "World_Topo_Map/MapServer/tile/{z}/{y}/{x}"
    ),
    "esri-natgeo": (
        "https://server.arcgisonline.com/ArcGIS/rest/services/"
        "NatGeo_World_Map/MapServer/tile/{z}/{y}/{x}"
    ),
    "esri-ocean": (
        "https://server.arcgisonline.com/ArcGIS/rest/services/"
        "Ocean/World_Ocean_Base/MapServer/tile/{z}/{y}/{x}"
    ),
    "esri-street": (
        "https://server.arcgisonline.com/ArcGIS/rest/services/"
        "World_Street_Map/MapServer/tile/{z}/{y}/{x}"
    ),
}

# (value, label) pairs for basemap pickers.
BASEMAP_STYLES: tuple[tuple[str, str], ...] = (
    ("open-street-map", "Open Street Map"),
    ("carto-positron", "Light (Carto)"),
    ("carto-darkmatter", "Dark (Carto)"),
    ("esri-satellite", "ESRI Satellite"),
    ("esri-topo", "ESRI Topographic"),
    ("esri-natgeo", "ESRI NatGeo"),
    ("esri-street", "ESRI Street"),
    ("esri-ocean", "ESRI Ocean"),
)


def basemap_style_and_layers(
    style: str | None,
    *,
    dark: bool = False,
) -> tuple[str, list]:
    """Resolve a style name into a (base style, raster layers) pair.

    ESRI styles are rendered as ``white-bg`` + a raster tile layer;
    native Plotly styles are returned as-is with no layers.
    """
    chosen = style or default_basemap(dark)
    if chosen in ESRI_TILES:
        return "white-bg", [
            {
                "below": "traces",
                "sourcetype": "raster",
                "source": [ESRI_TILES[chosen]],
                "opacity": 1.0,
            }
        ]
    return chosen, []


def resolve_crs_info(
    mode: str = "geo",
    *,
    zone: int = 50,
    hemisphere: str = "N",
    epsg: int | str = 4326,
) -> str:
    """Return a human-readable CRS description."""
    if mode == "geo":
        return "EPSG:4326 Geographic lat/lon (WGS 84)"
    if mode == "utm":
        if hemisphere.upper() == "N":
            code = 32600 + int(zone)
        else:
            code = 32700 + int(zone)
        hem = hemisphere.upper()
        return f"EPSG:{code} UTM Zone {zone}{hem} (WGS 84)"
    return f"EPSG:{epsg}"


def normalize_epsg(epsg: int | str) -> str:
    """Return an EPSG authority string."""
    value = str(epsg).upper().strip()
    if value.startswith("EPSG:"):
        return value
    return f"EPSG:{value}"


def transform_xy(
    x: np.ndarray,
    y: np.ndarray,
    *,
    crs: CRSConfig,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform coordinates between two CRS definitions."""
    Transformer = require_pyproj_transformer()
    transformer = Transformer.from_crs(
        normalize_epsg(crs.source),
        normalize_epsg(crs.target),
        always_xy=crs.always_xy,
    )
    xt, yt = transformer.transform(x, y)
    return (
        np.asarray(xt, dtype=float),
        np.asarray(yt, dtype=float),
    )


def reproject_xy_to_lonlat(
    x: np.ndarray,
    y: np.ndarray,
    *,
    epsg: int | str,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproject *x*, *y* from *epsg* to WGS84 lon/lat."""
    return transform_xy(
        x,
        y,
        crs=CRSConfig(source=epsg, target=4326),
    )


def build_contour_overlay(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    *,
    levels: int = 12,
    cmap: str = "Viridis",
    opacity: float = 0.65,
    grid_size: int = 80,
    mode: str = "lines",
) -> Any:
    """Build an interpolated map contour overlay.

    Returns a Plotly ``Contour`` trace on a regular grid.
    """
    go = require_plotly()
    xi, yi, zz = interpolate_overlay_grid(
        x,
        y,
        values,
        grid_size=grid_size,
    )
    return go.Contour(
        x=xi,
        y=yi,
        z=zz,
        ncontours=max(2, int(levels)),
        colorscale=cmap,
        opacity=float(opacity),
        contours=dict(
            coloring=_contour_coloring(mode),
            showlines=mode in {"lines", "both"},
        ),
        showscale=True,
    )


def build_geo_contour_image(
    lons: np.ndarray,
    lats: np.ndarray,
    values: np.ndarray,
    *,
    cmap: str = "jet",
    n_levels: int = 12,
    opacity: float = 0.6,
    mode: str = "filled+lines",
    log_scale: bool = False,
    grid_res: int = 160,
    expand: float = 0.06,
    interp: str = "cubic",
    smooth_sigma: float = 0.0,
) -> dict | None:
    """Render a Surfer-style filled-contour PNG for a basemap image layer.

    Interpolates sparse station *values* onto a regular lon/lat grid and
    rasterises filled (and/or line) contours to a transparent PNG. The
    result is meant to be added to a Plotly map as an ``image`` layer.

    Returns
    -------
    dict | None
        ``{"image": "data:image/png;base64,...", "coordinates": [...],
        "vmin": float, "vmax": float}`` (coordinates are the TL, TR, BR,
        BL ``[lon, lat]`` corners), or ``None`` when it cannot be built.
    """
    import base64
    import io

    from matplotlib.backends.backend_agg import (
        FigureCanvasAgg,
    )
    from matplotlib.figure import Figure

    lon = np.asarray(lons, dtype=float).ravel()
    lat = np.asarray(lats, dtype=float).ravel()
    val = np.asarray(values, dtype=float).ravel()
    good = np.isfinite(lon) & np.isfinite(lat) & np.isfinite(val)
    lon, lat, val = lon[good], lat[good], val[good]
    if val.size < 3:
        return None
    if log_scale:
        pos = val > 0
        if pos.sum() < 3:
            return None
        lon, lat, val = lon[pos], lat[pos], np.log10(val[pos])

    lon_min, lon_max = float(lon.min()), float(lon.max())
    lat_min, lat_max = float(lat.min()), float(lat.max())
    dlon = (lon_max - lon_min) * expand or 1e-4
    dlat = (lat_max - lat_min) * expand or 1e-4
    lon_min, lon_max = lon_min - dlon, lon_max + dlon
    lat_min, lat_max = lat_min - dlat, lat_max + dlat

    xi, yi, zz = interpolate_overlay_grid(
        lon,
        lat,
        val,
        grid_size=int(grid_res),
        method=interp,
    )
    if not np.isfinite(zz).any():
        return None
    if smooth_sigma and smooth_sigma > 0:
        try:
            from scipy.ndimage import gaussian_filter

            mask = np.isfinite(zz)
            filled = np.where(mask, zz, np.nanmean(zz[mask]))
            zz = gaussian_filter(filled, sigma=float(smooth_sigma))
            zz[~mask] = np.nan
        except Exception:
            pass

    vmin = float(np.nanmin(zz))
    vmax = float(np.nanmax(zz))
    fig = Figure(figsize=(6, 6), dpi=100)
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(lon_min, lon_max)
    ax.set_ylim(lat_min, lat_max)
    levels = max(2, int(n_levels))
    if mode in ("filled", "filled+lines"):
        ax.contourf(
            xi,
            yi,
            zz,
            levels=levels,
            cmap=cmap,
            alpha=float(opacity),
        )
    if mode in ("lines", "filled+lines"):
        ax.contour(
            xi,
            yi,
            zz,
            levels=levels,
            colors="k",
            linewidths=0.5,
            alpha=0.5,
        )
    buf = io.BytesIO()
    canvas.print_png(buf)
    fig.clf()
    b64 = base64.b64encode(buf.getvalue()).decode()
    return {
        "image": f"data:image/png;base64,{b64}",
        "coordinates": [
            [lon_min, lat_max],  # top-left
            [lon_max, lat_max],  # top-right
            [lon_max, lat_min],  # bottom-right
            [lon_min, lat_min],  # bottom-left
        ],
        "vmin": vmin,
        "vmax": vmax,
    }


def interpolate_overlay_grid(
    x: np.ndarray,
    y: np.ndarray,
    values: np.ndarray,
    *,
    grid_size: int = 80,
    method: str = "linear",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Interpolate scattered values to a regular grid."""
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    val_arr = np.asarray(values, dtype=float)
    good = np.isfinite(x_arr) & np.isfinite(y_arr) & np.isfinite(val_arr)
    if good.sum() < 3:
        raise ValueError("At least three finite points are required.")
    xi = np.linspace(
        float(np.nanmin(x_arr[good])),
        float(np.nanmax(x_arr[good])),
        int(grid_size),
    )
    yi = np.linspace(
        float(np.nanmin(y_arr[good])),
        float(np.nanmax(y_arr[good])),
        int(grid_size),
    )
    xx, yy = np.meshgrid(xi, yi)
    griddata = scipy_griddata()
    if griddata is not None:
        try:
            zz = griddata(
                np.column_stack([x_arr[good], y_arr[good]]),
                val_arr[good],
                (xx, yy),
                method=method,
            )
            if np.isfinite(zz).any():
                return xi, yi, zz
        except Exception:
            pass
    zz = _nearest_grid(
        x_arr[good],
        y_arr[good],
        val_arr[good],
        xx,
        yy,
    )
    return xi, yi, zz


def build_topography_overlay(
    x: np.ndarray,
    y: np.ndarray,
    elevation: np.ndarray,
    *,
    opacity: float = 0.70,
    colorscale: str = "Earth",
) -> Any:
    """Build a Plotly 3-D terrain mesh or surface."""
    go = require_plotly()

    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    z_arr = np.asarray(elevation, dtype=float)
    if z_arr.ndim == 2:
        return go.Surface(
            x=x_arr,
            y=y_arr,
            z=z_arr,
            surfacecolor=z_arr,
            colorscale=colorscale,
            opacity=float(opacity),
            name="Topography",
            showscale=False,
        )
    return go.Mesh3d(
        x=x_arr,
        y=y_arr,
        z=z_arr,
        intensity=z_arr,
        colorscale=colorscale,
        opacity=float(opacity),
        name="Topography",
        showscale=False,
    )


def build_station_label_overlay(
    x: np.ndarray,
    y: np.ndarray,
    labels: list[str] | np.ndarray,
    *,
    geo: bool = False,
    name: str = "Station labels",
    color: str = "#111827",
) -> Any:
    """Build a reusable station-label overlay trace."""
    go = require_plotly()
    kwargs = dict(
        mode="text",
        text=[str(label) for label in labels],
        name=name,
        textfont=dict(color=color),
        showlegend=False,
        hoverinfo="skip",
    )
    if geo:
        return _scatter_map_trace(
            go,
            lon=np.asarray(x, dtype=float),
            lat=np.asarray(y, dtype=float),
            **kwargs,
        )
    return go.Scatter(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        **kwargs,
    )


def build_profile_line_overlay(
    x: np.ndarray,
    y: np.ndarray,
    *,
    geo: bool = False,
    name: str = "Profile",
    color: str = "#2563eb",
    width: float = 2.0,
) -> Any:
    """Build a reusable profile-line overlay trace."""
    go = require_plotly()
    kwargs = dict(
        mode="lines",
        name=name,
        line=dict(color=color, width=float(width)),
        hoverinfo="skip",
        showlegend=False,
    )
    if geo:
        return _scatter_map_trace(
            go,
            lon=np.asarray(x, dtype=float),
            lat=np.asarray(y, dtype=float),
            **kwargs,
        )
    return go.Scatter(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        **kwargs,
    )


def build_basemap_layout(
    lon: np.ndarray,
    lat: np.ndarray,
    *,
    dark: bool = False,
    style: str | None = None,
    bearing: float = 0.0,
) -> BasemapConfig:
    """Return layout settings for geographic station maps."""
    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    good = np.isfinite(lon_arr) & np.isfinite(lat_arr)
    if not good.any():
        center = {"lat": 0.0, "lon": 0.0}
        zoom = 1
    else:
        center = {
            "lat": float(np.nanmean(lat_arr[good])),
            "lon": float(np.nanmean(lon_arr[good])),
        }
        zoom = auto_map_zoom(lon_arr[good], lat_arr[good])
    base_style, layers = basemap_style_and_layers(style, dark=dark)
    return BasemapConfig(
        style=base_style,
        center=center,
        zoom=zoom,
        bearing=float(bearing),
        layers=tuple(layers),
    )


def default_basemap(dark: bool = False) -> str:
    """Return a default public tile style."""
    return "carto-darkmatter" if dark else "open-street-map"


def auto_map_zoom(
    lon: np.ndarray,
    lat: np.ndarray,
) -> int:
    """Estimate a reasonable Plotly map zoom."""
    lon_arr = np.asarray(lon, dtype=float)
    lat_arr = np.asarray(lat, dtype=float)
    if lon_arr.size == 0 or lat_arr.size == 0:
        return 1
    lat_span = float(np.nanmax(lat_arr) - np.nanmin(lat_arr))
    lon_span = float(np.nanmax(lon_arr) - np.nanmin(lon_arr))
    span = max(lat_span, lon_span, 1e-6)
    return int(max(2, min(14, 8 - np.log2(span + 0.001))))


def _scatter_map_trace(go, **kwargs):
    if hasattr(go, "Scattermap"):
        return go.Scattermap(**kwargs)
    return go.Scattermapbox(**kwargs)


def _contour_coloring(mode: str) -> str:
    if mode in {"fill", "filled", "heatmap"}:
        return "fill"
    if mode == "both":
        return "heatmap"
    return "none"


def _nearest_grid(x, y, values, xx, yy):
    zz = np.empty_like(xx, dtype=float)
    points = np.column_stack([x, y])
    grid = np.column_stack([xx.ravel(), yy.ravel()])
    for i, point in enumerate(grid):
        dist = np.sum((points - point) ** 2, axis=1)
        zz.ravel()[i] = values[int(np.argmin(dist))]
    return zz
