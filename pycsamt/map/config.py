# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Configuration dataclasses for :mod:`pycsamt.map`."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

MapTheme = Literal["light", "dark", "publication"]
MapBackend = Literal["plotly", "matplotlib"]
VolumeMode = Literal["fence", "block", "depth", "surface"]


@dataclass
class StationMapOptions:
    """Options for 2-D station maps."""

    overlay: str = "index"
    frequency: float | None = None
    frequency_tolerance: float | None = None
    # For inversion surveys: the depth (m below surface) at which
    # ``overlay="resistivity"`` slices the precomputed sections. ``None``
    # keeps the classic behaviour (apparent resistivity at ``frequency``).
    depth: float | None = None
    # Draw the per-station scatter markers. Turning them off leaves just
    # the filled-contour depth slice / value image over the basemap.
    show_markers: bool = True
    component: str = "xy"
    theme: MapTheme = "light"
    backend: MapBackend = "plotly"
    selected_id: str | None = None
    line_filter: tuple[str, ...] | None = None
    marker_size: int = 10
    basemap: str | None = None
    cmap: str = "plasma"
    value_range: tuple[float, float] | None = None
    log_color: bool = False
    opacity: float = 0.92
    show_labels: bool = True
    label_fontsize: float = 8.0
    label_rotation: float = 0.0
    show_profiles: bool = True
    elevation_mode: Literal["markers", "contours"] = "markers"
    show_contours: bool = False
    contour_image: bool = False
    contour_mode: str = "filled+lines"
    contour_levels: int = 12
    contour_opacity: float = 0.55
    contour_interp: str = "linear"
    contour_smooth: float = 1.0
    contour_grid_res: int = 150
    bearing: float = 0.0
    title: str = ""


@dataclass
class ProfileMapOptions:
    """Options for pseudosection and profile-view maps."""

    quantity: str = "rho"
    component: str = "xy"
    components: tuple[str, ...] = ("xy", "yx")
    theme: MapTheme = "light"
    backend: MapBackend = "plotly"
    period_range: tuple[float, float] | None = None
    phase_range: tuple[float, float] | None = None
    value_range: tuple[float, float] | None = None
    x_axis: Literal["station", "distance"] = "station"
    log_rho: bool = True
    height_per_panel: int = 260
    show_errbar: bool = False
    cmap: str | None = None
    by_line: bool = False
    line_cols: int | None = None


@dataclass
class VolumeMapOptions:
    """Options for 3-D fence, block, and depth-slice maps."""

    mode: VolumeMode = "fence"
    quantity: str = "resistivity"
    component: str = "xy"
    theme: MapTheme = "light"
    cmap: str = "RdYlBu_r"
    depth_range: tuple[float, float] | None = None
    period_range: tuple[float, float] | None = None
    rho_range: tuple[float, float] | None = None
    iso_range: tuple[float, float] | None = None
    value_range: tuple[float, float] | None = None
    # Percentile clip (low, high) for the *auto* colour range — ignored
    # when ``value_range`` pins it explicitly. Keeps a handful of
    # extreme cells (e.g. a ModEM volume's above-topography air fill at
    # 1e10-1e13 ohm.m) from stretching the whole colourscale so the real
    # earth reads as one flat band. ``None`` or ``(0, 100)`` restores the
    # raw min/max behaviour.
    crange_percentile: tuple[float, float] | None = (2.0, 98.0)
    # Hard "hide every cell above this resistivity (ohm.m)" cutoff. Unlike
    # ``rho_range`` (a visibility band that deliberately leaves the colour
    # scale put), this also drops the excluded cells from the colour-range
    # computation — the one knob that reliably removes air / overburden
    # fill from a fence/depth/block view and its colourbar at once.
    rho_display_max: float | None = None
    log_color: bool = True
    opacity: float = 0.85
    show_contours: bool = False
    show_labels: bool = True
    n_slices: int = 5
    surface_count: int = 12
    line_spacing: float = 1.0
    azimuth: float = 0.0
    topography: bool = True
    show_terrain: bool = True
    terrain_opacity: float = 0.7
    show_stations: bool = False
    # A native Plotly Scatter3d marker symbol, or "triangle-down"/
    # "triangle-down-open" -- pycsamt's own extension, rendered as real
    # geometry since Plotly's 3-D marker enum has no triangle at all
    # (see pycsamt.map.volume.TRIANGLE_DOWN_SYMBOLS).
    station_symbol: str = "diamond"
    station_size: int = 4
    station_color: str = "#1f2937"
    station_labels: bool = False
    # Screen-space rotation (degrees) for the per-station labels. Only
    # meaningful with ``station_labels`` -- the labels are drawn as
    # scene annotations (``go.Scatter3d`` text cannot rotate) so a
    # crowded line can tilt its labels to 45/90 the way a 2-D section
    # does. 0 = horizontal.
    station_label_angle: float = 0.0
    # Which markers get a *label* (markers themselves are unaffected --
    # you still see every station's position). ``station_label_names``
    # wins when set: only those station ids are labelled. Otherwise
    # ``station_label_fraction`` keeps an evenly-spaced fraction per
    # line (1.0 = all, 0.5 = every other, ...), first and last always
    # kept. Lets a crowded line stay readable without hiding stations.
    station_label_fraction: float = 1.0
    station_label_names: tuple[str, ...] | None = None
    # Cap on how many station markers/labels are drawn per line, evenly
    # spaced along that line's own station order (first and last are
    # always kept) -- lets a crowded fence/block with many stations
    # thin the overlay to, say, 5 or 10 per line instead of every one.
    # None (default) shows every station, unchanged from before this
    # option existed.
    max_stations: int | None = None
    aspectmode: str = "data"
    x_unit: str = "m"
    depth_unit: str = "m"
    # Fence panels: resample each per-line section onto a denser
    # spline grid so it reads as a smooth curtain, not raw stripes.
    smooth_sections: bool = True
    section_res: int = 100
    # Block / iso-surface / anomaly structure smoothing. 0.0 (default)
    # renders the volume exactly as reconstructed on the raw
    # station/depth lattice -- the current behaviour. A positive value
    # (the mapview "Light / Medium / Strong / Very strong" presets pass
    # 0.8 / 1.5 / 2.5 / 4.0) first resamples the volume onto a finer
    # 3-D lattice -- 2x to ~3.6x denser per axis, higher = denser --
    # so Plotly's marching-cubes reconstructs rounded iso-surfaces (the
    # Geosoft-voxel look) instead of faceting on the sparse grid, then
    # applies a NaN-aware Gaussian whose width scales with that lattice
    # and, for a resistivity-band body, feathers the in/out-of-band
    # cliff so the envelope is a smooth contour, not a voxel staircase.
    volume_smoothing: float = 0.0
    title: str = ""


@dataclass
class ExportOptions:
    """Options used by figure export helpers."""

    path: str | Path
    format: str | None = None
    width: int | None = None
    height: int | None = None
    scale: float = 2.0
    include_plotlyjs: str | bool = "cdn"
