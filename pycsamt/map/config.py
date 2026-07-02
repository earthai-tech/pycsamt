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
    show_profiles: bool = True
    show_contours: bool = False
    contour_image: bool = False
    contour_mode: str = "filled+lines"
    contour_levels: int = 12
    contour_opacity: float = 0.55
    contour_interp: str = "cubic"
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
    station_symbol: str = "diamond"
    station_size: int = 4
    station_color: str = "#1f2937"
    station_labels: bool = False
    aspectmode: str = "data"
    x_unit: str = "m"
    depth_unit: str = "m"
    smooth_sections: bool = True
    section_res: int = 100
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
