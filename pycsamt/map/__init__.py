# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Code-first mapping API for pyCSAMT surveys.

This package provides reusable station, profile, and 3-D map
builders for scripts, notebooks, and applications.

Functions accept the same flexible inputs used by
``pycsamt.emtools.ensure_sites``: an EDI path, a directory,
an iterable of EDI-like objects, or a ``Sites`` instance.
"""

from __future__ import annotations

from ._core import (
    ComponentSpec,
    FrequencySelection,
    FrequencyValue,
    MapData,
    ProfileLine,
    StationRecord,
    component_spec,
    ensure_map_data,
    normalize_component,
    select_frequency,
    value_at_frequency_details,
)
from .config import (
    ExportOptions,
    MapTheme,
    ProfileMapOptions,
    StationMapOptions,
    VolumeMapOptions,
)
from .export import (
    export_figure,
    figure_to_dict,
    save_png,
    write_dict,
    write_html,
    write_image,
    write_json,
)
from .loader import load_lines, resolve_line_groups
from .overlays import (
    BASEMAP_STYLES,
    BasemapConfig,
    CRSConfig,
    basemap_style_and_layers,
    build_basemap_layout,
    build_contour_overlay,
    build_geo_contour_image,
    build_profile_line_overlay,
    build_station_label_overlay,
    build_topography_overlay,
    interpolate_overlay_grid,
    normalize_epsg,
    reproject_xy_to_lonlat,
    resolve_crs_info,
    transform_xy,
)
from .profile import (
    ProfileMap,
    build_profile_map,
    build_pseudosection,
    plot_profile_map,
    plot_pseudosection,
)
from .station import (
    StationMap,
    build_station_map,
    plot_station_map,
)
from .topo import (
    apply_elevations,
    fetch_elevations,
    parse_elevation_file,
)
from .view import MapView, open_app
from .volume import (
    Map3D,
    VolumeMap,
    build_3d_map,
    plot_3d_map,
    plot_volume_map,
)

__all__ = [
    "BASEMAP_STYLES",
    "BasemapConfig",
    "ComponentSpec",
    "CRSConfig",
    "basemap_style_and_layers",
    "build_geo_contour_image",
    "ExportOptions",
    "FrequencySelection",
    "FrequencyValue",
    "Map3D",
    "MapData",
    "MapTheme",
    "MapView",
    "ProfileLine",
    "ProfileMap",
    "ProfileMapOptions",
    "StationMap",
    "StationMapOptions",
    "StationRecord",
    "VolumeMapOptions",
    "VolumeMap",
    "apply_elevations",
    "build_3d_map",
    "build_basemap_layout",
    "build_contour_overlay",
    "build_profile_line_overlay",
    "build_profile_map",
    "build_pseudosection",
    "build_station_map",
    "build_station_label_overlay",
    "build_topography_overlay",
    "component_spec",
    "ensure_map_data",
    "export_figure",
    "fetch_elevations",
    "figure_to_dict",
    "interpolate_overlay_grid",
    "load_lines",
    "normalize_component",
    "normalize_epsg",
    "open_app",
    "parse_elevation_file",
    "plot_3d_map",
    "plot_profile_map",
    "plot_pseudosection",
    "plot_station_map",
    "plot_volume_map",
    "reproject_xy_to_lonlat",
    "resolve_crs_info",
    "resolve_line_groups",
    "save_png",
    "select_frequency",
    "transform_xy",
    "value_at_frequency_details",
    "write_html",
    "write_dict",
    "write_image",
    "write_json",
]
