# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Topography subpackage for pycsamt.

Provides a global configuration singleton (:data:`PYCSAMT_TOPO`) and a
set of utilities to extract elevation data from station collections,
transform flat depth-grids into terrain-following coordinates, and
render terrain polygons on 2-D section plots.

Quick start
-----------
Enable topography for all 2-D plots in the current session::

    from pycsamt.topo import configure_topo
    configure_topo(enabled=True)

Disable / reset::

    from pycsamt.topo import reset_topo
    reset_topo()

Temporarily override with a context manager::

    from pycsamt.topo import PYCSAMT_TOPO
    with PYCSAMT_TOPO.context(enabled=True, exaggeration=2.0):
        fig = section.plot()

Public API
----------
:data:`PYCSAMT_TOPO`
    Global :class:`~pycsamt.topo.config.TopoConfig` singleton.
:func:`configure_topo`
    Set attributes on the global singleton.
:func:`reset_topo`
    Restore all attributes to package defaults.
:func:`extract_elevation`
    Pull per-station elevation arrays from Sites / EDI collections.
:func:`extract_chainage`
    Compute along-profile cumulative distance (km).
:func:`has_elevation`
    Quick boolean check for meaningful elevation data.
:func:`extract_station_names`
    Return station name strings in collection order.
:func:`interp_elev`
    Interpolate station elevations to arbitrary x positions.
:func:`drape_section`
    Build terrain-following 2-D node arrays for pcolormesh.
:func:`mask_above_topo`
    NaN-mask data cells above the terrain surface.
:func:`station_surface_z`
    Draped z-coordinate for station marker pins.
:func:`draw_topo_section`
    Overlay terrain fill + station pins on a depth-section axes.
:func:`draw_topo_strip`
    Add a topo elevation strip above a pseudosection image axes.
:func:`add_station_labels`
    Draw rotated station name labels at marker positions.
"""

from .config import (           # noqa: F401
    TopoConfig,
    PYCSAMT_TOPO,
    configure_topo,
    reset_topo,
)
from .extract import (          # noqa: F401
    extract_elevation,
    extract_chainage,
    has_elevation,
    extract_station_names,
)
from .drape import (            # noqa: F401
    interp_elev,
    drape_section,
    mask_above_topo,
    station_surface_z,
)
from .overlay import (          # noqa: F401
    draw_topo_section,
    draw_topo_strip,
    add_station_labels,
)

__all__ = [
    # config
    "TopoConfig",
    "PYCSAMT_TOPO",
    "configure_topo",
    "reset_topo",
    # extract
    "extract_elevation",
    "extract_chainage",
    "has_elevation",
    "extract_station_names",
    # drape
    "interp_elev",
    "drape_section",
    "mask_above_topo",
    "station_surface_z",
    # overlay
    "draw_topo_section",
    "draw_topo_strip",
    "add_station_labels",
]
