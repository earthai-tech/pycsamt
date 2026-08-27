# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Topography helpers — ``per_station`` and ``raster`` PCSF kinds.

``per_station`` wires :class:`~pycsamt.format.schema.TopographyPerStation`
through the existing :mod:`pycsamt.map.topo` elevation sources instead
of adding a third, independent CSV/HDF5/NPZ parser next to the two
that already existed and had quietly drifted apart:
:func:`pycsamt.map.topo.parse_elevation_file` (``_ID_KEYS`` included
``"station_names"``) and
``pycsamt.app.web.callbacks.map3d._parse_topo_upload`` (a hand-rolled
duplicate whose own id-column list did not). This module reuses the
former; ``map3d.py`` was updated to delegate to it (see its own
docstring), so both live-app uploads and PCSF conversion now agree on
exactly one set of recognised column/dataset names.

``raster`` (:class:`~pycsamt.format.schema.TopographyRaster`) is a
standalone gridded-DEM surface, independent of any station table.
:func:`topography_from_grid` builds one from plain ``x``/``y``/
``elevation`` arrays a caller already has — it never parses a
georeferenced raster *file format* (GeoTIFF, ASCII grid, ...) itself,
so no GDAL/rasterio dependency is introduced here: reading such a file
is the caller's own responsibility, done by whatever means it likes,
outside pyCSAMT's dependency chain (e.g. with ``rasterio`` in a
one-off script that then hands PCSF three plain arrays).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from .schema import TopographyPerStation, TopographyRaster

if TYPE_CHECKING:
    from pycsamt.map._core import MapData

__all__ = [
    "topography_from_map_data",
    "topography_from_elevation_file",
    "topography_to_elev_map",
    "topography_from_grid",
    "topography_raster_to_grid",
]


def topography_from_map_data(data: "MapData") -> TopographyPerStation | None:
    """Build topography from a :class:`~pycsamt.map._core.MapData`'s
    own station elevations (typically real, EDI-derived values).

    Parameters
    ----------
    data : MapData
        Survey data, e.g. from :func:`pycsamt.map.load_lines`.

    Returns
    -------
    TopographyPerStation or None
        ``None`` when no station carries a finite elevation, so
        callers can leave :attr:`PCSFModel.topography` unset rather
        than persisting an all-``nan`` table.

    Examples
    --------
    >>> from pycsamt.map import load_lines
    >>> from pycsamt.format.topography import topography_from_map_data
    >>> data = load_lines("data/AMT/WILLY_DATA", detect="folder")  # doctest: +SKIP
    >>> topo = topography_from_map_data(data)  # doctest: +SKIP
    """
    station_id = [s.id for s in data.stations if s.elevation is not None]
    elevation = [
        float(s.elevation) for s in data.stations if s.elevation is not None
    ]
    if not station_id:
        return None
    return TopographyPerStation(
        station_id=station_id, elevation=np.asarray(elevation, dtype=float)
    )


def topography_from_elevation_file(
    content: str | bytes, filename: str
) -> TopographyPerStation | None:
    """Build topography from an uploaded elevation file.

    Thin wrapper around :func:`pycsamt.map.topo.parse_elevation_file`
    (CSV / HDF5 / NPZ, flexible station-id and elevation column/array
    names) — the same parser the "Upload file" elevation source in
    ``pycsamt.app.web`` uses, so a file that works there also works
    here.

    Parameters
    ----------
    content : str or bytes
        A Dash ``dcc.Upload``-style data URI (``"data:...;base64,..."``),
        raw base64 text, or raw bytes — anything
        :func:`~pycsamt.map.topo.parse_elevation_file` already accepts.
    filename : str
        Used only for its extension (``.csv``/``.h5``/``.hdf5``/``.npz``).

    Returns
    -------
    TopographyPerStation or None
        ``None`` when the file cannot be parsed (unrecognised format,
        missing id/elevation column) — matches
        :func:`~pycsamt.map.topo.parse_elevation_file`'s own
        best-effort, non-raising contract.
    """
    from ..map.topo import parse_elevation_file

    elev_map = parse_elevation_file(content, filename)
    if not elev_map:
        return None
    station_id = list(elev_map)
    elevation = np.asarray([elev_map[sid] for sid in station_id], dtype=float)
    return TopographyPerStation(station_id=station_id, elevation=elevation)


def topography_to_elev_map(topo: TopographyPerStation) -> dict[str, float]:
    """Return ``{station_id: elevation}``, the inverse of both builders.

    The same shape :func:`pycsamt.map.topo.apply_elevations` and
    :func:`pycsamt.map.topo.parse_elevation_file` already use, so a
    PCSF file's topography can be applied straight back onto a
    :class:`~pycsamt.map._core.MapData` with no extra conversion.
    """
    return {
        str(sid): float(elev)
        for sid, elev in zip(topo.station_id, topo.elevation)
    }


def topography_from_grid(
    x: np.ndarray, y: np.ndarray, elevation: np.ndarray
) -> TopographyRaster:
    """Build a gridded-DEM :class:`~pycsamt.format.schema.TopographyRaster`.

    Parameters
    ----------
    x : ndarray (n_x,)
        Grid x-coordinates (or longitude), increasing.
    y : ndarray (n_y,)
        Grid y-coordinates (or latitude), increasing.
    elevation : ndarray (n_y, n_x)
        Elevation surface, metres, on the ``(y, x)`` meshgrid implied
        by *x*/*y* (``numpy.meshgrid(x, y)``'s default row-major
        convention).

    Returns
    -------
    TopographyRaster

    Examples
    --------
    >>> import numpy as np
    >>> from pycsamt.format.topography import topography_from_grid
    >>> x = np.linspace(0.0, 500.0, 6)
    >>> y = np.linspace(0.0, 300.0, 4)
    >>> elevation = 100.0 + 0.01 * np.add.outer(y, x)
    >>> topo = topography_from_grid(x, y, elevation)
    >>> topo.elevation.shape
    (4, 6)
    """
    return TopographyRaster(
        x=np.asarray(x, dtype=float),
        y=np.asarray(y, dtype=float),
        elevation=np.asarray(elevation, dtype=float),
    )


def topography_raster_to_grid(
    topo: TopographyRaster,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ``(x, y, elevation)``, the inverse of :func:`topography_from_grid`."""
    return topo.x, topo.y, topo.elevation
