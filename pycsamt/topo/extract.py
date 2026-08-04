# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Extract elevation and chainage arrays from Sites / EDI collections.

These utilities are the bridge between pycsamt's data model and the
topography rendering pipeline.  They operate on any object that contains
EDI-like items with HEAD lat/lon/elev attributes — including
:class:`~pycsamt.site.base.Sites`, plain lists of
:class:`~pycsamt.seg.edi.EDIFile`, and
:class:`~pycsamt.stratagem.io.EDIBatch` objects.

Examples
--------
>>> from pycsamt.topo.extract import extract_elevation, extract_chainage
>>> elev = extract_elevation(sites)  # (n_stations,)  m a.s.l.
>>> chain = extract_chainage(sites)  # (n_stations,)  km
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

__all__ = [
    "extract_elevation",
    "extract_chainage",
    "has_elevation",
    "extract_station_names",
]


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def extract_elevation(sites: Any) -> np.ndarray:
    """Extract per-station elevation (m a.s.l.) from a Sites or EDI collection.

    Reads the ``.elev`` (or ``.elevation`` / ``.alt``) field from each
    station's HEAD coordinate block.  Returns a zero array and emits a
    :class:`UserWarning` when no valid non-zero elevation is found.

    Parameters
    ----------
    sites : Sites, EDICollection, list[EDIFile], or single EDIFile
        Any object that contains station data with HEAD lat/lon/elev.

    Returns
    -------
    numpy.ndarray, shape (n_stations,)
        Elevation in metres above sea level, one value per station in
        the order they appear in the collection.
    """
    elevs = [_read_elev(ed) for ed in _iter_edis(sites)]
    if not elevs:
        return np.zeros(0)
    arr = np.array(
        [v if (v is not None and np.isfinite(v)) else 0.0 for v in elevs],
        dtype=float,
    )
    if not np.any(arr != 0.0):
        warnings.warn(
            "No non-zero elevation data found in the station collection. "
            "Topography will appear flat (all zeros).  "
            "Set topography via Sites.with_topography() or "
            "configure_topo(source='array', elev_array=...).",
            UserWarning,
            stacklevel=2,
        )
    return arr


def extract_chainage(sites: Any) -> np.ndarray:
    """Compute along-profile cumulative distance (km) for each station.

    Uses a flat-Earth approximation: cumulative Euclidean distance in
    lat/lon space scaled to metres, converted to km.  Stations are
    assumed to be in profile order.

    Parameters
    ----------
    sites : Sites, EDICollection, list[EDIFile]

    Returns
    -------
    numpy.ndarray, shape (n_stations,)
        Cumulative chainage in **km** from the first station (starts at 0).
    """
    coords = _read_latlon(sites)
    n = len(coords)
    if n == 0:
        return np.zeros(0)
    if n == 1:
        return np.zeros(1)

    lats = np.array([c[0] for c in coords], dtype=float)
    lons = np.array([c[1] for c in coords], dtype=float)

    R_lat = 111_000.0  # m / deg latitude
    R_lon = 111_000.0 * np.cos(np.radians(np.mean(lats)))  # m / deg longitude

    dlat = np.diff(lats) * R_lat
    dlon = np.diff(lons) * R_lon
    segs = np.sqrt(dlat**2 + dlon**2)  # m
    return np.concatenate([[0.0], np.cumsum(segs)]) / 1000.0  # km


def has_elevation(sites: Any) -> bool:
    """Return True if any station carries a meaningful non-zero elevation.

    Parameters
    ----------
    sites : any station container

    Returns
    -------
    bool
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        arr = extract_elevation(sites)
    return bool(np.any(arr != 0.0))


def extract_station_names(sites: Any) -> list[str]:
    """Return station name / ID strings in collection order.

    Parameters
    ----------
    sites : any station container

    Returns
    -------
    list[str]
    """
    names = []
    for i, ed in enumerate(_iter_edis(sites)):
        name = _read_name(ed)
        names.append(name if name else f"S{i:03d}")
    return names


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _iter_edis(sites: Any):
    """Yield individual EDI-like objects from any container type."""
    # pycsamt Sites — has ._items list of Site objects with .edi attribute
    if hasattr(sites, "_items"):
        for site in sites._items:
            ed = getattr(site, "edi", site)
            yield ed
        return

    # EDIBatch / EDICollection — has .edi_objects_ or ._edis or similar
    for attr in ("edi_objects_", "_edis", "edis", "edi_files"):
        container = getattr(sites, attr, None)
        if container is not None:
            try:
                for item in container:
                    yield item
                return
            except TypeError:
                pass

    # Generic iterable (list / tuple of EDIFile / Site)
    try:
        items = list(sites)
        for item in items:
            yield getattr(item, "edi", item)
        return
    except TypeError:
        pass

    # Single EDI-like object
    yield sites


def _get_head(ed: Any) -> Any | None:
    """Return the head section object from any EDI-like object.

    Handles three HEAD storage patterns:
    - attribute access  ``ed.Head`` / ``ed.head`` (old pycsamt / MTpy style)
    - ``sections`` dict  ``ed.sections['head']`` (``seg.edi.EDIFile``)
    - ``DEFINEMEAS`` blocks (some legacy parsers)
    """
    for attr in ("Head", "head", "DEFINEMEAS", "definemeas"):
        h = getattr(ed, attr, None)
        if h is not None:
            return h
    # seg.edi.EDIFile stores sections in a dict
    sections = getattr(ed, "sections", None)
    if isinstance(sections, dict):
        for key in ("head", "Head", "MTSECT", "mtsect"):
            h = sections.get(key)
            if h is not None:
                return h
    return None


def _get_info(ed: Any) -> Any | None:
    """Return the ``>INFO`` section object from a ``seg.edi.EDIFile``.

    Some real deliveries (e.g. non-standard CSAMT exports) record
    ``LATITUDE=``/``LONGITUDE=``/``ELEVATION=`` under the ``>INFO`` block
    instead of the standard SEG-EDI ``>HEAD`` keys :func:`_get_head`
    reads -- confirmed directly on a real survey where ``Site.coords``
    (and this module's own HEAD-only lookup) came back empty even though
    the values are right there in the file, already parsed onto
    :class:`pycsamt.seg.heads.Info`'s own ``.latitude``/``.longitude``/
    ``.elevation`` attributes. This is a real, standards-compliant
    fallback (the SEG-EDI spec permits location metadata in ``>INFO``),
    not a delivery-specific hack.
    """
    get_section = getattr(ed, "get_section", None)
    if callable(get_section):
        try:
            return get_section("info")
        except (KeyError, TypeError, ValueError):
            return None
    sections = getattr(ed, "sections", None)
    if isinstance(sections, dict):
        for key in ("info", "Info", "INFO"):
            info = sections.get(key)
            if info is not None:
                return info
    return None


def _first_finite(obj: Any, names: tuple[str, ...]) -> float | None:
    """Return the first finite float among ``obj``'s ``names`` attributes."""
    for name in names:
        v = getattr(obj, name, None)
        if v is not None:
            try:
                f = float(v)
                if np.isfinite(f):
                    return f
            except (TypeError, ValueError):
                pass
    return None


def _read_elev(ed: Any) -> float | None:
    """Read elevation in metres from an EDI-like object's HEAD or INFO.

    A parsed field that is present but exactly ``0.0`` is treated the
    same as absent and superseded by a later, non-zero source -- the
    standard ``>HEAD`` block's ``Location``/elevation fields default to
    ``0.0`` when the file simply never set them (as opposed to a
    genuine sea-level station), which previously masked a real, non-zero
    value recorded under ``>INFO`` instead (see :func:`_get_info`).
    """
    candidate: float | None = None
    h = _get_head(ed)
    if h is not None:
        loc = getattr(h, "Location", None)
        if loc is not None:
            v = _first_finite(loc, ("elevation", "elev", "alt"))
            if v is not None and v != 0.0:
                return v
            candidate = candidate if candidate is not None else v
        v = _first_finite(h, ("elev", "elevation", "alt", "ALT", "z", "Z"))
        if v is not None and v != 0.0:
            return v
        candidate = candidate if candidate is not None else v
    v = _first_finite(ed, ("elev", "elevation", "alt"))
    if v is not None and v != 0.0:
        return v
    candidate = candidate if candidate is not None else v
    # >INFO fallback (see _get_info's own docstring) -- checked last among
    # non-zero sources, but preferred over an all-zero HEAD candidate.
    info = _get_info(ed)
    if info is not None:
        v = _first_finite(info, ("elevation", "elev", "alt"))
        if v is not None and v != 0.0:
            return v
        candidate = candidate if candidate is not None else v
    return candidate


def _read_latlon(sites: Any) -> list[tuple[float, float]]:
    """Extract (lat, lon) pairs in collection order."""
    coords = []
    for ed in _iter_edis(sites):
        lat = lon = None
        h = _get_head(ed)
        if h is not None:
            # seg.edi.EDIFile — Location sub-object
            loc = getattr(h, "Location", None)
            if loc is not None:
                lat = getattr(loc, "latitude", None) or getattr(
                    loc, "lat", None
                )
                lon = (
                    getattr(loc, "longitude", None)
                    or getattr(loc, "lon", None)
                    or getattr(loc, "long", None)
                )
            # Fall back to direct attributes on head
            if lat is None or lon is None:
                lat = getattr(h, "lat", None)
                lon = getattr(h, "lon", None) or getattr(h, "long", None)
        if lat is None or lon is None:
            # Fall back to the >INFO section (see _get_info's docstring).
            info = _get_info(ed)
            if info is not None:
                lat = lat if lat is not None else getattr(info, "latitude", None)
                lon = lon if lon is not None else getattr(info, "longitude", None)
        if lat is not None and lon is not None:
            try:
                coords.append((float(lat), float(lon)))
            except (TypeError, ValueError):
                pass
    return coords


def _read_name(ed: Any) -> str | None:
    """Read station name / dataid from an EDI-like object."""
    h = _get_head(ed)
    if h is not None:
        for name_attr in ("station", "dataid", "sitename", "name"):
            v = getattr(h, name_attr, None)
            if v is not None:
                return str(v).strip()
    for name_attr in ("station", "station_id", "name", "id"):
        v = getattr(ed, name_attr, None)
        if v is not None:
            return str(v).strip()
    return None
