# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from ..constants import _M_PER_DEG
from .location import Coord, chainage_along
from .utils import get_coords, station_name

__all__ = ["infer_line_orientation", "Profile"]


@dataclass
class Profile:
    r"""
    Profile(origin, azimuth, chainages=None, spacing_stats=None,
            gaps=None)

    Describe a 1-D survey line (profile) built from site
    locations. Stores the origin, line azimuth, per-site
    chainages, spacing statistics, and detected large gaps.

    The azimuth follows the convention 0 deg north, 90 deg east.
    Chainages are in meters and are computed consistently with
    :func:`pycsamt.site.location.chainage_along`.

    Parameters
    ----------
    origin : pycsamt.site.location.Coord
        Origin coordinate used for chainage computations.
    azimuth : float
        Line azimuth in degrees, 0 deg north, 90 deg east.
    chainages : dict, optional
        Mapping ``{station_name: chainage_m}``. Usually filled
        via :meth:`Profile.from_sites`.
    spacing_stats : dict, optional
        Precomputed spacing metrics. Filled automatically when
        chainages are set or updated.
    gaps : list of tuple, optional
        Large spacing gaps as ``[(s_left, s_right), ...]`` in
        meters.

    Attributes
    ----------
    origin : Coord
        Profile origin coordinate.
    azimuth : float
        Profile azimuth (deg).
    chainages : dict
        Per-site chainages in meters.
    spacing_stats : dict
        Keys include ``spacing_mean``, ``spacing_med``,
        ``spacing_min``, ``spacing_max`` (meters).
    gaps : list of tuple
        Detected large gaps as chainage intervals (meters).

    Notes
    -----
    Chainage for a site with local offsets :math:`(x,y)` relative
    to the origin and profile azimuth :math:`A` is

    .. math::

       s = x * cos(A) + y * sin(A)

    Examples
    --------
    >>> from pycsamt.site.profile import Profile
    >>> from pycsamt.site.location import Coord
    >>> class Head:
    ...     def __init__(self, lat, lon, name):
    ...         self.lat, self.lon, self.dataid = lat, lon, name
    ...
    >>> class EDI:
    ...     def __init__(self, name, lat, lon):
    ...         self._h = Head(lat, lon, name)
    ...     def get_section(self, key):
    ...         return self._h if key == "head" else None
    ...
    >>> sites = [EDI("A", 0.0, 0.00),
    ...          EDI("B", 0.0, 0.01),
    ...          EDI("C", 0.0, 0.02)]
    >>> prof = Profile.from_sites(sites)
    >>> round(prof.azimuth) in (89, 90, 91)
    True

    See Also
    --------
    pycsamt.site.profile.infer_line_orientation
    pycsamt.site.location.Coord
    pycsamt.site.location.chainage_along
    """

    origin: Coord
    azimuth: float
    chainages: dict[str, float] = field(default_factory=dict)
    spacing_stats: dict[str, float] = field(default_factory=dict)
    gaps: list[tuple[float, float]] = field(default_factory=list)

    @classmethod
    def from_sites(
        cls,
        sites: Iterable[Any],
        *,
        origin: Coord | None = None,
        azimuth: float | None = None,
    ) -> Profile:
        r"""
        from_sites(sites, *, origin=None, azimuth=None)

        Build a profile from an iterable of sites. If ``origin`` is
        omitted, the first site with finite coordinates is used. If
        ``azimuth`` is omitted, it is inferred with
        :func:`infer_line_orientation`.

        Parameters
        ----------
        sites : iterable
            EDI-like objects or wrappers with an ``.edi`` attribute.
        origin : Coord, optional
            Origin coordinate. If omitted, inferred from the first
            finite site.
        azimuth : float, optional
            Profile azimuth in degrees. If omitted, inferred from
            site positions.

        Returns
        -------
        Profile
            Profile with per-site chainages and spacing statistics
            computed.

        Notes
        -----
        Coordinates and names are obtained via
        :func:`pycsamt.site.utils.get_coords` and
        :func:`pycsamt.site.utils.station_name`. Chainages are
        computed with
        :func:`pycsamt.site.location.chainage_along`.

        Examples
        --------
        >>> from pycsamt.site.profile import Profile
        >>> prof = Profile.from_sites([])
        >>> isinstance(prof, Profile)
        True
        """

        items = _iter_sites(sites)
        if not items:
            o = origin or Coord(0.0, 0.0, 0.0)
            return cls(o, float(azimuth or 0.0))

        # origin default = first finite coordinate
        if origin is None:
            for _, la, lo, _ in items:
                if _finite(la) and _finite(lo):
                    origin = Coord(la, lo, 0.0)
                    break
        if origin is None:
            origin = Coord(0.0, 0.0, 0.0)

        # auto azimuth from points if not given
        if azimuth is None:
            azimuth = infer_line_orientation([it[3] for it in items])

        # chainages along the line
        s: dict[str, float] = {}
        for name, la, lo, _ in items:
            if _finite(la) and _finite(lo):
                s[name] = float(
                    chainage_along(
                        (origin.lat, origin.lon),
                        float(azimuth),
                        (la, lo),
                    )
                )

        prof = cls(origin, float(azimuth), s)
        prof._update_stats()
        return prof

    def sort_sites(self, sites: Iterable[Any]) -> list[Any]:
        r"""
        sort_sites(sites)

        Return the input sites ordered by chainage along the
        profile. Sites without finite chainage are dropped.

        Parameters
        ----------
        sites : iterable
            Same accepted types as :meth:`Profile.from_sites`.

        Returns
        -------
        list
            The subset of input sites sorted by increasing chainage.

        Examples
        --------
        >>> sorted_sites = prof.sort_sites(sites)
        >>> isinstance(sorted_sites, list)
        True
        """

        items = _iter_sites(sites)
        order = []
        for it in items:
            si = self.chainages.get(it[0], float("nan"))
            order.append((si, it[3]))
        order.sort(key=lambda t: (not _finite(t[0]), t[0]))
        return [o[1] for o in order if _finite(o[0])]

    def slice(self, s_min: float, s_max: float) -> dict[str, float]:
        r"""
        slice(s_min, s_max)

        Return chainages within the window ``s_min <= s <= s_max``
        as a dict ordered by chainage.

        Parameters
        ----------
        s_min : float
            Lower bound in meters.
        s_max : float
            Upper bound in meters.

        Returns
        -------
        dict
            Mapping ``{station_name: chainage_m}``, ordered by
            chainage.

        Examples
        --------
        >>> win = prof.slice(500.0, 1500.0)
        >>> isinstance(win, dict)
        True
        """
        out = {
            k: v
            for k, v in self.chainages.items()
            if _finite(v) and (s_min <= v <= s_max)
        }
        return dict(sorted(out.items(), key=lambda t: t[1]))

    def resample(self, step: float) -> np.ndarray:
        r"""
        resample(step)

        Build a regular chainage grid between the current minimum
        and maximum chainage (inclusive of the minimum). If
        ``step <= 0``, an empty array is returned.

        Parameters
        ----------
        step : float
            Grid spacing in meters.

        Returns
        -------
        numpy.ndarray
            1-D array of chainage locations in meters.

        Examples
        --------
        >>> grid = prof.resample(250.0)
        >>> (grid.ndim, grid.dtype.kind) == (1, 'f')
        True
        """

        if not self.chainages:
            return np.asarray([], float)
        s = np.asarray(list(self.chainages.values()), float)
        s = s[np.isfinite(s)]
        if s.size == 0 or step <= 0:
            return np.asarray([], float)
        a = float(np.nanmin(s))
        b = float(np.nanmax(s))
        n = int(max(1, math.floor((b - a) / step) + 1))
        return a + step * np.arange(n, dtype=float)

    def summary(self) -> dict[str, float]:
        r"""
        summary()

        Return a compact summary of the profile and its spacing
        statistics.

        Returns
        -------
        dict
            Keys include ``n_sites``, ``s_min``, ``s_max``,
            ``n_gaps``, and entries from ``spacing_stats``
            (``spacing_mean``, ``spacing_med``, ``spacing_min``,
            ``spacing_max``).

        Examples
        --------
        >>> info = prof.summary()
        >>> set(["n_sites", "n_gaps"]).issubset(info.keys())
        True
        """
        d: dict[str, float] = {}
        d.update(self.spacing_stats)
        d["n_sites"] = float(len(self.chainages))
        if self.chainages:
            vals = np.asarray(list(self.chainages.values()), float)
            m = np.isfinite(vals)
            if np.any(m):
                d["s_min"] = float(np.min(vals[m]))
                d["s_max"] = float(np.max(vals[m]))
        d["n_gaps"] = float(len(self.gaps))
        return d

    def _update_stats(self) -> None:
        if not self.chainages:
            self.spacing_stats = {}
            self.gaps = []
            return
        s = np.asarray(list(self.chainages.values()), float)
        s = s[np.isfinite(s)]
        s.sort()
        if s.size < 2:
            self.spacing_stats = {
                "spacing_mean": float("nan"),
                "spacing_med": float("nan"),
                "spacing_min": float("nan"),
                "spacing_max": float("nan"),
            }
            self.gaps = []
            return
        d = np.diff(s)
        self.spacing_stats = {
            "spacing_mean": float(np.mean(d)),
            "spacing_med": float(np.median(d)),
            "spacing_min": float(np.min(d)),
            "spacing_max": float(np.max(d)),
        }
        med = float(np.median(d))
        thr = 1.5 * med if med > 0 else float("inf")
        gaps: list[tuple[float, float]] = []
        for i in range(d.size):
            if d[i] > thr:
                gaps.append((float(s[i]), float(s[i + 1])))
        self.gaps = gaps


# ------------------------ helpers (minimal) ------------------------

def _finite(x: float) -> bool:
    return math.isfinite(float(x))


def _iter_sites(sites: Iterable[Any]) -> list[tuple[str, float, float, Any]]:
    """Return (name, lat, lon, edi_like) for each input item.

    Uses utils.station_name / utils.get_coords and supports objects
    that wrap an EDI as `.edi`.
    """
    out: list[tuple[str, float, float, Any]] = []
    for s in sites:
        ed = getattr(s, "edi", s)
        name = station_name(ed) or station_name(s)
        c = get_coords(ed)  # has fields: lat, lon, elev
        out.append((name, float(c.lat), float(c.lon), ed))
    return out


def _xy_local(lats: np.ndarray, lons: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    la0 = float(np.nanmean(lats))
    lo0 = float(np.nanmean(lons))
    dy = (lats - la0) * _M_PER_DEG
    dx = (lons - lo0) * _M_PER_DEG * math.cos(math.radians(la0))
    return dx, dy


def infer_line_orientation(sites: Iterable[Any]) -> float:
    r"""
    Infer the survey line azimuth from a collection of sites.

    This estimates the dominant line axis that best explains the
    site distribution. The estimate uses PCA on local Cartesian
    offsets derived from geographic coordinates. The result is an
    azimuth in degrees where 0 deg is north and 90 deg is east.

    Parameters
    ----------
    sites : iterable of objects
        Items may be EDI-like objects or wrappers exposing an
        ``.edi`` attribute. Site names and coordinates are
        resolved using :func:`pycsamt.site.utils.station_name`
        and :func:`pycsamt.site.utils.get_coords`.

    Returns
    -------
    float
        Azimuth in degrees, 0 deg north and 90 deg east. Because
        eigenvectors are sign-ambiguous, the orientation is
        defined modulo 180 deg (e.g., 45 deg and 225 deg are the
        same axis).

    Notes
    -----
    Local offsets :math:`(x, y)` are built with a small-extent
    flat-Earth approximation about the mean latitude
    :math:`lat_0` and longitude :math:`lon_0`:

    .. math::

       x = (lon - lon_0) * M\_PER\_DEG * cos(lat\_0)

       y = (lat - lat_0) * M\_PER\_DEG

    The principal component with the largest variance gives the
    line axis in the local frame :math:`(x=east, y=north)`. It is
    then converted to an azimuth (0 deg north, 90 deg east).

    Examples
    --------
    >>> from pycsamt.site.profile import infer_line_orientation
    >>> class Head:
    ...     def __init__(self, lat, lon, name):
    ...         self.lat, self.lon, self.dataid = lat, lon, name
    ...
    >>> class EDI:
    ...     def __init__(self, name, lat, lon):
    ...         self._h = Head(lat, lon, name)
    ...     def get_section(self, key):
    ...         return self._h if key == "head" else None
    ...
    >>> east = [EDI(f"S{i}", 0.0, i*0.01) for i in range(5)]
    >>> az = infer_line_orientation(east)
    >>> 80.0 <= az <= 100.0
    True

    See Also
    --------
    pycsamt.site.profile.Profile
    pycsamt.site.location.Coord
    pycsamt.site.location.chainage_along
    pycsamt.site.utils.get_coords
    pycsamt.site.utils.station_name

    References
    ----------
    .. [1] Jolliffe, I. T., "Principal Component Analysis",
           Springer, 2nd ed., 2002.
    """

    items = _iter_sites(sites)
    if not items:
        return 0.0
    la = np.asarray([i[1] for i in items], float)
    lo = np.asarray([i[2] for i in items], float)
    m = np.isfinite(la) & np.isfinite(lo)
    if not np.any(m):
        return 0.0
    x, y = _xy_local(la[m], lo[m])
    if x.size < 2:
        return 0.0

    # PCA major axis
    X = np.stack([x, y], axis=0)
    vals, vecs = np.linalg.eig(np.cov(X))
    k = int(np.argmax(vals))
    vx, vy = float(vecs[0, k]), float(vecs[1, k])

    # Convert to azimuth (0°=north, 90°=east).
    # The vector is in (x=east, y=north) => swap in atan2.
    ang = math.degrees(math.atan2(vx, vy))
    return (ang + 360.0) % 360.0
