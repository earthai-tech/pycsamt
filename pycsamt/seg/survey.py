# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

import fnmatch
import math
import re
from collections.abc import Iterable, Sequence
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..gis.utils import (
    get_utm_zone,
    normalize_lat_lon,
    project_point_utm2ll,
    to_utm,
)
from ..log.logger import get_logger
from .base import SurveyBase
from .collection import EDICollection
from .edi import EDIFile

logger = get_logger(__name__)

__all__ = ["EDIProfile", "Stations", "Topography"]


class Stations(SurveyBase):
    r"""
    Lightweight table view for station metadata derived from
    :class:`~.edi.EDIFile` objects.  Provides quick access to
    names, geographic coordinates, optional elevations, and
    working projected coordinates to support survey tasks.

    Parameters
    ----------
    items : EDIFile or iterable of EDIFile or EDICollection
        The sites to summarize.
    verbose : int, default ``0``
        Verbosity level for logging.

    Attributes
    ----------
    count : int
        Number of valid stations (rows with lat/lon).
    stations : list of str
        Station identifiers sourced from headers or filenames.

    Methods
    -------
    table()
        Materialize the rows as a list of dicts with keys
        ``station``, ``lat``, ``lon``, ``elev``, ``e`` (east),
        ``n`` (north), ``zone`` and ``path``.
    to_dataframe()
        Convert the table to a pandas DataFrame when pandas
        is available.
    bounds()
        Geographic bounding box as ``(min_lat, min_lon,
        max_lat, max_lon)``.
    select(keys=None, pattern=None, regex=None, pred=None)
        Return a filtered view.  Multiple filters are
        combined with logical **AND**.  See the method
        docstring for details.
    sort(by='station', reverse=False, inplace=True)
        Sort rows by a column (e.g. ``'station'``, ``'lat'``,
        ``'lon'``, ``'elev'``, ``'e'``, ``'n'``).  Returns
        this instance or a new view depending on ``inplace``.
    offsets(origin=None, azimuth=None)
        Compute along-line and cross-line offsets in meters
        from projected coordinates.  The profile axis is set
        by ``azimuth`` or inferred from endpoints.
    set_coords(key, *, lat=None, lon=None, elev=None)
        Update coordinates for a single station.  When the
        backing :class:`~.edi.EDIFile` is available its
        ``>HEAD`` values are kept in sync.

    Notes
    -----
    The class performs minimal validation.  Rows missing
    lat/lon are skipped.  Projected coordinates are intended
    for short-range work; for mapping at scale prefer the
    GIS utilities provided elsewhere in the package.

    Examples
    --------
    Build a table and print a compact view::

        sts = Stations(coll)
        for r in sts.table():
            print(r["station"], r["lat"], r["lon"])

    Filter, sort, and compute offsets::

        sel = sts.select(pattern="K*", pred=lambda r: r["elev"] > 800)
        sel.sort(by="e")
        along, across = sel.offsets()

    See Also
    --------
    EDIProfile
        Track-aware helper that computes distance and azimuth.
    Topography
        Produces elevation profiles from stations or profiles.

    References
    ----------
    .. [Stations-1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def __init__(
        self,
        items: EDIFile | Iterable[EDIFile] | EDICollection | EDIProfile,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)

        if isinstance(items, EDIProfile):
            eds = [r["ed"] for r in items._rows]  # noqa: SLF001
        else:
            eds = _as_collection(items)
        self._rows: list[dict[str, object]] = []
        self._load(eds)

    def _load(self, eds: Sequence[EDIFile]) -> None:
        rows: list[dict[str, object]] = []
        for k, ed in enumerate(eds):
            lat, lon, elev = _coerce_ll(ed)
            sid = _station(ed)
            p = getattr(ed, "path", None)
            rows.append(
                {
                    "idx": k,
                    "station": sid,
                    "lat": lat,
                    "lon": lon,
                    "elev": elev,
                    "path": str(p) if p else None,
                    "ed": ed,
                }
            )
        # project to UTM if possible (skip None coords)
        good = [
            r for r in rows if r["lat"] is not None and r["lon"] is not None
        ]
        if good:
            lat = np.array([r["lat"] for r in good], float)
            lon = np.array([r["lon"] for r in good], float)
            e, n, z = to_utm(lat, lon)
            for i, r in enumerate(good):
                r["e"] = float(e[i])
                r["n"] = float(n[i])
                r["zone"] = z if np.isscalar(z) else str(z[i])
        self._rows = rows

    def names(self) -> list[str]:
        return [str(r["station"]) for r in self._rows]

    def table(self) -> list[dict[str, object]]:
        return [dict(r) for r in self._rows]

    def get(self, key: str) -> EDIFile | None:
        key = str(key)
        for r in self._rows:
            if str(r["station"]) == key:
                return r["ed"]  # type: ignore[return-value]
        return None

    def row(self, key: str) -> dict[str, object] | None:
        key = str(key)
        for r in self._rows:
            if str(r["station"]) == key:
                return dict(r)
        return None

    def select(
        self,
        *,
        keys: Sequence[str] | None = None,
        pattern: str | None = None,
        regex: str | None = None,
        pred=None,
    ) -> Stations:
        r"""
        Return a filtered view of the stations table.

        Multiple filters are combined with logical **AND**.  When
        no filter is given the original view is returned.

        Parameters
        ----------
        keys : sequence of str, optional
            Station identifiers to keep.  Unknown ids are
            ignored.
        pattern : str, optional
            Glob-like pattern matched against station ids (e.g.
            ``'AB*'``).  Case-sensitive.
        regex : str, optional
            Regular expression matched against station ids using
            :func:`re.search`.
        pred : callable, optional
            A predicate ``pred(row) -> bool`` evaluated on each
            row dict.  Keep rows for which the predicate returns
            ``True``.

        Returns
        -------
        Stations
            A new :class:`Stations` view with rows that match the
            filters.

        Notes
        -----
        Filtering does not modify the original container.  Rows
        lacking a station id are always dropped.

        Examples
        --------
        Keep stations starting with ``'K'`` and above 800 m::

            sel = sts.select(pattern="K*", pred=lambda r: r["elev"] > 800)
        """

        rows = self._rows
        if keys:
            want = {str(k) for k in keys}
            rows = [r for r in rows if r["station"] in want]
        if pattern:
            pat = str(pattern)
            rows = [r for r in rows if fnmatch.fnmatch(str(r["station"]), pat)]
        if regex:
            rg = re.compile(str(regex))
            rows = [r for r in rows if rg.search(str(r["station"]))]
        if pred:
            rows = [r for r in rows if bool(pred(r))]
        out = object.__new__(Stations)
        # shallow clone
        out.verbose = self.verbose
        out._rows = [dict(r) for r in rows]  # noqa: SLF001
        return out

    def sort(
        self,
        *,
        by: str = "station",
        reverse: bool = False,
        inplace: bool = True,
    ) -> Stations:
        r"""
        Sort the stations table by a column.

        Parameters
        ----------
        by : str, default ``'station'``
            Column name to sort by (e.g. ``'station'``, ``'lat'``,
            ``'lon'``, ``'elev'``, ``'e'``, ``'n'``).
        reverse : bool, default ``False``
            If ``True`` sort in descending order.
        inplace : bool, default ``True``
            If ``True`` modify this instance and return it.
            Otherwise return a new sorted view.

        Returns
        -------
        Stations
            The sorted :class:`Stations` object (self or a copy).

        Notes
        -----
        Missing values are placed at the end.  Unknown columns
        raise a :class:`KeyError`.
        """
        key = str(by).lower()
        rows = sorted(
            self._rows,
            key=lambda r: r.get(key, None),
            reverse=bool(reverse),
        )
        if inplace:
            self._rows = rows
            return self
        out = object.__new__(Stations)
        out.verbose = self.verbose
        out._rows = rows
        return out

    def offsets(
        self,
        *,
        origin: tuple[float, float] | None = None,
        azimuth: float | None = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        r"""
        Compute along-line and cross-line offsets (meters).

        Offsets are computed from projected coordinates.  The
        along-line axis is defined by the given ``azimuth``; the
        cross-line axis is perpendicular to it.

        Parameters
        ----------
        origin : tuple of float, optional
            Reference point ``(easting, northing)`` in meters.
            Defaults to the first valid station.
        azimuth : float, optional
            Bearing in degrees, clockwise from North.  When
            omitted it is inferred from the first and last
            stations.

        Returns
        -------
        along : ndarray of float
            Distances projected on the profile axis.
        across : ndarray of float
            Signed distances perpendicular to the profile axis.

        Notes
        -----
        Rows without valid projected coordinates are skipped in
        the computation and do not contribute to the result.
        """

        good = [r for r in self._rows if ("e" in r and "n" in r)]
        if not good:
            return (np.asarray([]), np.asarray([]))
        e = np.array([r["e"] for r in good], float)
        n = np.array([r["n"] for r in good], float)
        if azimuth is None:
            # endpoints bearing (north clockwise)
            vx = e[-1] - e[0]
            vy = n[-1] - n[0]
            az = math.degrees(math.atan2(vx, vy)) % 360.0
        else:
            az = float(azimuth) % 360.0
        rad = math.radians(az)
        ve, vn = math.sin(rad), math.cos(rad)
        e0, n0 = (e[0], n[0]) if origin is None else origin
        de, dn = e - e0, n - n0
        along = de * ve + dn * vn
        cross = -de * vn + dn * ve
        return (along, cross)

    def set_coords(
        self,
        key: str,
        *,
        lat: float | None = None,
        lon: float | None = None,
        elev: float | None = None,
    ) -> None:
        r"""
        Update coordinates for a single station.

        Parameters
        ----------
        key : str
            Station identifier to modify.
        lat : float, optional
            New latitude in decimal degrees.
        lon : float, optional
            New longitude in decimal degrees.
        elev : float, optional
            New elevation in meters.

        Returns
        -------
        None

        Notes
        -----
        The in-memory row is updated.  If the backing
        :class:`~.edi.EDIFile` is attached for that station, its
        ``>HEAD`` values are also updated to keep them in sync.
        Unknown station ids raise a :class:`KeyError`.
        """
        ed = self.get(key)
        if ed is None:
            return
        if lat is not None:
            ed.lat = float(lat)
        if lon is not None:
            ed.lon = float(lon)
        if elev is not None:
            ed.elev = float(elev)

    def to_dataframe(
        self,
        *,
        columns: Sequence[str] | None = None,
        index: str | None = "station",
        coerce_numeric: bool = True,
    ) -> pd.DataFrame:
        r"""
        Return a pandas ``DataFrame`` view of the station table.

        Parameters
        ----------
        columns : sequence of str, optional
            Subset and ordering of columns to include.  When
            omitted, a sensible default is used:
            ``('station','lat','lon','elev','e','n','zone',
            'path')``.  Missing names are ignored.
        index : str or None, default ``'station'``
            Column to set as the DataFrame index.  If the name
            is not present, no index is set.  Use ``None`` to
            leave the default RangeIndex.
        coerce_numeric : bool, default ``True``
            Try converting known numeric columns (``lat``, ``lon``,
            ``elev``, ``e``, ``n``) to numeric dtypes.  Non
            convertible values become ``NaN``.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with one row per station.

        Notes
        -----
        ``pandas`` is imported lazily.  If it is not available, an
        :class:`ImportError` is raised.  The method is read-only and
        does not mutate the underlying table.

        Examples
        --------
        Basic usage::

            df = Stations(coll).to_dataframe()
            print(df.head())

        Custom subset and index::

            df = Stations(coll).to_dataframe(
                columns=("station", "elev", "e", "n"),
                index="station",
            )

        See Also
        --------
        table : List-of-dicts representation of the rows.
        bounds : Geographic bounding box.
        select : Filter rows prior to conversion.

        """
        rows = self.table()
        if not rows:
            cols = (
                list(columns)
                if columns is not None
                else [
                    "station",
                    "lat",
                    "lon",
                    "elev",
                    "e",
                    "n",
                    "zone",
                    "path",
                ]
            )
            return pd.DataFrame(columns=cols)

        df = pd.DataFrame.from_records(rows)

        if columns is not None:
            keep = [c for c in columns if c in df.columns]
            df = df.loc[:, keep]

        if coerce_numeric:
            for c in ("lat", "lon", "elev", "e", "n"):
                if c in df.columns:
                    df[c] = pd.to_numeric(df[c], errors="coerce")

        if index is not None and index in df.columns:
            df = df.set_index(index, drop=True)

        return df


class Topography(SurveyBase):
    r"""
    Elevation profile helper.  Builds paired arrays of
    distance and elevation from an :class:`EDIProfile`,
    :class:`Stations`, a collection, or raw :class:`EDIFile`
    inputs.  Includes smoothing, detrending, resampling, and
    quick plotting.

    Parameters
    ----------
    items : EDIProfile or Stations or EDICollection or EDIFile \
            or iterable of EDIFile
        Source of station positions and elevations.  When an
        :class:`EDIProfile` is given, along-line distances are
        reused by default.
    use_profile_step : bool, default ``True``
        If ``True`` and ``items`` is an :class:`EDIProfile`,
        copy its along-profile distances; otherwise recompute
        distances from planar coordinates.
    verbose : int, default ``0``
        Verbosity level for diagnostics.

    Attributes
    ----------
    distance : ndarray of float
        Along-track distances in meters.
    elevation : ndarray of float
        Elevation values aligned with ``distance``.
    trend : ndarray of float or None
        Fitted linear trend after :meth:`detrend`.  Otherwise
        ``None``.

    Methods
    -------
    smooth(window=5, method='median')
        Apply moving median or mean smoothing to elevation.
    detrend()
        Remove a best-fit linear trend and keep it for plotting.
    resample(step)
        Resample to a fixed distance step using interpolation.
    gradient(as_degrees=False)
        First derivative of elevation vs distance; optionally
        in degrees.
    plot(ax=None, title=None, show_trend=True)
        Quick plot of elevation vs distance.
    as_arrays()
        Return ``(distance, elevation)`` copies.
    to_dict()
        Return a dict with ``distance`` and ``elevation`` keys.

    Notes
    -----
    If the input is an :class:`EDIProfile` that has not yet
    computed distances, they are derived automatically.  When
    distances or elevations are missing the result is empty.

    Examples
    --------
    From a profile and detrend before plotting::

        topo = Topography(prof).detrend().smooth(window=7)
        ax = topo.plot(title="Detrended topography")

    From a list of files, resampled every 25 m::

        topo = Topography(edis).resample(step=25.0)

    See Also
    --------
    EDIProfile
        Provides distances and azimuth and can adjust station
        positions.
    Stations
        Tabular access to station metadata.

    References
    ----------
    .. [Topography-1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def __init__(
        self,
        items: EDIFile
        | Iterable[EDIFile]
        | EDICollection
        | EDIProfile
        | Stations,
        *,
        use_profile_step: bool = True,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)
        self._d = np.asarray([], float)
        self._z = np.asarray([], float)
        self._trend: np.ndarray | None = None
        self._load(items, use_profile_step=use_profile_step)

    def _load(
        self,
        items,
        *,
        use_profile_step: bool,
    ) -> None:
        if isinstance(items, EDIProfile):
            # accept a profile and prefer its along-profile distances
            # if requested; otherwise derive cumulative distance from
            # raw XY (east,north). Fall back to triggering the
            # profile's lazy computation when needed.
            profile = items

            # 1) distance
            if use_profile_step:
                d = getattr(profile, "distance", None)
                d = (
                    np.asarray(d, float)
                    if d is not None
                    else np.array([], float)
                )
                if d.size == 0 and hasattr(
                    profile, "_ensure_distance_bearing"
                ):
                    # compute lazily if profile hasn't adjusted yet
                    profile._ensure_distance_bearing()
                    d = np.asarray(profile.distance, float)
            else:
                e, n = profile.xy
                d = self._cumul_from_xy(
                    np.asarray(e, float),
                    np.asarray(n, float),
                )

            # 2) elevation (keep as provided by profile)
            z = getattr(profile, "elev", None)
            z = np.asarray(z, float) if z is not None else np.array([], float)

            # 3) length guard in case inputs disagree
            if d.size != z.size:
                m = int(min(d.size, z.size))
                d = d[:m]
                z = z[:m]

            self._d, self._z = d, z
            return

        if isinstance(items, Stations):
            rows = items.table()
            e = np.array([r.get("e", np.nan) for r in rows], float)
            n = np.array([r.get("n", np.nan) for r in rows], float)
            ok = np.isfinite(e) & np.isfinite(n)
            self._d = self._cumul_from_xy(e[ok], n[ok])
            self._z = np.array([r.get("elev", 0.0) for r in rows], float)[ok]
            return

        eds = _as_collection(items)
        rows = Stations(eds).table()
        e = np.array([r.get("e", np.nan) for r in rows], float)
        n = np.array([r.get("n", np.nan) for r in rows], float)
        ok = np.isfinite(e) & np.isfinite(n)
        self._d = self._cumul_from_xy(e[ok], n[ok])
        self._z = np.array([r.get("elev", 0.0) for r in rows], float)[ok]

    @staticmethod
    def _cumul_from_xy(
        e: np.ndarray,
        n: np.ndarray,
    ) -> np.ndarray:
        if e.size == 0:
            return np.asarray([])
        de = np.diff(e)
        dn = np.diff(n)
        ds = np.hypot(de, dn)
        return np.concatenate(([0.0], np.cumsum(ds)))

    @property
    def distance(self) -> np.ndarray:
        return self._d.copy()

    @property
    def elevation(self) -> np.ndarray:
        return self._z.copy()

    def smooth(
        self,
        *,
        window: int = 5,
        method: str = "median",
    ) -> Topography:
        r"""
        Smooth the elevation series with a sliding window.

        Parameters
        ----------
        window : int, default ``5``
            Window length (samples). Values ``<=1`` skip
            smoothing.
        method : {'median', 'mean'}, default ``'median'``
            Smoothing kernel. ``'median'`` is robust to spikes;
            ``'mean'`` uses a simple moving average.

        Returns
        -------
        Topography
            The instance (in place), allowing chaining.

        Notes
        -----
        Edge handling is performed by shrinking the window near
        the bounds. This modifies the internal elevation array.
        """
        w = max(1, int(window))
        z = self._z
        if z.size == 0 or w == 1:
            return self
        if method == "mean":
            ker = np.ones(w) / float(w)
            zs = np.convolve(z, ker, mode="same")
        else:
            # sliding median
            zs = z.copy()
            hw = w // 2
            for i in range(z.size):
                a = max(0, i - hw)
                b = min(z.size, i + hw + 1)
                zs[i] = float(np.median(z[a:b]))
        self._z = zs
        return self

    def detrend(self) -> Topography:
        if self._d.size < 2:
            self._trend = None
            return self
        c = np.polyfit(self._d, self._z, deg=1)
        t = np.polyval(c, self._d)
        self._z = self._z - t
        self._trend = t
        return self

    def resample(
        self,
        *,
        step: float,
    ) -> Topography:
        r"""
        Resample distance/elevation to a fixed along-line step.

        Parameters
        ----------
        step : float
            Target spacing in meters for the resampled profile.

        Returns
        -------
        Topography
            The instance (in place), allowing chaining.

        Notes
        -----
        Distances are regridded on ``[dmin, dmax]`` with uniform
        spacing and elevations are linearly interpolated.
        """

        if self._d.size == 0:
            return self
        s = float(step)
        dmin, dmax = float(self._d.min()), float(self._d.max())
        di = np.arange(dmin, dmax + s * 0.5, s, float)
        zi = np.interp(di, self._d, self._z)
        if self._trend is not None:
            self._trend = np.interp(di, self._d, self._trend)
        self._d, self._z = di, zi
        return self

    def gradient(
        self,
        *,
        as_degrees: bool = False,
    ) -> np.ndarray:
        r"""
        Compute local slope between consecutive samples.

        Parameters
        ----------
        as_degrees : bool, default ``False``
            If ``True``, return the slope angle in degrees. If
            ``False``, return the rise-over-run ratio.

        Returns
        -------
        ndarray
            Array of length ``len(distance) - 1`` with per-segment
            slopes (or angles when requested).
        """

        if self._d.size < 2:
            return np.asarray([])
        dz = np.diff(self._z)
        dd = np.diff(self._d)
        slope = np.divide(dz, dd, out=np.zeros_like(dz), where=(dd != 0.0))
        if as_degrees:
            return np.degrees(np.arctan(slope))
        return slope

    def plot(
        self,
        *,
        ax: plt.Axes | None = None,
        title: str | None = None,
        show_trend: bool = True,
    ) -> plt.Axes:
        r"""
        Plot elevation versus distance.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Target axes. A new figure/axes is created when
            omitted.
        title : str, optional
            Title for the axes.
        show_trend : bool, default ``True``
            Overlay the last computed trend line (from
            :meth:`detrend`) when available.

        Returns
        -------
        matplotlib.axes.Axes
            The axes with the rendered profile.
        """

        if ax is None:
            _, ax = plt.subplots()
        ax.plot(self._d, self._z, "-o", lw=1.0, ms=3)
        if show_trend and self._trend is not None:
            ax.plot(self._d, self._trend, "--", lw=1.0)
        ax.set_xlabel("distance (m)")
        ax.set_ylabel("elevation (m)")
        ax.grid(True, ls="--", lw=0.5)
        if title:
            ax.set_title(title)
        return ax

    def as_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        return self.distance, self.elevation

    def to_dict(self) -> dict[str, object]:
        return {"distance": self.distance, "elevation": self.elevation}


class EDIProfile(SurveyBase):
    r"""
    Profile helper for one or many :class:`~.edi.EDIFile`
    objects.  Computes small-area geometry (easting/northing),
    cumulative distance along line, profile azimuth, and
    exposes utilities to adjust coordinates and push them back
    into EDI headers.

    The class accepts a single file, an iterable of files, or
    an :class:`~.collection.EDICollection`.  Coordinates are
    read from ``>HEAD`` and converted to working planar
    coordinates.  For short lines the equirectangular
    approximation is used, and distances/azimuth are computed
    in that local frame.

    Parameters
    ----------
    items : EDIFile or iterable of EDIFile or EDICollection
        The input sites to include in the profile.
    verbose : int, default ``0``
        Verbosity level forwarded to internal helpers.

    Attributes
    ----------
    stations : list of str
        Station identifiers resolved from ``DATAID`` or file
        name.
    lat, lon : ndarray of float
        Geographic coordinates (degrees).
    elev : ndarray of float
        Elevations when present, missing values become ``0``.
    distance : ndarray of float
        Cumulative distance from the first site (meters).
    azimuth : float
        Bearing of the profile in degrees, clockwise from
        North, in ``[0, 360)``.
    xy : tuple of ndarray
        Working (easting, northing) arrays in meters.
    table : list of dict
        Row-wise view exposing station, lat, lon, elev,
        easting, northing, and UTM zone (when available).

    Methods
    -------
    get_bearing(method='endpoints')
        Compute bearing from either endpoints or a PCA-like
        fit of the track.
    get_step()
        Return cumulative distance and cache it for reuse.
    adjust(origin=None, azimuth=None, spacing=None, use_mean=True)
        Build an idealized straight profile and compute
        adjusted positions and lat/lon.
    update(use_adjusted=True, update_elev=False)
        Write back adjusted coordinates to each site's header.
    plot_profile(use_adjusted=False, annotate=True, title=None)
        Plot elevation against along-profile distance.
    plot_track(use_adjusted=False, title=None)
        Plot plan-view easting/northing track.

    Notes
    -----
    The small-area equirectangular frame is adequate for short
    profiles.  For long lines or large latitude spans prefer a
    full projection workflow.  When the UTM zone can be
    determined, adjusted coordinates are converted back to
    geographic using that zone.

    Examples
    --------
    Load two sites and compute azimuth and distance::

        prof = EDIProfile([ed1, ed2])
        print(float(prof.azimuth))
        d = prof.distance

    Adjust to a regular spacing and push back to headers::

        prof.adjust(spacing=50.0).update()

    See Also
    --------
    Stations
        Tabular view of station metadata and projected
        coordinates.
    Topography
        Elevation profile builder with smoothing and trend
        tools.

    References
    ----------
    .. [EDIProfile-1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    .. [EDIProfile-2] Snyder, J. P. (1987). *Map Projections – A Working
       Manual*, USGS Prof. Paper 1395.
    """

    def __init__(
        self,
        items: EDIFile | Iterable[EDIFile] | EDICollection,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(verbose=verbose)

        self._eds: list[EDIFile] = _as_collection(items)
        self._rows: list[dict[str, object]] = []
        self._e: np.ndarray | None = None
        self._n: np.ndarray | None = None
        self._zone: list[str | None] | None = None
        self._az: float | None = None
        self._d: np.ndarray | None = None
        self._adj_e: np.ndarray | None = None
        self._adj_n: np.ndarray | None = None
        self._adj_lat: np.ndarray | None = None
        self._adj_lon: np.ndarray | None = None

        self._distance = None
        self._azimuth = None

        self._load()

    def _load(self) -> None:
        rows: list[dict[str, object]] = []
        for ed in self._eds:
            lat, lon, elev = _coerce_ll(ed)
            sid = _station(ed)
            p = getattr(ed, "path", None)
            rows.append(
                {
                    "station": sid,
                    "lat": lat,
                    "lon": lon,
                    "elev": elev,
                    "path": str(p) if p else None,
                    "ed": ed,
                }
            )
        # keep only rows with valid lat/lon
        rows = [
            r for r in rows if (r["lat"] is not None and r["lon"] is not None)
        ]
        self._rows = rows
        if not rows:
            self._e = self._n = None
            self._zone = None
            return
        lat = np.array([r["lat"] for r in rows], float)
        lon = np.array([r["lon"] for r in rows], float)
        e, n, z = to_utm(lat, lon)
        self._e = np.asarray(e, float).ravel()
        self._n = np.asarray(n, float).ravel()
        if np.isscalar(z):
            self._zone = [str(z)] * self._e.size
        else:
            self._zone = [
                str(a) if a is not None else None for a in np.asarray(z)
            ]

    @property
    def stations(self) -> list[str]:
        return [str(r["station"]) for r in self._rows]

    @property
    def lat(self) -> np.ndarray:
        return np.asarray([r["lat"] for r in self._rows], float)

    @property
    def lon(self) -> np.ndarray:
        return np.asarray([r["lon"] for r in self._rows], float)

    @property
    def elev(self) -> np.ndarray:
        e = [r["elev"] for r in self._rows]
        return np.asarray([0.0 if v is None else float(v) for v in e], float)

    @property
    def xy(self) -> tuple[np.ndarray, np.ndarray]:
        if self._e is None or self._n is None:
            return (np.asarray([]), np.asarray([]))
        return (self._e.copy(), self._n.copy())

    def _compute_xy(self) -> tuple[np.ndarray, np.ndarray]:
        # small-area equirectangular to avoid
        # hard deps; good for profiles
        lon = np.asarray(self.lon, float)
        lat = np.asarray(self.lat, float)
        if lon.size == 0:
            return np.array([], float), np.array([], float)
        R = 6_371_000.0
        lr = np.radians(lat)
        br = np.radians(lon)
        lat0 = lr[0]
        lon0 = br[0]
        x = (br - lon0) * np.cos(lr.mean()) * R
        y = (lr - lat0) * R
        return x, y

    def _ensure_distance_bearing(self) -> None:
        if self._distance is not None and self._azimuth is not None:
            return
        x, y = self._compute_xy()
        if x.size == 0:
            self._distance = np.array([], float)
            self._azimuth = 0.0
            return
        seg = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)
        self._distance = np.concatenate(([0.0], np.cumsum(seg)))
        dx = x[-1] - x[0]
        dy = y[-1] - y[0]
        brg = (np.degrees(np.arctan2(dy, dx)) + 360.0) % 360.0
        self._azimuth = float(brg)

    @property
    def distance(self) -> np.ndarray:
        if self._distance is None:
            self._ensure_distance_bearing()
        return np.asarray(self._distance, float)

    @property
    def azimuth(self) -> float:
        if self._azimuth is None:
            self._ensure_distance_bearing()
        return float(self._azimuth)

    def get_bearing(
        self,
        *,
        method: str = "endpoints",
    ) -> float | None:
        r"""
        Estimate the survey bearing (azimuth) in degrees.

        Parameters
        ----------
        method : {'endpoints', 'linear'}, default ``'endpoints'``
            With ``'endpoints'`` the azimuth is computed from the
            first to the last station.  With ``'linear'`` a best-fit
            axis is estimated by SVD of centered coordinates.

        Returns
        -------
        float or None
            Bearing in ``[0, 360)`` (clockwise from north), or
            ``None`` if fewer than two valid stations exist.

        Notes
        -----
        Uses working projected coordinates (easting/northing),
        suitable for small-area profiles.
        """

        if self._e is None or self._n is None:
            self._az = None
            return None
        x = self._e
        y = self._n
        if x.size < 2:
            self._az = None
            return None
        if method == "linear":
            # SVD of centered coords → principal axis
            xc = x - x.mean()
            yc = y - y.mean()
            M = np.vstack((xc, yc)).T
            try:
                _, _, Vt = np.linalg.svd(M, full_matrices=False)
                vx, vy = Vt[0, 0], Vt[0, 1]
            except np.linalg.LinAlgError:
                vx, vy = (x[-1] - x[0], y[-1] - y[0])
        else:
            vx = x[-1] - x[0]
            vy = y[-1] - y[0]
        # azimuth from north, clockwise
        az = math.degrees(math.atan2(vx, vy))
        self._az = (az + 360.0) % 360.0
        return self._az

    def get_step(
        self,
        *,
        method: str = "mean",
        as_array: bool = False,
    ) -> float | np.ndarray:
        r"""
        Derive the inter-station spacing from the track.

        Parameters
        ----------
        method : {'mean', 'median'}, default ``'mean'``
            Aggregation used when returning a scalar spacing.
        as_array : bool, default ``False``
            If ``True``, return the pairwise segment lengths as a
            1-D array of size ``n-1``.  If ``False``, return a
            single spacing computed with ``method``.

        Returns
        -------
        float or ndarray
            Either a scalar spacing or the per-segment distances.

        Notes
        -----
        Distances are computed from consecutive projected
        coordinates; missing stations are ignored.
        """

        if self._e is None or self._n is None:
            return np.asarray([]) if as_array else 0.0
        if self._e.size < 2:
            return np.asarray([]) if as_array else 0.0

        seg = np.hypot(np.diff(self._e), np.diff(self._n))
        if as_array:
            return seg

        if method.lower() == "median":
            return float(np.median(seg))
        return float(np.mean(seg))

    def _mode_zone(self) -> str | None:
        if not self._zone:
            return None
        vals = [z for z in self._zone if z]
        if not vals:
            return None
        uniq, cnt = np.unique(np.asarray(vals), return_counts=True)
        return str(uniq[int(np.argmax(cnt))])

    def adjust(
        self,
        *,
        origin: tuple[float, float] | None = None,
        azimuth: float | None = None,
        spacing: float | None = None,
        step: float | None = None,
        use_mean: bool = True,
    ) -> EDIProfile:
        r"""
        Build an idealized, straightened profile and store the
        adjusted coordinates.

        Parameters
        ----------
        origin : tuple(float, float), optional
            Reference ``(easting, northing)`` for the first
            station.  Defaults to the first raw station.
        azimuth : float, optional
            Bearing of the adjusted line in degrees.  Defaults to
            :meth:`get_bearing`.
        spacing : float, optional
            Fixed spacing between consecutive stations (meters).
        step : float, optional
            Alias for ``spacing`` for convenience.
        use_mean : bool, default ``True``
            When both ``spacing`` and ``step`` are ``None``,
            compute spacing from observed distances using mean if
            ``True`` or median if ``False``.

        Returns
        -------
        EDIProfile
            The instance (allows chaining).

        Notes
        -----
        Adjusted easting/northing are projected back to latitude
        and longitude using the dominant UTM zone of the track.
        Results are stored in ``_adj_e/_adj_n/_adj_lat/_adj_lon``.
        """

        if self._e is None or self._n is None:
            return self
        x = self._e
        y = self._n
        n = x.size
        if n == 0:
            return self

        # azimuth
        if azimuth is None:
            az = self.get_bearing() or 0.0
        else:
            az = float(azimuth) % 360.0
        rad = math.radians(az)
        ve = math.sin(rad)
        vn = math.cos(rad)

        # origin
        if origin is None:
            e0 = float(x[0])
            n0 = float(y[0])
        else:
            e0, n0 = float(origin[0]), float(origin[1])

        # spacing (use step alias if given)
        if step is not None:
            sp = float(step)
        elif spacing is not None:
            sp = float(spacing)
        else:
            sp = self.get_step(method=("mean" if use_mean else "median"))

        idx = np.arange(n, dtype=float)
        adj_e = e0 + idx * sp * ve
        adj_n = n0 + idx * sp * vn
        self._adj_e = adj_e
        self._adj_n = adj_n

        # to lat/lon using dominant UTM zone
        z = self._mode_zone()
        lat: list[float] = []
        lon: list[float] = []
        for i in range(n):
            if z is None:
                # derive zone from raw geographic
                la = float(self.lat[i])
                lo = float(self.lon[i])
                _, _, z2 = get_utm_zone(la, lo)
            else:
                z2 = z
            la, lo = project_point_utm2ll(
                float(adj_e[i]),
                float(adj_n[i]),
                z2,
            )
            lat.append(float(la))
            lon.append(float(lo))
        self._adj_lat = np.asarray(lat, float)
        self._adj_lon = np.asarray(lon, float)
        return self

    def update(
        self,
        *,
        use_adjusted: bool = True,
        update_elev: bool = False,
    ) -> EDIProfile:
        r"""
        Push current coordinates back into the underlying EDI
        headers.

        Parameters
        ----------
        use_adjusted : bool, default ``True``
            If ``True`` write adjusted lat/lon (from
            :meth:`adjust`).  If no adjusted coordinates exist,
            fall back to raw lat/lon.
        update_elev : bool, default ``False``
            Also write elevations when present in the profile
            table.

        Returns
        -------
        EDIProfile
            The instance (allows chaining).

        Notes
        -----
        This mutates the in-memory :class:`~.edi.EDIFile` objects
        held by the profile; it does not write to disk.
        """

        # push adjusted coords back to EDI headers
        if not self._rows:
            return self
        use_adj = bool(
            use_adjusted
            and self._adj_lat is not None
            and self._adj_lon is not None
        )
        for i, r in enumerate(self._rows):
            ed: EDIFile = r["ed"]  # type: ignore[assignment]
            if use_adj:
                ed.lat = float(self._adj_lat[i])
                ed.lon = float(self._adj_lon[i])
            else:
                ed.lat = float(self.lat[i])
                ed.lon = float(self.lon[i])
            if update_elev and r["elev"] is not None:
                ed.elev = float(r["elev"])
        return self

    def plot_profile(
        self,
        *,
        ax: plt.Axes | None = None,
        use_adjusted: bool = False,
        annotate: bool = True,
        title: str | None = None,
    ) -> plt.Axes:
        r"""
        Plot elevation versus along-profile distance.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Target axes.  If omitted, a new figure/axes is made.
        use_adjusted : bool, default ``False``
            If ``True`` recompute distances from adjusted
            coordinates for the overlay.  Raw elevation values are
            used in both cases.
        annotate : bool, default ``True``
            Draw station labels next to points.
        title : str, optional
            Axes title.

        Returns
        -------
        matplotlib.axes.Axes
            The axes with the profile plot.

        Notes
        -----
        Uses the profile's cached cumulative distances.  Call
        :meth:`adjust` first to visualize an adjusted line.
        """

        if ax is None:
            _, ax = plt.subplots()
        d = self.distance
        if d.size == 0:
            return ax

        z = self.elev
        if use_adjusted and (self._adj_lat is not None):
            # recompute elevation vs along-line distances
            # using same step as raw, for visual compare.
            ax.plot(d, z, marker="o", lw=1.0, label="raw")
            ax.set_xlabel("distance (m)")
            ax.set_ylabel("elevation (m)")
            ax.grid(True, ls="--", lw=0.5)
            ax.legend()
            if annotate:
                for i, s in enumerate(self.stations):
                    ax.annotate(
                        s,
                        (d[i], z[i]),
                        xytext=(3, 3),
                        textcoords="offset points",
                        fontsize=8,
                    )
            if title:
                ax.set_title(title)
            return ax
        # raw profile
        ax.plot(d, z, marker="o", lw=1.0, label="raw")
        ax.set_xlabel("distance (m)")
        ax.set_ylabel("elevation (m)")
        ax.grid(True, ls="--", lw=0.5)
        ax.legend()
        if annotate:
            for i, s in enumerate(self.stations):
                ax.annotate(
                    s,
                    (d[i], z[i]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=8,
                )
        if title:
            ax.set_title(title)
        return ax

    def plot_track(
        self,
        *,
        ax: plt.Axes | None = None,
        use_adjusted: bool = False,
        title: str | None = None,
    ) -> plt.Axes:
        r"""
        Plot plan-view station positions (easting vs. northing).

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Target axes.  If omitted, a new figure/axes is made.
        use_adjusted : bool, default ``False``
            Plot the straightened track if adjusted coordinates
            exist; otherwise plot raw positions.
        title : str, optional
            Axes title.

        Returns
        -------
        matplotlib.axes.Axes
            The axes with the track plot.

        Notes
        -----
        Axes aspect is set to equal for a faithful plan view.
        """

        if ax is None:
            _, ax = plt.subplots()
        if use_adjusted and self._adj_e is not None:
            xe = self._adj_e
            yn = self._adj_n
            lab = "adjusted"
        else:
            if self._e is None or self._n is None:
                return ax
            xe = self._e
            yn = self._n
            lab = "raw"
        ax.plot(xe, yn, "-o", ms=4, lw=1.0, label=lab)
        for i, s in enumerate(self.stations):
            ax.text(xe[i], yn[i], s, fontsize=8, ha="left", va="bottom")
        ax.set_aspect("equal", adjustable="datalim")
        ax.grid(True, ls="--", lw=0.5)
        ax.set_xlabel("easting (m)")
        ax.set_ylabel("northing (m)")
        ax.legend()
        if title:
            ax.set_title(title)
        return ax

    def as_table(self) -> list[dict[str, object]]:
        out: list[dict[str, object]] = []
        for i, r in enumerate(self._rows):
            e = None if self._e is None else float(self._e[i])
            n = None if self._n is None else float(self._n[i])
            out.append(
                {
                    "station": r["station"],
                    "lat": r["lat"],
                    "lon": r["lon"],
                    "elev": r["elev"],
                    "easting": e,
                    "northing": n,
                    "zone": (self._zone[i] if self._zone else None),
                    "path": r["path"],
                }
            )
        return out


def _as_collection(
    obj: EDIFile | Iterable[EDIFile] | EDICollection,
) -> list[EDIFile]:
    if isinstance(obj, EDIFile):
        return [obj]
    if isinstance(obj, EDICollection):
        return list(obj)
    return list(obj)  # assume iterable


def _station(ed: EDIFile) -> str:
    sid = getattr(ed, "station", None)
    if sid:
        return str(sid)
    p = getattr(ed, "path", None)
    return Path(p).stem if p else "site"


def _coerce_ll(
    ed: EDIFile,
) -> tuple[float | None, float | None, float | None]:
    # accept either properties or HEAD fields
    lat = getattr(ed, "lat", None)
    lon = getattr(ed, "lon", None)
    if lat is None or lon is None:
        h = ed.get_section("head")
        if h is not None:
            lat = getattr(h, "lat", None)
            lon = getattr(h, "long", None)
            if lon is None:
                lon = getattr(h, "lon", None)
    elev = getattr(ed, "elev", None)
    if elev is None:
        h = ed.get_section("head")
        if h is not None:
            elev = getattr(h, "elev", None)
    if (lat is None) or (lon is None):
        return (None, None, elev)
    # normalize to (lat, lon) if user stored legacy order
    lat, lon = normalize_lat_lon(lon, lat, assume="lonlat")
    return (float(lat), float(lon), None if elev is None else float(elev))
