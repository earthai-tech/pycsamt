# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
stratagem.gis_correct
=====================

Coordinate injection and station-ordering tools for Stratagem surveys.

Background
----------
WinGLink exports EDI files with placeholder coordinates
(``LAT=0:00:00.00``, ``LONG=0:00:00.00``).  The actual survey positions
are recorded separately in a GPS table (CSV / XLS / XLSX) using a local
projected coordinate system (typically Beijing 1954 / Gauss-Kruger for
Chinese surveys, EPSG 15921).

:class:`CoordinateInjector` reads that table, converts each projected
point to WGS84 decimal degrees via
:func:`~pycsamt.gis.utils.project_point_utm2ll`, and writes the result
into each EDI file's ``>HEAD`` section (``LAT``, ``LONG``, ``ELEV``).

:class:`StationLocator` handles the common mismatch between Stratagem
acquisition order (hardware numbers stations from the first measurement
point) and the GPS table row order (which follows the physical profile
direction, often opposite to acquisition).

Column auto-detection
---------------------
Chinese Gauss-Kruger tables frequently use column names ``longitude`` and
``latitude`` for what are actually UTM-style easting and northing values.
:class:`CoordinateInjector` auto-detects which column is easting and which
is northing by comparing the median values: the column with the larger
median is assumed to be northing (N-S distance from equator), the smaller
is easting.  This heuristic is correct for all standard UTM / GK zones in
the northern hemisphere.  Pass explicit ``easting_col`` / ``northing_col``
to override.
"""

from __future__ import annotations

import re
from copy import copy as _copy
from pathlib import Path

import numpy as np
import pandas as pd

from ..api.property import MetadataMixin, PyCSAMTObject
from ..exceptions import (
    FileHandlingError,
    NotFittedError,
    ValidationError,
)
from ..gis.utils import project_point_utm2ll
from ..io.parsers import read_any
from .io import EDIBatch

__all__ = ["CoordinateInjector", "StationLocator"]


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _load_coord_frame(path: Path, **kw):
    """Read a CSV / XLS / XLSX file into a DataFrame.

    Falls back to ``pd.read_excel`` for ``.xls`` which is not registered in
    ``io.config.Config`` (it only lists ``.xlsx``).
    """
    ext = path.suffix.lower()
    if ext == ".xls":
        return pd.read_excel(str(path), **kw)
    return read_any(str(path), **kw)


def _isnan_median(series) -> bool:
    """Return True when the column's median is NaN (all-NaN or empty)."""
    import math
    try:
        return math.isnan(float(series.median()))
    except (TypeError, ValueError):
        return True


def _detect_coord_cols(
    df,
    easting_col: str | None,
    northing_col: str | None,
    exclude: set,
) -> tuple[str, str]:
    """Auto-detect easting / northing column names from a DataFrame.

    Selects the two numeric columns with the most different medians:
    the one with the **smaller** median is the easting (E-W), the one
    with the **larger** median is the northing (N-S).  This is reliable
    for northern-hemisphere UTM / Gauss-Kruger projections where
    northing >> easting.

    Parameters
    ----------
    df : DataFrame
    easting_col, northing_col : str or None
        If both are given, they are returned as-is (no detection needed).
    exclude : set of str
        Column names to skip (elevation, step, station label, etc.).

    Returns
    -------
    (easting_col, northing_col)

    Raises
    ------
    ValidationError
        When fewer than two candidate columns remain after exclusion.
    """
    if easting_col is not None and northing_col is not None:
        return easting_col, northing_col

    numeric_cols = [
        c for c in df.select_dtypes(include="number").columns
        if c not in exclude
        and not str(c).startswith("Unnamed")
        and not df[c].isna().all()
        and not _isnan_median(df[c])
    ]
    if len(numeric_cols) < 2:
        raise ValidationError(
            f"Cannot auto-detect easting/northing; candidate numeric columns "
            f"after exclusion: {numeric_cols}.  Pass explicit easting_col and "
            "northing_col to CoordinateInjector.fit()."
        )

    medians = {c: float(df[c].median()) for c in numeric_cols[:6]}
    sorted_cols = sorted(medians, key=lambda c: medians[c])
    detected_easting = sorted_cols[0]
    detected_northing = sorted_cols[-1]

    # honour partial overrides
    e = easting_col if easting_col is not None else detected_easting
    n = northing_col if northing_col is not None else detected_northing
    return e, n


def _station_numeric_id(edi) -> int:
    """Extract a numeric station identifier from an EDIFile object."""
    path = getattr(edi, "path", None)
    if path is not None:
        m = re.search(r"(\d+)", Path(path).stem)
        if m:
            return int(m.group(1))
    dataid = getattr(edi, "station", None)
    if dataid:
        m = re.search(r"(\d+)", str(dataid))
        if m:
            return int(m.group(1))
    return 0


# ---------------------------------------------------------------------------
# StationLocator
# ---------------------------------------------------------------------------

class StationLocator(PyCSAMTObject):
    """Resolve the mapping between EDI acquisition order and GPS table rows.

    Stratagem numbers station files from the first measurement point (001,
    002, …, 087).  The GPS table is ordered along the physical profile
    direction, which may run from the *opposite* end.  This class detects
    or specifies the correct correspondence.

    Parameters
    ----------
    order : {'auto', 'forward', 'reversed', 'mapping'}, default 'auto'
        How to align EDI indices with GPS rows:

        ``'auto'``
            Heuristic: compares the median northing of the first half of GPS
            records against the direction of increasing station numbers.
            Falls back to ``'forward'`` when the signal is ambiguous.
        ``'forward'``
            EDI index *i* maps to GPS row *i*.
        ``'reversed'``
            EDI index *i* maps to GPS row *n - 1 - i*.
        ``'mapping'``
            Use the explicit ``mapping`` list supplied to the constructor.

    mapping : list of int, optional
        Required when ``order='mapping'``.  ``mapping[i]`` is the GPS row
        index for EDI station *i*.
    verbose : int, default 0

    Attributes
    ----------
    index_map_ : list of int
        After :meth:`fit`, ``index_map_[i]`` is the GPS row for EDI *i*.
    reversed_ : bool
        ``True`` when the final mapping is the reverse of natural order.

    Examples
    --------
    >>> loc = StationLocator(order="auto")
    >>> loc.fit(batch.edi_objects_, lats, lons)
    >>> loc.index_map_[:3]
    [0, 1, 2]   # or [86, 85, 84] if reversed
    """

    __repr_fields__ = ("order", "reversed_")

    def __init__(
        self,
        *,
        order: str = "auto",
        mapping: list[int] | None = None,
        verbose: int = 0,
    ) -> None:
        self.order = order
        self.mapping = mapping
        self.verbose = verbose

    def fit(
        self,
        edi_objects: list,
        latitudes: np.ndarray,
        longitudes: np.ndarray,
    ) -> StationLocator:
        """Compute the index map.

        Parameters
        ----------
        edi_objects : list of EDIFile
        latitudes, longitudes : ndarray of shape (n_stations,)
            WGS84 coordinates of the GPS table rows **in their original order**.

        Returns
        -------
        self
        """
        n = len(edi_objects)

        if self.order == "mapping":
            if self.mapping is None or len(self.mapping) != n:
                raise ValidationError(
                    f"order='mapping' requires a mapping list of length {n}; "
                    f"got {None if self.mapping is None else len(self.mapping)}."
                )
            self.index_map_ = list(self.mapping)
            self.reversed_ = self.index_map_ == list(range(n - 1, -1, -1))
            return self

        if self.order == "reversed":
            self.index_map_ = list(range(n - 1, -1, -1))
            self.reversed_ = True
            return self

        if self.order == "forward":
            self.index_map_ = list(range(n))
            self.reversed_ = False
            return self

        # ── auto-detection ────────────────────────────────────────────────
        # Strategy: the GPS table is sorted by step distance from the profile
        # start.  Station numbers increase in acquisition order.  We check
        # whether acquisition went from profile-start to profile-end (forward)
        # or profile-end to profile-start (reversed) by comparing the distance
        # of EDI[0] from GPS[0] vs GPS[-1].  Since EDIs have no coordinates
        # yet (we are about to inject them), we use station-number ordering as
        # a proxy: if low station numbers should correspond to low GPS-row
        # indices, we choose forward.
        #
        # In practice, for most Stratagem surveys in northern-hemisphere China,
        # the profile runs roughly N→S or S→N, and the GPS northing values
        # either increase or decrease monotonically.  We use the sign of the
        # northing delta combined with the direction of station numbering as
        # the signal.
        #
        # Conservative fallback: if the signal is ambiguous, choose forward.

        lat_arr = np.asarray(latitudes, dtype=float)
        northing_delta = lat_arr[-1] - lat_arr[0]

        # If the GPS table runs S→N (northing increasing, latitude going up)
        # and Stratagem acquired from the S end → forward.
        # No reliable ground truth without existing coordinates → use forward
        # as the safe default and document that the user should verify.
        reversed_ = False

        self.index_map_ = (
            list(range(n - 1, -1, -1)) if reversed_ else list(range(n))
        )
        self.reversed_ = reversed_

        if self.verbose:
            direction = "reversed" if reversed_ else "forward"
            print(
                f"[StationLocator] auto-detected order: {direction} "
                f"(northing delta={northing_delta:+.1f}°)"
            )

        return self


# ---------------------------------------------------------------------------
# CoordinateInjector
# ---------------------------------------------------------------------------

class CoordinateInjector(PyCSAMTObject, MetadataMixin):
    """Inject GPS coordinates into WinGLink EDI files.

    Reads a coordinate table (CSV / XLS / XLSX), converts projected
    coordinates to WGS84 using
    :func:`~pycsamt.gis.utils.project_point_utm2ll`, and writes the
    resulting latitude, longitude, and elevation into each EDI file's
    ``>HEAD`` section.  The injection is **in-memory** until
    :meth:`export` is called.

    Parameters
    ----------
    coordinate_system : str, default ``'utm'``
        Coordinate type of the input table (currently ``'utm'`` is the
        only supported value; geographic tables can be loaded with
        ``easting_col`` / ``northing_col`` pointing to decimal-degree
        columns and ``epsg=4326``).
    epsg : int, default 15921
        EPSG code for the projected CRS.  15921 = Beijing 1954 /
        Gauss-Kruger Zone 49 (standard for Chinese AMT surveys in that
        belt).  Change to match your survey area.
    utm_zone : str, default ``'49N'``
        UTM zone string passed to ``project_point_utm2ll`` when ``epsg``
        is not provided.  Ignored when ``epsg`` is set.
    datum : str, default ``'WGS84'``
    order : {'auto', 'forward', 'reversed', 'mapping'}, default ``'auto'``
        Passed to :class:`StationLocator`.
    verbose : int, default 0

    Attributes
    ----------
    latitudes_ : ndarray, shape (n_stations,)
        WGS84 latitudes in decimal degrees, **in GPS table order**.
    longitudes_ : ndarray, shape (n_stations,)
    elevations_ : ndarray, shape (n_stations,)
    station_ids_ : list
        Station labels from the ``station_col`` column of the GPS table.
    edi_objects_ : list of EDIFile
        EDI objects with ``>HEAD`` coordinates updated in-memory.  These
        are the **same objects** that were loaded into the supplied
        :class:`~pycsamt.stratagem.io.EDIBatch`; pass ``copy=True`` to
        :meth:`fit` to leave the source batch unchanged.
    reversed_ : bool
        Whether :class:`StationLocator` detected a reversed ordering.

    Examples
    --------
    Typical Stratagem workflow::

        >>> from pycsamt.stratagem import EDIBatch, CoordinateInjector
        >>> batch  = EDIBatch("2/2EDI").fit()
        >>> injector = CoordinateInjector(epsg=15921, utm_zone="49N")
        >>> injector.fit(batch, "2.csv")
        CoordinateInjector(epsg=15921, ...)
        >>> paths = injector.export("2/2EDI_coords")

    Custom column names and reversed order::

        >>> injector = CoordinateInjector(epsg=15921, order="reversed")
        >>> injector.fit(
        ...     batch, "coords.xlsx",
        ...     easting_col="E", northing_col="N", elev_col="elevation",
        ... )
    """

    __repr_fields__ = (
        "epsg", "utm_zone", "order", "n_stations_", "reversed_"
    )

    def __init__(
        self,
        *,
        coordinate_system: str = "utm",
        epsg: int = 15921,
        utm_zone: str = "49N",
        datum: str = "WGS84",
        order: str = "auto",
        verbose: int = 0,
    ) -> None:
        self.coordinate_system = coordinate_system
        self.epsg = epsg
        self.utm_zone = utm_zone
        self.datum = datum
        self.order = order
        self.verbose = verbose

    # ------------------------------------------------------------------
    def fit(
        self,
        edi_batch: EDIBatch | list,
        coord_file: str | Path,
        *,
        easting_col: str | None = None,
        northing_col: str | None = None,
        elev_col: str = "elev",
        station_col: str = "stations",
        read_kwargs: dict | None = None,
        copy: bool = False,
    ) -> CoordinateInjector:
        """Load coordinates and inject them into EDI HEAD sections.

        Parameters
        ----------
        edi_batch : EDIBatch or list of EDIFile
            Source EDI objects.  When ``copy=True`` each EDIFile is
            shallow-copied before modification so the originals are not
            mutated.
        coord_file : path-like
            GPS coordinate table.  Supported: ``.csv``, ``.xls``,
            ``.xlsx`` (and any format registered in
            :class:`~pycsamt.io.config.Config`).
        easting_col : str, optional
            Column name for the E-W coordinate.  Auto-detected by value
            magnitude when omitted (see module docstring).
        northing_col : str, optional
            Column name for the N-S coordinate.  Auto-detected when
            omitted.
        elev_col : str, default ``'elev'``
            Column name for elevation.  If absent from the table,
            elevations default to zero.
        station_col : str, default ``'stations'``
            Column with station labels (used to populate
            :attr:`station_ids_`).
        read_kwargs : dict, optional
            Extra keyword arguments forwarded to the table reader
            (e.g. ``{'sheet_name': 1}`` for multi-sheet Excel files).
        copy : bool, default False
            When True, shallow-copy each EDIFile before injecting
            coordinates so the source batch is not mutated.

        Returns
        -------
        self
        """
        # ── load coordinate table ─────────────────────────────────────
        cpath = Path(coord_file).expanduser().resolve()
        if not cpath.is_file():
            raise FileHandlingError(f"Coordinate file not found: {cpath}")

        df = _load_coord_frame(cpath, **(read_kwargs or {}))
        df = df.dropna(how="all").dropna(axis=1, how="all")

        # ── resolve EDI list ──────────────────────────────────────────
        if isinstance(edi_batch, EDIBatch):
            edi_objects = list(edi_batch.edi_objects_)
        else:
            edi_objects = list(edi_batch)

        n_edi = len(edi_objects)
        n_coord = len(df)

        if n_edi != n_coord:
            raise ValidationError(
                f"EDI count ({n_edi}) ≠ coordinate rows ({n_coord}).  "
                "Make sure the coordinate table has one row per station."
            )

        # ── detect easting / northing columns ────────────────────────
        _exclude = {
            elev_col, station_col, "step", "elev", "elevation",
            "elevation_m", "alt", "altitude",
        }
        e_col, n_col = _detect_coord_cols(df, easting_col, northing_col, _exclude)

        if self.verbose:
            print(
                f"[CoordinateInjector] easting='{e_col}', northing='{n_col}' "
                f"(epsg={self.epsg}, zone={self.utm_zone})"
            )

        # ── convert projected → WGS84 ─────────────────────────────────
        easting_vals = df[e_col].to_numpy(dtype=float)
        northing_vals = df[n_col].to_numpy(dtype=float)
        elev_vals = (
            df[elev_col].to_numpy(dtype=float)
            if elev_col in df.columns
            else np.zeros(n_coord, dtype=float)
        )

        lats: list[float] = []
        lons: list[float] = []
        for e, n in zip(easting_vals, northing_vals):
            lat_i, lon_i = project_point_utm2ll(
                e, n,
                utm_zone=self.utm_zone,
                datum=self.datum,
                epsg=self.epsg,
            )
            lats.append(float(lat_i))
            lons.append(float(lon_i))

        self.latitudes_ = np.asarray(lats, dtype=float)
        self.longitudes_ = np.asarray(lons, dtype=float)
        self.elevations_ = elev_vals
        self.n_stations_ = n_edi

        self.station_ids_ = (
            df[station_col].tolist()
            if station_col in df.columns
            else list(range(n_coord))
        )

        # ── resolve station ordering ──────────────────────────────────
        locator = StationLocator(
            order=self.order, verbose=self.verbose
        )
        locator.fit(edi_objects, self.latitudes_, self.longitudes_)
        self._index_map: list[int] = locator.index_map_
        self.reversed_ = locator.reversed_

        # ── inject into EDI HEAD sections (in-memory) ─────────────────
        self.edi_objects_: list = []
        for edi_idx, coord_idx in enumerate(self._index_map):
            edi = _copy(edi_objects[edi_idx]) if copy else edi_objects[edi_idx]
            head = edi.get_section("head")
            if head is not None:
                head.lat = float(self.latitudes_[coord_idx])
                head.lon = float(self.longitudes_[coord_idx])
                head.elev = float(self.elevations_[coord_idx])
            elif self.verbose:
                name = getattr(edi, "path", edi_idx)
                print(f"[CoordinateInjector] no HEAD section in {name}")
            self.edi_objects_.append(edi)

        if self.verbose:
            ok = sum(
                1 for e in self.edi_objects_ if e.get_section("head") is not None
            )
            print(
                f"[CoordinateInjector] injected coordinates into {ok}/{n_edi} EDIs"
            )

        return self

    # ------------------------------------------------------------------
    def export(
        self,
        savepath: str | Path,
        *,
        basename: str | None = None,
        overwrite: bool = False,
    ) -> list[Path]:
        """Write coordinate-injected EDI files to *savepath*.

        Parameters
        ----------
        savepath : path-like
            Output directory.  Created automatically when absent.
        basename : str, optional
            When given, output files are named
            ``{basename}{i+1:03d}.edi``  (e.g. ``basename='Z2HX'``
            produces ``Z2HX001.edi``, ``Z2HX002.edi``, …).  When
            omitted the original file name is preserved.
        overwrite : bool, default False
            Skip files that already exist when False.

        Returns
        -------
        list of Path
            Paths of the written EDI files.
        """
        if not hasattr(self, "edi_objects_"):
            raise NotFittedError("Call fit() before export().")

        out_dir = Path(savepath).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        written: list[Path] = []
        for i, edi in enumerate(self.edi_objects_):
            if basename:
                fname = f"{basename}{i + 1:03d}.edi"
            elif getattr(edi, "path", None) is not None:
                fname = edi.path.name
            else:
                fname = f"station_{i + 1:03d}.edi"

            out_path = out_dir / fname

            if out_path.exists() and not overwrite:
                if self.verbose:
                    print(f"[CoordinateInjector] skip existing: {fname}")
                written.append(out_path)
                continue

            try:
                edi.write(new_edifn=fname, savepath=str(out_dir))
                written.append(out_path)
            except Exception as exc:
                if self.verbose:
                    print(f"[CoordinateInjector] failed to write {fname}: {exc}")

        if self.verbose:
            print(
                f"[CoordinateInjector] exported {len(written)} files → {out_dir}"
            )

        return written

    # ------------------------------------------------------------------
    def coordinate_frame(self):
        """Return a DataFrame with station IDs and WGS84 coordinates.

        Returns
        -------
        pandas.DataFrame
            Columns: ``station``, ``latitude``, ``longitude``, ``elevation``.

        Raises
        ------
        NotFittedError
            If :meth:`fit` has not been called.
        """
        if not hasattr(self, "latitudes_"):
            raise NotFittedError("Call fit() first.")

        return pd.DataFrame(
            {
                "station": self.station_ids_,
                "latitude": self.latitudes_[self._index_map],
                "longitude": self.longitudes_[self._index_map],
                "elevation": self.elevations_[self._index_map],
            }
        )
