# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import (
    Any,
    Callable,
    Union,
)

import numpy as np

from ..core.config import StationNamePolicy, get_config
from ..exceptions import SiteError
from ..log.logger import get_logger
from ..utils.handlers import columns_manager
from .cbase import CBBase, CoreParser, ParseMixin
from .edi import EDIFile

logger = get_logger(__name__)

__all__ = ["CollectionMixin", "EDICollection"]

Pathish = Union[str, Path]


class CollectionMixin(ParseMixin):
    r"""
    Utilities for working with collections of :class:`EDIFile`
    objects.

    The mixin adds non-IO helpers that operate on a set of
    items already loaded in memory.  Typical operations
    include filtering, mapping, grouping, and summarizing.
    It is designed to be combined with a small container
    base, such as :class:`CBBase`.

    Parameters
    ----------
    None
        This is a mixin.  It does not define its own public
        constructor parameters.

    Attributes
    ----------
    items : dict
        Mapping from collection key (for example the station
        id) to :class:`EDIFile`.  Expected to be provided by
        the host collection.
    order : list of str
        Insertion order of keys.  Used to preserve stable
        iteration and deterministic summaries.

    Notes
    -----
    The mixin focuses on pure-Python transformations.  The
    typical methods provided by implementations include:

    * ``select(keys=None, predicate=None)`` : return a new
      view with items that match the keys or satisfy the
      boolean predicate.
    * ``map(func)`` : apply a callable to each
      :class:`EDIFile` and collect the results in a list.
    * ``groupby(keyfunc)`` : bucket items by a user-defined
      key and return a mapping to lists of :class:`EDIFile`.
    * ``summary(fields=('station', 'n_freq', ...))`` :
      build a compact table of common attributes for quick
      inspection.

    These helpers do not read from disk and do not mutate
    :class:`EDIFile` content unless explicitly documented by
    the host class.

    Examples
    --------
    Filter by station prefix and compute sample counts::

        view = coll.select(
            predicate=lambda ed: ed.station.startswith('A')
        )
        counts = view.map(lambda ed: ed.Z.n_freq)

    Group by presence of tipper and list stations::

        groups = coll.groupby(
            lambda ed: ed.Tip.tipper is not None
        )
        stations_with_tip = [ed.station for ed in
                             groups.get(True, [])]

    See Also
    --------
    EDICollection
        Concrete collection that mixes this class with a
        container and a loader.
    CBBase
        Minimal container used by collections.
    CoreParser
        Parser used to build :class:`EDIFile` objects.
    EDIFile
        Single-file reader and writer.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def add_from(
        self,
        sources: Pathish | Sequence[Pathish],
        *,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int | None = None,
    ) -> CollectionMixin:
        v = self.verbose if verbose is None else int(verbose)
        pr = CoreParser(
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=v,
        )
        items = pr.parse(sources)
        for ed in items:
            self.add(ed)  # type: ignore[attr-defined]
        errs = pr.errors()
        if errs:
            logger.info("Added with %d errors.", len(errs))
        return self

    def select(
        self, stations: Sequence[str]
    ) -> EDICollection:
        stations = columns_manager( stations, empty_as_none=False )

        keep = {str(s) for s in stations}
        out = EDICollection(verbose=self.verbose)
        for ed in self:  # type: ignore[operator]
            sid = getattr(ed, "station", None) or "-"
            if sid in keep:
                out.add(ed)
        return out

    def where(
        self, fn: Callable[[EDIFile], bool]
    ) -> EDICollection:
        out = EDICollection(verbose=self.verbose)
        for ed in self:  # type: ignore[operator]
            if fn(ed):
                out.add(ed)
        return out

    def sort(
        self,
        key: str | Callable[[EDIFile], object] = "station",
        *,
        reverse: bool = False,
    ) -> EDICollection:
        def _key(ed: EDIFile) -> object:
            if callable(key):
                return key(ed)
            if key == "station":
                return getattr(ed, "station", "")
            if key == "n_freq":
                return int(getattr(ed.Z, "n_freq", 0) or 0)
            if key == "path":
                return str(ed.path) if ed.path else ""
            return getattr(ed, key, None)

        out = EDICollection(
            items=sorted(
                list(self),  # type: ignore[arg-type]
                key=_key,
                reverse=reverse,
            ),
            verbose=self.verbose,
        )
        return out

    @property
    def paths(self) -> list[str]:
        p: list[str] = []
        for ed in self:  # type: ignore[operator]
            p.append(str(ed.path) if ed.path else "-")
        return p

    def nf_stats(self) -> dict:
        vals: list[int] = []
        for ed in self:  # type: ignore[operator]
            vals.append(int(getattr(ed.Z, "n_freq", 0) or 0))
        a = np.asarray(vals, int)
        if a.size == 0:
            return {"min": 0, "max": 0, "mean": 0.0}
        return {
            "min": int(a.min()),
            "max": int(a.max()),
            "mean": float(a.mean()),
        }


class EDICollection(CBBase, CollectionMixin):
    r"""
    High-level collection of :class:`EDIFile` objects with
    integrated loading and summarization.

    The class couples the small container
    :class:`CBBase` with :class:`CoreParser` to discover,
    read, and manage many EDI files at once.  Items are
    indexed by station id or by absolute path, with a
    configurable duplicate policy.

    Parameters
    ----------
    sources : path-like or sequence of path-like, optional
        Files, folders, or glob patterns to load.  When
        omitted the collection is created empty.  Use
        :meth:`load` later to populate it.
    recursive : bool, default ``True``
        Recurse into subdirectories when scanning folders.
    on_dup : {'replace', 'keep'}, default ``'replace'``
        Duplicate policy for items that share the same key.
        With ``'replace'`` the last item wins.  With
        ``'keep'`` the first item is preserved.
    index_by : {'station', 'path'}, default ``'station'``
        Key used to index items.  When ``'station'`` the
        ``EDIFile.station`` is preferred with a fallback to
        file path if unspecified.
    strict : bool, default ``False``
        If ``True`` read errors are raised.  If ``False``
        they are captured and available via
        :meth:`errors`.
    verbose : int, default ``0``
        Verbosity forwarded to the underlying
        :class:`EDIFile` readers.

    Attributes
    ----------
    parser : CoreParser
        Internal parser used by :meth:`load`.  Exposed for
        advanced scenarios or diagnostics.
    items : dict
        Mapping from key to :class:`EDIFile`.
    order : list of str
        Insertion order of collection keys.
    meta : dict
        Free-form metadata that client code can populate.

    Methods
    -------
    load(sources)
        Discover and read items.  Merges into the current
        collection while respecting the duplicate policy.
    select(keys=None, predicate=None)
        Return a filtered view.  Provided by
        :class:`CollectionMixin`.
    map(func)
        Apply a callable to each item and collect results.
    groupby(keyfunc)
        Group items by a user function, returning a dict of
        lists.
    summary(fields=('station', 'n_freq', 'has_tipper'))
        Build a compact per-item table suitable for logging
        or quick inspection.
    errors()
        Return a list of ``(path, exception)`` pairs from the
        last :meth:`load` call.

    Notes
    -----
    The class aims to be pragmatic.  Loading is done with a
    single pass of :class:`CoreParser`.  Summaries are built
    from lightweight attributes that are cheap to access,
    such as ``station``, number of frequencies, and section
    presence.  Heavy analysis should be implemented outside
    the collection to keep the API small.

    Examples
    --------
    Load a folder, keep first duplicates, and summarize::

        coll = EDICollection(
            sources=['data/edi'],
            on_dup='keep', recursive=True
        )
        tbl = coll.summary(
            fields=('station', 'n_freq', 'tipper')
        )
        for row in tbl:
            print(row)

    Add more files later and filter by predicate::

        coll.load(['more/*.edi'])
        mt_only = coll.select(
            predicate=lambda ed: ed.has_section('mtsect')
        )
        stations = mt_only.map(lambda ed: ed.station)

    See Also
    --------
    CoreParser
        Multi-source reader used by :meth:`load`.
    CBBase
        Base container that stores items and preserves order.
    ParseMixin
        Discovery utilities employed by the parser.
    EDIFile
        Single EDI reader and writer stored in the
        collection.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def __init__(
        self,
        items: Iterable[EDIFile] | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        super().__init__(items=items, verbose=verbose)


    @classmethod
    def from_sources(
        cls,
        sources: Pathish | Sequence[Pathish],
        *,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> EDICollection:
        pr = CoreParser(
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=verbose,
        )
        items = pr.parse(sources)
        col = cls(items=items, verbose=verbose)
        errs = pr.errors()
        if errs:
            logger.info(
                "Loaded %d items with %d errors.",
                len(items),
                len(errs),
            )
        return col


    def merge(
        self,
        other: EDICollection,
        *,
        on_dup: str = "replace",
    ) -> EDICollection:
        out = EDICollection(
            items=list(self),
            verbose=self.verbose,
        )
        if on_dup not in {"replace", "keep"}:
            raise ValueError("on_dup must be keep|replace")
        for ed in other:
            sid = getattr(ed, "station", None)
            sid = sid or (str(ed.path) if ed.path else None)
            if sid in out._index and on_dup == "keep":
                continue
            out.add(ed)
        return out


    def _resolve(self, site: str):
        """Find item by key, station, stem or filename
        (case-insensitive)."""
        want = str(site).upper()
        cands = set(_site_candidates(want))

        # fast path: index keys (case-insensitive)
        for key, idx in self._index.items():
            if str(key).upper() in cands:
                return self._items[idx]

        # scan stations
        for ed in self:
            sid = getattr(ed, "station", None)
            if sid and str(sid).upper() in cands:
                return ed

        # scan paths
        low = want.lower()
        for ed in self:
            p = getattr(ed, "path", None)
            if p is None:
                continue
            if (
                p.stem.lower() in {low} |
                {c.lower() for c in cands}
                or p.name.lower() == low
                or str(p).lower() == low
            ):
                return ed

        all_sites = ", ".join(self.stations())
        raise SiteError(
            f"site not found: {site!r}. Try one of: "
            f"{all_sites}\nHint: call "
            "`<collection>.stations()` to list."
        )

    def _head(self, ed: EDIFile):
        return ed.get_section("head")


    def get(
        self,
        site: str,
        what: str,
        default: object | None = None,
    ) -> object | None:
        """
        Quick extractor for a single field/array.
        'what' supports:
          - 'freq', 'z', 'zxx','zxy','zyx','zyy'
          - 'tip','tx','ty'
          - 'station','dataid','lat','lon','elev'
          - 'path','filename'
          - 'spectra','timeseries'
        """
        try:
            ed = self._resolve(site)
        except KeyError:
            return default

        w = str(what).lower()

        # transfer functions
        if w == "freq":
            return getattr(ed.Z, "freq", None)
        if w == "z":
            return getattr(ed.Z, "z", None)
        if w in {"zxx", "zxy", "zyx", "zyy"}:
            m = {"zxx": (0, 0), "zxy": (0, 1),
                 "zyx": (1, 0), "zyy": (1, 1)}[w]
            z = getattr(ed.Z, "z", None)
            if z is None:
                return default
            return z[:, m[0], m[1]]

        if w in {"tip", "tx", "ty"}:
            t = getattr(ed.Tip, "tipper", None)
            if t is None or getattr(t, "size", 0) == 0:
                return default
            if w == "tip":
                return t
            return t[:, 0, 0] if w == "tx" else t[:, 0, 1]

        # headers / identity
        head = self._head(ed)
        if w == "station":
            return getattr(ed, "station", None)
        if w == "dataid":
            return getattr(head, "dataid", None) if head else None
        if w in {"lat", "lon", "elev"}:
            if not head:
                return default
            if w == "lat":
                return getattr(head, "lat", None)
            if w == "lon":
                # support both LONG and LON
                v = getattr(head, "long", None)
                return v if v is not None else getattr(head, "lon", None)
            return getattr(head, "elev", None)

        # path-ish
        if w in {"path", "filename"}:
            p = getattr(ed, "path", None)
            if p is None:
                return default
            return str(p) if w == "path" else p.name

        # higher-level blocks
        if w == "spectra":
            return ed.get_section("spectra")
        if w in {"ts", "timeseries"}:
            return ed.get_section("timeseries")

        return default

    def set(
        self,
        site: str,
        *,
        edi: EDIFile | None = None,
        update: dict[str, object] | None = None,
    ) -> EDIFile:
        """
        Replace or mutate a site's EDI.

        - If 'edi' is given, replace the stored object.
        - Else apply 'update' fields on the existing object:
          keys: 'station','dataid','lat','lon','elev','z',
                'freq','tip','spectra','timeseries'
        """
        ed = self._resolve(site)

        # replace whole object
        if edi is not None:
            # try a collection-level replace
            try:
                self.remove(site)          # optional on CBBase
            except Exception:
                pass
            try:
                self.add(edi)              # CBBase.add
            except Exception:
                # fallback: mutate mapping directly
                idx = getattr(self, "_index", None)
                if isinstance(idx, dict):
                    idx[site] = edi
            return edi

        # mutate in place
        upd = update or {}
        head = self._head(ed)

        if "station" in upd:
            try:
                ed.station = str(upd["station"])
            except Exception:
                if head:
                    head.dataid = str(upd["station"])

        if "dataid" in upd and head:
            head.dataid = str(upd["dataid"])

        if "lat" in upd and head:
            try:
                head.lat = float(upd["lat"])  # type: ignore[arg-type]
            except Exception:
                pass

        if "lon" in upd and head:
            try:
                # support LONG or LON
                if hasattr(head, "long"):
                    head.long = float(upd["lon"])  # type: ignore[arg-type]
                else:
                    head.lon = float(upd["lon"])  # type: ignore[attr-defined]
            except Exception:
                pass

        if "elev" in upd and head:
            try:
                head.elev = float(upd["elev"])  # type: ignore[arg-type]
            except Exception:
                pass

        if "freq" in upd and hasattr(ed.Z, "freq"):
            ed.Z._freq = np.asarray(upd["freq"], float)

        if "z" in upd and hasattr(ed.Z, "z"):
            ed.Z._z = np.asarray(upd["z"], complex)

        if "tip" in upd and hasattr(ed.Tip, "tipper"):
            ed.Tip._tipper = np.asarray(upd["tip"], complex)

        if "spectra" in upd:
            ed.add_section("spectra", upd["spectra"])

        if "timeseries" in upd or "ts" in upd:
            val = upd.get("timeseries", upd.get("ts"))
            ed.add_section("timeseries", val)

        return ed


    def adjust(
        self,
        site: str,
        *,
        dlat: float | None = None,
        dlon: float | None = None,
        lat: float | None = None,
        lon: float | None = None,
        elev: float | None = None,
        rename: str | None = None,
    ) -> EDIFile:
        """
        Shift or set position and/or rename station.
        """
        ed = self._resolve(site)
        head = self._head(ed)
        if head is None:
            raise ValueError("no >HEAD to adjust")

        # rename first so index logic can react later if needed
        if rename:
            try:
                ed.station = str(rename)
            except Exception:
                head.dataid = str(rename)

        if lat is not None:
            try:
                head.lat = float(lat)
            except Exception:
                pass
        if lon is not None:
            try:
                if hasattr(head, "long"):
                    head.long = float(lon)
                else:
                    head.lon = float(lon)   # type: ignore[attr-defined]
            except Exception:
                pass
        if elev is not None:
            try:
                head.elev = float(elev)
            except Exception:
                pass

        # relative shifts
        if dlat is not None:
            try:
                head.lat = float(head.lat) + float(dlat)  # type: ignore[arg-type]
            except Exception:
                pass
        if dlon is not None:
            try:
                if hasattr(head, "long"):
                    head.long = float(head.long) + float(dlon)
                else:
                    head.lon = float(head.lon) + float(dlon)  # type: ignore[attr-defined]
            except Exception:
                pass

        return ed

    def export(
        self,
        output_dir: Pathish,
        *,
        file_pattern: str = "{station}.edi",
        **edi_write_kwargs,
    ) -> dict:
        """
        Exports all EDIFile items to a directory with advanced options.

        This method orchestrates the writing of each EDIFile in the
        collection to a specified directory, with flexible naming and
        clear error reporting.

        Parameters
        ----------
        output_dir : Pathish
            The path to the directory where EDI files will be saved.
            It will be created if it does not exist.
        file_pattern : str, default="{station}.edi"
            A format string for output filenames. Can use attributes
            like `{station}`.
        **edi_write_kwargs
            Keyword arguments to be passed directly to each
            `EDIFile.write()` call (e.g., `datatype="mt"`,
            `synthesize_spectra=True`).

        Returns
        -------
        dict
            A dictionary with two keys:
            - 'successful': A list of paths to successfully written files.
            - 'failed': A list of tuples, where each tuple contains the
              station name and the error that occurred during writing.
        """
        out_dir = Path(str(output_dir)).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        successful_paths = []
        failed_items = []

        items_iterator = self._items

        # Optional: Use tqdm for a progress bar if available
        try:
            from tqdm import tqdm
        except ImportError:
            tqdm = None

        if tqdm:
            items_iterator = tqdm(self._items, desc="Exporting EDI Files")

        for ed in items_iterator:
            sid = self._site_of(ed) or "unknown_station"
            try:
                filename = file_pattern.format(station=sid)
                output_path = out_dir / filename

                # Delegate the actual writing to the EDIFile instance
                written_path = ed.write(
                    new_edifn=str(output_path), **edi_write_kwargs
                )
                successful_paths.append(written_path)
            except Exception as e:
                failed_items.append((sid, e))
                logger.error(f"Failed to write EDI for station {sid}: {e}")

        return {"successful": successful_paths, "failed": failed_items}

    def fetch(
        self,
        site: str | None = None,
        lat: float | None = None,
        lon: float | None = None,
        tol: float = 0.001,
        first: bool = False,
        strict: bool = False,
        policy: StationNamePolicy | None = None,
        **kwargs,
    ) -> EDIFile | list[EDIFile] | None:

        r"""
        Fetch EDI files by station, coordinates, or other metadata.

        This method searches the collection for
        :class:`~pycsamt.seg.edi.EDIFile` objects by
        station name, geographic coordinates, or any
        attribute exposed on the EDI or its ``>HEAD``
        section.

        Parameters
        ----------
        site : str or int or float, optional
            Station selector. Case-insensitive. Numeric-like
            values are normalized using the active
            :class:`~pycsamt.core.config.StationNamePolicy`
            and station-variant expansion. Variants tried
            include the raw token, integer form, zero-padded
            prefix form (e.g. ``S150``), and common ``×10``
            slips (``1500`` <-> ``150``).
        lat : float, optional
            Latitude in decimal degrees. Compared against the
            ``>HEAD`` field ``lat`` using the tolerance
            ``tol``.
        lon : float, optional
            Longitude in decimal degrees. Compared against the
            ``>HEAD`` field ``long`` using the tolerance
            ``tol``.
        tol : float, default 0.001
            Absolute tolerance (decimal degrees) used when
            matching ``lat`` and ``lon``.
        first : bool, default False
            If ``True``, return only the first matching object
            or ``None`` if no match (unless ``strict=True``,
            see *Raises*). If ``False``, return a list of all
            matches (possibly empty).
        strict : bool, default False
            If ``True`` and no match is found, raise
            :class:`~pycsamt.seg.exceptions.SiteError` with a
            helpful message suggesting
            ``<collection>.stations()``.
        policy : StationNamePolicy, optional
            Custom naming policy. Defaults to
            :func:`~pycsamt.core.config.get_config`
            ``.station_policy``.
        **kwargs : Any
            Additional attribute filters. Keys are matched
            (case-insensitive) first on the EDI object, then
            on its ``head`` section. String values compare in
            a case-insensitive way; non-strings compare by
            equality.

        Returns
        -------
        EDIFile or list of EDIFile or None
            If ``first=True``, an :class:`EDIFile` or ``None``
            when no match and ``strict=False``. If
            ``first=False``, a list of matches (possibly
            empty).

        Raises
        ------
        SiteError
            Raised when no match is found and ``strict=True``.
            The error message lists available stations and
            suggests calling ``<collection>.stations()``.

        Notes
        -----
        - Station matching expands candidate variants using the
          active :class:`StationNamePolicy` and a small set of
          robust transforms (integer form, zero-padded with
          prefix, and common ``×10`` slips).
        - Geographic matching uses absolute tolerance on decimal
          degrees and compares against ``>HEAD`` fields ``lat``
          and ``long``.
        - Arbitrary filters in ``**kwargs`` are applied after
          station and coordinate tests. String comparisons are
          case-insensitive; other types use equality.

        Examples
        --------
        >>> # List available stations
        >>> coll.stations()
        ['S150', 'S200', 'S250', 'S300', ...]

        >>> # Fetch a single station by numeric label
        >>> ed = coll.fetch(site=150, first=True)

        >>> # Enforce a helpful error when not found
        >>> coll.fetch(site='S999', first=True, strict=True)
        Traceback (most recent call last):
            ...
        SiteError: No station found. Try one of: ...

        >>> # Fetch by coordinates within tolerance
        >>> ed = coll.fetch(lat=-22.37, lon=139.19, first=True)

        >>> # Filter by metadata on >HEAD (case-insensitive)
        >>> zs = coll.fetch(acqby='Zonge')

        See Also
        --------
        pycsamt.core.config.StationNamePolicy
        pycsamt.core.config.ensure_station
        pycsamt.seg.collection.EDICollection.stations
        pycsamt.seg.collection.EDICollection._resolve

        References
        ----------
        .. [1] SEG EDI 1.0 Electromagnetic Data Interchange
               specification.
        """


        matches: list[EDIFile] = []

        site_cands = (
            set(_site_candidates(site, policy))
            if site is not None else set()
        )

        for edi in self:
            head = edi.get_section("head")
            ok = True

            # --- station match (case-insensitive, with variants)
            if site is not None:
                est = getattr(edi, "station", None)
                if not (
                    est and est.upper() in site_cands
                ):
                    ok = False

            # --- geo match with tolerance
            if lat is not None and ok:
                hlat = getattr(head, "lat", None)
                if hlat is None or abs(hlat - lat) > tol:
                    ok = False

            if lon is not None and ok:
                hlon = getattr(head, "long", None)
                if hlon is None or abs(hlon - lon) > tol:
                    ok = False

            # --- arbitrary attribute matches
            for key, val in kwargs.items():
                if not ok:
                    break
                aval = getattr(edi, key.lower(), None)
                if aval is None and head is not None:
                    aval = getattr(head, key.lower(), None)
                if isinstance(aval, str) and isinstance(val, str):
                    if aval.upper() != val.upper():
                        ok = False
                elif aval != val:
                    ok = False

            if ok:
                matches.append(edi)

        if first:
            if matches:
                return matches[0]
            if strict:
                all_sites = ", ".join(self.stations())
                raise SiteError(
                    "No station found. Try one of: "
                    f"{all_sites}\nHint: call "
                    "`<collection>.stations()` to list."
                )
            return None

        # first=False
        if not matches and strict:
            all_sites = ", ".join(self.stations())
            raise SiteError(
                "No stations matched your query. Try one of: "
                f"{all_sites}\nHint: call "
                "`<collection>.stations()` to list."
            )
        return matches

    @staticmethod
    def _site_of(it: EDIFile) -> str | None:
        """Safely extract the site/station
        name from an EDIFile."""
        return getattr(it, "station", None)

    @staticmethod
    def _lat_of(it: EDIFile) -> float | None:
        """Safely extract the latitude from
        an EDIFile's HEAD section."""
        head = it.get_section("head")
        v = getattr(head, "lat", None)
        return float(v) if isinstance(v, (int, float)) else None

    @staticmethod
    def _lon_of(it: EDIFile) -> float | None:
        """Safely extract the longitude from
        an EDIFile's HEAD section."""
        head = it.get_section("head")
        v = getattr(head, "long", None)
        return float(v) if isinstance(v, (int, float)) else None

    @property
    def sites(self) -> list[str]:
        """Returns a list of best-effort
        station names for each item."""
        out: list[str] = []
        for it in self._items:
            s = self._site_of(it)
            out.append(s if s is not None else "")
        return out

    @property
    def latitude(self) -> np.ndarray:
        """Returns a numpy array of latitudes
        (np.nan for missing values)."""
        vals: list[float] = []
        for it in self._items:
            v = self._lat_of(it)
            vals.append(np.nan if v is None else float(v))
        return np.asarray(vals, dtype=float)

    @property
    def longitude(self) -> np.ndarray:
        """Returns a numpy array of longitudes
        (np.nan for missing values)."""
        vals: list[float] = []
        for it in self._items:
            v = self._lon_of(it)
            vals.append(np.nan if v is None else float(v))
        return np.asarray(vals, dtype=float)

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"EDICollection(n={len(self)}, "
            f"stations={self.stations()!r})"
        )


    def _summary_stats(self, summary_data: list[dict]) -> str:
        """Creates a statistical summary block from summary data."""
        if not summary_data:
            return "  (No statistics available for an empty collection)\n"

        total_files = len(summary_data)
        with_tipper = sum(1 for r in summary_data if r.get("tipper"))
        with_spectra = sum(1 for r in summary_data if r.get("spectra"))
        with_ts = sum(1 for r in summary_data if r.get("ts"))

        # Get all frequencies from all files to find the true min/max
        all_freqs = np.concatenate([
            edi.Z.freq for edi in self if edi.Z.freq is not None
            and edi.Z.freq.size > 0
        ]) if total_files > 0 else np.array([])

        lats = self.latitude[~np.isnan(self.latitude)]
        lons = self.longitude[~np.isnan(self.longitude)]

        freq_range = (
            f"Min={np.min(all_freqs):.2E}, Max={np.max(all_freqs):.2E}"
            if all_freqs.size > 0 else "N/A"
        )
        lat_range =(
            f"{min(lats):.4f} to {max(lats):.4f}"
            if len(lats) > 0 else "N/A"
        )
        lon_range = (
            f"{min(lons):.4f} to {max(lons):.4f}"
            if len(lons) > 0 else "N/A"
        )

        lines = [
            "  " + "-"*65,
            "  Statistical Summary:",
            f"    Total Sites: {total_files}",
            f"    Content Counts: Tipper={with_tipper},"
            f" Spectra={with_spectra}, TimeSeries={with_ts}",
            f"    Frequency Range (Hz): {freq_range}",
            f"    Latitude Range:       {lat_range}",
            f"    Longitude Range:      {lon_range}",
            "  " + "-"*65,
        ]
        return "\n".join(lines)


    def __str__(self) -> str:  # pragma: no cover
        """Provides a detailed and statistical summary of the collection."""
        summary_data = self.summary()
        if not summary_data:
            return "EDICollection(n=0, stations=[])"

        # --- Header ---
        title = f" EDICollection Summary (Total Sites: {len(summary_data)}) "
        width = 70
        header = [
            "=" * width,
            title.center(width),
            "=" * width
        ]

        # --- Per-Site Details ---
        details = ["\nSite Details:"]
        max_station_len = max(
            [len(r['station']) for r in summary_data] + [len("Station")])

        table_header = (
            f"  {'Station'.ljust(max_station_len)} | Freqs"
            " | Tipper | Spectra | TimeSeries"
        )
        table_width = len(table_header) - 2
        details.append(table_header)
        details.append("  " + "-" * table_width)

        for r in summary_data:
            tip = "Y" if r["tipper"] else "N"
            sp = "Y" if r["spectra"] else "N"
            ts = "Y" if r["ts"] else "N"

            details.append(
                f"  {r['station'].ljust(max_station_len)} | "
                f"{r['n_freq']:<5} | {tip:^6} | {sp:^7} | {ts:^10}"
            )

        details.append("  " + "-" * table_width)
        details.append(
            "  *Y = Yes (data present); *N = No (data not present)"
            )

        # --- Statistical Summary ---
        stats_str = self._summary_stats(summary_data)

        # --- Combine and Return ---
        return "\n".join(header + [stats_str] + details)


def _site_candidates(
    site: Any,
    policy: StationNamePolicy | None = None,
) -> list[str]:
    pol = policy or get_config().station_policy

    def _try_float(v):
        try:
            return float(v)
        except Exception:
            return None

    s = str(site).strip()
    cands = {s.upper()}

    x = _try_float(s)
    if x is not None:
        n = int(round(x))
        cands.add(str(n).upper())
        cands.add(
            f"{pol.prefix}{n:0{pol.pad}d}".upper()
        )
        # common “×10” slip: try n*10 and n/10
        cands.add(str(n * 10).upper())
        cands.add(
            f"{pol.prefix}{n * 10:0{pol.pad}d}".upper()
        )
        if n % 10 == 0:
            m = n // 10
            cands.add(str(m).upper())
            cands.add(
                f"{pol.prefix}{m:0{pol.pad}d}".upper()
            )
    return list(cands)
