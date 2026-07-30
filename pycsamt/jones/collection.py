# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import (
    Callable,
    Union,
)

import numpy as np
import pandas as pd

from ..log.logger import get_logger
from .cbase import JCBBase, JCoreParser, JParseMixin
from .j import JFile

logger = get_logger(__name__)

__all__ = ["JCollectionMixin", "JCollection"]

Pathish = Union[str, Path]


class JCollectionMixin(JParseMixin):
    r"""
    Mixin that provides folder/glob expansion and robust
    parsing orchestration for Jones J-format collections.

    The mixin augments a stateful base (e.g.,
    :class:`JCBBase`) with helpers to normalize user input
    paths, expand folders and wildcards, deduplicate
    candidates, and parse items using a tolerant strategy.

    Implementations typically use
    :class:`~pycsamt.jones.j.JFile` as the per-item reader,
    but the mixin is agnostic to the concrete class as long
    as it exposes a compatible ``from_file`` constructor.

    Notes
    -----
    The design keeps discovery and parsing separate.  First
    the mixin expands and filters path candidates; then the
    owning class decides how to instantiate items and how to
    record failures.  This separation avoids tight coupling
    and keeps error handling clear.

    The mixin favors **permissive** heuristics while still
    catching obvious mistakes (non-files, unreadable paths,
    extensions that are unlikely to be J-format, etc.).

    Typical Methods
    ---------------
    from_paths(paths, *, strict=False, verbose=0)
        Expand and parse ``paths`` into items.  Returns the
        created items (usually a list).
    _iter_paths(paths)
        Yield :class:`pathlib.Path` objects from strings or
        paths, ignoring duplicates.
    _iter_j_files(paths, *, recursive=True)
        Yield J-candidates by expanding folders and glob
        patterns (e.g., ``"*.j"``, ``"**/*.txt"``).
    _is_j_candidate(path)
        Light heuristic that recognizes J-format candidates.

    Examples
    --------
    Expand a folder and parse all candidates:

    >>> mix = JCollectionMixin()
    >>> # usually inherited together with JCBBase
    >>> list(mix._iter_j_files(["data/j"]))[:1]  # doctest: +ELLIPSIS
    [PosixPath('...kb0-s001.txt')]

    See Also
    --------
    JCBBase
        Minimal stateful base for collections (stores items).
    JCollection
        High-level collection that mixes in this helper and
        provides user-facing APIs.
    pycsamt.jones.j.JFile
        High-level J reader used as default item type.

    References
    ----------
    .. [1] A. G. Jones (1994). *J-format v2.0*. MTNet notes.
    """

    def add_from(
        self,
        sources: Pathish | Sequence[Pathish],
        *,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int | None = None,
    ) -> JCollectionMixin:
        v = self.verbose if verbose is None else int(verbose)
        pr = JCoreParser(
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=v,
        )
        items = pr.parse(sources)
        for jf in items:
            self.add(jf)  # type: ignore[attr-defined]
        errs = pr.errors()
        if errs:
            logger.info("Added with %d errors.", len(errs))
        return self

    def select(self, stations: Sequence[str]) -> JCollection:
        keep = {str(s) for s in stations}
        out = JCollection(verbose=self.verbose)
        for jf in self:  # type: ignore[operator]
            sid = getattr(jf, "site", None) or "-"
            if sid in keep:
                out.add(jf)
        return out

    def where(self, fn: Callable[[JFile], bool]) -> JCollection:
        out = JCollection(verbose=self.verbose)
        for jf in self:  # type: ignore[operator]
            if fn(jf):
                out.add(jf)
        return out

    def sort(
        self,
        key: str | Callable[[JFile], object] = "station",
        *,
        reverse: bool = False,
    ) -> JCollection:
        def _key(jf: JFile) -> object:
            if callable(key):
                return key(jf)
            if key in {"station", "site", "name"}:
                return (
                    getattr(jf, "site", None)
                    or getattr(jf, "name", None)
                    or ""
                )
            if key == "n_freq":
                f = getattr(jf, "freq", None)
                return int(getattr(f, "size", len(f or [])))
            if key == "path":
                return str(jf.path) if jf.path else ""
            if key == "lat":
                return float(getattr(jf, "lat", 9e9) or 9e9)
            if key == "lon":
                return float(getattr(jf, "lon", 9e9) or 9e9)
            return getattr(jf, key, None)

        out = JCollection(
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
        for jf in self:  # type: ignore[operator]
            p.append(str(jf.path) if jf.path else "-")
        return p

    def nf_stats(self) -> dict:
        vals: list[int] = []
        for jf in self:  # type: ignore[operator]
            f = getattr(jf, "freq", None)
            vals.append(int(getattr(f, "size", 0)))
        a = np.asarray(vals, int)
        if a.size == 0:
            return {"min": 0, "max": 0, "mean": 0.0}
        return {
            "min": int(a.min()),
            "max": int(a.max()),
            "mean": float(a.mean()),
        }


class JCollection(JCBBase, JCollectionMixin):
    r"""
    High-level collection of Jones J-format files.

    Combines :class:`JCBBase` (state container) with
    :class:`JCollectionMixin` (discovery/orchestration) to
    offer a convenient interface for batch loading, quick
    metadata browsing, and simple filtering over many J files.

    The default per-item object is
    :class:`~pycsamt.jones.j.JFile`, which exposes parsed
    headers (:class:`~pycsamt.jones.heads.Heads`) and data
    blocks (:class:`~pycsamt.jones.blocks.JBlocks`) as well
    as MT objects (:class:`pycsamt.z.z.Z`,
    :class:`pycsamt.z.tipper.Tipper`,
    :class:`pycsamt.z.resphase.ResPhase`) when available.

    Parameters
    ----------
    verbose : int, default 0
        Verbosity for logging and warnings during discovery
        and parsing.

    Attributes
    ----------
    items : list of JFile
        Parsed items in insertion order.  The class follows
        Python's container protocol (``__len__``, ``__iter__``).
    n : int
        Number of items in the collection (``len(self)``).
    paths : list of Path
        Convenience view of the associated file paths.
    stations : list of str
        Station codes when available on items (best-effort).

    Methods
    -------
    parse(paths, *, strict=False, verbose=None)
        Expand ``paths`` (file, folder, glob), parse J files,
        and append items to the collection.  Returns the new
        items.
    where(**query)
        Optional, return a filtered view (e.g., by ``station``
        or by component family ``has_Z/has_R/has_T``).
    summary()
        Optional, return a compact text or dataframe summary.
    write_index(path, *, overwrite=True)
        Optional, serialize a CSV/TSV inventory for the
        collection.

    Notes
    -----
    The collection is intentionally lightweight.  It does not
    enforce a database schema, and it avoids implicitly
    reading large arrays until needed.  This makes it suitable
    for quick exploration and CI tests.

    Error handling is conservative: unreadable paths or parse
    errors are usually skipped with a warning (unless
    ``strict=True`` is requested).

    Examples
    --------
    Parse a folder and browse basic metadata:

    >>> col = JCollection(verbose=0)
    >>> _ = col.parse(["data/j"])
    >>> len(col) >= 1
    True
    >>> sorted(set(getattr(x, "site", None) for x in col))[
    ...     :3
    ... ]  # doctest: +ELLIPSIS
    [...]

    Filter by station (if ``where`` is provided):

    >>> sub = getattr(col, "where", lambda **k: col)(station="KB0001")
    >>> len(sub) >= 1
    True

    See Also
    --------
    JCBBase
        Underlying state container.
    JCollectionMixin
        Path expansion and tolerant orchestration.
    pycsamt.jones.j.JFile
        Per-item high-level J reader.

    References
    ----------
    .. [1] A. G. Jones (1994). *J-format v2.0*. MTNet notes.
    """

    def __init__(
        self,
        items: Iterable[JFile] | None = None,
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
    ) -> JCollection:
        pr = JCoreParser(
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
        other: JCollection,
        *,
        on_dup: str = "replace",
    ) -> JCollection:
        out = JCollection(
            items=list(self),
            verbose=self.verbose,
        )
        if on_dup not in {"replace", "keep"}:
            raise ValueError("on_dup must be keep|replace")
        for jf in other:
            sid = getattr(jf, "site", None)
            sid = sid or (str(jf.path) if jf.path else None)
            if sid in out._index and on_dup == "keep":
                continue
            out.add(jf)
        return out

    def _resolve(self, site: str) -> JFile:
        """Find item by key, site/station, stem,
        or filename (case-insensitive)."""
        site_upper = str(site).upper()
        # Case-insensitive check on the
        # index dictionary first for speed
        for key, idx in self._index.items():
            if str(key).upper() == site_upper:
                return self._items[idx]

        # Fallback to iterating through all items if not in index
        for jf in self:
            sid = getattr(jf, "site", None)
            if sid and str(sid).upper() == site_upper:
                return jf

        # Fallback to path matching (case-insensitive)
        site_lower = str(site).lower()
        for jf in self:
            p = getattr(jf, "path", None)
            if p is None:
                continue
            if (
                p.stem.lower() == site_lower
                or p.name.lower() == site_lower
                or str(p).lower() == site_lower
            ):
                return jf
        raise KeyError(f"site not found: {site!r}")

    def _heads(self, jf: JFile):
        return getattr(jf, "heads", None)

    def get(
        self,
        site: str,
        what: str,
        default: object | None = None,
    ) -> object | None:
        """
        Quick extractor.  'what' supports:
          - 'freq'
          - 'z','zxx','zxy','zyx','zyy'
          - 'tip','tx','ty'
          - 'rxy','ryx','rxx','ryy' (rho)
          - 'phixy','phiyx','phixx','phiyy'
          - 'station'/'site','lat','lon','elev','az','name'
          - 'path','filename'
        """
        try:
            jf = self._resolve(site)
        except KeyError:
            return default

        w = str(what).lower()

        # --------------- freq ---------------
        if w == "freq":
            return getattr(jf, "freq", None)

        # --------------- Z ------------------
        if w == "z":
            return getattr(getattr(jf, "Z", None), "z", None)
        if w in {"zxx", "zxy", "zyx", "zyy"}:
            m = {"zxx": (0, 0), "zxy": (0, 1), "zyx": (1, 0), "zyy": (1, 1)}[w]
            z = getattr(getattr(jf, "Z", None), "z", None)
            if z is None:
                return default
            a = np.asarray(z)
            return a[:, m[0], m[1]]

        # ------------- tipper ---------------
        if w in {"tip", "tx", "ty"}:
            tp = getattr(jf, "Tip", None)
            if tp is None:
                return default
            arr = getattr(tp, "tipper", None)
            if arr is None:
                return default
            a = np.asarray(arr)
            if a.ndim == 3 and a.shape[1:] == (1, 2):
                a = a[:, 0, :]
            if w == "tip":
                return a
            return a[:, 0] if w == "tx" else a[:, 1]

        # -------- resistivity/phase --------
        if w in {
            "rxx",
            "rxy",
            "ryx",
            "ryy",
            "phixx",
            "phixy",
            "phiyx",
            "phiyy",
        }:
            rp = getattr(jf, "Res", None)
            if rp is None:
                return default
            if w.startswith("phi"):
                key = "phase_" + w[3:]
            else:
                key = "res_" + w
            return getattr(rp, key, None)

        # ----------- headers/site ----------
        if w in {"station", "site"}:
            return getattr(jf, "site", None)
        if w == "name":
            return getattr(jf, "name", None)
        if w == "lat":
            return getattr(jf, "lat", None)
        if w == "lon":
            return getattr(jf, "lon", None)
        if w == "elev":
            return getattr(jf, "elev", None)
        if w == "az":
            return getattr(jf, "azimuth", None)

        # -------------- path-ish ------------
        if w in {"path", "filename"}:
            p = getattr(jf, "path", None)
            if p is None:
                return default
            return str(p) if w == "path" else p.name

        return default

    def set(
        self,
        site: str,
        *,
        jfile: JFile | None = None,
        update: dict[str, object] | None = None,
    ) -> JFile:
        """
        Replace or mutate a site's :class:`JFile`.

        - If 'jfile' is given, replace the stored object.
        - Else apply 'update' keys on the existing object:
          'station'/'site','lat','lon','elev','az',
          'freq','z','tip','resphase'
        """
        jf = self._resolve(site)

        if jfile is not None:
            try:
                self.remove(site)  # optional on CBBase
            except Exception:
                pass
            try:
                self.add(jfile)
            except Exception:
                idx = getattr(self, "_index", None)
                if isinstance(idx, dict):
                    idx[site] = jfile
            return jfile

        upd = update or {}
        heads = self._heads(jf)
        info = getattr(heads, "info", None)

        if "station" in upd or "site" in upd:
            sid = str(upd.get("station", upd.get("site")))
            try:
                # prefer JFile mutable name
                jf.station = sid  # type: ignore[attr-defined]
            except Exception:
                # fallback: set on head if present
                try:
                    heads.head.station = sid  # type: ignore[attr-defined]
                except Exception:
                    pass

        def _set_info(key: str, val: object) -> None:
            if info is None:
                return
            try:
                info.items[key] = str(val)
                # drop cached parsed site if any
                info._site_cache = None
            except Exception:
                pass

        if "lat" in upd:
            _set_info("LATITUDE", upd["lat"])
        if "lon" in upd:
            _set_info("LONGITUDE", upd["lon"])
        if "elev" in upd:
            _set_info("ELEVATION", upd["elev"])
        if "az" in upd:
            _set_info("AZIMUTH", upd["az"])

        if "freq" in upd and getattr(jf, "Z", None) is not None:
            jf.Z._freq = np.asarray(upd["freq"], float)
        if "z" in upd and getattr(jf, "Z", None) is not None:
            jf.Z._z = np.asarray(upd["z"], complex)
        if "tip" in upd and getattr(jf, "Tip", None) is not None:
            jf.Tip._tipper = np.asarray(upd["tip"], complex)

        if "resphase" in upd and getattr(jf, "Res", None) is not None:
            rp = np.asarray(upd["resphase"], float)
            try:
                jf.Res.resistivity = rp  # type: ignore[attr-defined]
            except Exception:
                pass

        return jf

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
    ) -> JFile:
        """
        Shift or set position and/or rename station.
        """
        jf = self._resolve(site)
        heads = self._heads(jf)
        info = getattr(heads, "info", None)
        if info is None:
            raise ValueError("no >INFO to adjust")

        if rename:
            try:
                jf.station = str(rename)  # type: ignore[attr-defined]
            except Exception:
                try:
                    heads.head.station = str(rename)  # type: ignore[attr-defined]
                except Exception:
                    pass

        def _getf(k: str, default: float) -> float:
            v = info.items.get(k)
            try:
                return float(v) if v is not None else default
            except Exception:
                return default

        if lat is not None:
            info.items["LATITUDE"] = str(float(lat))
        if lon is not None:
            info.items["LONGITUDE"] = str(float(lon))
        if elev is not None:
            info.items["ELEVATION"] = str(float(elev))

        if dlat is not None:
            cur = _getf("LATITUDE", 0.0)
            info.items["LATITUDE"] = str(cur + float(dlat))
        if dlon is not None:
            cur = _getf("LONGITUDE", 0.0)
            info.items["LONGITUDE"] = str(cur + float(dlon))

        info._site_cache = None
        return jf

    # ----------------- summary / repr ----------------
    def summary(
        self,
        fields: Sequence[str] = (
            "station",
            "n_freq",
            "has_z",
            "has_r",
            "has_t",
            "lat",
            "lon",
            "az",
        ),
    ) -> list[dict]:
        rows: list[dict] = []
        for jf in self:
            f = getattr(jf, "freq", None)
            n = int(getattr(f, "size", 0))
            row = {
                "station": getattr(jf, "site", None)
                or getattr(jf, "name", None)
                or "-",
                "n_freq": n,
                "has_z": bool(getattr(jf, "Z", None)),
                "has_r": bool(getattr(jf, "Res", None)),
                "has_t": bool(getattr(jf, "Tip", None)),
                "lat": getattr(jf, "lat", None),
                "lon": getattr(jf, "lon", None),
                "az": getattr(jf, "azimuth", None),
            }
            rows.append({k: row.get(k) for k in fields})
        return rows

    @staticmethod
    def _site_of(it) -> str | None:
        # Prefer common single-file attributes
        for nm in ("site", "station", "name"):
            v = getattr(it, nm, None)
            if isinstance(v, str) and v.strip():
                return v
        # Fallback to nested headers if present
        hd = getattr(it, "heads", None)
        sta = getattr(getattr(hd, "head", None), "station", None)
        if isinstance(sta, str) and sta.strip():
            return sta
        return None

    @staticmethod
    def _lat_of(it) -> float | None:
        # JFile exposes .lat; some objects use .latitude
        v = getattr(it, "lat", None)
        if v is None:
            v = getattr(it, "latitude", None)
        if v is None:
            hd = getattr(it, "heads", None)
            v = getattr(hd, "latitude", None)
        return float(v) if isinstance(v, (int, float)) else None

    @staticmethod
    def _lon_of(it) -> float | None:
        # JFile exposes .lon; some objects use .longitude
        v = getattr(it, "lon", None)
        if v is None:
            v = getattr(it, "longitude", None)
        if v is None:
            hd = getattr(it, "heads", None)
            v = getattr(hd, "longitude", None)
        return float(v) if isinstance(v, (int, float)) else None

    @property
    def sites(self) -> list[str]:
        """Best-effort station names (one per item)."""
        out: list[str] = []
        for it in self._items:
            s = self._site_of(it)
            out.append(s if s is not None else "")
        return out

    @property
    def latitude(self) -> np.ndarray:
        """Vector of latitudes (np.nan for missing)."""
        vals: list[float] = []
        for it in self._items:
            v = self._lat_of(it)
            vals.append(np.nan if v is None else float(v))
        return np.asarray(vals, dtype=float)

    @property
    def longitude(self) -> np.ndarray:
        """Vector of longitudes (np.nan for missing)."""
        vals: list[float] = []
        for it in self._items:
            v = self._lon_of(it)
            vals.append(np.nan if v is None else float(v))
        return np.asarray(vals, dtype=float)

    def export(
        self,
        output_dir: Pathish,
        *,
        file_pattern: str = "{station}.j",
        export_summary: bool = False,
        summary_filename: str = "summary.csv",
        **jfile_write_kwargs,
    ) -> dict:
        """
        Exports all JFile items to a directory with advanced options.

        This method orchestrates the writing of each JFile in the collection
        to a specified directory, with flexible naming, error handling, and
        an optional summary CSV file.

        Parameters
        ----------
        output_dir : Pathish
            The path to the directory where files will be saved. It will
            be created if it does not exist.
        file_pattern : str, default="{station}.j"
            A format string for output filenames. Can use attributes of
            JFile like `{station}`, `{name}`, `{site}`.
        export_summary : bool, default=False
            If True, a summary of the collection will also be saved as a
            CSV file in the output directory.
        summary_filename : str, default="summary.csv"
            The name of the summary file if `export_summary` is True.
        **jfile_write_kwargs
            Keyword arguments to be passed directly to each `JFile.write()`
            call (e.g., `datatype="ZRT"`, `overwrite=True`).

        Returns
        -------
        dict
            A dictionary with two keys:
            - 'successful': A list of paths to successfully written files.
            - 'failed': A list of tuples, where each tuple contains the
              station name and the error that occurred.
        """
        out_dir = Path(str(output_dir)).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        successful_paths = []
        failed_items = []

        # Use tqdm for a progress bar if available
        items_iterator = self._items
        try:
            from tqdm import tqdm

            items_iterator = tqdm(self._items, desc="Exporting J-Files")
        except (NameError, ImportError):
            pass  # tqdm not installed

        for jf in items_iterator:
            sid = jf.site or "unknown_station"
            try:
                filename = file_pattern.format(
                    station=sid, site=sid, name=jf.name
                )
                output_path = out_dir / filename

                # Delegate the actual writing to the JFile instance
                written_path = jf.write(
                    new_jfn=str(output_path), **jfile_write_kwargs
                )
                successful_paths.append(written_path)
            except Exception as e:
                failed_items.append((sid, e))
                logger.error(f"Failed to write J-file for station {sid}: {e}")

        # Export the summary CSV if requested
        if export_summary:
            try:
                summary_df = pd.DataFrame(self.summary())
                summary_path = out_dir / summary_filename
                summary_df.to_csv(summary_path, index=False)
                successful_paths.append(str(summary_path))
            except Exception as e:
                failed_items.append(("summary.csv", e))
                logger.error(f"Failed to write summary CSV: {e}")

        return {"successful": successful_paths, "failed": failed_items}

    def fetch(
        self,
        site: str | None = None,
        lat: float | None = None,
        lon: float | None = None,
        tol: float = 0.001,
        first: bool = False,
        **kwargs,
    ) -> JFile | list[JFile] | None:
        """
        Fetches JFile objects from the collection based on specified criteria.

        This method provides a flexible way to search for J-format files
        by site name, geographic coordinates, or any other attribute of
        the JFile object or its nested Heads sections.

        Parameters
        ----------
        site : str, optional
            The site or station name to search for. The comparison is
            case-insensitive.
        lat : float, optional
            The latitude to search for, in decimal degrees.
        lon : float, optional
            The longitude to search for, in decimal degrees.
        tol : float, default=0.001
            The tolerance in decimal degrees for geographic coordinate
            searches. A match is found if the absolute difference is
            within this tolerance.
        first : bool, default=False
            If True, returns only the first matching JFile object found,
            or None if no match is found. If False, returns a list of
            all matching objects.
        **kwargs : Any
            Additional keyword arguments to match against attributes of
            the JFile object or its nested Heads object (e.g.,
            `acqby='Contractor'`). The attribute name is case-insensitive.

        Returns
        -------
        JFile or list of JFile or None
            - If `first=True`, returns the first matching JFile or None.
            - If `first=False`, returns a list of all matching JFile objects.
              An empty list is returned if no matches are found.

        Examples
        --------
        >>> # Fetch a single site by its name
        >>> jfile_obj = jcollection.fetch(site="S01", first=True)

        >>> # Fetch all sites where the azimuth is 0
        >>> zero_az_files = jcollection.fetch(azimuth=0)

        >>> # Fetch all sites within a small geographic area
        >>> area_files = jcollection.fetch(lat=26.05, lon=-10.33, tol=0.1)
        """
        matches = []

        for jf in self:
            is_match = True

            # --- Match by site name (case-insensitive) ---
            if site is not None:
                jf_site = getattr(jf, "site", None)
                if not (jf_site and jf_site.upper() == site.upper()):
                    is_match = False

            # --- Match by geographic coordinates with tolerance ---
            if lat is not None and is_match:
                jf_lat = getattr(jf, "lat", None)
                if jf_lat is None or abs(jf_lat - lat) > tol:
                    is_match = False

            if lon is not None and is_match:
                jf_lon = getattr(jf, "lon", None)
                if jf_lon is None or abs(jf_lon - lon) > tol:
                    is_match = False

            # --- Match by other arbitrary attributes (case-insensitive) ---
            for key, value in kwargs.items():
                if not is_match:
                    break

                # Check for attribute directly on JFile object first
                attr_val = getattr(jf, key.lower(), None)

                # If not found, check on the nested head and info objects
                if attr_val is None and jf.heads and jf.heads.head:
                    attr_val = getattr(jf.heads.head, key.lower(), None)
                if attr_val is None and jf.heads and jf.heads.info:
                    attr_val = getattr(jf.heads.info, key.lower(), None)

                # Perform comparison
                if isinstance(attr_val, str) and isinstance(value, str):
                    if attr_val.upper() != value.upper():
                        is_match = False
                elif attr_val != value:
                    is_match = False

            if is_match:
                matches.append(jf)

        if first:
            return matches[0] if matches else None

        return matches

    def _summary_stats(self, summary_data: list[dict]) -> str:
        """Creates a statistical summary block from summary data."""
        if not summary_data:
            return "  (No statistics available for an empty collection)\n"

        total_files = len(summary_data)
        with_z = sum(1 for r in summary_data if r.get("has_z"))
        with_r = sum(1 for r in summary_data if r.get("has_r"))
        with_t = sum(1 for r in summary_data if r.get("has_t"))

        # Get all frequencies from all files to find the true min/max
        all_freqs = (
            np.concatenate(
                [
                    jf.freq
                    for jf in self
                    if jf.freq is not None and jf.freq.size > 0
                ]
            )
            if total_files > 0
            else np.array([])
        )

        lats = self.latitude[~np.isnan(self.latitude)]
        lons = self.longitude[~np.isnan(self.longitude)]

        freq_range = (
            f"Min={np.min(all_freqs):.2E}, Max={np.max(all_freqs):.2E}"
            if all_freqs.size > 0
            else "N/A"
        )
        lat_range = (
            f"{min(lats):.4f} to {max(lats):.4f}" if len(lats) > 0 else "N/A"
        )
        lon_range = (
            f"{min(lons):.4f} to {max(lons):.4f}" if len(lons) > 0 else "N/A"
        )

        lines = [
            "  " + "-" * 68,
            "  Statistical Summary:",
            f"    Total Sites: {total_files}",
            f"    Component Counts: Z={with_z}, R={with_r}, T={with_t}",
            f"    Frequency Range (Hz): {freq_range}",
            f"    Latitude Range:       {lat_range}",
            f"    Longitude Range:      {lon_range}",
            "  " + "-" * 68,
        ]
        return "\n".join(lines)

    def __str__(self) -> str:  # pragma: no cover
        """Provides a detailed and statistical summary of the collection."""
        summary_data = self.summary()
        if not summary_data:
            return "JCollection(n=0, stations=[])"

        # --- Header ---
        title = f" JCollection Summary (Total Sites: {len(summary_data)}) "
        width = 72
        header = ["=" * width, title.center(width), "=" * width]

        # --- Per-Site Details ---
        details = ["\nSite Details:"]
        # Calculate max station length, ensuring
        # it's at least as wide as the header
        max_station_len = max(
            [len(r["station"]) for r in summary_data] + [len("Station")]
        )

        table_header = (
            f"  {'Station'.ljust(max_station_len)} | Freqs"
            " | Z | R | T | Lat       | Lon        | Az    "
        )
        table_width = len(table_header) - 2
        details.append(table_header)
        details.append("  " + "-" * table_width)

        for r in summary_data:
            tip = "Y" if r["has_t"] else "N"
            zz = "Y" if r["has_z"] else "N"
            rr = "Y" if r["has_r"] else "N"
            lat_str = f"{r['lat']:.4f}" if r["lat"] is not None else "N/A"
            lon_str = f"{r['lon']:.4f}" if r["lon"] is not None else "N/A"
            az_str = f"{r['az']:.1f}" if r["az"] is not None else "N/A"

            details.append(
                f"  {r['station'].ljust(max_station_len)} | "
                f"{r['n_freq']:<5} | {zz:^1} | {rr:^1} | {tip:^1} | "
                f"{lat_str:<9} | {lon_str:<10} | {az_str}"
            )

        details.append("  " + "-" * table_width)
        details.append("  *Y = Yes (data present); *N = No (data not present)")

        # --- Statistical Summary ---
        stats_str = self._summary_stats(summary_data)

        # --- Combine and Return ---
        return "\n".join(header + [stats_str] + details)

    def __repr__(self) -> str:  # pragma: no cover
        return f"JCollection(n={len(self)}, stations={self.stations()!r})"
