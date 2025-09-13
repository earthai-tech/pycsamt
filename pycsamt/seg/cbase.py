# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Tuple,
    Union,
    Callable,
)

import numpy as np

from ..log.logger import get_logger
from .validation import IsEdi
from .edi import EDIFile

logger = get_logger(__name__)

__all__ = ["ParseMixin", "CoreParser", "CBBase"]


Pathish = Union[str, Path]
SrcType = Union[Pathish, Sequence[Pathish]]


class ParseMixin:
    r"""
    Helpers for discovering and normalizing EDI sources.

    The mixin accepts files, folders, and glob patterns,
    and yields only paths that point to ``.edi`` files.
    All path inputs are normalized with user-home and
    relative segments resolved.

    Parameters
    ----------
    None
        This is a mixin.  It does not define its own public
        constructor parameters.

    Attributes
    ----------
    EDI_SUFFIXES : set of str
        Allowed filename suffixes.  Defaults to
        ``{'.edi'}``.
    recursive : bool
        Expected to be provided by the host class.  If
        ``True``, directory searches use ``rglob``.
    _errors : list of tuple(Path, BaseException)
        Optional sink for discovery issues.  When present,
        helpers record unmatched patterns or missing files.

    Notes
    -----
    The mixin exposes small helpers that higher-level
    parsers can reuse:

    * ``_as_path`` converts any pathish value into an
      absolute :class:`pathlib.Path`.
    * ``_is_edi_path`` returns ``True`` for files whose
      suffix belongs to ``EDI_SUFFIXES``.
    * ``_iter_paths`` normalizes a single source or a
      sequence of sources into absolute paths.
    * ``_iter_edi_files`` walks files, directories, and
      glob patterns and yields only existing EDI files.
    * ``_fast_station`` performs a cheap scan of ``>HEAD``
      to locate ``DATAID``, when available.
    * ``_push_error`` records discovery problems into the
      host ``_errors`` list if present.

    Examples
    --------
    Basic file enumeration::

        class Finder(ParseMixin):
            recursive = True

        f = Finder()
        paths = list(f._iter_edi_files(
            ['data', 'logs/*.edi'], root=Path('.')
        ))

    See Also
    --------
    CoreParser
        High-level parser that builds :class:`EDIFile`
        objects and aggregates errors.
    EDIFile
        Reader and writer for single EDI files.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """


    EDI_SUFFIXES = {".edi"}

    def _as_path(self, p: Pathish) -> Path:
        return Path(str(p)).expanduser().resolve()

    def _is_edi_path(self, p: Path) -> bool:
        return p.is_file() and (p.suffix.lower() in self.EDI_SUFFIXES)

    def _iter_paths(self, src: SrcType) -> Iterator[Path]:
        if isinstance(src, (str, Path)):
            yield self._as_path(src)  # type: ignore[arg-type]
            return
        for s in src:
            yield self._as_path(s)

    def _push_error(self, src: Pathish, msg: str) -> None:
        """
        Record a parsing discovery error on the host if it
        exposes `_errors`. Create the store if missing.
        """
        store = getattr(self, "_errors", None)
        if store is None:
            setattr(self, "_errors", [])  # type: ignore[attr-defined]
            store = getattr(self, "_errors")
        err = FileNotFoundError(msg)
        store.append((self._as_path(src), err))

    def _iter_edi_files(
        self,
        sources: list[str | Path],
        *,
        root: Path | None = None,
    ) -> Iterator[Path]:
        base = Path(root) if root is not None else Path.cwd()

        def _is_glob(s: str) -> bool:
            return any(ch in s for ch in "*?[]")

        for src in sources:
            p = Path(str(src))

            # direct file
            if p.exists() and p.is_file():
                yield p
                continue

            # directory -> (r)glob for *.edi
            if p.exists() and p.is_dir():
                it = p.rglob("*.edi") if getattr(
                    self, "recursive", True
                ) else p.glob("*.edi")
                for m in it:
                    if m.is_file():
                        yield m
                continue

            # glob pattern (absolute or relative)
            if _is_glob(str(p)):
                if p.is_absolute():
                    parent, name = p.parent, p.name
                    it = parent.rglob(name) if "**" in name else \
                         parent.glob(name)
                else:
                    pat = str(p)
                    it = base.rglob(pat) if "**" in pat else \
                         base.glob(pat)
                any_yielded = False
                for m in it:
                    if m.is_file():
                        any_yielded = True
                        yield m
                if not any_yielded:
                    self._push_error(
                        src, f"No match for pattern: {src}"
                    )
                continue

            # anything else: record as not found
            self._push_error(
                src, f"Not found or unsupported: {src}"
            )

    def _fast_station(self, p: Path) -> Optional[str]:
        """
        Quick scan for DATAID inside >HEAD.  Best effort.
        """
        try:
            with p.open("r", encoding="utf-8") as f:
                in_head = False
                for raw in f:
                    s = raw.strip()
                    if not s:
                        continue
                    if s.startswith(">=") or (s.startswith(">") and
                                              not s.upper().startswith(
                                                  ">HEAD"
                                              )):
                        if in_head:
                            break
                    if s.upper().startswith(">HEAD"):
                        in_head = True
                        continue
                    if in_head and s.upper().startswith("DATAID="):
                        return s.split("=", 1)[1].strip()
        except Exception:
            return None
        return None

@dataclass
class _ParseResult:
    path: Path
    edi: Optional[EDIFile]
    error: Optional[BaseException]


class CoreParser(ParseMixin):
    r"""
    Robust multi-source EDI parser with error tracking and
    duplicate policies.

    The parser consumes any combination of files, folders,
    or glob patterns, builds :class:`EDIFile` objects, and
    records failures.  Duplicates can be handled by station
    id or kept as they appear.

    Parameters
    ----------
    recursive : bool, default ``True``
        Recurse into subdirectories when scanning folders.
    strict : bool, default ``False``
        If ``True``, any read error is raised.  If
        ``False``, errors are collected and the parse
        continues.
    on_dup : {'replace', 'keep'}, default ``'replace'``
        Duplicate policy by station id.  With ``'replace'``,
        the last item wins.  With ``'keep'``, the first seen
        item is preserved.
    verbose : int, default ``0``
        Verbosity passed to :class:`EDIFile`.

    Attributes
    ----------
    results : list of _ParseResult
        Structured outcomes for each discovered input.  Each
        entry stores the path, the optional :class:`EDIFile`,
        and the read error if one occurred.
    _errors : list of tuple(Path, BaseException)
        Discovery issues, such as unmatched glob patterns or
        missing files.  Filled by :meth:`ParseMixin._push_error`.

    Returns
    -------
    list of EDIFile
        The parsed and de-duplicated EDI objects.

    Notes
    -----
    Station identity is taken from ``EDIFile.station`` when
    available.  If missing, a fast scan of ``>HEAD`` is used
    as a fallback.  When no station can be resolved, the file
    path string is used as the key.

    Examples
    --------
    Parse a folder and a glob, keep the first copy::

        cp = CoreParser(on_dup='keep', recursive=True)
        edis = cp.parse(['data/edi', 'more/*.edi'])
        errs = cp.errors()

    See Also
    --------
    ParseMixin
        Source discovery utilities used by the parser.
    EDIFile
        Single-file reader used to load each path.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """

    def __init__(
        self,
        *,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> None:
        self.recursive = bool(recursive)
        self.strict = bool(strict)
        self.on_dup = str(on_dup).lower()
        self.verbose = int(verbose)

        if self.on_dup not in {"replace", "keep"}:
            raise ValueError("on_dup must be keep|replace")

        self.results: List[_ParseResult] = []
        self._errors: List[Tuple[Path, BaseException]] = []

    # parse a single file (with error policy)
    def _read_one(self, p: Path) -> _ParseResult:
        try:
            IsEdi._assert_edi(p, deep=False)
            edi = EDIFile(p, verbose=self.verbose)
            return _ParseResult(path=p, edi=edi, error=None)
        except BaseException as exc:  # noqa: BLE001
            if self.strict:
                raise
            logger.debug("Skip %s: %s", p, exc)
            return _ParseResult(path=p, edi=None, error=exc)

    def parse(self, sources: SrcType) -> List[EDIFile]:
        edis: List[EDIFile] = []
        self.results.clear()
   
        self._errors.clear()

        for p in self._iter_edi_files(
            sources if isinstance(sources, list)
            else list(self._iter_paths(sources))
        ):
            res = self._read_one(p)
            self.results.append(res)
            if res.edi is None:
                continue
            edis.append(res.edi)
            
        # handle duplicates by station (if any)
        by_station: Dict[str, EDIFile] = {}
        seen: Dict[str, int] = {}
        out: List[EDIFile] = []

        for ed in edis:
            sid = getattr(ed, "station", None)
            if not sid:
                sid = ed._fast_station(ed.path) if ed.path else None 
            sid = sid or (str(ed.path) if ed.path else "-")

            if sid in by_station and self.on_dup == "keep":
                continue
            by_station[sid] = ed
            seen[sid] = 1

        out.extend(by_station.values())
        return out

    def errors(self) -> List[Tuple[Path, BaseException]]:
        out = [
            (r.path, r.error)  # type: ignore[misc]
            for r in self.results
            if r.error is not None
        ]
        # ---------- NEW ----------
        out.extend(self._errors)
        return out
    
class CBBase:
    r"""
    Lightweight base for EDI collections.

    The class provides a minimal container over multiple
    :class:`EDIFile` objects.  It focuses on indexing by
    station id, fast lookup by path, and simple iteration.
    Subclasses can add project-specific logic, caching, or
    derived computations.

    Parameters
    ----------
    edis : sequence of EDIFile, optional
        Initial items to populate the collection.  When
        omitted the container starts empty.
    index_by : {'station', 'path'}, default ``'station'``
        Key to index items.  With ``'station'`` the
        ``EDIFile.station`` is used, with fallback to file
        path if missing.  With ``'path'`` the absolute path
        string is used as the key.
    on_dup : {'replace', 'keep'}, default ``'replace'``
        Duplicate policy when inserting items that share the
        same key.

    Attributes
    ----------
    items : dict
        Mapping from key to :class:`EDIFile`.
    order : list of str
        Insertion order of keys, which preserves a stable
        iteration sequence.
    meta : dict
        Arbitrary metadata for client code.  This is not used
        internally.

    Methods
    -------
    add(ed)
        Insert one :class:`EDIFile` respecting the duplicate
        policy.
    get(key)
        Return the stored :class:`EDIFile` by its key.
    __iter__()
        Iterate over stored :class:`EDIFile` objects in
        insertion order.
    keys()
        Yield collection keys in insertion order.
    values()
        Yield :class:`EDIFile` objects in insertion order.
    items()
        Yield ``(key, EDIFile)`` pairs in insertion order.

    Notes
    -----
    The base does not perform I/O or parsing.  Use
    :class:`CoreParser` to build items, then feed them into
    the collection.  The class aims to be small and easy to
    extend rather than fully featured.

    Examples
    --------
    Build a collection indexed by station::

        edis = CoreParser().parse(['data/edi'])
        coll = CBBase(edis, index_by='station')
        for key, ed in coll.items():
            print(key, ed.Z.n_freq)

    See Also
    --------
    CoreParser
        Parser that produces :class:`EDIFile` objects to be
        stored in collections.
    EDIFile
        Single-file container stored in the collection.

    References
    ----------
    .. [1] SEG EDI MT/EMAP standard (1987), MTNet.
       https://www.mtnet.info/docs/seg_mt_emap_1987.pdf
    """


    def __init__(
        self,
        items: Optional[Iterable[EDIFile]] = None,
        *,
        verbose: int = 0,
    ) -> None:
        self.verbose = int(verbose)
        self._items: List[EDIFile] = []
        self._index: Dict[str, int] = {}

        if items is not None:
            for ed in items:
                self.add(ed)

    # --------------- core collection protocol ----------------

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[EDIFile]:
        return iter(self._items)

    def __getitem__(self, key: Union[int, str]) -> EDIFile:
        if isinstance(key, int):
            return self._items[key]
        idx = self._index.get(str(key), None)
        if idx is None:
            raise KeyError(f"unknown station: {key!r}")
        return self._items[idx]

    def add(self, ed: EDIFile) -> None:
        sid = getattr(ed, "station", None)
        sid = sid or (str(ed.path) if ed.path else f"#{len(self)}")
        if sid in self._index:
            # replace semantics: keep last one
            i = self._index[sid]
            self._items[i] = ed
            return
        self._index[sid] = len(self._items)
        self._items.append(ed)

    def stations(self) -> List[str]:
        return list(self._index.keys())

    # --------------- ingest and bulk utilities ---------------

    @classmethod
    def load(
        cls,
        sources: SrcType,
        *,
        parser: Optional[CoreParser] = None,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> "CBBase":
        pr = parser or CoreParser(
            recursive=recursive,
            strict=strict,
            on_dup=on_dup,
            verbose=verbose,
        )
        items = pr.parse(sources)
        col = cls(items=items, verbose=verbose)
        errs = pr.errors()
        if errs:
            logger.info("Load completed with %d errors.", len(errs))
        return col

    def map(
        self,
        fn: Callable[[EDIFile], object],
    ) -> List[object]:
        out: List[object] = []
        for ed in self._items:
            out.append(fn(ed))
        return out

    def interpolate(
        self,
        new_freq: Union[np.ndarray, List[float]],
        *,
        kind: str = "slinear",
        bounds_error: bool = True,
        period_buffer: Optional[float] = None,
    ) -> "CBBase":
        out_items: List[EDIFile] = []
        for ed in self._items:
            z2 = ed.interpolate(
                new_freq,
                kind=kind,
                bounds_error=bounds_error,
                period_buffer=period_buffer,
            )
            fresh = EDIFile(ed.path, verbose=ed.verbose)
            fresh.Z = z2
            fresh.Tip = ed.Tip
            # carry over parsed aux sections if any
            for k in (
                "head",
                "info",
                "definemeasurement",
                "mtsect",
                "spectra_sect",
                "spectra_io",
                "spectra",
                "timeseries_sect",
                "timeseries_io",
                "timeseries",
                "other",
                "otherio",
            ):
                sec = getattr(ed, "get_section")(k)  # type: ignore[misc]  # noqa: E501
                if sec is not None:
                    fresh.add_section(k, sec)
            out_items.append(fresh)
        return CBBase(items=out_items, verbose=self.verbose)

    def write(
        self,
        savepath: Pathish,
        *,
        pattern: str = "{station}.edi",
        **kwargs,
    ) -> List[str]:
        out_dir = Path(str(savepath)).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        paths: List[str] = []
        for ed in self._items:
            sid = getattr(ed, "station", None) or "site"
            name = pattern.format(station=sid)
            p = out_dir / name
            s = ed.write(new_edifn=str(p), **kwargs)
            paths.append(s)
        return paths


    def summary(self) -> List[Dict[str, object]]:
        rows: List[Dict[str, object]] = []
        for ed in self._items:
            sid = getattr(ed, "station", None) or "-"
            nf = int(getattr(ed.Z, "n_freq", 0) or 0)
            tip = getattr(ed.Tip, "tipper", None)
            tip_ok = False
            if tip is not None:
                tip_ok = bool(tip.size and not np.all(tip == 0.0))
            rows.append(
                {
                    "station": sid,
                    "path": str(ed.path) if ed.path else "-",
                    "n_freq": nf,
                    "tipper": tip_ok,
                    "spectra": ed.get_section("spectra") is not None,
                    "ts": ed.get_section("timeseries") is not None,
                }
            )
        return rows
    
    @property 
    def items(self):
        """Internal: unified iterator over stored items."""
        return getattr(self, "_items", [])

    def __repr__(self) -> str:  # pragma: no cover
        return (
            f"CBBase(n={len(self)}, stations={self.stations()!r})"
        )

    def __str__(self) -> str:  # pragma: no cover
        lines = ["CBBase"]
        for r in self.summary():
            lines.append(
                f"  {r['station']}: nf={r['n_freq']} "
                f"tip={'Y' if r['tipper'] else 'N'} "
                f"sp={'Y' if r['spectra'] else 'N'} "
                f"ts={'Y' if r['ts'] else 'N'}"
            )
        return "\n".join(lines)
