# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Callable,
    Union,
)

import numpy as np

from ..log.logger import get_logger
from .config import (
    RE_BLANK,
    RE_COMMENT,
    RE_DATATYPE_UNITS,
    RE_INFO,
    RE_NPOINTS,
    RE_STATION,
)
from .j import JFile
from .validation import IsJ

logger = get_logger(__name__)

__all__ = ["JParseMixin", "JCoreParser", "JCBBase"]

Pathish = Union[str, Path]
SrcType = Union[Pathish, Sequence[Pathish]]


class JParseMixin:
    r"""
    Lightweight helpers for scanning Jones J-format text.

    The mixin implements tolerant, file-level utilities used by
    higher-level parsers and collections.  It focuses on small,
    dependency-free pieces such as path coercion, candidate file
    discovery, banner and header probing, and quick content
    reads.

    The mixin does **not** keep state.  Methods are small and
    side-effect free so they can be reused by different classes
    (e.g., :class:`JCoreParser`, :class:`JCBBase`,
    :class:`~pycsamt.jones.collection.JCollection`).

    Notes
    -----
    Utilities are designed to be permissive.  They handle odd
    encodings, mixed line endings, and common filename patterns.
    They also accept both :class:`pathlib.Path` and strings.

    The intent is to keep the heavy parsing in dedicated
    readers, while providing robust file discovery and quick
    checks here.

    Examples
    --------
    Discover candidate files in a folder:

    >>> m = JParseMixin()
    >>> root = "data/j"
    >>> list(m._iter_j_files([root]))[:1]  # doctest: +ELLIPSIS
    [PosixPath('...kb0-s001.txt')]

    Coerce user inputs to :class:`Path`:

    >>> p = m._as_path("data/j/kb0-s001.txt")
    >>> p.exists()
    True

    See Also
    --------
    JCoreParser :
        Adds content probing (station, dtype, counts).
    JCBBase :
        State-holding base for collections.
    pycsamt.jones.validation.IsJ :
        File-level validator for J candidates.
    pycsamt.jones.j.JFile :
        High-level reader that builds Z/Tipper/R blocks.

    References
    ----------
    .. [JParseMixin-1] A. G. Jones (1994). *J-format v2.0*. MTNet notes.
    """

    J_SUFFIXES = {".j", ".jones", ".txt", ".dat"}

    def _as_path(self, p: Pathish) -> Path:
        return Path(str(p)).expanduser().resolve()

    def _is_j_path(self, p: Path) -> bool:
        return p.is_file() and (p.suffix.lower() in self.J_SUFFIXES)

    def _iter_paths(self, src: SrcType) -> Iterator[Path]:
        if isinstance(src, (str, Path)):
            yield self._as_path(src)  # type: ignore[arg-type]
            return
        for s in src:
            yield self._as_path(s)

    def _push_error(self, src: Pathish, msg: str) -> None:
        store = getattr(self, "_errors", None)
        if store is None:
            self._errors = []  # type: ignore[attr-defined]
            store = self._errors
        err = FileNotFoundError(msg)
        try:
            # src may be an unmatched glob pattern; resolving it can raise
            # on Windows (GetFinalPathNameByHandle rejects "*"/"?"/"[]").
            p = self._as_path(src)
        except OSError:
            p = Path(str(src))
        store.append((p, err))

    def _iter_j_files(
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
                if self._is_j_path(p):
                    yield p
                else:
                    self._push_error(src, f"Not a J file: {src}")
                continue

            # directory -> (r)glob, filtered to J_SUFFIXES (.j/.jones/.txt/.dat)
            if p.exists() and p.is_dir():
                it = (
                    p.rglob("*")
                    if getattr(self, "recursive", True)
                    else p.glob("*")
                )
                any_yielded = False
                # Directory enumeration order (os.scandir/readdir) is
                # filesystem-defined, not insertion order -- sort so
                # discovery is reproducible across platforms instead of
                # an accident of the OS.
                for m in sorted(it):
                    if self._is_j_path(m):
                        any_yielded = True
                        yield m
                if not any_yielded:
                    self._push_error(src, f"No .j under: {src}")
                continue

            # glob pattern (absolute or relative)
            if _is_glob(str(p)):
                if p.is_absolute():
                    parent, name = p.parent, p.name
                    it = (
                        parent.rglob(name)
                        if "**" in name
                        else parent.glob(name)
                    )
                else:
                    pat = str(p)
                    it = base.rglob(pat) if "**" in pat else base.glob(pat)
                any_yielded = False
                for m in sorted(it):
                    if self._is_j_path(m):
                        any_yielded = True
                        yield m
                if not any_yielded:
                    self._push_error(src, f"No match: {src}")
                continue

            # anything else
            self._push_error(src, f"Not found: {src}")

    def _fast_station(self, p: Path) -> str | None:
        """Quick scan for the first station line."""
        try:
            with p.open("r", encoding="utf-8") as f:
                for raw in f:
                    s = raw.rstrip("\n")
                    if not s or RE_BLANK.match(s):
                        continue
                    if RE_COMMENT.match(s) or RE_INFO.match(s):
                        continue
                    m = RE_STATION.match(s)
                    if m:
                        return m.group("station").strip().upper()
        except Exception:
            return None
        return None

    def is_j_like(self, src: Pathish, *, deep: bool = True) -> bool:
        """
        Light heuristic to decide if ``src`` looks like a Jones file.

        Parameters
        ----------
        src : path-like
            File path.
        deep : bool, default=True
            If ``False``, only extension + existence is checked.

        Returns
        -------
        bool
            ``True`` if the file looks like J-format.
        """
        p = self._as_path(src)
        if not self._is_j_path(p):
            return False
        if not deep:
            return True

        try:
            with p.open("r", encoding="utf-8", errors="replace") as f:
                lines = []
                # limit reads to keep it snappy on large files
                for _ in range(300):
                    ln = f.readline()
                    if not ln:
                        break
                    lines.append(ln.rstrip("\n"))
        except Exception:
            return False

        found_banner = False
        found_info = False
        found_triple = False

        for i, s in enumerate(lines):
            if not s or RE_BLANK.match(s):
                continue

            # banner (very permissive)
            if s.lstrip().upper().startswith("#WRITTEN BY"):
                found_banner = True

            # any info line (>KEY=VAL)
            if RE_INFO.match(s):
                found_info = True

            # station + dtype/count nearby
            m = RE_STATION.match(s)
            if not m:
                continue

            # scan a few lines ahead to find dtype/count
            j = i + 1
            n = min(len(lines), i + 10)
            saw_dtype = False
            saw_count = False
            while j < n:
                t = lines[j]
                if RE_BLANK.match(t) or RE_COMMENT.match(t):
                    j += 1
                    continue
                if RE_STATION.match(t):
                    break
                if RE_DATATYPE_UNITS.match(t):
                    saw_dtype = True
                elif RE_NPOINTS.match(t):
                    saw_count = True
                # short-circuit if both seen in small window
                if saw_dtype and saw_count:
                    found_triple = True
                    break
                j += 1
            if found_triple:
                break
        # Require at least a header triple; banner/info are helpful hints
        return bool(found_triple or (found_banner and found_info))

    def is_j_file(self, src: Pathish, *, deep: bool = True) -> bool:
        """
        Spec-level check via :class:`IsJ`. Returns ``True`` if
        the file satisfies the structural requirements; ``False``
        otherwise. Never raises.
        """
        try:
            return IsJ._assert_j(src, deep=deep)
        except Exception:
            return False

    def is_j_candidate(self, src: Pathish, *, deep: bool = True) -> bool:
        """
        Alias of :meth:`is_j_like`. Some call-sites prefer this name.
        """
        return self.is_j_like(src, deep=deep)


@dataclass
class _JParseResult:
    path: Path
    jf: JFile | None
    error: BaseException | None


class JCoreParser(JParseMixin):
    r"""
    Core scanner for J files that extracts light metadata.

    Extends :class:`JParseMixin` with small, read-only probes
    that inspect the top header and first data head triple.
    This class can reveal the station code, presence of info
    keys, encountered data kinds, and the declared row count of
    the first block.  It avoids loading the whole file.

    The implementation aims to be fast and safe for directory
    crawls, where many files are filtered before full parse.

    Attributes
    ----------
    encoding : str, default ``'utf-8'``
        Encoding used when reading text.  Implementations may
        choose a tolerant variant such as ``'utf-8-sig'`` with
        ``errors='replace'``.

    Methods
    -------
    _read_text(path)
        Return the file content as a single text string.
    _scan_one(path)
        Return a small dict of hints (e.g., ``station``,
        ``kinds``, ``nrows``) obtained without full parsing.
    _looks_like_j(path)
        Heuristic check that a path is a J candidate.
    _iter_info_keys(text)
        Yield info keys (``>KEY=VALUE``) from the header area.

    Notes
    -----
    The scanner is structural.  It does not validate numeric
    content or cross-block consistency.  It is designed to be
    called many times in batch workflows.

    Examples
    --------
    Quick peek of a single file:

    >>> p = "data/j/kb0-s001.txt"
    >>> core = JCoreParser()
    >>> hints = core._scan_one(p)
    >>> "station" in hints and "kinds" in hints
    True

    Filter a folder to J candidates:

    >>> paths = list(core._iter_j_files(["data/j"]))
    >>> all(core._looks_like_j(p) for p in paths)
    True

    See Also
    --------
    JParseMixin :
        Path and discovery helpers used by this scanner.
    JCBBase :
        Collection base that can cache these hints.
    pycsamt.jones.heads.Heads :
        Full header reader (Info + Head + Banner).
    pycsamt.jones.blocks.JBlocks :
        Full data block reader.

    References
    ----------
    .. [JCoreParser-1] A. G. Jones (1994). *J-format v2.0*. MTNet notes.
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

        self.results: list[_JParseResult] = []
        self._errors: list[tuple[Path, BaseException]] = []

    def _read_one(self, p: Path) -> _JParseResult:
        try:
            IsJ._assert_j(p, deep=False)
            jf = JFile.from_file(p, verbose=self.verbose)
            return _JParseResult(path=p, jf=jf, error=None)
        except BaseException as exc:  # noqa: BLE001
            if self.strict:
                raise
            logger.debug("Skip %s: %s", p, exc)
            return _JParseResult(path=p, jf=None, error=exc)

    def parse(self, sources: SrcType) -> list[JFile]:
        self.results.clear()
        self._errors.clear()

        jfs: list[JFile] = []
        src_list = (
            sources
            if isinstance(sources, list)
            else list(self._iter_paths(sources))
        )
        for p in self._iter_j_files(src_list):
            res = self._read_one(p)
            self.results.append(res)
            if res.jf is None:
                continue
            jfs.append(res.jf)

        by_key: dict[str, JFile] = {}
        for jf in jfs:
            sid = jf.site or self._fast_station(jf.path)  # type: ignore[arg-type]  # noqa: E501
            sid = sid or (str(jf.path) if jf.path else "-")
            if sid in by_key and self.on_dup == "keep":
                continue
            by_key[sid] = jf

        return list(by_key.values())

    def errors(self) -> list[tuple[Path, BaseException]]:
        out = [(r.path, r.error) for r in self.results if r.error]
        out.extend(self._errors)
        return out


class JCBBase:
    r"""
    Minimal stateful base for collections of J files.

    The class manages a list of parsed items (usually
    :class:`~pycsamt.jones.j.JFile`), keeps light metadata for
    indexing and filtering, and offers a stable surface for
    higher-level collection types.  It is intentionally small so
    projects can extend it without tight coupling.

    Parameters
    ----------
    verbose : int, default 0
        Verbosity for logs and warnings.  Implementations are
        free to ignore it or to pass it to child readers.

    Attributes
    ----------
    items : list
        Stored per-file objects.  The class does not enforce a
        specific type, but ``JFile`` is the natural choice.
    n : int
        Number of stored items.
    paths : list of Path
        Convenience view of the underlying file paths (when
        available on items).

    Methods
    -------
    add(obj)
        Append an item.  Returns the item for chaining.
    parse(paths)
        Build and add items from paths.  Uses the light probes
        from :class:`JCoreParser` before full parse.
    filter(**kws)
        Optional helper that returns a new instance filtered by
        station, component family, or other metadata.
    __len__(), __iter__()
        Python container protocol for lists of items.

    Notes
    -----
    Subclasses are encouraged to keep the API stable while
    augmenting with domain-specific helpers, like saving index
    CSV files or computing per-site quality metrics.

    The base class uses tolerant discovery, so it can be used on
    mixed folders where some files are not Jones J-format.

    Examples
    --------
    Build a collection from a folder:

    >>> col = JCBBase(verbose=0)
    >>> items = col.parse(["data/j"])
    >>> len(items) == col.n
    True

    Iterate over sites and write quick summaries:

    >>> for it in col:
    ...     print(getattr(it, "site", None))  # doctest: +SKIP

    See Also
    --------
    pycsamt.jones.collection.JCollection :
        Higher-level collection with convenience features.
    pycsamt.jones.j.JFile :
        High-level reader used as per-item object.
    JCoreParser :
        Provides the fast scanning used during ``parse``.

    References
    ----------
    .. [JCBBase-1] A. G. Jones (1994). *J-format v2.0*. MTNet notes.
    """

    def __init__(
        self,
        items: Iterable[JFile] | None = None,
        *,
        verbose: int = 0,
    ) -> None:
        self.verbose = int(verbose)
        self._items: list[JFile] = []
        self._index: dict[str, int] = {}

        if items is not None:
            for jf in items:
                self.add(jf)

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self) -> Iterator[JFile]:
        return iter(self._items)

    def __getitem__(self, key: int | str) -> JFile:
        if isinstance(key, int):
            return self._items[key]
        idx = self._index.get(str(key), None)
        if idx is None:
            raise KeyError(f"unknown station: {key!r}")
        return self._items[idx]

    def add(self, jf: JFile) -> None:
        sid = jf.site or (str(jf.path) if jf.path else f"#{len(self)}")
        if sid in self._index:
            i = self._index[sid]
            self._items[i] = jf
            return
        self._index[sid] = len(self._items)
        self._items.append(jf)

    def stations(self) -> list[str]:
        return list(self._index.keys())

    @classmethod
    def load(
        cls,
        sources: SrcType,
        *,
        parser: JCoreParser | None = None,
        recursive: bool = True,
        strict: bool = False,
        on_dup: str = "replace",
        verbose: int = 0,
    ) -> JCBBase:
        pr = parser or JCoreParser(
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

    def map(self, fn: Callable[[JFile], object]) -> list[object]:
        out: list[object] = []
        for jf in self._items:
            out.append(fn(jf))
        return out

    def write(
        self,
        savepath: Pathish,
        *,
        pattern: str = "{station}.j",
        **kwargs,
    ) -> list[str]:
        out_dir = Path(str(savepath)).expanduser().resolve()
        out_dir.mkdir(parents=True, exist_ok=True)

        paths: list[str] = []
        for jf in self._items:
            sid = jf.site or "site"
            name = pattern.format(station=sid)
            p = out_dir / name
            s = jf.write(new_jfn=str(p), **kwargs)
            paths.append(s)
        return paths

    def summary(self) -> list[dict[str, object]]:
        rows: list[dict[str, object]] = []
        for jf in self._items:
            sid = jf.site or "-"
            f = jf.freq
            nf = int(getattr(f, "size", len(f))) if f is not None else 0
            tip_ok = False
            if jf.Tip is not None:
                tarr = getattr(jf.Tip, "tipper_array", None)
                if tarr is not None:
                    a = np.asarray(tarr)
                    tip_ok = bool(a.size and not np.all(a == 0.0))
            rows.append(
                {
                    "station": sid,
                    "path": str(jf.path) if jf.path else "-",
                    "n_freq": nf,
                    "tipper": tip_ok,
                    "has_r": bool(jf.Res is not None),
                    "has_z": bool(jf.Z is not None),
                }
            )
        return rows

    @property
    def items(self):
        """Internal: unified iterator over stored items."""
        return getattr(self, "_items", [])

    def __repr__(self) -> str:  # pragma: no cover
        return f"JCBBase(n={len(self)}, stations={self.stations()!r})"

    def __str__(self) -> str:  # pragma: no cover
        lines = ["JCBBase"]
        for r in self.summary():
            lines.append(
                "  {s}: nf={n} tip={t} R={hr} Z={hz}".format(
                    s=r["station"],
                    n=r["n_freq"],
                    t="Y" if r["tipper"] else "N",
                    hr="Y" if r["has_r"] else "N",
                    hz="Y" if r["has_z"] else "N",
                )
            )
        return "\n".join(lines)
