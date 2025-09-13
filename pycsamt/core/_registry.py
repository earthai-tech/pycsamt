# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, Iterable, List, Optional
from pathlib import Path
import hashlib
import json
import time
import uuid

from .base import CoreObject 

try:  # 3.11+
    import tomllib  # type: ignore[attr-defined]
except Exception:  # pragma: no cover
    tomllib = None  # type: ignore

__all__ = [
    "RegistryError",
    "Record",
    "Manifest",
    "ManifestStore",
    "FileManifestStore",
    "Registry",
    "guess_kind",
]


class RegistryError(RuntimeError):
    r"""
    Error raised for registry and manifest operations.
    
    Use this for signaling inconsistencies, corruption, or
    unsupported formats detected while loading or saving
    manifests or while mutating registry content.
    """

    pass

def _now() -> float:
    r"""
    Return the current POSIX time as ``float`` seconds.
    
    Returns
    -------
    float
        Seconds since the Unix epoch as returned by
        :func:`time.time`.
    
    Notes
    -----
    This is wall-clock time and may jump if the system clock
    changes.  For measuring durations prefer
    :func:`time.monotonic` (not used here to keep values
    serializable and human-interpretable).
    
    Examples
    --------
    >>> ts = _now()
    >>> isinstance(ts, float)
    True
    """

    return float(time.time())


def _uuid() -> str:
    r"""
    Return a random UUID4 hex string without dashes.
    
    Returns
    -------
    str
        32 lowercase hexadecimal characters, suitable as a
        compact record id.
    
    Notes
    -----
    Implemented via :func:`uuid.uuid4().hex`.  The value is not
    cryptographically stable across processes, but is unique with
    very high probability.
    
    Examples
    --------
    >>> rid = _uuid()
    >>> len(rid), rid.isalnum()
    (32, True)
    """

    return uuid.uuid4().hex


def _sha256(path: Path, block: int = 1 << 20) -> str:
    r"""
    Compute a streaming SHA-256 hex digest for a file.
    
    Parameters
    ----------
    path : pathlib.Path
        Path to the file to hash.  The file is opened in binary
        mode and read in blocks.
    block : int, optional
        Block size in bytes.  Defaults to ``1 << 20`` (1 MiB).
        Larger blocks may be faster on fast disks.
    
    Returns
    -------
    str
        Lowercase hexadecimal SHA-256 digest.
    
    Raises
    ------
    FileNotFoundError
        If ``path`` does not exist.
    PermissionError
        If the file cannot be opened.
    OSError
        For other I/O errors during reading.
    
    Notes
    -----
    Reads the file incrementally to keep memory usage bounded.
    
    Examples
    --------
    >>> from pathlib import Path
    >>> digest = _sha256(Path(__file__).parent / "__init__.py")
    >>> len(digest)
    64
    
    See Also
    --------
    hashlib.sha256
    """

    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def guess_kind(obj: Any) -> str:
    r"""
    Best-effort classification for common MT/EM artefacts.
    
    The function inspects ``obj.__class__.__module__`` and class
    name to assign a coarse category usable by higher-level
    tools.  It recognizes common terms such as ``"seg"``,
    ``"edi"``, ``"zonge"``, ``"avg"``, and ``"jones"``.
    
    Parameters
    ----------
    obj : Any
        Any Python object.  Only its class name and module path
        are inspected.  Iterables are handled specially.
    
    Returns
    -------
    str
        One of the following labels:
    
        * ``"edi"``       — an EDI item.
        * ``"edi_col"``   — an EDI collection.
        * ``"avg"``       — a Zonge AVG artefact.
        * ``"j"``         — a Jones J-file object.
        * ``"j_col"``     — a collection of J-files.
        * ``"list"``      — a list or tuple of items.
        * ``"meta"``      — a fallback for other metadata.
    
    Notes
    -----
    This is heuristic and string-based.  It is intentionally
    liberal to remain dependency-free.  Extend or override at the
    call site if stricter typing is needed.
    
    Examples
    --------
    >>> class Dummy: ...
    >>> guess_kind(Dummy())
    'meta'
    >>> guess_kind([])
    'list'
    
    See Also
    --------
    Record
    Manifest
    """

    try:
        mod = obj.__class__.__module__.lower()
        cls = obj.__class__.__name__.lower()
    except Exception:
        mod, cls = "", ""
    if ("seg" in mod) and ("edi" in (mod + cls)):
        if ("collection" in mod) or ("collection" in cls):
            return "edi_col"
        return "edi"
    if ("zonge" in mod) or ("avg" in (mod + cls)):
        return "avg"
    if ("jones" in mod) or (cls in {"jfile", "jcollection"}):
        if cls == "jcollection":
            return "j_col"
        return "j"
    if isinstance(obj, (list, tuple)):
        return "list"
    return "meta"


@dataclass
class Record(CoreObject):
    r"""
    A single registry item with identifiers and metadata.
    
    Parameters
    ----------
    rid : str
        Unique record id.  Often a value from :func:`_uuid`.
    kind : str
        Logical category.  See :func:`guess_kind` for common
        labels (e.g. ``"edi"``, ``"avg"``).
    path : str or None, optional
        Path (relative to :attr:`Manifest.root` or absolute) to
        the primary file associated with this record.
    fmt : str or None, optional
        Short format tag (e.g. ``"edi"`` or ``"json"``).  Free
        form and optional.
    dataid : str or None, optional
        External data id or survey id if applicable.
    station_id : str or None, optional
        Station or site identifier for field data.
    tags : list of str, optional
        Free-form labels for grouping or filtering.
    meta : dict, optional
        Arbitrary JSON-serializable metadata.
    checksum : str or None, optional
        SHA-256 hex digest of the file at :attr:`path`.  Use
        :func:`_sha256` to compute it.  Optional.
    created : float, optional
        Creation timestamp in seconds since epoch.  Defaults to
        :func:`_now`.
    updated : float, optional
        Last update timestamp in seconds since epoch.  Defaults
        to :func:`_now`.
    
    Attributes
    ----------
    rid : str
    kind : str
    path : str or None
    fmt : str or None
    dataid : str or None
    station_id : str or None
    tags : list of str
    meta : dict
    checksum : str or None
    created : float
    updated : float
    
    Notes
    -----
    The class is a :mod:`dataclasses` dataclass and is designed
    to round-trip with :meth:`to_dict` and :meth:`from_dict`.
    
    Examples
    --------
    Create a record with a checksum::
    
        from pathlib import Path
        p = Path("data/site001.edi")
        r = Record(rid=_uuid(), kind="edi", path=str(p),
                   checksum=_sha256(p))
    
    See Also
    --------
    Manifest
    guess_kind
    """

    rid: str
    kind: str
    path: Optional[str] = None
    fmt: Optional[str] = None
    dataid: Optional[str] = None
    station_id: Optional[str] = None
    tags: List[str] = field(default_factory=list)
    meta: Dict[str, Any] = field(default_factory=dict)
    checksum: Optional[str] = None
    created: float = field(default_factory=_now)
    updated: float = field(default_factory=_now)

    def touch(self) -> None:
        r"""
        Bump :attr:`updated` to the current time.
        
        Notes
        -----
        Uses :func:`_now`.  Does not mutate any other fields.
        """
        self.updated = _now()

    def to_dict(self) -> Dict[str, Any]:
        r"""
        Return a JSON-ready mapping for this record.
        
        Returns
        -------
        dict
            A plain ``dict`` produced via :func:`dataclasses.asdict`.
        """
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Record":
        r"""
        Create a :class:`Record` from a mapping.
        
        Parameters
        ----------
        d : dict
            Mapping with the same keys as :meth:`to_dict`.
        
        Returns
        -------
        Record
            New instance populated from the mapping.
        """
        return cls(**d)


@dataclass
class Manifest(CoreObject):
    r"""
    Container mapping record ids to :class:`Record` objects.
    
    Parameters
    ----------
    root : str
        Base directory for paths stored in records.  Callers may
        join this with :attr:`Record.path`.
    version : int, optional
        Schema or persistence version.  Defaults to ``1``.
    records : dict[str, Record], optional
        Mapping of record id to :class:`Record`.  Defaults to an
        empty mapping.
    
    Attributes
    ----------
    root : str
    version : int
    records : dict[str, Record]
    
    Notes
    -----
    The manifest is serializable through :meth:`to_dict` and
    :meth:`from_dict`.  Storage format (JSON, TOML, etc.) is
    left to :class:`ManifestStore` implementations.
    
    Examples
    --------
    Build a manifest programmatically::
    
        man = Manifest(root="data")
        rec = Record(rid=_uuid(), kind="edi", path="a.edi")
        man.records[rec.rid] = rec
    
    See Also
    --------
    Record
    ManifestStore
    json
    """
    root: str
    version: int = 1
    records: Dict[str, Record] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        r"""
        Return a nested mapping suitable for JSON or TOML.
        
        A ``dict`` with ``root``, ``version``, and a nested
        ``records`` mapping of ids to plain dicts.
        """

        return {
            "root": self.root,
            "version": self.version,
            "records": {k: r.to_dict() for k, r in self.records.items()},
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Manifest":
        r"""
        Create a :class:`Manifest` from a nested mapping.
        
        Parameters
        ----------
        d : dict
            Mapping previously produced by :meth:`to_dict`.
        
        Returns
        -------
        Manifest
            New manifest with :class:`Record` instances rebuilt.
        """

        recs = {
            k: Record.from_dict(v) 
            for k, v in d.get("records", {}).items()
        }
        return cls(
            root=d.get("root", "."), version=int(
                d.get("version", 1)), records=recs)


class ManifestStore(CoreObject):
    r"""
    Interface for loading and saving :class:`Manifest` objects.
    
    Subclasses provide concrete persistence (e.g. JSON on disk,
    TOML, or a database).  The interface is intentionally small.
    
    Notes
    -----
    Implementations should validate the manifest version and any
    schema expectations relevant to their storage.
    
    Examples
    --------
    A minimal JSON store sketch::
    
        class JsonStore(ManifestStore):
            def load(self, path: Path) -> Manifest:
                data = json.loads(path.read_text())
                return Manifest.from_dict(data)
    
            def save(self, man: Manifest, path: Path) -> None:
                path.write_text(json.dumps(man.to_dict(), indent=2))
    
    See Also
    --------
    Manifest
    Record
    json
    tomllib
    """
    def load(self, path: Path) -> Manifest:  # noqa: D401
        raise NotImplementedError

    def save(self, man: Manifest, path: Path) -> None:  # noqa: D401
        raise NotImplementedError


class FileManifestStore(ManifestStore):
    r"""
    Disk-backed :class:`Manifest` persistence utilities.
    
    This store loads and saves manifests from files on disk.
    Reading supports JSON by default, and can parse TOML when
    :mod:`tomllib` is available. Saving is JSON-only in this
    implementation.
    
    Behavior
    --------
    * If the given ``path`` does not exist on load, a new empty
      :class:`Manifest` is returned with ``root`` set to the
      parent directory of ``path``.
    * Files are read and written with ``utf-8`` encoding.
    * TOML saving is not implemented and raises
      :class:`RegistryError`.
    
    Notes
    -----
    This class is intentionally minimal and does not implement
    atomic writes or file locking. For concurrent writers use a
    safer store or wrap calls with your own locking.
    
    Examples
    --------
    >>> from pathlib import Path
    >>> from pycsamt.core._registry import (  # doctest: +SKIP
    ...     FileManifestStore, Manifest)
    >>> store = FileManifestStore()            # doctest: +SKIP
    >>> mpath = Path("manifest.json")          # doctest: +SKIP
    >>> man = store.load(mpath)                # doctest: +SKIP
    >>> isinstance(man, Manifest)              # doctest: +SKIP
    True
    >>> store.save(man, mpath)                 # doctest: +SKIP
    
    See Also
    --------
    Manifest
    ManifestStore
    json
    tomllib
    
    References
    ----------
    .. [1] Python Standard Library. *json* module.
    .. [2] Python 3.11+. *tomllib* — TOML parser.
    """

    def load(self, path: Path) -> Manifest:
        r"""
        Load a :class:`Manifest` from ``path``.
        
        Parameters
        ----------
        path : pathlib.Path
            File to read. If it does not exist, return a new empty
            manifest with ``root`` set to ``path.parent``.
        
        Returns
        -------
        Manifest
            Deserialized manifest.
        
        Raises
        ------
        json.JSONDecodeError
            If JSON parsing fails.
        OSError
            For I/O errors.
        
        Notes
        -----
        When the suffix is ``.toml`` and :mod:`tomllib` is missing,
        the content is parsed as JSON, which will likely fail. Install
        ``tomllib`` (Python 3.11+) to read TOML files.
        """

        if not path.exists():
            return Manifest(root=str(path.parent))
        txt = path.read_text(encoding="utf-8")
        if path.suffix.lower() == ".toml" and tomllib is not None:
            d = tomllib.loads(txt)
        else:
            d = json.loads(txt)
        return Manifest.from_dict(d)

    def save(self, man: Manifest, path: Path) -> None:
        r"""
        Serialize ``man`` to JSON at ``path`` using ``utf-8``.
        
        Parameters
        ----------
        man : Manifest
            Manifest to write.
        path : pathlib.Path
            Destination file. The suffix determines support.
        
        Raises
        ------
        RegistryError
            If a TOML file is requested (``.toml`` suffix).
        OSError
            For write or permission errors.
        
        Notes
        -----
        This method writes pretty-printed JSON (``indent=2``). TOML
        saving is not implemented here by design.
        """

        d = man.to_dict()
        if path.suffix.lower() == ".toml":
            raise RegistryError("TOML save not supported here")
        path.write_text(json.dumps(d, indent=2), encoding="utf-8")


class Registry(CoreObject):
    r"""
    High-level helper to manage a file-backed registry.
    
    A :class:`Registry` binds a filesystem ``root`` with a
    manifest file and provides convenience methods to add files
    or objects, query records, and persist changes.
    
    Parameters
    ----------
    root : path-like
        Base directory for data and the manifest file.
    manifest_name : str, optional
        Manifest filename within ``root``. Defaults to
        ``"manifest.json"``.
    store : ManifestStore or None, optional
        Custom persistence backend. Defaults to
        :class:`FileManifestStore`.
    
    Attributes
    ----------
    root : pathlib.Path
        Absolute path to the registry root.
    manifest_path : pathlib.Path
        Full path of the manifest file.
    store : ManifestStore
        Persistence backend in use.
    manifest : Manifest
        In-memory manifest loaded at construction.
    
    Notes
    -----
    On initialization, the manifest is loaded (or created if
    missing) and its ``root`` is set to the absolute ``root``
    path. All relative paths passed later are resolved under
    ``root``.
    
    Examples
    --------
    Create a registry and register an EDI file::
    
        >>> from pycsamt.core._registry import Registry  # doctest: +SKIP
        >>> reg = Registry("data")                       # doctest: +SKIP
        >>> rec = reg.add_path("site001.edi", kind="edi",
        ...                   fmt="edi")                 # doctest: +SKIP
        >>> reg.get(rec.rid).kind                        # doctest: +SKIP
        'edi'
    
    Register a Python object and tag it::
    
        >>> class Dummy:                                 # doctest: +SKIP
        ...     station = "S01"
        ...     station_id = 1
        >>> d = Dummy()                                  # doctest: +SKIP
        >>> rec = reg.add_obj(d, tags=["demo"])          # doctest: +SKIP
        >>> "demo" in rec.tags                           # doctest: +SKIP
        True
    
    Query and update metadata::
    
        >>> regs = reg.list(kind="edi")                  # doctest: +SKIP
        >>> found = reg.find(kind="edi")                 # doctest: +SKIP
        >>> _ = reg.update_meta(rec.rid, site="S001")    # doctest: +SKIP
    
    See Also
    --------
    Record
    Manifest
    ManifestStore
    FileManifestStore
    guess_kind
    
    References
    ----------
    .. [1] D. Wight (1991). *SEG MT/EMAP EDI Standard*.
       Society of Exploration Geophysicists.
    .. [2] Python Standard Library. *pathlib* module.
    """

    def __init__(
        self,
        root: Path | str,
        *,
        manifest_name: str = "manifest.json",
        store: Optional[ManifestStore] = None,
    ) -> None:
        self.root = Path(root)
        self.manifest_path = self.root / manifest_name
        self.store = store or FileManifestStore()
        self.manifest = self.store.load(self.manifest_path)
        self.manifest.root = str(self.root)

    #  basic ops 
    def save(self) -> None:
        r"""
        Persist the in-memory manifest to :attr:`manifest_path`.
        
        Notes
        -----
        Delegates to :attr:`store`. Use this after batch updates if
        you have modified :attr:`manifest` directly.
        """

        self.store.save(self.manifest, self.manifest_path)

    def add_path(
        self,
        p: Path | str,
        *,
        kind: Optional[str] = None,
        fmt: Optional[str] = None,
        dataid: Optional[str] = None,
        station_id: Optional[str] = None,
        tags: Optional[Iterable[str]] = None,
        meta: Optional[Dict[str, Any]] = None,
        with_hash: bool = True,
    ) -> Record:
        r"""
        Add a file-based record into the registry.
        
        Parameters
        ----------
        p : path-like
            File path. If relative, it is resolved under
            :attr:`root`.
        kind : str or None, optional
            Logical kind for the record. Defaults to ``"meta"`` when
            not provided.
        fmt : str or None, optional
            Short format tag, such as ``"edi"``.
        dataid : str or None, optional
            External dataset id or survey id.
        station_id : str or int or None, optional
            Station identifier. Converted to string if provided.
        tags : Iterable[str] or None, optional
            User-defined labels.
        meta : dict or None, optional
            Extra JSON-serializable metadata.
        with_hash : bool, optional
            If ``True`` and the file exists, compute SHA-256 and set
            :attr:`Record.checksum`. Defaults to ``True``.
        
        Returns
        -------
        Record
            The created record.
        
        Notes
        -----
        When the file does not exist, the checksum is left unset.
        The manifest is saved immediately after insertion.
        
        Examples
        --------
        >>> from pycsamt.core._registry import Registry  # doctest: +SKIP
        >>> reg = Registry("data")                       # doctest: +SKIP
        >>> rec = reg.add_path("a.edi", kind="edi")     # doctest: +SKIP
        """

        pth = Path(p)
        if not pth.is_absolute():
            pth = (self.root / pth).resolve()
        rid = _uuid()
        k = kind or "meta"
        cs = _sha256(pth) if with_hash and pth.exists() else None
        rec = Record(
            rid=rid,
            kind=k,
            path=str(pth),
            fmt=fmt,
            dataid=dataid,
            station_id=str(station_id) if station_id is not None else None,
            tags=list(tags or ()),
            meta=dict(meta or {}),
            checksum=cs,
        )
        self.manifest.records[rid] = rec
        self.save()
        return rec

    def add_obj(
        self,
        obj: Any,
        *,
        rid: Optional[str] = None,
        kind: Optional[str] = None,
        tags: Optional[Iterable[str]] = None,
        meta: Optional[Dict[str, Any]] = None,
        path: Optional[str | Path] = None,
    ) -> Record:
        r"""
        Add an object-backed record using heuristics for ``kind``.
        
        Parameters
        ----------
        obj : Any
            Object to classify. :func:`guess_kind` is used when
            ``kind`` is not provided.
        rid : str or None, optional
            Optional explicit record id. Defaults to a random id.
        kind : str or None, optional
            Override classification if given.
        tags : Iterable[str] or None, optional
            User-defined labels.
        meta : dict or None, optional
            Extra JSON-serializable metadata.
        path : str or pathlib.Path or None, optional
            Optional path reference to associate with the object.
        
        Returns
        -------
        Record
            The created record.
        
        Notes
        -----
        If the object exposes attributes ``station`` or
        ``station_id``, they are copied into the record when
        available. The manifest is saved immediately.
        
        Examples
        --------
        >>> from pycsamt.core._registry import Registry  # doctest: +SKIP
        >>> reg = Registry("data")                       # doctest: +SKIP
        >>> class Obj:                                   # doctest: +SKIP
        ...     station = "S02"; station_id = 2
        >>> rec = reg.add_obj(Obj(), tags=["x"])        # doctest: +SKIP
        """

        rid = rid or _uuid()
        k = kind or guess_kind(obj)
        rec = Record(
            rid=rid,
            kind=k,
            path=str(path) if path else None,
            fmt=None,
            dataid=getattr(obj, "station", None),
            station_id=str(getattr(obj, "station_id", "")) or None,
            tags=list(tags or ()),
            meta=dict(meta or {}),
        )
        self.manifest.records[rid] = rec
        self.save()
        return rec

    def get(self, rid: str) -> Record:
        r"""
        Return the :class:`Record` for ``rid`` or raise an error.
        """

        try:
            return self.manifest.records[rid]
        except KeyError as exc:
            raise RegistryError(f"unknown rid: {rid}") from exc

    def list(self, *, kind: Optional[str] = None) -> List[Record]:
        r"""
        List records, optionally filtering by ``kind``.
        
        Parameters
        ----------
        kind : str or None, optional
            If provided, only records whose ``kind`` matches are
            returned.
        
        Returns
        -------
        list of Record
            Records in insertion order as stored in the manifest.
        """

        vals = list(self.manifest.records.values())
        if kind:
            return [r for r in vals if r.kind == kind]
        return vals

    def find(
        self,
        *,
        tag: Optional[str] = None,
        kind: Optional[str] = None,
        dataid: Optional[str] = None,
    ) -> List[Record]:
        r"""
        Filter records by any combination of tag, kind, or dataid.
        
        Parameters
        ----------
        tag : str or None, optional
            Require that this tag is present in the record tags.
        kind : str or None, optional
            Require that the record kind matches this value.
        dataid : str or None, optional
            Require that the record data id matches this value.
        
        Returns
        -------
        list of Record
            Records satisfying all provided predicates.
        
        Examples
        --------
        >>> from pycsamt.core._registry import Registry  # doctest: +SKIP
        >>> reg = Registry("data")                       # doctest: +SKIP
        >>> reg.find(kind="edi")                         # doctest: +SKIP
        """

        out: List[Record] = []
        for r in self.manifest.records.values():
            if kind and r.kind != kind:
                continue
            if tag and tag not in (r.tags or ()): 
                continue
            if dataid and r.dataid != dataid:
                continue
            out.append(r)
        return out

    def update_meta(
        self,
        rid: str,
        **fields: Any,
    ) -> Record:
        r"""
        Update ``meta`` fields of the record with id ``rid``.
        
        Parameters
        ----------
        rid : str
            Record identifier.
        **fields : Any
            Key-value pairs merged into :attr:`Record.meta`.
        
        Returns
        -------
        Record
            The updated record.
        
        Notes
        -----
        Updates the ``updated`` timestamp and saves the manifest.
        """

        r = self.get(rid)
        r.meta.update(fields)
        r.touch()
        self.save()
        return r

    def remove(self, rid: str) -> None:
        r"""
        Remove the record with id ``rid`` and persist the manifest.
        
        Parameters
        ----------
        rid : str
            Record identifier.
        
        Notes
        -----
        The removal is silent when the id is unknown. The manifest
        is saved after the operation.
        """

        self.manifest.records.pop(rid, None)
        self.save()
