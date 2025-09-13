# -*- coding: utf-8 -*- 
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0

from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Tuple
from pathlib import Path
import json

import numpy as np

from ._registry  import Registry, guess_kind 
from .base import TFBundle, CoreObject
from .base import to_edi

__all__ = [
    "Packer",
    "register_packer",
    "get_packer",
    "list_packers",
    "pack_to_file",
    "unpack_from_file",
    "RegistryAPI",
]


Packer = Tuple[
    Callable[[Any], Dict[str, Any]],
    Callable[[Dict[str, Any]], Any],
]

_PACKERS: Dict[str, Packer] = {}


def register_packer(kind: str, packer: Packer) -> None:
    r"""
    Register a serializer/deserializer pair for a given ``kind``.
    
    Parameters
    ----------
    kind : str
        Case-insensitive key that identifies the packer.
    packer : tuple(callable, callable)
        Pair ``(pack, unpack)``. The first callable takes an
        object and returns a ``dict``; the second takes that
        mapping and returns the reconstructed object.
    
    Raises
    ------
    ValueError
        If ``kind`` is empty or not a ``str``.
    
    Notes
    -----
    Registration is global to the current process via an
    internal dictionary. Re-registering the same ``kind`` will
    overwrite the previous pair.
    
    Examples
    --------
    >>> def pack(x): return {"v": x}
    >>> def unpack(d): return d["v"]
    >>> register_packer("toy", (pack, unpack))
    >>> get_packer("toy") is not None
    True
    
    See Also
    --------
    get_packer
    list_packers
    pack_to_file
    unpack_from_file
    
    References
    ----------
    .. [1] NumPy. ``np.savez_compressed`` for NPZ payloads.
    """

    if not kind or not isinstance(kind, str):
        raise ValueError("kind must be non-empty str")
    _PACKERS[kind.lower()] = packer


def get_packer(kind: str) -> Optional[Packer]:
    r"""
    Return the registered packer for ``kind`` or ``None``.
    
    Parameters
    ----------
    kind : str
        Case-insensitive key used at registration time.
    
    Returns
    -------
    tuple(callable, callable) or None
        The ``(pack, unpack)`` pair if found, else ``None``.
    
    Examples
    --------
    >>> _ = get_packer("missing") is None
    >>> # After registering:
    >>> # register_packer("toy", (pack, unpack))
    >>> pk = get_packer("toy")
    >>> callable(pk[0]) and callable(pk[1])
    True
    
    See Also
    --------
    register_packer
    list_packers
    """

    return _PACKERS.get(kind.lower())


def list_packers() -> Dict[str, str]:
    r"""
    List available packers as a mapping of kind to signatures.
    
    Returns
    -------
    dict[str, str]
        Mapping ``kind -> "pack_name | unpack_name"`` for quick
        inspection.
    
    Notes
    -----
    Only names available via ``__name__`` are shown. For
    callables without a name, ``repr`` is used.
    
    Examples
    --------
    >>> isinstance(list_packers(), dict)
    True
    
    See Also
    --------
    register_packer
    get_packer
    """

    out: Dict[str, str] = {}
    for k, (pk, up) in _PACKERS.items():
        a = getattr(pk, "__name__", repr(pk))
        b = getattr(up, "__name__", repr(up))
        out[k] = f"{a} | {b}"
    return out


def _bundle_pack(obj: Any) -> Dict[str, Any]:
    r"""
    Pack a :class:`TFBundle` into an NPZ-friendly mapping.
    Validates type and converts arrays/scalars.
    """

    if not isinstance(obj, TFBundle):
        raise TypeError("expected TFBundle")

    def _arr(x):
        return None if x is None else np.asarray(x)

    payload = {
        "kind": "bundle",
        "freq": _arr(obj.freq),
        "z": _arr(obj.z),
        "z_err": _arr(obj.z_err),
        "tipper": _arr(obj.tipper),
        "tipper_err": _arr(obj.tipper_err),
        "rho": _arr(obj.rho),
        "phase": _arr(obj.phase),
        "station": obj.station,
        "station_id": obj.station_id,
        "lat": obj.lat,
        "lon": obj.lon,
        "elev": obj.elev,
        "azimuth": obj.azimuth,
        "attrs": json.dumps(obj.attrs or {}),
    }
    # drop None to avoid object arrays in NPZ
    return {k: v for k, v in payload.items() if v is not None}


def _bundle_unpack(d: Dict[str, Any]) -> TFBundle:
    r"""
    Rebuild a :class:`TFBundle` from a mapping produced by
    ``_bundle_pack``.
    
    Unboxes 0-D arrays and parses JSON attrs.
    """
    def _sc(x):
        # unbox 0-D numpy scalars to Python types
        if isinstance(x, np.ndarray) and x.shape == ():
            return x.item()
        return x

    attrs_txt = _sc(d.get("attrs"))
    attrs = json.loads(attrs_txt) if attrs_txt else {}

    return TFBundle(
        freq=d.get("freq"),
        z=d.get("z"),
        z_err=d.get("z_err"),
        tipper=d.get("tipper"),
        tipper_err=d.get("tipper_err"),
        rho=d.get("rho"),
        phase=d.get("phase"),
        station=_sc(d.get("station")),
        station_id=_sc(d.get("station_id")),
        lat=_sc(d.get("lat")),
        lon=_sc(d.get("lon")),
        elev=_sc(d.get("elev")),
        azimuth=_sc(d.get("azimuth")),
        attrs=attrs,
    )


register_packer("bundle", (_bundle_pack, _bundle_unpack))


def pack_to_file(
    obj: Any,
    path: Path | str,
    *,
    kind: str = "bundle",
) -> Path:
    r"""
    Serialize an object to ``.npz`` using a registered packer.
    
    Parameters
    ----------
    obj : Any
        Object to serialize. Must be supported by the chosen
        ``kind`` (i.e., its packer).
    path : path-like
        Output file path. ``.npz`` is recommended.
    kind : str, optional
        Packer name. Defaults to ``"bundle"``.
    
    Returns
    -------
    pathlib.Path
        The path actually written.
    
    Raises
    ------
    ValueError
        If no packer is registered for ``kind``.
    OSError
        For filesystem write errors.
    
    Notes
    -----
    Data are written with :func:`numpy.savez_compressed`, which
    stores arrays losslessly and compresses the archive.
    
    Examples
    --------
    Pack a bundle to disk::
    
        >>> # Assume a TFBundle instance: bndl
        >>> out = pack_to_file(bndl, "site001.npz", kind="bundle")
        >>> out.name.endswith(".npz")
        True
    
    See Also
    --------
    unpack_from_file
    register_packer
    get_packer
    numpy.savez_compressed
    
    References
    ----------
    .. [1] NumPy Reference. *np.savez_compressed*.
    """

    p = Path(path)
    pk = get_packer(kind)
    if pk is None:
        raise ValueError(f"no packer for kind: {kind}")
    payload = pk[0](obj)
    np.savez_compressed(p, **payload)
    return p


def unpack_from_file(path: Path | str, *, kind: str = "bundle") -> Any:
    r"""
    Deserialize an object from ``.npz`` using a registered
    packer.
    
    Parameters
    ----------
    path : path-like
        Input NPZ path created by :func:`pack_to_file`.
    kind : str, optional
        Packer name. Defaults to ``"bundle"``. Must match the
        packer used during serialization.
    
    Returns
    -------
    Any
        The reconstructed object.
    
    Raises
    ------
    ValueError
        If no packer is registered for ``kind``.
    OSError
        For filesystem read errors.
    KeyError
        If expected payload keys are missing for the packer.
    
    Notes
    -----
    Loading uses ``allow_pickle=False`` for safety; payloads must
    be composed of basic NumPy/JSON types.
    
    Examples
    --------
    Load a previously saved bundle::
    
        >>> obj = unpack_from_file("site001.npz", kind="bundle")
        >>> hasattr(obj, "z") and hasattr(obj, "freq")
        True
    
    See Also
    --------
    pack_to_file
    register_packer
    get_packer
    
    """

    p = Path(path)
    pk = get_packer(kind)
    if pk is None:
        raise ValueError(f"no packer for kind: {kind}")
    with np.load(p, allow_pickle=False) as npz:
        d = {k: npz[k] for k in npz.files}
    return pk[1](d)


class RegistryAPI(CoreObject):
    r"""
    Convenience façade over :class:`~pycsamt.core.registry.Registry`.
    
    This high-level helper wires a registry root, ensures a
    ``packs/`` folder exists under that root, and exposes both
    thin proxies to low-level registry calls and ergonomic
    helpers for packing objects to NPZ bundles and restoring
    them later.
    
    Parameters
    ----------
    root : path-like
        Base directory that will contain the manifest and, if
        needed, a ``packs/`` subdirectory for NPZ bundles.
    manifest_name : str, optional
        Manifest filename within ``root``. Defaults to
        ``"manifest.json"``.
    
    Attributes
    ----------
    low : Registry
        The underlying registry instance that performs the core
        operations (load/save, add, query).
    root : pathlib.Path
        Absolute path to the registry root directory.
    pack_dir : pathlib.Path
        Directory where NPZ bundles are stored (``root/packs``).
    
    Notes
    -----
    ``RegistryAPI`` aims to be pragmatic and file-system first.
    It does not do locking or concurrency control. For multi-
    process scenarios, add your own synchronization.
    
    The packing flow defaults to the ``"bundle"`` packer, which
    serializes :class:`TFBundle`-compatible objects into NPZ
    using :func:`numpy.savez_compressed`.
    
    Examples
    --------
    Create an API, add an EDI file, and list records::
    
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> _ = api.add_file("site001.edi", kind="edi", fmt="edi")
        >>> [r.kind for r in api.list(kind="edi")]
        ['edi']
    
    Pack an object to an NPZ bundle and store it in ``packs``::
    
        >>> class Dummy:
        ...     station = "S01"; station_id = 1
        >>> rec = api.add_object(Dummy(), pack=True,
        ...                      pack_name="S01.npz")
        >>> rec.kind
        'bundle'
    
    Materialize a record back into a Python object::
    
        >>> obj = api.materialize(rec.rid)  # may be Dummy or path
    
    Convert a record to an EDI representation::
    
        >>> edi_like = api.to_edi(rec.rid)
    
    See Also
    --------
    pycsamt.core.registry.Registry
    pycsamt.core.registry.pack_to_file
    pycsamt.core.registry.unpack_from_file
    pycsamt.core.registry.guess_kind
    pycsamt.seg.edi.EDIFile
    pycsamt.zonge.avg.AVG
    pycsamt.jones.j.JFile
    pycsamt.core.base.to_edi
    
    References
    ----------
    .. [1] D. Wight (1991). *SEG MT/EMAP EDI Standard*.
    .. [2] NumPy Reference. *np.savez_compressed*, *np.load*.
    """

    def __init__(
        self,
        root: Path | str,
        *,
        manifest_name: str = "manifest.json",
    ) -> None:
        self.low = Registry(
            root,
            manifest_name=manifest_name,
        )
        self.root = Path(self.low.root)
        self.pack_dir = self.root / "packs"
        self.pack_dir.mkdir(parents=True, exist_ok=True)

    # direct proxies 
    def save(self) -> None:
        r"""
        Persist the manifest to disk.
        
        Thin proxy to :meth:`Registry.save`. Use after mutating the
        underlying manifest directly, though typical helpers call
        save for you.
        """
        self.low.save()

    def list(self, *, kind: Optional[str] = None):
        r"""
        List all records, optionally filtered by kind.
        
        Parameters
        ----------
        kind : str or None, optional
            If provided, only records whose ``kind`` matches are
            returned.
        
        Returns
        -------
        list of Record
            Shallow views of stored :class:`Record` objects.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> isinstance(api.list(), list)
        True
        """
        return self.low.list(kind=kind)

    def get(self, rid: str):
        r"""
        Return the record by id or raise an error.
        
        Parameters
        ----------
        rid : str
            Record identifier returned at insertion time.
        
        Returns
        -------
        Record
            The matching record.
        
        Raises
        ------
        RegistryError
            If the id does not exist in the manifest.
        """

        return self.low.get(rid)

    def find(
        self,
        *,
        tag: Optional[str] = None,
        kind: Optional[str] = None,
        dataid: Optional[str] = None,
    ):
        r"""
        Filter records by tag, kind, and/or data id.
        
        Parameters
        ----------
        tag : str or None, optional
            Require that the tag is present in ``Record.tags``.
        kind : str or None, optional
            Require that the record kind matches the value.
        dataid : str or None, optional
            Require equality with ``Record.dataid``.
        
        Returns
        -------
        list of Record
            Records satisfying all given predicates.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> api.find(kind="edi")  # doctest: +ELLIPSIS
        [...]  # list of records
        """

        return self.low.find(tag=tag, kind=kind, dataid=dataid)

    # add helpers 
    def add_file(
        self,
        path: Path | str,
        *,
        kind: Optional[str] = None,
        fmt: Optional[str] = None,
        dataid: Optional[str] = None,
        station_id: Optional[str] = None,
        tags: Optional[list[str]] = None,
        meta: Optional[dict[str, Any]] = None,
        with_hash: bool = True,
    ):
        r"""
        Register a file path as a new record.
        
        Parameters
        ----------
        path : path-like
            File path. Relative paths are resolved under ``root``.
        kind : str or None, optional
            Logical kind (e.g. ``"edi"``, ``"avg"``). If omitted,
            defaults to ``"meta"``.
        fmt : str or None, optional
            Short format hint (e.g. ``"edi"``).
        dataid : str or None, optional
            External data or survey id.
        station_id : str or int or None, optional
            Station identifier. Cast to string if provided.
        tags : list[str] or None, optional
            User-defined labels to ease later queries.
        meta : dict[str, Any] or None, optional
            Extra JSON-serializable metadata.
        with_hash : bool, optional
            If ``True`` and the file exists, compute and store a
            SHA-256 checksum. Defaults to ``True``.
        
        Returns
        -------
        Record
            The created and persisted record.
        
        Notes
        -----
        This method saves the manifest on success.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> rec = api.add_file("site001.edi",
        ...                    kind="edi", fmt="edi")
        >>> rec.kind, rec.fmt
        ('edi', 'edi')
        """

        return self.low.add_path(
            path,
            kind=kind,
            fmt=fmt,
            dataid=dataid,
            station_id=station_id,
            tags=tags,
            meta=meta,
            with_hash=with_hash,
        )

    def add_object(
        self,
        obj: Any,
        *,
        kind: Optional[str] = None,
        tags: Optional[list[str]] = None,
        meta: Optional[dict[str, Any]] = None,
        pack: bool = False,
        pack_name: Optional[str] = None,
    ):
        r"""
        Register an object and, optionally, pack it to an NPZ bundle.
        
        Parameters
        ----------
        obj : Any
            Object to register. If ``kind`` is ``None``, a coarse
            label is inferred via :func:`guess_kind`.
        kind : str or None, optional
            Explicit kind override (e.g. ``"edi"``, ``"bundle"``).
        tags : list[str] or None, optional
            User-defined labels.
        meta : dict[str, Any] or None, optional
            Extra JSON-serializable metadata for the record.
        pack : bool, optional
            If ``True``, serialize ``obj`` with the ``"bundle"``
            packer into ``pack_dir`` and rewrite the record to refer
            to that NPZ on disk. Defaults to ``False``.
        pack_name : str or None, optional
            Output NPZ filename when ``pack=True``. Defaults to
            ``"<rid>.npz"``.
        
        Returns
        -------
        Record
            The created record. When ``pack=True``, the record
            ``kind`` becomes ``"bundle"`` and ``path`` points to the
            NPZ file inside ``pack_dir``.
        
        Notes
        -----
        Packing uses :func:`pack_to_file` with ``kind="bundle"``.
        The manifest is saved after packing.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> class Obj: station="S02"; station_id=2
        >>> rec = api.add_object(Obj(), tags=["demo"])
        >>> isinstance(rec.rid, str)
        True
        
        Pack immediately to NPZ::
        
            >>> rec2 = api.add_object(Obj(), pack=True,
            ...                       pack_name="S02.npz")
            >>> rec2.kind
            'bundle'
        
        See Also
        --------
        pycsamt.core.registry.pack_to_file
        pycsamt.core.registry.guess_kind
        """

        r = self.low.add_obj(
            obj,
            kind=kind or guess_kind(obj),
            tags=tags,
            meta=meta,
        )
        if pack:
            name = pack_name or f"{r.rid}.npz"
            out = self.pack_dir / name
            pack_to_file(obj, out, kind="bundle")
            r.path = str(out)
            r.kind = "bundle"
            self.save()
        return r

    # materializers 
    def materialize(self, rid: str) -> Any:
        r"""
        Load a record into a concrete Python object when possible.
        
        Parameters
        ----------
        rid : str
            Record id to resolve.
        
        Returns
        -------
        Any
            If the record is a ``"bundle"`` with a valid path,
            returns the unpacked object via
            :func:`unpack_from_file`. For known kinds with readable
            files this returns:
        
            * ``"edi"``  -> :class:`pycsamt.seg.edi.EDIFile`
            * ``"avg"``  -> :class:`pycsamt.zonge.avg.AVG`
            * ``"j"`` or ``"j_col"`` -> :class:`pycsamt.jones.j.JFile`
        
            Otherwise, returns the string path (or ``None``) of the
            record.
        
        Notes
        -----
        Errors during import or read are swallowed and the path is
        returned instead. This keeps the method side-effect free
        with respect to the registry while remaining robust.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> # Assuming a valid record id:
        >>> # obj = api.materialize(rid)
        """

        r = self.low.get(rid)
        if r.kind == "bundle" and r.path:
            return unpack_from_file(r.path)
        p = Path(r.path) if r.path else None
        try:
            if r.kind == "edi" and p and p.exists():
                from pycsamt.seg.edi import EDIFile
                return EDIFile.from_file(p)
            if r.kind == "avg" and p and p.exists():
                from pycsamt.zonge.avg import AVG
                return AVG.from_file(p)
            if r.kind in {"j", "j_col"} and p and p.exists():
                from pycsamt.jones.j import JFile
                return JFile.from_file(p)
        except Exception:
            return r.path
        return r.path

    def to_edi(self, rid: str, **kw: Any) -> Any:
        r"""
        Convert a record's materialized object into an EDI-like form.
        
        Parameters
        ----------
        rid : str
            Record id to resolve and convert.
        **kw : Any
            Forwarded to :func:`to_edi`, which performs the actual
            conversion.
        
        Returns
        -------
        Any
            An object produced by :func:`to_edi`. The concrete type
            depends on the available converters and inputs.
        
        Notes
        -----
        This is a convenience wrapper that first calls
        :meth:`materialize`, then delegates to
        :func:`pycsamt.core.base.to_edi`.
        
        Examples
        --------
        >>> from pycsamt.core.registry import RegistryAPI
        >>> api = RegistryAPI("data")
        >>> # Convert an existing record to an EDI representation
        >>> edi_like = api.to_edi("some-record-id")  # doctest: +SKIP
        
        See Also
        --------
        pycsamt.core.base.to_edi
        pycsamt.seg.edi.EDIFile
        """

        obj = self.materialize(rid)
        return to_edi(obj, **kw)
