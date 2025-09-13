# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Callable, Dict, Optional, Tuple
from pathlib import Path
import json

import numpy as np

from . import _registry as low
from .base import TFBundle
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
    if not kind or not isinstance(kind, str):
        raise ValueError("kind must be non-empty str")
    _PACKERS[kind.lower()] = packer


def get_packer(kind: str) -> Optional[Packer]:
    return _PACKERS.get(kind.lower())


def list_packers() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for k, (pk, up) in _PACKERS.items():
        a = getattr(pk, "__name__", repr(pk))
        b = getattr(up, "__name__", repr(up))
        out[k] = f"{a} | {b}"
    return out


def _bundle_pack(obj: Any) -> Dict[str, Any]:
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
    p = Path(path)
    pk = get_packer(kind)
    if pk is None:
        raise ValueError(f"no packer for kind: {kind}")
    payload = pk[0](obj)
    np.savez_compressed(p, **payload)
    return p


def unpack_from_file(path: Path | str, *, kind: str = "bundle") -> Any:
    p = Path(path)
    pk = get_packer(kind)
    if pk is None:
        raise ValueError(f"no packer for kind: {kind}")
    with np.load(p, allow_pickle=False) as npz:
        d = {k: npz[k] for k in npz.files}
    return pk[1](d)


class RegistryAPI:
    def __init__(
        self,
        root: Path | str,
        *,
        manifest_name: str = "manifest.json",
    ) -> None:
        self.low = low.Registry(
            root,
            manifest_name=manifest_name,
        )
        self.root = Path(self.low.root)
        self.pack_dir = self.root / "packs"
        self.pack_dir.mkdir(parents=True, exist_ok=True)

    # direct proxies -------------------------------------------------------
    def save(self) -> None:
        self.low.save()

    def list(self, *, kind: Optional[str] = None):
        return self.low.list(kind=kind)

    def get(self, rid: str):
        return self.low.get(rid)

    def find(
        self,
        *,
        tag: Optional[str] = None,
        kind: Optional[str] = None,
        dataid: Optional[str] = None,
    ):
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
        r = self.low.add_obj(
            obj,
            kind=kind or low.guess_kind(obj),
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
        obj = self.materialize(rid)
        return to_edi(obj, **kw)
