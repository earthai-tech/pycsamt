# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from typing import Any, Dict, Iterable, List, Optional
from pathlib import Path
import hashlib
import json
import time
import uuid

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
    pass

def _now() -> float:
    return float(time.time())


def _uuid() -> str:
    return uuid.uuid4().hex


def _sha256(path: Path, block: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def guess_kind(obj: Any) -> str:
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
class Record:
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
        self.updated = _now()

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Record":
        return cls(**d)


@dataclass
class Manifest:
    root: str
    version: int = 1
    records: Dict[str, Record] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "root": self.root,
            "version": self.version,
            "records": {k: r.to_dict() for k, r in self.records.items()},
        }

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Manifest":
        recs = {
            k: Record.from_dict(v) 
            for k, v in d.get("records", {}).items()
        }
        return cls(
            root=d.get("root", "."), version=int(
                d.get("version", 1)), records=recs)


class ManifestStore:
    def load(self, path: Path) -> Manifest:  # noqa: D401
        raise NotImplementedError

    def save(self, man: Manifest, path: Path) -> None:  # noqa: D401
        raise NotImplementedError


class FileManifestStore(ManifestStore):
    def load(self, path: Path) -> Manifest:
        if not path.exists():
            return Manifest(root=str(path.parent))
        txt = path.read_text(encoding="utf-8")
        if path.suffix.lower() == ".toml" and tomllib is not None:
            d = tomllib.loads(txt)
        else:
            d = json.loads(txt)
        return Manifest.from_dict(d)

    def save(self, man: Manifest, path: Path) -> None:
        d = man.to_dict()
        if path.suffix.lower() == ".toml":
            raise RegistryError("TOML save not supported here")
        path.write_text(json.dumps(d, indent=2), encoding="utf-8")


class Registry:
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

    # ---- basic ops -------------------------------------------------------
    def save(self) -> None:
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
        try:
            return self.manifest.records[rid]
        except KeyError as exc:
            raise RegistryError(f"unknown rid: {rid}") from exc

    def list(self, *, kind: Optional[str] = None) -> List[Record]:
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
        r = self.get(rid)
        r.meta.update(fields)
        r.touch()
        self.save()
        return r

    def remove(self, rid: str) -> None:
        self.manifest.records.pop(rid, None)
        self.save()
