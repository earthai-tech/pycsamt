# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Protocol


from .config import StationNamePolicy, get_config, get_adapter

__all__ = [
    "TFBundle",
    "SupportsToBundle",
    "SupportsFromBundle",
    "ensure_station",
    "pick_adapter_key",
    "to_edi",
]


@dataclass
class TFBundle:
    # neutral payload used across backends
    freq: Any | None = None
    z: Any | None = None
    z_err: Any | None = None
    tipper: Any | None = None
    tipper_err: Any | None = None
    rho: Any | None = None
    phase: Any | None = None
    station: str | None = None
    station_id: str | int | None = None
    lat: float | None = None
    lon: float | None = None
    elev: float | None = None
    azimuth: float | None = None
    attrs: Dict[str, Any] = field(default_factory=dict)

    def is_empty(self) -> bool:
        have_z = self.z is not None
        have_rp = (self.rho is not None) and (self.phase is not None)
        return not (have_z or have_rp)


class SupportsToBundle(Protocol):
    def to_bundle(self) -> TFBundle: ...  # noqa: E701


class SupportsFromBundle(Protocol):
    @classmethod
    def from_bundle(cls, bundle: TFBundle): ...  # noqa: E701


def ensure_station(
    name: Optional[str],
    station_id: Optional[str | int],
    *,
    policy: Optional[StationNamePolicy] = None,
) -> str:
    pol = policy or get_config().station_policy
    return pol.ensure(name, station_id)


def pick_adapter_key(
    obj: Any,
    *,
    hint: Optional[str] = None,
) -> Optional[str]:
    if hint:
        return hint.lower()

    try:
        mod = obj.__class__.__module__.lower()
        cls = obj.__class__.__name__.lower()
    except Exception:
        return None

    if ("zonge" in mod) or ("avg" in cls):
        return "avg"
    if ("jones" in mod) or (cls in {"jfile", "jcollection"}):
        return "j"
    if ("seg" in mod) and ("edi" in (mod + cls)):
        return "edi"
    return None


def to_edi(
    source: Any,
    *,
    key: Optional[str] = None,
    **kwargs: Any,
) -> Any:
    k = key or pick_adapter_key(source)
    if not k:
        raise RuntimeError("Cannot infer adapter key")

    factory = get_adapter(k)
    if not factory:
        raise RuntimeError(f"No adapter registered for: {k}")

    return factory(source, **kwargs)
