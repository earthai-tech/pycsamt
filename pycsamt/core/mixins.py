# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Iterator, Optional

from .base import TFBundle, ensure_station, to_edi

__all__ = [
    "bundle_from_edi",
    "BundleMixin",
    "BundleContainerMixin",
]


def _get(obj: Any, *names: str) -> Any:
    for n in names:
        if hasattr(obj, n):
            return getattr(obj, n)
    return None


def bundle_from_edi(obj: Any) -> TFBundle:
    freq = _get(obj, "freq", "frequency")
    z = _get(obj, "z", "Z")
    z_err = _get(obj, "z_err", "Z_err", "z_error")
    tip = _get(obj, "tipper", "Tipper", "tip")
    tip_err = _get(obj, "tipper_err", "tip_err")
    rho = _get(obj, "resistivity", "rho")
    pha = _get(obj, "phase", "phi", "phase_deg")
    name = _get(obj, "station", "site", "name")
    lat = _get(obj, "lat", "latitude")
    lon = _get(obj, "lon", "longitude")
    elev = _get(obj, "elev", "elevation")
    az = _get(obj, "azimuth", "az")

    # normalize tipper shape to (n, 2)
    try:
        if tip is not None and getattr(tip, "ndim", 1) == 3:
            if tip.shape[1:] == (1, 2):
                tip = tip[:, 0, :]
    except Exception:
        pass

    return TFBundle(
        freq=freq,
        z=z,
        z_err=z_err,
        tipper=tip,
        tipper_err=tip_err,
        rho=rho,
        phase=pha,
        station=name,
        lat=lat,
        lon=lon,
        elev=elev,
        azimuth=az,
    )


def _looks_collection(obj: Any) -> bool:
    if isinstance(obj, (list, tuple, set)):
        return True
    if isinstance(obj, (str, bytes)):
        return False
    if hasattr(obj, "__iter__") and not hasattr(obj, "z"):
        return True
    return False


class BundleMixin:
    def to_bundle(self) -> TFBundle:  # noqa: D401
        raise NotImplementedError

    @classmethod
    def from_bundle(cls, bundle: TFBundle):  # noqa: D401
        raise NotImplementedError

    @staticmethod
    def ensure_station_name(
        name: Optional[str],
        station_id: Optional[str | int],
    ) -> str:
        return ensure_station(name, station_id)

    def to_edi(self, *, key: Optional[str] = None, **k: Any) -> Any:
        return to_edi(self, key=key, **k)

    @classmethod
    def from_edi(cls, edi_obj: Any):
        if _looks_collection(edi_obj):
            out = []
            for it in edi_obj:  # type: ignore
                b = bundle_from_edi(it)
                out.append(cls.from_bundle(b))
            return out
        b = bundle_from_edi(edi_obj)
        return cls.from_bundle(b)


class BundleContainerMixin(BundleMixin):
    def iter_bundles(self) -> Iterator[TFBundle]:
        if hasattr(self, "items"):
            for _, v in self.items():  # type: ignore
                if hasattr(v, "to_bundle"):
                    yield v.to_bundle()
        elif hasattr(self, "__iter__"):
            for v in self:  # type: ignore
                if hasattr(v, "to_bundle"):
                    yield v.to_bundle()

    def to_edi_collection(
        self,
        *,
        key: Optional[str] = None,
        **k: Any,
    ) -> list[Any]:
        edis: list[Any] = []
        for b in self.iter_bundles():
            proxy = _ProxyFromBundle(b)
            edis.append(to_edi(proxy, key=key, **k))
        return edis


class _ProxyFromBundle:
    def __init__(self, b: TFBundle) -> None:
        self._b = b

    def to_bundle(self) -> TFBundle:
        return self._b
