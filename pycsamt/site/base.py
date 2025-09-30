# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Callable, Optional
from typing import List, Tuple, Sequence, Dict, Union

import math
import copy
from pathlib import Path

import numpy as np
import pandas as pd

from ..constants import _EARTH_R
from ..core.base import CoreObject
from ..seg.edi import EDIFile
from ..seg.collection import EDICollection
from ..session import normalize_session

from ..gis.utils import (
    assert_lat_value,
    assert_lon_value,
    assert_elevation_value,
    normalize_lat_lon,
)


__all__ = ["SiteMixin", "Site", "Sites"]


def _safe_get(d: Any, *names: str, default: Any = None) -> Any:
    for n in names:
        try:
            v = getattr(d, n)
            return v
        except Exception:
            pass
        try:
            v = d[n]  # type: ignore[index]
            return v
        except Exception:
            pass
    return default


def _clone_edi(edi: EDIFile) -> EDIFile:
    try:
        return copy.deepcopy(edi)
    except Exception:
        # best effort shallow clone
        c = EDIFile()
        for k, v in getattr(edi, "__dict__", {}).items():
            try:
                setattr(c, k, copy.copy(v))
            except Exception:
                setattr(c, k, v)
        return c


def _station_name(edi: EDIFile) -> str:
    head = _safe_get(edi, "HEAD", default=None)
    for k in ("station", "sitename", "name", "STATION"):
        v = _safe_get(head, k, default=None)
        if v:
            return str(v)
    return str(_safe_get(edi, "name", default="")) or "site"


def _set_station_name(edi: EDIFile, name: str) -> None:
    head = _safe_get(edi, "HEAD", default=None)
    if head is None:
        try:
            edi.HEAD = type("HEAD", (), {})()  # type: ignore
            head = edi.HEAD  # type: ignore[attr-defined]
        except Exception:
            pass
    for k in ("station", "sitename", "name", "STATION"):
        try:
            setattr(head, k, name)
            return
        except Exception:
            continue
    try:
        edi.name = name
    except Exception:
        pass


def _coords_tuple(edi: EDIFile) -> Tuple[float, float, float]:
    head = _safe_get(edi, "HEAD", default=None)
    lat = _safe_get(head, "lat", "LAT", default=None)
    lon = _safe_get(head, "lon", "LON", default=None)
    elev = _safe_get(head, "elev", "elevation", "ELEV", default=0)
    try:
        if lat is None or lon is None:
            return (np.nan, np.nan, float(elev or 0.0))
        latv = assert_lat_value(lat)  # type: ignore[arg-type]
        lonv = assert_lon_value(lon)  # type: ignore[arg-type]
        elevv = assert_elevation_value(elev)
        return (float(latv), float(lonv), float(elevv))
    except Exception:
        return (np.nan, np.nan, 0.0)


def _set_coords_on_head(
    edi: EDIFile,
    lat: float,
    lon: float,
    elev: Optional[float] = None,
) -> None:
    head = _safe_get(edi, "HEAD", default=None)
    if head is None:
        try:
            edi.HEAD = type("HEAD", (), {})()  # type: ignore
            head = edi.HEAD  # type: ignore[attr-defined]
        except Exception:
            return
    try:
        setattr(head, "lat", float(lat))
        setattr(head, "lon", float(lon))
        if elev is not None:
            setattr(head, "elev", float(elev))
    except Exception:
        pass


def _component_names() -> List[str]:
    return ["Zxx", "Zxy", "Zyx", "Zyy"]


def _extract_z_arrays(edi: EDIFile) -> Dict[str, Any]:
    Z = _safe_get(edi, "Z", default=None)
    out: Dict[str, Any] = {}
    out["freq"] = _safe_get(Z, "freq", "frequency", default=None)
    out["z"] = _safe_get(Z, "z", "impedance", default=None)
    out["z_err"] = _safe_get(
        Z, "z_error", "z_err", "impedance_err", default=None
    )
    out["rho"] = _safe_get(Z, "rho", "res", "resistivity", default=None)
    out["phase"] = _safe_get(Z, "phase", "phi", default=None)
    # tipper may live elsewhere; try common spots
    tip = _safe_get(edi, "T", "TIP", default=None)
    if tip is None:
        tip = _safe_get(Z, "tipper", "tip", default=None)
    out["tipper"] = tip
    return out


class SiteMixin(CoreObject):
    def __init__(self, edi: EDIFile) -> None:
        self.edi = edi

    # ---------- basic properties ----------

    @property
    def name(self) -> str:
        return _station_name(self.edi)

    @property
    def freq(self) -> Any:
        return _extract_z_arrays(self.edi)["freq"]

    @property
    def z(self) -> Any:
        return _extract_z_arrays(self.edi)["z"]

    @property
    def z_err(self) -> Any:
        return _extract_z_arrays(self.edi)["z_err"]

    @property
    def rho(self) -> Any:
        return _extract_z_arrays(self.edi)["rho"]

    @property
    def phase(self) -> Any:
        return _extract_z_arrays(self.edi)["phase"]

    @property
    def tipper(self) -> Any:
        return _extract_z_arrays(self.edi)["tipper"]

    @property
    def coords(self) -> Tuple[float, float, float]:
        return _coords_tuple(self.edi)

    @property
    def meta(self) -> Dict[str, Any]:
        head = _safe_get(self.edi, "HEAD", default=None)
        info = _safe_get(self.edi, "INFO", default=None)
        out: Dict[str, Any] = {}
        for k in ("station", "name", "lat", "lon", "elev"):
            v = _safe_get(head, k, default=None)
            if v is not None:
                out[k] = v
        if info is not None:
            out["INFO"] = dict(getattr(info, "__dict__", {}))
        return out

    # ---------- reads ----------

    def to_dataframe(self, kind: str = "z") -> pd.DataFrame:
        arrs = _extract_z_arrays(self.edi)
        f = np.asarray(arrs["freq"])
        if kind.lower() in ("z", "imp", "impedance"):
            z = np.asarray(arrs["z"])
            if z is None or z.size == 0:
                return pd.DataFrame(index=pd.Index([], name="f"))
            cols = _component_names()
            df = {}
            for i, c in enumerate(cols):
                try:
                    df[c] = z[:, i]
                except Exception:
                    df[c] = np.full_like(f, np.nan, dtype=float)
            return pd.DataFrame(df, index=pd.Index(f, name="f"))
        if kind.lower() in ("resphase", "rp", "rho_phase"):
            rho = np.asarray(arrs["rho"])
            ph = np.asarray(arrs["phase"])
            cols = _component_names()
            df: Dict[str, Any] = {}
            for i, c in enumerate(cols):
                rr = np.full_like(f, np.nan, dtype=float)
                pp = np.full_like(f, np.nan, dtype=float)
                try:
                    rr = rho[:, i]
                except Exception:
                    pass
                try:
                    pp = ph[:, i]
                except Exception:
                    pass
                df[f"rho_{c.lower()}"] = rr
                df[f"phi_{c.lower()}"] = pp
            return pd.DataFrame(df, index=pd.Index(f, name="f"))
        if kind.lower() in ("tip", "tipper", "t"):
            tip = arrs["tipper"]
            if tip is None:
                return pd.DataFrame(index=pd.Index(f, name="f"))
            tip = np.asarray(tip)
            df = {}
            for i, c in enumerate(("Tx", "Ty")):
                try:
                    df[c] = tip[:, i]
                except Exception:
                    df[c] = np.full_like(f, np.nan, dtype=float)
            return pd.DataFrame(df, index=pd.Index(f, name="f"))
        raise ValueError("Unknown kind: {!r}".format(kind))

    def quality_flags(self) -> Dict[str, bool]:
        arrs = _extract_z_arrays(self.edi)
        f_ok = arrs["freq"] is not None
        z = np.asarray(arrs["z"]) if arrs["z"] is not None else None
        e = (
            np.asarray(arrs["z_err"])
            if arrs["z_err"] is not None
            else None
        )
        ok = lambda a: a is not None and np.all(np.isfinite(a))
        return {
            "has_freq": bool(f_ok),
            "has_z": bool(ok(z)),
            "has_z_err": bool(ok(e)),
            "has_rho": bool(ok(arrs["rho"])),
            "has_phase": bool(ok(arrs["phase"])),
            "has_tipper": bool(ok(arrs["tipper"])),
        }

    def has_component(self, comp: str) -> bool:
        comp = comp.strip().lower()
        if comp in ("tip", "tx", "ty", "tipper"):
            tip = _extract_z_arrays(self.edi)["tipper"]
            if tip is None:
                return False
            a = np.asarray(tip)
            return a.size > 0 and np.any(np.isfinite(a))
        z = _extract_z_arrays(self.edi)["z"]
        if z is None:
            return False
        z = np.asarray(z)
        names = _component_names()
        try:
            i = [n.lower() for n in names].index(comp)
        except ValueError:
            return False
        a = z[:, i]
        return a.size > 0 and np.any(np.isfinite(a))

    def summary(self) -> Dict[str, Any]:
        lat, lon, elev = self.coords
        arrs = _extract_z_arrays(self.edi)
        nfreq = int(np.size(arrs["freq"] or []))
        present = [c for c in _component_names() if self.has_component(c)]
        return {
            "name": self.name,
            "nfreq": nfreq,
            "lat": lat,
            "lon": lon,
            "elev": elev,
            "components": present,
            "tipper": bool(self.has_component("tipper")),
        }

    # ---------- light writes (copy-on-write) ----------

    def rename(
        self,
        new: str,
        *,
        inplace: bool = False,
    ) -> "Site":
        if inplace:
            _set_station_name(self.edi, new)
            return self  # type: ignore[return-value]
        edi = _clone_edi(self.edi)
        _set_station_name(edi, new)
        return Site(edi)

    def set_coords(
        self,
        lat: Union[float, str],
        lon: Union[float, str],
        elev: Optional[Union[float, str]] = None,
        *,
        inplace: bool = False,
    ) -> "Site":
        la, lo = normalize_lat_lon(lon, lat, assume="lonlat")
        la = float(assert_lat_value(la))  # type: ignore[arg-type]
        lo = float(assert_lon_value(lo))  # type: ignore[arg-type]
        ev = (
            None
            if elev is None
            else float(assert_elevation_value(elev))
        )
        if inplace:
            _set_coords_on_head(self.edi, la, lo, ev)
            return self  # type: ignore[return-value]
        edi = _clone_edi(self.edi)
        _set_coords_on_head(edi, la, lo, ev)
        return Site(edi)

    def set_empty(self, *, inplace: bool = False) -> "Site":
        target = self.edi if inplace else _clone_edi(self.edi)
        Z = _safe_get(target, "Z", default=None)
        try:
            if Z is None:
                target.Z = type("Z", (), {})()  # type: ignore
                Z = target.Z  # type: ignore[attr-defined]
            setattr(Z, "freq", np.asarray([], dtype=float))
            setattr(Z, "z", np.asarray([], dtype=float))
            setattr(Z, "z_error", np.asarray([], dtype=float))
            setattr(Z, "rho", np.asarray([], dtype=float))
            setattr(Z, "phase", np.asarray([], dtype=float))
        except Exception:
            pass
        return self if inplace else Site(target)


class Site(SiteMixin):
    def __init__(self, edi: EDIFile) -> None:
        super().__init__(edi)

    def __repr__(self) -> str:
        s = self.summary()
        return (
            f"Site(name={s['name']!r}, nfreq={s['nfreq']}, "
            f"coords=({s['lat']:.5f},{s['lon']:.5f},{s['elev']:.1f}))"
        )


class Sites(CoreObject):
    def __init__(
        self,
        edic: Union[EDICollection, Sequence[EDIFile]],
    ) -> None:
        if isinstance(edic, EDICollection):
            items = list(edic)
        else:
            items = list(edic)
        self._items: List[Site] = [Site(e) for e in items]

    # ---------- basic container ----------

    def __len__(self) -> int:
        return len(self._items)

    def __iter__(self):
        return iter(self._items)

    def __getitem__(self, key: Union[int, str]) -> Site:
        if isinstance(key, int):
            return self._items[key]
        nm = str(key).lower()
        for s in self._items:
            if s.name.lower() == nm:
                return s
        raise KeyError(key)

    def by_index(self, i: int) -> Site:
        return self._items[i]

    def get(self, name: str) -> Optional[Site]:
        try:
            return self[name]
        except Exception:
            return None

    def as_list(self) -> List[EDIFile]:
        return [s.edi for s in self._items]

    # ---------- lookup helpers ----------

    def closest(
        self,
        lat: float,
        lon: float,
        tol: Optional[float] = None,
    ) -> Optional[Site]:
        la = float(assert_lat_value(lat))  # type: ignore[arg-type]
        lo = float(assert_lon_value(lon))  # type: ignore[arg-type]

        def hav(a, b, c, d):
            d1 = math.radians(c - a)
            d2 = math.radians(d - b)
            A = (
                math.sin(d1 / 2) ** 2
                + math.cos(math.radians(a))
                * math.cos(math.radians(c))
                * math.sin(d2 / 2) ** 2
            )
            return 2 * _EARTH_R * math.asin(min(1.0, math.sqrt(A)))

        best = (None, float("inf"))  # type: ignore[var-annotated]
        for s in self._items:
            sla, slo, _ = s.coords
            if not np.isfinite(sla) or not np.isfinite(slo):
                continue
            d = hav(la, lo, float(sla), float(slo))
            if d < best[1]:
                best = (s, d)
        if best[0] is None:
            return None
        if tol is not None and best[1] > tol:
            return None
        return best[0]

    # ---------- broadcast/edit ----------

    def map(self, fn: Callable[[Site], Any]) -> List[Any]:
        return [fn(s) for s in self._items]

    def edit_all(
        self,
        *,
        rename: Optional[Callable[[str], str]] = None,
        freq_slice: Optional[slice] = None,
        mask: Optional[Callable[[pd.DataFrame], pd.Series]] = None,
        inplace: bool = False,
    ) -> "Sites":
        items: List[EDIFile] = []
        for s in self._items:
            t = s if inplace else Site(_clone_edi(s.edi))
            if rename:
                t = t.rename(rename(t.name), inplace=True)
            if freq_slice is not None:
                arrs = _extract_z_arrays(t.edi)
                try:
                    f = np.asarray(arrs["freq"]) # #noqa VARIABLE defined not used. 
                    sl = freq_slice
                    Z = _safe_get(t.edi, "Z", default=None)
                    if Z is not None:
                        for fld in ("freq", "z", "z_error",
                                    "rho", "phase"):
                            a = _safe_get(Z, fld, default=None)
                            if a is not None:
                                a = np.asarray(a)
                                setattr(Z, fld, a[sl])
                except Exception:
                    pass
            if mask is not None:
                try:
                    df = t.to_dataframe("z")
                    m = mask(df)
                    Z = _safe_get(t.edi, "Z", default=None)
                    if Z is not None and hasattr(Z, "z"):
                        z = np.asarray(getattr(Z, "z"))
                        z[~np.asarray(m)] = np.nan
                        setattr(Z, "z", z)
                except Exception:
                    pass
            items.append(t.edi)
        return self if inplace else Sites(items)

    def with_topography(self, frame: Any) -> "Sites":
        # placeholder for integration with survey topo
        _ = frame
        return self

    # ---------- selection ----------

    def select(
        self,
        names: Optional[Sequence[str]] = None,
        predicate: Optional[Callable[[Site], bool]] = None,
    ) -> "Sites":
        out: List[EDIFile] = []
        if names:
            wanted = {str(n).lower() for n in names}
            out.extend(
                s.edi for s in self._items
                if s.name.lower() in wanted
            )
        elif predicate:
            out.extend(s.edi for s in self._items if predicate(s))
        else:
            out = self.as_list()
        return Sites(out)

    # ---------- import/export/profile ----------

    @classmethod
    def from_any(
        cls,
        source: Any,
        topo_src: Any | None = None,
    ) -> "Sites":
        with normalize_session(".tmp", topo_src=topo_src) as nz:
            obj = nz.load(source)
        if isinstance(obj, EDICollection):
            return cls(obj)
        if isinstance(obj, list) and all(
            isinstance(x, EDIFile) for x in obj
        ):
            return cls(obj)  # type: ignore[arg-type]
        if isinstance(obj, EDIFile):
            return cls([obj])
        # last resort: try to wrap iterable
        try:
            return cls(list(obj))  # type: ignore[arg-type]
        except Exception:
            return cls([])

    def write(
        self,
        outdir: Union[str, Path],
        *,
        template: str = "{station}.edi",
        exist_ok: bool = False,
    ) -> List[Path]:
        outp = Path(outdir)
        outp.mkdir(parents=True, exist_ok=True)
        paths: List[Path] = []
        for s in self._items:
            name = s.name or "site"
            fn = template.format(station=name)
            p = outp / fn
            if p.exists() and not exist_ok:
                raise FileExistsError(str(p))
            edi = s.edi
            tf = getattr(edi, "to_file", None)
            if callable(tf):
                try:
                    tf(str(p))
                    paths.append(p)
                    continue
                except Exception:
                    pass
            # fallback dumb writer
            with open(p, "w", encoding="utf-8") as f:
                f.write(f"# EDI placeholder for {name}\n")
            paths.append(p)
        return paths

    def to_profile(
        self,
        origin: Tuple[float, float],
        azimuth: float,
        *,
        crs: Optional[int] = None,
    ) -> Any:
        try:
            from .profile import to_profile as _to_prof  # type: ignore  # noqa: E501
            return _to_prof(self, origin, azimuth, crs=crs)
        except Exception:
            # lightweight chainage sort along azimuth
            la0, lo0 = origin
            la0 = float(assert_lat_value(la0))  # type: ignore[arg-type]
            lo0 = float(assert_lon_value(lo0))  # type: ignore[arg-type]
            az = math.radians(float(azimuth))
            # simple projected chainage (not CRS-true)
            items = []
            for s in self._items:
                la, lo, _ = s.coords
                if not np.isfinite(la) or not np.isfinite(lo):
                    continue
                dy = (la - la0) * 111_000.0
                dx = (lo - lo0) * 111_000.0 * math.cos(
                    math.radians(la0)
                )
                ch = dx * math.cos(az) + dy * math.sin(az)
                items.append((ch, s))
            items.sort(key=lambda t: t[0])
            return {"origin": origin, "azimuth": azimuth,
                    "sites": [s for _, s in items]}
