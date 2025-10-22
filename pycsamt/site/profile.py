# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Iterable, List, Tuple
import math
import numpy as np

from ..constants import _M_PER_DEG
from .location import (
    Coord,
    chainage_along,
)


__all__ = [
    "infer_line_orientation",
    "Profile",
]


def _finite(x: float) -> bool:
    return math.isfinite(float(x))


def _station_name(obj: Any) -> str:
    for k in ("station", "name", "dataid", "site", "sitename"):
        try:
            v = getattr(obj, k, None)
            if v:
                return str(v)
        except Exception:
            pass
        try:
            sec = obj.get_section("head")  # type: ignore
            v = getattr(sec, k, None)
            if v:
                return str(v)
        except Exception:
            pass
    return ""


def _coords(obj: Any) -> Tuple[float, float]:
    # EDI HEAD typically has lat, long (or lon)
    la = None
    lo = None
    try:
        h = obj.get_section("head")  # type: ignore
    except Exception:
        h = getattr(obj, "HEAD", None)
    for k in ("lat",):
        try:
            la = getattr(h, k, la)
        except Exception:
            pass
    for k in ("lon", "long"):
        try:
            lo = getattr(h, k, lo)
        except Exception:
            pass
    try:
        la = float(la)
    except Exception:
        la = float("nan")
    try:
        lo = float(lo)
    except Exception:
        lo = float("nan")
    return la, lo


def _iter_sites(
    sites: Iterable[Any],
) -> List[Tuple[str, float, float, Any]]:
    out: List[Tuple[str, float, float, Any]] = []
    for s in sites:
        ed = getattr(s, "edi", s)
        name = _station_name(ed) or _station_name(s)
        la, lo = _coords(ed)
        out.append((name, la, lo, ed))
    return out


def _xy_local(
    lats: np.ndarray, lons: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    la0 = float(np.nanmean(lats))
    lo0 = float(np.nanmean(lons))
    dy = (lats - la0) * _M_PER_DEG
    dx = (lons - lo0) * _M_PER_DEG * math.cos(
        math.radians(la0)
    )
    return dx, dy


def infer_line_orientation(sites: Iterable[Any]) -> float:
    items = _iter_sites(sites)
    if not items:
        return 0.0
    la = np.asarray([i[1] for i in items], float)
    lo = np.asarray([i[2] for i in items], float)
    m = np.isfinite(la) & np.isfinite(lo)
    if not np.any(m):
        return 0.0
    x, y = _xy_local(la[m], lo[m])
    if x.size < 2:
        return 0.0
    X = np.stack([x, y], axis=0)
    C = np.cov(X)
    vals, vecs = np.linalg.eig(C)
    k = int(np.argmax(vals))
    vx, vy = float(vecs[0, k]), float(vecs[1, k])
    ang = math.degrees(math.atan2(vy, vx))
    return (ang + 360.0) % 360.0

@dataclass
class Profile:
    origin: Coord
    azimuth: float
    chainages: Dict[str, float] = field(default_factory=dict)
    spacing_stats: Dict[str, float] = field(default_factory=dict)
    gaps: List[Tuple[float, float]] = field(default_factory=list)

    # build from sites
    @classmethod
    def from_sites(
        cls,
        sites: Iterable[Any],
        *,
        origin: Coord | None = None,
        azimuth: float | None = None,
    ) -> "Profile":
        items = _iter_sites(sites)
        if not items:
            o = origin or Coord(0.0, 0.0, 0.0)
            return cls(o, float(azimuth or 0.0))
        # origin default = first finite
        if origin is None:
            for _, la, lo, _ in items:
                if _finite(la) and _finite(lo):
                    origin = Coord(la, lo, 0.0)
                    break
        if origin is None:
            origin = Coord(0.0, 0.0, 0.0)
        if azimuth is None:
            azimuth = infer_line_orientation(
                [it[3] for it in items]
            )
        s: Dict[str, float] = {}
        for name, la, lo, _ in items:
            if _finite(la) and _finite(lo):
                s[name] = float(
                    chainage_along(
                        (origin.lat, origin.lon),
                        float(azimuth),
                        (la, lo),
                    )
                )
        prof = cls(origin, float(azimuth), s)
        prof._update_stats()
        return prof

    # sort by chainage
    def sort_sites(
        self,
        sites: Iterable[Any],
    ) -> List[Any]:
        items = _iter_sites(sites)
        order = []
        for it in items:
            si = self.chainages.get(it[0], float("nan"))
            order.append((si, it[3]))
        order.sort(key=lambda t: (not _finite(t[0]), t[0]))
        return [o[1] for o in order if _finite(o[0])]

    # window chainages
    def slice(
        self,
        s_min: float,
        s_max: float,
    ) -> Dict[str, float]:
        out = {
            k: v for k, v in self.chainages.items()
            if _finite(v) and (s_min <= v <= s_max)
        }
        return dict(sorted(out.items(), key=lambda t: t[1]))

    # regular grid along min..max
    def resample(self, step: float) -> np.ndarray:
        if not self.chainages:
            return np.asarray([], float)
        s = np.asarray(list(self.chainages.values()), float)
        s = s[np.isfinite(s)]
        if s.size == 0:
            return np.asarray([], float)
        a = float(np.nanmin(s))
        b = float(np.nanmax(s))
        if step <= 0:
            return np.asarray([], float)
        n = int(max(1, math.floor((b - a) / step) + 1))
        return a + step * np.arange(n, dtype=float)

    # short text summary
    def summary(self) -> Dict[str, float]:
        d: Dict[str, float] = {}
        d.update(self.spacing_stats)
        d["n_sites"] = float(len(self.chainages))
        if self.chainages:
            vals = np.asarray(
                list(self.chainages.values()), float
            )
            m = np.isfinite(vals)
            if np.any(m):
                d["s_min"] = float(np.min(vals[m]))
                d["s_max"] = float(np.max(vals[m]))
        d["n_gaps"] = float(len(self.gaps))
        return d

    # internal stats
    def _update_stats(self) -> None:
        if not self.chainages:
            self.spacing_stats = {}
            self.gaps = []
            return
        s = np.asarray(list(self.chainages.values()), float)
        s = s[np.isfinite(s)]
        s.sort()
        if s.size < 2:
            self.spacing_stats = {
                "spacing_mean": float("nan"),
                "spacing_med": float("nan"),
                "spacing_min": float("nan"),
                "spacing_max": float("nan"),
            }
            self.gaps = []
            return
        d = np.diff(s)
        self.spacing_stats = {
            "spacing_mean": float(np.mean(d)),
            "spacing_med": float(np.median(d)),
            "spacing_min": float(np.min(d)),
            "spacing_max": float(np.max(d)),
        }
        med = float(np.median(d))
        thr = 1.5 * med if med > 0 else float("inf")
        gaps: List[Tuple[float, float]] = []
        for i in range(d.size):
            if d[i] > thr:
                gaps.append((float(s[i]), float(s[i + 1])))
        self.gaps = gaps
