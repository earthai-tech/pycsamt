# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import (
    Any, Callable, Iterator, List, Optional, Sequence,
    Tuple, Union,
)
from collections import namedtuple
from pathlib import Path
import copy
import math
import re

import numpy as np

from ..seg.edi import EDIFile  
from ..seg.collection import EDICollection  
from ..seg.heads import Head  

from ..core.base import ensure_station  

__all__ = [
    "is_pathlike", "is_edi_file", "is_edi_collection",
    "iter_edifiles", "as_edicollection",
    "station_name", "set_station_name",
    "get_coords", "set_coords",
    "maybe_copy", "apply_inplace",
    "get_freq", "freq_match", "freq_select",
    "select_by_name", "match_name",
    "wrap_azimuth", "deg_to_mrad", "mrad_to_deg",
]

_Coord = namedtuple("_Coord", ["lat", "lon", "elev"])


def is_pathlike(x: Any) -> bool:
    return isinstance(x, (str, Path))


def is_edi_file(x: Any) -> bool:
    if x is None:
        return False
    if hasattr(x, "get_section") and hasattr(x, "Z"):
        return True
    return False


def is_edi_collection(x: Any) -> bool:
    if x is None:
        return False
    if isinstance(x, EDICollection):
        return True
    if hasattr(x, "__iter__") and not is_pathlike(x):
        try:
            it = iter(x)  # type: ignore
            first = next(it, None)
        except Exception:
            return False
        return is_edi_file(first)
    return False


def iter_edifiles(edic: Any) -> Iterator[EDIFile]:
    if is_edi_file(edic):
        yield edic  # type: ignore
        return
    if hasattr(edic, "__iter__") and not is_pathlike(edic):
        for it in edic:  # type: ignore
            if is_edi_file(it):
                yield it  # type: ignore


def as_edicollection(edic: Any) -> Optional[EDICollection]:
    if isinstance(edic, EDICollection):
        return edic
    items = list(iter_edifiles(edic))
    if not items:
        return None
    try:
        return EDICollection(items=items, verbose=0)  
    except TypeError:
        return EDICollection(items, verbose=0)  


def _get_head(ed: EDIFile) -> Any:
    try:
        return ed.get_section("head")  # type: ignore
    except Exception:
        return None


def _ensure_head(ed: EDIFile) -> Any:
    h = _get_head(ed)
    if h is not None:
        return h
    try:
        h = Head()
    except Exception:
        return None
    try:
        ed.set_section("head", h)  # type: ignore
    except Exception:
        try:
            setattr(ed, "Head", h)
        except Exception:
            pass
    return h


def station_name(ed: Any) -> str:
    if hasattr(ed, "station") and getattr(ed, "station"):
        return str(getattr(ed, "station"))
    h = _get_head(ed)
    if h is not None and getattr(h, "dataid", None):
        return str(getattr(h, "dataid"))
    for k in ("name", "site", "dataid"):
        v = getattr(ed, k, None)
        if v:
            return str(v)
    return ""


def set_station_name(
    ed: Any,
    name: Optional[str] = None,
    *,
    station_id: Optional[Union[str, int]] = None,
    policy: Any = None,
    inplace: bool = True,
) -> Any:
    tgt = ed if inplace else maybe_copy(ed)
    nm = ensure_station(name, station_id, policy=policy)
    try:
        setattr(tgt, "station", nm)
    except Exception:
        pass
    try:
        h = _ensure_head(tgt)
        if h is not None:
            h.dataid = nm
    except Exception:
        pass
    return tgt


def get_coords(ed: Any) -> _Coord:
    h = _get_head(ed)
    la = getattr(h, "lat", float("nan")) if h else float("nan")
    lo = getattr(h, "long", float("nan")) if h else float("nan")
    ev = getattr(h, "elev", float("nan")) if h else float("nan")
    try:
        return _Coord(float(la), float(lo), float(ev))
    except Exception:
        return _Coord(float("nan"), float("nan"), float("nan"))


def set_coords(
    ed: Any,
    *,
    lat: Optional[float] = None,
    lon: Optional[float] = None,
    elev: Optional[float] = None,
    inplace: bool = True,
) -> Any:
    tgt = ed if inplace else maybe_copy(ed)
    h = _ensure_head(tgt)
    if h is None:
        return tgt
    if lat is not None:
        try:
            h.lat = float(lat)
        except Exception:
            pass
    if lon is not None:
        try:
            h.long = float(lon)
        except Exception:
            pass
    if elev is not None:
        try:
            h.elev = float(elev)
        except Exception:
            pass
    return tgt


def maybe_copy(x: Any) -> Any:
    try:
        return copy.deepcopy(x)
    except:
        return x


def apply_inplace(
    x: Any,
    fn: Callable[[Any], Any],
    *,
    inplace: bool = False,
) -> Any:
    if inplace:
        return fn(x)
    y = maybe_copy(x)
    return fn(y)


def get_freq(ed: Any) -> np.ndarray:
    z = getattr(ed, "Z", None)
    if z is None:
        return np.asarray([], float)
    for k in ("freq", "_freq"):
        try:
            f = getattr(z, k)
            if f is None:
                continue
            a = np.asarray(f, float)
            return a
        except Exception:
            pass
    return np.asarray([], float)


def freq_match(
    f: np.ndarray,
    target: Union[float, Sequence[float]],
    *,
    tol: float = 1e-6,
) -> np.ndarray:
    a = np.asarray(f, float)
    t = np.atleast_1d(np.asarray(target, float))
    mask = np.zeros(a.shape, dtype=bool)
    for v in t:
        mask |= np.isfinite(a) & (np.abs(a - v) <= tol)
    idx = np.where(mask)[0]
    return idx.astype(int)


def freq_select(
    f: np.ndarray,
    sel: Union[
        slice, float, Sequence[float], Tuple[float, float],
    ],
    *,
    tol: float = 1e-6,
) -> np.ndarray:
    a = np.asarray(f, float)
    if isinstance(sel, slice):
        lo = float(sel.start) if sel.start else -np.inf
        hi = float(sel.stop) if sel.stop else np.inf
        m = (a >= lo - tol) & (a <= hi + tol)
        return np.where(m)[0].astype(int)
    if isinstance(sel, (tuple, list)) and len(sel) == 2:
        lo, hi = float(sel[0]), float(sel[1])
        m = (a >= lo - tol) & (a <= hi + tol)
        return np.where(m)[0].astype(int)
    if isinstance(sel, (float, int)):
        return freq_match(a, float(sel), tol=tol)
    try:
        return freq_match(a, list(sel), tol=tol)  # type: ignore
    except Exception:
        return np.asarray([], int)


def match_name(
    pat: Union[str, re.Pattern[str], Callable[[str], bool]],
    name: str,
) -> bool:
    try:
        if callable(pat):
            return bool(pat(name))
        if isinstance(pat, re.Pattern):
            return bool(pat.search(name))
        if "*" in pat or "?" in pat or "[" in pat:
            rx = re.compile(
                "^" + re.escape(pat)
                .replace("\\*", ".*")
                .replace("\\?", ".") + "$",
                flags=re.IGNORECASE,
            )
            return bool(rx.match(name))
        return name.upper() == str(pat).upper()
    except Exception:
        return False


def select_by_name(
    edic: Any,
    pat: Union[str, re.Pattern[str], Callable[[str], bool]],
) -> List[EDIFile]:
    out: List[EDIFile] = []
    for ed in iter_edifiles(edic):
        nm = station_name(ed)
        if match_name(pat, nm):
            out.append(ed)
    return out


def wrap_azimuth(az: float) -> float:
    a = float(az) % 360.0
    return a if a >= 0 else a + 360.0


def deg_to_mrad(x: Union[float, np.ndarray]) -> np.ndarray:
    a = np.asarray(x, float)
    return a * (math.pi / 180.0) * 1000.0


def mrad_to_deg(x: Union[float, np.ndarray]) -> np.ndarray:
    a = np.asarray(x, float)
    return a * (180.0 / math.pi) / 1000.0
