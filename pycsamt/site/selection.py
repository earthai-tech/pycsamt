# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Callable, Iterable, List
import re
import math

import numpy as np

from .utils import (
    as_edicollection, iter_edifiles, station_name,
    get_freq, get_coords,
)

__all__ = [
    "by_names", "by_index", "by_chainage", "by_freq",
    "by_bbox", "by_predicate",
    "keep_finite_z", "mask_large_phase_err", "drop_empty",
]


def _to_sites(x: Any):
    """Coerce any edi-like into a Sites wrapper."""
    from .base import Sites  # lazy to avoid cycles
    if isinstance(x, Sites):
        return x
    coll = as_edicollection(x) or []
    try:
        return Sites(coll)
    except TypeError:
        return Sites(edic=coll)


def _new_sites(src: Any, items: List[Any]):
    """Build a new Sites preserving wrapper semantics."""
    from .base import Sites  # lazy import
    try:
        return Sites(items)
    except TypeError:
        return Sites(edic=items)


def _in_box(lat: float, lon: float,
            mnla: float, mnlo: float,
            mxla: float, mxlo: float) -> bool:
    if not (math.isfinite(lat) and math.isfinite(lon)):
        return False
    return (mnla <= lat <= mxla) and (mnlo <= lon <= mxlo)


def _name_matches(nm: str, pat, case: bool) -> bool:
    if callable(pat):
        try:
            return bool(pat(nm))
        except Exception:
            return False
    if isinstance(pat, re.Pattern):
        return bool(pat.search(nm))
    s = str(pat)
    if not case:
        nm = nm.upper()
        s = s.upper()
    if ("*" in s) or ("?" in s) or ("[" in s):
        rx = re.compile(
            "^" + re.escape(s).replace("\\*", ".*")
            .replace("\\?", ".") + "$"
        )
        return bool(rx.match(nm))
    return nm == s


def _any_finite_z(ed: Any) -> bool:
    z = getattr(ed, "Z", None)
    if z is None:
        return False
    arr = getattr(z, "_z", None) or getattr(z, "z", None)
    if arr is None:
        # fallback to resistivity if present
        r = getattr(z, "_resistivity", None)
        if r is None:
            return False
        a = np.asarray(r, float)
        return np.isfinite(a).any()
    a = np.asarray(arr, complex)
    return np.isfinite(a.real).any() or np.isfinite(a.imag).any()


def _max_phase_err(ed: Any) -> float:
    z = getattr(ed, "Z", None)
    if z is None:
        return float("nan")
    pe = getattr(z, "_phase_err", None)
    if pe is None:
        return float("nan")
    a = np.asarray(pe, float)
    if a.size == 0:
        return float("nan")
    with np.errstate(all="ignore"):
        return float(np.nanmax(a))


def _is_empty_site(ed: Any) -> bool:
    f = get_freq(ed)
    if f.size == 0:
        return True
    z = getattr(ed, "Z", None)
    if z is None:
        return True
    has = any(
        getattr(z, k, None) is not None
        for k in ("_z", "z", "_resistivity")
    )
    return not has


def by_names(
    sites: Any,
    patterns: Iterable[Any] | Any,
    *,
    case: bool = False,
):
    s = _to_sites(sites)
    pats = (patterns if isinstance(patterns, (list, tuple))
            else [patterns])
    keep = []
    for ed in iter_edifiles(s.edic):
        nm = station_name(ed)
        ok = any(_name_matches(nm, p, case) for p in pats)
        if ok:
            keep.append(ed)
    return _new_sites(s, keep)


def by_index(sites: Any, indices: Iterable[int] | int):
    s = _to_sites(sites)
    idxs = (indices if isinstance(indices, (list, tuple))
            else [indices])
    n = len(list(iter_edifiles(s.edic)))
    idxs = [i for i in idxs if isinstance(i, int) and -n <= i < n]
    if not idxs:
        return _new_sites(s, [])
    items = [e for i, e in enumerate(iter_edifiles(s.edic))
             if i in set(idxs)]
    return _new_sites(s, items)


def by_chainage(
    sites: Any,
    smin: float,
    smax: float,
):
    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        h = getattr(ed, "get_section", lambda *_: None)("head")
        ch = getattr(h, "chainage", None) if h else None
        if ch is None:
            ch = getattr(ed, "chainage", None)
        try:
            v = float(ch)
        except Exception:
            continue
        if smin <= v <= smax:
            keep.append(ed)
    return _new_sites(s, keep)


def by_freq(
    sites: Any,
    fmin: float,
    fmax: float,
):
    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        f = get_freq(ed)
        if f.size == 0:
            continue
        m = np.isfinite(f) & (f >= fmin) & (f <= fmax)
        if bool(m.any()):
            keep.append(ed)
    return _new_sites(s, keep)


def by_bbox(
    sites: Any,
    minlat: float,
    minlon: float,
    maxlat: float,
    maxlon: float,
):
    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        c = get_coords(ed)
        if _in_box(c.lat, c.lon, minlat, minlon, maxlat, maxlon):
            keep.append(ed)
    return _new_sites(s, keep)


def by_predicate(
    sites: Any,
    pred: Callable[[Any], bool],
):
    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        try:
            if pred(ed):
                keep.append(ed)
        except Exception:
            pass
    return _new_sites(s, keep)


def keep_finite_z(sites: Any):
    s = _to_sites(sites)
    keep = [ed for ed in iter_edifiles(s.edic)
            if _any_finite_z(ed)]
    return _new_sites(s, keep)


def mask_large_phase_err(
    sites: Any,
    thresh: float,
):
    s = _to_sites(sites)
    keep = []
    for ed in iter_edifiles(s.edic):
        m = _max_phase_err(ed)
        if not math.isfinite(m):
            # no err → keep
            keep.append(ed)
        elif m <= float(thresh):
            keep.append(ed)
    return _new_sites(s, keep)


def drop_empty(sites: Any):
    s = _to_sites(sites)
    keep = [ed for ed in iter_edifiles(s.edic)
            if not _is_empty_site(ed)]
    return _new_sites(s, keep)
