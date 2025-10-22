# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import Any, Callable, Iterable, Optional, Tuple
import math
import numpy as np

from .utils import (
    iter_edifiles, copy_edi, station_name, 
    set_station_name, get_freq,
    ensure_head, set_coords as _set_coords,
    safe_compute_res_phase,
)


__all__ = [
    # single-site mutations
    "rotate", "select_freq", "rename", "set_coords",
    "fill_missing", "recompute_res_phase",
    # broadcast over Sites/collections
    "rotate_all", "select_freq_all", "rename_all",
    "set_coords_all",
]


def _to_site(ed: Any, *, inplace: bool) -> Any:
    """Return the object we will mutate."""
    if inplace:
        return ed
    return copy_edi(ed)


def _rotm(angle_deg: float) -> Tuple[np.ndarray, np.ndarray]:
    """2D rotation and its inverse for MT tensors."""
    th = float(angle_deg) * (math.pi / 180.0)
    c = math.cos(th)
    s = math.sin(th)
    r = np.array([[c, s], [-s, c]], dtype=float)
    rinv = np.array([[c, -s], [s, c]], dtype=float)
    return r, rinv


def _apply_freq_subset(ed: Any, sel) -> None:
    """
    Apply row selection to all Z/Tip arrays present.

    `sel` can be a boolean mask or integer indices.
    """
    z = getattr(ed, "Z", None)
    if z is not None:
        for k in ("_freq", "_z", "_z_err",
                  "_resistivity", "_phase", "_phase_err"):
            a = getattr(z, k, None)
            if a is None:
                continue
            arr = np.asarray(a)
            try:
                setattr(z, k, arr[sel])
            except Exception:
                pass

    tip = getattr(ed, "Tip", None)
    if tip is not None:
        for k in ("_freq", "_tipper", "_tipper_err",
                  "_amp", "_phase", "_azimuth"):
            a = getattr(tip, k, None)
            if a is None:
                continue
            arr = np.asarray(a)
            try:
                setattr(tip, k, arr[sel])
            except Exception:
                pass


def _fill_array(a: Optional[np.ndarray],
                shape: Tuple[int, ...],
                how: str) -> np.ndarray:
    if a is None:
        if how == "zero":
            return np.zeros(shape, float)
        return np.full(shape, np.nan, float)
    out = np.array(a, copy=True)
    if how == "zero":
        out = np.nan_to_num(out, nan=0.0,
                            posinf=0.0, neginf=0.0)
    else:
        out[~np.isfinite(out)] = np.nan
    return out


def rotate(site: Any, angle_deg: float,
           *, inplace: bool = False) -> Any:
    """Rotate Z (and Tipper if present) by `angle_deg`."""
    ed = _to_site(site, inplace=inplace)

    z = getattr(ed, "Z", None)
    if z is not None and getattr(z, "_z", None) is not None:
        r, rinv = _rotm(angle_deg)
        zz = np.asarray(z._z, complex)
        z._z = r[None, :, :] @ zz @ rinv[None, :, :]
        # naive err rotation (best-effort)
        if getattr(z, "_z_err", None) is not None:
            ze = np.asarray(z._z_err, float)
            ar = np.abs(r)
            ari = np.abs(rinv)
            z._z_err = (
                ar[None, :, :] * ze * ari[None, :, :]
            )

    tip = getattr(ed, "Tip", None)
    if tip is not None and getattr(tip, "_tipper", None) is not None:
        _, rinv = _rotm(angle_deg)
        t = np.asarray(tip._tipper, complex)
        if t.ndim == 2 and t.shape[1] == 2:
            tip._tipper = (t @ rinv.T)

    return ed


def select_freq(
    site: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    keep: Iterable[int] | np.ndarray | None = None,
    inplace: bool = False,
) -> Any:
    """
    Subset rows by frequency range or explicit indices.
    """
    ed = _to_site(site, inplace=inplace)
    f = get_freq(ed)

    if keep is not None:
        sel = np.asarray(list(keep))
        _apply_freq_subset(ed, sel)
        return ed

    if f.size == 0:
        return ed
    m = np.isfinite(f)
    if fmin is not None:
        m &= (f >= float(fmin))
    if fmax is not None:
        m &= (f <= float(fmax))
    _apply_freq_subset(ed, m)
    return ed


def rename(
    site: Any,
    name: Optional[str] = None,
    policy: Optional[Callable[[str], str]] = None,
    *,
    inplace: bool = False,
) -> Any:
    """
    Rename station using explicit `name` or a `policy`.
    """
    ed = _to_site(site, inplace=inplace)
    cur = station_name(ed)

    if name is None and policy is None:
        return ed

    new = name if name is not None else cur
    if policy is not None:
        try:
            new = str(policy(cur))
        except Exception:
            new = str(cur)

    set_station_name(ed, new)
    # keep head.dataid in sync
    h = ensure_head(ed)
    h.dataid = str(new or "")
    return ed


def set_coords(
    site: Any,
    *,
    lat: float | None = None,
    lon: float | None = None,
    elev: float | None = None,
    inplace: bool = False,
) -> Any:
    """Set HEAD lat/lon/elev (only values provided)."""
    ed = _to_site(site, inplace=inplace)
    _set_coords(ed, lat=lat, lon=lon, elev=elev)
    return ed


def fill_missing(
    site: Any,
    *,
    how: str = "zero",
    components: Iterable[str] = ("Z", "Tip"),
    inplace: bool = False,
) -> Any:
    """
    Replace NaN/Inf or absent arrays with zeros/NaNs.
    """
    how = str(how).lower()
    if how not in {"zero", "nan"}:
        raise ValueError("how must be 'zero' or 'nan'")
    ed = _to_site(site, inplace=inplace)
    f = get_freq(ed)
    n = int(f.size)

    if "Z" in set(map(str.upper, components)):
        z = getattr(ed, "Z", None)
        if z is not None:
            if getattr(z, "_z", None) is None:
                z._z = np.zeros((n, 2, 2), complex)
            else:
                z._z = _fill_array(z._z, (n, 2, 2), how)
            if getattr(z, "_z_err", None) is not None:
                z._z_err = _fill_array(z._z_err, (n, 2, 2),
                                       how)
            if getattr(z, "_resistivity", None) is not None:
                z._resistivity = _fill_array(
                    z._resistivity, (n, 2, 2), how
                )
            if getattr(z, "_phase", None) is not None:
                z._phase = _fill_array(z._phase, (n, 2, 2),
                                       how)
            if getattr(z, "_phase_err", None) is not None:
                z._phase_err = _fill_array(
                    z._phase_err, (n, 2, 2), how
                )

    if "TIP" in set(map(str.upper, components)):
        tp = getattr(ed, "Tip", None)
        if tp is not None:
            if getattr(tp, "_tipper", None) is None:
                tp._tipper = np.zeros((n, 2), complex)
            else:
                tp._tipper = _fill_array(tp._tipper, (n, 2),
                                         how)
            if getattr(tp, "_tipper_err", None) is not None:
                tp._tipper_err = _fill_array(
                    tp._tipper_err, (n, 2), how
                )
    return ed


def recompute_res_phase(site: Any,
                        *, inplace: bool = False) -> Any:
    """
    Recompute ρ, φ from Z using backend helper.
    """
    ed = _to_site(site, inplace=inplace)
    safe_compute_res_phase(ed)
    return ed


def _each_site(edic_or_sites: Any):
    """Yield (edi, index) from Sites or any edi-like."""
    try:
        from .base import Sites  # lazy
        if isinstance(edic_or_sites, Sites):
            for i, ed in enumerate(iter_edifiles(edic_or_sites.edic)):
                yield ed, i
            return
    except Exception:
        pass
    for i, ed in enumerate(iter_edifiles(edic_or_sites)):
        yield ed, i


def rotate_all(sites: Any, angle_deg: float,
               *, inplace: bool = False):
    items = []
    for ed, _ in _each_site(sites):
        items.append(rotate(ed, angle_deg, inplace=False))
    if inplace:
        # copy back into wrapper order if possible
        try:
            for (ed, i), ne in zip(_each_site(sites), items):
                sites.edic[i] = ne
            return sites
        except Exception:
            pass
    from .base import Sites  # lazy
    return Sites(items)


def select_freq_all(
    sites: Any,
    *,
    fmin: float | None = None,
    fmax: float | None = None,
    keep: Iterable[int] | np.ndarray | None = None,
    inplace: bool = False,
):
    items = []
    for ed, _ in _each_site(sites):
        items.append(
            select_freq(ed, fmin=fmin, fmax=fmax,
                        keep=keep, inplace=False)
        )
    if inplace:
        try:
            for (ed, i), ne in zip(_each_site(sites), items):
                sites.edic[i] = ne
            return sites
        except Exception:
            pass
    from .base import Sites
    return Sites(items)


def rename_all(
    sites: Any,
    *,
    policy: Callable[[str], str] | None = None,
    name_fn: Callable[[Any], str] | None = None,
    inplace: bool = False,
):
    """
    Rename every site. Provide `policy(name)->newname`
    or `name_fn(edi)->newname`.
    """
    items = []
    for ed, _ in _each_site(sites):
        cur = station_name(ed)
        new = (name_fn(ed) if name_fn else
               (policy(cur) if policy else cur))
        items.append(rename(ed, name=new, inplace=False))

    if inplace:
        try:
            for (ed, i), ne in zip(_each_site(sites), items):
                sites.edic[i] = ne
            return sites
        except Exception:
            pass
    from .base import Sites
    return Sites(items)


def set_coords_all(
    sites: Any,
    src: Any,
    *,
    inplace: bool = False,
):
    """
    Set coordinates from a frame/mapping/callable.

    `src` may be:
      - callable(edi)->(lat, lon, elev)
      - mapping[name]->(lat, lon, elev)
      - object with `.frame` DataFrame-like:
        columns: 'station','lat','lon','elev'
    """
    def _lookup_by_frame(name: str):
        fr = getattr(src, "frame", None)
        if fr is None:
            return None
        try:
            if "station" in fr.columns:
                row = fr[fr["station"] == name]
            else:
                row = None
            if row is None or len(row) == 0:
                return None
            lat = float(row["lat"].iloc[0])
            lon = float(row["lon"].iloc[0])
            elv = float(row.get("elev", 0.0).iloc[0])
            return lat, lon, elv
        except Exception:
            return None

    items = []
    for ed, _ in _each_site(sites):
        nm = station_name(ed)
        lat = lon = elv = None
        if callable(src):
            try:
                v = src(ed)
                if v is not None:
                    lat, lon, elv = v
            except Exception:
                pass
        elif hasattr(src, "get"):
            try:
                v = src.get(nm)
                if v is not None:
                    lat, lon, elv = v
            except Exception:
                pass
        else:
            v = _lookup_by_frame(nm)
            if v is not None:
                lat, lon, elv = v
        items.append(
            set_coords(ed, lat=lat, lon=lon, elev=elv,
                       inplace=False)
        )

    if inplace:
        try:
            for (ed, i), ne in zip(_each_site(sites), items):
                sites.edic[i] = ne
            return sites
        except Exception:
            pass
    from .base import Sites
    return Sites(items)
