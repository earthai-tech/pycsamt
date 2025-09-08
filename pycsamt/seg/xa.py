# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0-or-later

from __future__ import annotations

from pathlib import Path
from typing import Iterable, List, Dict

import numpy as np
import xarray as xr

from ..log.logger import get_logger
from .edi import EDIFile

logger = get_logger(__name__)

__all__ = ["XAMixin",  "EDIAcc", "build_dataset"]


def _site_id(ed: EDIFile) -> str:
    # prefer station; else file stem; else "site"
    sid = getattr(ed, "station", None)
    if sid:
        return str(sid)
    p = getattr(ed, "path", None)
    if isinstance(p, Path):
        return p.stem
    return "site"

def _meta(ed: EDIFile) -> Dict[str, object]:
    head = ed.get_section("head")
    info = ed.get_section("info")
    p = getattr(ed, "path", None)

    def _get(obj, name, dv=None):
        return getattr(obj, name, dv) if obj else dv

    # --- robust software extraction ---
    software = None
    proc = getattr(info, "Processing", None)
    if proc is not None and hasattr(proc, "ProcessingSoftware"):
        sw = getattr(proc, "ProcessingSoftware")
        # accept either an object with .name or a plain string
        software = getattr(sw, "name", sw)

    return {
        "path": str(p) if isinstance(p, Path) else None,
        "filename": p.name if isinstance(p, Path) else None,
        "dataid": _get(head, "dataid", None),
        "lat": _get(head, "lat", None),
        "lon": _get(head, "long", None) or _get(head, "lon", None),
        "elev": _get(head, "elev", None),
        "has_tip": bool(
            getattr(ed.Tip, "tipper", None) is not None
            and getattr(ed.Tip.tipper, "size", 0) > 0
        ),
        "nfreq": int(getattr(ed.Z, "n_freq", 0) or 0),
        "has_spec": ed.get_section("spectra") is not None,
        "has_ts": ed.get_section("timeseries") is not None,
        "software": software,
    }


def _pad2d(
    seqs: List[np.ndarray],
    *,
    pad: float = np.nan,
) -> np.ndarray:
    n = len(seqs)
    m = max((s.size for s in seqs), default=0)
    out = np.full((n, m), pad, float)
    for i, a in enumerate(seqs):
        if a.size:
            out[i, : a.size] = a
    return out

def _spec_pack(ed: EDIFile) -> Dict[str, xr.DataArray]:
    spec = ed.get_section("spectra")
    if spec is None:
        return {}

    # Try direct attributes first
    freq = getattr(spec, "freq", None)
    vals = getattr(spec, "values", None)
    bw = getattr(spec, "bw", None)
    avgt = getattr(spec, "avgt", None)
    rs = getattr(spec, "rotspec", None)

    if freq is None or vals is None:
        try:
            sect, io = spec.to_io()  # may not exist on all impls
            blks = getattr(io, "blocks", [])
            # prefer block.options['freq'] if present
            freq = [
                (getattr(b, "freq", None)
                 if getattr(b, "freq", None) is not None
                 else b.options.get("freq", None))
                for b in blks
            ]
            vals = [
                np.asarray(getattr(b, "values", []), float)
                for b in blks
            ]
            bw = [getattr(b, "bw", None) for b in blks]
            avgt = [getattr(b, "avgt", None) for b in blks]
            rs = [getattr(b, "rotspec", None) for b in blks]
        except Exception:
            return {}  # don’t break the whole dataset

    # Sanitize freq
    try:
        freq_clean = [float(x) for x in freq if x is not None]
    except Exception:
        return {}  # invalid freq list

    f = np.asarray(freq_clean, float)
    if f.size == 0:
        return {}

    # Values: ragged → padded
    try:
        vlist = [np.asarray(v, float) for v in vals]
    except Exception:
        return {}

    spec_vals = _pad2d(vlist, pad=np.nan)
    spec_len = np.array([v.size for v in vlist], int)
    sidx = np.arange(spec_vals.shape[1], dtype=int)

    out: Dict[str, xr.DataArray] = {
        "spec_vals": xr.DataArray(
            spec_vals, dims=("freq", "sidx"),
            coords={"freq": f, "sidx": sidx},
        ),
        "spec_len": xr.DataArray(
            spec_len, dims=("freq",), coords={"freq": f}
        ),
    }

    # Optional per-freq metadata (only if lengths match)
    def _ok_1d(a):
        try:
            a = np.asarray(a, float)
            return a.size == f.size, a
        except Exception:
            return False, None

    for name, arr in (
        ("spec_bw", bw),
        ("spec_avgt", avgt),
        ("spec_rotspec", rs),
    ):
        ok, arrv = _ok_1d(arr)
        if ok:
            out[name] = xr.DataArray(
                arrv, dims=("freq",), coords={"freq": f}
            )

    return out


def _ts_pack(ed: EDIFile) -> Dict[str, xr.DataArray]:
    ts = ed.get_section("timeseries")
    if ts is None:
        return {}

    # channels, dt per channel, sequences
    try:
        ch = list(ts.channels())
    except Exception:
        return {}

    if not ch:
        return {}

    series = [np.asarray(ts.get(c), float) for c in ch]
    npts = [a.size for a in series]
    mat = _pad2d(series, pad=np.nan)  # (ch_len→columns)
    # we want (sample, ch)
    mat = mat.T

    dt = np.array([float(ts.dt_map.get(c, np.nan)) for c in ch])
    nmax = mat.shape[0]
    # time grid per column with padding as nan
    T = np.full((nmax, len(ch)), np.nan, float)
    for j, (a, d) in enumerate(zip(series, dt)):
        if np.isfinite(d) and a.size:
            T[: a.size, j] = np.arange(a.size, float) * d

    out: Dict[str, xr.DataArray] = {}
    out["ts"] = xr.DataArray(
        mat,
        dims=("sample", "ch"),
        coords={"sample": np.arange(nmax), "ch": ch},
    )
    out["time"] = xr.DataArray(
        T,
        dims=("sample", "ch"),
        coords={"sample": np.arange(nmax), "ch": ch},
    )
    out["dt"] = xr.DataArray(
        dt, dims=("ch",), coords={"ch": ch}
    )
    out["npts"] = xr.DataArray(
        np.array(npts, int), dims=("ch",), coords={"ch": ch}
    )
    return out


# def _ds_from_edi(ed: EDIFile) -> "xr.Dataset":
#     sid = _site_id(ed)
#     f = np.asarray(ed.Z.freq, float)
#     if f.ndim != 1 or f.size == 0:
#         f = np.asarray([], float)

def _ds_from_edi(ed: EDIFile) -> xr.Dataset:
    sid = _site_id(ed)

    # prefer Z freq; if empty, fall back to Spectra freq
    f = np.asarray(getattr(ed.Z, "freq", []), float)
    if f.ndim != 1 or f.size == 0:
        sp = ed.get_section("spectra")
        sf = np.asarray(
            getattr(sp, "freq", []), float
            ) if sp is not None else np.asarray([], float)
        f = sf if (sf.ndim == 1 and sf.size > 0) else np.asarray([], float)

    # Z as complex (freq, i, j)
    z = np.asarray(getattr(ed.Z, "z", None))
    if z is None or z.size == 0:
        z = np.zeros((f.size, 2, 2), complex)

    ze = np.asarray(getattr(ed.Z, "z_err", None))
    if ze is None or ze.size == 0:
        ze = np.zeros((f.size, 2, 2), float)

    zrot = np.asarray(
        getattr(ed.Z, "rotation_angle", np.zeros(f.size))
    ).astype(float)
    if zrot.size != f.size:
        zrot = np.zeros(f.size, float)

    # tipper (freq, 1, 2) -> (freq, comp)
    tip = getattr(ed.Tip, "tipper", None)
    if tip is None or getattr(tip, "size", 0) == 0:
        tip = np.zeros((f.size, 1, 2), complex)
    tip = np.asarray(tip)
    tip_da = tip[:, 0, :] if tip.ndim == 3 else \
        np.zeros((f.size, 2), complex)

    terr = getattr(ed.Tip, "_tipper_err", None)
    if terr is None:
        terr = np.zeros_like(tip, float)
    terr = np.asarray(terr)
    terr_da = terr[:, 0, :] if terr.ndim == 3 else \
        np.zeros((f.size, 2), float)

    # optional scalar forms
    def _safe(a) -> np.ndarray:
        v = np.asarray(a) if a is not None else None
        if v is None or v.size == 0:
            return np.zeros((f.size, 2, 2), float)
        return v

    rho = np.stack(
        [
            _safe(ed.Z.res_xx),
            _safe(ed.Z.res_xy),
            _safe(ed.Z.res_yx),
            _safe(ed.Z.res_yy),
        ],
        axis=1,
    ).reshape(f.size, 2, 2)

    phi = np.stack(
        [
            _safe(ed.Z.phase_xx),
            _safe(ed.Z.phase_xy),
            _safe(ed.Z.phase_yx),
            _safe(ed.Z.phase_yy),
        ],
        axis=1,
    ).reshape(f.size, 2, 2)

    # base TF dataset
    ds = xr.Dataset(
        data_vars={
            "z": (("freq", "i", "j"), z),
            "z_err": (("freq", "i", "j"), ze),
            "zrot": (("freq",), zrot),
            "tip": (("freq", "tcomp"), tip_da),
            "tip_err": (("freq", "tcomp"), terr_da),
            "rho": (("freq", "i", "j"), rho),
            "phi": (("freq", "i", "j"), phi),
        },
        coords={
            "freq": f,
            "i": [0, 1],
            "j": [0, 1],
            "tcomp": ["tx", "ty"],
        },
        attrs=_meta(ed),
    ).expand_dims(site=[sid])

    # spectra/time-series if present
    sp = _spec_pack(ed)
    if sp:
        for k, v in sp.items():
            ds[k] = v

    ts = _ts_pack(ed)
    if ts:
        for k, v in ts.items():
            ds[k] = v

    return ds


def build_dataset(
    edis: Iterable[EDIFile],
    *,
    drop_empty: bool = True,
) -> "xr.Dataset":
    """
    Build a multi-site Dataset from an iterable of EDIFile.
    """
    parts: List[xr.Dataset] = []
    for ed in edis:
        
        try:
            ds = _ds_from_edi(ed)
            # consider dataset "empty" only if it has 
            # neither TF freq nor spectra/ts vars
            is_empty_tf = ds.sizes.get("freq", 0) == 0
            has_sp_or_ts = any(k in ds.data_vars for k in ("spec_vals", "ts"))
            if drop_empty and is_empty_tf and not has_sp_or_ts:
                continue  # truly empty
            parts.append(ds)
        except Exception as exc:
            logger.debug("skip %s: %s", ed, exc)
        
        # try:
        #     ds = _ds_from_edi(ed)
        #     if drop_empty and ds.sizes.get("freq", 0) == 0:
        #         parts.append(ds)
        #     else:
        #         parts.append(ds)
        # except Exception as exc:
        #     logger.debug("skip %s: %s", ed, exc)

    if not parts:
        # empty shell for downstream code
        base = xr.Dataset(
            data_vars={},
            coords={
                "site": [],
                "freq": [],
                "i": [0, 1],
                "j": [0, 1],
                "tcomp": ["tx", "ty"],
            },
            attrs={},
        )
        return base

    # align by coords, join outer on freq
    ds = xr.concat(parts, dim="site", join="outer")

    # component labels for convenience
    ds = ds.assign_coords(
        c=("c", ["zxx", "zxy", "zyx", "zyy"])
    )
    return ds


class XAMixin:
    """
    Add xarray export to a host collection-like object.
    """

    def to_xarray(
        self,
        *,
        drop_empty: bool = True,
    ) -> "xr.Dataset":
        # host must be iterable over EDIFile
        return build_dataset(self, drop_empty=drop_empty)

    def meta_table(self) -> "xr.Dataset":
        """
        Return a tiny Dataset with per-site metadata only.
        """
        rows: List[Dict[str, object]] = []
        for ed in self:
            m = _meta(ed)
            m["site"] = _site_id(ed)
            rows.append(m)
        if not rows:
            return xr.Dataset(coords={"site": []})
        cols = {}
        for k in rows[0].keys():
            cols[k] = ("site", [r.get(k) for r in rows])
        ds = xr.Dataset(coords={"site": cols.pop("site")[1]})
        for k, v in cols.items():
            ds[k] = v
        return ds


@xr.register_dataset_accessor("edi")
class EDIAcc:
    def __init__(self, ds: "xr.Dataset") -> None:
        self._ds = ds

    @property
    def stations(self) -> List[str]:
        s = self._ds.coords.get("site", None)
        return [] if s is None else [str(v) for v in s.data]

    def get(self, site: str) -> "xr.Dataset":
        return self._ds.sel(site=str(site))

    def components(self) -> List[str]:
        if "c" in self._ds.coords:
            return [str(v) for v in self._ds["c"].data]
        return ["zxx", "zxy", "zyx", "zyy"]

    def z_as_comp(self) -> "xr.DataArray":
        # reshape (i,j)->comp with mapping
        z = self._ds["z"]
        zz = xr.concat(
            [
                z.sel(i=0, j=0),
                z.sel(i=0, j=1),
                z.sel(i=1, j=0),
                z.sel(i=1, j=1),
            ],
            dim="c",
        )
        return zz.assign_coords(
            c=["zxx", "zxy", "zyx", "zyy"]
        )

    def band(
        self,
        fmin: float | None = None,
        fmax: float | None = None,
    ) -> "xr.Dataset":
        ds = self._ds
        if "freq" not in ds.coords:
            return ds
        cond = np.ones(ds.sizes["freq"], bool)
        fv = ds["freq"].data
        if fmin is not None:
            cond &= fv >= float(fmin)
        if fmax is not None:
            cond &= fv <= float(fmax)
        return ds.isel(freq=np.where(cond)[0])

    # --- spectra helpers ---
    def has_spectra(self) -> bool:
        return "spec_vals" in self._ds

    def spectra(self) -> "xr.Dataset":
        keep = [k for k in self._ds.data_vars
                if k.startswith("spec_")]
        return self._ds[keep] if keep else self._ds

    # --- time-series helpers ---
    def has_timeseries(self) -> bool:
        return "ts" in self._ds

    def timeseries(self) -> "xr.Dataset":
        keep = [k for k in ("ts", "time", "dt", "npts")]
        keep = [k for k in keep if k in self._ds]
        return self._ds[keep] if keep else self._ds

    def attrs(self) -> Dict[str, object]:
        return dict(self._ds.attrs)
