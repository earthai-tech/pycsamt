# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Point edits for the Edit ▸ Frequencies ▸ Point Editor (Qt-free).

Operations on one station's impedance, per frequency row and component
(``xy``, ``yx``):

* :func:`mask` -- the value becomes NaN (Z and its error): hidden from
  processing and inversion, kept as a row in the file;
* :func:`delete` -- the whole frequency row is removed (all components,
  tipper included) through :func:`pycsamt.site.edit.select_freq`;
* :func:`interpolate` -- the value is replaced by one interpolated from the
  component's valid neighbours (log-frequency; log |Z| and unwrapped
  phase), its error from theirs; clearly an estimate, and logged as one;
* :func:`restore` -- the original value of that row comes back;
* :func:`static_shift` -- a whole apparent-resistivity curve is scaled by
  a factor (Z by its square root, so phase is unchanged), the standard
  static-shift correction.

Every operation works on a copy and appends a line to the station's EDI
``INFO`` block (``PYCSAMT_POINT_EDITS=...``), so the exported file says
what was changed.  No operation moves a single measured value to an
arbitrary position.
"""

from __future__ import annotations

import copy
from typing import Any

import numpy as np

__all__ = [
    "COMPONENTS",
    "curves",
    "delete",
    "interpolate",
    "mask",
    "nearest_rows",
    "restore",
    "static_shift",
]

COMPONENTS = {"xy": (0, 1), "yx": (1, 0)}
_INFO_KEY = "PYCSAMT_POINT_EDITS"


def _Z(site):
    edi = site.edi
    return getattr(edi, "Z", None)


def _arrays(site) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    f = np.asarray(site.freq, dtype=float)
    z = np.array(site.z, dtype=complex).reshape(-1, 2, 2)
    ze = site.z_err
    ze = None if ze is None else np.array(ze, dtype=float).reshape(-1, 2, 2)
    return f, z, ze


def curves(site) -> dict:
    """Apparent resistivity and phase of each component, with errors.

    ``{"freq": f, "xy": {"rho", "phi", "rho_err", "phi_err", "valid"},
    "yx": {...}}``; rows in file order.
    """
    f, z, ze = _arrays(site)
    out = {"freq": f}
    for comp, (i, j) in COMPONENTS.items():
        zc = z[:, i, j]
        amp = np.abs(zc)
        with np.errstate(all="ignore"):
            rho = 0.2 / f * amp ** 2
            phi = np.degrees(np.angle(zc)) % 180.0
            if ze is not None:
                rel = ze[:, i, j] / amp
                rho_err = 2.0 * rho * rel
                phi_err = np.degrees(rel)
            else:
                rho_err = phi_err = np.full(f.shape, np.nan)
        out[comp] = {"rho": rho, "phi": phi, "rho_err": rho_err,
                     "phi_err": phi_err, "valid": np.isfinite(zc)}
    return out


def _write(site, z, ze) -> None:
    Z = _Z(site)
    Z.z = z
    if ze is not None:
        try:
            Z.z_err = ze
        except Exception:
            setattr(Z, "z_error", ze)


def _note(site, text: str) -> Any:
    """Append *text* to the station's INFO edit log (on a copy)."""
    from pycsamt.site.metadata import update_metadata

    old = ""
    try:
        info = site.edi.get_section("info")
        for line in getattr(info, "info_text", []) or []:
            if str(line).strip().upper().startswith(_INFO_KEY + "="):
                old = str(line).split("=", 1)[1].strip()
    except Exception:
        pass
    value = f"{old}; {text}" if old else text
    try:
        return update_metadata(site, {"info": {_INFO_KEY: value}},
                               validate_coordinates=False)
    except Exception:
        return site


def _rows_text(f, rows) -> str:
    fr = sorted(float(f[r]) for r in rows)
    shown = ", ".join(f"{x:.4g}" for x in fr[:6])
    return f"{shown}{' …' if len(fr) > 6 else ''} Hz"


def mask(site, rows, comps=("xy", "yx")):
    """Mask *rows* of *comps* (Z and error -> NaN)."""
    new = copy.deepcopy(site)
    f, z, ze = _arrays(new)
    rows = np.asarray(sorted(set(int(r) for r in rows)), dtype=int)
    for comp in comps:
        i, j = COMPONENTS[comp]
        z[rows, i, j] = np.nan
        if ze is not None:
            ze[rows, i, j] = np.nan
    _write(new, z, ze)
    return _note(new, f"masked {'/'.join(c.upper() for c in comps)} at "
                 f"{_rows_text(f, rows)}")


def delete(site, rows):
    """Remove whole frequency rows."""
    from pycsamt.site.edit import select_freq

    f, _z, _ze = _arrays(site)
    drop = set(int(r) for r in rows)
    keep = np.array([i not in drop for i in range(f.size)])
    if keep.all():
        return site
    if not keep.any():
        raise ValueError("cannot delete every frequency of a station")
    new = copy.deepcopy(site)
    edi = select_freq(new.edi, keep=keep, inplace=True)
    try:
        new.edi = edi
    except Exception:
        pass
    return _note(new, f"deleted {_rows_text(f, drop)}")


def interpolate(site, rows, comps=("xy", "yx")):
    """Replace *rows* of *comps* with values interpolated from the valid
    neighbours (log f; log |Z| and unwrapped phase; log error)."""
    new = copy.deepcopy(site)
    f, z, ze = _arrays(new)
    rows = sorted(set(int(r) for r in rows))
    lf = np.log10(f)
    done = []
    for comp in comps:
        i, j = COMPONENTS[comp]
        zc = z[:, i, j]
        good = np.isfinite(zc)
        good[rows] = False
        if good.sum() < 2:
            raise ValueError(f"{comp.upper()}: fewer than two valid points "
                             "to interpolate from")
        order = np.argsort(lf[good])
        xs = lf[good][order]
        amp = np.log10(np.abs(zc[good]))[order]
        ph = np.unwrap(np.angle(zc[good])[order])
        for r in rows:
            a = np.interp(lf[r], xs, amp)
            p = np.interp(lf[r], xs, ph)
            z[r, i, j] = 10 ** a * np.exp(1j * p)
            if ze is not None:
                eg = ze[:, i, j][good][order]
                ok = np.isfinite(eg) & (eg > 0)
                if ok.sum() >= 2:
                    ze[r, i, j] = 10 ** np.interp(lf[r], xs[ok],
                                                  np.log10(eg[ok]))
        done.append(comp.upper())
    _write(new, z, ze)
    return _note(new, f"interpolated {'/'.join(done)} at "
                 f"{_rows_text(f, rows)}")


def nearest_rows(freqs, targets, rtol: float = 0.01) -> list[int]:
    """Rows of *freqs* within *rtol* of each target frequency."""
    freqs = np.asarray(freqs, dtype=float)
    out = []
    for t in np.atleast_1d(targets):
        if not freqs.size:
            break
        k = int(np.argmin(np.abs(freqs - t)))
        if abs(freqs[k] - t) <= rtol * abs(t):
            out.append(k)
    return out


def restore(site, original, rows, comps=("xy", "yx")):
    """Put back the *original* values of *rows* (matched by frequency)."""
    new = copy.deepcopy(site)
    f, z, ze = _arrays(new)
    fo, zo, zeo = _arrays(original)
    back = []
    for r in rows:
        k = nearest_rows(fo, [f[int(r)]], rtol=1e-6)
        if not k:
            continue
        for comp in comps:
            i, j = COMPONENTS[comp]
            z[int(r), i, j] = zo[k[0], i, j]
            if ze is not None and zeo is not None:
                ze[int(r), i, j] = zeo[k[0], i, j]
        back.append(int(r))
    if not back:
        return site
    _write(new, z, ze)
    return _note(new, f"restored {'/'.join(c.upper() for c in comps)} at "
                 f"{_rows_text(f, back)}")


def static_shift(site, comp: str, factor: float):
    """Scale the apparent resistivity of *comp* by *factor* (phase kept)."""
    factor = float(factor)
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError("the static-shift factor must be positive")
    new = copy.deepcopy(site)
    _f, z, ze = _arrays(new)
    i, j = COMPONENTS[comp]
    s = np.sqrt(factor)
    z[:, i, j] = z[:, i, j] * s
    if ze is not None:
        ze[:, i, j] = ze[:, i, j] * s
    _write(new, z, ze)
    return _note(new, f"static shift {comp.upper()} rho x{factor:.4g}")


def info_log(site) -> str:
    """The station's accumulated point-edit log ("" when none)."""
    try:
        info = site.edi.get_section("info")
        for line in getattr(info, "info_text", []) or []:
            if str(line).strip().upper().startswith(_INFO_KEY + "="):
                return str(line).split("=", 1)[1].strip()
    except Exception:
        pass
    return ""
