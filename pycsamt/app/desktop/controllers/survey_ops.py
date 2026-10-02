# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Whole-survey frequency and tensor operations for Edit ▸ Frequencies /
Tensor (Qt-free).

Each :class:`SurveyOp` declares its settings as
:class:`~pycsamt.app.desktop.controllers.inversion_engines.Field` (the
dialog builds its form from them), runs on a copy of the survey through
the library function that already implements it, and
:func:`preview_table` compares every station before and after (number of
frequencies, band, masked points, rotation) so the user sees the effect
before applying.

Library functions used: :func:`pycsamt.emtools.frequency.select_band`,
:func:`~pycsamt.emtools.frequency.regrid_logspace`,
:func:`~pycsamt.emtools.frequency.decimate_step`,
:func:`~pycsamt.emtools.frequency.align_grid`,
:func:`~pycsamt.emtools.frequency.drop_duplicates`,
:func:`pycsamt.site.edit.rotate_all`,
:func:`pycsamt.emtools.tensor.rotate_to_strike`; gap filling reuses the
Point Editor's interpolation (:func:`~pycsamt.app.desktop.controllers.
point_edits.interpolate`), only *between* measured frequencies.
"""

from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd

from pycsamt.app.desktop.controllers.inversion_engines import Field

__all__ = ["OPS", "SurveyOp", "fill_gaps", "op", "preview_table", "run"]


@dataclass(frozen=True)
class SurveyOp:
    key: str
    label: str
    menu: str  # "frequencies" | "tensor"
    help: str
    fields: tuple = field(default_factory=tuple)


def _band(values: dict):
    lo = float(values.get("fmin") or 0.0)
    hi = float(values.get("fmax") or 0.0)
    return (lo if lo > 0 else None), (hi if hi > 0 else None)


def _trim(sites, v):
    from pycsamt.emtools.frequency import select_band

    lo, hi = _band(v)
    if lo is None and hi is None:
        raise ValueError("give a lowest and/or a highest frequency")
    if lo and hi and lo >= hi:
        raise ValueError("the lowest frequency must be below the highest")
    return select_band(sites, fmin=lo, fmax=hi, inplace=False)


def _regrid(sites, v):
    from pycsamt.emtools.frequency import regrid_logspace

    lo, hi = _band(v)
    return regrid_logspace(sites, fmin=lo, fmax=hi,
                           per_decade=int(v["per_decade"]),
                           method=v["method"], inplace=False)


def _decimate(sites, v):
    from pycsamt.emtools.frequency import decimate_step

    return decimate_step(sites, step=int(v["step"]), inplace=False)


def _align(sites, v):
    from pycsamt.emtools.frequency import align_grid

    mode = v["mode"]
    ref = v.get("ref") or None
    if mode == "ref" and not ref:
        raise ValueError("choose the reference station")
    return align_grid(sites, mode=mode, ref_station=ref,
                      method=v["method"], inplace=False)


def _dedupe(sites, v):
    from pycsamt.emtools.frequency import drop_duplicates

    return drop_duplicates(sites, tol=float(v["tol"]), inplace=False)


def fill_gaps(sites, comps=("xy", "yx")):
    """Interpolate masked (NaN) values that lie *between* measured
    frequencies of each station; band edges are never extrapolated."""
    from pycsamt.app.desktop.controllers import point_edits as pe
    from pycsamt.site.base import Sites

    out = []
    for site in sites:
        cur = site
        cv = pe.curves(cur)
        f = cv["freq"]
        for comp in comps:
            valid = cv[comp]["valid"]
            if valid.sum() < 2:
                continue
            fv = f[valid]
            inside = (f > fv.min()) & (f < fv.max())
            rows = list(np.nonzero(~valid & inside)[0])
            if rows:
                cur = pe.interpolate(cur, rows, (comp,))
                cv = pe.curves(cur)
        out.append(cur)
    return Sites(out)


def _fill(sites, v):
    comps = {"both": ("xy", "yx"), "xy": ("xy",),
             "yx": ("yx",)}[v["comps"]]
    return fill_gaps(sites, comps)


def _rotate(sites, v):
    from pycsamt.site.edit import rotate_all

    angle = float(v["angle"])
    if abs(angle) < 1e-12:
        raise ValueError("a rotation of 0° changes nothing")
    return rotate_all(sites, angle, inplace=False)


def _strike(sites, v):
    from pycsamt.emtools.tensor import rotate_to_strike

    return rotate_to_strike(copy.deepcopy(sites), method=v["method"],
                            inplace=False)


_METHOD = Field("method", "Interpolation", "choice", "nearest",
                choices=(("nearest", "Nearest measured frequency"),
                         ("linear", "Linear (in log f)")),
                help="Values at the new frequencies. Frequencies outside "
                     "a station's own band stay empty (masked).")
_FMIN = Field("fmin", "Lowest", "float", 0.0, 0.0, 1e6, 0.1, 5, unit="Hz",
              auto_zero=True, help="Auto = no lower limit")
_FMAX = Field("fmax", "Highest", "float", 0.0, 0.0, 1e6, 10.0, 5, unit="Hz",
              auto_zero=True, help="Auto = no upper limit")

OPS: tuple[SurveyOp, ...] = (
    SurveyOp("trim", "Trim Band", "frequencies",
             "Keep only the frequencies inside a band.", (_FMIN, _FMAX)),
    SurveyOp("regrid", "Regrid (log-spaced)", "frequencies",
             "Resample every station onto one log-spaced frequency grid.",
             (Field("per_decade", "Per decade", "int", 6, 1, 40, 1,
                    help="Frequencies per decade"), _FMIN, _FMAX, _METHOD)),
    SurveyOp("decimate", "Decimate", "frequencies",
             "Keep every n-th frequency.",
             (Field("step", "Keep every", "int", 2, 2, 20, 1,
                    unit="th"),)),
    SurveyOp("align", "Align to a Common Grid", "frequencies",
             "Put every station on the same frequencies: all of them "
             "(union), only the shared ones (intersection) or one "
             "station's.",
             (Field("mode", "Grid", "choice", "intersection",
                    choices=(("intersection", "Shared frequencies"),
                             ("union", "All frequencies"),
                             ("ref", "A reference station's"))),
              Field("ref", "Reference", "choice", "", choices=(("", "—"),)),
              _METHOD)),
    SurveyOp("dedupe", "Remove Duplicate Frequencies", "frequencies",
             "Drop repeated frequency rows.",
             (Field("tol", "Tolerance", "float", 1e-10, 0.0, 1.0, 1e-10, 12,
                    help="Relative difference below which two "
                         "frequencies are the same"),)),
    SurveyOp("fill", "Fill Gaps", "frequencies",
             "Interpolate masked values between measured frequencies "
             "(never beyond a station's band).",
             (Field("comps", "Components", "choice", "both",
                    choices=(("both", "XY and YX"), ("xy", "XY"),
                             ("yx", "YX"))),)),
    SurveyOp("rotate", "Rotate", "tensor",
             "Rotate every station's tensor by an angle (clockwise from "
             "north).",
             (Field("angle", "Angle", "float", 0.0, -180.0, 180.0, 5.0, 2,
                    unit="°"),)),
    SurveyOp("strike", "Rotate to Strike", "tensor",
             "Rotate each station to its own geoelectric strike.",
             (Field("method", "Method", "choice", "swift",
                    choices=(("swift", "Swift"),
                             ("phase_diff", "Phase difference"))),)),
)

_RUN: dict[str, Callable[[Any, dict], Any]] = {
    "trim": _trim, "regrid": _regrid, "decimate": _decimate,
    "align": _align, "dedupe": _dedupe, "fill": _fill, "rotate": _rotate,
    "strike": _strike,
}


def op(key: str) -> SurveyOp:
    return next(o for o in OPS if o.key == key)


def run(key: str, sites, values: dict):
    """Run operation *key* on a copy of *sites* with *values*."""
    from pycsamt.site.base import to_sites

    out = _RUN[key](copy.deepcopy(sites), dict(values))
    return to_sites(out)


def _stats(site) -> dict:
    from pycsamt.app.desktop.controllers import point_edits as pe

    try:
        cv = pe.curves(site)
    except Exception:
        return {"n": 0, "fmin": np.nan, "fmax": np.nan, "masked": 0}
    f = cv["freq"]
    ok = np.isfinite(f)
    masked = int((~cv["xy"]["valid"]).sum() + (~cv["yx"]["valid"]).sum())
    return {"n": int(ok.sum()), "fmin": float(np.nanmin(f)) if ok.any()
            else np.nan, "fmax": float(np.nanmax(f)) if ok.any() else np.nan,
            "masked": masked}


def _strike_angle(site) -> float:
    try:
        from pycsamt.site.compute import strike_estimate

        return float(strike_estimate(site.edi))
    except Exception:
        return float("nan")


def preview_table(before, after, *, with_strike: bool = False) -> pd.DataFrame:
    """One row per station: frequencies, band and masked points before ->
    after (and the strike angle when rotating to strike)."""
    rows = []
    after_by = {str(s.name): s for s in after}
    for s in before:
        name = str(s.name)
        a = after_by.get(name)
        b0, a0 = _stats(s), _stats(a) if a is not None else None
        row = {"station": name, "n_before": b0["n"],
               "n_after": a0["n"] if a0 else 0,
               "band_before": _fmt_band(b0), "band_after": _fmt_band(a0),
               "masked_before": b0["masked"],
               "masked_after": a0["masked"] if a0 else 0}
        if with_strike:
            row["strike"] = _strike_angle(s)
        rows.append(row)
    return pd.DataFrame(rows)


def _fmt_band(st) -> str:
    if not st or not np.isfinite(st["fmin"]):
        return "—"
    return f"{st['fmin']:.4g} – {st['fmax']:.4g} Hz"
