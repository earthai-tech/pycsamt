# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Airborne EM studio catalogue (Qt-free).

Views follow the *technology* of the loaded data -- ``ztem``,
``afmag_original`` (scalar tilt comparator), ``afmag_airmt`` (tensor
AirMt) or ``mobilemt`` -- because each ``pycsamt.emtools`` function reads
a different quantity (tipper, tilt angles, interstation tensor,
admittance) and draws an empty figure for the others.  The mapping below
was established by running every function on every technology.

Each view's options are built from its function signature
(:func:`options_for`), with the choices taken from the data: the
frequencies it was measured at, its stations, its flight lines.  The
geomagnetic / aircraft-attitude settings of the motion views are shared
(:data:`GEOMETRY_FIELDS`) so they are typed once.

:func:`run_view` returns a Figure for plots and a list of
``(title, DataFrame)`` for tables (diagnostics dictionaries are split into
a summary table plus their profile tables).
"""

from __future__ import annotations

from pycsamt.api.colormaps import colormap_choices

import importlib
import inspect
from dataclasses import dataclass, replace

import numpy as np
import pandas as pd

from pycsamt.app.desktop.controllers.inversion_engines import Field

__all__ = [
    "GEOMETRY_FIELDS",
    "TECH_LABELS",
    "VIEWS",
    "AirView",
    "data_context",
    "options_for",
    "run_view",
    "views_for",
]

TECH_LABELS = {
    "ztem": "ZTEM (tipper)",
    "afmag_original": "AFMAG (original tilt comparator)",
    "afmag_airmt": "AirMt / tensor AFMAG",
    "mobilemt": "MobileMT (admittance)",
}

_Z, _A, _M = "pycsamt.emtools.ztem", "pycsamt.emtools.afmag", \
    "pycsamt.emtools.mobilemt"
ALL = frozenset(TECH_LABELS)


@dataclass(frozen=True)
class AirView:
    key: str
    label: str
    group: str  # Profiles | Sections | Maps | Motion | Tables
    module: str
    fn: str
    techs: frozenset
    help: str = ""
    kind: str = "plot"  # plot | table
    polar: bool = False


V = AirView
VIEWS: tuple[AirView, ...] = (
    # ── ZTEM ─────────────────────────────────────────────────────────
    V("ztem_tipper", "Tipper profile", "Profiles", _Z,
      "plot_ztem_tipper_profile", frozenset({"ztem"}),
      "In-phase / quadrature tipper along the line at one frequency."),
    V("ztem_div", "Divergence profile", "Profiles", _Z,
      "plot_ztem_divergence_profile", frozenset({"ztem"}),
      "Total-divergence (Peaker) profile at one frequency."),
    V("ztem_rot", "Phase-rotation profile", "Profiles", _Z,
      "plot_ztem_phase_rotation_profile", frozenset({"ztem"}),
      "Raw vs Hilbert phase-rotated response: crossover → peak."),
    V("ztem_div_ps", "Divergence pseudosection", "Sections", _Z,
      "plot_ztem_divergence_psection", frozenset({"ztem"}),
      "Station × frequency total divergence."),
    V("ztem_div_grid", "Divergence sections, all lines", "Sections", _Z,
      "plot_ztem_divergence_psection_grid", frozenset({"ztem"}),
      "Every flight line's divergence section on one colour scale."),
    V("ztem_band", "Usable-band mask", "Sections", _Z,
      "plot_ztem_band_mask_psection", frozenset({"ztem"}),
      "|T| before / after masking outside the ZTEM band."),
    V("ztem_map", "Tipper / divergence map", "Maps", _Z, "plot_ztem_map",
      frozenset({"ztem"}), "Plan-view map at one frequency."),
    V("tilt_profile", "Tilt profile (from tipper)", "Profiles", _A,
      "plot_afmag_tilt_profile", frozenset({"ztem"}),
      "AFMAG-style tilt angles computed from the tipper."),
    V("tilt_ps", "Tilt pseudosection (from tipper)", "Sections", _A,
      "plot_afmag_tilt_psection", frozenset({"ztem"}),
      "Tilt angle, station × period."),
    V("tilt_polar", "Tilt polar (one station)", "Profiles", _A,
      "plot_afmag_tilt_polar", frozenset({"ztem"}),
      "Tilt azimuth and magnitude against period.", polar=True),
    # ── AFMAG original / AirMt ───────────────────────────────────────
    V("orig_tilt", "Tilt profile", "Profiles", _A,
      "plot_original_afmag_tilt_profile", frozenset({"afmag_original"}),
      "Classic AFMAG tilt-angle profile at one frequency."),
    V("orig_dual", "Dual-frequency profile", "Profiles", _A,
      "plot_original_afmag_dual_frequency_profile",
      frozenset({"afmag_original"}),
      "Low vs high frequency tilt: crossover shift → conductor depth."),
    V("airmt_tilt", "Tilt profile", "Profiles", _A,
      "plot_airmt_tilt_profile", frozenset({"afmag_airmt"}),
      "AirMt tilt from the interstation tensor, one frequency."),
    V("airmt_ps", "Tilt pseudosection", "Sections", _A,
      "plot_airmt_tilt_psection", frozenset({"afmag_airmt"}),
      "AirMt tilt, station × period."),
    # ── MobileMT ─────────────────────────────────────────────────────
    V("mmt_adm", "Admittance profile", "Profiles", _M,
      "plot_mobilemt_admittance_profile", frozenset({"mobilemt"}),
      "One admittance component along the line."),
    V("mmt_skew", "Skew profile", "Profiles", _M,
      "plot_mobilemt_skew_profile", frozenset({"mobilemt"}),
      "Admittance skew (3-D indicator) along the line."),
    V("mmt_cond", "Conductivity pseudosection", "Sections", _M,
      "plot_mobilemt_conductivity_psection", frozenset({"mobilemt"}),
      "Apparent conductivity, station × frequency."),
    # ── every technology ─────────────────────────────────────────────
    V("lines", "Flight lines", "Maps", _Z, "plot_ztem_flight_lines", ALL,
      "Plan view of the stations, one colour per flight line."),
    V("motion_map", "Motion susceptibility map", "Motion", _A,
      "plot_motion_susceptibility_map", ALL,
      "How much aircraft attitude couples into each station (needs the "
      "survey geometry below)."),
    V("motion_cmp", "Motion correction: before / after", "Motion", _A,
      "plot_afmag_correction_comparison", frozenset({"ztem"}),
      "Tilt before and after masking motion-susceptible frequencies."),
    # ── tables ───────────────────────────────────────────────────────
    V("t_div", "Total divergence", "Tables", _Z, "total_divergence_table",
      frozenset({"ztem"}), kind="table"),
    V("t_rot", "Phase rotation", "Tables", _Z, "phase_rotate_table",
      frozenset({"ztem"}), kind="table"),
    V("t_cross", "Crossover diagnostics", "Tables", _Z,
      "ztem_crossover_diagnostics", frozenset({"ztem"}), kind="table"),
    V("t_orig", "Tilt table", "Tables", _A, "original_afmag_tilt_table",
      frozenset({"afmag_original"}), kind="table"),
    V("t_cond", "Conductor diagnostics", "Tables", _A,
      "original_afmag_conductor_diagnostics",
      frozenset({"afmag_original"}), kind="table"),
    V("t_adm", "Admittance", "Tables", _M, "admittance_table",
      frozenset({"mobilemt"}), kind="table"),
    V("t_det", "Admittance determinant", "Tables", _M,
      "admittance_determinant_table", frozenset({"mobilemt"}),
      kind="table"),
    V("t_skew", "Admittance skew", "Tables", _M, "admittance_skew_table",
      frozenset({"mobilemt"}), kind="table"),
    V("t_motion", "Motion susceptibility", "Tables", _A,
      "motion_susceptibility_table", ALL, kind="table"),
)

GROUPS = ("Profiles", "Sections", "Maps", "Motion", "Tables")

GEOMETRY_FIELDS: tuple[Field, ...] = (
    Field("inclination", "Inclination", "float", 60.0, -90.0, 90.0, 1.0, 1,
          unit="°", help="Geomagnetic field inclination at the survey"),
    Field("declination", "Declination", "float", 0.0, -180.0, 180.0, 1.0, 1,
          unit="°"),
    Field("roll_amplitude_deg", "Roll", "float", 5.0, 0.0, 45.0, 0.5, 1,
          unit="°", help="Typical aircraft roll amplitude"),
    Field("pitch_amplitude_deg", "Pitch", "float", 5.0, 0.0, 45.0, 0.5, 1,
          unit="°"),
    Field("yaw_amplitude_deg", "Yaw", "float", 0.0, 0.0, 45.0, 0.5, 1,
          unit="°"),
)
_GEOMETRY = {f.key for f in GEOMETRY_FIELDS}


def views_for(techs) -> list[AirView]:
    """Views that can draw the loaded technologies."""
    techs = set(techs or ())
    return [v for v in VIEWS if v.techs & techs]


# ── data context ──────────────────────────────────────────────────────────
def data_context(asites) -> dict:
    """Frequencies, stations and lines of the loaded data."""
    freqs: set[float] = set()
    for s in asites:
        f = getattr(s, "freq", None)
        if f is not None:
            freqs.update(float(x) for x in np.asarray(f, float).ravel()
                         if np.isfinite(x))
    return {
        "freqs": sorted(freqs, reverse=True),
        "stations": [str(s.name) for s in asites],
        "lines": list(getattr(asites, "line_ids", ()) or ()),
        "techs": list(getattr(asites, "technologies", ()) or ()),
    }


def _freq_choices(ctx, auto_label="Auto (the function's choice)"):
    return (("", auto_label),) + tuple(
        (f"{f:.6g}", f"{f:.4g} Hz") for f in ctx.get("freqs", []))


_CMAPS = tuple(colormap_choices(first=("RdBu_r",)))
_COMPONENTS = {
    ("tzx", "tzy"): (("tzx", "Tzx"), ("tzy", "Tzy")),
    ("abs",): (("abs", "|T|"), ("tzx", "Tzx"), ("tzy", "Tzy")),
    ("real", "imag", "resultant"): (("resultant", "Resultant"),
                                    ("real", "In-phase"),
                                    ("imag", "Quadrature")),
    ("det",): tuple((c, c.upper() if len(c) > 2 else c) for c in
                    ("det", "xx", "xy", "yx", "yy", "hzx", "hzy")),
}


def _component_field(default) -> Field:
    for keys, choices in _COMPONENTS.items():
        if default in keys:
            return Field("component", "Component", "choice", default,
                         choices=choices)
    return Field("component", "Component", "choice", default,
                 choices=((default, str(default)),))


def options_for(view: AirView, ctx: dict) -> list[Field]:
    """Option fields for *view*'s function, choices from *ctx* (data)."""
    fn = getattr(importlib.import_module(view.module), view.fn)
    params = inspect.signature(fn).parameters
    out: list[Field] = []

    def d(name):
        v = params[name].default
        return None if v is inspect.Parameter.empty else v

    for name in params:
        if name in ("sites", "dataset", "before_sites", "ax", "axes",
                    "figsize", "panel_size", "station_style", "clim",
                    "after_sites", "system_spec", "flag_kws", "period_s",
                    "station_preset", "kwargs"):
            continue
        if name in _GEOMETRY:
            continue  # shared survey-geometry group
        if name == "frequency_hz":
            out.append(Field("frequency_hz", "Frequency", "choice", "",
                             choices=_freq_choices(ctx)))
        elif name in ("freq_low_hz", "freq_high_hz"):
            label = "Low frequency" if "low" in name else "High frequency"
            out.append(Field(name, label, "choice", "",
                             choices=_freq_choices(ctx, "Auto")))
        elif name == "component":
            out.append(_component_field(d(name)))
        elif name == "part":
            opts = (("real", "In-phase"), ("imag", "Quadrature"))
            if d(name) == "abs":
                opts = (("abs", "|value|"),) + opts
            out.append(Field("part", "Part", "choice", d(name), choices=opts))
        elif name == "quantity":
            out.append(Field("quantity", "Quantity", "choice", d(name),
                             choices=(("tipper", "Tipper"),
                                      ("divergence", "Total divergence"))))
        elif name == "source":
            out.append(Field("source", "Conductivity", "choice", d(name),
                             choices=(("theoretical", "Theoretical"),
                                      ("native", "As delivered"))))
        elif name == "station":
            out.append(Field("station", "Station", "choice", "",
                             choices=(("", "First station"),) + tuple(
                                 (s, s) for s in ctx.get("stations", []))))
        elif name == "line_id":
            out.append(Field("line_id", "Flight line", "choice", "",
                             choices=(("", "First line"),) + tuple(
                                 (s, s) for s in ctx.get("lines", []))))
        elif name in ("cmap", "delta_cmap"):
            dv = d(name)
            opts = (dv,) + tuple(c for c in _CMAPS if c != dv)
            out.append(Field(name, "Colour map" if name == "cmap" else
                             "Change colour map", "choice", dv,
                             choices=tuple((c, c) for c in opts)))
        elif name == "clim_pct" and isinstance(d(name), float):
            out.append(Field("clim_pct", "Colour clip", "float", d(name),
                             50.0, 100.0, 1.0, 1, unit="%",
                             help="Percentile of |value| setting the "
                                  "colour limits"))
        elif name in ("show_grid", "show_contour", "show_stations",
                      "as_percent"):
            label = {"show_grid": "Grid", "show_contour": "Contours",
                     "show_stations": "Stations",
                     "as_percent": "As percent"}[name]
            out.append(Field(name, label, "bool", bool(d(name))))
        elif name in ("n_contour_levels", "n_grid", "max_lines", "n_cols",
                      "station_label_step", "n_resample"):
            spec = {"n_contour_levels": ("Contour levels", 1, 20),
                    "n_grid": ("Grid size", 20, 400),
                    "max_lines": ("Max lines", 1, 30),
                    "n_cols": ("Columns", 1, 6),
                    "station_label_step": ("Label every", 1, 50),
                    "n_resample": ("Resample points", 0, 2000)}[name]
            dv = d(name)
            out.append(Field(name, spec[0], "int", int(dv or 0), spec[1],
                             spec[2], 1, auto_zero=dv is None))
        elif name == "spacing_m":
            out.append(Field("spacing_m", "Station spacing", "float",
                             d(name), 1.0, 10000.0, 10.0, 1, unit="m"))
        elif name == "band_hz":
            out.append(Field("band_lo", "Band from", "float", 0.0, 0.0, 1e6,
                             1.0, 1, unit="Hz", auto_zero=True,
                             help="Auto = the system's usable band"))
            out.append(Field("band_hi", "Band to", "float", 0.0, 0.0, 1e6,
                             10.0, 1, unit="Hz", auto_zero=True))
    return out


def needs_geometry(view: AirView) -> bool:
    fn = getattr(importlib.import_module(view.module), view.fn)
    return bool(_GEOMETRY & set(inspect.signature(fn).parameters))


# ── running ───────────────────────────────────────────────────────────────
def _kwargs(view: AirView, values: dict, geometry: dict) -> dict:
    fn = getattr(importlib.import_module(view.module), view.fn)
    params = inspect.signature(fn).parameters
    kw: dict = {}
    for k, v in values.items():
        if k in ("band_lo", "band_hi"):
            continue
        if k not in params:
            continue
        if k in ("frequency_hz", "freq_low_hz", "freq_high_hz",
                 "station", "line_id"):
            if v in ("", None):
                continue
            kw[k] = float(v) if k.endswith("_hz") else v
        elif k in ("n_resample",) and not v:
            continue
        else:
            kw[k] = v
    lo, hi = float(values.get("band_lo") or 0), float(values.get("band_hi")
                                                       or 0)
    if "band_hz" in params and lo > 0 and hi > lo:
        kw["band_hz"] = (lo, hi)
    for k in _GEOMETRY & set(params):
        kw[k] = float(geometry.get(k, 0.0))
    return kw


def run_view(view: AirView, asites, values: dict | None = None,
             geometry: dict | None = None):
    """A Figure (plot views) or ``[(title, DataFrame), ...]`` (tables)."""
    import matplotlib.pyplot as plt

    fn = getattr(importlib.import_module(view.module), view.fn)
    params = inspect.signature(fn).parameters
    kw = _kwargs(view, values or {}, geometry or {})
    if view.kind == "table":
        return _as_tables(view.label, fn(asites, **kw))
    if "ax" in params:
        subplot = {"projection": "polar"} if view.polar else {}
        fig, ax = plt.subplots(figsize=(9.5, 5.2) if not view.polar
                               else (6.5, 6.5), subplot_kw=subplot)
        out = fn(asites, ax=ax, **kw)
        return getattr(out, "figure", None) or fig
    out = fn(asites, **kw)
    return out if hasattr(out, "savefig") else plt.gcf()


def _as_tables(title: str, out) -> list[tuple[str, pd.DataFrame]]:
    if isinstance(out, pd.DataFrame):
        return [(title, out)]
    if isinstance(out, dict):
        scalars = {k: v for k, v in out.items()
                   if not isinstance(v, (pd.DataFrame, dict, list,
                                         np.ndarray))}
        tables = [(title, pd.DataFrame({"quantity": list(scalars),
                                        "value": list(scalars.values())}))]
        for k, v in out.items():
            if isinstance(v, pd.DataFrame):
                tables.append((f"{title}: {k}", v))
        return tables
    return [(title, pd.DataFrame({"value": [out]}))]


def with_choices(field_: Field, **changes) -> Field:
    return replace(field_, **changes)
