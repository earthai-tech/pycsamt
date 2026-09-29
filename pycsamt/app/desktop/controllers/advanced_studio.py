# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Advanced Tools studio: sections, per-plot options and inputs (Qt-free).

The analysis sections are :data:`~pycsamt.app.desktop.controllers.
advanced_controller.ADVANCED_GROUPS`.  For each emtools plot,
:func:`plot_options` reads the function signature and offers only the
options that function really takes (colour map, period window, period,
colour bar, strike method, ...).  :func:`plot_kwargs` turns the form values
back into keyword arguments, leaving the function's own default wherever
the user kept "Default"/"Auto".

Some plots need an input the survey has to supply:

* ``station`` -- one station (phase-tensor strip);
* ``lines`` -- stations grouped into profile lines (strip grid, strike rose
  by line), built by :func:`line_groups`;
* ``model`` -- a trained dictionary (ATOM pseudosection).
"""

from __future__ import annotations

import inspect
import re
from dataclasses import dataclass

from pycsamt.app.desktop.controllers.advanced_controller import (
    ADVANCED_GROUP_ICONS,
    ADVANCED_GROUPS,
    _POLAR_FNS,
)
from pycsamt.app.desktop.controllers.inversion_engines import Field

__all__ = [
    "SECTIONS",
    "Section",
    "line_groups",
    "plot_inputs",
    "plot_kwargs",
    "plot_options",
    "site_names",
]

_UTILITIES = {"Topography": "topo", "Conversion": "conv"}


@dataclass(frozen=True)
class Section:
    key: str  # "plots" | "topo" | "conv"
    label: str
    icon: str
    plots: tuple  # (label, fn_name, has_ax) for plot sections
    help: str = ""


_HELP = {
    "Strike Analysis": "Geoelectric strike: roses, stability, map sticks.",
    "Phase Tensor": "Ellipses, roses, skew and dimensionality.",
    "Induction / Tipper": "Induction arrows, tipper maps and hodograms.",
    "Impedance / Z": "Mohr circles, Argand, invariants, anisotropy.",
    "Depth Imaging": "Bostick, gradient and composite depth sections.",
    "Survey Tools": "Survey-wide fingerprints and comparisons.",
    "Topography": "Station elevations: source, style, preview.",
    "Conversion": "AVG / J / spectra files → EDI.",
}

SECTIONS: tuple[Section, ...] = tuple(
    Section(_UTILITIES.get(label, "plots"), label,
            ADVANCED_GROUP_ICONS.get(label, "advanced-tools"),
            tuple(plots), _HELP.get(label, ""))
    for label, plots in ADVANCED_GROUPS
)

# common colour maps offered for every function taking ``cmap``
_CMAPS = (("", "Default"), ("viridis", "viridis"), ("plasma", "plasma"),
          ("magma", "magma"), ("cividis", "cividis"),
          ("RdBu_r", "RdBu (diverging)"), ("coolwarm", "coolwarm"),
          ("Spectral_r", "Spectral"), ("jet_r", "jet (legacy)"),
          ("turbo", "turbo"))

_LINE_MODES = (("prefix", "By name prefix (e.g. L1-, L2-)"),
               ("single", "All stations on one line"))


def _fn(fn_name: str):
    import pycsamt.emtools as et

    return getattr(et, fn_name, None)


def _params(fn_name: str) -> dict:
    fn = _fn(fn_name)
    if fn is None:
        return {}
    try:
        return dict(inspect.signature(fn).parameters)
    except (TypeError, ValueError):
        return {}


def _plain(default):
    """A usable default (library sentinels / None -> None)."""
    if default is inspect.Parameter.empty or default is None:
        return None
    if isinstance(default, (bool, int, float, str, tuple)):
        return default
    return None


def plot_options(fn_name: str) -> list[Field]:
    """The option fields *fn_name* accepts, with its own defaults."""
    p = _params(fn_name)
    out: list[Field] = []
    if "cmap" in p:
        out.append(Field("cmap", "Colour map", "choice", "", choices=_CMAPS,
                         help="Default keeps the plot's own colour map."))
    if "period" in p:
        d = _plain(p["period"].default)
        out.append(Field("period", "Period", "float",
                         float(d) if isinstance(d, (int, float)) else 1.0,
                         1e-5, 1e5, 1.0, 4, unit="s",
                         help="Period shown (snapped to the nearest one "
                              "the data has)."))
    if "period_range" in p:
        out.append(Field("pmin", "Period from", "float", 0.0, 0.0, 1e5, 0.1,
                         4, unit="s", auto_zero=True,
                         help="Auto = the whole band"))
        out.append(Field("pmax", "Period to", "float", 0.0, 0.0, 1e5, 10.0,
                         4, unit="s", auto_zero=True,
                         help="Auto = the whole band"))
    if "method" in p and _plain(p["method"].default) == "consensus":
        out.append(Field("method", "Strike method", "choice", "consensus",
                         choices=(("consensus", "Consensus"),
                                  ("sweep", "Rotation sweep"),
                                  ("pt", "Phase tensor"))))
    if "bins" in p and isinstance(_plain(p["bins"].default), int):
        out.append(Field("bins", "Rose bins", "int", p["bins"].default, 8,
                         180, 4))
    for key, label in (("show_colorbar", "Colour bar"),
                       ("station_labels", "Station labels"),
                       ("normalize", "Normalise")):
        if key in p and isinstance(_plain(p[key].default), bool):
            out.append(Field(key, label, "bool", p[key].default))
    if "scale" in p:
        out.append(Field("scale", "Symbol scale", "float", 0.0, 0.0, 100.0,
                         0.1, 2, auto_zero=True, advanced=True,
                         help="Auto = the plot's own sizing"))
    return out


def plot_inputs(fn_name: str) -> tuple[str, ...]:
    """Inputs the survey must supply: "station", "lines", "model"."""
    p = _params(fn_name)
    need = []
    if fn_name == "plot_phase_tensor_strip" or (
            "station" in p and fn_name.endswith("_strip")):
        need.append("station")
    if "profiles" in p or "groups" in p:
        need.append("lines")
    if fn_name == "plot_atom_psection":
        need.append("model")
    return tuple(need)


def plot_kwargs(fn_name: str, values: dict, *, station: str = "",
                lines: dict | None = None) -> dict:
    """Keyword arguments for *fn_name* from the options form."""
    p = _params(fn_name)
    kw: dict = {}
    if values.get("cmap"):
        kw["cmap"] = values["cmap"]
    if "period" in values and "period" in p:
        kw["period"] = float(values["period"])
    lo, hi = float(values.get("pmin") or 0), float(values.get("pmax") or 0)
    if "period_range" in p and lo > 0 and hi > 0:
        kw["period_range"] = (min(lo, hi), max(lo, hi))
    for key in ("method", "bins", "show_colorbar", "station_labels",
                "normalize"):
        if key in values and key in p:
            kw[key] = values[key]
    if float(values.get("scale") or 0) > 0 and "scale" in p:
        kw["scale"] = float(values["scale"])
    if station and "station" in p:
        kw["station"] = station
    if lines:
        if "profiles" in p:
            kw["profiles"] = lines
        elif "groups" in p:
            kw["groups"] = lines
    return kw


def is_polar(fn_name: str) -> bool:
    return fn_name in _POLAR_FNS


def site_names(sites) -> list[str]:
    if sites is None:
        return []
    try:
        from pycsamt.emtools._core import _iter_items, _name

        return [str(_name(ed, i)) for i, ed in enumerate(_iter_items(sites))]
    except Exception:
        return []


def line_groups(names: list[str], mode: str = "prefix") -> dict:
    """Group station names into profile lines.

    ``"prefix"`` groups by the name with its trailing station number
    removed (``L1-03`` / ``L1-04`` -> ``L1``; ``gv100`` -> ``gv``);
    ``"single"`` puts every station on one line.
    """
    names = list(names)
    if not names:
        return {}
    if mode == "single":
        return {"all stations": names}
    groups: dict[str, list[str]] = {}
    for n in names:
        key = re.sub(r"[-_ .]*\d+[A-Za-z]?$", "", n) or "line"
        groups.setdefault(key, []).append(n)
    return groups


LINE_MODES = _LINE_MODES
