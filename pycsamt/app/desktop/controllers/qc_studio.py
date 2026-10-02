# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
QC studio (Qt-free): diagnostics catalogue, options, line scope, summary.

Diagnostics
    :data:`VIEWS` -- the transfer-function QC plots of
    :mod:`pycsamt.emtools`, grouped as Confidence, Coverage, Noise / SNR,
    Dimensionality & skew, Static shift, Distortion & source and Strike.
Options
    :func:`options_for` turns each function's parameters (from
    :func:`~pycsamt.app.desktop.controllers.qc_controller.qc_parameter_specs`)
    into form fields: choices as drop-downs, colour maps from the shared
    catalogue (``jet_r`` and ``jet`` included), stations from the survey,
    number pairs and lists as text.  :func:`to_kwargs` checks and converts
    them back.
Line scope
    QC profiles and pseudo-sections assume one profile.  With several lines
    :func:`render` draws one panel per line (or one line only) instead of
    projecting every station onto one bearing.
Summary
    :func:`scorecard` -- one row per station: composite confidence, data
    coverage, median SNR, median |skew|, flags and a pass / warn / fail
    status.
"""

from __future__ import annotations

import inspect
import math
from dataclasses import dataclass

import numpy as np

from pycsamt.api.colormaps import colormap_choices, colormap_label
from pycsamt.app.desktop.controllers.inversion_engines import Field
from pycsamt.app.desktop.controllers.qc_controller import (
    QC_STATIC_SHIFT_METHODS,
    QC_STATIC_SHIFT_PLOTS,
    QCController,
    describe_plot,
    qc_parameter_specs,
)

__all__ = [
    "GROUPS",
    "QCView",
    "VIEWS",
    "line_groups_of",
    "options_for",
    "render",
    "scorecard",
    "subset",
    "to_kwargs",
    "view",
]


@dataclass(frozen=True)
class QCView:
    fn: str
    label: str
    group: str
    per_line: bool = False  # a profile / pseudo-section along one line
    station: bool = False  # draws one station
    lines_arg: str = ""  # "groups" | "lines": takes the line grouping
    help: str = ""

    @property
    def key(self) -> str:
        return self.fn

    @property
    def description(self) -> str:
        return self.help or describe_plot(self.fn)


GROUPS = ("Confidence", "Coverage", "Noise / SNR", "Dimensionality & skew",
          "Static shift", "Distortion & source", "Strike")

_V = QCView
VIEWS: tuple[QCView, ...] = (
    # confidence
    _V("plot_confidence_profile", "Confidence profile", "Confidence",
       per_line=True, lines_arg="lines"),
    _V("plot_frequency_confidence_psection", "Confidence pseudo-section",
       "Confidence", per_line=True),
    _V("plot_confidence_heatmap", "Confidence heat map", "Confidence",
       per_line=True, help="Station × frequency confidence."),
    _V("plot_confidence_band_summary", "Confidence by band", "Confidence"),
    _V("plot_confidence_distribution", "Confidence distribution",
       "Confidence", help="Histogram and cumulative share of stations."),
    _V("plot_confidence_rank", "Station ranking", "Confidence",
       help="Stations ordered from least to most trusted."),
    _V("plot_confidence_map", "Confidence map", "Confidence",
       help="Station confidence in plan view."),
    _V("plot_confidence_grid_map", "Confidence grid map", "Confidence",
       help="Interpolated confidence surface over the survey area."),
    _V("plot_confidence_risk_map", "Confidence risk map", "Confidence",
       help="Where low confidence clusters in plan view."),
    _V("plot_confidence_component_map", "Confidence components",
       "Confidence", help="Each score component (coverage, uncertainty, "
                          "off-diagonal, diagonal, phase) mapped."),
    _V("plot_confidence_coverage_curve", "Confidence vs threshold",
       "Confidence", help="Share of stations kept for each confidence "
                          "threshold."),
    _V("plot_confidence_method_comparison", "Presence vs composite",
       "Confidence", help="The two confidence criteria side by side."),
    _V("plot_station_confidence_dashboard", "Station dashboard",
       "Confidence", station=True,
       help="One station: its score components and spectrum."),
    _V("plot_station_confidence_spectrum", "Station confidence spectrum",
       "Confidence", station=True),
    # coverage
    _V("plot_survey_inventory_overview", "Survey inventory", "Coverage",
       help="Frequencies and components recorded at each station."),
    _V("plot_coverage", "Coverage per station", "Coverage"),
    _V("plot_coverage_psection", "Coverage pseudo-section", "Coverage",
       per_line=True),
    _V("plot_coverage_quality_heatmap", "Coverage quality heat map",
       "Coverage", per_line=True),
    _V("plot_band_microstrips", "Band micro-strips", "Coverage",
       help="Per-band data presence strips along the line."),
    _V("plot_polar_coverage", "Polar coverage", "Coverage"),
    # noise
    _V("plot_snr_hist", "SNR histogram", "Noise / SNR"),
    _V("plot_snr_section", "SNR section", "Noise / SNR", per_line=True,
       help="Signal-to-noise ratio by station and period."),
    _V("plot_xyyx_crossover_map", "xy / yx crossover", "Noise / SNR",
       per_line=True, help="Where the two off-diagonal responses cross."),
    _V("plot_offdiag_antisym_residual", "Off-diagonal residual",
       "Noise / SNR", per_line=True,
       help="Departure of Zxy from −Zyx (1-D antisymmetry)."),
    _V("nr_qc_delta_offdiag_psection", "ΔZ off-diagonal (denoising)",
       "Noise / SNR", per_line=True),
    _V("nr_qc_snr_gain_profile", "SNR gain (denoising)", "Noise / SNR",
       per_line=True),
    _V("nr_qc_harmonic_waterfall", "Mains harmonics (denoising)",
       "Noise / SNR"),
    _V("nr_qc_station_offdiag_curves", "Station curves (denoising)",
       "Noise / SNR", station=True),
    # dimensionality
    _V("plot_skew_traffic_psection", "Skew traffic light",
       "Dimensionality & skew", per_line=True),
    _V("plot_skew_percentile_ribbon", "Skew percentile ribbon",
       "Dimensionality & skew"),
    _V("plot_skew_vote_band", "Skew vote by band", "Dimensionality & skew",
       help="Share of stations above the skew threshold in each band."),
    _V("plot_dimensionality_psection", "Dimensionality section",
       "Dimensionality & skew", per_line=True),
    _V("plot_dimensionality_grid", "Dimensionality grid",
       "Dimensionality & skew", per_line=True),
    _V("plot_dim_confidence_grid", "Dimensionality confidence",
       "Dimensionality & skew", per_line=True),
    _V("plot_phase_tensor_skewmap", "Phase-tensor skew section",
       "Dimensionality & skew", per_line=True),
    _V("plot_dim_map", "Dimensionality map", "Dimensionality & skew",
       help="1-D / 2-D / 3-D class of each station at one period."),
    _V("plot_dim_occupancy_area", "Dimensionality by period",
       "Dimensionality & skew",
       help="Share of 1-D, 2-D and 3-D stations through the periods."),
    # static shift
    _V("ss_qc_psection", "Static-shift section", "Static shift",
       per_line=True),
    _V("ss_qc_profile", "Static-shift profile", "Static shift",
       per_line=True),
    _V("ss_qc_station_curves", "Static-shift station curves",
       "Static shift", station=True),
    _V("plot_ss_radar", "Static-shift radar", "Static shift", station=True),
    # distortion
    _V("plot_ns_detection", "Near-surface detection", "Distortion & source"),
    _V("plot_distortion_radar", "Distortion radar", "Distortion & source",
       help="Galvanic-distortion indicators for a few stations."),
    _V("plot_consistency_fan", "Consistency fan", "Distortion & source",
       station=True, help="Bootstrap spread of one station's response."),
    _V("plot_overprint_section", "Source overprint section",
       "Distortion & source", per_line=True),
    _V("plot_field_zones", "Field zones", "Distortion & source",
       per_line=True),
    # strike
    _V("plot_strike_profile", "Strike profile", "Strike", per_line=True),
    _V("plot_strike_ribbon", "Strike ribbon", "Strike", per_line=True),
    _V("plot_strike_stability_bands", "Strike stability", "Strike",
       help="Strike estimates by period band and method."),
    _V("plot_strike_mapsticks", "Strike map sticks", "Strike"),
    _V("plot_strike_rose_by_line", "Strike rose by line", "Strike",
       lines_arg="groups", help="One strike rose per survey line."),
)


def view(key: str) -> QCView:
    return next(v for v in VIEWS if v.fn == key)


# ── options ───────────────────────────────────────────────────────────────
_HIDE = {"station", "station_order", "line_labels", "groups", "group_key",
         "lines", "stations", "figsize", "title", "weights"}
_CHOICE_LABELS = {"xy": "xy", "yx": "yx", "xx": "xx", "yy": "yy",
                  "det": "det", "pt": "PT", "rpca": "RPCA", "emap": "EMAP",
                  "logperiod": "log period", "inv_iqr": "1 / IQR"}

# string options the functions validate against a fixed set (their own
# error messages); keyed by name, or (name, default) where the same name
# means different things
_TEXT_CHOICES: dict = {
    "annotate": ("auto", "true", "false"),
    "annotate_stations": ("auto", "true", "false"),
    "coordinate_system": ("auto", "geographic", "projected"),
    "interpolation": ("linear", "cubic"),
    "line_statistic": ("median", "mean"),
    "map_aspect": ("auto", "equal", "geographic"),
    "metric": ("confidence", "coverage", "uncertainty", "offdiag",
               "diagonal", "phase", "spatial"),
    "order": ("worst", "best"),
    "shade_mode": ("score", "full", "none"),
    "section": ("dynamic", "pseudosection", "compact", "publication",
                "dashboard", "inversion"),
    "station_preset": ("pseudosection", "inversion", "survey"),
    "theta_axis": ("logperiod", "period"),
    ("method", "composite"): ("composite", "presence"),
    ("method", "consensus"): ("consensus", "sweep", "pt"),
    ("method", "pipeline"): ("pipeline", "notch", "smooth", "rpca",
                             "spatial", "hampel", "emap"),
    # the harmonic waterfall passes mains_hz: only these two accept it
    ("method", "notch"): ("notch", "pipeline"),
    ("mode", "auto"): ("auto", "scatter", "route", "contour"),
    ("mode", "contour"): ("contour", "scatter"),
}
_BOOL_OR_AUTO = {"annotate", "annotate_stations"}
_COSMETIC = ("color", "colors", "_fmt", "linestyle", "linestyles",
             "linewidths", "_kws", "colorbar_label", "edgecolor",
             "grid_shape", "count_kws")


def _label(name: str, spec_label: str) -> str:
    return {"ci_hi": "Safe threshold", "ci_lo": "Recoverable threshold",
            "snr_thresh": "SNR threshold", "skew_th": "Skew threshold",
            "vmax": "Colour max", "spacing_m": "Station spacing",
            "source_offset": "Source offset"}.get(name, spec_label)


def _choice_label(c) -> str:
    c = str(c)
    return _CHOICE_LABELS.get(c, c.replace("_", " ").capitalize())


def _float_field(name, label, default, section):
    d = float(default)
    lo = -1e9 if name.endswith(("_deg", "rotate")) or d < 0 else 0.0
    unit = {"spacing_m": "m", "source_offset": "m"}.get(name, "")
    if 0.0 <= d <= 1.0 and name.startswith(("ci_", "alpha", "q", "frac",
                                             "min_frac")):
        return Field(name, label, "float", d, 0.0, 1.0, 0.01, 3,
                     section=section)
    step = 10 ** math.floor(math.log10(abs(d))) if d else 0.1
    return Field(name, label, "float", d, lo, 1e9, step,
                 max(2, 4 if abs(d) < 1 else 2), section=section, unit=unit)


def options_for(v: QCView, stations=(), method: str | None = None
                ) -> list[Field]:
    """Form fields for *v* (``section`` is "Analysis" or "View")."""
    out: list[Field] = []
    for spec in qc_parameter_specs(v.fn, method):
        name = spec.name
        if name in _HIDE:
            continue
        section = "View" if spec.view else "Analysis"
        label = _label(name, spec.label)
        tip = f"{name} (default: {spec.default!r})"
        cosmetic = name.endswith(_COSMETIC) or "color" in name
        if cosmetic:
            section = "View"
        fixed = _TEXT_CHOICES.get((name, spec.default),
                                  _TEXT_CHOICES.get(name))
        if (fixed and isinstance(spec.default, str) and not spec.choices
                and spec.kind in {"str", "optional_str"}):
            opts = tuple(dict.fromkeys((spec.default, *fixed)))
            out.append(Field(name, label, "choice", spec.default,
                             choices=tuple((c, _choice_label(c))
                                           for c in opts),
                             section=section, help=tip))
            continue
        if name == "cmap" or name.endswith("_cmap"):
            cur = spec.default if isinstance(spec.default, str) else ""
            names = colormap_choices(cur or None, first=(cur,) if cur else ())
            out.append(Field(name, label.replace("Cmap", "colour map"),
                             "choice", cur, choices=tuple(
                                 (c, colormap_label(c)) for c in names)
                             if cur else (("", "Default"),) + tuple(
                                 (c, colormap_label(c)) for c in names),
                             section=section, help=tip))
        elif spec.kind == "bool":
            out.append(Field(name, label, "bool", bool(spec.default),
                             section=section, help=tip))
        elif spec.kind == "int":
            out.append(Field(name, label, "int", int(spec.default),
                             1 if name in {"half_window", "bins", "n_bins",
                                           "n_bands"} else 0,
                             100000, 1, section=section, help=tip))
        elif spec.kind == "float":
            out.append(_float_field(name, label, spec.default, section))
        elif spec.choices:
            labels = (QC_STATIC_SHIFT_METHODS if name == "method"
                      and v.fn in QC_STATIC_SHIFT_PLOTS else {})
            choices = tuple((c, labels.get(c, _choice_label(c)))
                            for c in spec.choices)
            if spec.default is None:
                choices = (("", "Default"),) + tuple(
                    c for c in choices if c[0] != "auto")
            default = "" if spec.default is None else spec.default
            out.append(Field(name, label, "choice", default, choices=choices,
                             section=section, help=tip))
        else:  # text: optional numbers, lists, pairs, free strings
            text = ("" if spec.default is None else
                    ", ".join(str(x) for x in spec.default)
                    if isinstance(spec.default, (tuple, list))
                    else str(spec.default))
            hint = {"sequence": "comma-separated", "optional_sequence":
                    "comma-separated; blank = automatic",
                    "optional_float": "number; blank = automatic",
                    "optional_int": "whole number; blank = automatic",
                    "optional_str": "blank = automatic"}.get(spec.kind, "")
            out.append(Field(name, label, "text", text, section=section,
                             advanced=cosmetic,
                             help=f"{tip} — {hint}" if hint else tip))
    return out


_KINDS: dict[tuple[str, str], str] = {}


def _kind(fn: str, name: str, method: str | None) -> str:
    key = (fn, name)
    if key not in _KINDS:
        for spec in qc_parameter_specs(fn, method):
            _KINDS[(fn, spec.name)] = spec.kind
    return _KINDS.get(key, "str")


def to_kwargs(v: QCView, values: dict, method: str | None = None
              ) -> tuple[dict, str | None]:
    """Form values -> keyword arguments; ``(kwargs, error message)``."""
    kw: dict = {}
    for name, raw in values.items():
        kind = _kind(v.fn, name, method)
        try:
            if name in _BOOL_OR_AUTO and raw in ("true", "false"):
                kw[name] = raw == "true"
                continue
            if isinstance(raw, str):
                text = raw.strip().strip("()[]")
                if not text:
                    if kind in {"sequence"}:
                        return {}, f"{name.replace('_', ' ')} cannot be empty."
                    continue  # blank = the function's own default
                if kind in {"optional_float", "float"}:
                    kw[name] = float(text)
                elif kind in {"optional_int", "int"}:
                    kw[name] = int(float(text))
                elif kind in {"sequence", "optional_sequence"}:
                    parts = [p.strip() for p in text.split(",") if p.strip()]
                    try:
                        kw[name] = tuple(float(p) for p in parts)
                    except ValueError:
                        kw[name] = tuple(parts)
                else:
                    kw[name] = text
            else:
                kw[name] = raw
        except ValueError:
            return {}, f"“{name.replace('_', ' ')}” needs a number."
    if "source_offset" in kw and kw["source_offset"] <= 0:
        return {}, "Source offset must be greater than zero metres."
    for a, b, what in (("ci_lo", "ci_hi", "The recoverable threshold must "
                        "be below the safe threshold."),
                       ("near_threshold", "far_threshold",
                        "Near threshold must be below far threshold.")):
        if a in kw and b in kw and kw[a] >= kw[b]:
            return {}, what
    return kw, None


# ── survey helpers ────────────────────────────────────────────────────────
def _items(sites):
    from pycsamt.emtools._core import _iter_items, _name

    return [(str(_name(ed, i)), ed) for i, ed in enumerate(_iter_items(sites))]


def station_names(sites) -> list[str]:
    try:
        return [n for n, _ed in _items(sites)] if sites is not None else []
    except Exception:
        return []


def line_groups_of(sites, lines: dict | None) -> dict[str, list[str]]:
    """``{line: [stations]}`` in survey order; one group without lines."""
    names = station_names(sites)
    if not lines:
        return {"": names} if names else {}
    out: dict[str, list[str]] = {}
    for n in names:
        out.setdefault(str(lines.get(n) or "unassigned"), []).append(n)
    return out


def subset(sites, names):
    """The stations *names* of *sites* (order kept) as a Sites-like list."""
    from pycsamt.emtools._core import ensure_sites

    keep = set(names)
    return ensure_sites([ed for n, ed in _items(sites) if n in keep],
                        recursive=False, verbose=0)


# ── rendering ─────────────────────────────────────────────────────────────
def _has_ax(fn_name: str) -> bool:
    import pycsamt.emtools as et

    try:
        return "ax" in inspect.signature(getattr(et, fn_name)).parameters
    except (TypeError, ValueError, AttributeError):
        return False


def render(v: QCView, sites, kwargs: dict, *, lines: dict | None = None,
           line: str = "", layout: str = "panels", station: str = "",
           ctrl: QCController | None = None):
    """Draw *v*; returns ``(figure, unavailable)``.

    *line* restricts to one line ("" = every line in scope).  For a
    per-line view over several lines, *layout* ``"panels"`` draws one
    stacked panel per line; ``"together"`` draws them on one axes (the
    confidence profile colours each line).
    """
    import matplotlib.pyplot as plt

    ctrl = ctrl or QCController()
    ctrl.dark = False  # figures are for publication: white
    groups = line_groups_of(sites, lines)
    if line and line in groups:
        sites = subset(sites, groups[line])
        groups = {line: groups[line]}
    kw = dict(kwargs)
    if v.station:
        kw["station"] = station or None
    if v.lines_arg == "groups":
        kw["groups"] = {(k or "survey"): list(n) for k, n in groups.items()}
    if v.fn == "plot_confidence_profile":
        kw.setdefault("annotate_low", False)
        if lines and len(groups) > 1 and layout == "together":
            kw["lines"] = {n: k for k, ns in groups.items() for n in ns}
    many = (v.per_line and len(groups) > 1 and layout == "panels"
            and _has_ax(v.fn))
    if many:
        ctrl.set_sites(sites)
        pre = ctrl._preflight(v.fn, source_offset=kw.get("source_offset"))
        if pre is not None:
            return plt.figure(), pre
    if not many:
        ctrl.set_sites(sites)
        fig = plt.figure(figsize=(11, 6))
        out = ctrl.draw(v.fn, _has_ax(v.fn), fig, **kw)
        if out is not None:
            plt.close(fig)
            fig = out
        return fig, ctrl.last_unavailable
    # one panel per line
    names = list(groups)
    # a grid, not a tall stack: the canvas fits the figure to the window,
    # so five stacked panels were each a few pixels high
    ncols = 1 if len(names) <= 2 else 2
    nrows = math.ceil(len(names) / ncols)
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(11 if ncols == 1 else 14,
                                      max(3.4 * nrows, 5)))
    for ax in axes.flat[len(names):]:
        ax.set_visible(False)
    import pycsamt.emtools as et

    fn = getattr(et, v.fn)
    params = inspect.signature(fn).parameters
    title = ""
    for k, (ax, name) in enumerate(zip(axes.flat, names)):
        try:
            fn(subset(sites, groups[name]), ax=ax, verbose=0,
               **{key: val for key, val in kw.items() if key in params})
            # one title for the figure, a line tag in each panel, x label
            # and legend once: stacked panels otherwise collide
            title = title or ax.get_title()
            ax.set_title("")
            ax.text(0.995, 0.03, f"Line {name}", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=9, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white",
                              ec="0.6", alpha=0.9), zorder=10)
            if k < len(names) - ncols:  # x label on the bottom row only
                ax.set_xlabel("")
            if k > 0 and ax.get_legend() is not None:
                ax.get_legend().remove()
        except Exception as exc:
            ax.cla()
            ax.axis("off")
            ax.text(0.5, 0.5, f"Line {name}: {exc}", ha="center",
                    va="center", transform=ax.transAxes, fontsize=9)
    if title:
        fig.suptitle(title, fontsize=11)
    try:
        fig.tight_layout()
    except Exception:
        pass
    return fig, None


# ── summary ───────────────────────────────────────────────────────────────
def scorecard(sites, lines: dict | None = None, *, ci_hi: float = 0.95,
              ci_lo: float = 0.85, min_frac_ok: float = 0.6,
              min_snr: float = 2.0, max_skew: float = 6.0):
    """One row per station with a pass / warn / fail status.

    The status follows the confidence ratio of Kouadio et al. (2024) --
    the share of frequencies with a valid tensor, which is what the
    *ci_hi* / *ci_lo* thresholds were defined for -- plus the flags of
    :func:`~pycsamt.emtools.qc.qc_flags`:

    * ``fail`` -- ratio below *ci_lo*, or coverage below *min_frac_ok*;
    * ``warn`` -- ratio below *ci_hi*, or low SNR / high skew;
    * ``pass`` -- otherwise.

    The composite score (coverage, uncertainty, off-diagonal consistency,
    diagonal leakage, phase smoothness) is reported beside it: it runs
    lower by construction (about 0.7 on clean data), so it ranks stations
    rather than gating them.
    """
    import pandas as pd

    from pycsamt.emtools.qc import qc_flags, station_confidence_table

    groups = line_groups_of(sites, lines)
    frames = []
    for line, names in groups.items():
        sub = subset(sites, names) if lines else sites
        conf = station_confidence_table(sub, method="composite", api=False)
        pres = station_confidence_table(sub, method="presence", api=False)
        pres = (pres.to_pandas() if hasattr(pres, "to_pandas")
                else pd.DataFrame(pres))
        flags = qc_flags(sub, min_frac_ok=min_frac_ok, min_snr_med=min_snr,
                         max_skew_med=max_skew)
        flags = (flags.to_pandas() if hasattr(flags, "to_pandas")
                 else pd.DataFrame(flags))
        conf = (conf.to_pandas() if hasattr(conf, "to_pandas")
                else pd.DataFrame(conf))
        conf = pd.merge(
            pres[["station", "confidence"]].rename(
                columns={"confidence": "ratio"}),
            conf[["station", "confidence", "uncertainty", "offdiag",
                  "phase"]], on="station", how="outer")
        df = pd.merge(conf,
                      flags[["station", "frac_ok", "snr_med", "skew_med",
                             "n_freq", "flags"]],
                      on="station", how="outer")
        df.insert(0, "line", line)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)

    def status(r) -> str:
        c = r["ratio"]
        f = str(r.get("flags") or "")
        if (np.isfinite(c) and c < ci_lo) or "low_coverage" in f:
            return "fail"
        if (np.isfinite(c) and c < ci_hi) or f:
            return "warn"
        return "pass"

    df["status"] = df.apply(status, axis=1)
    df["flags"] = df["flags"].fillna("").str.replace(",", ", ")
    order = ["line", "station", "status", "ratio", "confidence", "frac_ok",
             "snr_med", "skew_med", "uncertainty", "offdiag", "phase",
             "n_freq", "flags"]
    df = df[order].rename(columns={
        "line": "Line", "station": "Station", "status": "Status",
        "ratio": "Confidence ratio", "confidence": "Composite score",
        "frac_ok": "Coverage",
        "snr_med": "SNR (median)", "skew_med": "|Skew| median (°)",
        "uncertainty": "Uncertainty score", "offdiag": "Off-diagonal score",
        "phase": "Phase score", "n_freq": "Frequencies", "flags": "Flags"})
    if not lines:
        df = df.drop(columns="Line")
    return df
