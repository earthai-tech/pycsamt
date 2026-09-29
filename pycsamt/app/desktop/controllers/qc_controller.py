# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
QCController — drives all QC panel redraws.

No Qt imports.  Each draw method accepts a matplotlib Figure or Axes,
delegates to pycsamt.emtools, and applies the active theme style.
Functions that create their own multi-axes figure are handled via
_render_to_fig(); single-axes functions use the provided axes directly.
"""

from __future__ import annotations

import inspect
from dataclasses import dataclass, replace

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class QCUnavailableResult:
    """User-facing explanation for a QC result that cannot be produced."""

    title: str
    reason: str
    guidance: str = ""
    kind: str = "unavailable"


@dataclass(frozen=True)
class QCParameterSpec:
    """Description used to build one dynamic QC parameter editor."""

    name: str
    label: str
    default: object
    kind: str
    view: bool = False
    choices: tuple[str, ...] = ()
    method_specific: bool = False


QC_STATIC_SHIFT_PLOTS = {
    "ss_qc_psection", "ss_qc_profile", "ss_qc_station_curves",
    "ss_comparison_psection",
}
QC_STATIC_SHIFT_METHODS = {
    "ama": "AMA — adaptive moving average",
    "loess": "LOESS — local regression",
    "bilateral": "Bilateral — edge-preserving filter",
    "refmedian": "Reference median",
}
_SS_ESTIMATORS = {
    "ama": "correct_ss_ama",
    "loess": "estimate_ss_loess",
    "bilateral": "estimate_ss_bilateral",
    "refmedian": "estimate_ss_refmedian",
}


_HIDDEN_PARAMETERS = {
    "sites",
    "ax",
    "axes",
    "recursive",
    "on_dup",
    "strict",
    "verbose",
    "return_sites",
    "station_style",
    "colors",
    "marker",
    "ms",
    "lw",
    "ls",
    "eps",
}

_VIEW_PARAMETERS = {
    "figsize",
    "cmap",
    "title",
    "axis",
    "axis_y",
    "theta_axis",
    "period_axis",
    "log_y",
    "log_radius",
    "log_color",
    "colorbar",
    "show_colorbar",
    "show_mask",
    "show_ss",
    "show_errorbars",
    "shade_recoverable",
    "shade_mode",
    "annotate_low",
    "annotate_low_step",
    "station_labels",
    "station_label_step",
    "station_preset",
    "section",
    "alpha_by",
    "n_freq_ticks",
    "cbar_ticks",
    "smart_ylim",
    "ylim",
    "vlim",
    "contour_kr",
    "kr_levels",
    "contour_beta",
    "beta_levels",
    "segmented_colors",
    "fill_between",
}

_PARAMETER_CHOICES: dict[tuple[str, str] | str, tuple[str, ...]] = {
    "axis": ("period", "frequency"),
    "axis_y": ("logperiod", "period", "frequency"),
    "sort_by": ("auto", "chainage", "input", "station", "latitude", "longitude"),
    "rho_comp": ("xy", "yx", "xx", "yy", "det"),
    "agg": ("median", "mean"),
    "robust": ("median", "mean"),
    "robust_freq": ("median", "mean"),
    "robust_overall": ("median", "mean"),
    "summary": ("median", "mean"),
    "rotate": ("pt", "none"),
    "rotate_stat": ("median", "mean"),
    "radial": ("log10rho", "rho", "phase"),
    "weights": ("tri", "uniform"),
    ("plot_frequency_confidence_psection", "method"): (
        "composite", "presence", "error", "snr"
    ),
    ("plot_confidence_band_summary", "method"): (
        "composite", "presence", "error", "snr"
    ),
    ("plot_confidence_profile", "method"): (
        "presence", "composite", "error", "snr"
    ),
    ("plot_coverage_psection", "metric"): ("presence", "error", "snr"),
    ("plot_coverage_psection", "alpha_by"): ("none", "confidence"),
    ("plot_strike_profile", "method"): ("consensus", "pt", "swift"),
    ("plot_strike_ribbon", "method"): ("sweep", "consensus", "pt"),
}


def qc_parameter_specs(
    fn_name: str, method: str | None = None,
) -> list[QCParameterSpec]:
    """Return safe, editable public parameters for an emtools QC function."""
    import pycsamt.emtools as et

    fn = getattr(et, fn_name, None)
    if fn is None:
        return []
    specs: list[QCParameterSpec] = []
    for parameter in inspect.signature(fn).parameters.values():
        name = parameter.name
        if (
            name in _HIDDEN_PARAMETERS
            or parameter.kind
            in (inspect.Parameter.VAR_POSITIONAL, inspect.Parameter.VAR_KEYWORD)
            or parameter.default is inspect.Parameter.empty
        ):
            continue
        default = parameter.default
        annotation = str(parameter.annotation).lower()
        if isinstance(default, bool):
            kind = "bool"
        elif isinstance(default, int):
            kind = "int"
        elif isinstance(default, float):
            kind = "float"
        elif isinstance(default, str):
            kind = "str"
        elif isinstance(default, (tuple, list)):
            kind = "sequence"
        elif default is None:
            if "float" in annotation:
                kind = "optional_float"
            elif "int" in annotation:
                kind = "optional_int"
            elif "tuple" in annotation or "list" in annotation:
                kind = "optional_sequence"
            elif "str" in annotation or name in {"station", "sort_by"}:
                kind = "optional_str"
            elif name == "source_offset":
                kind = "optional_float"
            else:
                continue
        else:
            continue
        choices = _PARAMETER_CHOICES.get(
            (fn_name, name), _PARAMETER_CHOICES.get(name, ())
        )
        if fn_name in QC_STATIC_SHIFT_PLOTS and name == "method":
            choices = tuple(QC_STATIC_SHIFT_METHODS)
        label = name.replace("_", " ").strip().title()
        if name == "source_offset":
            label = "Source offset (m)"
        label = {
            "half_window": "Neighbours per side",
            "poly": "Polynomial degree",
            "it": "Robust iterations",
            "sig_dist": "Distance bandwidth",
            "sig_val": "Value bandwidth",
            "smooth_sites": "Station smoothing window",
            "pband": "Period band (s)",
            "max_skew": "Maximum skew (°)",
            "robust_freq": "Frequency statistic",
            "robust_overall": "Overall statistic",
            "summary": "Summary statistic",
        }.get(name, label)
        specs.append(
            QCParameterSpec(
                name=name,
                label=label,
                default=default,
                kind=kind,
                view=name in _VIEW_PARAMETERS,
                choices=choices,
            )
        )
    if fn_name in QC_STATIC_SHIFT_PLOTS:
        from pycsamt.emtools.ss import _ALLOWED

        method = method or "ama"
        if method in _SS_ESTIMATORS:
            existing = {spec.name for spec in specs}
            specs.extend(
                replace(spec, method_specific=True)
                for spec in qc_parameter_specs(_SS_ESTIMATORS[method])
                if spec.name in _ALLOWED[method]
                and spec.name not in existing
            )
    return specs

# ── Theme style dicts (same palette as PlotController) ────────────────────────

_DARK = dict(
    facecolor="#181825",
    labelcolor="#cdd6f4",
    title_color="#cdd6f4",
    tick_params={"colors": "#a6adc8", "labelsize": 8},
    spines_color="#45475a",
    grid_color="#313244",
    grid_alpha=0.35,
    grid_linestyle="--",
)
_LIGHT = dict(
    facecolor="#ffffff",
    labelcolor="#4c4f69",
    title_color="#4c4f69",
    tick_params={"colors": "#6c6f85", "labelsize": 8},
    spines_color="#bcc0cc",
    grid_color="#ccd0da",
    grid_alpha=0.5,
    grid_linestyle="--",
)


def _style_ax(ax, dark: bool = True) -> None:
    s = _DARK if dark else _LIGHT
    ax.set_facecolor(s["facecolor"])
    ax.title.set_color(s["title_color"])
    ax.xaxis.label.set_color(s["labelcolor"])
    ax.yaxis.label.set_color(s["labelcolor"])
    ax.tick_params(**s["tick_params"])
    for sp in ax.spines.values():
        sp.set_edgecolor(s["spines_color"])
    ax.grid(
        True,
        color=s["grid_color"],
        alpha=s["grid_alpha"],
        linestyle=s["grid_linestyle"],
    )
    fig = ax.get_figure()
    if fig is not None:
        fig.patch.set_facecolor("#1e1e2e" if dark else "#ffffff")


def _annotate_empty(ax, msg: str = "No data") -> None:
    """Clear *ax* and show a centred placeholder message (see the
    identical helper in plot_controller.py for the full rationale)."""
    ax.cla()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.5,
        0.5,
        msg,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        color="#585b70",
    )


def _style_all_axes(fig, dark: bool) -> None:
    for ax in fig.get_axes():
        try:
            _style_ax(ax, dark)
        except Exception:
            pass


# ── Plot-function catalogue ────────────────────────────────────────────────────
# Each entry: (label, fn_name, has_ax_param)
#   has_ax_param=True  → called as fn(sites, ax=ax, verbose=0)
#   has_ax_param=False → called as fn(sites, verbose=0), returns Figure

OVERVIEW_PLOTS: list = [
    ("QC Quicklook (multi-panel)", "plot_qc_quicklook", False),
    ("Coverage quality heatmap", "plot_coverage_quality_heatmap", True),
    (
        "Frequency confidence section",
        "plot_frequency_confidence_psection",
        True,
    ),
    ("Confidence band summary", "plot_confidence_band_summary", True),
    ("Confidence profile", "plot_confidence_profile", True),
]

COVERAGE_PLOTS: list = [
    ("Coverage (per-site presence)", "plot_coverage", True),
    ("Coverage pseudosection", "plot_coverage_psection", True),
    ("Polar coverage", "plot_polar_coverage", True),
]

NOISE_PLOTS: list = [
    ("SNR histogram", "plot_snr_hist", True),
    ("ΔZ off-diagonal psection", "nr_qc_delta_offdiag_psection", True),
    ("Harmonic waterfall", "nr_qc_harmonic_waterfall", True),
    ("SNR gain profile", "nr_qc_snr_gain_profile", True),
    ("Station off-diagonal curves", "nr_qc_station_offdiag_curves", True),
]

SKEW_DIM_PLOTS: list = [
    ("Skew traffic psection", "plot_skew_traffic_psection", True),
    ("Skew percentile ribbon", "plot_skew_percentile_ribbon", True),
    ("Dimensionality psection", "plot_dimensionality_psection", True),
    ("Dimensionality grid", "plot_dimensionality_grid", True),
    ("Phase tensor skew map", "plot_phase_tensor_skewmap", True),
    ("Dim confidence grid", "plot_dim_confidence_grid", True),
]

STATIC_SHIFT_PLOTS: list = [
    ("SS QC pseudosection", "ss_qc_psection", True),
    ("SS QC profile", "ss_qc_profile", True),
    ("SS QC station curves", "ss_qc_station_curves", True),
    ("SS radar", "plot_ss_radar", True),
]

DISTORTION_PLOTS: list = [
    ("Near-surface detection", "plot_ns_detection", True),
    ("Near-surface overprint section", "plot_overprint_section", True),
    ("Field zones", "plot_field_zones", True),
    ("Strike profile", "plot_strike_profile", True),
    ("Strike ribbon", "plot_strike_ribbon", True),
]

ALL_GROUPS: list = [
    ("Overview", OVERVIEW_PLOTS),
    ("Coverage", COVERAGE_PLOTS),
    ("Noise / SNR", NOISE_PLOTS),
    ("Skew / Dim", SKEW_DIM_PLOTS),
    ("Static Shift", STATIC_SHIFT_PLOTS),
    ("Distortion", DISTORTION_PLOTS),
]

GROUP_ICONS: dict[str, str] = {
    "Overview": "overview",
    "Coverage": "coverage",
    "Noise / SNR": "frequency-editor",
    "Skew / Dim": "skew",
    "Static Shift": "static-shift",
    "Distortion": "distorsion",
}

PLOT_DESCRIPTIONS: dict[str, str] = {
    "plot_qc_quicklook": (
        "Compact survey overview with station coverage, signal quality, "
        "skew, and broad data-health indicators."
    ),
    "plot_coverage_quality_heatmap": (
        "Shows per-station frequency coverage weighted by impedance error; "
        "bright cells indicate reliable sampled bands."
    ),
    "plot_frequency_confidence_psection": (
        "Pseudosection of frequency confidence across the line, useful for "
        "spotting weak bands before editing or inversion."
    ),
    "plot_confidence_band_summary": (
        "Summarizes confidence by period/frequency band so noisy windows are "
        "easy to compare."
    ),
    "plot_confidence_profile": (
        "Station-by-station confidence profile for locating weak or unstable "
        "sites along the survey line."
    ),
    "plot_coverage": (
        "Per-site data presence view showing how complete each station is "
        "across available frequencies."
    ),
    "plot_coverage_psection": (
        "Station-period coverage pseudosection for checking missing bands, "
        "gaps, and uneven acquisition density."
    ),
    "plot_polar_coverage": (
        "Polar summary of coverage distribution, highlighting whether data "
        "support is balanced over the selected bands."
    ),
    "plot_snr_hist": (
        "Histogram of signal-to-noise ratios for a quick view of survey-wide "
        "noise conditions."
    ),
    "nr_qc_delta_offdiag_psection": (
        "Compares off-diagonal impedance before and after denoising as a "
        "station-period pseudosection."
    ),
    "nr_qc_harmonic_waterfall": (
        "Tracks mains-harmonic reduction by station and harmonic index after "
        "noise filtering."
    ),
    "nr_qc_snr_gain_profile": (
        "Station profile of SNR improvement after the selected noise-removal "
        "pipeline."
    ),
    "nr_qc_station_offdiag_curves": (
        "Detailed off-diagonal impedance curves for one station, with guides "
        "for mains-related contamination."
    ),
    "plot_skew_traffic_psection": (
        "Traffic-light dimensionality section from skew indicators: green "
        "near 1-D, amber transitional, red more complex."
    ),
    "plot_skew_percentile_ribbon": (
        "Line-level skew envelope showing how dimensionality indicators vary "
        "through period bands."
    ),
    "plot_dimensionality_psection": (
        "Classifies station-period cells as 1-D, 2-D, or 3-D using phase "
        "tensor/skew criteria."
    ),
    "plot_dimensionality_grid": (
        "Confidence-aware dimensionality grid that combines structural class "
        "with data reliability."
    ),
    "plot_phase_tensor_skewmap": (
        "Phase-tensor skew pseudosection for locating periods and stations "
        "with strong 3-D behavior."
    ),
    "plot_dim_confidence_grid": (
        "Dimensionality confidence grid; transparent or weak cells warn where "
        "classification is less stable."
    ),
    "ss_qc_psection": (
        "Static-shift correction pseudosection showing how apparent "
        "resistivity changes after correction."
    ),
    "ss_qc_profile": (
        "Station profile of static-shift factors, useful for finding local "
        "near-surface offsets."
    ),
    "ss_qc_station_curves": (
        "Before/after station curves for inspecting static-shift correction "
        "at selected sites."
    ),
    "plot_ss_radar": (
        "Polar radar summary of static-shift magnitude and direction across "
        "stations."
    ),
    "plot_ns_detection": (
        "Near-surface distortion detection summary from static-shift and "
        "response-shape diagnostics."
    ),
    "plot_overprint_section": (
        "Source-overprint pseudosection highlighting station-period cells "
        "where source effects may bias CSAMT responses."
    ),
    "plot_field_zones": (
        "Field-zone diagnostic separating near-field, transition, and "
        "far-field behavior."
    ),
    "plot_strike_profile": (
        "Strike-angle profile by station for checking geoelectric orientation "
        "and rotation consistency."
    ),
    "plot_strike_ribbon": (
        "Period-band strike ribbon showing how preferred strike evolves "
        "across the line."
    ),
}


def describe_plot(fn_name: str) -> str:
    """Return the user-facing description for a QC dashboard plot."""
    return PLOT_DESCRIPTIONS.get(fn_name, "Render this QC diagnostic plot.")


# ── Controller ────────────────────────────────────────────────────────────────


class QCController:
    """
    Pure-Python controller that renders QC plots into matplotlib Figure objects.

    Attributes
    ----------
    dark : bool
        Apply dark-mode styling (default True).
    """

    def __init__(self) -> None:
        self._sites = None
        self.dark: bool = True
        self.last_unavailable: QCUnavailableResult | None = None

    def set_sites(self, sites) -> None:
        self._sites = sites

    def filter_sites(self, station_ids) -> None:
        """Restrict rendering to specific station IDs (empty = use all)."""
        if not station_ids or self._sites is None:
            return
        try:
            from pycsamt.emtools._core import (
                _iter_items,
                ensure_sites,
            )

            all_edis = list(_iter_items(self._sites))
            sel = set(station_ids)
            filtered = [
                ed
                for ed in all_edis
                if getattr(ed, "station", None) in sel
                or getattr(ed, "dataid", None) in sel
            ]
            if filtered:
                self._sites = ensure_sites(filtered, recursive=False)
        except Exception:
            pass

    def clear(self) -> None:
        self._sites = None
        self.last_unavailable = None

    def _unavailable(
        self,
        title: str,
        reason: str,
        guidance: str = "",
        *,
        kind: str = "unavailable",
    ) -> None:
        self.last_unavailable = QCUnavailableResult(
            title=title,
            reason=reason,
            guidance=guidance,
            kind=kind,
        )

    def _preflight(
        self, fn_name: str, source_offset=None
    ) -> QCUnavailableResult | None:
        """Return a precise prerequisite failure before plotting, if known."""
        geometry_plots = {
            "plot_field_zones": ("Field-zone analysis", "|k·r|"),
            "plot_overprint_section": (
                "Source-overprint analysis",
                "the β source-effect ratio",
            ),
        }
        if fn_name not in geometry_plots:
            return None
        try:
            from pycsamt.emtools._core import _iter_items, _name

            if fn_name == "plot_field_zones":
                from pycsamt.emtools.fieldzone import (
                    _resolve_offset,
                    _unwrap,
                )
            else:
                from pycsamt.emtools.source_effects import (
                    _resolve_offset,
                    _unwrap,
                )

            items = list(_iter_items(self._sites))
            analysis_name, quantity = geometry_plots[fn_name]
            if not items:
                return QCUnavailableResult(
                    f"{analysis_name} is unavailable",
                    "The loaded survey does not contain usable stations.",
                    "Load stations with impedance and frequency data.",
                    "missing-data",
                )
            missing_offsets: list[str] = []
            for index, item in enumerate(items):
                edi = _unwrap(item)
                station = _name(edi, index)
                offset = _resolve_offset(edi, source_offset, station)
                if offset is None or not np.isfinite(offset) or offset <= 0:
                    missing_offsets.append(station)
            if missing_offsets:
                missing_count = len(missing_offsets)
                total = len(items)
                verb = "is" if missing_count == 1 else "are"
                return QCUnavailableResult(
                    f"{analysis_name} needs source geometry",
                    (
                        f"{missing_count} of {total} loaded stations {verb} "
                        "missing the "
                        f"source–receiver offset required to compute {quantity}."
                    ),
                    (
                        "Enter Source offset (m) in Controls, or add measured "
                        "source_offset metadata for every station. The plot "
                        "will update automatically without terminal warnings."
                    ),
                    "missing-prerequisite",
                )
        except Exception:
            # Let the normal plotting path provide a more specific result.
            return None
        return None

    @staticmethod
    def _looks_empty(fig) -> bool:
        """Detect a placeholder-only single plot returned by an emtools helper."""
        axes = fig.get_axes()
        if not axes:
            return True
        if len(axes) == 1:
            ax = axes[0]
            messages = [ax.get_title(), *(text.get_text() for text in ax.texts)]
            if any("no data" in message.strip().lower() for message in messages):
                return True
        for ax in axes:
            if ax.lines or ax.images or ax.collections or ax.patches:
                return False
        # A figure with axes but no data artists is still an empty result,
        # whether or not the underlying helper added a "no data" label.
        return True

    # ── Main entry-point ──────────────────────────────────────────────────────

    def draw(
        self,
        fn_name: str,
        has_ax: bool,
        fig,
        figsize: tuple[float, float] | None = None,
        **draw_kwargs,
    ) -> object | None:
        """Render *fn_name* with optional kwargs forwarded to the emtools function.

        Only kwargs whose names appear in the function's signature are passed;
        unknown kwargs are silently discarded so callers can pass a broad dict.

        Parameters
        ----------
        fn_name : str
            Name of the emtools function to call.
        has_ax : bool
            True  → draw into *fig* in-place (single-axes); return None.
            False → function creates its own figure; return that Figure.
        fig : matplotlib.figure.Figure
            Target figure (resized when *figsize* is given).
        figsize : (float, float) or None
            If given, resize *fig* before drawing.
        **draw_kwargs
            Forwarded to the underlying emtools function (filtered by signature).
        """
        import inspect

        import pycsamt.emtools as et

        self.last_unavailable = None
        if figsize:
            fig.set_size_inches(*figsize)
        fig.clear()

        if self._sites is None:
            self._unavailable(
                "Load survey data to begin",
                "No stations are currently available to the QC dashboard.",
                "Load EDI or EMTF-XML files, then select a diagnostic.",
                kind="missing-data",
            )
            ax = fig.add_subplot(111)
            _annotate_empty(ax, "Load survey data first")
            _style_ax(ax, self.dark)
            return None

        fn = getattr(et, fn_name, None)
        if fn is None:
            self._unavailable(
                "This diagnostic is not available",
                f"The plotting function '{fn_name}' is not installed.",
                "Check the application installation or choose another diagnostic.",
                kind="unsupported",
            )
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Function not found: {fn_name}")
            _style_ax(ax, self.dark)
            return None

        preflight = self._preflight(
            fn_name, source_offset=draw_kwargs.get("source_offset")
        )
        if preflight is not None:
            self.last_unavailable = preflight
            # Keep the controller's Figure contract intact for non-Qt callers;
            # the desktop window replaces this fallback with the native view.
            ax = fig.add_subplot(111)
            _annotate_empty(ax, preflight.reason)
            _style_ax(ax, self.dark)
            return None

        # Filter draw_kwargs to only those the function actually accepts
        try:
            valid = set(inspect.signature(fn).parameters.keys())
            if fn_name in QC_STATIC_SHIFT_PLOTS:
                valid.update(
                    spec.name for spec in qc_parameter_specs(
                        fn_name, draw_kwargs.get("method")
                    )
                )
            kw = {
                k: v
                for k, v in draw_kwargs.items()
                if k in valid and v is not None
            }
            if figsize is not None and "figsize" in valid:
                kw["figsize"] = figsize
        except Exception:
            kw = {}

        try:
            if fn_name in _SPECIAL_DISPATCHERS:
                _SPECIAL_DISPATCHERS[fn_name](self, fn, fig, **kw)
            elif has_ax:
                ax = fig.add_subplot(111)
                fn(self._sites, ax=ax, verbose=0, **kw)
            else:
                src_fig = self._call_figure_fn(fn, **kw)
                if src_fig is None:
                    self._unavailable(
                        "No result was produced",
                        "The diagnostic completed without creating a figure.",
                        "Check that the selected stations contain the required data.",
                    )
                    ax = fig.add_subplot(111)
                    _annotate_empty(ax, "No figure produced")
                else:
                    _style_all_axes(src_fig, self.dark)
                    if self._looks_empty(src_fig):
                        self._unavailable(
                            "This diagnostic cannot be computed",
                            (
                                "The loaded stations do not provide usable "
                                "values for this plot."
                            ),
                            (
                                "Check the required channels, frequencies, "
                                "station count, and survey metadata, or choose "
                                "another diagnostic."
                            ),
                            kind="not-computable",
                        )
                        ax = fig.add_subplot(111)
                        _annotate_empty(ax, self.last_unavailable.reason)
                    else:
                        # Figure-returning QC helpers own their layouts and may
                        # contain inset/colorbar axes. Re-running tight_layout()
                        # is incompatible with those axes.
                        return src_fig
            if self._looks_empty(fig):
                self._unavailable(
                    "This diagnostic cannot be computed",
                    "The loaded stations do not provide usable values for this plot.",
                    (
                        "Check the required channels, frequencies, station count, "
                        "and survey metadata, or choose another diagnostic."
                    ),
                    kind="not-computable",
                )
        except Exception as exc:
            self._unavailable(
                "The diagnostic could not be computed",
                str(exc) or f"{fn_name} failed without an error message.",
                "Review the input data and prerequisites, then try again.",
                kind="error",
            )
            fig.clear()
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"{fn_name} error:\n{exc}")

        _style_all_axes(fig, self.dark)
        try:
            fig.tight_layout(pad=1.2)
        except Exception:
            pass
        return None

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _call_figure_fn(self, fn, **kwargs) -> object | None:
        """Call a multi-axes function and return the Figure it creates."""
        before = set(plt.get_fignums())
        result = fn(self._sites, verbose=0, **kwargs)
        after = set(plt.get_fignums())
        if hasattr(result, "get_axes"):
            return result
        new_nums = after - before
        if new_nums:
            return plt.figure(max(new_nums))
        return None


# ── Auto-bounds helper for polar coverage ─────────────────────────────────────


def _auto_rho_bounds(
    sites,
    lo_pct: float = 5.0,
    hi_pct: float = 95.0,
) -> tuple[float, float]:
    """
    Compute data-driven (q_lo, q_hi) as percentiles of observed ρ_a(xy)
    across all loaded stations.  Used when the QC panel calls
    plot_polar_coverage without external model bounds.
    """
    try:
        from pycsamt.emtools.diag import (
            _get_z_block,
            _iter_items,
            _rho_a_from_z,
            _unwrap,
        )

        rho_vals = []
        for _i, ed in enumerate(_iter_items(sites)):
            ed = _unwrap(ed)
            _, z_block, freqs = _get_z_block(ed)
            if z_block is None or freqs is None or freqs.size == 0:
                continue
            try:
                rho = _rho_a_from_z(z_block, freqs, "xy")
                rho_vals.extend(
                    float(r) for r in rho if np.isfinite(r) and r > 0
                )
            except Exception:
                continue
        if rho_vals:
            arr = np.asarray(rho_vals)
            return float(np.percentile(arr, lo_pct)), float(
                np.percentile(arr, hi_pct)
            )
    except Exception:
        pass
    return 0.5, 5000.0


def _dispatch_polar_coverage(ctrl: QCController, fn, fig, **kw) -> None:
    """Special dispatcher for plot_polar_coverage — adds polar axes and auto-bounds."""
    q_lo, q_hi = _auto_rho_bounds(ctrl._sites)
    ax = fig.add_subplot(111, projection="polar")
    fn(ctrl._sites, q_lo, q_hi, ax=ax, verbose=0, **kw)


def _dispatch_polar_ax(ctrl: QCController, fn, fig, **kw) -> None:
    """Special dispatcher for single-panel polar plots."""
    ax = fig.add_subplot(111, projection="polar")
    fn(ctrl._sites, ax=ax, verbose=0, **kw)


# ── Registry of functions that need non-standard dispatch ─────────────────────
# Maps fn_name → callable(ctrl, fn, fig)

_SPECIAL_DISPATCHERS: dict = {
    "plot_polar_coverage": _dispatch_polar_coverage,
    "plot_ss_radar": _dispatch_polar_ax,
}
