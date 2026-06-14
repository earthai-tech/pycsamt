# -*- coding: utf-8 -*-
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

from typing import Optional, Tuple  # Tuple kept for _auto_rho_bounds

import matplotlib.pyplot as plt
import numpy as np


# ── Theme style dicts (same palette as PlotController) ────────────────────────

_DARK = dict(
    facecolor="#181825", labelcolor="#cdd6f4", title_color="#cdd6f4",
    tick_params={"colors": "#a6adc8", "labelsize": 8},
    spines_color="#45475a", grid_color="#313244",
    grid_alpha=0.35, grid_linestyle="--",
)
_LIGHT = dict(
    facecolor="#eff1f5", labelcolor="#4c4f69", title_color="#4c4f69",
    tick_params={"colors": "#6c6f85", "labelsize": 8},
    spines_color="#bcc0cc", grid_color="#ccd0da",
    grid_alpha=0.5, grid_linestyle="--",
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
    ax.grid(True, color=s["grid_color"], alpha=s["grid_alpha"],
            linestyle=s["grid_linestyle"])
    fig = ax.get_figure()
    if fig is not None:
        fig.patch.set_facecolor("#1e1e2e" if dark else "#e6e9ef")


def _annotate_empty(ax, msg: str = "No data") -> None:
    ax.text(0.5, 0.5, msg, transform=ax.transAxes,
            ha="center", va="center", fontsize=11, color="#585b70")


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
    ("QC Quicklook (multi-panel)",    "plot_qc_quicklook",               False),
    ("Coverage quality heatmap",      "plot_coverage_quality_heatmap",   True),
    ("Frequency confidence section",  "plot_frequency_confidence_psection", True),
    ("Confidence band summary",       "plot_confidence_band_summary",    True),
    ("Confidence profile",            "plot_confidence_profile",         True),
]

COVERAGE_PLOTS: list = [
    ("Coverage (per-site presence)",  "plot_coverage",                   True),
    ("Coverage pseudosection",        "plot_coverage_psection",          True),
    ("Polar coverage",                "plot_polar_coverage",             True),
]

NOISE_PLOTS: list = [
    ("SNR histogram",                 "plot_snr_hist",                   True),
    ("ΔZ off-diagonal psection",      "nr_qc_delta_offdiag_psection",   True),
    ("Harmonic waterfall",            "nr_qc_harmonic_waterfall",        True),
    ("SNR gain profile",              "nr_qc_snr_gain_profile",          True),
    ("Station off-diagonal curves",   "nr_qc_station_offdiag_curves",   True),
]

SKEW_DIM_PLOTS: list = [
    ("Skew traffic psection",         "plot_skew_traffic_psection",      True),
    ("Skew percentile ribbon",        "plot_skew_percentile_ribbon",     True),
    ("Dimensionality psection",       "plot_dimensionality_psection",    True),
    ("Dimensionality grid",           "plot_dimensionality_grid",        True),
    ("Phase tensor skew map",         "plot_phase_tensor_skewmap",       True),
    ("Dim confidence grid",           "plot_dim_confidence_grid",        True),
]

STATIC_SHIFT_PLOTS: list = [
    ("SS QC pseudosection",           "ss_qc_psection",                  True),
    ("SS QC profile",                 "ss_qc_profile",                   True),
    ("SS QC station curves",          "ss_qc_station_curves",            True),
    ("SS radar",                      "plot_ss_radar",                   True),
]

DISTORTION_PLOTS: list = [
    ("Near-surface detection",        "plot_ns_detection",               True),
    ("Near-surface overprint section","plot_overprint_section",          True),
    ("Field zones",                   "plot_field_zones",                True),
    ("Strike profile",                "plot_strike_profile",             True),
    ("Strike ribbon",                 "plot_strike_ribbon",              True),
]

ALL_GROUPS: list = [
    ("Overview",      OVERVIEW_PLOTS),
    ("Coverage",      COVERAGE_PLOTS),
    ("Noise / SNR",   NOISE_PLOTS),
    ("Skew / Dim",    SKEW_DIM_PLOTS),
    ("Static Shift",  STATIC_SHIFT_PLOTS),
    ("Distortion",    DISTORTION_PLOTS),
]


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

    def set_sites(self, sites) -> None:
        self._sites = sites

    def clear(self) -> None:
        self._sites = None

    # ── Main entry-point ──────────────────────────────────────────────────────

    def draw(self, fn_name: str, has_ax: bool, fig) -> Optional[object]:
        """
        Render *fn_name*.

        Parameters
        ----------
        fn_name : str
            Name of the emtools function to call.
        has_ax : bool
            True  → draw into *fig* in-place; return None.
            False → function creates its own multi-axes figure; that figure
                    is styled and returned so the caller can hand it to
                    MplCanvas.show_figure(), preserving all axes content.
        fig : matplotlib.figure.Figure
            Target figure (mutated when has_ax=True or on error;
            ignored when has_ax=False and a figure is returned).

        Returns
        -------
        Figure or None
            New Figure to display (has_ax=False), or None (has_ax=True).
        """
        import pycsamt.emtools as et

        fig.clear()

        if self._sites is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, "Load survey data first")
            _style_ax(ax, self.dark)
            return None

        fn = getattr(et, fn_name, None)
        if fn is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Function not found: {fn_name}")
            _style_ax(ax, self.dark)
            return None

        try:
            if fn_name in _SPECIAL_DISPATCHERS:
                _SPECIAL_DISPATCHERS[fn_name](self, fn, fig)
            elif has_ax:
                ax = fig.add_subplot(111)
                fn(self._sites, ax=ax, verbose=0)
            else:
                src_fig = self._call_figure_fn(fn)
                if src_fig is None:
                    ax = fig.add_subplot(111)
                    _annotate_empty(ax, "No figure produced")
                else:
                    _style_all_axes(src_fig, self.dark)
                    try:
                        src_fig.tight_layout(pad=1.2)
                    except Exception:
                        pass
                    return src_fig
        except Exception as exc:
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

    def _call_figure_fn(self, fn) -> Optional[object]:
        """Call a multi-axes function and return the Figure it creates."""
        before = set(plt.get_fignums())
        result = fn(self._sites, verbose=0)
        after  = set(plt.get_fignums())
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
) -> Tuple[float, float]:
    """
    Compute data-driven (q_lo, q_hi) as percentiles of observed ρ_a(xy)
    across all loaded stations.  Used when the QC panel calls
    plot_polar_coverage without external model bounds.
    """
    try:
        from pycsamt.emtools.diag import (
            _iter_items, _unwrap, _get_z_block, _rho_a_from_z
        )
        rho_vals = []
        for i, ed in enumerate(_iter_items(sites)):
            ed = _unwrap(ed)
            _, z_block, freqs = _get_z_block(ed)
            if z_block is None or freqs is None or freqs.size == 0:
                continue
            try:
                rho = _rho_a_from_z(z_block, freqs, "xy")
                rho_vals.extend(float(r) for r in rho if np.isfinite(r) and r > 0)
            except Exception:
                continue
        if rho_vals:
            arr = np.asarray(rho_vals)
            return float(np.percentile(arr, lo_pct)), float(np.percentile(arr, hi_pct))
    except Exception:
        pass
    return 0.5, 5000.0


def _dispatch_polar_coverage(ctrl: "QCController", fn, fig) -> None:
    """Special dispatcher for plot_polar_coverage — adds polar axes and auto-bounds."""
    q_lo, q_hi = _auto_rho_bounds(ctrl._sites)
    ax = fig.add_subplot(111, projection="polar")
    fn(ctrl._sites, q_lo, q_hi, ax=ax, verbose=0)


# ── Registry of functions that need non-standard dispatch ─────────────────────
# Maps fn_name → callable(ctrl, fn, fig)

_SPECIAL_DISPATCHERS: dict = {
    "plot_polar_coverage": _dispatch_polar_coverage,
}
