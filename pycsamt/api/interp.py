"""Global visual controls for pycsamt.interp hydrogeophysical plots.

Centralises every colour, line-style, figure-geometry, and colourbar
choice used by the five hydro-geophysical plot classes:

* :class:`~pycsamt.interp.plot.PlotHydroSection`
* :class:`~pycsamt.interp.plot.PlotWaterTableProfile`
* :class:`~pycsamt.interp.plot.PlotTimeLapseSection`
* :class:`~pycsamt.interp.plot.PlotUncertaintySection`
* :class:`~pycsamt.interp.plot.PlotUncertaintyProfile`

This module follows the same singleton + preset + context-manager pattern
used by :mod:`pycsamt.api.section`, :mod:`pycsamt.api.style`, and
:mod:`pycsamt.api.station`.

Quick start
-----------
Apply a named preset globally::

    from pycsamt.api import use_interp
    use_interp("publication")           # all subsequent interp plots use it
    use_interp("accessible")            # colorblind-safe palette

Temporarily override for one block::

    from pycsamt.api.interp import PYCSAMT_INTERP
    with PYCSAMT_INTERP.context("dark"):
        PlotHydroSection(result).plot()
        PlotTimeLapseSection(tl).plot()
    # reverts automatically

Fine-tune one field::

    from pycsamt.api import configure_interp
    configure_interp(
        section__cmap_K="plasma",
        section__wt_color="white",
        profile__color_wt="royalblue",
    )

Pass a style directly to one plot (overrides global state)::

    from pycsamt.api.interp import PYCSAMT_INTERP
    fig = PlotHydroSection(result,
                           style=PYCSAMT_INTERP.publication).plot()
"""

from __future__ import annotations

import copy
from collections.abc import Generator
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any

__all__ = [
    # dataclasses
    "HydroSectionStyle",
    "HydroProfileStyle",
    "InterpStyle",
    # container + singleton
    "PyCSAMTInterp",
    "PYCSAMT_INTERP",
    # convenience functions
    "configure_interp",
    "reset_interp",
    "use_interp",
]

_PRESETS = {"default", "publication", "dark", "accessible"}

# ── quantity → cmap_* attribute name ─────────────────────────────────────────
_QUANTITY_CMAP_ATTR: dict[str, str] = {
    "K":          "cmap_K",
    "saturation": "cmap_Sw",
    "porosity":   "cmap_phi",
    "timelapse":  "cmap_timelapse",
    "spread":     "cmap_spread",
    "p50":        "cmap_p50",
    "crossplot":  "cmap_crossplot",
}


# ─────────────────────────────────────────────────────────────────────────────
# Leaf dataclass — 2-D colour sections
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class HydroSectionStyle:
    """Visual parameters for 2-D hydrogeophysical colour sections.

    Governs :class:`~pycsamt.interp.plot.PlotHydroSection`,
    :class:`~pycsamt.interp.plot.PlotTimeLapseSection`, and both panels of
    :class:`~pycsamt.interp.plot.PlotUncertaintySection`.

    Parameters
    ----------
    cmap_K : str
        Colourmap for the hydraulic-conductivity K section (default
        ``'viridis'``).
    cmap_Sw : str
        Colourmap for the water-saturation Sw section (``'RdYlBu'``).
    cmap_phi : str
        Colourmap for the porosity φ section (``'YlOrRd'``).
    cmap_timelapse : str
        Diverging colourmap for Δlog₁₀ρ and ΔSw time-lapse sections
        (``'RdBu_r'``).
    cmap_spread : str
        Colourmap for uncertainty spread (CV or P90–P10) panels
        (``'hot_r'`` — dark = high uncertainty).
    cmap_p50 : str
        Colourmap for the P50 estimate panel of
        :class:`~pycsamt.interp.plot.PlotUncertaintySection`.
    wt_color, wt_lw, wt_ls : str / float / str
        Water-table overlay: colour, line width, line style.
    station_color, station_alpha, station_lw : str / float / float
        Station-tick axvline appearance.
    cb_fraction, cb_pad, cb_fontsize : float / float / int
        Colourbar geometry and tick label size.
    """

    cmap_K:         str   = "viridis"
    cmap_Sw:        str   = "RdYlBu"
    cmap_phi:       str   = "YlOrRd"
    cmap_timelapse: str   = "RdBu_r"
    cmap_spread:    str   = "hot_r"
    cmap_p50:       str   = "viridis"

    wt_color:       str   = "deepskyblue"
    wt_lw:          float = 2.5
    wt_ls:          str   = "--"

    station_color:  str   = "k"
    station_alpha:  float = 0.45
    station_lw:     float = 0.40

    cb_fraction:    float = 0.025
    cb_pad:         float = 0.01
    cb_fontsize:    int   = 8
    nan_color:      str   = "0.88"   # colour for NaN / masked cells in imshow

    # ── petrophysical cross-plot ──────────────────────────────────────────
    cmap_crossplot:      str   = "viridis"   # PlotPetrophysicalCrossPlot
    scatter_alpha:       float = 0.60
    scatter_size:        float = 12.0
    hs_alpha:            float = 0.12        # Hashin-Shtrikman fill opacity
    model_curve_color:   str   = "k"
    model_curve_ls:      str   = "--"
    model_curve_lw:      float = 1.6

    # ── 1-D depth profile ─────────────────────────────────────────────────
    rho_curve_color:     str   = "0.15"      # PlotResistivityDepthProfile
    rho_fill_color:      str   = "steelblue"
    rho_fill_alpha:      float = 0.12
    zone_boundary_color: str   = "0.55"
    zone_boundary_ls:    str   = "--"
    zone_label_fontsize: int   = 7

    # ── convenience methods ───────────────────────────────────────────────

    def cmap_for(self, quantity: str) -> str:
        """Return the colourmap name for *quantity*.

        Parameters
        ----------
        quantity : str
            One of ``'K'``, ``'saturation'``, ``'porosity'``,
            ``'timelapse'``, ``'spread'``, ``'p50'``.
        """
        attr = _QUANTITY_CMAP_ATTR.get(quantity)
        if attr is None:
            raise ValueError(
                f"Unknown quantity {quantity!r}. "
                f"Expected one of {sorted(_QUANTITY_CMAP_ATTR)}."
            )
        return str(getattr(self, attr))

    def wt_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for the water-table ``ax.plot`` call."""
        d: dict[str, Any] = dict(
            color=self.wt_color,
            linewidth=self.wt_lw,
            linestyle=self.wt_ls,
            zorder=5,
        )
        d.update(overrides)
        return d

    def station_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for the station-tick ``ax.axvline`` calls."""
        d: dict[str, Any] = dict(
            color=self.station_color,
            linewidth=self.station_lw,
            alpha=self.station_alpha,
        )
        d.update(overrides)
        return d

    def cb_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments forwarded to ``fig.colorbar``."""
        d: dict[str, Any] = dict(
            fraction=self.cb_fraction,
            pad=self.cb_pad,
        )
        d.update(overrides)
        return d

    def crossplot_scatter_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for the ``ax.scatter`` call in a petrophysical cross-plot."""
        d: dict[str, Any] = dict(
            cmap=self.cmap_crossplot,
            alpha=self.scatter_alpha,
            s=self.scatter_size,
            edgecolors="none",
            zorder=3,
        )
        d.update(overrides)
        return d

    def hs_fill_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for ``ax.fill_between`` of Hashin-Shtrikman bounds."""
        d: dict[str, Any] = dict(alpha=self.hs_alpha, color="0.65", zorder=1)
        d.update(overrides)
        return d

    def model_curve_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for the model-curve ``ax.plot`` call."""
        d: dict[str, Any] = dict(
            color=self.model_curve_color,
            linestyle=self.model_curve_ls,
            linewidth=self.model_curve_lw,
            zorder=5,
        )
        d.update(overrides)
        return d

    def rho_curve_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for the ρ(z) depth-profile curve."""
        d: dict[str, Any] = dict(
            color=self.rho_curve_color,
            linewidth=1.4,
            zorder=3,
        )
        d.update(overrides)
        return d

    def zone_boundary_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for geological zone boundary ``ax.axhline`` calls."""
        d: dict[str, Any] = dict(
            color=self.zone_boundary_color,
            linestyle=self.zone_boundary_ls,
            linewidth=0.7,
        )
        d.update(overrides)
        return d

    def copy(self, **kw: Any) -> HydroSectionStyle:
        """Return a shallow copy with optional field overrides."""
        new = copy.copy(self)
        for k, v in kw.items():
            setattr(new, k, v)
        return new


# ─────────────────────────────────────────────────────────────────────────────
# Leaf dataclass — 1-D profile plots
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class HydroProfileStyle:
    """Visual parameters for 1-D hydrogeophysical profile plots.

    Governs :class:`~pycsamt.interp.plot.PlotWaterTableProfile` and
    :class:`~pycsamt.interp.plot.PlotUncertaintyProfile`.

    Parameters
    ----------
    color_wt : str
        Colour for water-table bars/lines/envelopes (``'steelblue'``).
    color_T : str
        Colour for transmissivity bars/lines/envelopes (``'seagreen'``).
    envelope_alpha : float
        Fill-between transparency for P10–P90 bands (0–1).
    p50_lw : float
        Line width for the P50 central estimate.
    scatter_size : float
        Marker size for station-value scatter dots.
    ref_color, ref_ls, ref_lw : str / str / float
        Reference-depth line style (axhline).
    grid_alpha : float
        Background grid transparency.
    grid_axis : str
        ``'y'``, ``'x'``, or ``'both'``.
    station_color, station_alpha, station_lw : str / float / float
        Station-tick axvline appearance.
    """

    color_wt:       str   = "steelblue"
    color_T:        str   = "seagreen"

    envelope_alpha: float = 0.25
    p50_lw:         float = 2.0
    scatter_size:   float = 14.0

    ref_color:      str   = "tomato"
    ref_ls:         str   = "--"
    ref_lw:         float = 1.4

    grid_alpha:     float = 0.30
    grid_axis:      str   = "y"

    station_color:  str   = "0.6"
    station_alpha:  float = 0.50
    station_lw:     float = 0.35

    # ── Dar-Zarrouk / aquifer characterization ────────────────────────────
    color_TR:              str   = "darkorange"   # PlotAquiferCharacterization
    color_S:               str   = "mediumpurple"
    tr_threshold:          float = 1000.0         # Ω·m² — productive aquifer
    s_threshold_moderate:  float = 1.0            # S — poor/moderate boundary
    s_threshold_good:      float = 5.0            # S — moderate/good boundary

    # ── uncertainty histogram ─────────────────────────────────────────────
    hist_alpha:    float = 0.70                   # PlotUncertaintyHistogram
    hist_edgecolor:str   = "white"
    hist_bins:     int   = 30
    kde_color:     str   = "k"
    kde_lw:        float = 1.8

    # ── convenience methods ───────────────────────────────────────────────

    def envelope_kwargs(self, color: str, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for ``ax.fill_between`` (uncertainty band)."""
        d: dict[str, Any] = dict(color=color, alpha=self.envelope_alpha)
        d.update(overrides)
        return d

    def line_kwargs(self, color: str, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for ``ax.plot`` (P50 central line)."""
        d: dict[str, Any] = dict(color=color, linewidth=self.p50_lw)
        d.update(overrides)
        return d

    def scatter_kwargs(self, color: str, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for ``ax.scatter``."""
        d: dict[str, Any] = dict(color=color, s=self.scatter_size, zorder=5)
        d.update(overrides)
        return d

    def ref_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for the reference-depth ``ax.axhline`` call."""
        d: dict[str, Any] = dict(
            color=self.ref_color,
            linestyle=self.ref_ls,
            linewidth=self.ref_lw,
        )
        d.update(overrides)
        return d

    def station_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Keyword arguments for station-tick ``ax.axvline`` calls."""
        d: dict[str, Any] = dict(
            color=self.station_color,
            linewidth=self.station_lw,
            alpha=self.station_alpha,
        )
        d.update(overrides)
        return d

    def grid_kwargs(self) -> dict[str, Any]:
        """Keyword arguments for ``ax.grid``."""
        return dict(alpha=self.grid_alpha, axis=self.grid_axis)

    def bar_kwargs(self, color: str, **overrides: Any) -> dict[str, Any]:
        """Kwargs for Dar-Zarrouk bar charts."""
        d: dict[str, Any] = dict(color=color, alpha=0.80, edgecolor="none")
        d.update(overrides)
        return d

    def hist_kwargs(self, color: str, **overrides: Any) -> dict[str, Any]:
        """Kwargs for ``ax.hist`` in uncertainty histograms."""
        d: dict[str, Any] = dict(
            color=color,
            alpha=self.hist_alpha,
            edgecolor=self.hist_edgecolor,
            bins=self.hist_bins,
            density=True,
        )
        d.update(overrides)
        return d

    def kde_kwargs(self, **overrides: Any) -> dict[str, Any]:
        """Kwargs for the KDE curve ``ax.plot`` call."""
        d: dict[str, Any] = dict(color=self.kde_color, linewidth=self.kde_lw)
        d.update(overrides)
        return d

    def copy(self, **kw: Any) -> HydroProfileStyle:
        """Return a shallow copy with optional field overrides."""
        new = copy.copy(self)
        for k, v in kw.items():
            setattr(new, k, v)
        return new


# ─────────────────────────────────────────────────────────────────────────────
# Composition dataclass — one complete interp preset bundle
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class InterpStyle:
    """Complete visual bundle for one hydrogeophysical interpretation preset.

    Parameters
    ----------
    section : HydroSectionStyle
        Controls all 2-D colour-section plots.
    profile : HydroProfileStyle
        Controls all 1-D profile plots.
    figsize_section : tuple
        Default figure size (width, height) for single-panel 2-D sections.
    figsize_profile : tuple
        Default figure size for 2-panel 1-D profile plots.
    figsize_uncertainty : tuple
        Default figure size for 2-panel uncertainty sections.
    """

    section:             HydroSectionStyle = field(
        default_factory=HydroSectionStyle
    )
    profile:             HydroProfileStyle = field(
        default_factory=HydroProfileStyle
    )
    figsize_section:       tuple[float, float] = (13.0, 5.0)
    figsize_profile:       tuple[float, float] = (13.0, 6.0)
    figsize_uncertainty:   tuple[float, float] = (13.0, 8.0)
    figsize_crossplot:     tuple[float, float] = (7.0,  5.5)
    figsize_depth_profile: tuple[float, float] = (4.5,  9.0)
    figsize_aquifer_char:  tuple[float, float] = (13.0, 9.0)
    figsize_tl_panel:      tuple[float, float] = (3.5,  4.5)

    def copy(self, **kw: Any) -> InterpStyle:
        """Return a deep copy with optional top-level field overrides."""
        new = copy.deepcopy(self)
        for k, v in kw.items():
            setattr(new, k, v)
        return new


# ─────────────────────────────────────────────────────────────────────────────
# Preset catalogue
# ─────────────────────────────────────────────────────────────────────────────

def _make_default() -> InterpStyle:
    return InterpStyle(
        section=HydroSectionStyle(
            cmap_K="viridis", cmap_Sw="RdYlBu", cmap_phi="YlOrRd",
            cmap_timelapse="RdBu_r", cmap_spread="hot_r", cmap_p50="viridis",
            cmap_crossplot="viridis",
            wt_color="deepskyblue", wt_lw=2.5, wt_ls="--",
            station_color="k", station_alpha=0.45, station_lw=0.40,
            cb_fraction=0.025, cb_pad=0.01, cb_fontsize=8,
            scatter_alpha=0.60, scatter_size=12.0,
            hs_alpha=0.12, model_curve_color="k", model_curve_ls="--", model_curve_lw=1.6,
            rho_curve_color="0.15", rho_fill_color="steelblue", rho_fill_alpha=0.12,
            zone_boundary_color="0.55", zone_boundary_ls="--", zone_label_fontsize=7,
        ),
        profile=HydroProfileStyle(
            color_wt="steelblue", color_T="seagreen",
            envelope_alpha=0.25, p50_lw=2.0, scatter_size=14.0,
            ref_color="tomato", ref_ls="--", ref_lw=1.4,
            grid_alpha=0.30, grid_axis="y",
            station_color="0.6", station_alpha=0.50, station_lw=0.35,
            color_TR="darkorange", color_S="mediumpurple",
            tr_threshold=1000.0, s_threshold_moderate=1.0, s_threshold_good=5.0,
            hist_alpha=0.70, hist_edgecolor="white", hist_bins=30,
            kde_color="k", kde_lw=1.8,
        ),
        figsize_section=(13.0, 5.0), figsize_profile=(13.0, 6.0),
        figsize_uncertainty=(13.0, 8.0), figsize_crossplot=(7.0, 5.5),
        figsize_depth_profile=(4.5, 9.0), figsize_aquifer_char=(13.0, 9.0),
        figsize_tl_panel=(3.5, 4.5),
    )


def _make_publication() -> InterpStyle:
    """Clean, journal-ready style — grayscale-compatible, compact figures."""
    return InterpStyle(
        section=HydroSectionStyle(
            cmap_K="plasma", cmap_Sw="RdBu", cmap_phi="cividis",
            cmap_timelapse="PiYG", cmap_spread="Greys", cmap_p50="plasma",
            cmap_crossplot="cividis",
            wt_color="black", wt_lw=1.2, wt_ls="--",
            station_color="k", station_alpha=0.30, station_lw=0.30,
            cb_fraction=0.020, cb_pad=0.01, cb_fontsize=7, nan_color="0.93",
            scatter_alpha=0.50, scatter_size=10.0,
            hs_alpha=0.10, model_curve_color="0.2", model_curve_ls="--", model_curve_lw=1.4,
            rho_curve_color="k", rho_fill_color="0.6", rho_fill_alpha=0.10,
            zone_boundary_color="0.45", zone_boundary_ls=":", zone_label_fontsize=6,
        ),
        profile=HydroProfileStyle(
            color_wt="k", color_T="0.40",
            envelope_alpha=0.15, p50_lw=1.8, scatter_size=10.0,
            ref_color="0.45", ref_ls=":", ref_lw=1.0,
            grid_alpha=0.20, grid_axis="y",
            station_color="0.70", station_alpha=0.40, station_lw=0.30,
            color_TR="0.30", color_S="0.60",
            tr_threshold=1000.0, s_threshold_moderate=1.0, s_threshold_good=5.0,
            hist_alpha=0.60, hist_edgecolor="white", hist_bins=30,
            kde_color="k", kde_lw=1.6,
        ),
        figsize_section=(8.5, 3.8), figsize_profile=(8.5, 5.0),
        figsize_uncertainty=(8.5, 6.5), figsize_crossplot=(5.5, 4.5),
        figsize_depth_profile=(3.8, 8.0), figsize_aquifer_char=(8.5, 7.5),
        figsize_tl_panel=(2.8, 3.8),
    )


def _make_dark() -> InterpStyle:
    """Dark-background theme — high contrast for presentations."""
    return InterpStyle(
        section=HydroSectionStyle(
            cmap_K="inferno", cmap_Sw="cool", cmap_phi="magma",
            cmap_timelapse="RdBu_r", cmap_spread="hot", cmap_p50="inferno",
            cmap_crossplot="plasma",
            wt_color="cyan", wt_lw=2.0, wt_ls="--",
            station_color="w", station_alpha=0.50, station_lw=0.40,
            cb_fraction=0.025, cb_pad=0.01, cb_fontsize=8, nan_color="0.18",
            scatter_alpha=0.65, scatter_size=14.0,
            hs_alpha=0.18, model_curve_color="w", model_curve_ls="--", model_curve_lw=1.8,
            rho_curve_color="w", rho_fill_color="deepskyblue", rho_fill_alpha=0.15,
            zone_boundary_color="0.70", zone_boundary_ls="--", zone_label_fontsize=7,
        ),
        profile=HydroProfileStyle(
            color_wt="cyan", color_T="#00ff88",
            envelope_alpha=0.20, p50_lw=2.0, scatter_size=14.0,
            ref_color="orangered", ref_ls="--", ref_lw=1.4,
            grid_alpha=0.20, grid_axis="y",
            station_color="w", station_alpha=0.40, station_lw=0.35,
            color_TR="orange", color_S="violet",
            tr_threshold=1000.0, s_threshold_moderate=1.0, s_threshold_good=5.0,
            hist_alpha=0.65, hist_edgecolor="0.3", hist_bins=30,
            kde_color="w", kde_lw=1.8,
        ),
        figsize_section=(13.0, 5.0), figsize_profile=(13.0, 6.0),
        figsize_uncertainty=(13.0, 8.0), figsize_crossplot=(7.0, 5.5),
        figsize_depth_profile=(4.5, 9.0), figsize_aquifer_char=(13.0, 9.0),
        figsize_tl_panel=(3.5, 4.5),
    )


def _make_accessible() -> InterpStyle:
    """Colorblind-safe palette (Paul Tol / IBM accessible colours)."""
    return InterpStyle(
        section=HydroSectionStyle(
            cmap_K="cividis", cmap_Sw="BrBG", cmap_phi="viridis",
            cmap_timelapse="PuOr", cmap_spread="Greys", cmap_p50="cividis",
            cmap_crossplot="cividis",
            wt_color="#0077bb", wt_lw=2.0, wt_ls="--",
            station_color="k", station_alpha=0.45, station_lw=0.40,
            cb_fraction=0.025, cb_pad=0.01, cb_fontsize=8,
            scatter_alpha=0.60, scatter_size=12.0,
            hs_alpha=0.12, model_curve_color="#ee3377", model_curve_ls="--", model_curve_lw=1.6,
            rho_curve_color="0.15", rho_fill_color="#0077bb", rho_fill_alpha=0.12,
            zone_boundary_color="0.55", zone_boundary_ls="--", zone_label_fontsize=7,
        ),
        profile=HydroProfileStyle(
            color_wt="#0077bb", color_T="#009988",
            envelope_alpha=0.25, p50_lw=2.0, scatter_size=14.0,
            ref_color="#ee3377", ref_ls="--", ref_lw=1.4,
            grid_alpha=0.30, grid_axis="y",
            station_color="0.6", station_alpha=0.50, station_lw=0.35,
            color_TR="#ee7733", color_S="#aa3377",
            tr_threshold=1000.0, s_threshold_moderate=1.0, s_threshold_good=5.0,
            hist_alpha=0.70, hist_edgecolor="white", hist_bins=30,
            kde_color="#333333", kde_lw=1.8,
        ),
        figsize_section=(13.0, 5.0), figsize_profile=(13.0, 6.0),
        figsize_uncertainty=(13.0, 8.0), figsize_crossplot=(7.0, 5.5),
        figsize_depth_profile=(4.5, 9.0), figsize_aquifer_char=(13.0, 9.0),
        figsize_tl_panel=(3.5, 4.5),
    )


# ─────────────────────────────────────────────────────────────────────────────
# Container singleton
# ─────────────────────────────────────────────────────────────────────────────

class PyCSAMTInterp:
    """Package-wide visual control container for hydro-geophysical plots.

    Holds one :class:`InterpStyle` per named preset.  The *default* preset
    is what plot classes use when no explicit ``style=`` argument is given.

    Examples
    --------
    >>> from pycsamt.api.interp import PYCSAMT_INTERP
    >>> PYCSAMT_INTERP.use("publication")
    >>> # now all subsequent interp plots use the publication style
    >>> with PYCSAMT_INTERP.context("dark"):
    ...     fig = PlotHydroSection(result).plot()
    >>> # reverts
    """

    def __init__(self) -> None:
        self._init_defaults()

    def _init_defaults(self) -> None:
        self.default:    InterpStyle = _make_default()
        self.publication:InterpStyle = _make_publication()
        self.dark:       InterpStyle = _make_dark()
        self.accessible: InterpStyle = _make_accessible()

    # ── preset access ──────────────────────────────────────────────────────

    def style_for(self, preset: str) -> InterpStyle:
        """Return the :class:`InterpStyle` for *preset*.

        Parameters
        ----------
        preset : str
            One of ``'default'``, ``'publication'``, ``'dark'``,
            ``'accessible'``.
        """
        key = str(preset).lower().strip()
        if key not in _PRESETS:
            raise ValueError(
                f"interp preset must be one of {sorted(_PRESETS)}, "
                f"got {preset!r}."
            )
        return getattr(self, key)

    # ── configure / use / reset ────────────────────────────────────────────

    def use(self, preset: str) -> None:
        """Copy a named preset into the *default* slot.

        After this call, all plot classes that read from ``PYCSAMT_INTERP``
        will use the chosen preset's styles.

        Parameters
        ----------
        preset : str
            Name of the preset to activate.
        """
        self.default = copy.deepcopy(self.style_for(preset))

    def configure(self, **kw: Any) -> None:
        """Configure the *default* style using ``sub__field`` dotted paths.

        Examples
        --------
        >>> PYCSAMT_INTERP.configure(
        ...     section__cmap_K="plasma",
        ...     section__wt_color="white",
        ...     profile__color_wt="royalblue",
        ...     figsize_section=(10, 4),
        ... )
        """
        for path, value in kw.items():
            parts = path.split("__")
            obj = self.default
            for part in parts[:-1]:
                obj = getattr(obj, part)
            setattr(obj, parts[-1], value)

    def reset(self) -> None:
        """Reset all presets to package defaults."""
        self._init_defaults()

    # ── context manager ────────────────────────────────────────────────────

    @contextmanager
    def context(
        self,
        preset: str | None = None,
        **kw: Any,
    ) -> Generator[PyCSAMTInterp, None, None]:
        """Temporarily override the *default* style, then restore it.

        Parameters
        ----------
        preset : str, optional
            Named preset to activate for the duration of the block.
        **kw
            Additional dotted-path overrides applied after the preset
            (same semantics as :meth:`configure`).

        Examples
        --------
        >>> with PYCSAMT_INTERP.context("publication",
        ...                             section__wt_lw=0.8):
        ...     fig = PlotHydroSection(result).plot()
        """
        snapshot = copy.deepcopy(self.default)
        try:
            if preset is not None:
                self.use(preset)
            if kw:
                self.configure(**kw)
            yield self
        finally:
            self.default = snapshot

    # ── introspection ──────────────────────────────────────────────────────

    def summary(self) -> str:
        """Return a human-readable summary of all presets."""
        lines = ["PyCSAMTInterp — hydro-geophysical plot styles"]
        for name in sorted(_PRESETS):
            sty = self.style_for(name)
            active = "  ← active" if name == "default" else ""
            lines.append(
                f"  {name:<12} "
                f"K={sty.section.cmap_K!r}  "
                f"Sw={sty.section.cmap_Sw!r}  "
                f"wt_color={sty.section.wt_color!r}  "
                f"fig_sec={sty.figsize_section}"
                f"{active}"
            )
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()

    # ── snapshot / restore helpers ─────────────────────────────────────────

    def _snapshot(self) -> dict[str, InterpStyle]:
        return {name: copy.deepcopy(self.style_for(name)) for name in _PRESETS}

    def _restore(self, snapshot: dict[str, InterpStyle]) -> None:
        for name, value in snapshot.items():
            setattr(self, name, value)


# ─────────────────────────────────────────────────────────────────────────────
# Module-level singleton and convenience functions
# ─────────────────────────────────────────────────────────────────────────────

PYCSAMT_INTERP: PyCSAMTInterp = PyCSAMTInterp()


def configure_interp(**kw: Any) -> None:
    """Configure :data:`PYCSAMT_INTERP` using dotted-path keys.

    Equivalent to ``PYCSAMT_INTERP.configure(**kw)``.
    """
    PYCSAMT_INTERP.configure(**kw)


def reset_interp() -> None:
    """Reset :data:`PYCSAMT_INTERP` to package defaults."""
    PYCSAMT_INTERP.reset()


def use_interp(preset: str) -> None:
    """Activate a named preset in :data:`PYCSAMT_INTERP`.

    Parameters
    ----------
    preset : str
        One of ``'default'``, ``'publication'``, ``'dark'``,
        ``'accessible'``.
    """
    PYCSAMT_INTERP.use(preset)


# ─────────────────────────────────────────────────────────────────────────────
# Style resolution helpers (used by interp/plot.py)
# ─────────────────────────────────────────────────────────────────────────────

def resolve_section_style(
    style: str | InterpStyle | HydroSectionStyle | None,
) -> HydroSectionStyle:
    """Resolve a *style* argument to a :class:`HydroSectionStyle`.

    * ``None``          → ``PYCSAMT_INTERP.default.section``
    * ``str``           → named preset's section style
    * :class:`InterpStyle` → ``.section`` attribute
    * :class:`HydroSectionStyle` → passed through unchanged
    """
    if style is None:
        return PYCSAMT_INTERP.default.section
    if isinstance(style, str):
        return PYCSAMT_INTERP.style_for(style).section
    if isinstance(style, InterpStyle):
        return style.section
    if isinstance(style, HydroSectionStyle):
        return style
    raise TypeError(
        f"style must be None, str, InterpStyle, or HydroSectionStyle, "
        f"got {type(style).__name__}."
    )


def resolve_profile_style(
    style: str | InterpStyle | HydroProfileStyle | None,
) -> HydroProfileStyle:
    """Resolve a *style* argument to a :class:`HydroProfileStyle`.

    * ``None``          → ``PYCSAMT_INTERP.default.profile``
    * ``str``           → named preset's profile style
    * :class:`InterpStyle` → ``.profile`` attribute
    * :class:`HydroProfileStyle` → passed through unchanged
    """
    if style is None:
        return PYCSAMT_INTERP.default.profile
    if isinstance(style, str):
        return PYCSAMT_INTERP.style_for(style).profile
    if isinstance(style, InterpStyle):
        return style.profile
    if isinstance(style, HydroProfileStyle):
        return style
    raise TypeError(
        f"style must be None, str, InterpStyle, or HydroProfileStyle, "
        f"got {type(style).__name__}."
    )


def resolve_figsize(
    user_figsize: tuple[float, float] | None,
    style: str | InterpStyle | None,
    kind: str,
) -> tuple[float, float]:
    """Return the effective figure size.

    *user_figsize* takes priority; falls back to the active preset's
    size for *kind* (``'section'``, ``'profile'``, or ``'uncertainty'``).

    Leaf styles (:class:`HydroSectionStyle`, :class:`HydroProfileStyle`)
    carry no figsize; the active default bundle is used instead.
    """
    if user_figsize is not None:
        return user_figsize
    if isinstance(style, str):
        bundle = PYCSAMT_INTERP.style_for(style)
    elif isinstance(style, InterpStyle):
        bundle = style
    else:
        # None, HydroSectionStyle, HydroProfileStyle → active default
        bundle = PYCSAMT_INTERP.default
    return getattr(bundle, f"figsize_{kind}")
