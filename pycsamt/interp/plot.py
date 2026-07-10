# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""plot — visualization for pycsamt.interp results.

Three plot classes cover the main interpretation deliverables:

* :class:`PlotStratigraphicLog` — single-station pseudo-stratigraphic
  log in the style of Fig. 5d / Fig. 7 of Kouadio et al. (2022).
* :class:`PlotFenceDiagram` — multi-station panel with all logs
  arranged along the profile.
* :class:`PlotCalibratedModel` — side-by-side CRM vs NM with the
  misfit G (%) map overlaid.

All classes follow the same pattern::

    fig = PlotStratigraphicLog(log).plot()
    fig.savefig("S17_log.png", dpi=200)
"""
from __future__ import annotations

from collections.abc import Sequence
from typing import Optional, Union

import numpy as np

from ..api.interp import (
    HydroProfileStyle,
    HydroSectionStyle,
    InterpStyle,
    resolve_figsize,
    resolve_profile_style,
    resolve_section_style,
)
from ._base import ResistivityModel
from .lithology import StratigraphicLog

# Type alias accepted by the style= parameter of every hydro plot class.
_StyleArg = Optional[Union[str, InterpStyle, HydroSectionStyle, HydroProfileStyle]]

__all__ = [
    "PlotStratigraphicLog",
    "PlotFenceDiagram",
    "PlotCalibratedModel",
    # hydro-geophysics plots
    "PlotHydroSection",
    "PlotWaterTableProfile",
    "PlotTimeLapseSection",
    # uncertainty plots
    "PlotUncertaintySection",
    "PlotUncertaintyProfile",
    # new diagnostic / novel plots
    "PlotPetrophysicalCrossPlot",
    "PlotAquiferCharacterization",
    "PlotMultiTimeLapseGrid",
    "PlotResistivityDepthProfile",
    "PlotUncertaintyHistogram",
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _require_mpl():
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        return matplotlib, plt
    except ImportError as exc:
        raise ImportError(
            "matplotlib is required for pycsamt.interp.plot"
        ) from exc


def _cmap_with_bad(name: str, nan_color: str):
    """Return a copy of colourmap *name* with NaN cells rendered as *nan_color*."""
    import matplotlib.pyplot as plt
    cmap = plt.get_cmap(name).copy()
    cmap.set_bad(color=nan_color)
    return cmap


def _hatch_for(lithology: str) -> str:
    """Return a matplotlib hatch string keyed on lithology name."""
    name = lithology.lower()
    if "granite" in name or "igneous" in name or "basement" in name:
        return "+"
    if "fractured" in name or "fault" in name:
        return "x"
    if "aquifer" in name or "water" in name:
        return "o"
    if "clay" in name or "shale" in name:
        return "---"
    if "sand" in name or "alluvium" in name:
        return "..."
    if "basalt" in name or "gabbro" in name:
        return "///"
    if "limestone" in name or "dolomite" in name or "marble" in name:
        return r"\\\\"
    if "schist" in name or "gneiss" in name or "quartzite" in name:
        return "|||"
    if "ore" in name or "sulfide" in name or "graphite" in name:
        return "**"
    return ""


# ---------------------------------------------------------------------------
# PlotStratigraphicLog
# ---------------------------------------------------------------------------

class PlotStratigraphicLog:
    """Single-station pseudo-stratigraphic log.

    Reproduces the two-panel layout of Fig. 5d / Fig. 7 in
    Kouadio et al. (2022):

    * **Left panel** — colour / hatch blocks for each geological layer
      with lithology annotations and thickness values.
    * **Right panel** — log₁₀(ρ) depth curve overlaid on the same
      depth axis.

    Parameters
    ----------
    log : StratigraphicLog
    figsize : tuple
    depth_unit : str
        Label for the depth axis (default ``'m'``).
    title : str, optional
    annotation_kws : dict, optional
        Extra keyword arguments passed to ``ax.annotate``.
    """

    def __init__(
        self,
        log: StratigraphicLog,
        *,
        figsize: tuple[float, float] = (8, 10),
        depth_unit: str = "m",
        title: str | None = None,
        annotation_kws: dict | None = None,
    ) -> None:
        self.log = log
        self.figsize = figsize
        self.depth_unit = depth_unit
        self.title = title or f"Pseudo-Stratigraphic Log — {log.station_name}"
        self.annotation_kws = annotation_kws or {"fontsize": 8}

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        log = self.log
        fig, (ax_log, ax_rho) = plt.subplots(
            1, 2, figsize=self.figsize,
            sharey=True,
            gridspec_kw={"width_ratios": [1, 1.4]},
        )

        # ── left: stratigraphic column ──────────────────────────────
        for layer in log.layers:
            ax_log.barh(
                y=(layer.top + layer.bottom) / 2,
                width=1.0,
                height=layer.thickness,
                color=layer.color,
                alpha=0.75,
                hatch=_hatch_for(layer.lithology),
                edgecolor="0.3",
                linewidth=0.5,
            )
            mid = (layer.top + layer.bottom) / 2
            label = f"{layer.lithology}\n({layer.thickness:.1f} {self.depth_unit})"
            ax_log.annotate(
                label,
                xy=(0.5, mid),
                xycoords=("axes fraction", "data"),
                ha="center", va="center",
                **self.annotation_kws,
            )

        ax_log.set_xlim(0, 1)
        ax_log.set_ylim(log.z_centers[-1] + 10, -5)
        ax_log.set_xticks([])
        ax_log.set_ylabel(f"Depth ({self.depth_unit})")
        ax_log.set_title("Lithology", fontsize=9)
        ax_log.invert_yaxis()

        # ── right: resistivity curve ─────────────────────────────────
        valid = ~np.isnan(log.rho_log10)
        ax_rho.plot(log.rho_log10[valid], log.z_centers[valid],
                    color="0.2", linewidth=1.4, zorder=3)
        ax_rho.fill_betweenx(
            log.z_centers[valid], log.rho_log10[valid],
            alpha=0.12, color="steelblue"
        )

        for layer in log.layers:
            ax_rho.axhline(layer.top, color="0.55", linewidth=0.6, linestyle="--")

        ax_rho.set_xlabel(r"$\log_{10}(\rho\ /\ \Omega\mathrm{m})$")
        ax_rho.set_title("Resistivity", fontsize=9)
        ax_rho.grid(axis="x", alpha=0.3)

        fig.suptitle(self.title, fontweight="bold", y=1.01)
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotFenceDiagram
# ---------------------------------------------------------------------------

class PlotFenceDiagram:
    """Multi-station fence diagram of pseudo-stratigraphic logs.

    Plots every log as a vertical panel side by side, sharing the depth
    axis, so the lateral geological evolution along the profile is
    immediately visible.

    Parameters
    ----------
    logs : list of StratigraphicLog
        Ordered list of station logs (West → East, or South → North).
    figsize : tuple, optional
        Defaults to ``(2 * n_logs, 10)``.
    title : str, optional
    max_depth : float, optional
        Truncate display at this depth (metres).
    """

    def __init__(
        self,
        logs: Sequence[StratigraphicLog],
        *,
        figsize: tuple[float, float] | None = None,
        title: str = "Fence Diagram",
        max_depth: float | None = None,
    ) -> None:
        self.logs = list(logs)
        self.figsize = figsize or (2 * len(logs), 10)
        self.title = title
        self.max_depth = max_depth

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        n = len(self.logs)
        if n == 0:
            raise ValueError("No logs provided.")

        fig, axes = plt.subplots(
            1, n, figsize=self.figsize,
            sharey=True,
            gridspec_kw={"wspace": 0.05},
        )
        if n == 1:
            axes = [axes]

        z_max = self.max_depth or max(
            log.z_centers[-1] for log in self.logs
        )

        for ax, log in zip(axes, self.logs):
            for layer in log.layers:
                if layer.top > z_max:
                    continue
                bottom = min(layer.bottom, z_max)
                ax.barh(
                    y=(layer.top + bottom) / 2,
                    width=1.0,
                    height=bottom - layer.top,
                    color=layer.color,
                    alpha=0.80,
                    hatch=_hatch_for(layer.lithology),
                    edgecolor="0.3",
                    linewidth=0.4,
                )

            ax.set_xlim(0, 1)
            ax.set_ylim(z_max + 10, -5)
            ax.set_xticks([])
            ax.set_title(log.station_name, fontsize=7, pad=3)
            ax.invert_yaxis()

        axes[0].set_ylabel("Depth (m)")
        fig.suptitle(self.title, fontweight="bold")
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotCalibratedModel
# ---------------------------------------------------------------------------

class PlotCalibratedModel:
    """Compare CRM vs NM and display the G (%) misfit map.

    Three sub-plots stacked vertically:

    1. CRM — original inversion result (log₁₀ρ colour image)
    2. NM — calibrated New Model (same colour scale)
    3. Misfit G (%) — diverging colour scale highlighting where the
       model was corrected the most

    Parameters
    ----------
    crm : ResistivityModel
        Original CRM.
    nm : ResistivityModel
        Calibrated NM from :meth:`~ModelCalibrator.calibrated_model`.
    misfit_map : ndarray (n_z, n_x), optional
        G (%) array from :meth:`~ModelCalibrator.misfit_map`.
        If ``None``, computed from the difference between *nm* and *crm*.
    figsize : tuple
    cmap_rho : str
        Matplotlib colourmap for the resistivity panels.
    vmin_rho, vmax_rho : float
        Colour-scale limits for log₁₀(ρ).
    title : str, optional
    """

    def __init__(
        self,
        crm: ResistivityModel,
        nm: ResistivityModel,
        misfit_map: np.ndarray | None = None,
        *,
        figsize: tuple[float, float] = (12, 10),
        cmap_rho: str = "jet",
        vmin_rho: float = 1.0,
        vmax_rho: float = 5.0,
        title: str | None = None,
    ) -> None:
        self.crm = crm
        self.nm = nm
        if misfit_map is not None:
            self._misfit = np.asarray(misfit_map)
        else:
            diff = nm.rho_2d - crm.rho_2d
            with np.errstate(divide="ignore", invalid="ignore"):
                self._misfit = 100.0 * np.abs(diff) / np.maximum(np.abs(crm.rho_2d), 1e-12)
        self.figsize = figsize
        self.cmap_rho = cmap_rho
        self.vmin_rho = vmin_rho
        self.vmax_rho = vmax_rho
        self.title = title or "CRM vs Calibrated NM"

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        fig, axes = plt.subplots(3, 1, figsize=self.figsize, sharex=True)
        ax_crm, ax_nm, ax_g = axes

        x = self.crm.x_centers
        z = self.crm.z_centers

        extent = [x[0], x[-1], z[-1], z[0]]

        def _rho_im(ax, data, label):
            im = ax.imshow(
                data, aspect="auto", extent=extent,
                cmap=self.cmap_rho,
                vmin=self.vmin_rho, vmax=self.vmax_rho,
                origin="upper",
            )
            ax.set_ylabel("Depth (m)")
            ax.set_title(label, fontsize=9)
            plt.colorbar(im, ax=ax, label=r"$\log_{10}(\rho)$",
                         fraction=0.03, pad=0.01)

        _rho_im(ax_crm, self.crm.rho_2d, "CRM — Inversion result")
        _rho_im(ax_nm,  self.nm.rho_2d,  "NM — Calibrated model")

        # Misfit panel
        g_clip = np.clip(self._misfit, 0, 10)
        im_g = ax_g.imshow(
            g_clip, aspect="auto", extent=extent,
            cmap="RdYlBu_r", vmin=0, vmax=10, origin="upper",
        )
        ax_g.set_ylabel("Depth (m)")
        ax_g.set_xlabel("Profile distance (m)")
        ax_g.set_title("Misfit G (%)", fontsize=9)
        plt.colorbar(im_g, ax=ax_g, label="G (%)", fraction=0.03, pad=0.01)

        # Station markers on all panels
        if len(self.crm.station_x):
            for ax in axes:
                for sx in self.crm.station_x:
                    ax.axvline(sx, color="k", linewidth=0.4, alpha=0.5)

        fig.suptitle(self.title, fontweight="bold")
        fig.tight_layout()
        return fig


# ---------------------------------------------------------------------------
# PlotHydroSection
# ---------------------------------------------------------------------------

class PlotHydroSection:
    """2-D hydrogeological section from an :class:`~pycsamt.interp.hydromodel.EMHydroResult`.

    Renders a colour-image cross-section of one quantitative hydro map
    (hydraulic conductivity K, water saturation Sw, or porosity φ) with
    optional overlays:

    * **Water-table line** — dashed blue line at the estimated WT depth.
    * **Station markers** — thin vertical tick-marks at each profile station.

    Parameters
    ----------
    result : EMHydroResult
        Quantitative hydro output from :class:`~pycsamt.interp.hydromodel.EMHydroModel`.
    quantity : str
        Which map to display:
        ``'K'`` — hydraulic conductivity (log₁₀ scale),
        ``'saturation'`` — water saturation Sw,
        ``'porosity'`` — effective porosity φ.
    cmap : str
        Matplotlib colourmap.  Defaults: K → ``'viridis'``,
        Sw → ``'RdYlBu'``, φ → ``'YlOrRd'``.
    vmin, vmax : float, optional
        Colour-scale limits.  Auto-derived from the data if ``None``.
    show_water_table : bool
        Overlay the water-table depth profile (default ``True``).
    figsize : tuple
    title : str, optional
    depth_min : float, optional
        Start display at this depth (m).  Use to zoom past near-surface
        artefacts and push a shallow water-table line into view.
    depth_max : float, optional
        Truncate display at this depth (m).

    Examples
    --------
    >>> fig = PlotHydroSection(result, quantity='K').plot()
    >>> fig = PlotHydroSection(result, quantity='saturation',
    ...                        cmap='Blues', vmin=0, vmax=1).plot()
    """

    _LABELS: dict = {
        "K":          r"$\log_{10}(K\ /\ \mathrm{m\,s^{-1}})$",
        "saturation": r"Water saturation $S_w$",
        "porosity":   r"Porosity $\phi$",
    }

    def __init__(
        self,
        result,
        quantity: str = "K",
        *,
        style: _StyleArg = None,
        cmap: str | None = None,
        vmin: float | None = None,
        vmax: float | None = None,
        show_water_table: bool = True,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
        depth_min: float | None = None,
        depth_max: float | None = None,
    ) -> None:
        if quantity not in self._LABELS:
            raise ValueError(
                f"quantity must be one of {list(self._LABELS)}, got {quantity!r}."
            )
        self.result           = result
        self.quantity         = quantity
        self.style            = style
        self.cmap             = cmap
        self.vmin             = vmin
        self.vmax             = vmax
        self.show_water_table = show_water_table
        self.figsize          = figsize
        self.title            = title
        self.depth_min        = depth_min
        self.depth_max        = depth_max

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty    = resolve_section_style(self.style)
        result = self.result
        model  = result.resistivity_model

        if self.quantity == "K":
            raw  = result.hydraulic_K
            data = np.log10(np.where(raw > 0, raw, np.nan))
        elif self.quantity == "saturation":
            data = result.saturation
        else:
            data = result.porosity

        z, x   = model.z_centers, model.x_centers
        z0     = self.depth_min if self.depth_min is not None else float(z[0])
        dz     = self.depth_max if self.depth_max is not None else float(z[-1])
        z_mask = (z >= z0) & (z <= dz)
        data_plot, z_plot = data[z_mask, :], z[z_mask]

        vmin   = self.vmin if self.vmin is not None else float(np.nanpercentile(data_plot,  2))
        vmax   = self.vmax if self.vmax is not None else float(np.nanpercentile(data_plot, 98))
        cmap   = _cmap_with_bad(
            self.cmap if self.cmap is not None else sty.cmap_for(self.quantity),
            sty.nan_color,
        )
        fsz    = resolve_figsize(self.figsize, self.style, "section")
        extent = [x[0], x[-1], z_plot[-1], z_plot[0]]

        fig, ax = plt.subplots(figsize=fsz, layout="constrained")
        im = ax.imshow(
            data_plot, aspect="auto", extent=extent,
            cmap=cmap, vmin=vmin, vmax=vmax, origin="upper",
        )
        cb = fig.colorbar(im, ax=ax, label=self._LABELS[self.quantity],
                          **sty.cb_kwargs())
        cb.ax.tick_params(labelsize=sty.cb_fontsize)

        if self.show_water_table:
            wt_plot = np.where(np.isfinite(result.water_table),
                               result.water_table, np.nan)
            ax.plot(x, wt_plot, **sty.wt_kwargs(), label="Water table")
            ax.legend(fontsize=8, loc="lower right")
        # enforce depth clip on the y-axis (imshow extent sets the image bounds
        # but does not clip when depth_min > z[0], so we set limits explicitly)
        ax.set_ylim(dz, z0)

        if len(model.station_x):
            for sx in model.station_x:
                ax.axvline(sx, **sty.station_kwargs())

        ax.set_xlabel("Profile distance (m)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(
            self.title or f"Hydro section — {self.quantity}  [{result.method_tag}]",
            fontsize=10,
        )
        return fig


# ---------------------------------------------------------------------------
# PlotWaterTableProfile
# ---------------------------------------------------------------------------

class PlotWaterTableProfile:
    """Profile plot: water-table depth and transmissivity along the section.

    Two stacked panels share the same x-axis (profile distance):

    * **Top** — water-table depth (m) as a stem/bar plot.  Shallower is
      better; the y-axis is inverted so deep values plot downward.
    * **Bottom** — transmissivity T (m²/s) on a log₁₀ scale.

    Parameters
    ----------
    result : EMHydroResult
    figsize : tuple
    color_wt : str
        Colour for the water-table bars (default ``'steelblue'``).
    color_T : str
        Colour for the transmissivity bars (default ``'seagreen'``).
    reference_depth : float, optional
        Draw a horizontal dashed line at this depth on the WT panel
        (e.g. a known piezometric level).
    title : str, optional

    Examples
    --------
    >>> fig = PlotWaterTableProfile(result, reference_depth=20.0).plot()
    """

    def __init__(
        self,
        result,
        *,
        style: _StyleArg = None,
        color_wt: str | None = None,
        color_T:  str | None = None,
        reference_depth: float | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        self.result          = result
        self.style           = style
        self.color_wt        = color_wt
        self.color_T         = color_T
        self.reference_depth = reference_depth
        self.figsize         = figsize
        self.title           = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty    = resolve_profile_style(self.style)
        result = self.result
        model  = result.resistivity_model
        x      = model.x_centers
        cwt    = self.color_wt if self.color_wt is not None else sty.color_wt
        cT     = self.color_T  if self.color_T  is not None else sty.color_T
        fsz    = resolve_figsize(self.figsize, self.style, "profile")

        wt    = result.water_table
        T_log = np.log10(np.where(result.transmissivity > 0,
                                  result.transmissivity, np.nan))

        fig, (ax_wt, ax_T) = plt.subplots(
            2, 1, figsize=fsz, sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": [1, 1]},
        )

        # ── water table panel ─────────────────────────────────────────────
        valid_wt = np.isfinite(wt)
        ax_wt.bar(x[valid_wt], wt[valid_wt], width=np.diff(x).mean() * 0.7,
                  color=cwt, alpha=0.75, edgecolor="none")
        ax_wt.scatter(x[valid_wt], wt[valid_wt],
                      **sty.scatter_kwargs(cwt, s=18))

        if self.reference_depth is not None:
            ax_wt.axhline(self.reference_depth,
                          **sty.ref_kwargs(),
                          label=f"Ref. {self.reference_depth} m")
            ax_wt.legend(fontsize=8)

        ax_wt.set_ylabel("Water-table depth (m)")
        ax_wt.invert_yaxis()
        ax_wt.grid(**sty.grid_kwargs())
        ax_wt.set_title(
            self.title or f"Water table & transmissivity [{result.method_tag}]",
            fontsize=10,
        )

        # ── transmissivity panel ──────────────────────────────────────────
        valid_T = np.isfinite(T_log)
        ax_T.bar(x[valid_T], T_log[valid_T], width=np.diff(x).mean() * 0.7,
                 color=cT, alpha=0.75, edgecolor="none",
                 bottom=T_log[valid_T].min() - 0.5 if valid_T.any() else 0)
        ax_T.scatter(x[valid_T], T_log[valid_T],
                     **sty.scatter_kwargs(cT, s=18))
        ax_T.set_ylabel(r"$\log_{10}(T\ /\ \mathrm{m^2\,s^{-1}})$")
        ax_T.set_xlabel("Profile distance (m)")
        ax_T.grid(**sty.grid_kwargs())

        if len(model.station_x):
            for sx in model.station_x:
                ax_wt.axvline(sx, **sty.station_kwargs())
                ax_T.axvline(sx,  **sty.station_kwargs())

        return fig


# ---------------------------------------------------------------------------
# PlotTimeLapseSection
# ---------------------------------------------------------------------------

class PlotTimeLapseSection:
    """Difference section for time-lapse EM monitoring.

    Shows Δlog₁₀(ρ) or Δ*Sw* between a selected survey and the baseline
    as a diverging colour image, allowing rapid visual identification of
    zones that became more conductive (wetting) or more resistive (drying).

    Parameters
    ----------
    timelapse : TimeLapseEM
        Time-lapse container holding all surveys.
    quantity : str
        ``'rho'`` — raw resistivity change Δlog₁₀ρ (default).
        ``'saturation'`` — saturation change ΔSw (requires *petro* and *rho_w*).
    survey_idx : int
        Index of the comparison survey in ``timelapse.surveys`` relative to the
        baseline (default 0 → first non-baseline survey).
    baseline_idx : int
        Index of the baseline survey (default 0).
    petro : ArchieModel, optional
        Required when ``quantity='saturation'``.
    rho_w : float
        Pore-water resistivity for saturation inversion (default 20 Ω·m).
    phi : float or ndarray
        Porosity for Archie inversion (default 0.25).
    cmap : str
        Diverging colourmap (default ``'RdBu_r'``).
        Positive (blue) → resistivity increase / drying.
        Negative (red)  → resistivity decrease / wetting.
    vmax : float, optional
        Symmetric colour-scale limit.  Auto-derived from the 98th percentile
        of ``|data|`` if ``None``.
    show_water_table : bool
        Overlay the water-table profile from the *comparison* survey.
    figsize : tuple
    title : str, optional
    depth_max : float, optional

    Examples
    --------
    >>> from pycsamt.interp.timelapse import TimeLapseEM
    >>> tl  = TimeLapseEM([model_dry, model_wet], labels=['dry', 'wet'])
    >>> fig = PlotTimeLapseSection(tl, quantity='rho').plot()
    >>> fig = PlotTimeLapseSection(
    ...     tl, quantity='saturation',
    ...     petro=ArchieModel(), rho_w=20.0
    ... ).plot()
    """

    def __init__(
        self,
        timelapse,
        quantity: str = "rho",
        *,
        style: _StyleArg = None,
        survey_idx: int = 0,
        baseline_idx: int = 0,
        petro=None,
        rho_w: float = 20.0,
        phi: float | np.ndarray = 0.25,
        cmap: str | None = None,
        vmax: float | None = None,
        show_water_table: bool = True,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
        depth_max: float | None = None,
    ) -> None:
        if quantity not in ("rho", "saturation"):
            raise ValueError("quantity must be 'rho' or 'saturation'.")
        if quantity == "saturation" and petro is None:
            raise ValueError("petro (ArchieModel) is required for quantity='saturation'.")
        self.timelapse        = timelapse
        self.quantity         = quantity
        self.style            = style
        self.survey_idx       = survey_idx
        self.baseline_idx     = baseline_idx
        self.petro            = petro
        self.rho_w            = rho_w
        self.phi              = phi
        self.cmap             = cmap
        self.vmax             = vmax
        self.show_water_table = show_water_table
        self.figsize          = figsize
        self.title            = title
        self.depth_max        = depth_max

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty = resolve_section_style(self.style)
        tl  = self.timelapse

        if self.quantity == "rho":
            deltas   = tl.resistivity_change(baseline_idx=self.baseline_idx)
            cb_label = r"$\Delta\log_{10}(\rho\ /\ \Omega\mathrm{m})$"
            cmap     = self.cmap if self.cmap is not None else sty.cmap_timelapse
        else:
            deltas   = tl.saturation_change(
                self.petro, phi=self.phi, rho_w=self.rho_w,
                baseline_idx=self.baseline_idx,
            )
            cb_label = r"$\Delta S_w$"
            cmap     = self.cmap if self.cmap is not None else sty.cmap_timelapse

        if self.survey_idx >= len(deltas):
            raise IndexError(
                f"survey_idx={self.survey_idx} out of range; "
                f"only {len(deltas)} comparison surveys available."
            )
        delta      = deltas[self.survey_idx]
        others     = [i for i in range(tl.n_surveys) if i != self.baseline_idx]
        comp_label = tl.labels[others[self.survey_idx]]
        base_label = tl.labels[self.baseline_idx]

        ref_model  = tl.surveys[self.baseline_idx]
        x, z       = ref_model.x_centers, ref_model.z_centers
        dz         = self.depth_max or float(z[-1])
        z_mask     = z <= dz
        data_plot  = delta[z_mask, :]

        abs_max = float(np.nanpercentile(np.abs(data_plot), 98))
        vmax    = self.vmax if self.vmax is not None else max(abs_max, 1e-6)
        fsz     = resolve_figsize(self.figsize, self.style, "section")
        extent  = [x[0], x[-1], z[z_mask][-1], z[z_mask][0]]
        cmap    = _cmap_with_bad(cmap, sty.nan_color)

        fig, ax = plt.subplots(figsize=fsz, layout="constrained")
        im = ax.imshow(data_plot, aspect="auto", extent=extent,
                       cmap=cmap, vmin=-vmax, vmax=vmax, origin="upper")
        cb = fig.colorbar(im, ax=ax, label=cb_label, **sty.cb_kwargs())
        cb.ax.tick_params(labelsize=sty.cb_fontsize)

        if self.show_water_table:
            comp_survey = tl.surveys[others[self.survey_idx]]
            from .petrophysics import ArchieModel as _Archie
            from .petrophysics import water_table_from_profile
            _archie = self.petro if self.petro is not None else _Archie()
            wt = np.array([
                (lambda d: d if d is not None else np.nan)(
                    water_table_from_profile(
                        comp_survey.rho_2d[:, ix], comp_survey.z_centers,
                        _archie, rho_w=self.rho_w,
                    )
                )
                for ix in range(comp_survey.n_x)
            ])
            ax.plot(x, wt, **sty.wt_kwargs(), label=f"WT ({comp_label})")
            ax.legend(fontsize=8, loc="lower right")

        if len(ref_model.station_x):
            for sx in ref_model.station_x:
                ax.axvline(sx, **sty.station_kwargs())

        ax.set_xlabel("Profile distance (m)")
        ax.set_ylabel("Depth (m)")
        ax.set_title(
            self.title or (
                f"Time-lapse section — {self.quantity}:  "
                f"{comp_label} − {base_label}"
            ),
            fontsize=10,
        )
        return fig


# ---------------------------------------------------------------------------
# PlotUncertaintySection
# ---------------------------------------------------------------------------

class PlotUncertaintySection:
    """Two-panel section showing the P50 estimate and the uncertainty spread.

    Panel layout:

    * **Top** — P50 (median) of the selected quantity as a colour image.
    * **Bottom** — Uncertainty spread: either the coefficient of variation
      (CV = std/mean, for K) or the P90–P10 range (for Sw, porosity).

    Parameters
    ----------
    unc : UncertaintyResult
    quantity : str
        ``'K'`` (default) — hydraulic conductivity (log₁₀ scale for P50,
        CV for spread).
        ``'saturation'`` — Sw (P50 and P90–P10 range).
        ``'porosity'`` — φ (P50 and P90–P10 range).
    cmap_p50 : str
        Colourmap for the P50 panel (defaults mirror :class:`PlotHydroSection`).
    cmap_spread : str
        Colourmap for the spread panel (default ``'hot_r'`` — dark = high
        uncertainty).
    vmin_p50, vmax_p50 : float, optional
        Colour limits for P50 panel.
    vmax_spread : float, optional
        Upper colour limit for the spread panel.  Auto if ``None``.
    show_water_table : bool
        Overlay median water-table on both panels (default ``True``).
    figsize : tuple
    title : str, optional
    depth_min : float, optional
        Start display at this depth (m) — mirrors :class:`PlotHydroSection`.
    depth_max : float, optional

    Examples
    --------
    >>> fig = PlotUncertaintySection(unc_result, quantity='K').plot()
    >>> fig = PlotUncertaintySection(unc_result, quantity='saturation',
    ...                              depth_max=200.0).plot()
    """

    _META = {
        "K":          {"p50_cmap": "viridis", "spread_label": "CV  (std / mean K)",
                       "p50_label": r"$\log_{10}(K_{P50}\ /\ \mathrm{m\,s^{-1}})$"},
        "saturation": {"p50_cmap": "RdYlBu",  "spread_label": r"$S_w$ P90 − P10",
                       "p50_label": r"$S_{w,P50}$"},
        "porosity":   {"p50_cmap": "YlOrRd",  "spread_label": r"$\phi$ P90 − P10",
                       "p50_label": r"$\phi_{P50}$"},
    }

    def __init__(
        self,
        unc,
        quantity: str = "K",
        *,
        style: _StyleArg = None,
        cmap_p50:    str | None = None,
        cmap_spread: str | None = None,
        vmin_p50:    float | None = None,
        vmax_p50:    float | None = None,
        vmax_spread: float | None = None,
        show_water_table: bool = True,
        figsize: tuple[float, float] | None = None,
        title:   str | None = None,
        depth_min: float | None = None,
        depth_max: float | None = None,
    ) -> None:
        if quantity not in self._META:
            raise ValueError(f"quantity must be one of {list(self._META)}.")
        self.unc              = unc
        self.quantity         = quantity
        self.style            = style
        self.cmap_p50         = cmap_p50
        self.cmap_spread      = cmap_spread
        self.vmin_p50         = vmin_p50
        self.vmax_p50         = vmax_p50
        self.vmax_spread      = vmax_spread
        self.show_water_table = show_water_table
        self.figsize          = figsize
        self.title            = title
        self.depth_min        = depth_min
        self.depth_max        = depth_max

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty   = resolve_section_style(self.style)
        unc   = self.unc
        model = unc.resistivity_model
        meta  = self._META[self.quantity]
        x, z  = model.x_centers, model.z_centers
        z0    = self.depth_min if self.depth_min is not None else float(z[0])
        dz    = self.depth_max if self.depth_max is not None else float(z[-1])
        zm    = (z >= z0) & (z <= dz)

        cmap_p50    = _cmap_with_bad(self.cmap_p50    or sty.cmap_for(self.quantity), sty.nan_color)
        cmap_spread = _cmap_with_bad(self.cmap_spread or sty.cmap_spread, sty.nan_color)
        fsz         = resolve_figsize(self.figsize, self.style, "uncertainty")

        if self.quantity == "K":
            p50_raw  = unc.p50_K[zm, :]
            p50_plot = np.log10(np.where(p50_raw > 0, p50_raw, np.nan))
            spread   = unc.cv_K[zm, :]
        elif self.quantity == "saturation":
            p50_plot = unc.mean_Sw[zm, :]
            spread   = unc.p90_Sw[zm, :] - unc.p10_Sw[zm, :]
        else:
            p50_plot = unc.mean_phi[zm, :]
            spread   = unc.std_phi[zm, :] * 2.0

        vmin_p  = self.vmin_p50  if self.vmin_p50  is not None else float(np.nanpercentile(p50_plot,  2))
        vmax_p  = self.vmax_p50  if self.vmax_p50  is not None else float(np.nanpercentile(p50_plot, 98))
        vmax_sp = self.vmax_spread if self.vmax_spread is not None else float(np.nanpercentile(spread, 98))
        extent  = [x[0], x[-1], z[zm][-1], z[zm][0]]

        fig, (ax_p50, ax_sp) = plt.subplots(
            2, 1, figsize=fsz, sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": [1, 1]},
        )

        im1 = ax_p50.imshow(p50_plot, aspect="auto", extent=extent,
                            cmap=cmap_p50, vmin=vmin_p, vmax=vmax_p, origin="upper")
        cb1 = fig.colorbar(im1, ax=ax_p50, label=meta["p50_label"], **sty.cb_kwargs())
        cb1.ax.tick_params(labelsize=sty.cb_fontsize)
        ax_p50.set_ylabel("Depth (m)")
        ax_p50.set_title("P50 estimate", fontsize=9)

        im2 = ax_sp.imshow(spread, aspect="auto", extent=extent,
                           cmap=cmap_spread, vmin=0, vmax=vmax_sp, origin="upper")
        cb2 = fig.colorbar(im2, ax=ax_sp, label=meta["spread_label"], **sty.cb_kwargs())
        cb2.ax.tick_params(labelsize=sty.cb_fontsize)
        ax_sp.set_ylabel("Depth (m)")
        ax_sp.set_xlabel("Profile distance (m)")
        ax_sp.set_title("Uncertainty spread", fontsize=9)

        if self.show_water_table:
            wt = np.where(np.isfinite(unc.p50_wt), unc.p50_wt, np.nan)
            for ax in (ax_p50, ax_sp):
                ax.plot(x, wt, **sty.wt_kwargs(), label="WT P50")
            ax_p50.legend(fontsize=8, loc="lower right")

        for ax in (ax_p50, ax_sp):
            ax.set_ylim(dz, z0)

        if len(model.station_x):
            for sx in model.station_x:
                ax_p50.axvline(sx, **sty.station_kwargs())
                ax_sp.axvline(sx,  **sty.station_kwargs())

        fig.suptitle(
            self.title or (
                f"Uncertainty section — {self.quantity}  "
                f"[{unc.method_tag}  N={unc.n_samples}]"
            ),
            fontweight="bold",
        )
        return fig


# ---------------------------------------------------------------------------
# PlotUncertaintyProfile
# ---------------------------------------------------------------------------

class PlotUncertaintyProfile:
    """Profile plot: water-table depth and transmissivity with P10–P90 envelopes.

    Four data series are shown along the profile x-axis:

    * **Top panel** — water-table depth.  Shaded band = P10–P90 range;
      solid line = P50 median.  Optional reference depth.
    * **Bottom panel** — log₁₀(T).  Shaded band = P10–P90; solid line = P50.

    Parameters
    ----------
    unc : UncertaintyResult
    figsize : tuple
    color_wt : str
        Colour for the water-table envelope (default ``'steelblue'``).
    color_T : str
        Colour for the transmissivity envelope (default ``'seagreen'``).
    reference_depth : float, optional
        Horizontal reference line on the WT panel.
    title : str, optional

    Examples
    --------
    >>> fig = PlotUncertaintyProfile(unc_result, reference_depth=20.0).plot()
    """

    def __init__(
        self,
        unc,
        *,
        style: _StyleArg = None,
        color_wt: str | None = None,
        color_T:  str | None = None,
        reference_depth: float | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        self.unc             = unc
        self.style           = style
        self.color_wt        = color_wt
        self.color_T         = color_T
        self.reference_depth = reference_depth
        self.figsize         = figsize
        self.title           = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty   = resolve_profile_style(self.style)
        unc   = self.unc
        model = unc.resistivity_model
        x     = model.x_centers
        cwt   = self.color_wt if self.color_wt is not None else sty.color_wt
        cT    = self.color_T  if self.color_T  is not None else sty.color_T
        fsz   = resolve_figsize(self.figsize, self.style, "profile")

        p10_wt, p50_wt, p90_wt = unc.p10_wt, unc.p50_wt, unc.p90_wt
        p10_T  = np.log10(np.where(unc.p10_T > 0, unc.p10_T, np.nan))
        p50_T  = np.log10(np.where(unc.p50_T > 0, unc.p50_T, np.nan))
        p90_T  = np.log10(np.where(unc.p90_T > 0, unc.p90_T, np.nan))

        fig, (ax_wt, ax_T) = plt.subplots(
            2, 1, figsize=fsz, sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": [1, 1]},
        )

        # ── WT panel ──────────────────────────────────────────────────────
        ax_wt.fill_between(x, p10_wt, p90_wt,
                           **sty.envelope_kwargs(cwt), label="P10–P90")
        ax_wt.plot(x, p50_wt, **sty.line_kwargs(cwt), label="P50")
        ax_wt.scatter(x, p50_wt, **sty.scatter_kwargs(cwt))

        if self.reference_depth is not None:
            ax_wt.axhline(self.reference_depth, **sty.ref_kwargs(),
                          label=f"Ref. {self.reference_depth} m")

        ax_wt.set_ylabel("Water-table depth (m)")
        ax_wt.invert_yaxis()
        ax_wt.legend(fontsize=8, ncol=3)
        ax_wt.grid(**sty.grid_kwargs())
        ax_wt.set_title(
            self.title or (
                f"WT & T uncertainty profile  "
                f"[{unc.method_tag}  N={unc.n_samples}]"
            ),
            fontsize=10,
        )

        # ── T panel ───────────────────────────────────────────────────────
        ax_T.fill_between(x, p10_T, p90_T,
                          **sty.envelope_kwargs(cT), label="P10–P90")
        ax_T.plot(x, p50_T, **sty.line_kwargs(cT), label="P50")
        ax_T.scatter(x, p50_T, **sty.scatter_kwargs(cT))

        ax_T.set_ylabel(r"$\log_{10}(T_{P50}\ /\ \mathrm{m^2\,s^{-1}})$")
        ax_T.set_xlabel("Profile distance (m)")
        ax_T.legend(fontsize=8)
        ax_T.grid(**sty.grid_kwargs())

        if len(model.station_x):
            for sx in model.station_x:
                ax_wt.axvline(sx, **sty.station_kwargs())
                ax_T.axvline(sx,  **sty.station_kwargs())

        return fig


# ---------------------------------------------------------------------------
# PlotPetrophysicalCrossPlot
# ---------------------------------------------------------------------------

class PlotPetrophysicalCrossPlot:
    """ρ vs φ scatter colored by Sw with the fitted Archie/WS model curve.

    Reproduces the cross-plot of Fig. 3b/3c in Chen et al. (2026) for EM data.
    Each point is one model cell; the colour encodes water saturation (or
    depth). The petrophysical model curve is drawn at the mean observed
    saturation. Hashin-Shtrikman bounds are optionally shown as a shaded
    envelope — no existing tool combines all three in one panel.

    Parameters
    ----------
    result : EMHydroResult
    petro : ArchieModel or WaxmanSmitsModel, optional
        Defaults to ``result.config.petro``.
    color_by : str
        ``'saturation'`` (default) or ``'depth'``.
    show_hs_bounds : bool
        Overlay Hashin-Shtrikman bounds (default ``True``).
    rho_matrix : float
        Rock-matrix resistivity for HS bounds (Ω·m; default 5 000).
    depth_range : (float, float), optional
        Restrict to this depth window (m).
    Sw_for_curve : float, optional
        Sw value used to draw the model curve. Defaults to mean(Sw).
    log_rho : bool
        Log₁₀ scale on the y-axis (default ``True``).
    style, figsize, title
    """

    def __init__(
        self,
        result,
        *,
        style: _StyleArg = None,
        petro=None,
        color_by: str = "saturation",
        show_hs_bounds: bool = True,
        rho_matrix: float = 5000.0,
        depth_range: tuple[float, float] | None = None,
        Sw_for_curve: float | None = None,
        log_rho: bool = True,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        self.result         = result
        self.style          = style
        self.petro          = petro
        self.color_by       = color_by
        self.show_hs_bounds = show_hs_bounds
        self.rho_matrix     = rho_matrix
        self.depth_range    = depth_range
        self.Sw_for_curve   = Sw_for_curve
        self.log_rho        = log_rho
        self.figsize        = figsize
        self.title          = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()
        from .petrophysics import (
            ArchieModel,
            HashinShtrikmanBounds,
        )

        sty   = resolve_section_style(self.style)
        res   = self.result
        model = res.resistivity_model
        cfg   = res.config
        petro = self.petro if self.petro is not None else cfg.petro
        rho_w = cfg.rho_w
        fsz   = resolve_figsize(self.figsize, self.style, "crossplot")

        rho_lin = 10.0 ** model.rho_2d
        z_2d    = np.broadcast_to(model.z_centers[:, np.newaxis], rho_lin.shape).copy()

        mask = np.ones(rho_lin.shape, dtype=bool)
        if self.depth_range is not None:
            z_lo, z_hi = self.depth_range
            mask = (z_2d >= z_lo) & (z_2d <= z_hi)

        phi_f = res.porosity[mask].ravel()
        rho_f = rho_lin[mask].ravel()
        Sw_f  = res.saturation[mask].ravel()
        z_f   = z_2d[mask].ravel()

        ok = np.isfinite(rho_f) & np.isfinite(phi_f) & (rho_f > 0) & (phi_f > 0)
        phi_f, rho_f, Sw_f, z_f = phi_f[ok], rho_f[ok], Sw_f[ok], z_f[ok]

        c_vals  = Sw_f  if self.color_by == "saturation" else z_f
        c_label = r"$S_w$" if self.color_by == "saturation" else "Depth (m)"
        y_vals  = np.log10(rho_f) if self.log_rho else rho_f
        y_label = (r"$\log_{10}(\rho\ /\ \Omega\mathrm{m})$"
                   if self.log_rho else r"$\rho$ (Ω·m)")

        fig, ax = plt.subplots(figsize=fsz, layout="constrained")

        if not len(phi_f):
            ax.set_xlabel(r"Porosity $\phi$")
            ax.set_ylabel(y_label)
            ax.set_title(
                self.title or f"Petrophysical cross-plot [{res.method_tag}]",
                fontsize=10,
            )
            ax.text(0.5, 0.5, "No data in selected depth range",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=11, color="0.5")
            return fig

        if self.show_hs_bounds:
            phi_r = np.linspace(max(float(phi_f.min()), 0.01),
                                min(float(phi_f.max()), 0.70), 120)
            hs    = HashinShtrikmanBounds(rho_matrix=self.rho_matrix, rho_fluid=rho_w)
            r_lo, r_hi = hs.bounds(phi_r)
            y_lo = np.log10(r_lo) if self.log_rho else r_lo
            y_hi = np.log10(r_hi) if self.log_rho else r_hi
            ax.fill_between(phi_r, y_lo, y_hi, label="HS bounds",
                            **sty.hs_fill_kwargs())

        sc = ax.scatter(phi_f, y_vals, c=c_vals, **sty.crossplot_scatter_kwargs())
        cb = fig.colorbar(sc, ax=ax, label=c_label, **sty.cb_kwargs())
        cb.ax.tick_params(labelsize=sty.cb_fontsize)

        Sw_c  = (self.Sw_for_curve if self.Sw_for_curve is not None
                 else float(np.nanmean(Sw_f)))
        phi_c = np.linspace(max(float(phi_f.min()), 0.02),
                            min(float(phi_f.max()), 0.70), 200)
        if isinstance(petro, ArchieModel):
            rho_c = petro.forward(phi=phi_c, Sw=Sw_c, rho_w=rho_w)
            lbl_c = f"Archie (Sw={Sw_c:.2f})"
        else:
            sigma_w = 1e3 / rho_w
            rho_c   = petro.forward(phi=phi_c, Sw=Sw_c, sigma_w=sigma_w)
            lbl_c   = f"WS (Sw={Sw_c:.2f})"
        y_c = np.log10(np.clip(rho_c, 1e-2, 1e7)) if self.log_rho else rho_c
        ax.plot(phi_c, y_c, label=lbl_c, **sty.model_curve_kwargs())

        ax.set_xlabel(r"Porosity $\phi$")
        ax.set_ylabel(y_label)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.25)
        ax.set_title(
            self.title or f"Petrophysical cross-plot [{res.method_tag}]",
            fontsize=10,
        )
        return fig


# ---------------------------------------------------------------------------
# PlotAquiferCharacterization
# ---------------------------------------------------------------------------

class PlotAquiferCharacterization:
    """Dar-Zarrouk profile — TR, S, water table, and transmissivity.

    Three or four stacked panels with a shared profile-distance x-axis:

    1. **TR** = Σρᵢhᵢ (Ω·m²) — aquifer productivity indicator.
       Threshold line at ``sty.tr_threshold``.
    2. **S**  = Σhᵢ/ρᵢ (siemens) — clay-seal protective capacity.
       Narain-Mehrotra class lines at ``sty.s_threshold_moderate`` and
       ``sty.s_threshold_good``.
    3. **WT** — water-table depth (m).
    4. **T**  — log₁₀(T) (m²/s), optional.

    Parameters
    ----------
    result : EMHydroResult
    show_transmissivity : bool
        Add a fourth T panel (default ``True``).
    log_TR : bool
        Use a log₁₀ y-axis for the TR panel (default ``True``).
        Recommended when a resistive basement dominates TR and compresses
        the productive-aquifer bars to near-zero on a linear scale.
    reference_depth : float, optional
        Horizontal dashed line on the WT panel.
    style, figsize, title
    """

    def __init__(
        self,
        result,
        *,
        style: _StyleArg = None,
        show_transmissivity: bool = True,
        log_TR: bool = True,
        reference_depth: float | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        self.result              = result
        self.style               = style
        self.show_transmissivity = show_transmissivity
        self.log_TR              = log_TR
        self.reference_depth     = reference_depth
        self.figsize             = figsize
        self.title               = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty    = resolve_profile_style(self.style)
        result = self.result
        model  = result.resistivity_model
        x      = model.x_centers
        fsz    = resolve_figsize(self.figsize, self.style, "aquifer_char")
        n_pan  = 4 if self.show_transmissivity else 3

        fig, axes = plt.subplots(
            n_pan, 1, figsize=fsz, sharex=True, layout="constrained",
            gridspec_kw={"height_ratios": [1] * n_pan},
        )
        ax_TR, ax_S, ax_wt = axes[0], axes[1], axes[2]
        ax_T = axes[3] if self.show_transmissivity else None
        bar_w = np.diff(x).mean() * 0.7 if len(x) > 1 else 50.0

        # 1 — Transverse Resistance
        TR = result.dar_zarrouk_TR
        if self.log_TR:
            tr_plot = np.log10(np.where(TR > 0, TR, np.nan))
            tr_thr  = np.log10(sty.tr_threshold)
            tr_ylabel = r"$\log_{10}(\mathrm{TR}\ /\ \Omega\mathrm{m}^2)$"
        else:
            tr_plot = TR
            tr_thr  = sty.tr_threshold
            tr_ylabel = "TR (Ω·m²)"
        vld_tr = np.isfinite(tr_plot)
        ax_TR.bar(x[vld_tr], tr_plot[vld_tr], width=bar_w,
                  **sty.bar_kwargs(sty.color_TR))
        ax_TR.axhline(tr_thr, color="crimson", lw=1.2, ls="--",
                      label=f">{sty.tr_threshold:.0f} Ω·m² (productive)")
        ax_TR.set_ylabel(tr_ylabel)
        ax_TR.legend(fontsize=7)
        ax_TR.grid(**sty.grid_kwargs())
        ax_TR.set_title(
            self.title or f"Aquifer characterization [{result.method_tag}]",
            fontsize=10,
        )

        # 2 — Longitudinal Conductance + protection classes
        ax_S.bar(x, result.dar_zarrouk_S, width=bar_w,
                 **sty.bar_kwargs(sty.color_S))
        for thr, lbl, ls in [
            (sty.s_threshold_moderate, "Poor | Moderate", ":"),
            (sty.s_threshold_good,     "Moderate | Good", "--"),
        ]:
            ax_S.axhline(thr, color="0.35", lw=1.0, ls=ls, label=lbl)
        ax_S.set_ylabel("S (siemens)")
        ax_S.legend(fontsize=7, ncol=2)
        ax_S.grid(**sty.grid_kwargs())

        # 3 — Water Table
        wt      = result.water_table
        valid   = np.isfinite(wt)
        ax_wt.bar(x[valid], wt[valid], width=bar_w,
                  **sty.bar_kwargs(sty.color_wt))
        ax_wt.scatter(x[valid], wt[valid], **sty.scatter_kwargs(sty.color_wt, s=16))
        if self.reference_depth is not None:
            ax_wt.axhline(self.reference_depth, **sty.ref_kwargs(),
                          label=f"Ref. {self.reference_depth} m")
            ax_wt.legend(fontsize=7)
        ax_wt.set_ylabel("WT depth (m)")
        ax_wt.invert_yaxis()
        ax_wt.grid(**sty.grid_kwargs())

        # 4 — Transmissivity
        if ax_T is not None:
            T_log = np.log10(np.where(result.transmissivity > 0,
                                      result.transmissivity, np.nan))
            vld = np.isfinite(T_log)
            ax_T.bar(x[vld], T_log[vld], width=bar_w, **sty.bar_kwargs(sty.color_T))
            ax_T.scatter(x[vld], T_log[vld], **sty.scatter_kwargs(sty.color_T, s=16))
            ax_T.set_ylabel(r"$\log_{10}(T)$")
            ax_T.set_xlabel("Profile distance (m)")
            ax_T.grid(**sty.grid_kwargs())
        else:
            ax_wt.set_xlabel("Profile distance (m)")

        if len(model.station_x):
            for sx in model.station_x:
                for ax in axes:
                    ax.axvline(sx, **sty.station_kwargs())

        return fig


# ---------------------------------------------------------------------------
# PlotMultiTimeLapseGrid
# ---------------------------------------------------------------------------

class PlotMultiTimeLapseGrid:
    """Grid of EM sections at successive time steps — Fig. 5c/5d equivalent.

    N mini-sections in one row, shared colourbar on the right.

    Parameters
    ----------
    timelapse : TimeLapseEM
    quantity : str
        ``'rho'`` — absolute log₁₀ρ.
        ``'delta_rho'`` — Δlog₁₀ρ from baseline (diverging).
        ``'delta_saturation'`` — ΔSw from baseline (requires *petro*).
    surveys : list of int, optional
        Survey indices to show (default: all).
    baseline_idx : int
        Baseline for delta quantities (default 0).
    petro, rho_w, phi
        Petrophysics for saturation conversion.
    vmin, vmax : float, optional
    depth_max : float, optional
    figsize_panel : (w, h), optional
        Size of each mini panel.
    style, title
    """

    def __init__(
        self,
        timelapse,
        quantity: str = "rho",
        *,
        style: _StyleArg = None,
        surveys: list[int] | None = None,
        baseline_idx: int = 0,
        petro=None,
        rho_w: float = 20.0,
        phi: float = 0.25,
        vmin: float | None = None,
        vmax: float | None = None,
        depth_max: float | None = None,
        figsize_panel: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        if quantity not in ("rho", "delta_rho", "delta_saturation"):
            raise ValueError("quantity must be 'rho', 'delta_rho', or 'delta_saturation'.")
        if quantity == "delta_saturation" and petro is None:
            raise ValueError("petro (ArchieModel) required for 'delta_saturation'.")
        self.timelapse    = timelapse
        self.quantity     = quantity
        self.style        = style
        self.surveys      = surveys
        self.baseline_idx = baseline_idx
        self.petro        = petro
        self.rho_w        = rho_w
        self.phi          = phi
        self.vmin         = vmin
        self.vmax         = vmax
        self.depth_max    = depth_max
        self.figsize_panel = figsize_panel
        self.title        = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty      = resolve_section_style(self.style)
        tl       = self.timelapse
        idx_list = self.surveys if self.surveys is not None else list(range(tl.n_surveys))
        n_panels = len(idx_list)
        if not n_panels:
            raise ValueError("surveys list is empty.")

        pw, ph = (self.figsize_panel if self.figsize_panel is not None
                  else resolve_figsize(None, self.style, "tl_panel"))
        fsz    = (pw * n_panels + 1.2, ph)

        is_delta = self.quantity != "rho"
        ref_m    = tl.surveys[self.baseline_idx]
        z, x     = ref_m.z_centers, ref_m.x_centers
        dz       = self.depth_max or float(z[-1])
        z_mask   = z <= dz

        panels: list[np.ndarray] = []
        for si in idx_list:
            surv = tl.surveys[si]
            if self.quantity == "rho":
                data = surv.rho_2d[z_mask, :]
            elif self.quantity == "delta_rho":
                data = (surv.rho_2d - ref_m.rho_2d)[z_mask, :]
            else:
                from .petrophysics import (
                    ArchieModel as _Archie,
                )
                archie = self.petro if self.petro is not None else _Archie()
                data   = (archie.saturation(10.0 ** surv.rho_2d, self.phi, self.rho_w) -
                          archie.saturation(10.0 ** ref_m.rho_2d, self.phi, self.rho_w))[z_mask, :]
            panels.append(data)

        stack = np.stack(panels, axis=0)
        if is_delta:
            amax = float(np.nanpercentile(np.abs(stack), 98))
            vmax = self.vmax if self.vmax is not None else max(amax, 1e-6)
            vmin = self.vmin if self.vmin is not None else -vmax
            cmap = _cmap_with_bad(sty.cmap_timelapse, sty.nan_color)
        else:
            vmin = self.vmin if self.vmin is not None else float(np.nanpercentile(stack,  2))
            vmax = self.vmax if self.vmax is not None else float(np.nanpercentile(stack, 98))
            cmap = _cmap_with_bad(sty.cmap_K, sty.nan_color)

        extent = [x[0], x[-1], z[z_mask][-1], z[z_mask][0]]

        fig, axes = plt.subplots(
            1, n_panels, figsize=fsz, sharey=True, layout="constrained",
            gridspec_kw={"wspace": 0.04},
        )
        if n_panels == 1:
            axes = [axes]

        im_last = None
        for ax, data, si in zip(axes, panels, idx_list):
            im_last = ax.imshow(data, aspect="auto", extent=extent,
                                cmap=cmap, vmin=vmin, vmax=vmax, origin="upper")
            ax.set_title(tl.labels[si], fontsize=8)
            ax.set_xlabel("Distance (m)", fontsize=7)
            ax.tick_params(labelsize=7)
        axes[0].set_ylabel("Depth (m)")

        cb_label = {"rho": r"$\log_{10}(\rho)$",
                    "delta_rho": r"$\Delta\log_{10}(\rho)$",
                    "delta_saturation": r"$\Delta S_w$"}[self.quantity]
        fig.colorbar(im_last, ax=axes, label=cb_label, fraction=0.015, pad=0.01)

        if len(ref_m.station_x):
            for ax in axes:
                for sx in ref_m.station_x:
                    ax.axvline(sx, **sty.station_kwargs())

        fig.suptitle(self.title or f"Time-lapse grid — {self.quantity}",
                     fontsize=10, fontweight="bold")
        return fig


# ---------------------------------------------------------------------------
# PlotResistivityDepthProfile
# ---------------------------------------------------------------------------

class PlotResistivityDepthProfile:
    """1-D resistivity depth profile with zone shading — Fig. 3a equivalent.

    Plots the EM inversion ρ(z) curve at one station with a fill and
    optional hydraulic zone shading derived from an
    :class:`~pycsamt.interp.hydromodel.EMHydroResult`.

    Parameters
    ----------
    source : ResistivityModel or EMHydroResult
    station : str or int
        Station name or column index.
    depth_max : float, optional
    show_zones : bool
        Shade aquifer / vadose / basement zones from EMHydroResult (``True``).
    borehole : Borehole, optional
        If given, a narrow borehole panel is added on the right.
    log_rho : bool
        Log-scale x-axis (default ``True``).
    style, figsize, title
    """

    _ZONE_COLORS = {
        "aquifer":              "#aaddff",
        "fractured/weathered":  "#d0eedd",
        "vadose/weathered":     "#ffffcc",
        "resistive basement":   "#d3d3d3",
        "clay":                 "#f5deb3",
        "saline":               "#f4a460",
    }

    def __init__(
        self,
        source,
        station: str | int = 0,
        *,
        style: _StyleArg = None,
        depth_max: float | None = None,
        show_zones: bool = True,
        borehole=None,
        log_rho: bool = True,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        self.source    = source
        self.station   = station
        self.style     = style
        self.depth_max = depth_max
        self.show_zones= show_zones
        self.borehole  = borehole
        self.log_rho   = log_rho
        self.figsize   = figsize
        self.title     = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty = resolve_section_style(self.style)
        fsz = resolve_figsize(self.figsize, self.style, "depth_profile")

        from .hydromodel import EMHydroResult
        has_hydro = isinstance(self.source, EMHydroResult)
        model     = self.source.resistivity_model if has_hydro else self.source

        n_panels = 2 if self.borehole is not None else 1
        fig, axes = plt.subplots(
            1, n_panels, figsize=fsz, sharey=True, layout="constrained",
            gridspec_kw={"width_ratios": [3, 1]} if n_panels == 2 else None,
        )
        ax = axes[0] if n_panels == 2 else axes

        # ── column resolution ──────────────────────────────────────
        if isinstance(self.station, str):
            ix = model.station_names.index(self.station)
        else:
            ix = int(self.station)
        name = (model.station_names[ix]
                if model.station_names and ix < len(model.station_names)
                else f"S{ix:03d}")

        z_all  = model.z_centers
        dz     = self.depth_max if self.depth_max is not None else float(z_all[-1])
        z_mask = z_all <= dz
        z      = z_all[z_mask]
        rho_l  = model.rho_2d[z_mask, ix]

        # ── zone shading from EMHydroResult ────────────────────────
        if has_hydro and self.show_zones:
            Sw_col = self.source.saturation[z_mask, ix]
            for iz in range(len(z) - 1):
                sw = Sw_col[iz]
                if sw >= 0.85:
                    color = self._ZONE_COLORS["aquifer"]
                elif sw < 0.25:
                    color = self._ZONE_COLORS["resistive basement"]
                else:
                    color = self._ZONE_COLORS["vadose/weathered"]
                ax.axhspan(z[iz], z[iz + 1], color=color, alpha=0.28, zorder=0)

        # ── ρ curve ────────────────────────────────────────────────
        x_vals = rho_l if self.log_rho else 10.0 ** rho_l
        ax.plot(x_vals, z, **sty.rho_curve_kwargs())
        ax.fill_betweenx(z, x_vals, color=sty.rho_fill_color,
                         alpha=sty.rho_fill_alpha)

        ax.invert_yaxis()
        ax.set_xlabel(r"$\log_{10}(\rho\ /\ \Omega\mathrm{m})$"
                      if self.log_rho else r"$\rho$ (Ω·m)")
        ax.set_ylabel("Depth (m)")
        ax.grid(alpha=0.25, axis="x")
        ax.set_title(self.title or f"ρ depth profile — {name}", fontsize=10)

        # ── optional borehole panel ────────────────────────────────
        if self.borehole is not None and n_panels == 2:
            ax_bh = axes[1]
            for intv in getattr(self.borehole, "intervals", []):
                top    = getattr(intv, "top",      None)
                bottom = getattr(intv, "bottom",   None)
                color  = getattr(intv, "color",    "0.7")
                label  = getattr(intv, "lithology","")
                if top is not None and bottom is not None:
                    ax_bh.barh(y=(top + bottom) / 2, width=1.0,
                               height=bottom - top, color=color,
                               alpha=0.80, edgecolor="0.4")
                    ax_bh.annotate(label,
                                   xy=(0.5, (top + bottom) / 2),
                                   xycoords=("axes fraction", "data"),
                                   ha="center", va="center", fontsize=6)
            ax_bh.set_xlim(0, 1)
            ax_bh.set_xticks([])
            ax_bh.set_title("Borehole", fontsize=8)

        return fig


# ---------------------------------------------------------------------------
# PlotUncertaintyHistogram
# ---------------------------------------------------------------------------

class PlotUncertaintyHistogram:
    """Posterior histograms of key hydro quantities per station.

    Shows the full MC posterior for water-table depth or transmissivity at
    up to 6 stations, with histogram bars, optional KDE, and P10/P50/P90
    vertical lines.  When raw ensemble arrays are not passed, a Gaussian
    approximation is drawn from the stored statistics.

    Parameters
    ----------
    unc : UncertaintyResult
    quantity : str
        ``'water_table'`` (default) or ``'transmissivity'``.
    stations : list of str or int, optional
        Stations to display (default: up to 6 evenly spaced).
    wt_ensemble : ndarray (n_samples, n_x), optional
    T_ensemble  : ndarray (n_samples, n_x), optional
    show_kde : bool
    show_percentiles : bool
    log_x : bool, optional
        Log scale (default ``True`` for T, ``False`` for WT).
    ncols : int
        Subplot columns (default min(n_sta, 3)).
    style, figsize, title
    """

    def __init__(
        self,
        unc,
        quantity: str = "water_table",
        *,
        style: _StyleArg = None,
        stations: list[str | int] | None = None,
        wt_ensemble: np.ndarray | None = None,
        T_ensemble:  np.ndarray | None = None,
        show_kde: bool = True,
        show_percentiles: bool = True,
        log_x: bool | None = None,
        ncols: int | None = None,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
    ) -> None:
        if quantity not in ("water_table", "transmissivity"):
            raise ValueError("quantity must be 'water_table' or 'transmissivity'.")
        self.unc              = unc
        self.quantity         = quantity
        self.style            = style
        self.stations         = stations
        self.wt_ensemble      = wt_ensemble
        self.T_ensemble       = T_ensemble
        self.show_kde         = show_kde
        self.show_percentiles = show_percentiles
        self.log_x            = log_x
        self.ncols            = ncols
        self.figsize          = figsize
        self.title            = title

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()

        sty   = resolve_profile_style(self.style)
        unc   = self.unc
        model = unc.resistivity_model
        n_x   = model.n_x
        names = (model.station_names if model.station_names
                 else [f"S{i:03d}" for i in range(n_x)])

        if self.stations is None:
            step = max(1, n_x // 6)
            idx_list = list(range(0, n_x, step))[:6]
        else:
            idx_list = [
                names.index(s) if isinstance(s, str) else int(s)
                for s in self.stations
            ]

        n_sta  = len(idx_list)
        n_cols = self.ncols if self.ncols is not None else min(n_sta, 3)
        n_rows = int(np.ceil(n_sta / n_cols))
        fsz    = self.figsize or (4.0 * n_cols, 3.5 * n_rows)

        fig, axes = plt.subplots(n_rows, n_cols, figsize=fsz,
                                 layout="constrained")
        axes_flat = np.array(axes).ravel() if n_sta > 1 else [axes]

        is_T   = self.quantity == "transmissivity"
        is_log = self.log_x if self.log_x is not None else is_T
        ens    = self.T_ensemble if is_T else self.wt_ensemble
        mean_a = unc.mean_T  if is_T else unc.mean_wt
        std_a  = unc.std_T   if is_T else unc.std_wt
        p10_a  = unc.p10_T   if is_T else unc.p10_wt
        p50_a  = unc.p50_T   if is_T else unc.p50_wt
        p90_a  = unc.p90_T   if is_T else unc.p90_wt
        color  = sty.color_T if is_T else sty.color_wt

        for plot_i, (ax, ix) in enumerate(zip(axes_flat, idx_list)):
            name = names[ix] if ix < len(names) else f"S{ix:03d}"

            if ens is not None:
                raw = ens[:, ix]
                raw = raw[np.isfinite(raw)]
                if is_log:
                    raw = np.log10(np.where(raw > 0, raw, np.nan))
                    raw = raw[np.isfinite(raw)]
            else:
                mu  = np.log10(max(mean_a[ix], 1e-20)) if is_log else mean_a[ix]
                sig = max(std_a[ix], 1e-12)
                raw = np.random.default_rng(42 + ix).normal(mu, sig, 500)

            if len(raw) < 2:
                ax.set_visible(False)
                continue

            ax.hist(raw, **sty.hist_kwargs(color))

            if self.show_kde:
                try:
                    from scipy.stats import gaussian_kde
                    kde = gaussian_kde(raw)
                    xg  = np.linspace(raw.min(), raw.max(), 200)
                    ax.plot(xg, kde(xg), **sty.kde_kwargs())
                except Exception:
                    pass

            if self.show_percentiles:
                for val, ls, lbl in [
                    (p10_a[ix], ":", "P10"), (p50_a[ix], "-", "P50"), (p90_a[ix], ":", "P90")
                ]:
                    v = np.log10(max(val, 1e-20)) if is_log else val
                    if np.isfinite(v):
                        ax.axvline(v, color=sty.ref_color, lw=1.0,
                                   ls=ls, label=lbl)
                if plot_i == 0:
                    ax.legend(fontsize=6)

            ax.set_title(name, fontsize=8)
            ax.set_xlabel(r"$\log_{10}(T)$" if is_log and is_T
                          else "WT depth (m)", fontsize=7)
            ax.set_ylabel("Density", fontsize=7)
            ax.tick_params(labelsize=7)

        for ax in axes_flat[n_sta:]:
            ax.set_visible(False)

        fig.suptitle(
            self.title or (
                f"Posterior histogram — {self.quantity.replace('_', ' ')}  "
                f"[N={unc.n_samples}]"
            ),
            fontsize=10, fontweight="bold",
        )
        return fig
