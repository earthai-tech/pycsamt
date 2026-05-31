# -*- coding: utf-8 -*-
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

from typing import List, Optional, Sequence, Tuple, Union

import numpy as np

from .lithology import Layer, StratigraphicLog, RockDatabase
from ._base import ResistivityModel

__all__ = [
    "PlotStratigraphicLog",
    "PlotFenceDiagram",
    "PlotCalibratedModel",
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
        figsize: Tuple[float, float] = (8, 10),
        depth_unit: str = "m",
        title: Optional[str] = None,
        annotation_kws: Optional[dict] = None,
    ) -> None:
        self.log = log
        self.figsize = figsize
        self.depth_unit = depth_unit
        self.title = title or f"Pseudo-Stratigraphic Log — {log.station_name}"
        self.annotation_kws = annotation_kws or {"fontsize": 8}

    def plot(self):
        """Render and return the matplotlib Figure."""
        _, plt = _require_mpl()
        import matplotlib.patches as mpatches

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
        figsize: Optional[Tuple[float, float]] = None,
        title: str = "Fence Diagram",
        max_depth: Optional[float] = None,
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
        misfit_map: Optional[np.ndarray] = None,
        *,
        figsize: Tuple[float, float] = (12, 10),
        cmap_rho: str = "jet",
        vmin_rho: float = 1.0,
        vmax_rho: float = 5.0,
        title: Optional[str] = None,
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
        import matplotlib.colors as mcolors

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
