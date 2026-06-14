# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Plotting utilities for MARE2DEM data and results.

Port of:
  * ``plotMARE2DEM_SurveyLayout.m`` → :class:`PlotSurveyLayout`
  * ``m2d_plot_poly.m``             → :func:`plot_poly`
  * MT response visualization       → :class:`PlotConvergence` (already functional)
  * Resistivity section             → :class:`PlotModel` (stub until mesh port)
  * Response comparison             → :class:`PlotResponse` (stub)
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np

from .base import Mare2DEMBase

if TYPE_CHECKING:
    import matplotlib.axes
    import matplotlib.figure

__all__ = [
    "PlotConvergence",
    "PlotSurveyLayout",
    "PlotRxParams",
    "PlotTxParams",
    "plot_poly",
    "PlotModel",
    "PlotResponse",
]


# ---------------------------------------------------------------------------
# PlotConvergence (functional — already ported)
# ---------------------------------------------------------------------------

class PlotConvergence(Mare2DEMBase):
    """Plot RMS misfit convergence from a MARE2DEM log.

    Parameters
    ----------
    log_or_result : Mare2DEMLog or InversionResult
        Source of iteration records.

    Examples
    --------
    >>> from pycsamt.models.mare2dem import InversionResult, PlotConvergence
    >>> result = InversionResult("./run")
    >>> pc = PlotConvergence(result)
    >>> fig = pc.plot()
    """

    def __init__(self, log_or_result: Any, **kwargs):
        super().__init__(**kwargs)
        from .log import Mare2DEMLog
        from .results import InversionResult
        if isinstance(log_or_result, Mare2DEMLog):
            self._log = log_or_result
        elif isinstance(log_or_result, InversionResult):
            self._log = log_or_result.log
        else:
            raise TypeError(
                "log_or_result must be a Mare2DEMLog or InversionResult."
            )

    def plot(
        self,
        ax: "matplotlib.axes.Axes | None" = None,
        *,
        savefig: str | Path | None = None,
        dpi: int = 150,
        target_rms: float | None = None,
    ) -> "matplotlib.figure.Figure":
        """Draw the RMS-vs-iteration convergence curve."""
        import matplotlib.pyplot as plt

        if self._log is None or not self._log.iterations:
            raise ValueError("No iteration records to plot.")

        iters = [r.iteration for r in self._log.iterations]
        rms = self._log.rms_history()

        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 4))
        else:
            fig = ax.figure

        ax.plot(iters, rms, "-o", ms=4, label="Normalized RMS")
        if target_rms is not None:
            ax.axhline(target_rms, ls="--", color="gray", label="Target RMS")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("Normalized RMS misfit")
        ax.set_title("MARE2DEM convergence")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()

        if savefig is not None:
            fig.savefig(str(savefig), dpi=dpi, bbox_inches="tight")
        return fig


# ---------------------------------------------------------------------------
# PlotSurveyLayout — port of plotMARE2DEM_SurveyLayout.m
# ---------------------------------------------------------------------------

class PlotSurveyLayout(Mare2DEMBase):
    """Plot survey receiver and transmitter positions on a map.

    Port of ``plotMARE2DEM_SurveyLayout.m`` (map view only; GUI removed).

    Parameters
    ----------
    em : EMDataFile
        Data file supplying receiver / transmitter positions and the
        UTM origin metadata.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.

    Examples
    --------
    >>> from pycsamt.models.mare2dem import read_emdata
    >>> from pycsamt.models.mare2dem.plot import PlotSurveyLayout
    >>> em = read_emdata("survey.emdata")
    >>> fig = PlotSurveyLayout(em).plot()
    """

    def __init__(self, em: Any, **kwargs):
        super().__init__(**kwargs)
        self._em = em

    def plot(
        self,
        ax: "matplotlib.axes.Axes | None" = None,
        *,
        savefig: str | Path | None = None,
        dpi: int = 150,
        units: str = "km",
    ) -> "matplotlib.figure.Figure":
        """Draw the survey map: receivers and transmitters in UTM.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on.
        savefig : path-like, optional
            Save figure to this path.
        dpi : int, default 150
            DPI for the saved figure.
        units : {"m", "km"}, default "km"
            Display units for axes labels.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        em = self._em
        scale = 1e-3 if units == "km" else 1.0
        label = "km" if units == "km" else "m"

        utm = em.utm
        n0, e0 = utm.north0, utm.east0
        theta = utm.theta
        cc = np.cos(np.radians(theta))
        ss = np.sin(np.radians(theta))

        def _to_map(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            """Project profile (x, y) to UTM (E, N)."""
            N = cc * x + (-ss) * y + n0
            E = ss * x + cc * y + e0
            return E * scale, N * scale

        if ax is None:
            fig, ax = plt.subplots(figsize=(8, 6))
        else:
            fig = ax.figure

        handles, labels = [], []

        if em.csem is not None:
            if len(em.csem.receivers):
                E, N = _to_map(em.csem.receivers[:, 0], em.csem.receivers[:, 1])
                h = ax.plot(E, N, "bo", ms=4, label="CSEM Receivers")[0]
                handles.append(h); labels.append("CSEM Receivers")
            if len(em.csem.transmitters):
                E, N = _to_map(em.csem.transmitters[:, 0], em.csem.transmitters[:, 1])
                h = ax.plot(E, N, "rx", ms=8, mew=1.5, label="Transmitters")[0]
                handles.append(h); labels.append("Transmitters")

        if em.mt is not None and len(em.mt.receivers):
            E, N = _to_map(em.mt.receivers[:, 0], em.mt.receivers[:, 1])
            h = ax.plot(E, N, "go", ms=7, label="MT Receivers")[0]
            handles.append(h); labels.append("MT Receivers")

        ax.set_xlabel(f"Easting ({label})")
        ax.set_ylabel(f"Northing ({label})")
        title = getattr(em, "path", None)
        ax.set_title(str(title) if title else "MARE2DEM Survey Layout",
                     fontsize=9)
        ax.legend(handles, labels, fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.set_aspect("equal")
        fig.tight_layout()

        if savefig is not None:
            fig.savefig(str(savefig), dpi=dpi, bbox_inches="tight")
        return fig


# ---------------------------------------------------------------------------
# PlotRxParams — receiver geometry QC (port of plotRxParams sub-function)
# ---------------------------------------------------------------------------

class PlotRxParams(Mare2DEMBase):
    """Plot receiver geometry parameters (x, y, z, θ, α, β).

    Port of the ``plotRxParams`` sub-function in
    ``plotMARE2DEM_SurveyLayout.m``.

    Parameters
    ----------
    em : EMDataFile
        Survey data file.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.
    """

    def __init__(self, em: Any, **kwargs):
        super().__init__(**kwargs)
        self._em = em

    def plot(
        self,
        *,
        fig=None,
        savefig: str | Path | None = None,
        dpi: int = 150,
        units: str = "km",
    ) -> "matplotlib.figure.Figure":
        """Draw 6-panel receiver parameter overview."""
        import matplotlib.pyplot as plt

        em = self._em
        scale = 1e-3 if units == "km" else 1.0
        label = "km" if units == "km" else "m"

        datasets: list[tuple[np.ndarray, str, str]] = []
        if em.csem is not None and len(em.csem.receivers):
            datasets.append((em.csem.receivers, "CSEM", "b"))
        if em.mt is not None and len(em.mt.receivers):
            datasets.append((em.mt.receivers, "MT", "r"))

        if not datasets:
            raise ValueError("No receiver data found in EMDataFile.")

        param_labels = ["x", "y (km)", "z (km)", "θ (deg)", "α (deg)", "β (deg)"]
        n_panels = 6

        if fig is None:
            fig, axes = plt.subplots(n_panels, 1, figsize=(8, 12), sharex=True)
        else:
            axes = fig.axes

        for rx_arr, name, color in datasets:
            y_km = rx_arr[:, 1] * scale
            for panel_idx in range(n_panels):
                ax = axes[panel_idx]
                val = rx_arr[:, panel_idx]
                if panel_idx in (1, 2):
                    val = val * scale
                ax.plot(y_km, val, ".", color=color, ms=4, label=name)
                ax.set_ylabel(param_labels[panel_idx], fontsize=8)
                ax.grid(True, alpha=0.3)

        axes[0].legend(fontsize=8)
        axes[0].set_title("Receiver Parameters", fontsize=9)
        axes[-1].set_xlabel(f"Receiver y position ({label})", fontsize=9)
        if axes[2].get_ylim()[0] < axes[2].get_ylim()[1]:
            axes[2].invert_yaxis()

        fig.tight_layout()
        if savefig is not None:
            fig.savefig(str(savefig), dpi=dpi, bbox_inches="tight")
        return fig


# ---------------------------------------------------------------------------
# PlotTxParams
# ---------------------------------------------------------------------------

class PlotTxParams(Mare2DEMBase):
    """Plot CSEM transmitter geometry parameters (x, y, z, azimuth, dip).

    Port of the ``plotTxParams`` sub-function in
    ``plotMARE2DEM_SurveyLayout.m``.
    """

    def __init__(self, em: Any, **kwargs):
        super().__init__(**kwargs)
        self._em = em

    def plot(
        self,
        *,
        fig=None,
        savefig: str | Path | None = None,
        dpi: int = 150,
        units: str = "km",
    ) -> "matplotlib.figure.Figure":
        """Draw transmitter parameter overview."""
        import matplotlib.pyplot as plt

        em = self._em
        if em.csem is None or len(em.csem.transmitters) == 0:
            raise ValueError("No CSEM transmitter data in EMDataFile.")

        scale = 1e-3 if units == "km" else 1.0
        label = "km" if units == "km" else "m"

        tx = em.csem.transmitters
        n_panels = min(5, tx.shape[1])
        param_labels = ["x", "y", "z", "Azimuth (deg)", "Dip (deg)"]

        if fig is None:
            fig, axes = plt.subplots(n_panels, 1, figsize=(8, 10), sharex=True)
        else:
            axes = fig.axes

        y_km = tx[:, 1] * scale
        for p in range(n_panels):
            val = tx[:, p]
            if p in (0, 1, 2):
                val = val * scale
            axes[p].plot(y_km, val, "bo", ms=4)
            axes[p].set_ylabel(param_labels[p], fontsize=8)
            axes[p].grid(True, alpha=0.3)

        axes[0].set_title("Transmitter Parameters", fontsize=9)
        axes[-1].set_xlabel(f"Transmitter y position ({label})", fontsize=9)
        if axes[2].get_ylim()[0] < axes[2].get_ylim()[1]:
            axes[2].invert_yaxis()

        fig.tight_layout()
        if savefig is not None:
            fig.savefig(str(savefig), dpi=dpi, bbox_inches="tight")
        return fig


# ---------------------------------------------------------------------------
# plot_poly — port of m2d_plot_poly.m
# ---------------------------------------------------------------------------

def plot_poly(
    poly_file: str | Path,
    ax: "matplotlib.axes.Axes | None" = None,
    *,
    linewidth: float = 1.0,
    color: str = "k",
    savefig: str | Path | None = None,
    dpi: int = 150,
) -> "matplotlib.axes.Axes":
    """Plot a Triangle ``.poly`` PSLG mesh file.

    Port of ``m2d_plot_poly.m``.

    Parameters
    ----------
    poly_file : path-like
        Path to the ``.poly`` file.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.
    linewidth : float, default 1.0
    color : str, default "k"
    savefig : path-like, optional
    dpi : int, default 150

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the PSLG drawn.

    Examples
    --------
    >>> from pycsamt.models.mare2dem.plot import plot_poly
    >>> ax = plot_poly("mare2dem.poly")
    """
    import matplotlib.pyplot as plt
    from .iotools.poly import read_poly

    pf = read_poly(poly_file)
    nodes = pf.nodes
    segs = pf.segments

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    if len(nodes) == 0 or len(segs) == 0:
        return ax

    # fast multi-segment plot: insert NaN between segments
    y0 = nodes[segs[:, 0] - 1 if segs.min() >= 1 else segs[:, 0], 0]
    y1 = nodes[segs[:, 1] - 1 if segs.min() >= 1 else segs[:, 1], 0]
    z0 = nodes[segs[:, 0] - 1 if segs.min() >= 1 else segs[:, 0], 1]
    z1 = nodes[segs[:, 1] - 1 if segs.min() >= 1 else segs[:, 1], 1]

    # rebuild with NaN separators
    nan_col = np.full(len(segs), np.nan)
    Y = np.column_stack([y0, y1, nan_col]).ravel()
    Z = np.column_stack([z0, z1, nan_col]).ravel()

    ax.plot(Y, Z, "-", lw=linewidth, color=color)
    ax.set_xlabel("y (m)")
    ax.set_ylabel("z (m)")
    ax.invert_yaxis()
    ax.set_title(str(poly_file))
    ax.set_aspect("equal")

    if savefig is not None:
        ax.figure.savefig(str(savefig), dpi=dpi, bbox_inches="tight")
    return ax


# ---------------------------------------------------------------------------
# PlotModel and PlotResponse — stubs pending mesh parser completion
# ---------------------------------------------------------------------------

class PlotModel(Mare2DEMBase):
    """Plot the log10-resistivity section.

    .. note::
       Full FEM mesh visualization will be implemented when the full
       ``.resistivity`` + Triangle mesh parser is complete.
    """

    def __init__(self, model_or_result: Any, **kwargs):
        super().__init__(**kwargs)
        from .mesh import ResistivityModel
        from .results import InversionResult
        if isinstance(model_or_result, ResistivityModel):
            self._model = model_or_result
        elif isinstance(model_or_result, InversionResult):
            self._model = model_or_result.model
        else:
            raise TypeError(
                "model_or_result must be a ResistivityModel or InversionResult."
            )

    def plot(self, ax=None, *, savefig=None, dpi: int = 150):
        raise NotImplementedError(
            "PlotModel.plot requires the full FEM mesh parser, which will be "
            "added when the Mamba2D mesh builder port is complete."
        )


class PlotResponse(Mare2DEMBase):
    """Compare observed and predicted EM data.

    .. note::
       Full response visualization will be implemented when the
       ``.emdata`` response parser is complete.
    """

    def __init__(self, result: Any, **kwargs):
        super().__init__(**kwargs)
        self._result = result

    def plot(self, ax=None, *, savefig=None, dpi: int = 150):
        raise NotImplementedError(
            "PlotResponse.plot is not yet implemented.  It will be "
            "available after the plotMARE2DEM parser is fully ported."
        )
