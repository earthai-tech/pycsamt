# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Plotting utilities for MARE2DEM data and results.

Port of:
  * ``plotMARE2DEM_SurveyLayout.m`` -> :class:`PlotSurveyLayout`
  * ``m2d_plot_poly.m``             -> :func:`plot_poly`
  * MT response visualization       -> :class:`PlotConvergence`
  * Resistivity section             -> :class:`PlotModel`
  * Response comparison             -> :class:`PlotResponse`
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
        ax: matplotlib.axes.Axes | None = None,
        *,
        savefig: str | Path | None = None,
        dpi: int = 150,
        target_rms: float | None = None,
    ) -> matplotlib.figure.Figure:
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
        ax: matplotlib.axes.Axes | None = None,
        *,
        savefig: str | Path | None = None,
        dpi: int = 150,
        units: str = "km",
    ) -> matplotlib.figure.Figure:
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
    ) -> matplotlib.figure.Figure:
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
    ) -> matplotlib.figure.Figure:
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
    ax: matplotlib.axes.Axes | None = None,
    *,
    linewidth: float = 1.0,
    color: str = "k",
    savefig: str | Path | None = None,
    dpi: int = 150,
) -> matplotlib.axes.Axes:
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
# ---------------------------------------------------------------------------
# PlotModel — 2-D resistivity section
# ---------------------------------------------------------------------------

class PlotModel(Mare2DEMBase):
    """Plot the log10-resistivity 2-D section.

    Renders a colour-filled triangular mesh via
    :func:`matplotlib.axes.Axes.tripcolor` when Triangle
    mesh output files (``.node`` + ``.ele``) are found next
    to the ``.resistivity`` file. Falls back to a histogram
    of region resistivity values when no mesh is present.

    Parameters
    ----------
    model_or_result : ResistivityModel or InversionResult
        Source of resistivity values and mesh location.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.
    """

    def __init__(
        self,
        model_or_result: Any,
        **kwargs,
    ):
        super().__init__(**kwargs)
        from .iotools.resistivity import ResistivityFile
        from .mesh import ResistivityModel
        from .results import InversionResult

        self._rf: ResistivityFile | None = None
        self._workdir: Path | None = None
        self._title_stem: str = ""

        if isinstance(model_or_result, ResistivityModel):
            self._rf = model_or_result._rf
            if model_or_result.path is not None:
                p = Path(model_or_result.path)
                self._workdir = p.parent
                self._title_stem = p.stem
        elif isinstance(
            model_or_result, InversionResult
        ):
            rf_like = model_or_result.model
            if isinstance(rf_like, ResistivityModel):
                self._rf = rf_like._rf
                if rf_like.path is not None:
                    self._title_stem = (
                        Path(rf_like.path).stem
                    )
            elif isinstance(rf_like, ResistivityFile):
                self._rf = rf_like
                self._title_stem = Path(
                    rf_like.resistivity_file
                ).stem
            self._workdir = model_or_result.workdir
        else:
            raise TypeError(
                "model_or_result must be a "
                "ResistivityModel or InversionResult."
            )

    # -------------------------------------------------------

    def _load_mesh(self):
        """Return (nodes, tris_1based, rho_per_tri) or (None, None, None)."""
        if self._workdir is None or self._rf is None:
            return None, None, None

        rf = self._rf
        if rf is None:
            return None, None, None

        stem = Path(rf.poly_file).stem
        node_p = ele_p = None
        for infix in (".1", ""):
            np_ = self._workdir / f"{stem}{infix}.node"
            ep_ = self._workdir / f"{stem}{infix}.ele"
            if np_.exists() and ep_.exists():
                node_p, ele_p = np_, ep_
                break

        if node_p is None:
            return None, None, None

        def _first_data(lines):
            return next(
                (
                    i for i, ln in enumerate(lines)
                    if ln.strip()
                    and not ln.strip().startswith("#")
                ),
                None,
            )

        nlines = node_p.read_text(
            errors="replace"
        ).splitlines()
        ni = _first_data(nlines)
        if ni is None:
            return None, None, None

        n_nodes = int(nlines[ni].split()[0])
        ni += 1
        nrows: list = []
        while len(nrows) < n_nodes and ni < len(nlines):
            ln = nlines[ni].strip()
            ni += 1
            if ln and not ln.startswith("#"):
                nrows.append(
                    [float(v) for v in ln.split()]
                )
        if not nrows:
            return None, None, None
        nodes = np.array(nrows, dtype=float)[:, 1:3]

        elines = ele_p.read_text(
            errors="replace"
        ).splitlines()
        ei = _first_data(elines)
        if ei is None:
            return None, None, None

        ehdr = elines[ei].split()
        ei += 1
        n_tri = int(ehdr[0])
        n_dpt = int(ehdr[1]) if len(ehdr) > 1 else 3
        n_ea = int(ehdr[2]) if len(ehdr) > 2 else 0

        erows: list = []
        rattrs: list = []
        while len(erows) < n_tri and ei < len(elines):
            ln = elines[ei].strip()
            ei += 1
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split()
            erows.append(
                [int(v) for v in parts[1: 1 + n_dpt]]
            )
            if n_ea > 0 and len(parts) > 1 + n_dpt:
                rattrs.append(
                    int(float(parts[1 + n_dpt]))
                )
            else:
                rattrs.append(0)

        if not erows:
            return None, None, None

        tris = np.array(erows, dtype=int)
        ri = np.array(rattrs, dtype=int)

        rho_table = rf.resistivity
        n_reg = len(rho_table)
        rho_vals = np.zeros(len(tris))
        for k, r in enumerate(ri):
            if 0 < r <= n_reg:
                rho_vals[k] = float(rho_table[r - 1, 0])

        return nodes, tris, rho_vals

    # -------------------------------------------------------

    def _plot_hist(self, ax, savefig, dpi):
        """Fallback: histogram of resistivity values."""
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots(figsize=(7, 4))
        else:
            fig = ax.figure

        rf = self._rf
        rho = (
            rf.resistivity[:, 0]
            if (
                rf is not None
                and len(rf.resistivity)
            )
            else np.array([])
        )
        if len(rho):
            ax.hist(
                rho,
                bins=40,
                color="steelblue",
                edgecolor="none",
            )
        ax.set_xlabel("log10 rho (ohm-m)")
        ax.set_ylabel("Count")
        ax.set_title("Resistivity distribution")
        fig.tight_layout()
        if savefig is not None:
            fig.savefig(
                str(savefig), dpi=dpi,
                bbox_inches="tight",
            )
        return fig

    # -------------------------------------------------------

    def plot(
        self,
        ax=None,
        *,
        savefig=None,
        dpi: int = 150,
        cmap: str = "turbo_r",
        vmin: float | None = None,
        vmax: float | None = None,
    ) -> matplotlib.figure.Figure:
        r"""Plot the 2-D resistivity section.

        When Triangle mesh output files (``.node`` + ``.ele``)
        are found next to the ``.resistivity`` file, a
        colour-filled triangular section is rendered. Otherwise
        a histogram of resistivity values is shown.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Target axes. When ``None`` a new figure is created.
        savefig : path-like, optional
            Save the figure to this path.
        dpi : int, default 150
            Output resolution in dots per inch.
        cmap : str, default ``"turbo_r"``
            Colour map for the resistivity section.
        vmin : float, optional
            Lower log10-rho bound for the colour axis.
        vmax : float, optional
            Upper log10-rho bound for the colour axis.

        Returns
        -------
        matplotlib.figure.Figure
        """
        import matplotlib.pyplot as plt

        nodes, tris, rho_vals = self._load_mesh()
        if nodes is None:
            return self._plot_hist(ax, savefig, dpi)

        if ax is None:
            fig, ax = plt.subplots(figsize=(10, 5))
        else:
            fig = ax.figure

        y = nodes[:, 0]
        z = nodes[:, 1]
        els = tris - 1

        rf = self._rf
        gb = (
            rf.global_bounds
            if (
                rf is not None
                and rf.global_bounds is not None
            )
            else None
        )
        _vmin = vmin if vmin is not None else (
            float(gb[0]) if gb is not None else None
        )
        _vmax = vmax if vmax is not None else (
            float(gb[1]) if gb is not None else None
        )

        tc = ax.tripcolor(
            y, z, els,
            facecolors=rho_vals,
            cmap=cmap,
            vmin=_vmin,
            vmax=_vmax,
        )
        cb = plt.colorbar(tc, ax=ax)
        cb.set_label("log10 rho (ohm-m)")
        ax.invert_yaxis()
        ax.set_xlabel("y (m)")
        ax.set_ylabel("Depth (m)")
        title = "MARE2DEM resistivity model"
        if self._title_stem:
            title = f"{title} ({self._title_stem})"
        ax.set_title(title, fontsize=9)
        fig.tight_layout()

        if savefig is not None:
            fig.savefig(
                str(savefig), dpi=dpi,
                bbox_inches="tight",
            )
        return fig


# ---------------------------------------------------------------------------
# PlotResponse — observed vs predicted MT response
# ---------------------------------------------------------------------------

class PlotResponse(Mare2DEMBase):
    """Compare observed and predicted MT responses.

    Generates a per-receiver grid of subplots overlaying TE
    and TM apparent resistivity and phase versus period.

    Parameters
    ----------
    result : InversionResult
        Inversion output containing observed data and
        MARE2DEM predicted response.
    **kwargs :
        Forwarded to :class:`Mare2DEMBase`.
    """

    def __init__(self, result: Any, **kwargs):
        super().__init__(**kwargs)
        from .results import InversionResult

        if not isinstance(result, InversionResult):
            raise TypeError(
                "result must be an InversionResult."
            )
        self._result = result

    def plot(
        self,
        ax=None,
        *,
        savefig=None,
        dpi: int = 150,
        station: str | None = None,
        max_rx: int = 4,
        figsize: tuple | None = None,
    ) -> matplotlib.figure.Figure:
        r"""Plot observed vs predicted MT responses.

        One figure row per receiver, two columns: apparent
        resistivity on the left and phase on the right.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Ignored -- the method creates its own figure.
        savefig : path-like, optional
            Save the figure to this path.
        dpi : int, default 150
            Output resolution in dots per inch.
        station : str, optional
            Plot only this receiver by name. When ``None``,
            the first *max_rx* receivers are plotted.
        max_rx : int, default 4
            Maximum number of receivers when *station* is
            ``None``.
        figsize : (float, float), optional
            Figure size in inches.

        Returns
        -------
        matplotlib.figure.Figure

        Raises
        ------
        ValueError
            When the result contains no MT data.
        """
        import matplotlib.pyplot as plt

        em = self._result.data
        resp = self._result.response

        if em is None or em.mt is None:
            raise ValueError(
                "No MT data found in InversionResult."
            )

        freqs = em.mt.frequencies
        if len(freqs) == 0:
            raise ValueError(
                "MT section has no frequencies."
            )

        periods = 1.0 / freqs
        rx_names = em.mt.receiver_name
        n_rx = len(em.mt.receivers)

        if station is not None:
            rx_set = [
                i
                for i, nm in enumerate(rx_names)
                if nm == station
            ]
            if not rx_set:
                raise ValueError(
                    f"Station {station!r} not found."
                )
        else:
            rx_set = list(range(min(n_rx, max_rx)))

        n_plot = len(rx_set)
        if figsize is None:
            figsize = (12, 4 * n_plot)

        fig, axes = plt.subplots(
            n_plot, 2,
            figsize=figsize,
            squeeze=False,
        )

        obs = em.data
        pred = (
            resp.data if resp is not None else None
        )

        _MODES = [
            (123, 104, "tab:red",  "TE"),
            (125, 106, "tab:blue", "TM"),
        ]

        for row, irx in enumerate(rx_set):
            ax_r = axes[row, 0]
            ax_p = axes[row, 1]
            nm = (
                rx_names[irx]
                if irx < len(rx_names)
                else f"Rx{irx + 1}"
            )

            for c_rho, c_phi, color, lbl in _MODES:
                codes = obs[:, 0].astype(int)
                rxcol = obs[:, 3].astype(int)
                m_rho = (
                    (codes == c_rho)
                    & (rxcol == irx + 1)
                )
                m_phi = (
                    (codes == c_phi)
                    & (rxcol == irx + 1)
                )

                if m_rho.any():
                    fi = obs[m_rho, 1].astype(int) - 1
                    ax_r.errorbar(
                        periods[fi],
                        obs[m_rho, 4],
                        yerr=obs[m_rho, 5],
                        fmt="o",
                        ms=3,
                        color=color,
                        label=f"{lbl} obs",
                        zorder=3,
                    )

                if pred is not None:
                    pc = pred[:, 0].astype(int)
                    pr = pred[:, 3].astype(int)
                    mr = (
                        (pc == c_rho)
                        & (pr == irx + 1)
                    )
                    if mr.any():
                        fi_r = (
                            pred[mr, 1].astype(int) - 1
                        )
                        ax_r.plot(
                            periods[fi_r],
                            pred[mr, 6],
                            "-",
                            color=color,
                            lw=1.2,
                            label=f"{lbl} pred",
                            zorder=2,
                        )

                if m_phi.any():
                    fi = obs[m_phi, 1].astype(int) - 1
                    ax_p.errorbar(
                        periods[fi],
                        obs[m_phi, 4],
                        yerr=obs[m_phi, 5],
                        fmt="o",
                        ms=3,
                        color=color,
                        label=f"{lbl} obs",
                        zorder=3,
                    )

                if pred is not None:
                    pc = pred[:, 0].astype(int)
                    pr = pred[:, 3].astype(int)
                    mp = (
                        (pc == c_phi)
                        & (pr == irx + 1)
                    )
                    if mp.any():
                        fi_r = (
                            pred[mp, 1].astype(int) - 1
                        )
                        ax_p.plot(
                            periods[fi_r],
                            pred[mp, 6],
                            "-",
                            color=color,
                            lw=1.2,
                            label=f"{lbl} pred",
                            zorder=2,
                        )

            for axi, ylabel, ttl in [
                (
                    ax_r,
                    "log10 rho (ohm-m)",
                    f"App. resistivity ({nm})",
                ),
                (ax_p, "Phase (deg)", f"Phase ({nm})"),
            ]:
                axi.set_xscale("log")
                axi.set_xlabel("Period (s)")
                axi.set_ylabel(ylabel)
                axi.set_title(ttl, fontsize=9)
                axi.legend(fontsize=7)
                axi.grid(True, alpha=0.3)

        fig.tight_layout()
        if savefig is not None:
            fig.savefig(
                str(savefig), dpi=dpi,
                bbox_inches="tight",
            )
        return fig
