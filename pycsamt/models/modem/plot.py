# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Diagnostic plotting helpers for ModEM inversion results.

The classes in this module turn an
:class:`~pycsamt.models.modem.results.InversionResult` into
matplotlib figures for checking convergence, inspecting 2-D
and 3-D resistivity models, and viewing apparent-resistivity
or phase responses.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from ...api.section import PYCSAMT_SECTION, SectionStyle
from ...api.station import PYCSAMT_STATION_RENDERING, StationMarkerStyle

from .base import ModEmBase
from .results import InversionResult

__all__ = [
    "PlotMisfit",
    "PlotModel2D",
    "PlotModel3D",
    "PlotSection",
    "PlotResponse",
    "PlotPseudo",
]


# ----------------------------------------------------------------------
# Shared base
# ----------------------------------------------------------------------

class _ModEmPlotBase(ModEmBase):
    """Shared result handling for ModEM plotting helpers."""

    def __init__(self, result: InversionResult | None = None, **kwargs):
        """Initialize a plotter with an optional result object."""
        super().__init__(**kwargs)
        self.result = result

    def _check_result(self) -> InversionResult:
        """Return the attached result or raise a clear error."""
        if self.result is None:
            msg = (
                "No InversionResult attached. Pass result= "
                "to the constructor."
            )
            raise ValueError(msg)
        return self.result


def _resolve_section_style(section: str | SectionStyle) -> SectionStyle:
    """Return a copied section style for ModEM section plots."""
    if isinstance(section, SectionStyle):
        return section.copy()
    return PYCSAMT_SECTION.style_for(str(section)).copy()


# ----------------------------------------------------------------------
# PlotMisfit
# ----------------------------------------------------------------------

class PlotMisfit(_ModEmPlotBase):
    """Plot RMS misfit as a function of inversion iteration.

    ``PlotMisfit`` visualizes the convergence history parsed
    from :class:`~pycsamt.models.modem.log.ModEmLog`. The plot
    shows normalized RMS misfit against iteration number and,
    optionally, marks the lowest RMS value. A reference line at
    :math:`RMS=1` is included because values near one generally
    indicate a data fit comparable to the assigned errors.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded ModEM inversion result. The result must contain
        a parsed ``log`` attribute.
    show_best : bool, default True
        Whether to mark the iteration with the lowest parsed
        RMS value.

    Examples
    --------
    >>> from pycsamt.models.modem.results import InversionResult
    >>> from pycsamt.models.modem.plot import PlotMisfit
    >>> result = InversionResult("modem_run")
    >>> fig = PlotMisfit(result=result).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        show_best: bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.show_best = show_best

    def plot(self):
        """Return a matplotlib figure containing the RMS curve.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing one axes with RMS history.

        Raises
        ------
        ValueError
            If no result is attached or the result has no
            parsed log.
        """
        import matplotlib.pyplot as plt

        r = self._check_result()
        if r.log is None:
            raise ValueError("InversionResult has no log loaded.")

        iters = r.log.iterations
        rms = r.log.rms

        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(iters, rms, "k.-", lw=1.5, ms=5, label="RMS")

        if self.show_best and rms.size:
            bi = int(np.argmin(rms))
            ax.scatter(
                [iters[bi]],
                [rms[bi]],
                color="tab:red",
                zorder=5,
                label=f"best={rms[bi]:.3f}",
            )

        ax.axhline(1.0, color="gray", ls="--", lw=0.8, label="RMS=1")
        ax.set_xlabel("Iteration")
        ax.set_ylabel("RMS misfit")
        ax.set_title(f"ModEM inversion - mode={r.mode}")
        ax.legend(fontsize=8)
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotModel2D
# ----------------------------------------------------------------------

class PlotModel2D(_ModEmPlotBase):
    """Plot a 2-D ModEM resistivity cross-section.

    The plot displays either the initial or final 2-D model as
    a depth section. Resistivity is shown on a logarithmic
    colour scale because magnetotelluric models commonly span
    several orders of magnitude.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result containing ``model_initial`` or
        ``model_final``.
    which : {"final", "initial"}, default "final"
        Which model to display.
    depth_max : float, optional
        Maximum depth in metres to display. If omitted, all
        layers are plotted.
    rho_min, rho_max : float, default 1.0, 1000.0
        Resistivity colour-scale limits in ohm metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap name.

    Examples
    --------
    >>> from pycsamt.models.modem.plot import PlotModel2D
    >>> fig = PlotModel2D(result=result, depth_max=5000).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        which: str = "final",
        depth_max: float | None = None,
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        section: str | SectionStyle = "inversion",
        figsize: tuple[float, float] | None = None,
        show_stations: bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.which = which
        self.depth_max = depth_max
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap
        self.section_style = _resolve_section_style(section)
        self.figsize = figsize
        self.show_stations = show_stations

    def plot(self):
        """Return a matplotlib figure containing the model section."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        model = (
            r.model_final if self.which == "final" else r.model_initial
        )
        if model is None:
            msg = f"No {self.which} model in InversionResult."
            raise ValueError(msg)

        x_nodes = model.x_nodes
        z_nodes = model.z_nodes
        rho = model.rho_linear

        x_km = x_nodes / 1e3
        z_km = z_nodes / 1e3

        # depth crop
        nz = model.nz
        if self.depth_max is not None:
            iz_max = int(np.searchsorted(z_nodes, self.depth_max))
            iz_max = min(iz_max, nz)
        else:
            iz_max = nz

        fig, ax = plt.subplots(
            figsize=self.figsize or self.section_style.figsize_for(
                n_stations=model.nx,
                n_y=iz_max,
                colorbar=True,
            ),
        )
        norm = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)
        pm = ax.pcolormesh(
            x_km,
            z_km[:iz_max + 1],
            rho[:iz_max, :],
            norm=norm,
            cmap=self.cmap,
            shading="flat",
        )
        self.section_style.apply_axis(
            ax,
            xlabel="Offset (km)",
            ylabel="Depth (km)",
            title=f"ModEM 2D model - {self.which}",
        )
        self.section_style.add_colorbar(
            pm,
            ax,
            label="Resistivity (Ω·m)",
        )

        if self.show_stations and r.data_obs is not None:
            mk = PYCSAMT_STATION_RENDERING.style_for("inversion").marker
            sx = (
                r.data_obs.x_coords / 1e3
                + model.x_nodes[-1] / 2e3
            )
            ax.scatter(
                sx,
                np.zeros(len(sx)),
                **mk.kwargs(),
            )

        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotModel3D
# ----------------------------------------------------------------------

class PlotModel3D(_ModEmPlotBase):
    """Plot horizontal slices through a 3-D ModEM model.

    ``PlotModel3D`` extracts model layers nearest to requested
    depths and renders each as an x-y resistivity map. The
    selected depths are interpreted in metres below the model
    top, including the depth coordinate convention stored in
    the model object.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result containing a 3-D model.
    depths : sequence of float, optional
        Depths in metres at which to extract slices. If
        omitted, the first four active earth-layer centres are
        used.
    which : {"final", "initial"}, default "final"
        Which model to display.
    rho_min, rho_max : float, default 1.0, 1000.0
        Resistivity colour-scale limits in ohm metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap name.
    n_cols : int, default 2
        Number of columns in the subplot grid.
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        depths: Sequence[float] | None = None,
        which: str = "final",
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        n_cols: int = 2,
        section: str | SectionStyle = "inversion",
        show_stations: bool = True,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.depths = depths
        self.which = which
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap
        self.n_cols = n_cols
        self.section_style = _resolve_section_style(section)
        self.show_stations = show_stations

    def plot(self):
        """Return a matplotlib figure containing model slices."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        model = (
            r.model_final if self.which == "final" else r.model_initial
        )
        if model is None:
            msg = f"No {self.which} model in InversionResult."
            raise ValueError(msg)

        z_nodes = model.z_nodes
        z_centres = 0.5 * (z_nodes[:-1] + z_nodes[1:])
        n_air = model.n_air

        if self.depths is None:
            active = z_centres[n_air:n_air + 4]
        else:
            active = np.asarray(self.depths, dtype=float)

        # find layer index for each requested depth
        iz_list = [
            int(np.argmin(np.abs(z_centres - d))) for d in active
        ]
        n_slices = len(iz_list)
        n_cols = min(self.n_cols, n_slices)
        n_rows = int(np.ceil(n_slices / n_cols))

        x_km = model.x_nodes / 1e3
        y_km = model.y_nodes / 1e3
        norm = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)

        fig, axes = plt.subplots(
            n_rows,
            n_cols,
            figsize=(5 * n_cols, 4.5 * n_rows),
            squeeze=False,
        )

        # station positions shifted to model-grid coordinates
        sx_km = sy_km = None
        if self.show_stations and r.data_obs is not None:
            cx = model.x_nodes[-1] / 2e3
            cy = model.y_nodes[-1] / 2e3
            sx_km = r.data_obs.x_coords / 1e3 + cx
            sy_km = r.data_obs.y_coords / 1e3 + cy
            mk = PYCSAMT_STATION_RENDERING.style_for("inversion").marker

        for k, iz in enumerate(iz_list):
            row, col = divmod(k, n_cols)
            ax = axes[row][col]
            rho_slice = model.rho_linear[iz, :, :]
            pm = ax.pcolormesh(
                x_km,
                y_km,
                rho_slice,
                norm=norm,
                cmap=self.cmap,
                shading="flat",
            )
            depth_m = float(z_centres[iz])
            self.section_style.add_colorbar(
                pm, ax, label="Resistivity (Ω·m)"
            )
            ax.set_title(f"z = {depth_m/1e3:.2f} km", fontsize=9)
            ax.set_xlabel("X (km)", fontsize=8)
            ax.set_ylabel("Y (km)", fontsize=8)
            ax.set_aspect("equal")

            if sx_km is not None:
                ax.scatter(sx_km, sy_km, **mk.kwargs(s=12))

        # hide unused axes
        for k in range(n_slices, n_rows * n_cols):
            row, col = divmod(k, n_cols)
            axes[row][col].set_visible(False)

        fig.suptitle(f"ModEM 3D model slices - {self.which}", y=1.01)
        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotSection helpers
# ----------------------------------------------------------------------

_AIR_LOGE = 11.5   # ln(~100 000 Ω·m) — cells above this are treated as air


def _detect_terrain(rho_loge_2d: np.ndarray) -> np.ndarray:
    """Return the first non-air layer index for each column of a 2-D slice.

    Parameters
    ----------
    rho_loge_2d : ndarray, shape (nz, ncols)
        Log-natural resistivity values.

    Returns
    -------
    ndarray, shape (ncols,), int
        Per-column index of the first earth layer.
    """
    nz, ncols = rho_loge_2d.shape
    idxs = np.full(ncols, nz, dtype=int)
    for j in range(ncols):
        hits = np.where(rho_loge_2d[:, j] < _AIR_LOGE)[0]
        if hits.size:
            idxs[j] = hits[0]
    return idxs


def _calibrate_z_ref(
    rho_loge_3d: np.ndarray,
    z_widths: np.ndarray,
    x_widths: np.ndarray,
    y_col: int,
    site_coords: dict,
) -> float:
    """Estimate the elevation reference for the model z-axis.

    In the ModEM convention the station Z in the data file is the
    negative of the real elevation above sea level (Z = -100 m →
    100 m above sea level).  Model z-nodes are measured downward from
    the model top.  This function samples each station's N-S column
    at the fixed E-W index ``y_col``, detects the first non-air layer,
    and derives ``z_ref`` from ``elevation_m = z_ref - z_node_at_terrain``.

    Notes
    -----
    In the ModEM 3-D grid (Mod3Dxy convention):

    * ``x_widths`` (nx = 168) — **N-S** (northing) direction.
    * ``y_widths`` (ny =  48) — **E-W** (easting)  direction.
    * ``rho_loge`` shape: (nz, ny, nx).

    Returns
    -------
    float
        Estimated ``z_ref`` in metres, defaulting to 380.0 if no
        stations can be matched.
    """
    z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])
    x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])  # N-S axis
    cx = float(x_nodes[-1]) / 2.0   # N-S model centre

    refs: list[float] = []
    for _name, (x_m, _y_m, z_m) in site_coords.items():
        # station northing (x_m) → model N-S position → x index
        x_model = x_m + cx
        xi = int(np.argmin(np.abs(x_nodes - x_model)))
        xi = max(0, min(xi, rho_loge_3d.shape[2] - 1))  # shape (nz, ny, nx)
        col = rho_loge_3d[:, y_col, xi]
        hits = np.where(col < _AIR_LOGE)[0]
        if not hits.size:
            continue
        # z_ref = terrain_z_node - z_m  (z_m is negative of elevation)
        refs.append(float(z_nodes[hits[0]]) - float(z_m))

    return float(np.median(refs)) if refs else 380.0


# ----------------------------------------------------------------------
# PlotSection
# ----------------------------------------------------------------------

class PlotSection(_ModEmPlotBase):
    """Plot a vertical cross-section with topography and station markers.

    ``PlotSection`` cuts a North-South (or East-West) slice through the
    3-D model at a chosen easting (or northing) offset and displays the
    resistivity, the terrain surface, and the station positions along
    the profile.  It replicates the style of the reference BMP plots
    produced by the ModEM MATLAB toolbox (``pMod3Dlldepsec``), using a
    jet colourmap on a log₁₀ scale from 0 to 3 (1–1000 Ω·m).

    The terrain polygon is inferred directly from the model: cells with
    resistivity above the air threshold (~10⁵ Ω·m, ``rho_loge > 11.5``)
    are treated as air.  The upper boundary of the first earth layer in
    each column is drawn as the terrain surface.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result.
    profile_offset : float, default 0.0
        Offset in metres from the model centre along the axis
        perpendicular to the section.  For a North-South section
        (``direction="NS"``) this is the easting (Y coordinate of
        stations).  For an East-West section it is the northing.
    direction : {"NS", "EW"}, default "NS"
        Orientation of the section.
    which : {"final", "initial"}, default "final"
        Which model to display.
    depth_max : float, default 2000.0
        Maximum depth below sea level to display (metres).
    rho_min, rho_max : float, default 1.0, 1000.0
        Resistivity colour-scale limits (Ω·m, plotted as log₁₀).
    cmap : str, default "jet_r"
        Matplotlib colourmap.
    show_terrain : bool, default True
        Whether to draw the terrain surface and fill the air zone.
    show_stations : bool, default True
        Whether to plot station markers on the terrain surface.
    station_tol : float, default 300.0
        Distance in metres perpendicular to the profile within which
        stations are included.
    figsize : tuple of float, optional
        Figure size in inches.  Defaults to ``(11, 5)``.
    title : str, optional
        Figure title.  A default is derived from the model and profile.

    Examples
    --------
    >>> from pycsamt.models.modem.results import InversionResult
    >>> from pycsamt.models.modem.plot import PlotSection
    >>> result = InversionResult("modem_run")
    >>> # North-South profile at easting +360 m (profile L18)
    >>> fig = PlotSection(result=result, profile_offset=360).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        profile_offset: float = 0.0,
        direction: str = "NS",
        which: str = "final",
        depth_max: float = 2000.0,
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        show_terrain: bool = True,
        show_stations: bool = True,
        show_station_names: bool = True,
        station_tol: float = 300.0,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.profile_offset = float(profile_offset)
        self.direction = direction.upper()
        self.which = which
        self.depth_max = float(depth_max)
        self.rho_min = float(rho_min)
        self.rho_max = float(rho_max)
        self.cmap = cmap
        self.show_terrain = show_terrain
        self.show_stations = show_stations
        self.show_station_names = show_station_names
        self.station_tol = float(station_tol)
        self.figsize = figsize or (11.0, 5.0)
        self.title = title

    def plot(self):
        """Return a matplotlib figure with the resistivity cross-section.

        Returns
        -------
        matplotlib.figure.Figure
            Figure containing one axes with the depth section.

        Raises
        ------
        ValueError
            If no result is attached, or the result has no 3-D final/
            initial model.
        """
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        model = (
            r.model_final if self.which == "final" else r.model_initial
        )
        if model is None:
            raise ValueError(f"No {self.which} model in InversionResult.")

        # In the ModEM Mod3Dxy grid convention used by this dataset:
        #   x_widths (nx=168) → N-S (northing) direction
        #   y_widths (ny= 48) → E-W (easting)  direction
        #   rho_loge shape   → (nz, ny, nx)
        x_nodes = model.x_nodes   # N-S nodes (169 values)
        y_nodes = model.y_nodes   # E-W nodes  (49 values)
        z_nodes = model.z_nodes   # depth downward (66 values)

        cx = float(x_nodes[-1]) / 2.0   # N-S model centre (m)
        cy = float(y_nodes[-1]) / 2.0   # E-W model centre (m)

        if self.direction == "NS":
            # Fix easting (E-W = model y dimension) → nearest y column
            y_target = self.profile_offset + cy
            y_col = int(np.argmin(np.abs(y_nodes - y_target)))
            y_col = max(0, min(y_col, model.ny - 1))
            # N-S slice: rho[:, y_col, :] → shape (nz, nx=168)
            rho_2d = model.rho_loge[:, y_col, :]
            # Along-profile axis: northing from centre
            prof_km = (x_nodes - cx) / 1e3          # 169 nodes
            ax_label = "N-S distance from centre (km)"
        else:
            # Fix northing (N-S = model x dimension) → nearest x row
            x_target = self.profile_offset + cx
            x_row = int(np.argmin(np.abs(x_nodes - x_target)))
            x_row = max(0, min(x_row, model.nx - 1))
            # E-W slice: rho[:, :, x_row] → shape (nz, ny=48)
            rho_2d = model.rho_loge[:, :, x_row]
            # Along-profile axis: easting from centre
            prof_km = (y_nodes - cy) / 1e3          # 49 nodes
            ax_label = "E-W distance from centre (km)"

        # Calibrate elevation axis from data_obs station Z values
        if r.data_obs is not None and self.direction == "NS":
            z_ref = _calibrate_z_ref(
                model.rho_loge,
                model.z_widths,
                model.x_widths,   # N-S widths
                y_col,            # fixed E-W column
                r.data_obs.site_coords,
            )
        else:
            z_ref = 380.0

        # Elevation axis: positive = above sea level (m)
        elev_m = z_ref - z_nodes          # 66 values
        elev_km = elev_m / 1e3

        # Depth crop (keep z where elevation > -depth_max)
        iz_max = int(np.searchsorted(z_nodes, z_ref + self.depth_max))
        iz_max = min(iz_max, model.nz)

        rho_2d_crop = rho_2d[:iz_max, :]
        elev_km_crop = elev_km[:iz_max + 1]

        # Detect terrain per column
        terrain_idx = _detect_terrain(rho_2d_crop)

        # Mask air cells → NaN for display
        rho_disp = rho_2d_crop.copy().astype(float)
        for j in range(rho_disp.shape[1]):
            rho_disp[:terrain_idx[j], j] = np.nan
        log10_rho = rho_disp / np.log(10.0)

        # Profile axis cell centres
        prof_c = 0.5 * (prof_km[:-1] + prof_km[1:])   # ncols centres

        # Terrain profile in elevation (km)
        terrain_elev_km = np.array(
            [
                float(elev_km[terrain_idx[j]])
                for j in range(rho_2d_crop.shape[1])
            ]
        )

        # ------------------------------------------------------------------
        # Build figure
        # ------------------------------------------------------------------
        fig, ax = plt.subplots(figsize=self.figsize)

        norm = mcolors.Normalize(
            vmin=np.log10(self.rho_min), vmax=np.log10(self.rho_max)
        )
        pm = ax.pcolormesh(
            prof_km,
            elev_km_crop,
            log10_rho,
            norm=norm,
            cmap=self.cmap,
            shading="flat",
        )
        cb = fig.colorbar(pm, ax=ax, pad=0.02, shrink=0.85)
        cb.set_label("Resistivity (Ω·m)", fontsize=9)
        ticks = [1, 10, 100, 1000]
        tick_vals = [np.log10(t) for t in ticks if self.rho_min <= t <= self.rho_max]
        cb.set_ticks(tick_vals)
        cb.set_ticklabels([str(t) for t in ticks if self.rho_min <= t <= self.rho_max])

        # Terrain fill (white = air zone above terrain surface)
        if self.show_terrain and terrain_idx.any():
            top_y = float(elev_km_crop[0])
            # Clamp terrain values so padded edge columns don't vanish
            terrain_elev_km_plot = np.where(
                terrain_idx < len(elev_km_crop),
                terrain_elev_km,
                top_y,
            )
            terrain_x = np.concatenate([
                [prof_km[0]],
                prof_c,
                [prof_km[-1]],
            ])
            terrain_y = np.concatenate([
                [top_y],
                terrain_elev_km_plot,
                [top_y],
            ])
            ax.fill_between(
                terrain_x, terrain_y, top_y,
                step=None,
                color="white",
                zorder=3,
                linewidth=0,
            )
            ax.plot(prof_c, terrain_elev_km_plot, color="k", lw=1.0, zorder=4)

        # Station markers, names, and profile x-limits
        # Vertical pads (km) above the detected terrain surface:
        #   _MK_PAD  → raises marker centre so the inverted-triangle tip
        #              lands exactly on the topographic line
        #   _NM_PAD  → bottom of station-name text, comfortably in the
        #              white air zone above the marker
        _MK_PAD = 0.018   # km (~18 m)
        _NM_PAD = 0.050   # km (~50 m)

        st_xs: list[float] = []
        st_ys: list[float] = []   # terrain elevation at each station column
        st_names: list[str] = []
        if self.show_stations and r.data_obs is not None:
            mk = PYCSAMT_STATION_RENDERING.style_for("inversion").marker
            for nm, (x_m, y_m, z_m) in r.data_obs.site_coords.items():
                perp = y_m if self.direction == "NS" else x_m
                along = x_m if self.direction == "NS" else y_m
                if abs(perp - self.profile_offset) <= self.station_tol:
                    # Snap to detected terrain at the nearest profile column
                    # so the marker sits ON the topographic line, not inside
                    # the earth (station Z is placed 0.1 m below terrain).
                    col_idx = int(np.argmin(np.abs(prof_c - along / 1e3)))
                    col_idx = max(0, min(col_idx, len(terrain_elev_km) - 1))
                    st_xs.append(along / 1e3)
                    st_ys.append(float(terrain_elev_km[col_idx]))
                    st_names.append(nm)
            if st_xs:
                ax.scatter(
                    st_xs,
                    [y + _MK_PAD for y in st_ys],
                    **mk.kwargs(s=30, zorder=6, clip_on=False),
                )
                if self.show_station_names:
                    order = np.argsort(st_xs)
                    for i in order:
                        label = st_names[i].split("-")[-1]  # e.g. "001A"
                        ax.text(
                            st_xs[i],
                            st_ys[i] + _NM_PAD,
                            label,
                            rotation=90,
                            ha="center",
                            va="bottom",
                            fontsize=5,
                            color="k",
                            zorder=7,
                            clip_on=False,
                        )

        # Auto-crop horizontal extent to station region + margin
        if st_xs:
            margin = 0.5   # km
            xlo = min(st_xs) - margin
            xhi = max(st_xs) + margin
        else:
            xlo = float(prof_km[0])
            xhi = float(prof_km[-1])
        ax.set_xlim(xlo, xhi)

        # Axes formatting
        ax.set_ylim(float(elev_km_crop[-1]), float(elev_km_crop[0]))
        ax.set_xlabel(ax_label, fontsize=9)
        ax.set_ylabel("Elevation (km)", fontsize=9)

        title = self.title or (
            f"ModEM {self.direction} section — offset {self.profile_offset:+.0f} m "
            f"({self.which})"
        )
        ax.set_title(title, fontsize=10)

        fig.tight_layout()
        return fig


# ----------------------------------------------------------------------
# PlotResponse helpers
# ----------------------------------------------------------------------

_RESP_COMPS = ("ZXX", "ZXY", "ZYX", "ZYY")
_RESP_STYLE_KEY = {"ZXX": "xx", "ZXY": "xy", "ZYX": "yx", "ZYY": "yy"}
_RESP_LATEX = {
    "ZXX": r"$Z_{xx}$", "ZXY": r"$Z_{xy}$",
    "ZYX": r"$Z_{yx}$", "ZYY": r"$Z_{yy}$",
}
_MU0 = 4.0 * np.pi * 1e-7


_ERR_MASK = 1e10   # ModEM uses 2e15 to flag dead/masked data rows


def _collect_z_rows(
    data,
    site_name: str,
    comp: str,
    filter_masked: bool = True,
) -> list:
    """Return sorted [(period, real, imag, error)] for one site + component.

    Parameters
    ----------
    filter_masked : bool, default True
        When True, rows whose ``error`` exceeds ``_ERR_MASK`` (the
        ModEM sentinel ``2e15``) are dropped — appropriate for observed
        data.  Pass ``False`` for forward-response / predicted files
        where the sentinel is always present but the Z values are valid.
    """
    if site_name not in data.site_names:
        return []
    si = data.site_names.index(site_name)
    rows = []
    for blk in data.blocks:
        for row in blk["rows"]:
            if row[1] == si and row[5] == comp:
                err = float(row[8])
                if filter_masked and err >= _ERR_MASK:
                    continue
                rows.append((float(row[0]), float(row[6]),
                              float(row[7]), err))
    rows.sort(key=lambda x: x[0])
    return rows


def _rho_phase_from_rows(rows: list):
    """Convert impedance rows to (periods, rho_a, drho, phi_deg, dphi_deg)."""
    if not rows:
        return None
    p   = np.array([r[0] for r in rows])
    re  = np.array([r[1] for r in rows])
    im  = np.array([r[2] for r in rows])
    err = np.array([r[3] for r in rows])
    omega = 2.0 * np.pi / p
    z2    = re**2 + im**2
    rho_a = z2 / (_MU0 * omega)
    valid = (z2 > 0) & np.isfinite(z2)
    drho  = np.where(valid, 2.0 * np.sqrt(z2) * np.abs(err) / (_MU0 * omega), np.nan)
    phi   = np.degrees(np.arctan2(im, re))
    dphi  = np.where(valid, np.degrees(np.abs(err) / np.sqrt(z2)), np.nan)
    return p, rho_a, drho, phi, dphi


def _z_rms(obs_rows: list, pred_rows: list, rtol: float = 1e-4) -> float:
    """RMS misfit on real + imaginary parts, with period matching tolerance."""
    if not obs_rows or not pred_rows:
        return float("nan")
    pred_p = np.array([r[0] for r in pred_rows])
    sq: list[float] = []
    for p_o, re_o, im_o, err in obs_rows:
        if err <= 0 or not np.isfinite(err):
            continue
        di = np.argmin(np.abs(pred_p - p_o))
        if abs(pred_p[di] - p_o) / max(p_o, 1e-15) < rtol:
            re_p, im_p = pred_rows[di][1], pred_rows[di][2]
            sq.append(((re_o - re_p) / err) ** 2)
            sq.append(((im_o - im_p) / err) ** 2)
    return float(np.sqrt(np.mean(sq))) if sq else float("nan")


# ----------------------------------------------------------------------
# PlotResponse
# ----------------------------------------------------------------------

class PlotResponse(_ModEmPlotBase):
    """Per-station MT response in MTPy style.

    Plots apparent resistivity and phase for all four impedance
    components (Z_xx, Z_xy, Z_yx, Z_yy) of each selected station.
    Each station occupies 4 columns (one per component) built with
    ``GridSpecFromSubplotSpec`` so the ρ_a panel is exactly twice
    the height of the φ panel with zero whitespace between them.

    Observed data is drawn with error bars using component colours
    from :data:`~pycsamt.api.style.PYCSAMT_STYLE`.  When the result
    contains predicted data (``result.data_pred``), the modelled
    response is overlaid as a dotted line in the same colour.  The
    component RMS misfit is shown in each subplot title.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded inversion result.
    stations : sequence of str, optional
        Station names to plot.  Defaults to the first
        ``max_stations`` stations in ``result.data_obs``.
    max_stations : int, default 4
        Maximum number of stations to show.
    show_tipper : bool, default False
        Whether to add a third row for tipper (Re/Im Tx and Ty).
    period_min, period_max : float, optional
        Period range in seconds to display.
    figsize : tuple of float, optional
        Figure size in inches.  Derived automatically if omitted.

    Examples
    --------
    >>> from pycsamt.models.modem.results import InversionResult
    >>> from pycsamt.models.modem.plot import PlotResponse
    >>> result = InversionResult("modem_run")
    >>> fig = PlotResponse(result=result, stations=["23-18-010U"]).plot()
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        stations: Sequence[str] | None = None,
        max_stations: int = 4,
        show_tipper: bool = False,
        period_min: float | None = None,
        period_max: float | None = None,
        figsize: tuple[float, float] | None = None,
        style: str = "modem",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.stations     = stations
        self.max_stations = max_stations
        self.show_tipper  = show_tipper
        self.period_min   = period_min
        self.period_max   = period_max
        self.figsize      = figsize
        self.style        = style

    def plot(self):
        """Return a matplotlib figure with the per-station response panels."""
        import contextlib
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as mgridspec
        from ...api.style import PYCSAMT_STYLE

        r = self._check_result()
        if r.data_obs is None:
            raise ValueError("InversionResult has no data_obs loaded.")

        data_obs  = r.data_obs
        data_pred = r.data_pred
        names = list(self.stations or data_obs.site_names)[: self.max_stations]
        n_st  = len(names)

        inner_hr     = [2, 1, 1] if self.show_tipper else [2, 1]
        n_rows_inner = len(inner_hr)

        # Figure geometry
        col_w = 3.0;  rho_h = 2.4;  phs_h = 1.3
        tip_h = 1.3 if self.show_tipper else 0.0
        gap_h = 0.5
        st_h  = rho_h + phs_h + tip_h + 0.35
        fig_w = col_w * 4 + 1.0
        fig_h = st_h * n_st + gap_h * max(0, n_st - 1) + 0.4
        fig   = plt.figure(
            figsize=self.figsize or (fig_w, max(3.0, fig_h)),
            constrained_layout=False,
        )

        # Apply the named style preset only for the duration of this method.
        _ctx = (
            PYCSAMT_STYLE.context(self.style)
            if self.style and self.style.lower() != "pycsamt"
            else contextlib.nullcontext()
        )
        with _ctx:
            mt = PYCSAMT_STYLE.mt

            _TOP = 0.04 if n_st == 1 else 0.02
            _BOT = 0.06;  _L = 0.07;  _R = 0.01
            outer_gs = mgridspec.GridSpec(
                n_st, 1, figure=fig,
                hspace=gap_h / st_h,
                top=1.0 - _TOP, bottom=_BOT, left=_L, right=1.0 - _R,
            )

            _first_axes: list = []

            for st_i, name in enumerate(names):
                inner_gs = mgridspec.GridSpecFromSubplotSpec(
                    n_rows_inner, 4,
                    subplot_spec=outer_gs[st_i],
                    height_ratios=inner_hr,
                    hspace=0.0,   # ρ_a and φ share zero whitespace
                    wspace=0.10,
                )
                ax_r0 = None

                for ci, comp in enumerate(_RESP_COMPS):
                    cstyle   = getattr(mt, _RESP_STYLE_KEY[comp])
                    obs_rows = _collect_z_rows(data_obs, name, comp)
                    prd_rows = (
                        _collect_z_rows(data_pred, name, comp, filter_masked=False)
                        if data_pred else []
                    )

                    rms     = _z_rms(obs_rows, prd_rows)
                    rms_str = f"  rms={rms:.2f}" if np.isfinite(rms) else ""

                    ax_r = fig.add_subplot(inner_gs[0, ci])
                    ax_p = fig.add_subplot(inner_gs[1, ci], sharex=ax_r)
                    plt.setp(ax_r.get_xticklabels(), visible=False)
                    if ax_r0 is None:
                        ax_r0 = ax_r
                        _first_axes.append((ax_r, name))

                    # ── Observed: error bars in component colour ───────────
                    rp = _rho_phase_from_rows(obs_rows)
                    if rp is not None:
                        p, rho, drho, phi, dphi = rp
                        if self.period_min is not None:
                            m = p >= self.period_min
                            p, rho, drho, phi, dphi = (
                                p[m], rho[m], drho[m], phi[m], dphi[m])
                        if self.period_max is not None:
                            m = p <= self.period_max
                            p, rho, drho, phi, dphi = (
                                p[m], rho[m], drho[m], phi[m], dphi[m])
                        ekw = cstyle.errorbar_kwargs()
                        ekw["label"] = _RESP_LATEX[comp]
                        ax_r.errorbar(p, rho, yerr=drho, **ekw)
                        ekw_p = {k: v for k, v in ekw.items() if k != "label"}
                        ax_p.errorbar(p, phi, yerr=dphi, **ekw_p)

                    # ── Predicted: dotted line in predicted_color ─────────
                    rp2 = _rho_phase_from_rows(prd_rows)
                    if rp2 is not None:
                        pp, rho2, _, phi2, _ = rp2
                        if self.period_min is not None:
                            m2 = pp >= self.period_min
                            pp, rho2, phi2 = pp[m2], rho2[m2], phi2[m2]
                        if self.period_max is not None:
                            m2 = pp <= self.period_max
                            pp, rho2, phi2 = pp[m2], rho2[m2], phi2[m2]
                        # Use predicted_color from style; fall back to obs color
                        pred_c = (
                            cstyle.predicted_color
                            if getattr(cstyle, "predicted_color", "")
                            else cstyle.color
                        )
                        pred_ls = getattr(cstyle, "predicted_ls", ":")
                        label_pred = _RESP_LATEX[comp].rstrip("$") + r"^{m}$"
                        ax_r.plot(
                            pp, rho2,
                            color=pred_c, ls=pred_ls, lw=1.8, alpha=0.85,
                            label=label_pred,
                        )
                        ax_p.plot(
                            pp, phi2,
                            color=pred_c, ls=pred_ls, lw=1.8, alpha=0.85,
                        )

                    # ── Titles / axes ─────────────────────────────────────
                    ax_r.set_title(_RESP_LATEX[comp] + rms_str, fontsize=8, pad=3)
                    ax_r.set_xscale("log");  ax_r.set_yscale("log")
                    ax_p.set_xscale("log")
                    if ci == 0:
                        ax_r.set_ylabel(r"$\rho_a\ (\Omega{\cdot}m)$", fontsize=8)
                        ax_p.set_ylabel(r"$\phi\ (°)$", fontsize=8)
                    else:
                        ax_r.tick_params(labelleft=False)
                        ax_p.tick_params(labelleft=False)
                    ax_p.set_xlabel(r"$T\ (s)$", fontsize=7)
                    ax_r.tick_params(labelsize=6, which="both", top=True)
                    ax_p.tick_params(labelsize=6, which="both")
                    ax_r.legend(fontsize=6, loc="upper left",
                                framealpha=0.6, borderpad=0.3)

            # Station-name annotations (above first ρ_a axis per station)
            for _ax, _nm in _first_axes:
                _ax.annotate(
                    _nm,
                    xy=(0.0, 1.0), xycoords="axes fraction",
                    xytext=(0.0, 1.06), textcoords="axes fraction",
                    fontsize=8, fontweight="bold",
                    va="bottom", ha="left", annotation_clip=False,
                )

            # Small legend note at the top-right corner
            fig.text(
                1.0 - _R, 1.0 - _TOP + 0.005,
                r"obs $\bullet$ — pred $\cdots$",
                fontsize=7, ha="right", va="bottom",
                transform=fig.transFigure, color="#444444",
            )

        return fig


# ----------------------------------------------------------------------
# PlotPseudo
# ----------------------------------------------------------------------

class PlotPseudo(_ModEmPlotBase):
    """Plot apparent-resistivity and phase pseudo-sections.

    ``PlotPseudo`` selects one component from observed data and
    arranges apparent resistivity and phase on station offset
    versus period grids. This is a quick survey-scale view of
    lateral and period-dependent response variation.

    Parameters
    ----------
    result : InversionResult, optional
        Loaded result containing observed data.
    component : str
        Data component to display (e.g. ``'TE'``, ``'ZXY'``).
    rho_min, rho_max : float, default 1.0, 1000.0
        Apparent-resistivity colour-scale limits in ohm
        metres.
    cmap : str, default "jet_r"
        Matplotlib colourmap used for apparent resistivity.
    """

    def __init__(
        self,
        result: InversionResult | None = None,
        component: str = "TE",
        rho_min: float = 1.0,
        rho_max: float = 1000.0,
        cmap: str = "jet_r",
        **kwargs,
    ):
        super().__init__(result=result, **kwargs)
        self.component = component
        self.rho_min = rho_min
        self.rho_max = rho_max
        self.cmap = cmap

    def plot(self):
        """Return a matplotlib figure containing pseudo-sections."""
        import matplotlib.colors as mcolors
        import matplotlib.pyplot as plt

        r = self._check_result()
        if r.data_obs is None:
            raise ValueError("InversionResult has no data_obs loaded.")

        data = r.data_obs
        periods = data.periods
        site_names = data.site_names
        offsets = data.offsets / 1e3

        n_per = len(periods)
        n_site = len(site_names)
        rho_mat = np.full((n_per, n_site), np.nan)
        phs_mat = np.full((n_per, n_site), np.nan)
        mu0 = 4 * np.pi * 1e-7

        for blk in data.blocks:
            for row in blk["rows"]:
                comp = row[5]
                if comp != self.component:
                    continue
                period = row[0]
                si = row[1]
                real, imag = row[6], row[7]
                pi_idx = np.argmin(np.abs(periods - period))
                z2 = real ** 2 + imag ** 2
                rho_mat[pi_idx, si] = z2 / (mu0 * 2 * np.pi / period)
                phs_mat[pi_idx, si] = np.degrees(np.arctan2(imag, real))

        log_p = np.log10(periods)

        # sort stations by easting so pcolormesh gets monotonic x-axis
        sort_idx = np.argsort(offsets)
        offsets = offsets[sort_idx]
        rho_mat = rho_mat[:, sort_idx]
        phs_mat = phs_mat[:, sort_idx]

        fig, (ax_rho, ax_phs) = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
        )

        norm_rho = mcolors.LogNorm(vmin=self.rho_min, vmax=self.rho_max)
        pm_r = ax_rho.pcolormesh(
            offsets,
            log_p,
            rho_mat,
            norm=norm_rho,
            cmap=self.cmap,
            shading="nearest",
        )
        fig.colorbar(pm_r, ax=ax_rho, label="rho_a (ohm m)")
        ax_rho.set_ylabel("log10 Period (s)")
        ax_rho.set_title(
            f"ModEM pseudo-section - {self.component} rho_a"
        )

        pm_p = ax_phs.pcolormesh(
            offsets,
            log_p,
            phs_mat,
            vmin=0,
            vmax=90,
            cmap="RdBu_r",
            shading="nearest",
        )
        fig.colorbar(pm_p, ax=ax_phs, label="Phase (deg)")
        ax_phs.set_ylabel("log10 Period (s)")
        ax_phs.set_xlabel("Offset (km)")
        ax_phs.set_title(
            f"ModEM pseudo-section - {self.component} phase"
        )

        fig.tight_layout()
        return fig
