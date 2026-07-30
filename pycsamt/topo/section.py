# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""One-call topography-embedded 2-D resistivity section plots.

This module generalizes the terrain-following pipeline in
:mod:`pycsamt.topo.extract`/:mod:`pycsamt.topo.drape`/
:mod:`pycsamt.topo.overlay` into a single entry point that accepts
*any* pycsamt resistivity model or inversion result and *any* source
of station topography, then renders a colour section that drapes
correctly over real terrain.

Two rendering modes are supported, mirroring the two modes already
provided by :mod:`pycsamt.topo.overlay`:

``kind="pcolormesh"`` (default)
    A profile-distance vs. elevation section.  The flat depth grid is
    warped into terrain-following coordinates with
    :func:`~pycsamt.topo.drape.drape_section` and the terrain surface
    is drawn with :func:`~pycsamt.topo.overlay.draw_topo_section`.
``kind="imshow"``
    A station-index vs. depth pseudosection.  The colour grid is left
    flat (depth increasing downward) and a compact elevation profile
    strip is inserted above the axes with
    :func:`~pycsamt.topo.overlay.draw_topo_strip`.

Accepted ``model`` inputs
--------------------------
* ``(x_centers, z_centers, rho_2d)`` — a plain array triple, metres by
  default (see *model_unit*).
* :class:`pycsamt.interp.ResistivityModel` — the package's unified
  method-agnostic 2-D model container.
* :class:`pycsamt.inversion.results.InversionResult` — the
  backend-neutral result (converted via ``result.to_resistivity_model()``).
* A native :class:`pycsamt.models.occam2d.results.InversionResult`
  (``rho_2d`` + ``mesh``) — converted via
  :meth:`~pycsamt.interp.ResistivityModel.from_occam2d`.
* A native 2-D :class:`pycsamt.models.modem.results.InversionResult`
  (``mode == "2d"``) — built from ``model_final``/``model_initial``.
* An AI agent-style result exposing ``pred_rho`` (dict, or
  :class:`pycsamt.agents._base.AgentResult`), e.g. the output of
  :class:`pycsamt.agents.inv3d_agent.Inv3DAgent`.

Accepted topography inputs
---------------------------
* ``sites`` — any Sites / EDI collection accepted by
  :func:`~pycsamt.topo.extract.extract_elevation`.
* ``elevation`` (+ optional ``chainage``) — explicit arrays.
* Model-derived — inferred from air-like cells (``log10(rho)`` above
  *air_log10_threshold*) when neither of the above is given.  This is
  a *relative* terrain profile (no absolute vertical datum) unless
  ``sites``/``elevation`` supply one.

Examples
--------
Plot a backend-neutral inversion result with topography from its
source EDI collection::

    >>> from pycsamt.topo import plot_topo_section
    >>> ax = plot_topo_section(result, sites=sites, depth_max=1500.0)  # doctest: +SKIP

Plot an AI 3-D inversion result, depth-cropped to 1.5 km, as a
pseudosection with an elevation strip::

    >>> ax = plot_topo_section(
    ...     agent_result, sites=sites, kind="imshow",
    ...     model_unit="km", depth_max=1.5,
    ... )  # doctest: +SKIP

Build the data without plotting (e.g. for a custom figure)::

    >>> from pycsamt.topo import build_topo_section
    >>> section = build_topo_section(model, sites=sites)  # doctest: +SKIP
    >>> section.values.shape  # doctest: +SKIP
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import TopoConfig
from .drape import drape_section, interp_elev

__all__ = ["TopoSection", "build_topo_section", "plot_topo_section"]

# log10(rho) threshold above which a cell is treated as "air" when
# inferring terrain directly from a model grid.  Matches the ModEM
# convention (~10**5 ohm.m, i.e. ln(rho) > 11.5) expressed in log10.
_AIR_LOG10_THRESHOLD = 5.0


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------


@dataclass
class TopoSection:
    """Resolved, terrain-embedded 2-D resistivity section.

    Returned by :func:`build_topo_section`.  Carries both the
    terrain-draped grid (for ``pcolormesh``) and the flat cell-centre
    grid plus raw topography arrays (for ``imshow`` / custom plots).

    Attributes
    ----------
    x_nodes_km : ndarray, shape (n_x+1,)
        Profile-distance pcolormesh node positions (km).
    z_draped_km : ndarray, shape (n_z+1, n_x+1)
        Terrain-following elevation node grid (km a.s.l.), ready for
        ``ax.pcolormesh(x_nodes_km, z_draped_km, values)``.
    values : ndarray, shape (n_z, n_x)
        Cell values — ``log10(rho)`` when ``log_rho=True``, linear
        resistivity otherwise.  NaN-masked above terrain when
        ``clip_above_surface=True``.
    x_centers_km, z_centers_km : ndarray
        Flat (undraped) cell-centre coordinates (km), post depth-crop.
    z_nodes_km : ndarray, shape (n_z+1,)
        Flat (undraped) depth node positions (km), post depth-crop.
    chainage_km, elev_km : ndarray
        Resolved topography source arrays (one value per topography
        sample point — station count, which may differ from ``n_x``).
    surface_km : ndarray, shape (n_x+1,)
        Terrain elevation interpolated to ``x_nodes_km`` (km a.s.l.).
    station_x_km : ndarray
        Marker x-positions (km) for station pins.
    station_names : list of str
        Station labels, aligned with ``station_x_km``.
    depth_min_km, depth_max_km : float
        Effective (post-crop) depth range, km below the flat datum.
    exaggeration : float
        Vertical exaggeration applied while draping.
    log_rho : bool
        Whether ``values`` is log10(rho) (``True``) or linear rho.
    method : str
        Source model tag (``"occam2d"``, ``"modem"``, ``"ai"``, ...).
    rms : float
        Inversion RMS misfit, if available; ``nan`` otherwise.
    topo_source : str
        Which topography source was actually used: ``"sites"``,
        ``"array"``, ``"model"``, or ``"flat"``.
    """

    x_nodes_km: np.ndarray
    z_draped_km: np.ndarray
    values: np.ndarray
    x_centers_km: np.ndarray
    z_centers_km: np.ndarray
    z_nodes_km: np.ndarray
    chainage_km: np.ndarray
    elev_km: np.ndarray
    surface_km: np.ndarray
    station_x_km: np.ndarray
    station_names: list[str]
    depth_min_km: float
    depth_max_km: float
    exaggeration: float
    log_rho: bool
    method: str
    rms: float
    topo_source: str


@dataclass
class _GridInfo:
    """Internal adapter output — a flat cell-centre resistivity grid."""

    x_centers: np.ndarray
    z_centers: np.ndarray
    rho_log10: np.ndarray
    station_x: np.ndarray
    station_names: list[str]
    method: str
    rms: float
    unit: str  # "m" or "km"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def build_topo_section(
    model: Any,
    *,
    sites: Any = None,
    elevation: Any = None,
    chainage: Any = None,
    station_names: Sequence[str] | None = None,
    station_x: Any = None,
    topo_source: str = "auto",
    model_unit: str = "m",
    depth_min: float = 0.0,
    depth_max: float | None = None,
    exaggeration: float = 1.0,
    log_rho: bool = True,
    interp_method: str = "linear",
    clip_above_surface: bool = True,
    smooth_sigma: float | tuple[float, float] | None = None,
    air_log10_threshold: float = _AIR_LOG10_THRESHOLD,
) -> TopoSection:
    """Resolve a model + topography source into a terrain-draped section.

    This is the data-building half of :func:`plot_topo_section` — use
    it directly when you want the resolved arrays without a figure.

    Parameters
    ----------
    model : object
        Any of the input forms documented in the module docstring:
        a ``(x_centers, z_centers, rho_2d)`` tuple, a
        :class:`pycsamt.interp.ResistivityModel`, a backend-neutral or
        native Occam2D/ModEM ``InversionResult``, or an AI agent
        result exposing ``pred_rho``.
    sites : object, optional
        Station/EDI collection to extract chainage + elevation from
        (see :mod:`pycsamt.topo.extract`).  Takes priority over
        *elevation* when ``topo_source="auto"``.
    elevation : array_like, optional
        Explicit per-station elevation (m a.s.l.).  Paired with
        *chainage* (km); if *chainage* is omitted, the model's own
        station positions are used.
    chainage : array_like, optional
        Explicit per-station along-profile distance (km).  Only used
        together with *elevation*.
    station_names : sequence of str, optional
        Overrides the station labels carried by *model* / *sites*.
    station_x : array_like, optional
        Overrides the marker x-positions carried by *model*.
    topo_source : {"auto", "sites", "array", "model"}
        Which topography source to use.  ``"auto"`` prefers *sites*,
        then *elevation*, then model-derived (air-cell) inference,
        then falls back to a flat datum with a warning.
    model_unit : {"m", "km"}
        Unit of *model*'s coordinate arrays (and of *depth_min* /
        *depth_max*).  Ignored for AI-style results, whose
        ``depths_km``/coordinates are always treated as km.
    depth_min, depth_max : float, optional
        Depth window to display, in *model_unit* units.  ``None``
        keeps the full depth range.
    exaggeration : float, default 1.0
        Vertical exaggeration applied to the terrain relief.
    log_rho : bool, default True
        Keep values as log10(rho).  ``False`` converts to linear
        Ω·m (``10 ** log10_rho``).
    interp_method : {"linear", "cubic", "nearest"}
        Elevation interpolation method, forwarded to
        :func:`~pycsamt.topo.drape.interp_elev`.
    clip_above_surface : bool, default True
        NaN-mask cells that lie above the local terrain surface.
    smooth_sigma : float or (float, float), optional
        Gaussian-smoothing sigma (depth, distance), applied to
        *values* before draping.  Requires scipy; silently skipped
        with a warning when scipy is unavailable.
    air_log10_threshold : float, default 5.0
        log10(rho) threshold used for model-derived terrain inference
        (``topo_source in {"auto", "model"}``).

    Returns
    -------
    TopoSection

    Raises
    ------
    TypeError
        If *model* is not a recognised input form.
    ValueError
        If *topo_source* is invalid, or is forced to a source that
        cannot be resolved (e.g. ``"sites"`` without *sites*).
    """
    if model_unit not in ("m", "km"):
        raise ValueError("model_unit must be 'm' or 'km'.")

    grid = _extract_grid(
        model,
        station_x=station_x,
        station_names=station_names,
        unit=model_unit,
    )

    x_centers_km = _to_km(grid.x_centers, grid.unit)
    z_centers_km = _to_km(grid.z_centers, grid.unit)
    station_x_km = _to_km(grid.station_x, grid.unit)

    chain_km, elev_m, names, topo_src_used = _resolve_topography(
        grid,
        sites=sites,
        elevation=elevation,
        chainage=chainage,
        station_names=station_names,
        topo_source=topo_source,
        air_log10_threshold=air_log10_threshold,
    )
    elev_km = np.asarray(elev_m, dtype=float) / 1000.0

    values = np.asarray(grid.rho_log10, dtype=float).copy()

    depth_min_km = _to_km(depth_min, model_unit)
    if depth_max is None:
        depth_max_km = float(z_centers_km.max()) if z_centers_km.size else 0.0
    else:
        depth_max_km = _to_km(depth_max, model_unit)

    mask = (z_centers_km >= depth_min_km) & (z_centers_km <= depth_max_km)
    if not np.any(mask):
        warnings.warn(
            f"depth range [{depth_min_km}, {depth_max_km}] km selects no "
            "layers; showing the full depth range instead.",
            UserWarning,
            stacklevel=2,
        )
        mask = np.ones_like(z_centers_km, dtype=bool)
    z_centers_km = z_centers_km[mask]
    values = values[mask, :]

    if smooth_sigma is not None:
        try:
            from scipy.ndimage import gaussian_filter

            values = gaussian_filter(values, sigma=smooth_sigma)
        except ImportError:
            warnings.warn(
                "smooth_sigma was given but scipy is not installed; "
                "skipping smoothing.",
                UserWarning,
                stacklevel=2,
            )

    if not log_rho:
        values = 10.0**values

    x_nodes_km = _cell_edges(x_centers_km)
    z_nodes_km = _cell_edges(z_centers_km)
    if z_nodes_km.size:
        z_nodes_km[0] = max(0.0, float(z_nodes_km[0]))

    elev_at_centres_km = interp_elev(
        chain_km, elev_km, x_centers_km, method=interp_method
    )
    x_nodes_km, z_draped_km, values_draped = drape_section(
        x_nodes_km,
        z_nodes_km,
        values,
        elev_at_centres_km,
        exaggeration=exaggeration,
        clip_above_surface=clip_above_surface,
    )
    surface_km = interp_elev(
        chain_km, elev_km, x_nodes_km, method=interp_method
    )

    return TopoSection(
        x_nodes_km=x_nodes_km,
        z_draped_km=z_draped_km,
        values=values_draped,
        x_centers_km=x_centers_km,
        z_centers_km=z_centers_km,
        z_nodes_km=z_nodes_km,
        chainage_km=chain_km,
        elev_km=elev_km,
        surface_km=surface_km,
        station_x_km=station_x_km,
        station_names=names,
        depth_min_km=depth_min_km,
        depth_max_km=(
            float(z_centers_km.max()) if z_centers_km.size else depth_min_km
        ),
        exaggeration=float(exaggeration),
        log_rho=bool(log_rho),
        method=grid.method,
        rms=grid.rms,
        topo_source=topo_src_used,
    )


def plot_topo_section(
    model: Any,
    *,
    ax: Any = None,
    kind: str = "pcolormesh",
    sites: Any = None,
    elevation: Any = None,
    chainage: Any = None,
    station_names: Sequence[str] | None = None,
    station_x: Any = None,
    topo_source: str = "auto",
    model_unit: str = "m",
    depth_min: float = 0.0,
    depth_max: float | None = None,
    exaggeration: float = 1.0,
    log_rho: bool = True,
    interp_method: str = "linear",
    clip_above_surface: bool = True,
    smooth_sigma: float | tuple[float, float] | None = None,
    air_log10_threshold: float = _AIR_LOG10_THRESHOLD,
    cmap: str = "jet_r",
    vmin: float | None = None,
    vmax: float | None = None,
    vmin_percentile: float = 2.0,
    vmax_percentile: float = 98.0,
    colorbar: bool = True,
    section: Any = "inversion",
    show_stations: bool = True,
    show_station_names: bool = True,
    topo_cfg: TopoConfig | None = None,
    station_marker: Any = None,
    dark: bool = False,
    title: str | None = None,
    figsize: tuple[float, float] | None = None,
    savepath: str | None = None,
    savefig_kw: dict[str, Any] | None = None,
    return_data: bool = False,
):
    """Plot a resistivity model or inversion result draped over topography.

    A single entry point that accepts any of the model forms described
    in the :mod:`pycsamt.topo.section` module docstring (raw arrays,
    :class:`pycsamt.interp.ResistivityModel`, backend-neutral or native
    Occam2D/ModEM ``InversionResult``, AI agent results) together with
    any topography source (``sites``, explicit arrays, or model-derived
    inference), and renders a publication-style section using the
    shared :mod:`pycsamt.api` styling (:data:`pycsamt.api.section.PYCSAMT_SECTION`,
    :data:`pycsamt.api.station.PYCSAMT_STATION_RENDERING`).

    Parameters
    ----------
    model : object
        See the module docstring for accepted forms.
    ax : matplotlib Axes, optional
        Existing axes to draw into.  A new figure/axes is created when
        omitted, sized via the selected *section* style.
    kind : {"pcolormesh", "imshow"}
        ``"pcolormesh"`` drapes the grid over real terrain (profile
        distance vs. elevation).  ``"imshow"`` keeps a flat
        station-index vs. depth pseudosection with a compact elevation
        strip inserted above it.
    sites, elevation, chainage, station_names, station_x :
        Topography source, forwarded to :func:`build_topo_section`.
    topo_source, model_unit :
        Topography-resolution and unit controls, forwarded to
        :func:`build_topo_section`.
    depth_min, depth_max, exaggeration :
        Depth window and vertical exaggeration, forwarded to
        :func:`build_topo_section`.
    log_rho, interp_method, clip_above_surface, smooth_sigma, air_log10_threshold :
        Value scaling and grid-building controls, forwarded to
        :func:`build_topo_section` — see its docstring for details.
    cmap : str, default "jet_r"
        Matplotlib colormap.
    vmin, vmax : float, optional
        Explicit colour-scale limits.  When omitted, computed from
        *vmin_percentile* / *vmax_percentile* of the visible values.
    vmin_percentile, vmax_percentile : float, default 2.0, 98.0
        Percentile bounds used to auto-scale the colour map.
    colorbar : bool, default True
        Draw a colorbar using the shared section colorbar style.
    section : str or pycsamt.api.section.SectionStyle, default "inversion"
        Section style preset name or explicit style object.
    show_stations : bool, default True
        Draw station markers (pins for ``pcolormesh``, the elevation
        strip's markers for ``imshow``).
    show_station_names : bool, default True
        Draw station name labels alongside the markers.
    topo_cfg : pycsamt.topo.config.TopoConfig, optional
        Full terrain-rendering style override (fill colour/alpha, line
        style, marker pad, ...).  Defaults to a config that only
        toggles ``station_pins_at_surface`` from *show_stations* and
        turns off the above-surface fill (see *Notes*).
    station_marker : pycsamt.api.station.StationMarkerStyle, optional
        Station-pin style override, forwarded to
        :func:`~pycsamt.topo.overlay.draw_topo_section` /
        :func:`~pycsamt.topo.overlay.draw_topo_strip`.  Defaults to a
        white-filled, black-edged marker sized for legibility against
        a busy ``cmap`` background; pass a
        :class:`~pycsamt.api.station.StationMarkerStyle` to override.
    dark : bool, default False
        Use light-on-dark label colours for the terrain overlay.
    title : str, optional
        Axes title.  Defaults to a method/rms/topo-source summary.
    figsize : (float, float), optional
        Explicit figure size, overriding the section style's sizing.
    savepath : str, optional
        Save the figure via :func:`pycsamt.api.plot.save_fig`.
    savefig_kw : dict, optional
        Extra keyword arguments forwarded to ``save_fig``.
    return_data : bool, default False
        When ``True``, return ``(ax, TopoSection)`` instead of ``ax``.

    Returns
    -------
    matplotlib.axes.Axes
        Or ``(ax, TopoSection)`` when ``return_data=True``.

    Examples
    --------
    >>> from pycsamt.topo import plot_topo_section
    >>> ax = plot_topo_section(result, sites=sites)  # doctest: +SKIP

    Cropped to the shallow 1.5 km and rendered as a pseudosection:

    >>> ax = plot_topo_section(
    ...     result,
    ...     sites=sites,
    ...     kind="imshow",
    ...     depth_max=1500.0,
    ... )  # doctest: +SKIP
    """
    import matplotlib.pyplot as plt

    from ..api.plot import save_fig
    from ..api.section import PYCSAMT_SECTION, SectionStyle
    from .overlay import draw_topo_section, draw_topo_strip

    kind = str(kind).lower()
    if kind not in ("pcolormesh", "imshow"):
        raise ValueError("kind must be 'pcolormesh' or 'imshow'.")

    section_style = (
        section.copy()
        if isinstance(section, SectionStyle)
        else PYCSAMT_SECTION.style_for(str(section)).copy()
    )

    data = build_topo_section(
        model,
        sites=sites,
        elevation=elevation,
        chainage=chainage,
        station_names=station_names,
        station_x=station_x,
        topo_source=topo_source,
        model_unit=model_unit,
        depth_min=depth_min,
        depth_max=depth_max,
        exaggeration=exaggeration,
        log_rho=log_rho,
        interp_method=interp_method,
        clip_above_surface=clip_above_surface,
        smooth_sigma=smooth_sigma,
        air_log10_threshold=air_log10_threshold,
    )

    finite_values = data.values[np.isfinite(data.values)]
    if vmin is None:
        vmin = (
            float(np.nanpercentile(finite_values, vmin_percentile))
            if finite_values.size
            else 0.0
        )
    if vmax is None:
        vmax = (
            float(np.nanpercentile(finite_values, vmax_percentile))
            if finite_values.size
            else 1.0
        )

    if ax is None:
        fsize = figsize or section_style.figsize_for(
            n_stations=max(len(data.station_names), data.x_centers_km.size, 1),
            n_y=max(data.z_centers_km.size, 1),
            labels=data.station_names,
            colorbar=colorbar,
        )
        fig, ax = plt.subplots(
            figsize=fsize, constrained_layout=section_style.figure.constrained
        )
    else:
        fig = ax.get_figure()

    cbar_label = (
        r"$\log_{10}\rho$ ($\Omega\cdot$m)"
        if log_rho
        else r"$\rho$ ($\Omega\cdot$m)"
    )
    if topo_cfg is None:
        # Deliberately not forwarding `exaggeration` here: drape_section
        # leaves its surface row unscaled by exaggeration (only the depth
        # term is stretched), while draw_topo_section scales absolute
        # elevation by cfg.exaggeration. Leaving this at the TopoConfig
        # default (1.0) keeps the overlay's terrain line/pins aligned with
        # the mesh surface at any `exaggeration` value; wiring the two
        # together would reintroduce the mismatch.
        #
        # fill_alpha=0.0: the space above the terrain line has no data in
        # it (it is literally air), so it is left as plain white/background
        # rather than tinted with the fill colour. clip_below_surface stays
        # True so mismatched-column gaps (see docs/user_guide/topo/concepts)
        # blend into that same white instead of standing out separately.
        topo_cfg = TopoConfig(
            station_pins_at_surface=show_stations, fill_alpha=0.0
        )
    if station_marker is None:
        from ..api.station import PYCSAMT_STATION_RENDERING, StationMarkerStyle

        _base = PYCSAMT_STATION_RENDERING.inversion.marker
        station_marker = StationMarkerStyle(
            marker=_base.marker,
            size=_base.size,
            facecolor="white",
            edgecolor="black",
            linewidth=_base.linewidth,
            alpha=_base.alpha,
            offset=_base.offset,
            zorder=_base.zorder,
        )

    if kind == "pcolormesh":
        im = ax.pcolormesh(
            data.x_nodes_km,
            data.z_draped_km,
            data.values,
            shading="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_xlabel("Profile distance (km)")
        ax.set_ylabel("Elevation (km)")
        if section_style.figure.title:
            # Rotated station names occupy the band immediately above the
            # terrain.  Reserve that band so the section title cannot pass
            # through the labels on dense profiles.
            ax.set_title(
                _title(data, title),
                pad=62 if (show_stations and show_station_names) else None,
            )
        ax.grid(bool(section_style.axis.grid))
        if ax.yaxis_inverted():
            ax.invert_yaxis()

        top = float(np.nanmax(data.surface_km))
        bottom = float(np.nanmin(data.z_draped_km))
        margin = 0.04 * max(top - bottom, 1e-6)
        ax.set_ylim(bottom - margin * 0.25, top + margin)
        ax.set_xlim(float(data.x_nodes_km.min()), float(data.x_nodes_km.max()))

        if colorbar:
            section_style.add_colorbar(im, ax, label=cbar_label)
        draw_topo_section(
            ax,
            data.chainage_km,
            data.elev_km * 1000.0,
            data.station_names if show_station_names else None,
            station_x_km=data.station_x_km,
            cfg=topo_cfg,
            dark=dark,
            marker_style=station_marker,
        )
    else:  # imshow
        n_x = data.x_centers_km.size
        extent = [
            -0.5,
            n_x - 0.5,
            float(data.z_nodes_km[-1]) if data.z_nodes_km.size else 1.0,
            float(data.z_nodes_km[0]) if data.z_nodes_km.size else 0.0,
        ]
        im = ax.imshow(
            data.values,
            extent=extent,
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
        )
        section_style.apply_axis(
            ax,
            xlabel="Station",
            ylabel="Depth (km)",
            title=_title(data, title),
        )
        ax.set_xlim(-0.5, n_x - 0.5)
        if colorbar:
            section_style.add_colorbar(im, ax, label=cbar_label)
        if show_stations:
            elev_strip_m = (
                interp_elev(
                    data.chainage_km,
                    data.elev_km,
                    data.x_centers_km,
                    method=interp_method,
                )
                * 1000.0
            )
            strip_names = (
                data.station_names
                if show_station_names and len(data.station_names) == n_x
                else None
            )
            draw_topo_strip(
                fig,
                ax,
                np.arange(n_x, dtype=float),
                elev_strip_m,
                strip_names,
                cfg=topo_cfg,
                dark=dark,
                marker_style=station_marker,
            )

    if section_style.figure.tight:
        try:
            fig.tight_layout()
        except Exception:
            pass

    if savepath:
        save_fig(ax, savepath, **(savefig_kw or {}))

    if return_data:
        return ax, data
    return ax


# ---------------------------------------------------------------------------
# Internal: model adapters
# ---------------------------------------------------------------------------


def _to_km(value: Any, unit: str) -> Any:
    """Convert a scalar or array from *unit* ("m"/"km") to km."""
    if isinstance(value, np.ndarray):
        return (
            value / 1000.0 if unit == "m" else value.astype(float, copy=True)
        )
    return float(value) / 1000.0 if unit == "m" else float(value)


def _cell_edges(centres: np.ndarray) -> np.ndarray:
    """Return node/edge positions bracketing 1-D cell-centre values."""
    centres = np.asarray(centres, dtype=float).ravel()
    if centres.size == 0:
        return np.asarray([0.0, 1.0])
    if centres.size == 1:
        width = max(abs(float(centres[0])) * 0.1, 0.05)
        return np.asarray([centres[0] - width, centres[0] + width])
    edges = np.empty(centres.size + 1, dtype=float)
    edges[1:-1] = 0.5 * (centres[:-1] + centres[1:])
    edges[0] = centres[0] - (edges[1] - centres[0])
    edges[-1] = centres[-1] + (centres[-1] - edges[-2])
    return edges


def _infer_terrain_km(
    z_centers_km: np.ndarray,
    rho_log10_2d: np.ndarray,
    threshold: float,
) -> np.ndarray | None:
    """Infer a *relative* terrain profile from air-like cells in a grid.

    Returns ``None`` when no column contains a cell above *threshold*,
    signalling that model-derived topography is not available.
    """
    rho_log10_2d = np.asarray(rho_log10_2d, dtype=float)
    z_centers_km = np.asarray(z_centers_km, dtype=float)
    nz, nx = rho_log10_2d.shape
    if nz == 0 or nx == 0:
        return None
    is_earth = rho_log10_2d < threshold
    any_earth = is_earth.any(axis=0)
    if not np.any(any_earth):
        return None
    first_idx = np.argmax(is_earth, axis=0)
    first_idx = np.where(any_earth, first_idx, 0)
    if not np.any(first_idx > 0):
        # No column actually starts with an air cap above threshold —
        # there is no real air layer to infer terrain relief from.
        return None
    return -z_centers_km[first_idx]


def _extract_grid(
    model: Any,
    *,
    station_x: Any = None,
    station_names: Sequence[str] | None = None,
    unit: str = "m",
) -> _GridInfo:
    """Adapt any supported *model* input into a flat :class:`_GridInfo`."""
    # 1. Raw (x_centers, z_centers, rho_2d) triple.
    if isinstance(model, (tuple, list)) and len(model) == 3:
        x_c, z_c, rho = model
        x_c = np.asarray(x_c, dtype=float)
        z_c = np.asarray(z_c, dtype=float)
        rho = np.asarray(rho, dtype=float)
        sx = (
            np.asarray(station_x, dtype=float)
            if station_x is not None
            else x_c
        )
        names = (
            list(station_names)
            if station_names is not None
            else [f"S{i:03d}" for i in range(x_c.size)]
        )
        return _GridInfo(
            x_centers=x_c,
            z_centers=z_c,
            rho_log10=rho,
            station_x=sx,
            station_names=names,
            method="array",
            rms=float("nan"),
            unit=unit,
        )

    if isinstance(model, np.ndarray):
        raise TypeError(
            "A bare 2-D array must be paired with coordinates; pass "
            "model=(x_centers, z_centers, rho_2d) instead."
        )

    # 2. pycsamt.interp.ResistivityModel (or lookalike).
    if (
        hasattr(model, "rho_2d")
        and hasattr(model, "x_centers")
        and hasattr(model, "z_centers")
    ):
        x_c = np.asarray(model.x_centers, dtype=float)
        sx_default = getattr(model, "station_x", None)
        sx = (
            np.asarray(station_x, dtype=float)
            if station_x is not None
            else (
                np.asarray(sx_default, dtype=float)
                if sx_default is not None and len(sx_default)
                else x_c
            )
        )
        names_default = list(getattr(model, "station_names", None) or [])
        names = (
            list(station_names)
            if station_names is not None
            else (names_default or [f"S{i:03d}" for i in range(x_c.size)])
        )
        return _GridInfo(
            x_centers=x_c,
            z_centers=np.asarray(model.z_centers, dtype=float),
            rho_log10=np.asarray(model.rho_2d, dtype=float),
            station_x=sx,
            station_names=names,
            method=str(getattr(model, "method", "generic")),
            rms=float(getattr(model, "rms", float("nan"))),
            unit=unit,
        )

    # 3. Backend-neutral pycsamt.inversion.results.InversionResult.
    to_rm = getattr(model, "to_resistivity_model", None)
    if callable(to_rm):
        rm = to_rm()
        info = _extract_grid(
            rm, station_x=station_x, station_names=station_names, unit="m"
        )
        backend = getattr(model, "backend", "")
        method = getattr(model, "method", info.method)
        info.method = f"{backend}:{method}".strip(":") or info.method
        info.rms = float(getattr(model, "rms", info.rms))
        return info

    # 4. Native Occam2D InversionResult (rho_2d + mesh, no x_centers).
    if (
        hasattr(model, "rho_2d")
        and hasattr(model, "mesh")
        and model.rho_2d is not None
    ):
        from ..interp import ResistivityModel

        rm = ResistivityModel.from_occam2d(model)
        return _extract_grid(
            rm, station_x=station_x, station_names=station_names, unit="m"
        )

    # 5. Native 2-D ModEM InversionResult (mode + model_final/initial).
    if hasattr(model, "mode") and (
        getattr(model, "model_final", None) is not None
        or getattr(model, "model_initial", None) is not None
    ):
        mm = (
            model.model_final
            if getattr(model, "model_final", None) is not None
            else model.model_initial
        )
        if str(model.mode).lower() != "2d":
            raise ValueError(
                "Native 3-D ModEM InversionResult objects hold a full "
                "volume, not a single profile. Extract a 2-D cut first "
                "(e.g. pycsamt.models.modem.section.station_curtain or "
                "pycsamt.models.modem.plot.PlotSection) and pass the "
                "resulting (x_centers, z_centers, rho_2d) grid, or a "
                "pycsamt.interp.ResistivityModel."
            )
        x_widths = np.asarray(mm.x_widths, dtype=float)
        z_widths = np.asarray(mm.z_widths, dtype=float)
        x_nodes = np.concatenate([[0.0], np.cumsum(x_widths)])
        z_nodes = np.concatenate([[0.0], np.cumsum(z_widths)])
        x_c = 0.5 * (x_nodes[:-1] + x_nodes[1:])
        z_c = 0.5 * (z_nodes[:-1] + z_nodes[1:])
        rho_log10 = np.asarray(mm.rho_loge, dtype=float) / np.log(10.0)
        sx = (
            np.asarray(station_x, dtype=float)
            if station_x is not None
            else x_c
        )
        names = (
            list(station_names)
            if station_names is not None
            else [f"S{i:03d}" for i in range(x_c.size)]
        )
        return _GridInfo(
            x_centers=x_c,
            z_centers=z_c,
            rho_log10=rho_log10,
            station_x=sx,
            station_names=names,
            method="modem",
            rms=float("nan"),
            unit=unit,
        )

    # 6. AI agent-style result: dict, or AgentResult (both support
    #    `in` / `.get`), exposing `pred_rho` (n_sta, n_layers), log10 rho.
    getter = getattr(model, "get", None)
    if callable(getter):
        try:
            has_pred_rho = "pred_rho" in model
        except TypeError:
            has_pred_rho = False
        if has_pred_rho:
            pred_rho = np.asarray(model["pred_rho"], dtype=float)
            depths_km = model.get("depths_km")
            depths_km = (
                np.asarray(depths_km, dtype=float)
                if depths_km is not None
                else np.arange(pred_rho.shape[1], dtype=float)
            )
            names = list(
                station_names
                or model.get("station_names")
                or [f"S{i:03d}" for i in range(pred_rho.shape[0])]
            )
            coords = model.get("station_coords")
            if station_x is not None:
                x_c = np.asarray(station_x, dtype=float)
            elif coords is not None:
                coords = np.asarray(coords, dtype=float)
                if coords.shape[0] > 1:
                    seg = np.sqrt((np.diff(coords, axis=0) ** 2).sum(axis=1))
                    x_c = np.concatenate([[0.0], np.cumsum(seg)]) / 1000.0
                else:
                    x_c = np.zeros(coords.shape[0])
            else:
                x_c = np.arange(pred_rho.shape[0], dtype=float)
            rms = model.get("rms_global", float("nan"))
            return _GridInfo(
                x_centers=x_c,
                z_centers=depths_km,
                rho_log10=pred_rho.T,
                station_x=x_c,
                station_names=names,
                method="ai",
                rms=float(rms) if rms is not None else float("nan"),
                unit="km",
            )

    # 7. MARE2DEM triangular mesh — genuinely unsupported (no regular grid).
    mod_module = type(model).__module__
    if "mare2dem" in mod_module:
        raise NotImplementedError(
            "Native MARE2DEM results use an unstructured triangular mesh "
            "and are not natively supported by plot_topo_section(). "
            "Regrid the mesh onto a regular (x_centers, z_centers) grid "
            "first (e.g. with scipy.interpolate.griddata against the "
            ".poly node coordinates), or build a "
            "pycsamt.interp.ResistivityModel manually, then pass that "
            "instead."
        )

    raise TypeError(
        f"Unsupported model type: {type(model)!r}. Expected a "
        "(x_centers, z_centers, rho_2d) tuple, a "
        "pycsamt.interp.ResistivityModel, an object exposing "
        "to_resistivity_model() (pycsamt.inversion.results.InversionResult), "
        "a native Occam2D or 2-D ModEM InversionResult, or an AI agent "
        "result exposing 'pred_rho'."
    )


def _resolve_topography(
    grid: _GridInfo,
    *,
    sites: Any = None,
    elevation: Any = None,
    chainage: Any = None,
    station_names: Sequence[str] | None = None,
    topo_source: str = "auto",
    air_log10_threshold: float = _AIR_LOG10_THRESHOLD,
) -> tuple[np.ndarray, np.ndarray, list[str], str]:
    """Resolve (chainage_km, elev_m, names, source_used) for *grid*."""
    valid_sources = {"auto", "sites", "array", "model"}
    if topo_source not in valid_sources:
        raise ValueError(
            f"topo_source must be one of {sorted(valid_sources)}."
        )

    names = (
        list(station_names)
        if station_names is not None
        else list(grid.station_names)
    )

    if topo_source == "sites" or (topo_source == "auto" and sites is not None):
        if sites is None:
            raise ValueError(
                "topo_source='sites' requires the `sites` argument."
            )
        from .extract import (
            extract_chainage,
            extract_elevation,
            extract_station_names,
        )

        chain_km = extract_chainage(sites)
        elev_m = extract_elevation(sites)
        if not names:
            names = extract_station_names(sites)
        return chain_km, elev_m, names, "sites"

    if topo_source == "array" or (
        topo_source == "auto" and elevation is not None
    ):
        if elevation is None:
            raise ValueError(
                "topo_source='array' requires the `elevation` argument."
            )
        elev_m = np.asarray(elevation, dtype=float)
        if chainage is not None:
            chain_km = np.asarray(chainage, dtype=float)
        else:
            chain_km = _to_km(grid.station_x, grid.unit)
        return chain_km, elev_m, names, "array"

    if topo_source in ("model", "auto"):
        z_km = _to_km(grid.z_centers, grid.unit)
        terrain_km = _infer_terrain_km(
            z_km, grid.rho_log10, air_log10_threshold
        )
        if terrain_km is not None:
            x_km = _to_km(grid.x_centers, grid.unit)
            return x_km, terrain_km * 1000.0, names, "model"
        if topo_source == "model":
            raise ValueError(
                "topo_source='model' requested but no air-like cells "
                f"(log10 rho > {air_log10_threshold}) were found to "
                "infer terrain from."
            )

    warnings.warn(
        "No topography source resolved (no sites/elevation given, and "
        "the model carries no detectable air layer); rendering with a "
        "flat datum. Pass `sites=` or `elevation=` to embed real terrain.",
        UserWarning,
        stacklevel=3,
    )
    x_km = _to_km(grid.x_centers, grid.unit)
    return x_km, np.zeros_like(x_km), names, "flat"


def _title(data: TopoSection, title: str | None) -> str | None:
    if title is not None:
        return title
    rms = "" if not np.isfinite(data.rms) else f", rms={data.rms:.3g}"
    return f"{data.method} section (topo: {data.topo_source}){rms}"
