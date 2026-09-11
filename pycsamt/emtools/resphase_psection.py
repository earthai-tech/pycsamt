# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""MTPy-style apparent-resistivity and phase pseudo-sections.

:class:`PlotResPhasePseudoSection` draws the classic magnetotelluric
pseudo-section: apparent resistivity above phase, one column per
impedance component, and one stacked ``[resistivity, phase]`` row-pair
per site group.  By default every group shares one period (vertical)
axis, one resistivity colour scale, and one phase colour scale; the
station axis (names + inverted-triangle markers) is drawn on top, from
:data:`pycsamt.api.PYCSAMT_STATION_RENDERING`.
"""

from __future__ import annotations

import copy
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

from ..api.control import wrap_phase
from ..api.labels import LOG10_PERIOD_LABEL, PERIOD_LABEL, STATION_LABEL
from ..api.station import PYCSAMT_STATION_RENDERING
from ..core.base import CoreObject
from ._core import _iter_items, _name, ensure_sites
from .inspect import _df_resphase

__all__ = [
    "PlotResPhasePseudoSection",
    "plot_res_phase_pseudosection",
]

# Recognised impedance components, in canonical order.
_COMPONENTS: tuple[str, ...] = ("xy", "yx", "xx", "yy")

# Panel titles.
_Z_LABEL: dict[str, str] = {
    "xy": r"$Z_{xy}$",
    "yx": r"$Z_{yx}$",
    "xx": r"$Z_{xx}$",
    "yy": r"$Z_{yy}$",
}

# pyCSAMT default colormaps -- the same pair ``plot_mt_composite_section``
# and the apparent-resistivity helpers in ``pycsamt.emtools.advanced`` use.
_DEFAULT_RES_CMAP = "jet_r"
_DEFAULT_PHASE_CMAP = "plasma"

# Named phase windows: ``(vmin, vmax)`` in degrees.
_PHASE_PRESETS: dict[str, tuple[float, float]] = {
    "0-90": (0.0, 90.0),
    "-45-45": (-45.0, 45.0),
    "-90-90": (-90.0, 90.0),
    "-180-180": (-180.0, 180.0),
}

# Log-period tolerance for snapping a station's samples onto the shared
# grid -- the same 0.3-decade window used by ``plot_mt_composite_section``.
_LOG_SNAP = 0.30


class PlotResPhasePseudoSection(CoreObject):
    r"""Apparent-resistivity / phase pseudo-section, MTPy style.

    Each **column** is one impedance component; within a column the
    apparent-resistivity panel sits above the phase panel.  Passing
    several site groups stacks them as extra ``[resistivity, phase]``
    row-pairs down the figure.  Every panel shares the period axis; by
    default every resistivity panel shares one colour scale and every
    phase panel shares another, so the columns and groups are directly
    comparable.

    Parameters
    ----------
    sites : Sites-like, or mapping of str to Sites-like
        One survey (any input :func:`~pycsamt.emtools.ensure_sites`
        accepts -- a directory, an EDI list, a ``Sites``), or a mapping
        ``{group_label: survey}`` to stack several groups vertically.
    components : sequence of {"xy", "yx", "xx", "yy"}, optional
        Which components to draw, **in the order given**.  ``None``
        (default) auto-detects: every component that carries finite data
        in at least one group is plotted, in canonical ``xy, yx, xx, yy``
        order.
    res_phase_ratio : float, default ``2/3``
        Height of a resistivity panel relative to its phase panel.  The
        default makes resistivity two-thirds the height of phase; pass
        ``2.0`` to give resistivity two-thirds of the row-pair and phase
        the remaining third, or ``1.0`` for equal panels.
    station_side : {"top", "bottom", "none"}, default "top"
        Where the station axis (ticks, names, and inverted-triangle
        markers) is drawn.  ``"top"`` matches the pyCSAMT section
        convention (:data:`~pycsamt.api.PYCSAMT_STATION_RENDERING`) and
        gives every stacked line its own station axis above its
        resistivity panel; ``"bottom"`` places it under each phase panel;
        ``"none"`` omits it.
    station_markers : bool, default True
        Draw the triangle station markers from
        :data:`~pycsamt.api.PYCSAMT_STATION_RENDERING`.  ``False`` keeps
        only the tick names.
    panel_labels : bool or sequence of str, default False
        Add a per-group corner tag.  ``True`` uses ``(a)``, ``(b)``, …;
        a sequence supplies the tags verbatim.
    grid : bool, default False
        Draw faint horizontal period grid lines on every panel.
    show_colorbar : bool, default True
        Draw the two shared colourbars (resistivity, phase) on the right.
    phase_range : str or tuple of float, default ``"0-90"``
        Phase colour limits in degrees: one of the presets ``"0-90"``,
        ``"-45-45"``, ``"-90-90"``, ``"-180-180"``, the string
        ``"auto"`` (2nd/98th percentiles), or an explicit
        ``(vmin, vmax)`` tuple.  The ``yx`` and ``yy`` phases are always
        rotated by 180 degrees first, so a 1-D response plots near
        45 degrees in every component (the MTPy convention).
    res_range : tuple of float, optional
        Resistivity colour limits in :math:`\Omega\cdot\mathrm{m}`.
        ``None`` uses the 5th/95th percentiles.
    log_res : bool, default True
        Map :math:`\log_{10}\rho_a` rather than :math:`\rho_a` to colour.
    log_period : bool, default False
        Label the vertical axis with :math:`\log_{10}T` on a linear axis.
        ``False`` keeps :math:`T` in seconds on a logarithmic axis.
    res_cmap, phase_cmap : str or Colormap, optional
        Colormaps.  ``None`` uses the pyCSAMT defaults
        (``"jet_r"`` for resistivity, ``"plasma"`` for phase).
    share_scales : bool, default True
        Share one resistivity colour scale and one phase colour scale
        across every panel.  When ``False`` each panel is scaled on its
        own data.
    res_scale_from : {"offdiag", "all"}, default "offdiag"
        When *share_scales* is on, whether the shared resistivity limits
        are taken from the off-diagonal components only (``xy``/``yx`` --
        the physically meaningful ones) or from every drawn component.
        Ignored when only diagonal components are present.
    period_range : tuple of float, optional
        ``(T_min, T_max)`` in seconds; samples outside are dropped.
    station_order : sequence of str, optional
        Explicit left-to-right station order.  ``None`` keeps the order
        the sites are read in.
    n_grid : int, default 60
        Number of log-spaced period rows per group.
    share_period : bool, default True
        Resample every group onto one common period grid and share the
        vertical axis.  Set ``False`` when stacking lines from different
        instruments (e.g. an AMT line above an LMT line) so each keeps
        its own period window instead of a mostly-empty shared axis.
    station_label_step : int, optional
        Show every *n*-th station label.  ``None`` picks a step that
        keeps at most ~12 labels.
    figsize : tuple of float, optional
    title : str, optional
        Figure suptitle.
    axes : sequence of matplotlib.axes.Axes, optional
        A flat, row-major grid of ``2 * n_groups * n_components`` axes to
        draw into (resistivity then phase for each column, top group
        first).  When given, the figure and its two shared colourbars are
        the caller's responsibility.  ``None`` (default) builds a fresh
        figure.

    Examples
    --------
    All present components of one survey:

    >>> from pycsamt.emtools import PlotResPhasePseudoSection
    >>> fig = PlotResPhasePseudoSection(
    ...     "data/MT/broken-hill/edis"
    ... ).plot()

    Two lines stacked, resistivity twice the height of phase, station
    markers above each line, lettered panels:

    >>> fig = PlotResPhasePseudoSection(
    ...     {"Line A": sites_a, "Line B": sites_b},
    ...     components=["xy", "yx"],
    ...     res_phase_ratio=2.0,
    ...     station_side="top",
    ...     panel_labels=True,
    ... ).plot()
    """

    def __init__(
        self,
        sites: Any,
        *,
        components: list[str] | tuple[str, ...] | None = None,
        res_phase_ratio: float = 2.0 / 3.0,
        phase_range: str | tuple[float, float] = "0-90",
        res_range: tuple[float, float] | None = None,
        log_res: bool = True,
        log_period: bool = False,
        res_cmap: Any = None,
        phase_cmap: Any = None,
        share_scales: bool = True,
        res_scale_from: str = "offdiag",
        period_range: tuple[float, float] | None = None,
        station_order: list[str] | None = None,
        n_grid: int = 60,
        share_period: bool = True,
        station_label_step: int | None = None,
        station_side: str = "top",
        station_markers: bool = True,
        panel_labels: bool | list[str] | tuple[str, ...] = False,
        grid: bool = False,
        show_colorbar: bool = True,
        figsize: tuple[float, float] | None = None,
        title: str | None = None,
        axes: Any = None,
        **kwargs: Any,
    ) -> None:
        super().__init__(**kwargs)
        self._sites = sites
        self.components = (
            list(components) if components is not None else None
        )
        self.res_phase_ratio = max(float(res_phase_ratio), 1e-3)
        self.phase_range = phase_range
        self.res_range = res_range
        self.log_res = bool(log_res)
        self.log_period = bool(log_period)
        self.res_cmap = res_cmap or _DEFAULT_RES_CMAP
        self.phase_cmap = phase_cmap or _DEFAULT_PHASE_CMAP
        self.share_scales = bool(share_scales)
        self.res_scale_from = str(res_scale_from)
        self.period_range = period_range
        self.station_order = (
            list(station_order) if station_order is not None else None
        )
        self.n_grid = max(int(n_grid), 4)
        self.share_period = bool(share_period)
        self.station_label_step = station_label_step
        self.station_side = str(station_side).lower()
        self.station_markers = bool(station_markers)
        self.panel_labels = panel_labels
        self.grid = bool(grid)
        self.show_colorbar = bool(show_colorbar)
        self.figsize = figsize
        self.title = title
        self.axes = axes
        if self.station_side not in ("top", "bottom", "none"):
            raise ValueError(
                "station_side must be 'top', 'bottom', or 'none'; "
                f"got {station_side!r}."
            )

    # -- data -------------------------------------------------------

    def _groups(self) -> dict[str, Any]:
        """Return an ordered ``{label: sites-like}`` mapping."""
        if isinstance(self._sites, dict):
            return dict(self._sites)
        return {"": self._sites}

    def _survey_frame(self, sites: Any):
        """Return ``(stations, resphase DataFrame with a period column)``."""
        S = ensure_sites(sites, recursive=True)
        df = _df_resphase(S, kind="resphase")
        if df is None or df.empty:
            return [], None
        df = df.copy()
        df["period"] = 1.0 / df["freq"].replace(0, np.nan)
        if self.period_range is not None:
            lo, hi = self.period_range
            df = df[(df["period"] >= float(lo)) & (df["period"] <= float(hi))]
        if df.empty:
            return [], None
        seen = [
            _name(ed, i) for i, ed in enumerate(_iter_items(S))
        ]
        stations = [s for s in dict.fromkeys(seen) if s in set(df["station"])]
        if not stations:
            stations = list(dict.fromkeys(df["station"]))
        if self.station_order is not None:
            stations = [s for s in self.station_order if s in stations]
        return stations, df

    def _grid_image(self, df, stations, per_grid, column: str) -> np.ndarray:
        """Snap one column onto the shared ``(n_grid, n_station)`` grid."""
        lpg = np.log10(per_grid)
        out = np.full((len(per_grid), len(stations)), np.nan)
        if column not in df.columns:
            return out
        for si, st in enumerate(stations):
            sub = df[df["station"] == st]
            ps = sub["period"].to_numpy(float)
            vs = sub[column].to_numpy(float)
            ok = np.isfinite(ps) & np.isfinite(vs) & (ps > 0)
            ps, vs = ps[ok], vs[ok]
            if ps.size == 0:
                continue
            lps = np.log10(ps)
            for gi, lg in enumerate(lpg):
                j = int(np.argmin(np.abs(lps - lg)))
                if abs(lps[j] - lg) < _LOG_SNAP:
                    out[gi, si] = vs[j]
        return out

    @staticmethod
    def _robust_period_window(df) -> tuple[float, float] | None:
        """Return an outlier-tolerant ``(T_min, T_max)`` for one frame.

        Keeps only the half-decade period bins in which at least a
        fraction of the line's stations carry a finite apparent
        resistivity, so a lone station's stray high-frequency segment in
        an otherwise long-period line does not stretch the grid across
        empty decades.
        """
        per = df["period"].to_numpy(float)
        rho_cols = [c for c in df.columns if c.startswith("rho_")]
        has = np.zeros(len(df), dtype=bool)
        for c in rho_cols:
            has |= np.isfinite(df[c].to_numpy(float))
        ok = np.isfinite(per) & (per > 0) & has
        if not ok.any():
            return None
        lp = np.log10(per[ok])
        sta = df["station"].to_numpy()[ok]
        n_st = len(set(df["station"]))
        thr = max(2, int(np.ceil(0.15 * n_st)))
        bins: dict[float, set] = {}
        for b, s in zip(np.round(lp * 2.0) / 2.0, sta):
            bins.setdefault(float(b), set()).add(s)
        keep = sorted(b for b, ss in bins.items() if len(ss) >= thr)
        if not keep:
            keep = sorted(bins)
        return 10.0 ** min(keep), 10.0 ** (max(keep) + 0.25)

    def _build(self):
        """Assemble every panel image plus each group's period grid.

        Returns
        -------
        dict
            ``components`` and ``groups`` -- the latter a list of
            ``{"label", "stations", "res", "phase", "period_grid"}``
            records, ``res``/``phase`` keyed by component.  When
            *share_period* is on every group carries the same grid.
        """
        groups_in = self._groups()
        frames = {
            lab: self._survey_frame(s) for lab, s in groups_in.items()
        }
        frames = {
            lab: fr for lab, fr in frames.items() if fr[1] is not None
        }
        if not frames:
            raise ValueError(
                "no apparent-resistivity / phase data in the input"
            )

        # per-group period window, then one shared grid or one each
        windows = {
            lab: self._robust_period_window(df)
            for lab, (_, df) in frames.items()
        }
        windows = {k: v for k, v in windows.items() if v is not None}
        if not windows:
            raise ValueError("no finite periods in the input")
        if self.share_period:
            lo = min(w[0] for w in windows.values())
            hi = max(w[1] for w in windows.values())
            grids = {
                lab: np.logspace(
                    np.log10(lo), np.log10(hi), self.n_grid
                )
                for lab in frames
            }
        else:
            grids = {
                lab: np.logspace(
                    np.log10(windows[lab][0]),
                    np.log10(windows[lab][1]),
                    self.n_grid,
                )
                for lab in windows
            }

        # resolve components
        comps = self.components
        if comps is None:
            comps = [
                c
                for c in _COMPONENTS
                if any(
                    np.isfinite(df[f"rho_{c}"]).any()
                    for _, df in frames.values()
                    if f"rho_{c}" in df.columns
                )
            ]
        if not comps:
            raise ValueError("no usable impedance component found")

        groups = []
        for lab, (stations, df) in frames.items():
            per_grid = grids[lab]
            res: dict[str, np.ndarray] = {}
            phase: dict[str, np.ndarray] = {}
            for c in comps:
                res[c] = self._grid_image(
                    df, stations, per_grid, f"rho_{c}"
                )
                ph = self._grid_image(df, stations, per_grid, f"phi_{c}")
                # bring yx/yy out of the third quadrant so a 1-D response
                # sits near +45 deg in every component (MTPy convention)
                if c in ("yx", "yy"):
                    ph = ph + 180.0
                phase[c] = wrap_phase(ph, (-180.0, 180.0))
            groups.append(
                {
                    "label": lab,
                    "stations": stations,
                    "res": res,
                    "phase": phase,
                    "period_grid": per_grid,
                }
            )
        return {"components": comps, "groups": groups}

    # -- limits ----------------------------------------------------

    def _res_limits(self, built) -> tuple[float, float]:
        if self.res_range is not None:
            lo, hi = float(self.res_range[0]), float(self.res_range[1])
        else:
            comps = built["components"]
            pick = comps
            if self.res_scale_from == "offdiag":
                off = [c for c in comps if c in ("xy", "yx")]
                pick = off or comps
            vals = np.concatenate(
                [
                    g["res"][c].ravel()
                    for g in built["groups"]
                    for c in pick
                ]
            )
            vals = vals[np.isfinite(vals) & (vals > 0)]
            if vals.size == 0:
                lo, hi = 1.0, 1000.0
            else:
                lo, hi = np.percentile(vals, [5, 95])
        if self.log_res:
            lo = np.log10(max(lo, 1e-6))
            hi = np.log10(max(hi, 10 ** (lo + 0.5) if self.log_res else lo))
        return float(lo), float(hi)

    def _phase_limits(self, built) -> tuple[float, float]:
        pr = self.phase_range
        if isinstance(pr, str):
            if pr == "auto":
                vals = np.concatenate(
                    [
                        g["phase"][c].ravel()
                        for g in built["groups"]
                        for c in built["components"]
                    ]
                )
                vals = vals[np.isfinite(vals)]
                if vals.size == 0:
                    return 0.0, 90.0
                return tuple(
                    float(v) for v in np.percentile(vals, [2, 98])
                )
            if pr not in _PHASE_PRESETS:
                raise ValueError(
                    f"phase_range must be a (vmin, vmax) tuple, 'auto', "
                    f"or one of {sorted(_PHASE_PRESETS)}; got {pr!r}."
                )
            return _PHASE_PRESETS[pr]
        return float(pr[0]), float(pr[1])

    # -- decoration helpers --------------------------------------

    def _panel_label_list(self, n: int) -> list[str]:
        """Return the per-group corner tags, or ``[]`` when disabled."""
        pl = self.panel_labels
        if not pl:
            return []
        if pl is True:
            return [f"({chr(97 + i)})" for i in range(n)]
        return [str(x) for x in pl][:n]

    def _station_style(self, ci: int, n_col: int, n_st: int):
        """Return a configured copy of the pseudo-section station style."""
        st = copy.deepcopy(
            PYCSAMT_STATION_RENDERING.style_for("pseudosection")
        )
        st.show_markers = self.station_markers
        if self.station_side == "bottom":
            st.side = "bottom"
            # the "pseudosection" preset marker sits just above the axis
            # (axes-fraction 1.025); flip it below for a bottom axis
            st.marker.offset = -0.03
            st.xlabel = STATION_LABEL if ci == n_col // 2 else ""
        else:
            st.side = "top"
            # a repeated "Station" caption above every column only adds
            # noise -- the triangle markers already read as stations
            st.xlabel = ""
        # columns are narrow, so thin the names harder than the global
        # station-axis default would -- and prefer a step that divides
        # ``n_st - 1`` so the last tick lands on the final station
        # instead of being crammed next to the previous one
        if self.station_label_step:
            st.every = int(self.station_label_step)
        elif n_st > 3:
            target = n_st / (6.0 if n_col <= 2 else 4.5)
            divs = [d for d in range(1, n_st) if (n_st - 1) % d == 0]
            st.every = min(
                divs, key=lambda d: abs(d - target)
            ) if divs else max(1, round(target))
        else:
            st.every = 1
        return st

    # -- plot ----------------------------------------------------

    def plot(self):
        """Render the pseudo-section.

        Returns
        -------
        matplotlib.figure.Figure

        Raises
        ------
        ValueError
            If the input carries no resistivity/phase data, if no
            component is usable, or if *phase_range* is not recognised.
        """

        built = self._build()
        comps = built["components"]
        groups = built["groups"]
        n_col = len(comps)
        n_grp = len(groups)
        side = self.station_side
        top_axis = side == "top"
        bot_axis = side == "bottom"

        r_lo, r_hi = self._res_limits(built)
        p_lo, p_hi = self._phase_limits(built)
        ratio = self.res_phase_ratio
        plabels = self._panel_label_list(n_grp)

        ylabel = (
            LOG10_PERIOD_LABEL if self.log_period else PERIOD_LABEL
        )

        # ---- figure + nested grid (one inner grid per site group) -----
        # margins are sized in inches, then converted to figure
        # fractions, so a tall (high res_phase_ratio) figure does not
        # grow a huge blank band above the panels
        n_needed = 2 * n_grp * n_col
        show_cbar = self.show_colorbar and self.axes is None
        gs_left = 0.13 if n_grp > 1 else 0.10
        gs_right = 0.865 if show_cbar else 0.965
        row_h = 2.4 * (ratio + 1.0)  # one [res, phase] pair, inches
        top_in = (0.72 if self.title else 0.32) + (
            0.62 if top_axis else 0.0
        )
        bot_in = 0.85 if bot_axis else 0.48
        gap_in = (1.15 if (top_axis or bot_axis) else 0.2) * (n_grp - 1)

        if self.axes is not None:
            axlist = [
                a for a in np.asarray(self.axes, dtype=object).ravel()
            ]
            if len(axlist) < n_needed:
                raise ValueError(
                    f"axes must supply at least {n_needed} axes "
                    f"(2 x {n_grp} groups x {n_col} components); "
                    f"got {len(axlist)}."
                )
            fig = axlist[0].figure
            inners = None
            fig_h = fig.get_figheight()
        else:
            axlist = None
            fig_h = row_h * n_grp + top_in + bot_in + gap_in
            figsize = self.figsize or (2.2 * n_col + 1.9, fig_h)
            fig = plt.figure(figsize=figsize)
            fig_h = figsize[1]

        gs_top = 1.0 - top_in / fig_h
        gs_bottom = bot_in / fig_h

        if self.axes is None:
            group_gap = gap_in / row_h if n_grp > 1 else 0.0
            outer = fig.add_gridspec(
                n_grp,
                1,
                hspace=group_gap,
                left=gs_left,
                right=gs_right,
                top=gs_top,
                bottom=gs_bottom,
            )
            inners = [
                outer[gi].subgridspec(
                    2,
                    n_col,
                    height_ratios=[ratio, 1.0],
                    hspace=0.07,
                    wspace=0.09,
                )
                for gi in range(n_grp)
            ]

        title_pad = 22.0 if top_axis else 6.0

        im_r = im_p = None
        share_y = None
        all_axes: list[tuple[Any, np.ndarray]] = []
        for gi, g in enumerate(groups):
            stations = g["stations"]
            n_st = len(stations)
            x_edges = np.arange(n_st + 1) - 0.5
            per_grid = g["period_grid"]
            yc = np.log10(per_grid) if self.log_period else per_grid
            y_edges = _edges(yc)
            axr0 = axp0 = None
            for ci, c in enumerate(comps):
                link_y = share_y if self.share_period else None
                if inners is not None:
                    axr = fig.add_subplot(
                        inners[gi][0, ci], sharey=link_y
                    )
                    axp = fig.add_subplot(
                        inners[gi][1, ci], sharex=axr, sharey=link_y
                    )
                    if share_y is None:
                        share_y = axr
                else:
                    axr = axlist[(2 * gi) * n_col + ci]
                    axp = axlist[(2 * gi + 1) * n_col + ci]
                all_axes += [(axr, y_edges), (axp, y_edges)]
                if ci == 0:
                    axr0, axp0 = axr, axp

                R = g["res"][c]
                if self.log_res:
                    R = np.log10(np.clip(R, 1e-6, None))
                if self.share_scales:
                    rk = dict(vmin=r_lo, vmax=r_hi)
                    pk = dict(vmin=p_lo, vmax=p_hi)
                else:
                    rk = pk = {}
                im_r = axr.pcolormesh(
                    x_edges, y_edges, R, cmap=self.res_cmap,
                    shading="flat", **rk,
                )
                im_p = axp.pcolormesh(
                    x_edges, y_edges, g["phase"][c], cmap=self.phase_cmap,
                    shading="flat", **pk,
                )
                axr.tick_params(labelsize=8)
                axp.tick_params(labelsize=8)
                axr.tick_params(labelbottom=False)

                if gi == 0:
                    axr.set_title(
                        _Z_LABEL.get(c, c), fontsize=12, pad=title_pad
                    )
                if ci == 0:
                    axr.set_ylabel(ylabel, fontsize=9)
                    axp.set_ylabel(ylabel, fontsize=9)
                else:
                    axr.tick_params(labelleft=False)
                    axp.tick_params(labelleft=False)

                if self.grid:
                    for a in (axr, axp):
                        a.grid(
                            True, which="both", axis="y",
                            ls=":", lw=0.4, alpha=0.45,
                        )

                # station axis -- one per group, on the resistivity row
                # (top) or the phase row (bottom); columns align by
                # station index so each line gets its own axis
                if side == "none":
                    axp.tick_params(labelbottom=False, labeltop=False)
                else:
                    target = axr if top_axis else axp
                    other = axp if top_axis else axr
                    self._station_style(ci, n_col, n_st).apply(
                        target,
                        np.arange(n_st, dtype=float),
                        stations,
                        xlim=(-0.5, n_st - 0.5),
                    )
                    other.tick_params(labelbottom=False, labeltop=False)

            if plabels and axr0 is not None:
                axr0.annotate(
                    plabels[gi],
                    xy=(0.0, 1.0),
                    xycoords="axes fraction",
                    xytext=(-34 if n_grp > 1 else -30, 4),
                    textcoords="offset points",
                    fontsize=11,
                    fontweight="bold",
                    va="bottom",
                    ha="left",
                )
            if n_grp > 1 and g["label"] and axr0 is not None:
                y_top = axr0.get_position().y1
                y_bot = axp0.get_position().y0
                fig.text(
                    0.028,
                    0.5 * (y_top + y_bot),
                    g["label"],
                    rotation=90,
                    va="center",
                    ha="center",
                    fontsize=10,
                    fontweight="bold",
                )

        # period axis: log scale (unless the values are already log10),
        # period increasing downward, each group clipped to its own grid
        for a, yed in all_axes:
            if not self.log_period:
                a.set_yscale("log")
            a.set_ylim(yed[-1], yed[0])

        # two shared colourbars on the right -- only when we own the
        # figure; a caller passing *axes* manages its own colour key
        if show_cbar:
            rlabel = (
                r"$\log_{10}\rho_a$  ($\Omega\cdot$m)"
                if self.log_res
                else r"$\rho_a$  ($\Omega\cdot$m)"
            )
            mid = 0.5 * (gs_top + gs_bottom)
            cb_x = gs_right + 0.022
            cax_r = fig.add_axes(
                [cb_x, mid + 0.035, 0.016, gs_top - mid - 0.07]
            )
            cb_r = fig.colorbar(im_r, cax=cax_r)
            cb_r.set_label(rlabel, fontsize=9)
            cb_r.ax.tick_params(labelsize=8)
            cax_p = fig.add_axes(
                [cb_x, gs_bottom + 0.035, 0.016, mid - gs_bottom - 0.07]
            )
            cb_p = fig.colorbar(im_p, cax=cax_p)
            cb_p.set_label(r"$\varphi$  (deg)", fontsize=9)
            cb_p.ax.tick_params(labelsize=8)

        if self.title:
            fig.suptitle(self.title, fontsize=12, y=0.995, va="top")
        return fig


def _edges(centres: np.ndarray) -> np.ndarray:
    """Return ``n+1`` bin edges bracketing ``n`` monotonic centres."""
    c = np.asarray(centres, dtype=float)
    if c.size == 1:
        return np.array([c[0] - 0.5, c[0] + 0.5])
    mid = 0.5 * (c[:-1] + c[1:])
    return np.concatenate(
        [[c[0] - (mid[0] - c[0])], mid, [c[-1] + (c[-1] - mid[-1])]]
    )


def plot_res_phase_pseudosection(
    sites: Any, *, axes: Any = None, **kwargs: Any
):
    """Functional wrapper around :class:`PlotResPhasePseudoSection`.

    Equivalent to
    ``PlotResPhasePseudoSection(sites, axes=axes, **kwargs).plot()``.
    """
    return PlotResPhasePseudoSection(sites, axes=axes, **kwargs).plot()
