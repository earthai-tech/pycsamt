# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
PlotController — drives all scientific panel redraws.

No Qt imports: usable from the Dash web app.  Each draw method
accepts a matplotlib ``Figure`` or ``Axes`` and draws in-place,
delegating to the existing pycsamt.emtools plot functions.

Design contract
---------------
* Methods never raise — they catch exceptions and annotate the axes with
  an error message so the canvas still updates cleanly.
* All period-range filtering is done before handing data to emtools.
* Dark-style helper ``style_axes`` is applied to every axes after drawing.
"""

from __future__ import annotations

from typing import Optional, Tuple

import numpy as np


# ──────────────────────────────────────────────────────────────────────────────
# Dark-style helpers (applied to every axes we touch)
# ──────────────────────────────────────────────────────────────────────────────

_DARK: dict = {
    "facecolor":        "#181825",
    "labelcolor":       "#cdd6f4",
    "title_color":      "#cdd6f4",
    "tick_params":      {"colors": "#a6adc8", "labelsize": 8},
    "spines_color":     "#45475a",
    "grid_color":       "#313244",
    "grid_alpha":       0.35,
    "grid_linestyle":   "--",
}

_LIGHT: dict = {
    "facecolor":        "#eff1f5",
    "labelcolor":       "#4c4f69",
    "title_color":      "#4c4f69",
    "tick_params":      {"colors": "#6c6f85", "labelsize": 8},
    "spines_color":     "#bcc0cc",
    "grid_color":       "#ccd0da",
    "grid_alpha":       0.5,
    "grid_linestyle":   "--",
}


def style_axes(ax, dark: bool = True) -> None:
    """Apply dark or light style to a matplotlib Axes in-place."""
    s = _DARK if dark else _LIGHT
    ax.set_facecolor(s["facecolor"])
    ax.title.set_color(s["title_color"])
    ax.xaxis.label.set_color(s["labelcolor"])
    ax.yaxis.label.set_color(s["labelcolor"])
    ax.tick_params(**s["tick_params"])
    for spine in ax.spines.values():
        spine.set_edgecolor(s["spines_color"])
    ax.grid(True, color=s["grid_color"], alpha=s["grid_alpha"],
            linestyle=s["grid_linestyle"])
    fig = ax.get_figure()
    if fig is not None:
        fig.patch.set_facecolor("#1e1e2e" if dark else "#e6e9ef")


def _annotate_empty(ax, msg: str = "No data") -> None:
    ax.text(
        0.5, 0.5, msg,
        transform=ax.transAxes,
        ha="center", va="center",
        fontsize=11, color="#585b70",
    )


def _relabel_colorbar_log10(cb_ax) -> None:
    """Relabel colorbar tick values as their log₁₀ equivalents.

    The pseudosection data is plotted in linear (Ω·m) scale.  When the
    colorbar label indicates log₁₀ units, tick values like 70 000 should
    display as 4.85, 10 000 as 4.0, etc.  Ticks that are ≤ 0 are
    labelled with an empty string.
    """
    try:
        ticks = cb_ax.get_yticks()
        orient = "y"
        if len(ticks) == 0:
            ticks = cb_ax.get_xticks()
            orient = "x"
        labels = []
        for t in ticks:
            if t > 0:
                v = np.log10(t)
                labels.append(
                    f"{int(round(v))}" if abs(v - round(v)) < 0.08
                    else f"{v:.2f}"
                )
            else:
                labels.append("")
        if orient == "y":
            cb_ax.set_yticklabels(labels, fontsize=8)
        else:
            cb_ax.set_xticklabels(labels, fontsize=8)
    except Exception:
        pass


def _add_station_markers(ax, dark: bool = True) -> None:
    """Add triangular station markers at the top of a pseudosection axes.

    Called after et.pseudosection() which already placed tick labels at the top.
    Uses PYCSAMT_STATION_RENDERING so the marker style matches the package API.
    The dark flag adjusts edgecolor for visibility against the panel background.
    """
    try:
        from pycsamt.api.station import PYCSAMT_STATION_RENDERING
        xticks = np.asarray(ax.get_xticks(), dtype=float)
        xlim   = ax.get_xlim()
        valid  = xticks[(xticks >= xlim[0] - 0.5) & (xticks <= xlim[1] + 0.5)]
        if valid.size == 0:
            return
        m = PYCSAMT_STATION_RENDERING.style_for("pseudosection").marker
        kw = m.kwargs()
        if dark:
            kw["edgecolors"] = "#cccccc"
            kw["facecolors"] = "none"
        ax.scatter(
            valid,
            np.full(valid.size, m.offset),
            transform=ax.get_xaxis_transform(),
            **kw,
        )
    except Exception:
        pass


def _apply_log10_xfmt(ax) -> None:
    """Format a log-scale x-axis to show log₁₀ integer exponents.

    Converts tick positions like 10⁻³, 10⁰, 10² to the plain integers
    -3, 0, 2.  Correct when the axis label reads log₁₀(T) (s).
    Minor-tick labels are suppressed.
    """
    import matplotlib.ticker as mticker

    def _fmt(x, _):
        if x <= 0:
            return ""
        v = np.log10(x)
        iv = round(v)
        return f"${iv:d}$" if abs(v - iv) < 0.02 else f"${v:.1f}$"

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_fmt))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())


def _fix_psection_axes(ax, colorbar_label: str = "") -> None:
    """Post-process a pseudosection axes after et.pseudosection().

    et.pseudosection() uses imshow() whose Y-axis is in *row-index* space
    (0 … n_period-1), not period seconds.  Tick labels are the period values
    set explicitly by set_yticklabels inside et.pseudosection.  Reading those
    stored label strings and converting to log₁₀ notation is the only robust
    approach — do NOT call set_ylim with period values on this axis.
    """
    # ── y-axis: convert stored period tick-label strings to log₁₀(T) ──
    try:
        labels = [t.get_text() for t in ax.get_yticklabels()]
        new_labels = []
        for lbl in labels:
            try:
                val = float(lbl)
                if val > 0:
                    v = np.log10(val)
                    new_labels.append(
                        f"{int(round(v))}" if abs(v - round(v)) < 0.08
                        else f"{v:.1f}"
                    )
                else:
                    new_labels.append(lbl)
            except ValueError:
                new_labels.append(lbl)
        if any(s for s in new_labels):        # only set if we got something
            ax.set_yticklabels(new_labels, fontsize=8)
        ax.set_ylabel(r"$\log_{10}(T)\ \mathrm{(s)}$", fontsize=9)
    except Exception:
        pass

    # ── colorbar label + tick conversion ──────────────────────────────
    if colorbar_label:
        try:
            fig = ax.get_figure()
            caxes = [a for a in fig.axes if a is not ax]
            if caxes:
                cb_ax = caxes[-1]
                cb_ax.set_ylabel(colorbar_label, fontsize=8)
                # When label implies log10 units (e.g. log₁₀(ρₐ)), the data
                # was plotted in linear (Ω·m) scale but we want log10 tick
                # labels so values like 70000 → 4.85, 10000 → 4.0, etc.
                if "log" in colorbar_label.lower():
                    _relabel_colorbar_log10(cb_ax)
        except Exception:
            pass


# ──────────────────────────────────────────────────────────────────────────────
# PlotController
# ──────────────────────────────────────────────────────────────────────────────

class PlotController:
    """
    Pure-Python controller that drives all scientific panel redraws.

    Attributes
    ----------
    dark : bool
        Whether to apply dark-mode styling (default True).
    """

    def __init__(self) -> None:
        self._sites         = None
        self._station_id:   Optional[str]   = None
        self._period_range: Optional[Tuple[float, float]] = None
        self._components:   Tuple[str, ...] = ("xy", "yx")
        self._phase_ylim:   Optional[Tuple[float, float]] = None
        self._show_errbar:  bool = True
        self._bw_mode:      bool = False
        self.dark: bool     = True

    # ── State setters ──────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        self._sites = sites

    def set_station(self, station_id: Optional[str]) -> None:
        self._station_id = station_id

    def set_period_range(
        self,
        T_min: Optional[float],
        T_max: Optional[float],
    ) -> None:
        if T_min is None and T_max is None:
            self._period_range = None
        else:
            self._period_range = (T_min or 0.0, T_max or 1e9)

    def set_components(self, components: Tuple[str, ...]) -> None:
        """Set which impedance components to plot (e.g. ('xy','yx','xx','yy'))."""
        self._components = components if components else ("xy", "yx")

    def set_show_errbar(self, show: bool) -> None:
        """Toggle error bars on ρₐ / φ response curves."""
        self._show_errbar = show

    def set_bw_mode(self, on: bool) -> None:
        """Switch all response curves to black (publication B/W mode).

        Line styles and error bars are preserved; only the colour changes
        to black so the figure is suitable for greyscale publication.
        """
        self._bw_mode = on

    def set_phase_ylim(
        self,
        ymin: Optional[float],
        ymax: Optional[float],
    ) -> None:
        """Force the phase axes y-limits.  Pass (None, None) for auto."""
        if ymin is None and ymax is None:
            self._phase_ylim = None
        else:
            self._phase_ylim = (ymin, ymax)

    @property
    def _x_label(self) -> str:
        """X-axis label string derived from the API x-axis type.

        Currently always period-based.  When the API gains a frequency-mode
        flag this property is the single place to update.
        """
        return r"$\log_{10}(T)\ \mathrm{(s)}$"

    def clear(self) -> None:
        self._sites      = None
        self._station_id = None

    # ── Period-range filter helpers ────────────────────────────────────

    def _active_period_range(self):
        """Return (T_min, T_max) in seconds if a real filter is active, else None."""
        if self._period_range is None:
            return None
        T_min, T_max = self._period_range
        if T_min <= 1e-9 and T_max >= 1e8:
            return None          # full range → no filter needed
        return (T_min, T_max)

    def _clip_xaxis_period(self, *axes) -> None:
        """Set X-axis limits to the active period range on all supplied axes."""
        pr = self._active_period_range()
        if pr is None:
            return
        T_min, T_max = pr
        for ax in axes:
            try:
                ax.set_xlim(T_min, T_max)
            except Exception:
                pass

    def _clip_yaxis_period(self, *axes) -> None:
        """Set Y-axis limits to the active period range, respecting axis direction."""
        pr = self._active_period_range()
        if pr is None:
            return
        T_min, T_max = pr
        for ax in axes:
            try:
                yb, yt = ax.get_ylim()
                new_lo = max(min(yb, yt), T_min)
                new_hi = min(max(yb, yt), T_max)
                # Preserve inverted axes (some pseudosections put long T at bottom)
                ax.set_ylim(new_lo, new_hi) if yb <= yt else ax.set_ylim(new_hi, new_lo)
            except Exception:
                pass

    # ── Paired color / style tables ────────────────────────────────────

    # XY and XX share blue; YX and YY share red.
    # Within each pair the off-diagonal (xy/yx) uses a solid line while
    # the diagonal (xx/yy) uses a dashed line.
    _COMP_COLOR = {
        "xy": "#1f77b4",
        "xx": "#1f77b4",
        "yx": "#d62728",
        "yy": "#d62728",
    }
    _COMP_LS = {
        "xy": "-",
        "xx": "--",
        "yx": "-",
        "yy": "--",
    }

    # ── ρₐ / φ response curves (single station) ───────────────────────

    def draw_rho_phi(self, fig) -> None:
        """Draw apparent-resistivity and phase curves for the selected station.

        Paired colour scheme: xy/xx share blue (solid/dashed), yx/yy share
        red (solid/dashed).  Error bars appear on *both* ρₐ and φ when
        enabled.  Log scales are applied to both axes so tick labels stay
        clean (no ugly exponentials).
        """
        from pycsamt.emtools.plot import (
            _comp_slice,
            _phase_deg,
            _rhoa_from,
            _err_rhoa,
            _err_phase_deg,
            _zblk_flex,
        )

        fig.clear()
        # Tight 2:1 layout — rho gets the extra space freed by hiding
        # the duplicate x-axis on the top panel.
        gs   = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.0)
        ax_r = fig.add_subplot(gs[0])
        ax_p = fig.add_subplot(gs[1], sharex=ax_r)

        if self._sites is None:
            _annotate_empty(ax_r, "Load survey data first")
            _annotate_empty(ax_p)
            for ax in (ax_r, ax_p):
                style_axes(ax, self.dark)
            return

        if self._station_id is None:
            _annotate_empty(ax_r, "Select a station to view response curves")
            _annotate_empty(ax_p)
            for ax in (ax_r, ax_p):
                style_axes(ax, self.dark)
            return

        try:
            site = self._get_site(self._station_id)
            if site is None:
                raise ValueError(f"Station '{self._station_id}' not found")

            out = _zblk_flex(site)
            if len(out) == 4:
                _, z, fr, ze = out
            else:
                _, z, fr = out[:3]
                ze = None

            if z is None or fr is None:
                raise ValueError("No impedance data available")

            fr = np.asarray(fr, dtype=float)
            T  = 1.0 / (fr + 1e-24)

            for comp in self._components:
                color = "k" if self._bw_mode else self._COMP_COLOR.get(comp, "k")
                ls    = self._COMP_LS.get(comp, "-")
                label = f"Z$_{{{comp.upper()}}}$"

                zz = _comp_slice(z, comp)
                ee = None
                if ze is not None and isinstance(ze, np.ndarray) \
                        and ze.shape == z.shape:
                    ee = _comp_slice(ze, comp)

                rho     = _rhoa_from(zz, fr)
                rho_err = _err_rhoa(zz, ee, fr) if self._show_errbar else None
                ph      = _phase_deg(zz)
                ph_err  = _err_phase_deg(zz, ee) if self._show_errbar else None

                kw = dict(
                    color=color, ls=ls,
                    marker=".", ms=5, mfc="white", mec=color, mew=0.9,
                    elinewidth=0.65, capsize=2.5, lw=0.9,
                )
                ax_r.errorbar(T, rho, yerr=rho_err, label=label, **kw)
                ax_p.errorbar(T, ph,  yerr=ph_err,  label=label, **kw)

            ax_r.set_xscale("log")
            ax_r.set_yscale("log")
            ax_r.set_ylabel(
                r"$\rho_a\ (\Omega\cdot\mathrm{m})$", fontsize=9
            )
            ax_r.set_title(
                r"$\rho_a$"
                f" / $\\varphi$ — {self._station_id}",
                fontsize=10, pad=3,
            )
            # x ticks/labels only on the bottom (phase) panel
            ax_r.tick_params(labelbottom=False)

            ax_p.set_xscale("log")
            ax_p.set_ylabel(r"$\varphi\ (°)$", fontsize=9)
            ax_p.set_xlabel(self._x_label, fontsize=9)
            # Tick values match the label: show -3, -2, -1, 0... not 10^{-3}
            _apply_log10_xfmt(ax_p)

            if self._phase_ylim is not None:
                ax_p.set_ylim(self._phase_ylim)

            n_comps = len(self._components)
            ax_r.legend(ncols=n_comps, fontsize=8, loc="best",
                        framealpha=0.65)

            # Apply period/frequency range filter (axes share X via sharex)
            self._clip_xaxis_period(ax_p)

        except Exception as exc:
            _annotate_empty(ax_r, f"Plot error:\n{exc}")
            _annotate_empty(ax_p)

        for ax in (ax_r, ax_p):
            style_axes(ax, self.dark)

        fig.tight_layout(pad=1.0, h_pad=0.4)

    # ── Apparent-resistivity pseudosection ────────────────────────────

    def draw_rho_pseudosection(self, ax) -> None:
        """Draw log(ρₐ_xy) pseudosection on *ax*."""
        import pycsamt.emtools as et
        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            et.pseudosection(
                self._sites,
                quantity="rho_xy",
                ax=ax,
                verbose=0,
                period_range=self._active_period_range(),
            )
            _fix_psection_axes(ax, colorbar_label=r"$\log_{10}(\rho_a)\ (\Omega{\cdot}m)$")
            _add_station_markers(ax, self.dark)
            self._mark_station(ax)
            ax.set_title(r"$\rho_a$ (XY) — Pseudosection", fontsize=10,
                         loc="left", pad=4)
        except Exception as exc:
            _annotate_empty(ax, f"Pseudosection error:\n{exc}")
        style_axes(ax, self.dark)

    # ── Phase pseudosection ────────────────────────────────────────────

    def draw_phase_pseudosection(self, ax) -> None:
        """Draw phase_xy pseudosection on *ax*."""
        import pycsamt.emtools as et
        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            et.pseudosection(
                self._sites,
                quantity="phi_xy",
                ax=ax,
                verbose=0,
                period_range=self._active_period_range(),
            )
            _fix_psection_axes(ax, colorbar_label=r"$\varphi$ (°)")
            _add_station_markers(ax, self.dark)
            self._mark_station(ax)
            ax.set_title(r"$\varphi$ (XY) — Pseudosection", fontsize=10,
                         loc="left", pad=4)
        except Exception as exc:
            _annotate_empty(ax, f"Pseudosection error:\n{exc}")
        style_axes(ax, self.dark)

    # ── Tipper ────────────────────────────────────────────────────────

    def draw_tipper(self, ax) -> None:
        """Draw tipper components on *ax*."""
        import pycsamt.emtools as et
        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            target = (
                self._get_site(self._station_id)
                if self._station_id else self._sites
            )
            et.plot_tipper_components(target, ax=ax, verbose=0)
            ax.set_title(
                f"Tipper{' — ' + self._station_id if self._station_id else ''}",
                fontsize=10,
                pad=3,
            )
            self._clip_xaxis_period(ax)
        except Exception as exc:
            _annotate_empty(ax, f"Tipper error:\n{exc}")
        style_axes(ax, self.dark)

    # ── Phase tensor ──────────────────────────────────────────────────

    def draw_phase_tensor(self, ax) -> None:
        """Draw phase-tensor pseudosection — stations at TOP, short period at TOP."""
        import numpy as np
        import pycsamt.emtools as et
        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            et.plot_phase_tensor_psection(
                self._sites,
                ax=ax,
                period_range=self._active_period_range(),
                period_up=False,   # high freq (short T, near-surface) at TOP
                verbose=0,
            )

            # ── Station labels and markers at top ─────────────────────
            # Collect tick positions before repositioning
            xticks = ax.get_xticks()
            xlim   = ax.get_xlim()
            # Preserve existing label strings
            st_labels = [t.get_text() for t in ax.get_xticklabels()]

            # Move station name labels to top edge
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            # Re-apply labels with top-friendly alignment (vertical, centered)
            ax.set_xticklabels(st_labels, rotation=90, ha="center",
                               va="bottom", fontsize=8)

            # Station surface pins — hollow ▽ at y=1 (axes-fraction top edge).
            # mfc='none' keeps these unfilled: only inversion uses filled ▼.
            trans     = ax.get_xaxis_transform()
            valid     = xticks[(xticks >= xlim[0] - 0.1) & (xticks <= xlim[1] + 0.1)]
            pin_color = "#cccccc" if self.dark else "#444444"
            ax.plot(
                valid,
                np.ones(len(valid)),
                marker="v",
                markerfacecolor="none",
                markeredgecolor=pin_color,
                markeredgewidth=1.2,
                markersize=6,
                linestyle="none",
                transform=trans,
                clip_on=False,
                zorder=5,
            )

            ax.set_title(
                r"Phase Tensor $\beta$ — Pseudosection",
                fontsize=10, pad=3,
            )
        except Exception as exc:
            _annotate_empty(ax, f"Phase tensor error:\n{exc}")
        style_axes(ax, self.dark)

    # ── Publication-quality multi-panel view ──────────────────────────

    def draw_publication_view(self, fig) -> None:
        """Build an MTpy-style publication figure for the selected station.

        Layout (per checked component column):
            ┌─────────────────────┐
            │  App. Res.  (2×)    │
            ├─────────────────────┤
            │  Phase      (1×)    │
            ├─────────────────────┤
            │  T_x / T_y  (0.5×) │  ← only when tipper data present
            └─────────────────────┘
        Zero hspace between rows so the panels are tight.  The figure is
        drawn into *fig* and never creates its own window — the caller is
        responsible for displaying it (e.g. ``PublicationViewDialog``).
        """
        from pycsamt.emtools.plot import (
            _comp_slice, _phase_deg, _rhoa_from,
            _err_rhoa, _err_phase_deg, _zblk_flex,
        )
        from pycsamt.emtools._core import _get_t_block

        fig.clear()

        if self._sites is None or self._station_id is None:
            ax = fig.add_subplot(111)
            _annotate_empty(ax, "Select a station first")
            style_axes(ax, self.dark)
            return

        try:
            site = self._get_site(self._station_id)
            if site is None:
                raise ValueError(f"Station '{self._station_id}' not found")

            # ── impedance block ────────────────────────────────────────
            out = _zblk_flex(site)
            if len(out) == 4:
                _, z, fr, ze = out
            else:
                _, z, fr = out[:3]
                ze = None
            if z is None or fr is None:
                raise ValueError("No impedance data")
            fr = np.asarray(fr, dtype=float)
            T  = 1.0 / (fr + 1e-24)

            # ── tipper block (optional) ────────────────────────────────
            t_out = _get_t_block(site, with_errors=True)
            if len(t_out) == 4:
                _, tipper, tfr, terr = t_out
            else:
                _, tipper, tfr = t_out[:3]
                terr = None
            has_tipper = (
                tipper is not None and tfr is not None
                and np.asarray(tipper).size > 0
            )
            if has_tipper:
                tipper = np.asarray(tipper, complex)
                Tt = 1.0 / (np.asarray(tfr, float) + 1e-24)

            # ── GridSpec ───────────────────────────────────────────────
            comps   = sorted(self._components)   # xx, xy, yx, yy
            n_cols  = max(1, len(comps))
            if has_tipper:
                n_rows      = 4          # rho, phase, Tx, Ty
                h_ratios    = [4, 2, 1, 1]
            else:
                n_rows      = 2
                h_ratios    = [2, 1]

            gs = fig.add_gridspec(
                n_rows, n_cols,
                height_ratios=h_ratios,
                hspace=0.0,
                wspace=0.10,
            )

            ax_rho_row   = []
            ax_phase_row = []
            ax_tx_row    = []
            ax_ty_row    = []

            for j in range(n_cols):
                ar = fig.add_subplot(gs[0, j])
                ap = fig.add_subplot(gs[1, j], sharex=ar)
                ax_rho_row.append(ar)
                ax_phase_row.append(ap)
                if has_tipper:
                    atx = fig.add_subplot(gs[2, j], sharex=ar)
                    aty = fig.add_subplot(gs[3, j], sharex=ar)
                    ax_tx_row.append(atx)
                    ax_ty_row.append(aty)

            # ── draw per component column ──────────────────────────────
            for j, comp in enumerate(comps):
                ar = ax_rho_row[j]
                ap = ax_phase_row[j]
                color = "k" if self._bw_mode else self._COMP_COLOR.get(comp, "k")
                ls    = self._COMP_LS.get(comp, "-")

                zz = _comp_slice(z, comp)
                ee = None
                if ze is not None and isinstance(ze, np.ndarray) \
                        and ze.shape == z.shape:
                    ee = _comp_slice(ze, comp)

                rho     = _rhoa_from(zz, fr)
                rho_err = _err_rhoa(zz, ee, fr) if self._show_errbar else None
                ph      = _phase_deg(zz)
                ph_err  = _err_phase_deg(zz, ee) if self._show_errbar else None

                kw = dict(
                    color=color, ls=ls,
                    marker=".", ms=4, mfc="white", mec=color, mew=0.8,
                    elinewidth=0.6, capsize=2, lw=0.9,
                )
                ar.errorbar(T, rho, yerr=rho_err, **kw)
                ap.errorbar(T, ph,  yerr=ph_err,  **kw)

                ar.set_xscale("log")
                ar.set_yscale("log")
                ap.set_xscale("log")

                # Component label
                ar.set_title(
                    f"$Z_{{{comp.upper()}}}$",
                    fontsize=9, pad=2,
                )

                # Y labels only on leftmost column
                if j == 0:
                    ar.set_ylabel(
                        r"$\rho_a\ (\Omega\cdot\mathrm{m})$",
                        fontsize=8,
                    )
                    ap.set_ylabel(r"$\varphi\ (°)$", fontsize=8)
                else:
                    ar.set_yticklabels([])
                    ap.set_yticklabels([])

                # X tick labels only on bottom row
                if has_tipper:
                    ar.tick_params(labelbottom=False)
                    ap.tick_params(labelbottom=False)
                else:
                    ar.tick_params(labelbottom=False)
                    ap.tick_params(labelbottom=True, labelsize=7)
                    ap.set_xlabel(self._x_label, fontsize=8)
                    _apply_log10_xfmt(ap)

                if self._phase_ylim is not None:
                    ap.set_ylim(self._phase_ylim)

                # Tipper sub-panels
                if has_tipper:
                    atx = ax_tx_row[j]
                    aty = ax_ty_row[j]
                    tip_kw = dict(
                        marker=".", ms=3, mfc="white", lw=0.8,
                        elinewidth=0.5, capsize=1.5,
                    )
                    # Tx real/imag — same column style as impedance
                    tx_vals = tipper[:, 0]
                    atx.errorbar(
                        Tt, np.real(tx_vals), color=color, ls="-",
                        label="real", **tip_kw,
                    )
                    atx.errorbar(
                        Tt, np.imag(tx_vals), color=color, ls="--",
                        label="imag", **tip_kw,
                    )
                    atx.axhline(0, color="0.55", lw=0.5)
                    atx.set_xscale("log")

                    # Ty
                    ty_vals = tipper[:, 1]
                    aty.errorbar(
                        Tt, np.real(ty_vals), color=color, ls="-",
                        **tip_kw,
                    )
                    aty.errorbar(
                        Tt, np.imag(ty_vals), color=color, ls="--",
                        **tip_kw,
                    )
                    aty.axhline(0, color="0.55", lw=0.5)
                    aty.set_xscale("log")
                    aty.set_xlabel(self._x_label, fontsize=8)
                    _apply_log10_xfmt(aty)

                    if j == 0:
                        atx.set_ylabel(r"$T_x\ (-)$", fontsize=8)
                        aty.set_ylabel(r"$T_y\ (-)$", fontsize=8)
                    else:
                        atx.set_yticklabels([])
                        aty.set_yticklabels([])
                    atx.tick_params(labelbottom=False, labelsize=7)
                    aty.tick_params(labelbottom=True, labelsize=7)

                for ax in ([ar, ap]
                           + ([ax_tx_row[j], ax_ty_row[j]] if has_tipper else [])):
                    ax.grid(True, which="both", alpha=0.2, lw=0.4)
                    style_axes(ax, self.dark)

            # Apply period/frequency range filter to all component axes
            self._clip_xaxis_period(*ax_rho_row, *ax_phase_row)
            if has_tipper:
                self._clip_xaxis_period(*ax_tx_row, *ax_ty_row)

            # Global title
            fig.suptitle(
                f"Response — {self._station_id}",
                fontsize=10, y=1.01,
            )

        except Exception as exc:
            fig.clear()
            ax = fig.add_subplot(111)
            _annotate_empty(ax, f"Publication view error:\n{exc}")
            style_axes(ax, self.dark)

    # ── Helpers ───────────────────────────────────────────────────────

    def _get_site(self, station_id: str):
        """Return the Site object for *station_id*, or None."""
        if self._sites is None:
            return None
        try:
            return self._sites.get(station_id)
        except Exception:
            pass
        try:
            # as_list() yields EDIFile objects; use .station (not .name)
            for edi in self._sites.as_list():
                if edi.station == station_id:
                    return self._sites.get(edi.station)
        except Exception:
            pass
        return None

    def _mark_station(self, ax) -> None:
        """Add a vertical dashed marker at the selected station."""
        if self._station_id is None:
            return
        try:
            # Try to find the x-tick position of the selected station
            labels = [t.get_text() for t in ax.get_xticklabels()]
            if self._station_id in labels:
                idx = labels.index(self._station_id)
                ticks = ax.get_xticks()
                if idx < len(ticks):
                    ax.axvline(
                        ticks[idx],
                        color="#f5c542",
                        linestyle="--",
                        linewidth=1.0,
                        alpha=0.8,
                        label=self._station_id,
                    )
        except Exception:
            pass
