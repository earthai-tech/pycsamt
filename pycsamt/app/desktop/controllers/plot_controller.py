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

import numpy as np

# ──────────────────────────────────────────────────────────────────────────────
# Dark-style helpers (applied to every axes we touch)
# ──────────────────────────────────────────────────────────────────────────────

_DARK: dict = {
    "facecolor": "#181825",
    "labelcolor": "#cdd6f4",
    "title_color": "#cdd6f4",
    "tick_params": {"colors": "#a6adc8", "labelsize": 8},
    "spines_color": "#45475a",
    "grid_color": "#313244",
    "grid_alpha": 0.35,
    "grid_linestyle": "--",
}

_LIGHT: dict = {
    "facecolor": "#ffffff",
    "labelcolor": "#4c4f69",
    "title_color": "#4c4f69",
    "tick_params": {"colors": "#6c6f85", "labelsize": 8},
    "spines_color": "#bcc0cc",
    "grid_color": "#ccd0da",
    "grid_alpha": 0.5,
    "grid_linestyle": "--",
}


def style_axes(ax, dark: bool = True) -> None:
    """Apply dark or light style to a matplotlib Axes in-place."""
    s = _DARK if dark else _LIGHT
    ax.set_facecolor(s["facecolor"])
    ax.title.set_color(s["title_color"])
    ax.xaxis.label.set_color(s["labelcolor"])
    ax.yaxis.label.set_color(s["labelcolor"])
    # Use explicit labelcolor so set_xticklabels() Text objects are also covered
    tc = s["tick_params"]["colors"]
    ax.tick_params(
        axis="both",
        which="both",
        colors=tc,
        labelcolor=tc,
        labelsize=s["tick_params"]["labelsize"],
    )
    for spine in ax.spines.values():
        spine.set_edgecolor(s["spines_color"])
    ax.grid(
        True,
        color=s["grid_color"],
        alpha=s["grid_alpha"],
        linestyle=s["grid_linestyle"],
    )
    fig = ax.get_figure()
    if fig is not None:
        fig.patch.set_facecolor("#1e1e2e" if dark else "#ffffff")


_DARK_TEXT_MAP = {
    "#2166ac": "#89b4fa",  # dark blue -> Catppuccin blue
    "#b2182b": "#f38ba8",  # dark red  -> Catppuccin red
}


def _readable_on_dark(color, fallback: str) -> str:
    """*color* if it already reads on a dark background, else a lighter
    version of the same hue (near-greys become the theme's label colour)."""
    from matplotlib.colors import to_hex, to_rgb

    try:
        hx = to_hex(color).lower()
        r, g, b = to_rgb(color)
    except (ValueError, TypeError):
        return fallback
    if hx in _DARK_TEXT_MAP:
        return _DARK_TEXT_MAP[hx]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    if lum >= 0.45:
        return hx
    if max(r, g, b) - min(r, g, b) < 0.15:  # black / grey text
        return fallback
    mix = 0.6  # towards white, keeping the hue
    return to_hex((r + (1 - r) * mix, g + (1 - g) * mix, b + (1 - b) * mix))


def _style_figure_full(ref_ax, dark: bool) -> None:
    """
    Post-process every axes in the figure, then fix theme-unaware elements
    produced by emtools (hardcoded white bboxes, black reference ellipse, etc.).

    Call AFTER emtools functions + style_axes on the primary axes.
    """
    s = _DARK if dark else _LIGHT
    fig = ref_ax.get_figure()
    if fig is None:
        return

    # ── 1. Style ALL axes (catches colorbar / inset axes) ────────────────
    for ax in fig.axes:
        ax.set_facecolor(s["facecolor"])
        ax.xaxis.label.set_color(s["labelcolor"])
        ax.yaxis.label.set_color(s["labelcolor"])
        tc = s["tick_params"]["colors"]
        ax.tick_params(
            axis="both",
            which="both",
            colors=tc,
            labelcolor=tc,
            labelsize=s["tick_params"]["labelsize"],
        )
        for spine in ax.spines.values():
            spine.set_edgecolor(s["spines_color"])
        try:
            ax.title.set_color(s["title_color"])
        except Exception:
            pass

    # ── 2./3. Annotation text readable on the theme ──────────────────────
    # emtools draws labels on white boxes (fc="white").  On dark the box is
    # repainted dark, so the text must be light: every too-dark colour is
    # lightened keeping its hue (blue stays blue) -- the old two-colour map
    # left the default-black "size reference" / station labels black on
    # near-black.
    if dark:
        bbox_fc, bbox_ec = "#1a1a2e", "#6c7086"
    else:
        bbox_fc, bbox_ec = "white", "none"
    for ax in fig.axes:
        for txt in ax.texts:
            bb = txt.get_bbox_patch()
            if bb is not None:
                bb.set_facecolor(bbox_fc)
                bb.set_edgecolor(bbox_ec)
                bb.set_alpha(0.88)
            if dark:
                txt.set_color(_readable_on_dark(txt.get_color(),
                                                s["labelcolor"]))

    # ── 4. Fix unfilled patches: reference ellipse edgecolor "k" → visible
    ref_ec = "#cccccc" if dark else "#444444"
    for patch in ref_ax.patches:
        if not patch.get_fill():
            patch.set_edgecolor(ref_ec)


def _annotate_empty(ax, msg: str = "No data") -> None:
    """Clear *ax* and show a centred placeholder message.

    Used for "no data yet" / "select a station" / plot-error states.
    Clearing first guarantees a clean placeholder even when called from
    an ``except`` block after a plot function partially drew before
    failing. Hiding ticks and spines avoids a meaningless default 0..1
    grid on an otherwise-empty axes -- style_axes()'s later
    ``ax.grid(True, ...)`` call still runs but draws nothing since there
    are no tick positions left to align gridlines to.
    """
    ax.cla()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(
        0.5,
        0.5,
        msg,
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=11,
        color="#585b70",
    )


def _clear_figure(fig) -> None:
    """Clear a reusable figure without resetting logarithmic limits to zero."""
    for ax in tuple(fig.axes):
        try:
            if ax.get_xscale() != "linear":
                ax.set_xscale("linear")
            if ax.get_yscale() != "linear":
                ax.set_yscale("linear")
        except Exception:
            pass
    fig.clear()


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
                    f"{int(round(v))}"
                    if abs(v - round(v)) < 0.08
                    else f"{v:.2f}"
                )
            else:
                labels.append("")
        if orient == "y":
            cb_ax.set_yticks(ticks)
            cb_ax.set_yticklabels(labels, fontsize=8)
        else:
            cb_ax.set_xticks(ticks)
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
        from pycsamt.api.station import (
            PYCSAMT_STATION_RENDERING,
        )

        xticks = np.asarray(ax.get_xticks(), dtype=float)
        xlim = ax.get_xlim()
        valid = xticks[(xticks >= xlim[0] - 0.5) & (xticks <= xlim[1] + 0.5)]
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


def _thin_pseudosection_xlabels(ax, keep_label: str | None = None) -> None:
    """Thin cluttered per-station x-tick labels on a pseudosection axes.

    et.pseudosection() labels every station unconditionally, which turns
    into unreadable overlapping text once a profile has more than a
    handful of stations. This reuses the same density-aware
    PYCSAMT_STATION_RENDERING mechanism already applied to the Phase
    Tensor and Res/Phase Section tabs: tick *positions* are left
    untouched (so the triangle markers added by _add_station_markers,
    which reads those same positions, still mark every station) — only
    the label *text* at non-selected positions is blanked out, and the
    final station is always kept.

    Must be called AFTER _mark_station(), which locates the highlighted
    station by matching its rendered tick-label text — blanking that
    text first would make the highlight silently fail to draw.
    *keep_label*, when given, is force-kept visible even where thinning
    would otherwise drop it (keeps the selected station's name on-screen
    alongside its highlight line).
    """
    try:
        from pycsamt.api.station import PYCSAMT_STATION_RENDERING

        labels = [t.get_text() for t in ax.get_xticklabels()]
        if len(labels) <= 1:
            return
        style = PYCSAMT_STATION_RENDERING.style_for("pseudosection")
        figwidth = ax.figure.get_figwidth() if ax.figure is not None else 10.0
        idx = set(int(i) for i in style.label_indices(labels, figwidth))
        if keep_label is not None and keep_label in labels:
            idx.add(labels.index(keep_label))
        new_labels = [lbl if i in idx else "" for i, lbl in enumerate(labels)]
        ax.set_xticklabels(new_labels, rotation=90, fontsize=7, ha="center")
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
                        f"{int(round(v))}"
                        if abs(v - round(v)) < 0.08
                        else f"{v:.1f}"
                    )
                else:
                    new_labels.append(lbl)
            except ValueError:
                new_labels.append(lbl)
        if any(s for s in new_labels):  # only set if we got something
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
        self._sites = None
        self._station_id: str | None = None
        self._period_range: tuple[float, float] | None = None
        self._components: tuple[str, ...] = ("xy", "yx")
        self._phase_ylim: tuple[float, float] | None = None
        self._show_errbar: bool = True
        self._bw_mode: bool = False
        # Shared between the Phase Tensor pseudosection and PT Strip tabs:
        # one toggle controls whether both colour ellipses by signed skew
        # beta (range roughly -10..10 deg, hard to eyeball against the
        # 1-D/2-D vs 3-D structure threshold) or by |beta| (0..max, so the
        # skew_threshold boundary reads directly off the colour scale).
        self._pt_abs_skew: bool = False
        self._pt_annotations: bool = True
        self.dark: bool = True
        # Phase-tensor DataFrame cache: id(sites) → built DataFrame
        # Avoids re-running build_phase_tensor_table() on every tab switch.
        self._pt_df_cache: dict = {}

    # ── State setters ──────────────────────────────────────────────────

    def set_sites(self, sites) -> None:
        self._sites = sites
        self._pt_df_cache.clear()  # stale when new data loaded

    def set_station(self, station_id: str | None) -> None:
        self._station_id = station_id

    def set_period_range(
        self,
        T_min: float | None,
        T_max: float | None,
    ) -> None:
        if T_min is None and T_max is None:
            self._period_range = None
        else:
            self._period_range = (T_min or 0.0, T_max or 1e9)

    def set_components(self, components: tuple[str, ...]) -> None:
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

    def set_pt_annotations(self, show: bool) -> None:
        """Show the in-plot labels of the Phase Tensor (|β| legend, size
        reference) and PT Strip (station name) tabs."""
        self._pt_annotations = bool(show)

    def set_pt_absolute_skew(self, absolute: bool) -> None:
        """Use absolute (|beta|) rather than signed skew colouring.

        Shared by both the Phase Tensor pseudosection and PT Strip tabs.
        """
        self._pt_abs_skew = bool(absolute)

    def set_phase_ylim(
        self,
        ymin: float | None,
        ymax: float | None,
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

    def phase_tensor_key(self) -> tuple:
        """Return a tuple that uniquely identifies the current PT plot settings.

        Panel-level code compares this to the key stored after the last draw
        and skips the full redraw when they match.
        """
        return (
            id(self._sites),
            self._period_range,
            self.dark,
            self._pt_abs_skew,
            self._pt_annotations,  # the labels checkbox must force a redraw
        )

    def invalidate_phase_tensor(self) -> None:
        """Force the next draw_phase_tensor() call to recompute from scratch."""
        self._pt_df_cache.clear()

    def _get_or_build_pt_df(self):
        """Return the cached phase-tensor DataFrame, building it if needed."""
        import pandas as pd

        if self._sites is None:
            return pd.DataFrame()
        key = id(self._sites)
        if key not in self._pt_df_cache:
            from pycsamt.emtools.tensor import (
                build_phase_tensor_table,
            )

            self._pt_df_cache[key] = build_phase_tensor_table(
                self._sites, verbose=0
            )
        return self._pt_df_cache[key]

    def clear(self) -> None:
        self._sites = None
        self._station_id = None
        self._pt_df_cache.clear()

    # ── Period-range filter helpers ────────────────────────────────────

    def _active_period_range(self):
        """Return (T_min, T_max) in seconds if a real filter is active, else None."""
        if self._period_range is None:
            return None
        T_min, T_max = self._period_range
        if T_min <= 1e-9 and T_max >= 1e8:
            return None  # full range → no filter needed
        return (T_min, T_max)

    def _clip_xaxis_period(self, *axes, log: bool = False) -> None:
        """Set X-axis limits to the active period range on all supplied axes.

        *log* must be True when the axis holds pre-transformed
        log10(period) values on a linear scale (the "logperiod"
        convention) rather than raw period on a log-scale axis.
        """
        import math

        pr = self._active_period_range()
        if pr is None:
            return
        T_min, T_max = pr
        if log:
            T_min, T_max = math.log10(T_min), math.log10(T_max)
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
                ax.set_ylim(new_lo, new_hi) if yb <= yt else ax.set_ylim(
                    new_hi, new_lo
                )
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
            _err_phase_deg,
            _err_rhoa,
            _phase_deg,
            _rhoa_from,
            _zblk_flex,
        )

        _clear_figure(fig)
        # Tight 2:1 layout — rho gets the extra space freed by hiding
        # the duplicate x-axis on the top panel.
        gs = fig.add_gridspec(2, 1, height_ratios=[2, 1], hspace=0.0)
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
            T = 1.0 / (fr + 1e-24)

            for comp in self._components:
                color = (
                    "k" if self._bw_mode else self._COMP_COLOR.get(comp, "k")
                )
                ls = self._COMP_LS.get(comp, "-")
                label = f"Z$_{{{comp.upper()}}}$"

                zz = _comp_slice(z, comp)
                ee = None
                if (
                    ze is not None
                    and isinstance(ze, np.ndarray)
                    and ze.shape == z.shape
                ):
                    ee = _comp_slice(ze, comp)

                rho = _rhoa_from(zz, fr)
                rho_err = _err_rhoa(zz, ee, fr) if self._show_errbar else None
                ph = _phase_deg(zz)
                ph_err = _err_phase_deg(zz, ee) if self._show_errbar else None

                kw = dict(
                    color=color,
                    ls=ls,
                    marker=".",
                    ms=5,
                    mfc="white",
                    mec=color,
                    mew=0.9,
                    elinewidth=0.65,
                    capsize=2.5,
                    lw=0.9,
                )
                ax_r.errorbar(T, rho, yerr=rho_err, label=label, **kw)
                ax_p.errorbar(T, ph, yerr=ph_err, label=label, **kw)

            ax_r.set_xscale("log")
            ax_r.set_yscale("log")
            ax_r.set_ylabel(r"$\rho_a\ (\Omega\cdot\mathrm{m})$", fontsize=9)
            ax_r.set_title(
                r"$\rho_a$" f" / $\\varphi$ — {self._station_id}",
                fontsize=10,
                pad=3,
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
            ax_r.legend(ncols=n_comps, fontsize=8, loc="best", framealpha=0.65)

            # Apply period/frequency range filter (axes share X via sharex)
            self._clip_xaxis_period(ax_p)

        except Exception as exc:
            _annotate_empty(ax_r, f"Plot error:\n{exc}")
            _annotate_empty(ax_p)

        for ax in (ax_r, ax_p):
            style_axes(ax, self.dark)

        # Explicit margins are stable inside the resizable Qt canvas; calling
        # tight_layout here repeatedly warns when decorations cannot fit.
        fig.subplots_adjust(
            left=0.11, right=0.97, bottom=0.10, top=0.94, hspace=0.08
        )

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
            _fix_psection_axes(
                ax, colorbar_label=r"$\log_{10}(\rho_a)\ (\Omega{\cdot}m)$"
            )
            _add_station_markers(ax, self.dark)
            self._mark_station(ax)
            _thin_pseudosection_xlabels(ax, keep_label=self._station_id)
            ax.set_title(
                r"$\rho_a$ (XY) — Pseudosection",
                fontsize=10,
                loc="left",
                pad=4,
            )
        except Exception as exc:
            _annotate_empty(ax, f"Pseudosection error:\n{exc}")
        style_axes(ax, self.dark)
        _style_figure_full(ax, self.dark)

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
            _thin_pseudosection_xlabels(ax, keep_label=self._station_id)
            ax.set_title(
                r"$\varphi$ (XY) — Pseudosection",
                fontsize=10,
                loc="left",
                pad=4,
            )
        except Exception as exc:
            _annotate_empty(ax, f"Pseudosection error:\n{exc}")
        style_axes(ax, self.dark)
        _style_figure_full(ax, self.dark)

    # ── Tipper ────────────────────────────────────────────────────────

    def draw_tipper(self, ax) -> None:
        """Draw tipper components on *ax*."""
        import pycsamt.emtools as et

        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            selected = self._get_site(self._station_id) if self._station_id else None
            target = selected if selected is not None else self._sites
            et.plot_tipper_components(target, ax=ax, axis="logperiod", verbose=0)
            station_label = self._station_id if selected is not None else None
            ax.set_title(
                f"Tipper{' — ' + station_label if station_label else ''}",
                fontsize=10,
                pad=3,
            )
            self._clip_xaxis_period(ax, log=True)
        except Exception as exc:
            _annotate_empty(ax, f"Tipper error:\n{exc}")
        style_axes(ax, self.dark)
        _style_figure_full(ax, self.dark)

    # ── Phase tensor ──────────────────────────────────────────────────

    def draw_phase_tensor(self, ax) -> None:
        """Draw phase-tensor pseudosection — stations at TOP, short period at TOP.

        The phase-tensor DataFrame is built once and cached keyed by the Sites
        object identity.  Subsequent calls for the same survey (e.g. period
        range or dark-mode changes) reuse the cached DataFrame and skip the
        expensive build_phase_tensor_table() step.

        Dark mode: pass a light edgecolor so ellipses are visible against the
        dark background.  The rest of the contrast fixes are in _style_figure_full.
        """
        import numpy as np

        import pycsamt.emtools as et

        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        try:
            # Use cached DataFrame when available (avoids expensive recomputation)
            pt_df = self._get_or_build_pt_df()

            # In dark mode use a light edge so ellipses are visible on #181825.
            # In light mode keep the default "k" border.
            edge_kw = dict(
                edgecolor="#cccccc" if self.dark else "k",
                linewidth=0.35,
            )

            et.plot_phase_tensor_psection(
                pt_df,  # pre-built DataFrame → skips build step
                ax=ax,
                period_range=self._active_period_range(),
                period_up=False,  # high freq (short T, near-surface) at TOP
                verbose=0,
                c_by=("|beta|" if self._pt_abs_skew else "skew"),
                symmetric_clim=not self._pt_abs_skew,
                annotations=self._pt_annotations,
                **edge_kw,
            )

            # ── Station labels and markers at top ─────────────────────
            # Collect tick positions before repositioning
            xticks = ax.get_xticks()
            xlim = ax.get_xlim()
            # Preserve existing label strings
            st_labels = [t.get_text() for t in ax.get_xticklabels()]

            # Move station name labels to top edge
            ax.xaxis.tick_top()
            ax.xaxis.set_label_position("top")
            # Re-apply labels with top-friendly alignment (vertical, centered)
            ax.set_xticklabels(
                st_labels, rotation=90, ha="center", va="bottom", fontsize=8
            )

            # Station surface pins — hollow ▽ at y=1 (axes-fraction top edge).
            # mfc='none' keeps these unfilled: only inversion uses filled ▼.
            trans = ax.get_xaxis_transform()
            valid = xticks[
                (xticks >= xlim[0] - 0.1) & (xticks <= xlim[1] + 0.1)
            ]
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
                fontsize=10,
                pad=3,
            )

            # Explicitly recolor station tick labels — set_xticklabels()
            # creates Text objects whose color isn't always updated by the
            # later tick_params call when ticks have been moved to top.
            tick_col = (_DARK if self.dark else _LIGHT)["tick_params"][
                "colors"
            ]
            for lbl in ax.get_xticklabels():
                lbl.set_color(tick_col)

            # Reserve headroom for the rotated (90°) station labels above
            # the axes. plot_phase_tensor_psection's colorbar is attached
            # via make_axes_locatable (no subplotspec), so MplCanvas's own
            # fit_to_view() intentionally skips both tight_layout() and its
            # auto label-space reservation for this figure (see that
            # method's "manually positioned colorbars/insets ... retain
            # the layout defined by their plotting function" docstring) —
            # without this, a short/wide docked panel clips the labels
            # clean off the top of the canvas (confirmed by reproducing at
            # a realistic panel aspect ratio: labels present as Text
            # objects, invisible because the default ~12% top margin isn't
            # enough for 90°-rotated multi-character station names).
            #
            # Same root cause pushes the RIGHT edge too: _attach_cbar()
            # appends the colorbar via make_axes_locatable right at the
            # main axes' current right edge, so its own label ("skew",
            # "|beta|", ...) sits right at the figure's right border and
            # gets clipped by the canvas edge on a docked panel unless we
            # explicitly reserve room for it here (fit_to_view() skips its
            # own auto-margin pass for this same locator-colorbar reason).
            ax.figure.subplots_adjust(top=0.78, right=0.88)

        except Exception as exc:
            _annotate_empty(ax, f"Phase tensor error:\n{exc}")

        style_axes(ax, self.dark)
        _style_figure_full(ax, self.dark)

    def draw_phase_tensor_strip(self, ax) -> None:
        """Draw the single-station phase-tensor ellipse strip for the
        currently selected station (see the Station selector at the top
        of the Profile Viewer).

        Reuses the same cached phase-tensor DataFrame as
        :meth:`draw_phase_tensor` (built once per Sites object).
        """
        import pycsamt.emtools as et

        if self._sites is None:
            _annotate_empty(ax, "Load survey data first")
            style_axes(ax, self.dark)
            return
        if self._station_id is None:
            _annotate_empty(ax, "Select a station to view its ellipse strip")
            style_axes(ax, self.dark)
            return
        try:
            pt_df = self._get_or_build_pt_df()
            edge_kw = dict(
                edgecolor="#cccccc" if self.dark else "k",
                linewidth=0.35,
            )
            et.plot_phase_tensor_strip(
                pt_df,  # pre-built DataFrame → skips build step
                station=self._station_id,
                period_range=self._active_period_range(),
                ax=ax,
                verbose=0,
                c_by=("|beta|" if self._pt_abs_skew else "skew"),
                symmetric_clim=not self._pt_abs_skew,
                axis_style="logperiod",
                station_label=self._pt_annotations,
                **edge_kw,
            )
            ax.set_title(
                f"Phase-tensor ellipse strip — {self._station_id}",
                fontsize=10,
                pad=3,
            )

            # Same root cause as draw_phase_tensor()'s margin fix:
            # _attach_cbar() attaches this colorbar via
            # make_axes_locatable (no subplotspec), so MplCanvas's own
            # fit_to_view() skips tight_layout()/auto margin reservation
            # for this figure. On a short docked panel the default
            # margins aren't enough for the x-label ("log10(T) (s)"),
            # y-label ("Phase (deg)") and colorbar label ("Skew beta
            # (deg)") to all fit -- confirmed by reproducing at a
            # realistic panel aspect ratio (9.8x1.8in): the x-label was
            # dropped off the canvas entirely.
            ax.figure.subplots_adjust(bottom=0.32, right=0.88, top=0.82)

        except Exception as exc:
            _annotate_empty(ax, f"Phase tensor strip error:\n{exc}")

        style_axes(ax, self.dark)
        _style_figure_full(ax, self.dark)

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
        from pycsamt.emtools._core import _get_t_block
        from pycsamt.emtools.plot import (
            _comp_slice,
            _err_phase_deg,
            _err_rhoa,
            _phase_deg,
            _rhoa_from,
            _zblk_flex,
        )

        _clear_figure(fig)

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
            T = 1.0 / (fr + 1e-24)

            # ── tipper block (optional) ────────────────────────────────
            t_out = _get_t_block(site, with_errors=True)
            if len(t_out) == 4:
                _, tipper, tfr, terr = t_out
            else:
                _, tipper, tfr = t_out[:3]
            has_tipper = (
                tipper is not None
                and tfr is not None
                and np.asarray(tipper).size > 0
            )
            if has_tipper:
                tipper = np.asarray(tipper, complex)
                Tt = 1.0 / (np.asarray(tfr, float) + 1e-24)

            # ── GridSpec ───────────────────────────────────────────────
            comps = sorted(self._components)  # xx, xy, yx, yy
            n_cols = max(1, len(comps))
            if has_tipper:
                n_rows = 4  # rho, phase, Tx, Ty
                h_ratios = [4, 2, 1, 1]
            else:
                n_rows = 2
                h_ratios = [2, 1]

            gs = fig.add_gridspec(
                n_rows,
                n_cols,
                height_ratios=h_ratios,
                hspace=0.0,
                wspace=0.10,
            )

            ax_rho_row = []
            ax_phase_row = []
            ax_tx_row = []
            ax_ty_row = []

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
                color = (
                    "k" if self._bw_mode else self._COMP_COLOR.get(comp, "k")
                )
                ls = self._COMP_LS.get(comp, "-")

                zz = _comp_slice(z, comp)
                ee = None
                if (
                    ze is not None
                    and isinstance(ze, np.ndarray)
                    and ze.shape == z.shape
                ):
                    ee = _comp_slice(ze, comp)

                rho = _rhoa_from(zz, fr)
                rho_err = _err_rhoa(zz, ee, fr) if self._show_errbar else None
                ph = _phase_deg(zz)
                ph_err = _err_phase_deg(zz, ee) if self._show_errbar else None

                kw = dict(
                    color=color,
                    ls=ls,
                    marker=".",
                    ms=4,
                    mfc="white",
                    mec=color,
                    mew=0.8,
                    elinewidth=0.6,
                    capsize=2,
                    lw=0.9,
                )
                ar.errorbar(T, rho, yerr=rho_err, **kw)
                ap.errorbar(T, ph, yerr=ph_err, **kw)

                ar.set_xscale("log")
                ar.set_yscale("log")
                ap.set_xscale("log")

                # Component label
                ar.set_title(
                    f"$Z_{{{comp.upper()}}}$",
                    fontsize=9,
                    pad=2,
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
                        marker=".",
                        ms=3,
                        mfc="white",
                        lw=0.8,
                        elinewidth=0.5,
                        capsize=1.5,
                    )
                    # Tx real/imag — same column style as impedance
                    tx_vals = tipper[:, 0]
                    atx.errorbar(
                        Tt,
                        np.real(tx_vals),
                        color=color,
                        ls="-",
                        label="real",
                        **tip_kw,
                    )
                    atx.errorbar(
                        Tt,
                        np.imag(tx_vals),
                        color=color,
                        ls="--",
                        label="imag",
                        **tip_kw,
                    )
                    atx.axhline(0, color="0.55", lw=0.5)
                    atx.set_xscale("log")

                    # Ty
                    ty_vals = tipper[:, 1]
                    aty.errorbar(
                        Tt,
                        np.real(ty_vals),
                        color=color,
                        ls="-",
                        **tip_kw,
                    )
                    aty.errorbar(
                        Tt,
                        np.imag(ty_vals),
                        color=color,
                        ls="--",
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

                for ax in [ar, ap] + (
                    [ax_tx_row[j], ax_ty_row[j]] if has_tipper else []
                ):
                    ax.grid(True, which="both", alpha=0.2, lw=0.4)
                    style_axes(ax, self.dark)

            # Apply period/frequency range filter to all component axes
            self._clip_xaxis_period(*ax_rho_row, *ax_phase_row)
            if has_tipper:
                self._clip_xaxis_period(*ax_tx_row, *ax_ty_row)

            # Global title
            fig.suptitle(
                f"Response — {self._station_id}",
                fontsize=10,
                y=1.01,
            )

        except Exception as exc:
            _clear_figure(fig)
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
