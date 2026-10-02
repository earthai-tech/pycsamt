# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
correction_views — Before/After comparison figures for the Correction window.

Every impedance correction is judged the same way: how did the *response*
change?  This module renders that comparison for any pair of datasets
(raw vs. corrected/preview) in three complementary modes, for both 1-D
sounding curves and 2-D pseudosections:

=============  ==========================  ================================
Compare mode   1-D curves                  2-D pseudosection
=============  ==========================  ================================
Before/After   side-by-side columns,       Before row over After row, one
               shared axis limits          shared colour scale per quantity
Overlay        one axes; before faded /    After filled, Before drawn as
               dashed, after solid         contour lines on top
Diff           ρ ratio (after/before, ×)   log₁₀ ρ ratio and Δφ on a
               and Δφ (°) vs period        diverging, zero-centred scale
=============  ==========================  ================================

Apparent resistivity *and* phase are shown by default, so a correction that
(wrongly or rightly) alters phase is visible at a glance — e.g. a clean
static-shift correction shows a flat ρ ratio and Δφ ≡ 0.

Nothing here ever leaves an empty axes behind: when a view cannot be drawn
the functions raise :class:`PlotUnavailable` carrying a human-readable
title/reason/guidance, which the window shows in place of the canvas.

Qt-free on purpose (only numpy + matplotlib) so it is unit-testable without
a QApplication.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

# ── Public vocabulary (also used as combo-box labels by the window) ──────────

COMPARE_MODES = ("Before / After", "Overlay", "Diff")
QUANTITY_CHOICES = {"ρ_a + φ": ("rho", "phi"), "ρ_a": ("rho",), "φ": ("phi",)}
COMPONENT_CHOICES = {"XY": ("xy",), "YX": ("yx",), "XY + YX": ("xy", "yx")}
ALL_STATIONS = "All stations"

_COMP_LABEL = {"xy": "XY", "yx": "YX"}
# Component colours for single-station views (colour-blind-safe blue/orange).
_COMP_COLOR = {"xy": "#1f6fd1", "yx": "#e8590c"}
# Line style per component when colour already encodes the station.
_COMP_LS = {"xy": "-", "yx": ":"}
_BEFORE_GREY = "#8a94a6"


class PlotUnavailable(Exception):
    """Raised when a comparison view cannot be drawn, with the reason why.

    The window shows ``title``/``reason``/``guidance`` in a placeholder card
    instead of a canvas with empty axes.
    """

    def __init__(self, title: str, reason: str = "", guidance: str = "") -> None:
        super().__init__(f"{title}: {reason}" if reason else title)
        self.title = title
        self.reason = reason
        self.guidance = guidance


# ── Response extraction ───────────────────────────────────────────────────────


@dataclass
class StationResponse:
    """ρ_a / φ of one station, periods sorted ascending."""

    T: np.ndarray
    rho: dict = field(default_factory=dict)  # comp -> ρ_a [Ω·m]
    phi: dict = field(default_factory=dict)  # comp -> φ [deg]

    def get(self, quantity: str, comp: str) -> np.ndarray:
        return (self.rho if quantity == "rho" else self.phi)[comp]


def _wrap_deg(a: np.ndarray) -> np.ndarray:
    return (a + 180.0) % 360.0 - 180.0


def extract_responses(sites) -> dict[str, StationResponse | None]:
    """Ordered ``{station: StationResponse | None}`` for any Sites-like input.

    ρ_a uses the field-unit form 0.2·T·|Z|² (Z in mV/km/nT), matching the
    controller's other plots; φ_yx is shifted by +180° into the first
    quadrant so XY and YX phases share one axis, the usual MT convention.
    """
    from pycsamt.emtools._core import _get_z_block, _iter_items, _name

    out: dict[str, StationResponse | None] = {}
    if sites is None:
        return out
    for i, ed in enumerate(_iter_items(sites)):
        name = _name(ed, i)
        if name in out:  # keep duplicate station ids distinguishable
            name = f"{name} ({i})"
        try:
            _, z, freqs = _get_z_block(ed)
        except Exception:
            z, freqs = None, None
        if z is None or freqs is None or np.size(freqs) == 0:
            out[name] = None
            continue
        f = np.asarray(freqs, dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            T = 1.0 / f
        order = np.argsort(T)
        T, zz = T[order], np.asarray(z)[order]
        resp = StationResponse(T=T)
        for comp, (i0, j0), shift in (("xy", (0, 1), 0.0), ("yx", (1, 0), 180.0)):
            zc = zz[:, i0, j0]
            with np.errstate(invalid="ignore", over="ignore"):
                rho = 0.2 * T * np.abs(zc) ** 2
            rho = np.where(np.isfinite(rho) & (rho > 0), rho, np.nan)
            phi = _wrap_deg(np.degrees(np.angle(zc)) + shift)
            phi = np.where(np.isfinite(zc), phi, np.nan)
            resp.rho[comp] = rho
            resp.phi[comp] = phi
        out[name] = resp
    return out


def _has_values(resp: StationResponse | None, quantities, comps) -> bool:
    if resp is None:
        return False
    return any(
        np.isfinite(resp.get(q, c)).sum() >= 2 for q in quantities for c in comps
    )


def _aligned(before: StationResponse, after: StationResponse, quantity, comp):
    """Return (T, before_vals, after_vals) on the *before* period grid."""
    tb, vb = before.T, before.get(quantity, comp)
    ta, va = after.T, after.get(quantity, comp)
    if tb.shape == ta.shape and np.allclose(tb, ta, rtol=1e-6, equal_nan=True):
        return tb, vb, va
    ok = np.isfinite(ta) & (ta > 0) & np.isfinite(va)
    if ok.sum() < 2:
        return tb, vb, np.full_like(vb, np.nan)
    x = np.log10(ta[ok])
    y = np.log10(va[ok]) if quantity == "rho" else va[ok]
    yi = np.interp(np.log10(tb), x, y, left=np.nan, right=np.nan)
    return tb, vb, (10**yi if quantity == "rho" else yi)


def _same_data(b: dict, a: dict) -> bool:
    if b.keys() != a.keys():
        return False
    for k in b:
        rb, ra = b[k], a[k]
        if (rb is None) != (ra is None):
            return False
        if rb is None:
            continue
        if rb.T.shape != ra.T.shape:
            return False
        for c in ("xy", "yx"):
            if not np.allclose(rb.rho[c], ra.rho[c], equal_nan=True, rtol=1e-9):
                return False
            if not np.allclose(rb.phi[c], ra.phi[c], equal_nan=True, atol=1e-9):
                return False
    return True


# ── Shared axis cosmetics ─────────────────────────────────────────────────────


def _style_axes(ax, s: dict) -> None:
    ax.set_facecolor(s["bg"])
    ax.tick_params(which="both", colors=s["tick"], labelsize=7)
    for sp in ax.spines.values():
        sp.set_edgecolor(s["spine"])


def _grid(ax, s: dict) -> None:
    ax.grid(True, which="major", color=s["grid"], alpha=0.6, ls="-", lw=0.5)
    ax.grid(True, which="minor", color=s["grid"], alpha=0.3, ls=":", lw=0.4)


def _qlabel(quantity: str) -> str:
    return "ρ_a  [Ω·m]" if quantity == "rho" else "φ  [°]"


def _stations_for(before: dict, after: dict) -> list[str]:
    names = list(before.keys())
    names += [n for n in after.keys() if n not in before]
    return names


def _resolve_station(station, names: list[str]) -> str | None:
    if not station or station == ALL_STATIONS:
        return None
    if station not in names:
        raise PlotUnavailable(
            f"Station “{station}” is not in this dataset",
            "The selected station no longer exists after the current "
            "correction (or belongs to a previously loaded survey).",
            "Pick another station, or choose “All stations”.",
        )
    return station


# ── 1-D sounding curves ───────────────────────────────────────────────────────


def _curve_style(state: str, comp: str, single: bool, station_color):
    """Matplotlib kwargs encoding (state, component, station) unambiguously.

    single-station: colour = component; before = dashed + open markers +
    faded, after = solid + filled markers.
    all-stations:   colour = station (after) / neutral grey (before);
    line style = component (XY solid, YX dotted).
    """
    # "before" sits above "after" (zorder) so where a correction left the
    # data unchanged the dashed/grey line stays visible on top of the solid.
    if single:
        col = _COMP_COLOR[comp]
        if state == "before":
            return dict(color=col, ls="--", lw=1.0, alpha=0.65, marker="o",
                        ms=3.4, mfc="white", mew=0.9, zorder=3)
        return dict(color=col, ls="-", lw=1.5, alpha=0.95, marker="o",
                    ms=3.2, mew=0, zorder=2)
    if state == "before":
        return dict(color=_BEFORE_GREY, ls=_COMP_LS[comp], lw=0.8, alpha=0.55,
                    zorder=3)
    return dict(color=station_color, ls=_COMP_LS[comp], lw=1.0, alpha=0.85,
                zorder=2)


def _station_colors(names):
    import matplotlib.pyplot as plt

    cmap = plt.get_cmap("turbo")
    n = max(len(names) - 1, 1)
    return {nm: cmap(0.05 + 0.9 * i / n) for i, nm in enumerate(names)}


def _station_key(fig, axes, names, colors):
    """Discrete colour bar mapping curve colour -> station name."""
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap

    n = len(names)
    sm = ScalarMappable(norm=BoundaryNorm(np.arange(n + 1) - 0.5, n),
                        cmap=ListedColormap([colors[nm] for nm in names]))
    cb = fig.colorbar(sm, ax=axes, pad=0.01, fraction=0.03, aspect=40)
    step = max(1, int(np.ceil(n / 16)))
    ticks = list(range(0, n, step))
    cb.set_ticks(ticks)
    cb.set_ticklabels([names[i] for i in ticks])
    cb.ax.tick_params(labelsize=6)
    cb.set_label("Station", fontsize=7)
    return cb


def _plot_state_curves(ax, data, names, quantity, comps, state, single,
                       colors, force_state_style=None):
    """Draw one dataset's curves into *ax*; return number of curves drawn."""
    drawn = 0
    for nm in names:
        resp = data.get(nm)
        if resp is None:
            continue
        for c in comps:
            v = resp.get(quantity, c)
            ok = np.isfinite(v) & np.isfinite(resp.T) & (resp.T > 0)
            if ok.sum() < 2:
                continue
            style = _curve_style(force_state_style or state, c, single,
                                 colors.get(nm))
            ax.plot(resp.T[ok], v[ok], **style)
            drawn += 1
    return drawn


def _finish_curve_axes(ax, quantity, s, *, show_xlabel):
    ax.set_xscale("log")
    if quantity == "rho":
        ax.set_yscale("log")
    ax.set_ylabel(_qlabel(quantity), fontsize=8, color=s["fg"])
    if show_xlabel:
        ax.set_xlabel("Period  [s]", fontsize=8, color=s["fg"])
    else:
        ax.tick_params(labelbottom=False)
    _grid(ax, s)
    _style_axes(ax, s)


def _legend_curves(ax, comps, single, mode, s, n_stations):
    from matplotlib.lines import Line2D

    handles, labels = [], []
    if single:
        for c in comps:
            if mode in ("Overlay",):
                for st in ("before", "after"):
                    handles.append(Line2D([], [], **_curve_style(st, c, True, None)))
                    labels.append(f"{_COMP_LABEL[c]} {st}")
            else:
                handles.append(Line2D([], [], **_curve_style("after", c, True, None)))
                labels.append(_COMP_LABEL[c])
    else:
        if mode == "Overlay":
            handles.append(Line2D([], [], color=_BEFORE_GREY, lw=1.2, alpha=0.7))
            labels.append("Before (grey)")
            handles.append(Line2D([], [], color="#1f9e89", lw=1.4))
            labels.append(f"After (colour = station, n={n_stations})")
        if len(comps) > 1:
            for c in comps:
                handles.append(Line2D([], [], color=s["fg"], ls=_COMP_LS[c], lw=1.0))
                labels.append(_COMP_LABEL[c])
    if handles:
        ax.legend(handles, labels, fontsize=7, loc="best", framealpha=0.85,
                  ncol=2 if len(handles) > 3 else 1)


def render_curves(
    fig,
    before,
    after,
    *,
    mode: str = "Before / After",
    station: str | None = None,
    quantities=("rho", "phi"),
    components=("xy",),
    theme: dict,
    after_title: str = "After",
    before_title: str = "Before (raw)",
) -> None:
    """ρ_a / φ sounding curves comparing *before* and *after* in *fig*.

    Raises :class:`PlotUnavailable` if there is nothing drawable.
    """
    s = theme
    fig.clear()
    fig.patch.set_facecolor(s["fig_bg"])
    b, a = extract_responses(before), extract_responses(after)
    if not b:
        raise PlotUnavailable(
            "No data loaded",
            "There are no stations to plot.",
            "Load survey data from the main window.",
        )
    names = _stations_for(b, a)
    sel = _resolve_station(station, names)
    names = [sel] if sel else names
    single = sel is not None
    quantities, comps = tuple(quantities), tuple(components)

    if not any(_has_values(b.get(n), quantities, comps) for n in names):
        who = f"station {sel}" if single else "any station"
        raise PlotUnavailable(
            f"No valid {'/'.join(_COMP_LABEL[c] for c in comps)} impedance "
            f"for {who}",
            "The selected impedance component has no finite values, so no "
            "apparent resistivity or phase can be computed.",
            "Try the other component (XY / YX) or another station.",
        )
    colors = _station_colors(names)
    heights = [2 if q == "rho" else 1.3 for q in quantities]

    if mode == "Diff":
        _render_curves_diff(fig, b, a, names, quantities, comps, single,
                            colors, heights, s, sel)
        return

    ncols = 2 if mode == "Before / After" else 1
    gs = fig.add_gridspec(len(quantities), ncols, height_ratios=heights,
                          hspace=0.08, wspace=0.08)
    axes = np.empty((len(quantities), ncols), dtype=object)
    for r in range(len(quantities)):
        for c in range(ncols):
            share_x = axes[0, 0] if (r or c) else None
            share_y = axes[r, 0] if c else None
            axes[r, c] = fig.add_subplot(gs[r, c], sharex=share_x, sharey=share_y)

    for r, q in enumerate(quantities):
        last = r == len(quantities) - 1
        if mode == "Before / After":
            for c, (data, st) in enumerate(((b, "before"), (a, "after"))):
                ax = axes[r, c]
                # Both columns are drawn in the solid "after" style: each
                # column holds one state, so fading one would only mislead.
                n = _plot_state_curves(ax, data, names, q, comps, st, single,
                                       colors, force_state_style="after")
                _finish_curve_axes(ax, q, s, show_xlabel=last)
                if c:
                    ax.tick_params(labelleft=False)
                    ax.set_ylabel("")
                if n == 0:
                    ax.text(0.5, 0.5, "no valid values in this state",
                            transform=ax.transAxes, ha="center", va="center",
                            fontsize=8, color=s["muted"])
        else:  # Overlay
            ax = axes[r, 0]
            _plot_state_curves(ax, b, names, q, comps, "before", single, colors)
            _plot_state_curves(ax, a, names, q, comps, "after", single, colors)
            _finish_curve_axes(ax, q, s, show_xlabel=last)

    who = sel if single else f"{len(names)} stations"
    if mode == "Before / After":
        axes[0, 0].set_title(f"{before_title}  —  {who}", fontsize=9,
                             color=s["title"])
        axes[0, 1].set_title(f"{after_title}  —  {who}", fontsize=9,
                             color=s["title"])
        _legend_curves(axes[0, 1], comps, single, mode, s, len(names))
    else:
        axes[0, 0].set_title(f"Overlay: {before_title} vs {after_title}  —  {who}",
                             fontsize=9, color=s["title"])
        _legend_curves(axes[0, 0], comps, single, mode, s, len(names))
    if not single:
        _station_key(fig, axes.ravel().tolist(), names, colors)


def _render_curves_diff(fig, b, a, names, quantities, comps, single, colors,
                        heights, s, sel):
    if _same_data(b, a):
        raise PlotUnavailable(
            "Nothing to compare yet",
            "No correction has been previewed or applied, so the corrected "
            "data are identical to the raw data.",
            "Choose a correction and click Preview or Apply.",
        )
    gs = fig.add_gridspec(len(quantities), 1, height_ratios=heights, hspace=0.08)
    ax0 = None
    any_drawn = False
    for r, q in enumerate(quantities):
        ax = fig.add_subplot(gs[r, 0], sharex=ax0)
        ax0 = ax0 or ax
        drawn = 0
        stats = []
        for nm in names:
            rb, ra = b.get(nm), a.get(nm)
            if rb is None or ra is None:
                continue
            for c in comps:
                T, vb, va = _aligned(rb, ra, q, c)
                with np.errstate(divide="ignore", invalid="ignore"):
                    d = va / vb if q == "rho" else va - vb
                ok = np.isfinite(d) & np.isfinite(T) & (T > 0)
                if q == "rho":
                    ok &= d > 0
                if ok.sum() < 2:
                    continue
                style = _curve_style("after", c, single, colors.get(nm))
                ax.plot(T[ok], d[ok], **style)
                stats.append(np.nanmedian(d[ok]))
                drawn += 1
        any_drawn |= drawn > 0
        ax.axhline(1.0 if q == "rho" else 0.0, color=s["fg"], lw=0.9,
                   ls="--", alpha=0.6, zorder=0)
        ax.set_xscale("log")
        if q == "rho":
            ax.set_yscale("log")
            ax.set_ylabel("ρ_after / ρ_before  [×]", fontsize=8, color=s["fg"])
            _symmetric_log_ylim(ax)
        else:
            ax.set_ylabel("φ_after − φ_before  [°]", fontsize=8, color=s["fg"])
            lo, hi = ax.get_ylim()
            m = max(abs(lo), abs(hi), 2.0)
            ax.set_ylim(-m, m)
        if stats:
            med = float(np.nanmedian(stats))
            txt = (f"median ×{med:.3g}" if q == "rho" else f"median {med:+.2f}°")
            unchanged = (abs(np.log10(med)) < 1e-4) if q == "rho" else abs(med) < 1e-3
            if unchanged and _all_close_to_ref(ax, q):
                txt = "unchanged by this correction"
            ax.text(0.99, 0.95, txt, transform=ax.transAxes, ha="right",
                    va="top", fontsize=7.5, color=s["fg"],
                    bbox=dict(boxstyle="round,pad=0.25", fc="white",
                              ec=s["spine"], alpha=0.85))
        last = r == len(quantities) - 1
        if last:
            ax.set_xlabel("Period  [s]", fontsize=8, color=s["fg"])
        else:
            ax.tick_params(labelbottom=False)
        _grid(ax, s)
        _style_axes(ax, s)
    if not any_drawn:
        raise PlotUnavailable(
            "No overlapping values to difference",
            "Before and after share no station/period where both are finite.",
            "Check that the correction kept the stations and periods, or use "
            "Before / After instead.",
        )
    who = sel if single else f"{len(names)} stations"
    fig.axes[0].set_title(f"Change caused by the correction  —  {who}",
                          fontsize=9, color=s["title"])
    if single or len(comps) > 1:
        _legend_curves(fig.axes[0], comps, single, "Diff", s, len(names))
    if not single:
        _station_key(fig, list(fig.axes), names, colors)


def _symmetric_log_ylim(ax):
    lo, hi = ax.get_ylim()
    m = max(abs(np.log10(max(lo, 1e-12))), abs(np.log10(max(hi, 1e-12))), 0.05)
    ax.set_ylim(10**-m, 10**m)


def _all_close_to_ref(ax, q) -> bool:
    ref = 1.0 if q == "rho" else 0.0
    for ln in ax.get_lines():
        y = np.asarray(ln.get_ydata(), dtype=float)
        if y.size > 2:  # skip the reference line itself
            if q == "rho":
                if np.nanmax(np.abs(np.log10(y))) > 1e-4:
                    return False
            elif np.nanmax(np.abs(y - ref)) > 1e-3:
                return False
    return True


# ── 2-D pseudosections ────────────────────────────────────────────────────────


def _period_grid(datas, names, n=56):
    Ts = [d[nm].T for d in datas for nm in names if d.get(nm) is not None]
    Ts = [t[np.isfinite(t) & (t > 0)] for t in Ts]
    Ts = [t for t in Ts if t.size]
    if not Ts:
        return None
    all_T = np.concatenate(Ts)
    lo, hi = np.nanpercentile(all_T, 1), np.nanpercentile(all_T, 99)
    if not (np.isfinite(lo) and np.isfinite(hi) and hi > lo):
        return None
    return np.logspace(np.log10(lo), np.log10(hi), n)


def _section_grid(data, names, T_grid, quantity, comp):
    """(n_T, n_station) grid: log10 ρ_a or φ, NaN where undefined."""
    g = np.full((T_grid.size, len(names)), np.nan)
    for j, nm in enumerate(names):
        r = data.get(nm)
        if r is None:
            continue
        v = r.get(quantity, comp)
        ok = np.isfinite(v) & np.isfinite(r.T) & (r.T > 0)
        if quantity == "rho":
            ok &= v > 0
        if ok.sum() < 2:
            continue
        y = np.log10(v[ok]) if quantity == "rho" else v[ok]
        g[:, j] = np.interp(np.log10(T_grid), np.log10(r.T[ok]), y,
                            left=np.nan, right=np.nan)
    return g


def _edges(T_grid, n_st):
    lT = np.log10(T_grid)
    d = np.diff(lT)
    T_edges = 10 ** np.concatenate([[lT[0] - d[0] / 2], (lT[:-1] + lT[1:]) / 2,
                                    [lT[-1] + d[-1] / 2]])
    return np.arange(n_st + 1) - 0.5, T_edges


def _decorate_section(ax, names, T_edges, s, *, sel, affected, show_names,
                      show_ylabel):
    ax.set_yscale("log")
    ax.set_ylim(T_edges[-1], T_edges[0])  # short periods (shallow) on top
    ax.set_xlim(-0.5, len(names) - 0.5)
    step = max(1, int(np.ceil(len(names) / 14)))
    idx = list(range(0, len(names), step))
    ax.set_xticks(idx)
    if show_names:
        ax.set_xticklabels([names[i] for i in idx], rotation=45, ha="right",
                           fontsize=6)
        for lbl, i in zip(ax.get_xticklabels(), idx):
            if names[i] == sel:
                lbl.set_fontweight("bold")
    else:
        ax.tick_params(labelbottom=False)
    if show_ylabel:
        ax.set_ylabel("Period  [s]", fontsize=8, color=s["fg"])
    else:
        ax.tick_params(labelleft=False)
    for nm in affected or ():
        if nm in names:
            ax.axvline(names.index(nm), color="#d6336c", lw=0.9, ls="--",
                       alpha=0.75, zorder=4)
    if sel in names:
        j = names.index(sel)
        ax.axvspan(j - 0.5, j + 0.5, fc="none", ec="black", lw=1.4, zorder=6)
    _style_axes(ax, s)


def _sel_suffix(sel):
    return f"  ·  ▭ {sel}" if sel else ""


def _log_factor_formatter():
    from matplotlib.ticker import FuncFormatter

    return FuncFormatter(lambda v, _p: f"×{10**v:.2g}")


def render_section(
    fig,
    before,
    after,
    *,
    mode: str = "Before / After",
    station: str | None = None,
    quantities=("rho",),
    components=("xy",),
    affected_stations=None,
    theme: dict,
    after_title: str = "After",
    before_title: str = "Before (raw)",
) -> None:
    """Period × station pseudosections comparing *before* and *after*.

    One column per (quantity, component) panel.  Raises
    :class:`PlotUnavailable` when nothing can be drawn.
    """
    s = theme
    fig.clear()
    fig.patch.set_facecolor(s["fig_bg"])
    b, a = extract_responses(before), extract_responses(after)
    if not b:
        raise PlotUnavailable("No data loaded", "There are no stations to plot.",
                              "Load survey data from the main window.")
    names = _stations_for(b, a)
    if len(names) < 2:
        raise PlotUnavailable(
            "A pseudosection needs at least two stations",
            f"This dataset has {len(names)} station — a period × station "
            "section cannot be interpolated from a single sounding.",
            "Switch Display to “Curves (1-D)”.",
        )
    sel = _resolve_station(station, names)
    T_grid = _period_grid((b, a), names)
    if T_grid is None:
        raise PlotUnavailable("No valid periods",
                              "No station has usable frequencies.",
                              "Check the loaded EDI files.")
    panels = [(q, c) for q in quantities for c in components]
    X, T_edges = _edges(T_grid, len(names))
    grids_b = {p: _section_grid(b, names, T_grid, *p) for p in panels}
    grids_a = {p: _section_grid(a, names, T_grid, *p) for p in panels}
    if not any(np.isfinite(g).any() for g in grids_b.values()):
        raise PlotUnavailable(
            "No valid impedance for this component",
            "Every station is NaN for the selected quantity/component.",
            "Try the other component (XY / YX).",
        )

    import numpy.ma as ma

    ncol = len(panels)
    if mode == "Before / After":
        gs = fig.add_gridspec(2, ncol, hspace=0.12, wspace=0.10)
        axes = [[fig.add_subplot(gs[r, c]) for c in range(ncol)] for r in range(2)]
        for c, p in enumerate(panels):
            q, comp = p
            both = np.concatenate([grids_b[p].ravel(), grids_a[p].ravel()])
            both = both[np.isfinite(both)]
            vmin, vmax = ((np.percentile(both, 2), np.percentile(both, 98))
                          if both.size else (0, 1))
            cmap = "jet_r" if q == "rho" else "jet"
            im = None
            for r, (g, ttl) in enumerate(((grids_b[p], before_title),
                                          (grids_a[p], after_title))):
                ax = axes[r][c]
                im = ax.pcolormesh(X, T_edges, ma.masked_invalid(g), cmap=cmap,
                                   vmin=vmin, vmax=vmax, shading="flat")
                _decorate_section(ax, names, T_edges, s, sel=sel,
                                  affected=affected_stations, show_names=r == 1,
                                  show_ylabel=c == 0)
                ax.set_title(f"{ttl}  —  {_panel_name(q, comp)}{_sel_suffix(sel)}",
                             fontsize=8.5,
                             color=s["title"])
            cb = fig.colorbar(im, ax=[axes[0][c], axes[1][c]], pad=0.015,
                              fraction=0.04, aspect=30)
            cb.set_label(_cbar_label(q), fontsize=7)
            cb.ax.tick_params(labelsize=6)
        return

    if mode == "Overlay":
        gs = fig.add_gridspec(1, ncol, wspace=0.10)
        for c, p in enumerate(panels):
            q, comp = p
            ax = fig.add_subplot(gs[0, c])
            ga, gb = grids_a[p], grids_b[p]
            both = np.concatenate([gb.ravel(), ga.ravel()])
            both = both[np.isfinite(both)]
            vmin, vmax = ((np.percentile(both, 2), np.percentile(both, 98))
                          if both.size else (0, 1))
            im = ax.pcolormesh(X, T_edges, ma.masked_invalid(ga),
                               cmap="jet_r" if q == "rho" else "jet",
                               vmin=vmin, vmax=vmax, shading="flat", alpha=0.9)
            levels = np.linspace(vmin, vmax, 8)
            if np.isfinite(gb).sum() > 4 and vmax > vmin:
                cs = ax.contour(np.arange(len(names)), T_grid, gb, levels=levels,
                                colors="black", linewidths=0.8, linestyles="--")
                fmt = _ohm_label if q == "rho" else (lambda v: f"{v:.0f}°")
                ax.clabel(cs, cs.levels[::2], fontsize=6, fmt=fmt, inline=True)
            _decorate_section(ax, names, T_edges, s, sel=sel,
                              affected=affected_stations, show_names=True,
                              show_ylabel=c == 0)
            ax.set_title(f"{_panel_name(q, comp)}: {after_title} (colour) · "
                         f"{before_title} (dashed contours){_sel_suffix(sel)}", fontsize=8.5,
                         color=s["title"])
            cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.05, aspect=30)
            cb.set_label(f"After  {_cbar_label(q)}", fontsize=7)
            cb.ax.tick_params(labelsize=6)
        return

    # Diff
    if _same_data(b, a):
        raise PlotUnavailable(
            "Nothing to compare yet",
            "No correction has been previewed or applied, so the corrected "
            "data are identical to the raw data.",
            "Choose a correction and click Preview or Apply.",
        )
    gs = fig.add_gridspec(1, ncol, wspace=0.10)
    for c, p in enumerate(panels):
        q, comp = p
        ax = fig.add_subplot(gs[0, c])
        d = grids_a[p] - grids_b[p]  # log10 ratio for ρ, Δφ for φ
        fin = np.abs(d[np.isfinite(d)])
        floor = 0.02 if q == "rho" else 1.0
        m = max(float(np.percentile(fin, 98)) if fin.size else 0.0, floor)
        im = ax.pcolormesh(X, T_edges, ma.masked_invalid(d), cmap="RdBu_r",
                           vmin=-m, vmax=m, shading="flat")
        _decorate_section(ax, names, T_edges, s, sel=sel,
                          affected=affected_stations, show_names=True,
                          show_ylabel=c == 0)
        ttl = ("ρ_after / ρ_before" if q == "rho" else "φ_after − φ_before")
        ax.set_title(f"{ttl}  —  {_COMP_LABEL[comp]}{_sel_suffix(sel)}", fontsize=8.5,
                     color=s["title"])
        cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.05, aspect=30)
        if q == "rho":
            cb.formatter = _log_factor_formatter()
            cb.update_ticks()
            cb.set_label("resistivity factor (log scale, white = unchanged)",
                         fontsize=7)
        else:
            cb.set_label("Δφ  [°]  (white = unchanged)", fontsize=7)
        cb.ax.tick_params(labelsize=6)
        if fin.size and fin.max() < (1e-4 if q == "rho" else 1e-3):
            ax.text(0.5, 0.5, "unchanged by this correction",
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, bbox=dict(boxstyle="round", fc="white",
                                          ec=s["spine"], alpha=0.9))


def _ohm_label(v):
    r = 10**v
    if r >= 1000:
        return f"{r / 1000:.3g}k"
    return f"{r:.3g}"


def _panel_name(q, comp):
    return f"{'ρ_a' if q == 'rho' else 'φ'} {_COMP_LABEL[comp]}"


def _cbar_label(q):
    return "log₁₀ ρ_a  [Ω·m]" if q == "rho" else "φ  [°]"


# ── Blank-figure safety net for legacy ax-based plotters ──────────────────────


def figure_blank_reason(fig) -> str | None:
    """Return the message text if *fig* holds no data artists, else ``None``.

    Legacy controller plotters report problems by writing a sentence into an
    otherwise empty axes; the window uses this to swap such figures for the
    placeholder card instead of showing bare axes.
    """
    axes = [ax for ax in fig.axes if ax.get_visible()]
    if not axes:
        return "The plotting routine produced no figure."
    # ``ax.patches`` holds only user-added patches (bars, spans), never the
    # axes' own background, so this works for polar (rose) axes too.
    if any(ax.lines or ax.collections or ax.patches or ax.images for ax in axes):
        return None
    msgs = [t.get_text() for ax in axes for t in ax.texts if t.get_text()]
    return "; ".join(dict.fromkeys(msgs)) or "There is nothing to draw for this view."


__all__ = [
    "ALL_STATIONS",
    "COMPARE_MODES",
    "COMPONENT_CHOICES",
    "QUANTITY_CHOICES",
    "PlotUnavailable",
    "StationResponse",
    "extract_responses",
    "figure_blank_reason",
    "render_curves",
    "render_section",
]
