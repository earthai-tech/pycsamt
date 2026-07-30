# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
callbacks/correction.py — Data Correction page (v2 redesign).

Architecture
------------
No module-level controller singleton.  Every callback builds a fresh
CorrectionController from the session cache, replays the stored step list
(CORR_STORE), then executes the requested operation.

CORR_STORE holds {"steps": [{"fn_name": ..., "kwargs": {...}, "label": ...}]}
which is JSON-serialisable and lives in a client-side dcc.Store.

Auto-render chain
-----------------
  apply / undo / reset / drop-freq  →  write CORR_STORE
                                     →  auto_render fires (Input CORR_STORE)
                                     →  active panel(s) re-rendered

Tab system
----------
  CORR_ACTIVE_TAB (dcc.Store) → clientside callbacks toggle panel display
  Panel IDs:  corr-panel-{ba|rho-phi|overlay|diff}
  Button IDs: corr-tab-btn-{ba|rho-phi|overlay|diff}

Before/After distinction
------------------------
  auto_render always shows: left = raw_sites, right = current_sites
  preview_correction shows: left = current_sites (before step), right = preview result
  Labels are updated dynamically to make the distinction clear.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import dash_bootstrap_components as dbc
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from dash import (
    ALL,
    Input,
    Output,
    State,
    ctx,
    dcc,
    html,
    no_update,
)
from dash.exceptions import PreventUpdate

from pycsamt.app.desktop.controllers.correction_controller import (
    CATALOGUE,
    CorrectionController,
)
from pycsamt.app.web.cache import cache_get
from pycsamt.app.web.layout import IDs
from pycsamt.app.web.utils import (
    apply_web_dark_theme,
    apply_web_light_theme,
    empty_src,
    fig_to_src,
    filter_sites_by_lines,
    no_active_lines_src,
)

# ── Catalogue helpers ─────────────────────────────────────────────────────────


def _get_specs(cat: str, method: str) -> list:
    return CATALOGUE.get(cat or "", {}).get(method or "", {}).get("params", [])


def _get_fn_name(cat: str, method: str) -> str | None:
    return CATALOGUE.get(cat or "", {}).get(method or "", {}).get("fn")


def _get_desc(cat: str, method: str) -> str:
    return CATALOGUE.get(cat or "", {}).get(method or "", {}).get("desc", "")


def _build_kwargs(cat: str, method: str, param_vals: list) -> dict:
    """Map ordered param_vals (from ALL State) back to {name: value} kwargs."""
    specs = _get_specs(cat, method)
    kwargs = {}
    for spec, val in zip(specs, param_vals or []):
        if val is None:
            kwargs[spec.name] = spec.default
        elif spec.kind == "spin":
            try:
                kwargs[spec.name] = int(float(val))
            except (ValueError, TypeError):
                kwargs[spec.name] = spec.default
        elif spec.kind == "dspin":
            try:
                kwargs[spec.name] = float(val)
            except (ValueError, TypeError):
                kwargs[spec.name] = spec.default
        elif spec.kind == "check":
            kwargs[spec.name] = bool(val)
        else:
            kwargs[spec.name] = val
    for spec in specs[len(param_vals or []) :]:
        kwargs[spec.name] = spec.default
    return kwargs


# ── Controller builder ────────────────────────────────────────────────────────


def _build_ctrl(
    session_id: str, steps: list, dark: bool, sites_override=None
) -> CorrectionController | None:
    # Always fetch the FULL session cache for strat EDI objects so Stratagem
    # methods work even when sites_override is a line-filtered subset.
    all_sites = cache_get(session_id)
    sites = sites_override if sites_override is not None else all_sites
    if sites is None:
        return None
    ctrl = CorrectionController()
    ctrl.dark = dark
    # Populate _strat_edi_objects from the full (unfiltered) session data so
    # Stratagem corrections (static shift, freq filter, …) never raise
    # "No EDI directory loaded".  The individual step methods deepcopy the list
    # before modifying it, so sharing references here is safe.
    if all_sites is not None:
        raw_edis = list(_iter_edi(all_sites))
        ctrl._strat_edi_objects = raw_edis
        ctrl._strat_raw_edi_objects = raw_edis
    ctrl.set_raw_sites(sites)
    for step in steps or []:
        try:
            ctrl.apply(
                step["fn_name"],
                step.get("kwargs", {}),
                label=step.get("label", step["fn_name"]),
            )
        except Exception:
            pass
    return ctrl


# ── Active-lines filter helper ────────────────────────────────────────────────


def _apply_line_filter(session_id, active_lines_store, store_data):
    """Return (raw_sites, error_sentinel) where sentinel is 'muted' or None."""
    raw = cache_get(session_id)
    if raw is None:
        return None, None
    _als = active_lines_store or {}
    active = _als.get("active", _als.get("all"))
    if active is not None:
        records = (store_data or {}).get("station_records", [])
        filtered = filter_sites_by_lines(raw, records, active)
        if filtered is None:
            return None, "muted"
        raw = filtered
    return raw, None


# ── ρ/φ display-scope filter (line + station selectors) ──────────────────────


def _filter_display_sites(sites, sel_lines, sel_stations, records):
    """Return a Sites filtered to user-selected lines and/or stations for display.

    Both filters are display-only; the correction chain is unchanged.
    If neither selector has a value, returns *sites* unchanged.
    """
    if not sel_lines and not sel_stations:
        return sites
    try:
        target_ids: set | None = None
        if sel_lines and records:
            import pandas as pd

            df = pd.DataFrame(records)
            if "Line" in df.columns and "ID" in df.columns:
                target_ids = set(df[df["Line"].isin(sel_lines)]["ID"].tolist())
        if sel_stations:
            stn_set = set(sel_stations)
            target_ids = (
                (target_ids & stn_set) if target_ids is not None else stn_set
            )

        if (
            not target_ids
        ):  # empty intersection → show nothing, fall back to all
            return sites

        from pycsamt.emtools._core import ensure_sites

        edis = [ed for ed in _iter_edi(sites) if _edi_name(ed) in target_ids]
        if not edis:
            return sites
        return ensure_sites(
            edis, recursive=True, on_dup="replace", strict=False, verbose=0
        )
    except Exception:
        return sites


# ── Pseudosection rendering ───────────────────────────────────────────────────


def _render_ps(
    ctrl: CorrectionController, sites, title: str, dark: bool
) -> str:
    """ρ_a pseudosection as a data-URI PNG."""
    if sites is None:
        return empty_src(dark=dark)
    try:
        fig = plt.figure(figsize=(9, 4))
        ctrl.plot_rho_pseudosection(sites, fig, title=title)
        src = fig_to_src(fig)
        plt.close(fig)
        return src
    except Exception:
        plt.close("all")
        return empty_src(dark=dark)


# ── Overlay helpers ───────────────────────────────────────────────────────────


def _common_freq_grid(all_freqs, n_pts=60):
    """Build a common log-spaced frequency grid spanning the full range."""
    valid = [np.asarray(f).ravel() for f in all_freqs]
    valid = [f for f in valid if f.size > 0]
    if not valid:
        return None
    all_f = np.concatenate(valid)
    all_f = all_f[all_f > 0]
    if all_f.size < 2:
        return None
    return np.geomspace(all_f.min(), all_f.max(), n_pts)


def _interp_to_grid(freq, values, grid):
    """Log-linear interpolation of (freq, values) onto grid; NaN outside range."""
    freq = np.asarray(freq, float)
    values = np.asarray(values, float)
    valid = np.isfinite(freq) & np.isfinite(values) & (freq > 0) & (values > 0)
    if valid.sum() < 2:
        return np.full(len(grid), np.nan)
    log_f = np.log10(freq[valid])
    log_v = np.log10(values[valid])
    order = np.argsort(log_f)
    interped = np.interp(
        np.log10(grid), log_f[order], log_v[order], left=np.nan, right=np.nan
    )
    return 10.0**interped


def _interp_phi_to_grid(freq, values, grid):
    """Linear-frequency interpolation for phase (not log-valued)."""
    freq = np.asarray(freq, float)
    values = np.asarray(values, float)
    valid = np.isfinite(freq) & np.isfinite(values) & (freq > 0)
    if valid.sum() < 2:
        return np.full(len(grid), np.nan)
    log_f = np.log10(freq[valid])
    order = np.argsort(log_f)
    return np.interp(
        np.log10(grid),
        log_f[order],
        values[valid][order],
        left=np.nan,
        right=np.nan,
    )


def _collect_rho_phase_edi(edis, comp):
    """Yield (freq, rho, phase) numpy arrays from a list of EDI objects."""
    for ed in edis:
        z = getattr(ed, "Z", None)
        if z is None:
            continue
        freq = getattr(z, "freq", None)
        if freq is None:
            continue
        freq = np.asarray(freq, float)
        if freq.size == 0:
            continue
        rho = getattr(z, f"res_{comp}", None)
        phi = getattr(z, f"phase_{comp}", None)
        if rho is None:
            continue
        yield (
            freq,
            np.asarray(rho, float),
            (
                np.asarray(phi, float)
                if phi is not None
                else np.full(freq.size, np.nan)
            ),
        )


def _band_stats_rho(triplets, grid):
    """Return (med, p25, p75) for rho interpolated to grid (log space)."""
    mat = np.array([_interp_to_grid(f, r, grid) for f, r, _ in triplets])
    if mat.ndim < 2 or mat.shape[0] == 0:
        nans = np.full(len(grid), np.nan)
        return nans, nans, nans
    return (
        np.nanpercentile(mat, 50, axis=0),
        np.nanpercentile(mat, 25, axis=0),
        np.nanpercentile(mat, 75, axis=0),
    )


def _band_stats_phi(triplets, grid):
    """Return (med, p25, p75) for phase interpolated to grid (linear)."""
    mat = np.array([_interp_phi_to_grid(f, p, grid) for f, _, p in triplets])
    if mat.ndim < 2 or mat.shape[0] == 0:
        nans = np.full(len(grid), np.nan)
        return nans, nans, nans
    return (
        np.nanpercentile(mat, 50, axis=0),
        np.nanpercentile(mat, 25, axis=0),
        np.nanpercentile(mat, 75, axis=0),
    )


def _ovl_ax_style(ax, ax_bg, txt_col, grid_col, log_y=True, log_x=True):
    ax.set_facecolor(ax_bg)
    ax.tick_params(colors=txt_col, labelsize=7)
    for sp in ax.spines.values():
        sp.set_edgecolor(grid_col)
    ax.grid(True, color=grid_col, lw=0.35, ls="--", alpha=0.6)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")


def _render_overlay_median_band(raw_sites, cur_sites, dark, comp="both"):
    """Median ± IQR shaded band — most publishable overlay summary."""
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []
    comps = ["xy", "yx"] if comp == "both" else [comp]

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"

    fig, (ax_rho, ax_phi) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True, facecolor=bg
    )
    fig.subplots_adjust(
        hspace=0.07, left=0.10, right=0.97, top=0.92, bottom=0.10
    )

    any_data = False
    for c in comps:
        col, _ = _COMP_STYLE[c]
        raw_trips = list(_collect_rho_phase_edi(edis_raw, c))
        cur_trips = list(_collect_rho_phase_edi(edis_cur, c))
        if not raw_trips and not cur_trips:
            continue
        all_freqs = [t[0] for t in raw_trips + cur_trips]
        grid = _common_freq_grid(all_freqs)
        if grid is None:
            continue
        T = 1.0 / grid
        lbl = c.upper()
        # Raw — gray shaded + dashed median
        if raw_trips:
            med, p25, p75 = _band_stats_rho(raw_trips, grid)
            ok = np.isfinite(med) & (med > 0)
            if ok.any():
                ax_rho.fill_between(
                    T,
                    p25,
                    p75,
                    alpha=0.15,
                    color="gray",
                    where=np.isfinite(p25) & np.isfinite(p75),
                )
                ax_rho.plot(
                    T[ok],
                    med[ok],
                    ls="--",
                    color="gray",
                    lw=1.5,
                    label=f"Raw {lbl}",
                )
            med_p, p25_p, p75_p = _band_stats_phi(raw_trips, grid)
            okp = np.isfinite(med_p)
            if okp.any():
                ax_phi.fill_between(
                    T,
                    p25_p,
                    p75_p,
                    alpha=0.15,
                    color="gray",
                    where=np.isfinite(p25_p) & np.isfinite(p75_p),
                )
                ax_phi.plot(T[okp], med_p[okp], ls="--", color="gray", lw=1.5)
        # Corrected — component color shaded + solid median
        if cur_trips:
            med, p25, p75 = _band_stats_rho(cur_trips, grid)
            ok = np.isfinite(med) & (med > 0)
            if ok.any():
                ax_rho.fill_between(
                    T,
                    p25,
                    p75,
                    alpha=0.18,
                    color=col,
                    where=np.isfinite(p25) & np.isfinite(p75),
                )
                ax_rho.plot(
                    T[ok],
                    med[ok],
                    ls="-",
                    color=col,
                    lw=2.0,
                    label=f"Corrected {lbl}",
                )
                any_data = True
            med_p, p25_p, p75_p = _band_stats_phi(cur_trips, grid)
            okp = np.isfinite(med_p)
            if okp.any():
                ax_phi.fill_between(
                    T,
                    p25_p,
                    p75_p,
                    alpha=0.18,
                    color=col,
                    where=np.isfinite(p25_p) & np.isfinite(p75_p),
                )
                ax_phi.plot(T[okp], med_p[okp], ls="-", color=col, lw=2.0)

    if not any_data and not edis_raw:
        plt.close(fig)
        return empty_src(dark=dark)

    _ovl_ax_style(ax_rho, ax_bg, txt_col, grid_col, log_y=True)
    _ovl_ax_style(ax_phi, ax_bg, txt_col, grid_col, log_y=False)
    ax_rho.set_ylabel("ρ_a (Ω·m)", color=txt_col, fontsize=9)
    ax_phi.set_ylabel("Phase (°)", color=txt_col, fontsize=9)
    ax_phi.set_xlabel("Period (s)", color=txt_col, fontsize=9)
    ax_rho.set_title(
        f"Station ensemble — Median ± IQR  ({len(edis_raw)} stations)",
        color=txt_col,
        fontsize=10,
        pad=6,
    )
    ax_rho.legend(
        fontsize=8,
        facecolor=ax_bg,
        edgecolor=grid_col,
        labelcolor=txt_col,
        framealpha=0.85,
    )
    fig.patch.set_facecolor(bg)
    plt.setp(ax_rho.get_xticklabels(), visible=False)
    src = fig_to_src(fig)
    plt.close(fig)
    return src


def _render_overlay_spaghetti(
    raw_sites, cur_sites, dark, comp="both", n_max=40
):
    """Individual station curves (thin, semi-transparent) + thick median line."""
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []
    n_show = min(len(edis_raw), int(n_max or 40))
    comps = ["xy", "yx"] if comp == "both" else [comp]

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"

    fig, (ax_rho, ax_phi) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True, facecolor=bg
    )
    fig.subplots_adjust(
        hspace=0.07, left=0.10, right=0.97, top=0.92, bottom=0.10
    )

    for c in comps:
        col, _ = _COMP_STYLE[c]
        raw_trips = list(_collect_rho_phase_edi(edis_raw[:n_show], c))
        cur_trips = list(_collect_rho_phase_edi(edis_cur[:n_show], c))
        all_freqs = [t[0] for t in raw_trips + cur_trips]
        grid = _common_freq_grid(all_freqs)
        if grid is None:
            continue
        T = 1.0 / grid

        # Individual raw curves — gray, very thin, low alpha
        for f, r, _ in raw_trips:
            rr = _interp_to_grid(f, r, grid)
            ok = np.isfinite(rr) & (rr > 0)
            if ok.any():
                ax_rho.plot(T[ok], rr[ok], color="gray", lw=0.6, alpha=0.25)
        # Individual corrected curves — component color, very thin, low alpha
        for f, r, p in cur_trips:
            rr = _interp_to_grid(f, r, grid)
            ok = np.isfinite(rr) & (rr > 0)
            if ok.any():
                ax_rho.plot(T[ok], rr[ok], color=col, lw=0.6, alpha=0.30)
            pp = _interp_phi_to_grid(f, p, grid)
            okp = np.isfinite(pp)
            if okp.any():
                ax_phi.plot(T[okp], pp[okp], color=col, lw=0.6, alpha=0.30)
        # Raw phase — gray, thin
        for f, _, p in raw_trips:
            pp = _interp_phi_to_grid(f, p, grid)
            okp = np.isfinite(pp)
            if okp.any():
                ax_phi.plot(T[okp], pp[okp], color="gray", lw=0.6, alpha=0.25)

        # Thick median lines
        lbl = c.upper()
        if raw_trips:
            med, _, _ = _band_stats_rho(raw_trips, grid)
            ok = np.isfinite(med) & (med > 0)
            if ok.any():
                ax_rho.plot(
                    T[ok],
                    med[ok],
                    ls="--",
                    color="gray",
                    lw=2.2,
                    zorder=5,
                    label=f"Median raw {lbl}",
                )
            med_p, _, _ = _band_stats_phi(raw_trips, grid)
            okp = np.isfinite(med_p)
            if okp.any():
                ax_phi.plot(
                    T[okp],
                    med_p[okp],
                    ls="--",
                    color="gray",
                    lw=2.2,
                    zorder=5,
                )
        if cur_trips:
            med, _, _ = _band_stats_rho(cur_trips, grid)
            ok = np.isfinite(med) & (med > 0)
            if ok.any():
                ax_rho.plot(
                    T[ok],
                    med[ok],
                    ls="-",
                    color=col,
                    lw=2.5,
                    zorder=6,
                    label=f"Median corr {lbl}",
                )
            med_p, _, _ = _band_stats_phi(cur_trips, grid)
            okp = np.isfinite(med_p)
            if okp.any():
                ax_phi.plot(
                    T[okp], med_p[okp], ls="-", color=col, lw=2.5, zorder=6
                )

    _ovl_ax_style(ax_rho, ax_bg, txt_col, grid_col, log_y=True)
    _ovl_ax_style(ax_phi, ax_bg, txt_col, grid_col, log_y=False)
    ax_rho.set_ylabel("ρ_a (Ω·m)", color=txt_col, fontsize=9)
    ax_phi.set_ylabel("Phase (°)", color=txt_col, fontsize=9)
    ax_phi.set_xlabel("Period (s)", color=txt_col, fontsize=9)
    ax_rho.set_title(
        f"{n_show} stations — thin=individual, thick=median",
        color=txt_col,
        fontsize=10,
        pad=6,
    )
    ax_rho.legend(
        fontsize=8,
        facecolor=ax_bg,
        edgecolor=grid_col,
        labelcolor=txt_col,
        framealpha=0.85,
    )
    plt.setp(ax_rho.get_xticklabels(), visible=False)
    fig.patch.set_facecolor(bg)
    src = fig_to_src(fig)
    plt.close(fig)
    return src


def _render_overlay_per_line(
    raw_sites, cur_sites, dark, comp="both", records=None
):
    """One column per survey line, 2 rows (ρ_a, φ) — median before/after per line."""
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []
    records = records or []
    comps = ["xy", "yx"] if comp == "both" else [comp]

    # Build station-name → line mapping
    stn_to_line: dict[str, str] = {
        r["ID"]: r.get("Line", "L?") for r in records if "ID" in r
    }
    # Group EDI indices by line (fallback: all in one group)
    line_groups: dict[str, list[int]] = {}
    for i, ed in enumerate(edis_raw):
        stn = (
            getattr(ed, "station", None)
            or getattr(ed, "dataid", None)
            or f"S{i + 1}"
        )
        line = stn_to_line.get(stn, stn_to_line.get(stn.strip(), "All"))
        line_groups.setdefault(line, []).append(i)

    lines = sorted(line_groups.keys())
    n_lines = min(len(lines), 5)  # cap at 5 columns
    if n_lines == 0:
        return empty_src(dark=dark)

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"

    fig_w = max(5.0, 3.8 * n_lines)
    fig, axes = plt.subplots(
        2, n_lines, figsize=(fig_w, 6), sharex=False, facecolor=bg
    )
    fig.subplots_adjust(
        hspace=0.08, wspace=0.38, left=0.09, right=0.97, top=0.92, bottom=0.10
    )
    if n_lines == 1:
        axes = np.array(axes).reshape(2, 1)

    for col_i, line_name in enumerate(lines[:n_lines]):
        idxs = line_groups[line_name]
        ax_rho = axes[0, col_i]
        ax_phi = axes[1, col_i]
        _ovl_ax_style(ax_rho, ax_bg, txt_col, grid_col, log_y=True)
        _ovl_ax_style(ax_phi, ax_bg, txt_col, grid_col, log_y=False)
        n_stn = len(idxs)
        ax_rho.set_title(
            f"{line_name} ({n_stn} stn)", color=txt_col, fontsize=8.5, pad=3
        )
        if col_i == 0:
            ax_rho.set_ylabel("ρ_a (Ω·m)", color=txt_col, fontsize=8)
            ax_phi.set_ylabel("Phase (°)", color=txt_col, fontsize=8)
        ax_phi.set_xlabel("Period (s)", color=txt_col, fontsize=7.5)
        plt.setp(ax_rho.get_xticklabels(), visible=False)

        for c in comps:
            col_c, _ = _COMP_STYLE[c]
            raw_sub = [edis_raw[i] for i in idxs if i < len(edis_raw)]
            cur_sub = [edis_cur[i] for i in idxs if i < len(edis_cur)]
            raw_trips = list(_collect_rho_phase_edi(raw_sub, c))
            cur_trips = list(_collect_rho_phase_edi(cur_sub, c))
            all_freqs = [t[0] for t in raw_trips + cur_trips]
            grid = _common_freq_grid(all_freqs, n_pts=50)
            if grid is None:
                continue
            T = 1.0 / grid
            lbl = c.upper()
            if raw_trips:
                med, p25, p75 = _band_stats_rho(raw_trips, grid)
                ok = np.isfinite(med) & (med > 0)
                if ok.any():
                    ax_rho.fill_between(
                        T,
                        p25,
                        p75,
                        alpha=0.12,
                        color="gray",
                        where=np.isfinite(p25) & np.isfinite(p75),
                    )
                    ax_rho.plot(
                        T[ok],
                        med[ok],
                        ls="--",
                        color="gray",
                        lw=1.5,
                        label=f"Raw {lbl}" if col_i == 0 else None,
                    )
                med_p, p25_p, p75_p = _band_stats_phi(raw_trips, grid)
                okp = np.isfinite(med_p)
                if okp.any():
                    ax_phi.fill_between(
                        T,
                        p25_p,
                        p75_p,
                        alpha=0.12,
                        color="gray",
                        where=np.isfinite(p25_p) & np.isfinite(p75_p),
                    )
                    ax_phi.plot(
                        T[okp], med_p[okp], ls="--", color="gray", lw=1.5
                    )
            if cur_trips:
                med, p25, p75 = _band_stats_rho(cur_trips, grid)
                ok = np.isfinite(med) & (med > 0)
                if ok.any():
                    ax_rho.fill_between(
                        T,
                        p25,
                        p75,
                        alpha=0.15,
                        color=col_c,
                        where=np.isfinite(p25) & np.isfinite(p75),
                    )
                    ax_rho.plot(
                        T[ok],
                        med[ok],
                        ls="-",
                        color=col_c,
                        lw=1.8,
                        label=f"Corr {lbl}" if col_i == 0 else None,
                    )
                med_p, p25_p, p75_p = _band_stats_phi(cur_trips, grid)
                okp = np.isfinite(med_p)
                if okp.any():
                    ax_phi.fill_between(
                        T,
                        p25_p,
                        p75_p,
                        alpha=0.15,
                        color=col_c,
                        where=np.isfinite(p25_p) & np.isfinite(p75_p),
                    )
                    ax_phi.plot(
                        T[okp], med_p[okp], ls="-", color=col_c, lw=1.8
                    )

        if col_i == 0:
            ax_rho.legend(
                fontsize=7,
                facecolor=ax_bg,
                edgecolor=grid_col,
                labelcolor=txt_col,
                framealpha=0.85,
            )

    fig.patch.set_facecolor(bg)
    src = fig_to_src(fig)
    plt.close(fig)
    return src


def _render_overlay(
    ctrl: CorrectionController,
    before,
    after,
    dark: bool,
    ovl_style: str = "median-band",
    comp: str = "both",
    records: list | None = None,
) -> str:
    """Dispatch to one of three overlay renderers based on ovl_style."""
    try:
        if ovl_style == "spaghetti":
            return _render_overlay_spaghetti(before, after, dark, comp=comp)
        elif ovl_style == "per-line":
            return _render_overlay_per_line(
                before, after, dark, comp=comp, records=records
            )
        else:  # "median-band" (default)
            return _render_overlay_median_band(before, after, dark, comp=comp)
    except Exception:
        plt.close("all")
        return empty_src(dark=dark)


def _render_diff(ctrl: CorrectionController, before, after, dark: bool) -> str:
    try:
        fig, ax = plt.subplots(figsize=(9, 4))
        ctrl.plot_diff(before, after, ax)
        src = fig_to_src(fig)
        plt.close(fig)
        return src
    except Exception:
        plt.close("all")
        return empty_src(dark=dark)


# ── Per-station ρ/φ grid ──────────────────────────────────────────────────────


def _safe_attr(z_obj, *attrs, freq_arr=None):
    """Try multiple attribute names and return a 1-D float array or None."""
    for a in attrs:
        v = getattr(z_obj, a, None)
        if v is not None:
            try:
                arr = np.asarray(v, float).ravel()
                if freq_arr is not None:
                    if arr.size == freq_arr.size:
                        return arr
                else:
                    return arr
            except Exception:
                pass
    return None


# API-defined component colors and markers (from PYCSAMT_STYLE.mt)
_XY_COL = "#1f77b4"  # mt.xy.color
_YX_COL = "#d62728"  # mt.yx.color
_XY_MRK = "o"  # mt.xy.marker
_YX_MRK = "s"  # mt.yx.marker
_COMP_STYLE = {
    "xy": (_XY_COL, _XY_MRK),
    "yx": (_YX_COL, _YX_MRK),
}


def _get_comp(z_obj, comp: str):
    """Return (T_periods, rho, phi) for component 'xy' or 'yx'.

    Uses the correct Z attribute names: res_xy/res_yx and phase_xy/phase_yx.
    Returns (empty_array, None, None) when data is unavailable.
    """
    if z_obj is None:
        return np.array([]), None, None
    freq = getattr(z_obj, "freq", None)
    if freq is None:
        return np.array([]), None, None
    freq = np.asarray(freq, float)
    if freq.size == 0:
        return np.array([]), None, None
    T = 1.0 / np.where(freq > 0, freq, np.nan)
    rho = getattr(z_obj, f"res_{comp}", None)
    phi = getattr(z_obj, f"phase_{comp}", None)
    return (
        T,
        (np.asarray(rho, float) if rho is not None else None),
        (np.asarray(phi, float) if phi is not None else None),
    )


def _plot_comp(
    ax_rho,
    ax_phi,
    T,
    rho,
    phi,
    color,
    marker,
    ls,
    mfc,
    lw,
    ms,
    label=None,
    mew=1.0,
    zorder=2,
    alpha=1.0,
):
    """Plot one (rho, phi) trace on the given axes pair."""
    kw_rho = dict(
        color=color,
        lw=lw,
        ms=ms,
        marker=marker,
        mfc=mfc,
        ls=ls,
        mew=mew,
        zorder=zorder,
        alpha=alpha,
    )
    kw_phi = dict(**kw_rho)
    if rho is not None and T.size > 0:
        ok = np.isfinite(T) & np.isfinite(rho) & (rho > 0)
        if ok.any():
            if label:
                kw_rho["label"] = label
            ax_rho.loglog(T[ok], rho[ok], **kw_rho)
    if phi is not None and T.size > 0:
        ok = np.isfinite(T) & np.isfinite(phi)
        if ok.any():
            ax_phi.semilogx(T[ok], phi[ok], **kw_phi)


def _ax_style(ax, ax_bg, txt_col, grid_col):
    ax.set_facecolor(ax_bg)
    ax.tick_params(colors=txt_col, labelsize=6.0)
    for sp in ax.spines.values():
        sp.set_edgecolor(grid_col)
    ax.grid(True, color=grid_col, lw=0.35, ls="--", alpha=0.7)


def _legend_handles(comps, mode, ax_bg, grid_col, txt_col):
    """Build legend proxy artists for the active comp+mode combination."""
    import matplotlib.lines as mlines

    handles = []
    for c in comps:
        col, mrk = _COMP_STYLE[c]
        if mode in ("raw", "both"):
            handles.append(
                mlines.Line2D(
                    [],
                    [],
                    color=col,
                    marker=mrk,
                    ls="-",
                    mfc=col,
                    ms=4,
                    lw=1.2,
                    label=f"raw {c.upper()}",
                )
            )
        if mode in ("corrected", "both"):
            handles.append(
                mlines.Line2D(
                    [],
                    [],
                    color=col,
                    marker=mrk,
                    ls=":",
                    mfc=ax_bg,
                    mew=1.1,
                    ms=4,
                    lw=1.5,
                    label=f"corr {c.upper()}",
                )
            )
    return handles


def _render_rho_phi_grid(
    ctrl,
    raw_sites,
    cur_sites,
    mode: str,
    n_max: int,
    dark: bool,
    comp: str = "both",
) -> str:
    """Per-station ρ_a/φ grid — XY and YX overlaid, raw vs corrected overlaid.

    Encoding
    --------
    Color  → component: XY = blue (#1f77b4), YX = red (#d62728)
    Line   → source:    raw = solid,  corrected = dotted
    Marker → source:    raw = filled, corrected = open (bg-fill)
    """
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []

    n = min(len(edis_raw), max(1, int(n_max or 6)))
    if n == 0:
        return empty_src(dark=dark)

    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig_w = max(5.0, 3.8 * ncols)
    fig_h = max(4.5, 3.2 * nrows)

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"

    comps = ["xy", "yx"] if comp == "both" else [comp]

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=bg)
    outer = gridspec.GridSpec(
        nrows, ncols, figure=fig, hspace=0.60, wspace=0.42
    )

    for idx in range(n):
        row_g = idx // ncols
        col_g = idx % ncols
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row_g, col_g], hspace=0.08
        )
        ax_rho = fig.add_subplot(inner[0])
        ax_phi = fig.add_subplot(inner[1], sharex=ax_rho)
        for ax in (ax_rho, ax_phi):
            _ax_style(ax, ax_bg, txt_col, grid_col)

        ed_raw = edis_raw[idx]
        station = getattr(ed_raw, "station", None) or getattr(
            ed_raw, "dataid", f"S{idx + 1}"
        )
        z_raw = getattr(ed_raw, "Z", None)

        if z_raw is None:
            ax_rho.text(
                0.5,
                0.5,
                "no Z",
                ha="center",
                va="center",
                color=txt_col,
                fontsize=7,
                transform=ax_rho.transAxes,
            )
            ax_rho.set_title(station, fontsize=7.5, color=txt_col, pad=2)
            continue

        # Corrected EDI (matched by position)
        z_cur = None
        if idx < len(edis_cur):
            z_cur = getattr(edis_cur[idx], "Z", None)

        for c in comps:
            col, mrk = _COMP_STYLE[c]
            # Raw traces — solid, filled marker
            if mode in ("raw", "both"):
                T_r, rho_r, phi_r = _get_comp(z_raw, c)
                _plot_comp(
                    ax_rho,
                    ax_phi,
                    T_r,
                    rho_r,
                    phi_r,
                    color=col,
                    marker=mrk,
                    ls="-",
                    mfc=col,
                    lw=1.2,
                    ms=3.5,
                    zorder=2,
                )
            # Corrected traces — dotted, open marker
            if mode in ("corrected", "both") and z_cur is not None:
                T_c, rho_c, phi_c = _get_comp(z_cur, c)
                _plot_comp(
                    ax_rho,
                    ax_phi,
                    T_c,
                    rho_c,
                    phi_c,
                    color=col,
                    marker=mrk,
                    ls=":",
                    mfc=ax_bg,
                    lw=1.5,
                    ms=4.0,
                    mew=1.1,
                    zorder=3,
                )

        ax_rho.set_title(station, fontsize=7.5, color=txt_col, pad=2)
        ax_rho.set_ylabel("ρ_a\n(Ω·m)", fontsize=6.0, color=txt_col)
        ax_phi.set_ylabel("φ (°)", fontsize=6.0, color=txt_col)
        ax_phi.set_xlabel("T (s)", fontsize=6.0, color=txt_col)
        plt.setp(ax_rho.get_xticklabels(), visible=False)

        if idx == 0:
            handles = _legend_handles(comps, mode, ax_bg, grid_col, txt_col)
            if handles:
                ax_rho.legend(
                    handles=handles,
                    fontsize=5.5,
                    loc="upper right",
                    facecolor=ax_bg,
                    edgecolor=grid_col,
                    labelcolor=txt_col,
                    framealpha=0.85,
                )

    src = fig_to_src(fig)
    plt.close(fig)
    return src


# ── Station-Pairs renderer (static shift view) ───────────────────────────────


def _render_paired_stations(
    ctrl,
    raw_sites,
    cur_sites,
    mode: str,
    n_max: int,
    dark: bool,
    comp: str = "both",
) -> str:
    """Side-by-side raw vs corrected per-station ρ/φ (static-shift view).

    Each station block = two panels [raw | corrected].
    Each panel shows XY and YX components on shared rho/phase axes.
    Encoding: XY = blue (#1f77b4, ○), YX = red (#d62728, s)
    Useful for comparing the vertical shift in log(ρ_a) from static shift.
    """
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []

    n = min(len(edis_raw), max(1, int(n_max or 4)))
    if n == 0:
        return empty_src(dark=dark)

    comps = ["xy", "yx"] if comp == "both" else [comp]
    fig_w = max(6.0, 3.6 * n * 2)
    fig_h = 5.5

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"
    sep_col = "#45475a" if dark else "#c1c2d0"
    hdr_raw = "#6c7086" if dark else "#9ca0b0"
    hdr_cor = "#89b4fa" if dark else "#1e66f5"

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=bg)
    outer = gridspec.GridSpec(1, n, figure=fig, wspace=0.06)

    for s_idx in range(n):
        inner = gridspec.GridSpecFromSubplotSpec(
            2,
            2,
            subplot_spec=outer[0, s_idx],
            hspace=0.06,
            wspace=0.30,
        )
        ax_rr = fig.add_subplot(inner[0, 0])
        ax_cr = fig.add_subplot(inner[0, 1])
        ax_rp = fig.add_subplot(inner[1, 0], sharex=ax_rr)
        ax_cp = fig.add_subplot(inner[1, 1], sharex=ax_cr)

        for ax in (ax_rr, ax_cr, ax_rp, ax_cp):
            _ax_style(ax, ax_bg, txt_col, grid_col)

        ed_raw = edis_raw[s_idx]
        station = getattr(ed_raw, "station", None) or getattr(
            ed_raw, "dataid", f"S{s_idx + 1}"
        )
        z_raw = getattr(ed_raw, "Z", None)

        # Station title spanning both sub-cols
        mid_x = (
            inner[0, 0].get_position(fig).x0 + inner[0, 1].get_position(fig).x1
        ) / 2
        top_y = inner[0, 0].get_position(fig).y1 + 0.026
        fig.text(
            mid_x,
            top_y,
            station,
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=txt_col,
            fontweight="bold",
        )

        ax_rr.set_title("/raw", fontsize=7, color=hdr_raw, pad=2)
        ax_cr.set_title("/corrected", fontsize=7, color=hdr_cor, pad=2)

        if s_idx > 0:
            lx = inner[0, 0].get_position(fig).x0 - 0.01
            fig.add_artist(
                plt.Line2D(
                    [lx, lx],
                    [0.06, 0.95],
                    transform=fig.transFigure,
                    color=sep_col,
                    lw=0.8,
                    ls="--",
                )
            )

        if z_raw is None:
            for ax in (ax_rr, ax_rp):
                ax.text(
                    0.5,
                    0.5,
                    "no Z",
                    ha="center",
                    va="center",
                    color=txt_col,
                    fontsize=7,
                    transform=ax.transAxes,
                )
            for ax in (ax_cr, ax_cp):
                ax.text(
                    0.5,
                    0.5,
                    "—",
                    ha="center",
                    va="center",
                    color=txt_col,
                    fontsize=7,
                    transform=ax.transAxes,
                )
            plt.setp([ax_rr, ax_cr], xticklabels=[])
            continue

        z_cur = (
            getattr(edis_cur[s_idx], "Z", None)
            if s_idx < len(edis_cur)
            else None
        )

        # Draw raw panel
        for c in comps:
            col, mrk = _COMP_STYLE[c]
            T_r, rho_r, phi_r = _get_comp(z_raw, c)
            lbl = c.upper() if s_idx == 0 else None
            _plot_comp(
                ax_rr,
                ax_rp,
                T_r,
                rho_r,
                phi_r,
                color=col,
                marker=mrk,
                ls="-",
                mfc=col,
                lw=1.2,
                ms=3.5,
                label=lbl,
            )

        # Draw corrected panel
        if z_cur is not None:
            for c in comps:
                col, mrk = _COMP_STYLE[c]
                T_c, rho_c, phi_c = _get_comp(z_cur, c)
                _plot_comp(
                    ax_cr,
                    ax_cp,
                    T_c,
                    rho_c,
                    phi_c,
                    color=col,
                    marker=mrk,
                    ls="-",
                    mfc=col,
                    lw=1.2,
                    ms=3.5,
                )
        else:
            for ax in (ax_cr, ax_cp):
                ax.text(
                    0.5,
                    0.5,
                    "no corrected\ndata",
                    ha="center",
                    va="center",
                    color=txt_col,
                    fontsize=6.5,
                    transform=ax.transAxes,
                )

        if s_idx == 0:
            ax_rr.set_ylabel("ρ_a (Ω·m)", fontsize=6.0, color=txt_col)
            ax_rp.set_ylabel("φ (°)", fontsize=6.0, color=txt_col)
            handles = _legend_handles(comps, "raw", ax_bg, grid_col, txt_col)
            # relabel to just XY / YX
            for h, c in zip(handles, comps):
                h.set_label(c.upper())
            if handles:
                ax_rr.legend(
                    handles=handles,
                    fontsize=5.5,
                    loc="upper right",
                    facecolor=ax_bg,
                    edgecolor=grid_col,
                    labelcolor=txt_col,
                    framealpha=0.85,
                )

        ax_rp.set_xlabel("T (s)", fontsize=6.0, color=txt_col)
        ax_cp.set_xlabel("T (s)", fontsize=6.0, color=txt_col)
        plt.setp([ax_rr, ax_cr], xticklabels=[])

    fig.subplots_adjust(top=0.87, bottom=0.08)
    src = fig_to_src(fig)
    plt.close(fig)
    return src


# ── Freq-Bands renderer (Task 2 — drop-frequency view) ───────────────────────


def _removed_periods(ed_raw, ed_cur, tol_rel: float = 0.02):
    """Return 1-D period array of frequencies present in raw but absent in cur.

    Uses a relative tolerance to handle floating-point rounding when comparing
    frequency arrays from two EDI objects.
    """
    z_raw = getattr(ed_raw, "Z", None)
    z_cur = getattr(ed_cur, "Z", None) if ed_cur is not None else None
    if z_raw is None:
        return np.array([])
    freq_raw = np.asarray(
        _f if (_f := getattr(z_raw, "freq", None)) is not None else [], float
    )
    if freq_raw.size == 0:
        return np.array([])
    if z_cur is None:
        return 1.0 / freq_raw[freq_raw > 0]

    freq_cur = np.asarray(
        _f if (_f := getattr(z_cur, "freq", None)) is not None else [], float
    )
    if freq_cur.size == 0:
        return 1.0 / freq_raw[freq_raw > 0]

    removed = []
    for f in freq_raw:
        if f <= 0:
            continue
        if freq_cur.size == 0 or np.all(
            np.abs(freq_cur - f) / (f + 1e-12) > tol_rel
        ):
            removed.append(1.0 / f)
    return np.asarray(removed)


def _render_freq_bands(
    ctrl,
    raw_sites,
    cur_sites,
    mode: str,
    n_max: int,
    dark: bool,
    comp: str = "both",
    dropped_freqs: list | None = None,
) -> str:
    """Per-station ρ/φ grid with removed frequency bands highlighted.

    Removed bands shown as faint red axvspan, removed data points on the
    raw curve marked with open red circles/squares.
    XY and YX plotted with their API colors.
    """
    edis_raw = list(_iter_edi(raw_sites))
    edis_cur = list(_iter_edi(cur_sites)) if cur_sites is not None else []

    n = min(len(edis_raw), max(1, int(n_max or 6)))
    if n == 0:
        return empty_src(dark=dark)

    comps = ["xy", "yx"] if comp == "both" else [comp]
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig_w = max(5.0, 3.8 * ncols)
    fig_h = max(4.5, 3.2 * nrows)

    bg = "#1e1e2e" if dark else "#eff1f5"
    ax_bg = "#181825" if dark else "#ffffff"
    txt_col = "#cdd6f4" if dark else "#4c4f69"
    grid_col = "#313244" if dark else "#dce0e8"
    drop_col = "#f38ba8" if dark else "#d20f39"

    fig = plt.figure(figsize=(fig_w, fig_h), facecolor=bg)
    outer = gridspec.GridSpec(
        nrows, ncols, figure=fig, hspace=0.60, wspace=0.42
    )

    for idx in range(n):
        row_g = idx // ncols
        col_g = idx % ncols
        inner = gridspec.GridSpecFromSubplotSpec(
            2, 1, subplot_spec=outer[row_g, col_g], hspace=0.08
        )
        ax_rho = fig.add_subplot(inner[0])
        ax_phi = fig.add_subplot(inner[1], sharex=ax_rho)
        for ax in (ax_rho, ax_phi):
            _ax_style(ax, ax_bg, txt_col, grid_col)

        ed_raw = edis_raw[idx]
        ed_cur = edis_cur[idx] if idx < len(edis_cur) else None
        station = getattr(ed_raw, "station", None) or getattr(
            ed_raw, "dataid", f"S{idx + 1}"
        )
        z_raw = getattr(ed_raw, "Z", None)
        z_cur = getattr(ed_cur, "Z", None) if ed_cur is not None else None

        if z_raw is None:
            ax_rho.set_title(station, fontsize=7.5, color=txt_col, pad=2)
            continue

        # Removed periods: use explicit list when provided (manual drop),
        # otherwise auto-detect from raw vs cur freq array comparison.
        if dropped_freqs is not None:
            drop_arr = np.asarray(
                [float(f) for f in dropped_freqs if f > 0], float
            )
            T_removed = (1.0 / drop_arr) if drop_arr.size > 0 else np.array([])
        else:
            T_removed = _removed_periods(ed_raw, ed_cur)

        # Shade removed bands
        for T_r in np.sort(T_removed):
            ax_rho.axvspan(
                T_r * 0.70, T_r * 1.43, alpha=0.18, color=drop_col, zorder=1
            )
            ax_phi.axvspan(
                T_r * 0.70, T_r * 1.43, alpha=0.18, color=drop_col, zorder=1
            )

        for c in comps:
            col, mrk = _COMP_STYLE[c]
            T_r, rho_r, phi_r = _get_comp(z_raw, c)

            # Raw traces
            if mode in ("raw", "both") and T_r.size > 0:
                _plot_comp(
                    ax_rho,
                    ax_phi,
                    T_r,
                    rho_r,
                    phi_r,
                    color=col,
                    marker=mrk,
                    ls="-",
                    mfc=col,
                    lw=0.9,
                    ms=3.0,
                    alpha=0.6,
                    zorder=2,
                )
                # Mark removed data points on raw curve
                if T_removed.size > 0 and rho_r is not None:
                    for T_drop in T_removed:
                        diffs = np.abs(T_r - T_drop)
                        mi = np.argmin(diffs)
                        if diffs[mi] / (T_drop + 1e-12) < 0.1:
                            if (
                                rho_r is not None
                                and rho_r[mi] > 0
                                and np.isfinite(rho_r[mi])
                            ):
                                ax_rho.loglog(
                                    [T_r[mi]],
                                    [rho_r[mi]],
                                    mrk,
                                    color=drop_col,
                                    ms=5.5,
                                    mfc="none",
                                    mew=1.3,
                                    zorder=5,
                                )
                            if phi_r is not None and np.isfinite(phi_r[mi]):
                                ax_phi.semilogx(
                                    [T_r[mi]],
                                    [phi_r[mi]],
                                    mrk,
                                    color=drop_col,
                                    ms=5.5,
                                    mfc="none",
                                    mew=1.3,
                                    zorder=5,
                                )

            # Corrected traces (retained data)
            if mode in ("corrected", "both") and z_cur is not None:
                T_c, rho_c, phi_c = _get_comp(z_cur, c)
                _plot_comp(
                    ax_rho,
                    ax_phi,
                    T_c,
                    rho_c,
                    phi_c,
                    color=col,
                    marker=mrk,
                    ls=":",
                    mfc=ax_bg,
                    lw=1.4,
                    ms=4.0,
                    mew=1.1,
                    zorder=3,
                )

        n_drop = T_removed.size
        if n_drop > 0:
            ax_rho.text(
                0.98,
                0.97,
                f"−{n_drop} freq",
                ha="right",
                va="top",
                transform=ax_rho.transAxes,
                fontsize=5.5,
                color=drop_col,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    fc=ax_bg,
                    ec=drop_col,
                    lw=0.6,
                    alpha=0.85,
                ),
            )

        ax_rho.set_title(station, fontsize=7.5, color=txt_col, pad=2)
        ax_rho.set_ylabel("ρ_a\n(Ω·m)", fontsize=6.0, color=txt_col)
        ax_phi.set_ylabel("φ (°)", fontsize=6.0, color=txt_col)
        ax_phi.set_xlabel("T (s)", fontsize=6.0, color=txt_col)
        plt.setp(ax_rho.get_xticklabels(), visible=False)

        if idx == 0:
            handles = _legend_handles(comps, mode, ax_bg, grid_col, txt_col)
            if handles:
                ax_rho.legend(
                    handles=handles,
                    fontsize=5.5,
                    loc="upper right",
                    facecolor=ax_bg,
                    edgecolor=grid_col,
                    labelcolor=txt_col,
                    framealpha=0.85,
                )

    src = fig_to_src(fig)
    plt.close(fig)
    return src


# ── Data extraction helpers ───────────────────────────────────────────────────


def _collect_rho(sites):
    try:
        arr = []
        for ed in _iter_edi(sites):
            z = getattr(ed, "Z", None)
            if z is None:
                continue
            for attr in ("resistivity", "rho_xy", "resistivity_xy", "rho"):
                v = getattr(z, attr, None)
                if v is not None:
                    a = np.asarray(v, float).ravel()
                    arr.append(a[np.isfinite(a)])
                    break
        return np.concatenate(arr) if arr else None
    except Exception:
        return None


def _iter_edi(sites):
    """Yield EDIFile objects from a Sites, list-of-EDI, or similar container.

    Sites.__iter__ yields Site wrappers, not EDIFile objects.  Sites.as_list()
    returns the underlying EDIFile list; we use that when available.
    Site objects carry a `.edi` attribute with the EDIFile — we unwrap it.
    """
    if sites is None:
        return
    # Sites objects have as_list() → [EDIFile, ...]
    if hasattr(sites, "as_list"):
        try:
            yield from sites.as_list()
            return
        except Exception:
            pass
    # Explicit EDI-list attributes
    for attr in ("edifiles", "edi_list", "_edi_list"):
        v = getattr(sites, attr, None)
        if v is not None:
            try:
                for item in v:
                    yield getattr(item, "edi", item)
            except Exception:
                pass
            return
    # Generic iteration: Site wrappers → unwrap to .edi; plain EDIFiles pass through
    try:
        for item in sites:
            yield getattr(item, "edi", item)
    except Exception:
        pass


def _edi_name(ed) -> str:
    """Return the station name from an EDIFile (or Site) object."""
    if hasattr(ed, "edi"):  # Site wrapper
        ed = ed.edi
    # EDIFile.station reads dataid from the HEAD section
    return (
        getattr(ed, "station", None)
        or getattr(ed, "dataid", None)
        or getattr(ed, "name", "")
        or ""
    )


def _edi_to_rows(ed) -> list:
    rows = []
    try:
        station = getattr(ed, "station", None) or getattr(ed, "dataid", "")
        z = getattr(ed, "Z", None)
        if z is None:
            return rows
        freq = np.asarray(
            _f if (_f := getattr(z, "freq", None)) is not None else [], float
        )
        for attr, col in [
            ("resistivity_xy", "rho_xy"),
            ("resistivity_yx", "rho_yx"),
            ("phase_xy", "phase_xy"),
            ("phase_yx", "phase_yx"),
        ]:
            v = getattr(z, attr, None)
            if v is not None and len(np.asarray(v).ravel()) == len(freq):
                for i, f in enumerate(freq):
                    d = {"station": station, "freq_hz": float(f)}
                    d[col] = float(np.asarray(v).ravel()[i])
                    rows.append(d)
    except Exception:
        pass
    return rows


# ── UI builder helpers ────────────────────────────────────────────────────────


def _diff_stats_strip(mn, mx, mu, sd, n_neg, n_pos):
    def _kv(label, value, color="#cdd6f4"):
        return html.Div(
            [
                html.Span(
                    label,
                    style={
                        "fontSize": "9px",
                        "color": "#6c7086",
                        "display": "block",
                    },
                ),
                html.Span(
                    value,
                    style={
                        "fontSize": "12px",
                        "color": color,
                        "fontWeight": "600",
                    },
                ),
            ],
            className="corr-diff-kv",
        )

    return html.Div(
        [
            _kv("Δlog₁₀(ρₐ) mean", f"{mu:+.3f}", "#89b4fa"),
            _kv("Δlog₁₀(ρₐ) σ", f"{sd:.3f}", "#6c7086"),
            _kv("min Δ", f"{mn:+.3f}", "#a6e3a1" if mn >= 0 else "#f38ba8"),
            _kv("max Δ", f"{mx:+.3f}", "#a6e3a1" if mx >= 0 else "#f38ba8"),
            _kv("cells increased", str(n_pos), "#a6e3a1"),
            _kv("cells decreased", str(n_neg), "#f38ba8"),
        ],
        className="corr-diff-stats-row",
    )


def _stack_el(steps: list):
    if not steps:
        return html.Span(
            "No corrections applied yet.", className="text-muted small"
        )
    items = []
    for i, s in enumerate(steps):
        label = s.get("label", s["fn_name"])
        kw = s.get("kwargs", {})
        kw_str = ", ".join(f"{k}={v}" for k, v in kw.items()) if kw else ""
        items.append(
            html.Li(
                [
                    html.Span(f"{i + 1}", className="corr-stack-idx"),
                    html.Div(
                        [
                            html.Span(label, className="corr-stack-label"),
                            html.Span(kw_str, className="corr-stack-kw")
                            if kw_str
                            else None,
                        ],
                        className="corr-stack-body",
                    ),
                    dbc.Button(
                        html.I(className="bi bi-x"),
                        id={"type": "corr-step-del", "index": i},
                        color="danger",
                        outline=True,
                        size="sm",
                        className="corr-step-del-btn",
                        n_clicks=0,
                    ),
                ],
                className="corr-stack-item",
            ),
        )
    return html.Ul(items, className="corr-stack-list")


def _chain_badge(n: int) -> str:
    return f"{n} step{'s' if n != 1 else ''}" if n > 0 else ""


def _ctx_bar_content(steps: list) -> list:
    n = len(steps)
    if n == 0:
        return [
            html.I(className="bi bi-info-circle me-1"),
            "Chain is empty — raw data shown.",
        ]
    last = steps[-1].get("label", steps[-1]["fn_name"])
    return [
        html.I(className="bi bi-layers me-1", style={"color": "var(--green)"}),
        html.Span(
            f"{n} correction{'s' if n != 1 else ''} applied",
            style={"fontWeight": "600", "color": "var(--green)"},
        ),
        html.Span(" · ", style={"color": "var(--sub0)"}),
        html.Span(f"Last: {last}", style={"color": "var(--sub0)"}),
    ]


# ── ρ/φ view hint builder ────────────────────────────────────────────────────

_STATIC_SHIFT_FNS = {
    "correct_ss_ama",
    "_correct_ss_loess",
    "_correct_ss_bilateral",
    "_correct_ss_refmedian",
    "correct_static_shift",
    "_strat_static_shift",
}
_DROP_FREQ_FNS = {"_strat_freq_filter", "drop_freqs_manual"}


def _rho_phi_hint(steps: list, current_style: str) -> str:
    """Return a context-aware hint for the ρ/φ view bar."""
    fn_names = {s.get("fn_name", "") for s in steps}
    has_ss = bool(fn_names & _STATIC_SHIFT_FNS)
    has_df = bool(fn_names & _DROP_FREQ_FNS)

    if current_style == "pairs":
        return (
            "Station Pairs: left = raw station, right = after correction. "
            "Ideal for inspecting vertical log(ρ_a) offsets from static shift."
        )
    if current_style == "freq-bands":
        if has_df:
            return (
                "Freq Bands: red-shaded regions = removed frequencies, "
                "open red circles = removed data points on raw curve, "
                "blue line = retained corrected data."
            )
        return (
            "Freq Bands: highlights frequencies removed by drop-bad-freq. "
            "Apply 'Drop Bad Frequencies' first to see removed bands."
        )
    # Default grid style
    if has_ss:
        return (
            "Tip: Switch to 'Station Pairs' to inspect static shift — "
            "shows raw vs corrected side-by-side per station."
        )
    if has_df:
        return (
            "Tip: Switch to 'Freq Bands' to see which frequencies were removed "
            "and their effect on the ρ_a/φ curves."
        )
    return (
        "ρ/φ Grid: per-station ρ_a/φ vs period. "
        "Use 'Station Pairs' for static-shift view, 'Freq Bands' for drop-freq view."
    )


# ── Params panel builder ──────────────────────────────────────────────────────


def _params_panel(cat: str, method: str):
    specs = _get_specs(cat, method)
    if not specs:
        return html.Div(
            "No parameters for this correction.",
            className="text-muted small mt-1",
        )
    rows = []
    for spec in specs:
        lbl = html.Div(
            spec.label, className="param-label", title=(spec.tip or "")
        )
        if spec.kind in ("spin", "dspin"):
            lo, hi, step_ = spec.opts
            widget = dbc.Input(
                id={"type": "corr-param", "name": spec.name},
                type="number",
                min=lo,
                max=hi,
                step=step_,
                value=spec.default,
                size="sm",
                className="param-input",
                debounce=True,
            )
        elif spec.kind == "combo":
            widget = dbc.Select(
                id={"type": "corr-param", "name": spec.name},
                options=[{"label": o, "value": o} for o in spec.opts],
                value=spec.default,
                size="sm",
                className="param-input",
            )
        elif spec.kind == "check":
            widget = dbc.Switch(
                id={"type": "corr-param", "name": spec.name},
                label="",
                value=bool(spec.default),
                className="param-check ms-1",
            )
        else:
            widget = dbc.Input(
                id={"type": "corr-param", "name": spec.name},
                type="text",
                value=str(spec.default),
                size="sm",
                className="param-input",
                debounce=True,
            )
        rows.append(html.Div([lbl, widget], className="param-row"))
    return html.Div(rows, className="params-list")


# ── Callback registration ─────────────────────────────────────────────────────


def register_correction(app) -> None:

    # ── T1. Clientside: switch active tab (store + button classes) ─────────
    app.clientside_callback(
        """
        function(n0, n1, n2, n3, current) {
            var slugs = ['ba', 'rho-phi', 'overlay', 'diff'];
            var ctx   = window.dash_clientside.callback_context;
            if (!ctx || !ctx.triggered || ctx.triggered.length === 0)
                return Array(4).fill(window.dash_clientside.no_update).concat([current]);
            var prop_id = ctx.triggered[0].prop_id;
            var active  = current || 'ba';
            for (var i = 0; i < slugs.length; i++) {
                if (prop_id.indexOf('corr-tab-btn-' + slugs[i]) !== -1) {
                    active = slugs[i];
                    break;
                }
            }
            var classes = slugs.map(function(s) {
                return 'corr-tab-btn' + (s === active ? ' active' : '');
            });
            return [classes[0], classes[1], classes[2], classes[3], active];
        }
        """,
        Output("corr-tab-btn-ba", "className"),
        Output("corr-tab-btn-rho-phi", "className"),
        Output("corr-tab-btn-overlay", "className"),
        Output("corr-tab-btn-diff", "className"),
        Output(IDs.CORR_ACTIVE_TAB, "data"),
        Input("corr-tab-btn-ba", "n_clicks"),
        Input("corr-tab-btn-rho-phi", "n_clicks"),
        Input("corr-tab-btn-overlay", "n_clicks"),
        Input("corr-tab-btn-diff", "n_clicks"),
        State(IDs.CORR_ACTIVE_TAB, "data"),
        prevent_initial_call=True,
    )

    # ── T2. Clientside: show/hide panels driven by store ───────────────────
    app.clientside_callback(
        """
        function(active) {
            var a = active || 'ba';
            var styles = {
                'ba':      {display: 'flex'},
                'rho-phi': {display: 'flex', flexDirection: 'column'},
                'overlay': {display: 'flex', flexDirection: 'column'},
                'diff':    {display: 'flex', flexDirection: 'column'}
            };
            var none = {display: 'none'};
            return ['ba','rho-phi','overlay','diff'].map(function(s){
                return s === a ? styles[s] : none;
            });
        }
        """,
        Output("corr-panel-ba", "style"),
        Output("corr-panel-rho-phi", "style"),
        Output("corr-panel-overlay", "style"),
        Output("corr-panel-diff", "style"),
        Input(IDs.CORR_ACTIVE_TAB, "data"),
        prevent_initial_call=False,
    )

    # ── L1. Populate line selector from loaded data ────────────────────────
    @app.callback(
        Output(IDs.CORR_LINE_SEL, "options"),
        Input(IDs.STORE_DATA, "data"),
    )
    def populate_corr_line_opts(store_data):
        records = (store_data or {}).get("station_records", [])
        lines = sorted({r.get("Line", "") for r in records if r.get("Line")})
        return [{"label": ln, "value": ln} for ln in lines]

    # ── L2. Populate station selector based on selected lines ──────────────
    @app.callback(
        Output(IDs.CORR_STN_SEL, "options"),
        Input(IDs.CORR_LINE_SEL, "value"),
        State(IDs.STORE_DATA, "data"),
    )
    def populate_corr_stn_opts(sel_lines, store_data):
        records = (store_data or {}).get("station_records", [])
        if sel_lines:
            records = [r for r in records if r.get("Line") in sel_lines]
        stations = sorted({r.get("ID", "") for r in records if r.get("ID")})
        return [{"label": s, "value": s} for s in stations]

    # ── 1. Sync method dropdown when category changes ──────────────────────
    @app.callback(
        Output(IDs.CORR_METHOD, "options"),
        Output(IDs.CORR_METHOD, "value"),
        Input(IDs.CORR_CAT, "value"),
    )
    def sync_methods(cat):
        if not cat:
            raise PreventUpdate
        methods = CATALOGUE.get(cat, {})
        opts = [{"label": k, "value": k} for k in methods]
        return opts, (opts[0]["value"] if opts else None)

    # ── 2. Sync params panel + description when cat/method changes ─────────
    @app.callback(
        Output(IDs.CORR_PARAMS_PANEL, "children"),
        Output(IDs.CORR_DESC, "children"),
        Input(IDs.CORR_CAT, "value"),
        Input(IDs.CORR_METHOD, "value"),
    )
    def sync_params_panel(cat, method):
        desc = _get_desc(cat or "", method or "")
        panel = _params_panel(cat or "", method or "")
        desc_el = (
            html.Div(desc, className="corr-desc-text") if desc else html.Div()
        )
        return panel, desc_el

    # ── 3. Auto-render: fires on nav, data load, store change, tab switch ──
    @app.callback(
        Output(IDs.IMG_CORR_LEFT, "src"),
        Output(IDs.IMG_CORR_RIGHT, "src"),
        Output(IDs.IMG_CORR_OVERLAY, "src"),
        Output(IDs.IMG_CORR_DIFF, "src"),
        Output(IDs.IMG_CORR_RHO_PHI, "src"),
        Output(IDs.CORR_FEEDBACK, "children"),
        Output("corr-ctx-bar", "children"),
        Output("corr-label-left", "children"),
        Output("corr-label-right", "children"),
        Output("corr-rho-phi-hint", "children"),
        Output("corr-overlay-hint", "children"),
        Input(IDs.NAV_SECTION, "data"),
        Input(IDs.STORE_DATA, "data"),
        Input(IDs.CORR_STORE, "data"),
        Input(IDs.CORR_ACTIVE_TAB, "data"),
        Input(IDs.CORR_RHO_PHI_MODE, "value"),
        Input(IDs.CORR_RHO_PHI_STN, "value"),
        Input(IDs.CORR_RHO_PHI_STYLE, "value"),
        Input(IDs.CORR_LINE_SEL, "value"),
        Input(IDs.CORR_STN_SEL, "value"),
        Input(IDs.CORR_COMP_SEL, "value"),
        Input(IDs.CORR_OVERLAY_STYLE, "value"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        prevent_initial_call=True,
    )
    def auto_render(
        nav_section,
        store_data,
        corr_store,
        active_tab,
        rho_mode,
        rho_max_stn,
        rho_style,
        sel_lines,
        sel_stations,
        sel_comp,
        ovl_style,
        theme,
        session_id,
        active_lines_store,
    ):
        if nav_section != "correction":
            return (no_update,) * 11

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        raw_sites, err = _apply_line_filter(
            session_id, active_lines_store, store_data
        )
        if err == "muted":
            w = no_active_lines_src(dark)
            msg = "All lines muted — enable at least one line."
            return (
                w,
                w,
                w,
                w,
                w,
                msg,
                [msg],
                no_update,
                no_update,
                no_update,
                no_update,
            )

        steps = (corr_store or {}).get("steps", [])
        ctrl = _build_ctrl(session_id, steps, dark, sites_override=raw_sites)

        if ctrl is None:
            hint = [
                html.I(className="bi bi-info-circle me-1"),
                "Load survey data first to use corrections.",
            ]
            return (
                no_update,
                no_update,
                no_update,
                no_update,
                no_update,
                html.Div(hint, className="text-muted small"),
                hint,
                no_update,
                no_update,
                no_update,
                no_update,
            )

        raw = ctrl.raw_sites
        cur = ctrl.current_sites
        n = len(steps)
        tab = active_tab or "ba"

        ctx_children = _ctx_bar_content(steps)
        lbl_left = [html.I(className="bi bi-camera me-1"), "Raw / Original"]
        lbl_right = [
            html.I(className="bi bi-stars me-1"),
            (
                f"After {n} correction{'s' if n != 1 else ''}"
                if n
                else "Corrected (no chain yet)"
            ),
        ]

        style = rho_style or "grid"
        comp = sel_comp or "both"
        rho_hint = _rho_phi_hint(steps, style)

        left_src = right_src = ovl_src = diff_src = rho_phi_src = no_update
        fb = (
            f"{n} correction{'s' if n != 1 else ''} in chain"
            if n
            else "No corrections applied — showing raw data"
        )

        n_max = max(1, int(rho_max_stn or 6))
        mode = rho_mode or "both"
        records = (store_data or {}).get("station_records", [])

        _OVL_HINTS = {
            "median-band": (
                "Median ± IQR: statistical summary across all stations — "
                "gray=raw, color=corrected. Shaded region = IQR."
            ),
            "spaghetti": (
                "Spaghetti + median: thin curves = individual stations, "
                "thick dashed=raw median, thick solid=corrected median."
            ),
            "per-line": (
                "Per-line panels: one column per survey line — "
                "median ± IQR before (gray dashed) vs after (colored) correction."
            ),
        }
        ovl_style_val = ovl_style or "median-band"
        ovl_hint = _OVL_HINTS.get(ovl_style_val, _OVL_HINTS["median-band"])

        if tab == "ba":
            left_src = _render_ps(ctrl, raw, "Raw (original)", dark)
            right_src = _render_ps(
                ctrl,
                cur,
                f"After {n} correction{'s' if n != 1 else ''}"
                if n
                else "No corrections applied",
                dark,
            )
        elif tab == "rho-phi":
            rho_raw = _filter_display_sites(
                raw, sel_lines, sel_stations, records
            )
            rho_cur = _filter_display_sites(
                cur, sel_lines, sel_stations, records
            )
            n_shown = len(list(_iter_edi(rho_raw)))
            if (sel_lines or sel_stations) and n_shown:
                sel_desc = []
                if sel_lines:
                    sel_desc.append(
                        f"{len(sel_lines)} line{'s' if len(sel_lines) != 1 else ''}"
                    )
                if sel_stations:
                    sel_desc.append(
                        f"{len(sel_stations)} station{'s' if len(sel_stations) != 1 else ''}"
                    )
                rho_hint = f"Showing {n_shown} stations from {', '.join(sel_desc)}. {rho_hint}"
            kw = dict(mode=mode, n_max=n_max, dark=dark, comp=comp)
            if style == "pairs":
                rho_phi_src = _render_paired_stations(
                    ctrl, rho_raw, rho_cur, **kw
                )
            elif style == "freq-bands":
                rho_phi_src = _render_freq_bands(ctrl, rho_raw, rho_cur, **kw)
            else:
                rho_phi_src = _render_rho_phi_grid(
                    ctrl, rho_raw, rho_cur, **kw
                )
        elif tab == "overlay":
            ovl_src = _render_overlay(
                ctrl,
                raw,
                cur,
                dark,
                ovl_style=ovl_style_val,
                comp=comp,
                records=records,
            )
        elif tab == "diff":
            diff_src = _render_diff(ctrl, raw, cur, dark)

        return (
            left_src,
            right_src,
            ovl_src,
            diff_src,
            rho_phi_src,
            fb,
            ctx_children,
            lbl_left,
            lbl_right,
            rho_hint,
            ovl_hint,
        )

    # ── 4. Preview — shows current-state vs preview-result ────────────────
    @app.callback(
        Output(IDs.IMG_CORR_LEFT, "src", allow_duplicate=True),
        Output(IDs.IMG_CORR_RIGHT, "src", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.CORR_SPINNER, "children"),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Output("corr-label-left", "children", allow_duplicate=True),
        Output("corr-label-right", "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_PREVIEW, "n_clicks"),
        State(IDs.CORR_CAT, "value"),
        State(IDs.CORR_METHOD, "value"),
        State({"type": "corr-param", "name": ALL}, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def preview_correction(
        n,
        cat,
        method,
        param_vals,
        corr_store,
        theme,
        session_id,
        active_lines_store,
        store_data,
    ):
        if not n:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        raw_sites, err = _apply_line_filter(
            session_id, active_lines_store, store_data
        )
        if err == "muted":
            w = no_active_lines_src(dark)
            return (
                w,
                w,
                "All lines muted.",
                "",
                False,
                "",
                no_update,
                no_update,
            )

        if raw_sites is None and cache_get(session_id) is None:
            e = empty_src(dark=dark)
            return (
                e,
                e,
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
                no_update,
                no_update,
            )

        steps = (corr_store or {}).get("steps", [])
        ctrl = _build_ctrl(session_id, steps, dark, sites_override=raw_sites)

        if ctrl is None:
            e = empty_src(dark=dark)
            return (
                e,
                e,
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
                no_update,
                no_update,
            )

        fn_name = _get_fn_name(cat or "", method or "")
        if not fn_name:
            return (
                no_update,
                no_update,
                f"Unknown correction: {method}",
                "",
                True,
                f"No function for: {method}",
                no_update,
                no_update,
            )

        kwargs = _build_kwargs(cat or "", method or "", param_vals)

        try:
            n_chain = len(steps)
            before = ctrl.current_sites
            result = ctrl.preview(fn_name, kwargs)
            after = result if result is not None else before

            # Labels make the distinction explicit
            lbl_left = [
                html.I(className="bi bi-camera me-1"),
                (
                    f"Current ({n_chain} step{'s' if n_chain != 1 else ''} applied)"
                    if n_chain
                    else "Raw (no corrections yet)"
                ),
            ]
            lbl_right = [
                html.I(
                    className="bi bi-eye me-1",
                    style={"color": "var(--yellow)"},
                ),
                html.Span(
                    "Preview: ",
                    style={"color": "var(--yellow)", "fontWeight": "600"},
                ),
                method,
            ]

            left_title = (
                "Before preview"
                if n_chain == 0
                else f"Current ({n_chain} step{'s' if n_chain != 1 else ''})"
            )
            left_src = _render_ps(ctrl, before, left_title, dark)
            right_src = _render_ps(ctrl, after, f"Preview: {method}", dark)
            fb = f"Preview of '{method}' — click Apply to commit"
            return (
                left_src,
                right_src,
                fb,
                "",
                False,
                "",
                lbl_left,
                lbl_right,
            )

        except Exception as exc:
            e = empty_src(dark=dark)
            return (
                e,
                e,
                f"Preview failed: {exc}",
                "",
                True,
                str(exc),
                no_update,
                no_update,
            )

    # ── 5. Apply — add step to CORR_STORE ─────────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data"),
        Output(IDs.CORR_STACK, "children"),
        Output("corr-chain-badge", "children"),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_APPLY, "n_clicks"),
        State(IDs.CORR_CAT, "value"),
        State(IDs.CORR_METHOD, "value"),
        State({"type": "corr-param", "name": ALL}, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def apply_correction(n, cat, method, param_vals, corr_store, session_id):
        if not n:
            raise PreventUpdate
        if cache_get(session_id) is None:
            return (
                no_update,
                no_update,
                no_update,
                "No data loaded.",
                True,
                "Load survey data first.",
            )
        fn_name = _get_fn_name(cat or "", method or "")
        if not fn_name:
            return (
                no_update,
                no_update,
                no_update,
                f"Unknown: {method}",
                True,
                f"No function for: {method}",
            )
        kwargs = _build_kwargs(cat or "", method or "", param_vals)
        steps = list((corr_store or {}).get("steps", []))
        try:
            ctrl = _build_ctrl(session_id, steps, dark=True)
            ctrl.preview(fn_name, kwargs)
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                f"Validation failed: {exc}",
                True,
                str(exc),
            )
        steps.append({"fn_name": fn_name, "kwargs": kwargs, "label": method})
        n_s = len(steps)
        return (
            {"steps": steps},
            _stack_el(steps),
            _chain_badge(n_s),
            f"Applied: {method}  ({n_s} step{'s' if n_s != 1 else ''} total)",
            False,
            "",
        )

    # ── 6. Drop Bad Frequencies — quick-apply ─────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_DROP_FREQ, "n_clicks"),
        State(IDs.CORR_SNR_THRESH, "value"),
        State(IDs.CORR_MIN_FRAC, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def apply_drop_freq(n, snr_thresh, min_frac, corr_store, session_id):
        if not n:
            raise PreventUpdate
        if cache_get(session_id) is None:
            return (
                no_update,
                no_update,
                no_update,
                "No data loaded.",
                True,
                "Load survey data first.",
            )

        snr_v = float(snr_thresh or 2.5)
        frac_v = float(min_frac or 0.4)
        fn_name = "_strat_freq_filter"
        kwargs = {
            "fmin": 0.0,
            "fmax": 0.0,
            "snr_thresh": snr_v,
            "min_frac": frac_v,
        }
        label = f"Drop bad freq (SNR<{snr_v}, frac≥{frac_v})"

        steps = list((corr_store or {}).get("steps", []))
        try:
            ctrl = _build_ctrl(session_id, steps, dark=True)
            ctrl.preview(fn_name, kwargs)
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                f"Validation failed: {exc}",
                True,
                str(exc),
            )

        steps.append({"fn_name": fn_name, "kwargs": kwargs, "label": label})
        n_s = len(steps)
        return (
            {"steps": steps},
            _stack_el(steps),
            _chain_badge(n_s),
            f"Applied: {label}  ({n_s} step{'s' if n_s != 1 else ''} total)",
            False,
            "",
        )

    # ── 6b. Smooth ρ/φ — Preview ──────────────────────────────────────────
    @app.callback(
        Output(IDs.IMG_CORR_RHO_PHI, "src", allow_duplicate=True),
        Output(IDs.CORR_ACTIVE_TAB, "data", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.CORR_SPINNER, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_SMOOTH_PREVIEW, "n_clicks"),
        State(IDs.CORR_SMOOTH_DEGREE, "value"),
        State(IDs.CORR_SMOOTH_BLEND, "value"),
        State(IDs.CORR_SMOOTH_COMP, "value"),
        State(IDs.CORR_SMOOTH_ROBUST, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.CORR_RHO_PHI_STN, "value"),
        State(IDs.CORR_COMP_SEL, "value"),
        prevent_initial_call=True,
    )
    def smooth_preview(
        n,
        degree,
        blend,
        comp_val,
        robust_val,
        corr_store,
        theme,
        session_id,
        active_lines_store,
        store_data,
        rho_max_stn,
        sel_comp,
    ):
        if not n:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        raw_sites, err = _apply_line_filter(
            session_id, active_lines_store, store_data
        )
        if err == "muted":
            return (
                empty_src(dark=dark),
                "rho-phi",
                "All lines muted.",
                "",
                False,
                "",
            )

        if raw_sites is None and cache_get(session_id) is None:
            return (
                empty_src(dark=dark),
                "rho-phi",
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
            )

        steps = (corr_store or {}).get("steps", [])
        ctrl = _build_ctrl(session_id, steps, dark, sites_override=raw_sites)
        if ctrl is None:
            return (
                empty_src(dark=dark),
                "rho-phi",
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
            )

        kwargs = {
            "components": comp_val or "offdiag",
            "degree": max(0, int(degree or 3)),
            "blend": float(np.clip(float(blend or 1.0), 0.0, 1.0)),
            "smooth_rho": True,
            "smooth_phase": True,
            "robust": bool(robust_val),
        }
        try:
            before = ctrl.current_sites
            smoothed = ctrl.preview("smooth_rho_phase", kwargs)
            if smoothed is None:
                smoothed = before

            comp = sel_comp or "both"
            n_max = max(1, int(rho_max_stn or 6))
            # Show raw vs smoothed in ρ/φ grid (mode="both" = both overlaid)
            src = _render_rho_phi_grid(
                ctrl,
                before,
                smoothed,
                mode="both",
                n_max=n_max,
                dark=dark,
                comp=comp,
            )
            fb = (
                f"Preview: smooth ρ/φ (degree={kwargs['degree']}, "
                f"blend={kwargs['blend']:.2f}, comp={kwargs['components']}) "
                "— click Apply to commit"
            )
            return src, "rho-phi", fb, "", False, ""
        except Exception as exc:
            return (
                empty_src(dark=dark),
                "rho-phi",
                f"Smooth preview failed: {exc}",
                "",
                True,
                str(exc),
            )

    # ── 6c. Smooth ρ/φ — Apply to chain ───────────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_SMOOTH_APPLY, "n_clicks"),
        State(IDs.CORR_SMOOTH_DEGREE, "value"),
        State(IDs.CORR_SMOOTH_BLEND, "value"),
        State(IDs.CORR_SMOOTH_COMP, "value"),
        State(IDs.CORR_SMOOTH_ROBUST, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def smooth_apply(
        n, degree, blend, comp_val, robust_val, corr_store, session_id
    ):
        if not n:
            raise PreventUpdate
        if cache_get(session_id) is None:
            return (
                no_update,
                no_update,
                no_update,
                "No data loaded.",
                True,
                "Load survey data first.",
            )

        kwargs = {
            "components": comp_val or "offdiag",
            "degree": max(0, int(degree or 3)),
            "blend": float(np.clip(float(blend or 1.0), 0.0, 1.0)),
            "smooth_rho": True,
            "smooth_phase": True,
            "robust": bool(robust_val),
        }
        label = (
            f"Smooth ρ/φ (deg={kwargs['degree']}, "
            f"blend={kwargs['blend']:.2f}, {kwargs['components']})"
        )
        steps = list((corr_store or {}).get("steps", []))
        try:
            ctrl = _build_ctrl(session_id, steps, dark=True)
            ctrl.preview("smooth_rho_phase", kwargs)
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                f"Validation failed: {exc}",
                True,
                str(exc),
            )

        steps.append(
            {"fn_name": "smooth_rho_phase", "kwargs": kwargs, "label": label}
        )
        n_s = len(steps)
        return (
            {"steps": steps},
            _stack_el(steps),
            _chain_badge(n_s),
            f"Applied: {label}  ({n_s} step{'s' if n_s != 1 else ''} total)",
            False,
            "",
        )

    # ── 6d. Manual freq drop — Load frequency list ────────────────────────
    @app.callback(
        Output(IDs.CORR_MFREQ_LIST, "options"),
        Output(IDs.CORR_MFREQ_LIST, "value"),
        Output(IDs.CORR_MFREQ_STATUS, "children"),
        Input(IDs.BTN_CORR_MFREQ_LOAD, "n_clicks"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_DATA, "data"),
        prevent_initial_call=True,
    )
    def mfreq_load(n, session_id, active_lines_store, store_data):
        if not n:
            raise PreventUpdate
        raw_sites, err = _apply_line_filter(
            session_id, active_lines_store, store_data
        )
        sites = raw_sites if raw_sites is not None else cache_get(session_id)
        if sites is None:
            return [], [], "No data loaded — load a survey first."

        # Collect all unique frequencies across all stations
        all_freqs: set[float] = set()
        for ed in _iter_edi(sites):
            z = getattr(ed, "Z", None)
            if z is None:
                continue
            fr = getattr(z, "freq", None)
            if fr is None:
                continue
            for f in np.asarray(fr, float).ravel():
                if np.isfinite(f) and f > 0:
                    all_freqs.add(float(f))

        if not all_freqs:
            return [], [], "No valid frequencies found in loaded data."

        # Sort descending (high freq = short period first)
        sorted_freqs = sorted(all_freqs, reverse=True)
        opts = []
        for f in sorted_freqs:
            T = 1.0 / f
            if T < 0.001:
                t_str = f"{T * 1000:.3f} ms"
            elif T < 1.0:
                t_str = f"{T * 1000:.2f} ms"
            elif T < 1000.0:
                t_str = f"{T:.4f} s"
            else:
                t_str = f"{T:.1f} s"
            label = f"T = {t_str}  |  f = {f:.4g} Hz"
            opts.append({"label": label, "value": round(f, 8)})

        status = f"{len(opts)} frequencies loaded — check the ones to remove."
        return opts, [], status

    # ── 6e. Manual freq drop — Select All ─────────────────────────────────
    @app.callback(
        Output(IDs.CORR_MFREQ_LIST, "value", allow_duplicate=True),
        Input(IDs.BTN_CORR_MFREQ_ALL, "n_clicks"),
        State(IDs.CORR_MFREQ_LIST, "options"),
        prevent_initial_call=True,
    )
    def mfreq_sel_all(n, opts):
        if not n or not opts:
            raise PreventUpdate
        return [o["value"] for o in opts]

    # ── 6f. Manual freq drop — Select None ────────────────────────────────
    @app.callback(
        Output(IDs.CORR_MFREQ_LIST, "value", allow_duplicate=True),
        Input(IDs.BTN_CORR_MFREQ_NONE, "n_clicks"),
        prevent_initial_call=True,
    )
    def mfreq_sel_none(n):
        if not n:
            raise PreventUpdate
        return []

    # ── 6g. Manual freq drop — Preview ────────────────────────────────────
    @app.callback(
        Output(IDs.IMG_CORR_RHO_PHI, "src", allow_duplicate=True),
        Output(IDs.CORR_ACTIVE_TAB, "data", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.CORR_SPINNER, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_MFREQ_PREVIEW, "n_clicks"),
        State(IDs.CORR_MFREQ_LIST, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.STORE_THEME, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_ACTIVE_LINES, "data"),
        State(IDs.STORE_DATA, "data"),
        State(IDs.CORR_RHO_PHI_STN, "value"),
        State(IDs.CORR_COMP_SEL, "value"),
        prevent_initial_call=True,
    )
    def mfreq_preview(
        n,
        drop_vals,
        corr_store,
        theme,
        session_id,
        active_lines_store,
        store_data,
        rho_max_stn,
        sel_comp,
    ):
        if not n:
            raise PreventUpdate

        dark = (theme or "dark") == "dark"
        if dark:
            apply_web_dark_theme()
        else:
            apply_web_light_theme()

        if not drop_vals:
            return (
                empty_src(dark=dark),
                "rho-phi",
                "No frequencies selected — check at least one to preview.",
                "",
                False,
                "",
            )

        raw_sites, err = _apply_line_filter(
            session_id, active_lines_store, store_data
        )
        if err == "muted":
            return (
                empty_src(dark=dark),
                "rho-phi",
                "All lines muted.",
                "",
                False,
                "",
            )
        if raw_sites is None and cache_get(session_id) is None:
            return (
                empty_src(dark=dark),
                "rho-phi",
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
            )

        steps = (corr_store or {}).get("steps", [])
        ctrl = _build_ctrl(session_id, steps, dark, sites_override=raw_sites)
        if ctrl is None:
            return (
                empty_src(dark=dark),
                "rho-phi",
                "No data loaded.",
                "",
                True,
                "Load survey data first.",
            )

        kwargs = {"drop_freqs": list(drop_vals)}
        try:
            before = ctrl.current_sites
            dropped = ctrl.preview("drop_freqs_manual", kwargs)
            if dropped is None:
                dropped = before

            comp = sel_comp or "both"
            n_max = max(1, int(rho_max_stn or 6))
            # Use freq-bands style to highlight which periods were removed
            src = _render_freq_bands(
                ctrl,
                before,
                dropped,
                mode="both",
                n_max=n_max,
                dark=dark,
                comp=comp,
                dropped_freqs=list(drop_vals),
            )
            n_drop = len(drop_vals)
            fb = (
                f"Preview: {n_drop} frequenc{'ies' if n_drop != 1 else 'y'} will be removed "
                "— red bands mark dropped periods. Click 'Apply Drop' to commit."
            )
            return src, "rho-phi", fb, "", False, ""
        except Exception as exc:
            return (
                empty_src(dark=dark),
                "rho-phi",
                f"Manual freq preview failed: {exc}",
                "",
                True,
                str(exc),
            )

    # ── 6h. Manual freq drop — Apply to chain ─────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Output(IDs.TOAST_ERROR, "is_open", allow_duplicate=True),
        Output(IDs.TOAST_BODY, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_MFREQ_APPLY, "n_clicks"),
        State(IDs.CORR_MFREQ_LIST, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_ID, "data"),
        prevent_initial_call=True,
    )
    def mfreq_apply(n, drop_vals, corr_store, session_id):
        if not n:
            raise PreventUpdate
        if not drop_vals:
            return (
                no_update,
                no_update,
                no_update,
                "No frequencies selected — nothing to apply.",
                False,
                "",
            )
        if cache_get(session_id) is None:
            return (
                no_update,
                no_update,
                no_update,
                "No data loaded.",
                True,
                "Load survey data first.",
            )

        kwargs = {"drop_freqs": list(drop_vals)}
        n_drop = len(drop_vals)
        label = f"Manual freq drop ({n_drop} frequenc{'ies' if n_drop != 1 else 'y'})"
        steps = list((corr_store or {}).get("steps", []))
        try:
            ctrl = _build_ctrl(session_id, steps, dark=True)
            ctrl.preview("drop_freqs_manual", kwargs)
        except Exception as exc:
            return (
                no_update,
                no_update,
                no_update,
                f"Validation failed: {exc}",
                True,
                str(exc),
            )

        steps.append(
            {"fn_name": "drop_freqs_manual", "kwargs": kwargs, "label": label}
        )
        n_s = len(steps)
        return (
            {"steps": steps},
            _stack_el(steps),
            _chain_badge(n_s),
            f"Applied: {label}  ({n_s} step{'s' if n_s != 1 else ''} total)",
            False,
            "",
        )

    # ── 7. Undo ────────────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_UNDO, "n_clicks"),
        State(IDs.CORR_STORE, "data"),
        prevent_initial_call=True,
    )
    def undo_correction(n, corr_store):
        if not n:
            raise PreventUpdate
        steps = list((corr_store or {}).get("steps", []))
        if not steps:
            return no_update, no_update, no_update, "Nothing to undo."
        removed = steps.pop()
        n_left = len(steps)
        return (
            {"steps": steps},
            _stack_el(steps),
            _chain_badge(n_left),
            f"Undone '{removed.get('label', '')}' — "
            f"{n_left} step{'s' if n_left != 1 else ''} remaining",
        )

    # ── 8. Reset ───────────────────────────────────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Input(IDs.BTN_CORR_RESET, "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_corrections(n):
        if not n:
            raise PreventUpdate
        return (
            {"steps": []},
            html.Span(
                "No corrections applied yet.", className="text-muted small"
            ),
            "",
            "Reset to raw — all corrections cleared.",
        )

    # ── 9. Delete individual step from chain ───────────────────────────────
    @app.callback(
        Output(IDs.CORR_STORE, "data", allow_duplicate=True),
        Output(IDs.CORR_STACK, "children", allow_duplicate=True),
        Output("corr-chain-badge", "children", allow_duplicate=True),
        Output(IDs.CORR_FEEDBACK, "children", allow_duplicate=True),
        Input({"type": "corr-step-del", "index": ALL}, "n_clicks"),
        State(IDs.CORR_STORE, "data"),
        prevent_initial_call=True,
    )
    def delete_step(n_clicks_list, corr_store):
        triggered = ctx.triggered_id
        if not triggered or not isinstance(triggered, dict):
            raise PreventUpdate
        idx = triggered.get("index")
        if idx is None:
            raise PreventUpdate
        vals = ctx.inputs_list[0]
        clicked = [item for item in vals if item.get("value", 0) > 0]
        if not clicked:
            raise PreventUpdate

        steps = list((corr_store or {}).get("steps", []))
        if 0 <= idx < len(steps):
            removed = steps.pop(idx)
            fb = f"Removed step {idx + 1}: '{removed.get('label', '')}'"
        else:
            raise PreventUpdate

        n_left = len(steps)
        return ({"steps": steps}, _stack_el(steps), _chain_badge(n_left), fb)

    # ── 10. Difference stats strip ─────────────────────────────────────────
    @app.callback(
        Output(IDs.CORR_DIFF_STATS, "children"),
        Input(IDs.CORR_STORE, "data"),
        Input(IDs.CORR_ACTIVE_TAB, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def update_diff_stats(corr_store, active_tab, session_id, theme):
        if active_tab != "diff":
            return no_update
        dark = (theme or "dark") == "dark"
        steps = (corr_store or {}).get("steps", [])
        if not steps:
            return html.Span(
                "No corrections applied yet.", className="text-muted small"
            )
        try:
            ctrl = _build_ctrl(session_id, steps, dark=dark)
            if ctrl is None:
                return no_update
            rho_raw = _collect_rho(ctrl.raw_sites)
            rho_cur = _collect_rho(ctrl.current_sites)
            if rho_raw is None or rho_cur is None or rho_raw.size == 0:
                return html.Span(
                    "No data for stats.", className="text-muted small"
                )
            delta = np.log10(rho_cur + 1e-12) - np.log10(rho_raw + 1e-12)
            mn, mx = float(delta.min()), float(delta.max())
            mu, sd = float(delta.mean()), float(delta.std())
            n_neg = int((delta < 0).sum())
            n_pos = int((delta > 0).sum())
            return _diff_stats_strip(mn, mx, mu, sd, n_neg, n_pos)
        except Exception:
            return no_update

    # ── 11. Export corrected data ──────────────────────────────────────────
    @app.callback(
        Output(IDs.CORR_DOWNLOAD, "data"),
        Input(IDs.BTN_CORR_EXPORT, "n_clicks"),
        State(IDs.CORR_EXPORT_FMT, "value"),
        State(IDs.CORR_STORE, "data"),
        State(IDs.SESSION_ID, "data"),
        State(IDs.STORE_THEME, "data"),
        prevent_initial_call=True,
    )
    def export_corrected(n, fmt, corr_store, session_id, theme):
        if not n:
            raise PreventUpdate
        dark = (theme or "dark") == "dark"
        steps = (corr_store or {}).get("steps", [])
        ctrl = _build_ctrl(session_id, steps, dark=dark)
        if ctrl is None:
            raise PreventUpdate
        try:
            import pandas as pd

            sites = ctrl.current_sites
            rows = []
            for ed in _iter_edi(sites):
                rows.extend(_edi_to_rows(ed))
            if not rows:
                raise PreventUpdate
            df = pd.DataFrame(rows)
            if fmt == "csv":
                return dcc.send_data_frame(
                    df.to_csv, "corrected_data.csv", index=False
                )
            elif fmt == "xyz":
                cols = [
                    c
                    for c in ["station", "lon", "lat", "rho_xy", "phase_xy"]
                    if c in df.columns
                ]
                return dcc.send_data_frame(
                    df[cols].to_csv,
                    "corrected_profile.xyz",
                    index=False,
                    sep="\t",
                )
            else:
                return dcc.send_data_frame(
                    df.to_csv, "corrected_edi.csv", index=False
                )
        except Exception:
            raise PreventUpdate
