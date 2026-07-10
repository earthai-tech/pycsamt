# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Visualisation helpers for :class:`~pycsamt.ai.processing.qc.EMQCScorer` output.

Six public functions cover the main inspection use-cases:

* :func:`plot_qc_scores`              — per-station bar chart, multi-profile
* :func:`plot_qc_heatmap`             — station × frequency quality heat-map
* :func:`plot_qc_feature_heatmap`     — one panel per QC feature
* :func:`plot_qc_score_distribution`  — score histogram + KDE per profile
* :func:`plot_qc_score_spread`        — violin / box / strip per profile
* :func:`plot_qc_summary`             — 3-panel summary figure
"""
from __future__ import annotations

from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..plot._style import (
    EMStyle,
    StationTickConfig,
)

__all__ = [
    "plot_qc_scores",
    "plot_qc_heatmap",
    "plot_qc_feature_heatmap",
    "plot_qc_score_distribution",
    "plot_qc_score_spread",
    "plot_qc_summary",
]

# ── default profile colour cycle ──────────────────────────────────────────── #
_PROFILE_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c",
    "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22",
]

_FEATURE_LABELS: dict[str, str] = {
    "snr":        "SNR",
    "swift_skew": "Swift skew |β|",
    "asym":       r"Asymmetry  log$_{10}$(|Z$_{xy}$|/|Z$_{yx}$|)",
    "phase_xy":   "Phase XY (°)",
    "phase_yx":   "Phase YX (°)",
    "score":      "QC score",
}


# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_profile_colors(
    profile_names: list[str],
    user_colors: Any | None,
) -> dict[str, str]:
    if isinstance(user_colors, dict):
        return {p: user_colors.get(p, _PROFILE_COLORS[i % len(_PROFILE_COLORS)])
                for i, p in enumerate(profile_names)}
    if isinstance(user_colors, (list, tuple)):
        return {p: user_colors[i % len(user_colors)]
                for i, p in enumerate(profile_names)}
    return {p: _PROFILE_COLORS[i % len(_PROFILE_COLORS)]
            for i, p in enumerate(profile_names)}


def _normalise_qc_input(
    scores: Any,
    station_labels: list[str] | None,
) -> tuple[dict[str, np.ndarray], dict[str, list[str] | None]]:
    """
    Convert flexible input to ``{profile → 1-D score array}`` + label dicts.

    Accepted formats
    ----------------
    * ``ndarray`` (n_st,)          single profile
    * ``dict[str, ndarray]``       one entry per profile
    * ``pd.DataFrame``             with ``score`` column; optional
                                   ``profile`` and ``station`` columns
    """
    if isinstance(scores, np.ndarray):
        return {"_": scores.ravel().astype(float)}, {"_": station_labels}

    if isinstance(scores, dict):
        pscores = {k: np.asarray(v, dtype=float).ravel()
                   for k, v in scores.items()}
        return pscores, {k: None for k in pscores}

    if isinstance(scores, pd.DataFrame):
        if "score" not in scores.columns:
            raise ValueError(
                "DataFrame must contain a 'score' column "
                "(use EMQCScorer.score_table())."
            )
        has_prof = "profile" in scores.columns
        has_st   = "station" in scores.columns

        if has_prof:
            pscores, plabels = {}, {}
            for prof, grp in scores.groupby("profile"):
                key = str(prof)
                if has_st:
                    agg = grp.groupby("station")["score"].median()
                    pscores[key] = agg.values.astype(float)
                    plabels[key] = list(agg.index.astype(str))
                else:
                    pscores[key] = grp["score"].values.astype(float)
                    plabels[key] = None
        else:
            if has_st:
                agg = scores.groupby("station")["score"].median()
                pscores = {"_": agg.values.astype(float)}
                plabels = {"_": list(agg.index.astype(str))}
            else:
                pscores = {"_": scores["score"].values.astype(float)}
                plabels = {"_": station_labels}
        return pscores, plabels

    raise TypeError(
        f"'scores' must be ndarray, dict, or DataFrame; "
        f"got {type(scores).__name__}."
    )


def _make_tick_config(
    every: int | str,
    rotation: float,
    fontsize: int,
    override: StationTickConfig | None,
) -> StationTickConfig:
    """Build a StationTickConfig from individual params or use *override*."""
    if override is not None:
        return override
    return StationTickConfig(every=every, rotation=rotation, fontsize=fontsize)


def _logT_grid_from_df(
    df: pd.DataFrame,
    station_order: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    freqs = np.sort(df["freq"].unique())[::-1]
    n_f, n_st = freqs.size, len(station_order)
    mat = np.full((n_f, n_st), np.nan)
    st_idx = {s: i for i, s in enumerate(station_order)}
    fr_idx = {f: i for i, f in enumerate(freqs)}
    for _, row in df.iterrows():
        si = st_idx.get(row["station"])
        fi = fr_idx.get(row["freq"])
        if si is not None and fi is not None:
            mat[fi, si] = float(row["score"])
    return mat, freqs, np.arange(n_st, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# Internal bar-chart renderer (reused by plot_qc_scores and plot_qc_summary)
# ─────────────────────────────────────────────────────────────────────────────

def _render_bar_chart(
    ax: plt.Axes,
    pscores: dict[str, np.ndarray],
    plabels: dict[str, list[str] | None],
    pcolors: dict[str, str],
    *,
    score_threshold: float,
    show_scatter: bool,
    scatter_src: pd.DataFrame | None,
    tick_cfg: StationTickConfig,
    bar_width: float = 0.80,
    bar_alpha: float = 0.88,
    show_profile_labels: bool = True,
    profile_label_y: float = 1.04,
    profile_label_fontsize: int = 8,
    separator_color: str = "gray",
    separator_lw: float = 0.7,
    separator_ls: str = "--",
    show_zone_labels: bool = True,
    zone_label_fontsize: int = 8,
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.4,
    threshold_ls: str = "--",
    bad_zone_color: str = "#fde8e8",
    bad_zone_alpha: float = 0.30,
    xlabel: str = "Station",
    show_grid: bool = True,
) -> None:
    """Render bar chart onto *ax* (no figure creation, no return value)."""
    ax.axhspan(0, score_threshold, color=bad_zone_color,
               alpha=bad_zone_alpha, zorder=0)

    profiles = list(pscores.keys())
    x_offset = 0
    all_positions: list[float] = []
    all_labels:    list[str]   = []

    for i, prof in enumerate(profiles):
        sc   = pscores[prof]
        n    = len(sc)
        col  = pcolors[prof]
        labs = plabels[prof]
        xpos = np.arange(x_offset, x_offset + n, dtype=float)

        # ── optional per-frequency scatter ────────────────────────────────
        if show_scatter and scatter_src is not None \
                and "freq" in scatter_src.columns \
                and "station" in scatter_src.columns:
            sub = (scatter_src[scatter_src["profile"] == prof]
                   if "profile" in scatter_src.columns else scatter_src)
            if labs:
                rng = np.random.default_rng(i)
                for si, st in enumerate(labs):
                    rows = sub[sub["station"] == st]["score"].values
                    if rows.size:
                        jit = rng.uniform(-0.18, 0.18, rows.size)
                        ax.scatter(
                            np.full(rows.size, x_offset + si) + jit,
                            rows, color=col, alpha=0.22, s=5, zorder=1,
                        )

        # ── bars ──────────────────────────────────────────────────────────
        sc_plot = np.where(np.isfinite(sc), sc, 0.0)
        ax.bar(xpos, sc_plot, color=col, width=bar_width,
               edgecolor="none", alpha=bar_alpha, zorder=2)

        # thin top-edge line at score value
        for xi, yi in zip(xpos, sc_plot):
            ax.hlines(yi, xi - bar_width * 0.44, xi + bar_width * 0.44,
                      colors=col, lw=0.9, zorder=3, alpha=0.65)

        # ── profile separator ─────────────────────────────────────────────
        if i < len(profiles) - 1:
            ax.axvline(x_offset + n - 0.5,
                       color=separator_color, lw=separator_lw,
                       ls=separator_ls, alpha=0.6, zorder=1)

        # ── profile label above bars ──────────────────────────────────────
        if show_profile_labels and prof != "_":
            ax.text(
                x_offset + (n - 1) / 2.0, profile_label_y, prof,
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=profile_label_fontsize,
                fontweight="bold", color=col,
            )

        # ── accumulate tick info ──────────────────────────────────────────
        all_positions.extend(xpos.tolist())
        if labs is not None:
            all_labels.extend(labs)
        else:
            all_labels.extend([str(j) for j in range(x_offset, x_offset + n)])

        x_offset += n

    # ── threshold line ─────────────────────────────────────────────────────
    ax.axhline(score_threshold, color=threshold_color,
               lw=threshold_lw, ls=threshold_ls, zorder=4)

    # ── zone labels on right margin ────────────────────────────────────────
    if show_zone_labels:
        y_good   = (1.0 + score_threshold) / 2.0
        y_review = score_threshold / 2.0
        kw = dict(transform=ax.transAxes, va="center",
                  fontsize=zone_label_fontsize, fontweight="bold")
        ax.text(1.002, y_good,   "Good ▶",   ha="left",
                color="#2ca02c", **kw)
        ax.text(1.002, y_review, "◀ Review", ha="left",
                color=threshold_color, **kw)

    # ── x-axis via StationTickConfig ───────────────────────────────────────
    pos_arr = np.asarray(all_positions)
    tick_cfg.apply(
        ax, pos_arr, all_labels,
        xlabel=xlabel,
        xlim=(-0.8, x_offset - 0.2),
    )

    ax.set_ylim(0, 1.15)
    ax.set_ylabel("QC score", fontsize=9)
    if show_grid:
        ax.grid(True, axis="y", ls=":", lw=0.4, color="gray",
                alpha=0.5, zorder=0)
    ax.set_axisbelow(True)


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_scores
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_scores(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    score_threshold: float = 0.5,
    # threshold styling
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.4,
    threshold_ls: str = "--",
    # rejection zone
    bad_zone_color: str = "#fde8e8",
    bad_zone_alpha: float = 0.30,
    # bar styling
    bar_width: float = 0.80,
    bar_alpha: float = 0.88,
    # scatter overlay (per-frequency, when DataFrame with freq column)
    show_scatter: bool = True,
    # profile annotations
    show_profile_labels: bool = True,
    profile_label_fontsize: int = 8,
    # separator
    separator_color: str = "gray",
    separator_lw: float = 0.7,
    separator_ls: str = "--",
    # zone labels
    show_zone_labels: bool = True,
    zone_label_fontsize: int = 8,
    # ── station tick control ───────────────────────────────────────────────
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    # ── axes styling ──────────────────────────────────────────────────────
    xlabel: str = "Station",
    ylabel: str = "QC score",
    title: str = "",
    ylim: tuple[float, float] | None = None,
    show_grid: bool = True,
    figsize: tuple[float, float] = (10.0, 4.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Per-station QC-score bar chart, optionally grouped by profile.

    Handles label overlap automatically: ``tick_every="auto"`` (the
    default) picks the smallest *nice* step (1 → 2 → 5 → 10 → 25 → …)
    so that station labels never overlap, regardless of how many stations
    are shown.  Pass a :class:`~pycsamt.ai.plot._style.StationTickConfig`
    via *station_tick_config* for fine-grained reusable control.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Quality scores.  Accepted forms:

        * ``ndarray`` shape ``(n_st,)`` — single profile.
        * ``dict[str, ndarray]`` — one 1-D array per profile name.
        * :class:`pandas.DataFrame` from ``EMQCScorer.score_table()`` —
          requires a ``score`` column; optional ``station``, ``freq``,
          ``profile`` columns are used when present.

    station_labels : list of str or None
        X-axis tick labels.  Inferred from the ``station`` column when
        *scores* is a DataFrame.  Defaults to integer index strings.
    profile_colors : dict, list, or None
        ``{profile_name: color}`` mapping, ordered list, or ``None``
        for the built-in colour cycle.
    score_threshold : float, default ``0.5``
        Threshold line and rejection-zone boundary.
    threshold_color, threshold_lw, threshold_ls : str, float, str
    bad_zone_color, bad_zone_alpha : str, float
    bar_width : float, default ``0.80``
    bar_alpha : float, default ``0.88``
    show_scatter : bool, default ``True``
        When *scores* is a DataFrame with a ``freq`` column, draw
        semi-transparent per-frequency score dots behind each bar.
    show_profile_labels : bool, default ``True``
    profile_label_fontsize : int, default ``8``
    separator_color, separator_lw, separator_ls : str, float, str
    show_zone_labels : bool, default ``True``
        Add "Good ▶" / "◀ Review" margin annotations.
    zone_label_fontsize : int, default ``8``
    tick_every : int or ``"auto"``, default ``"auto"``
        Station tick step.  ``"auto"`` computes the step from the figure
        width and label length so labels never overlap.  Pass an integer
        to force a fixed step, e.g. ``tick_every=5`` shows every 5th label.
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    station_tick_config : :class:`~pycsamt.ai.plot._style.StationTickConfig` or None
        If provided, overrides *tick_every*, *tick_label_rotation*, and
        *tick_fontsize* entirely.
    xlabel, ylabel, title : str
    ylim : (ymin, ymax) or None
    show_grid : bool, default ``True``
    figsize : (w, h), default ``(10.0, 4.2)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`

    Examples
    --------
    >>> from pycsamt.ai.processing.plot import plot_qc_scores
    >>> plot_qc_scores({"L18": sc_18, "L22": sc_22}, tick_every=5)

    >>> from pycsamt.ai.plot import StationTickConfig
    >>> cfg = StationTickConfig(every=10, rotation=30, fontsize=8)
    >>> plot_qc_scores(scores_dict, station_tick_config=cfg)
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles         = list(pscores.keys())
    pcolors          = _resolve_profile_colors(profiles, profile_colors)
    tick_cfg         = _make_tick_config(tick_every, tick_label_rotation,
                                         tick_fontsize, station_tick_config)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _render_bar_chart(
        ax, pscores, plabels, pcolors,
        score_threshold=score_threshold,
        show_scatter=show_scatter,
        scatter_src=(scores if isinstance(scores, pd.DataFrame) else None),
        tick_cfg=tick_cfg,
        bar_width=bar_width,
        bar_alpha=bar_alpha,
        show_profile_labels=show_profile_labels,
        profile_label_fontsize=profile_label_fontsize,
        separator_color=separator_color,
        separator_lw=separator_lw,
        separator_ls=separator_ls,
        show_zone_labels=show_zone_labels,
        zone_label_fontsize=zone_label_fontsize,
        threshold_color=threshold_color,
        threshold_lw=threshold_lw,
        threshold_ls=threshold_ls,
        bad_zone_color=bad_zone_color,
        bad_zone_alpha=bad_zone_alpha,
        xlabel=xlabel,
        show_grid=show_grid,
    )

    ax.set_ylabel(ylabel, fontsize=9)
    if ylim is not None:
        ax.set_ylim(*ylim)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")

    # ── legend ─────────────────────────────────────────────────────────────
    if len(profiles) > 1 or profiles[0] != "_":
        handles = [
            mpatches.Patch(fc=pcolors[p], label=p, alpha=bar_alpha)
            for p in profiles if p != "_"
        ]
        handles.append(
            plt.Line2D([], [], color=threshold_color, ls=threshold_ls,
                       lw=threshold_lw,
                       label=f"Review threshold ({score_threshold:.2f})")
        )
        ax.legend(handles=handles, loc="lower right",
                  fontsize=7.5, framealpha=0.9)

    return ax


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_heatmap
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_heatmap(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    score_threshold: float = 0.5,
    cmap: str = "RdYlGn",
    vmin: float = 0.0,
    vmax: float = 1.0,
    period_up: bool = True,
    show_threshold_contour: bool = True,
    contour_color: str = "white",
    contour_lw: float = 0.8,
    n_yticks: int = 7,
    colorbar_label: str = "QC score",
    xlabel: str = "Station",
    ylabel: str = "Period (s)",
    title: str = "",
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    figsize: tuple[float, float] = (10.0, 5.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Station × frequency QC-score heat-map.

    Parameters
    ----------
    scores : DataFrame
        Must contain ``station``, ``freq``, and ``score`` columns
        (output of :meth:`~pycsamt.ai.processing.qc.EMQCScorer.score_table`).
    station_labels : list of str or None
    score_threshold : float, default ``0.5``
    cmap : str, default ``"RdYlGn"``
    vmin, vmax : float
    period_up : bool, default ``True``
    show_threshold_contour : bool, default ``True``
    contour_color, contour_lw : str, float
    n_yticks : int, default ``7``
    colorbar_label, xlabel, ylabel, title : str
    tick_every : int or ``"auto"``
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    station_tick_config : :class:`~pycsamt.ai.plot._style.StationTickConfig` or None
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not isinstance(scores, pd.DataFrame):
        raise TypeError(
            "plot_qc_heatmap expects a DataFrame with columns "
            "'station', 'freq', 'score' (from EMQCScorer.score_table())."
        )
    if not {"station", "freq", "score"}.issubset(scores.columns):
        raise ValueError("DataFrame must have columns 'station', 'freq', 'score'.")

    tick_cfg = _make_tick_config(tick_every, tick_label_rotation,
                                  tick_fontsize, station_tick_config)

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(scores["station"].unique(),
                          key=lambda x: (str(x).isdigit(), x))

    mat, freqs, x_cents = _logT_grid_from_df(scores, st_order)
    n_f, n_st = mat.shape

    log_T  = np.log10(1.0 / freqs)
    d_lT   = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0]    = log_T[0]   - 0.5 * abs(d_lT[0])  if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else np.array([])
    y_edge[-1]   = log_T[-1]  + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    x_edge = np.concatenate([
        [x_cents[0] - 0.5],
        0.5 * (x_cents[:-1] + x_cents[1:]),
        [x_cents[-1] + 0.5],
    ]) if n_st > 1 else np.array([-0.5, 0.5])

    X, Y = np.meshgrid(x_edge, y_edge)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    qm = ax.pcolormesh(X, Y, mat, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", rasterized=True)
    if period_up:
        ax.invert_yaxis()

    if show_threshold_contour:
        try:
            x_c = 0.5 * (x_edge[:-1] + x_edge[1:])
            y_c = 0.5 * (y_edge[:-1] + y_edge[1:])
            Xc, Yc = np.meshgrid(x_c, y_c)
            ax.contour(Xc, Yc, mat, levels=[score_threshold],
                       colors=[contour_color], linewidths=[contour_lw])
        except Exception:
            pass

    cb = plt.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(colorbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # y-axis (log-period)
    pos  = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(f"$10^{{{r}}}$" if abs(r - v) < 0.04
                    else f"$10^{{{v:.1f}}}$")
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    # x-axis via StationTickConfig
    tick_cfg.apply(ax, x_cents, st_order, xlabel=xlabel,
                   xlim=(x_edge[0], x_edge[-1]))

    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")

    return ax


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_feature_heatmap
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_feature_heatmap(
    df: pd.DataFrame,
    *,
    features: list[str] | None = None,
    station_labels: list[str] | None = None,
    cmaps: dict[str, str] | None = None,
    period_up: bool = True,
    n_yticks: int = 5,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 6,
    station_tick_config: StationTickConfig | None = None,
    clim_pct: tuple[float, float] = (2.0, 98.0),
    title: str = "",
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    One panel per QC feature — station × frequency heat-maps.

    Parameters
    ----------
    df : DataFrame
        Output of ``EMQCScorer.score_table()``.  Must contain
        ``station``, ``freq``, and the feature columns.
    features : list of str or None
        Feature columns to plot.  Defaults to all five QC features.
    station_labels : list of str or None
    cmaps : dict or None
        ``{feature: cmap_name}`` overrides.
    period_up : bool, default ``True``
    n_yticks : int, default ``5``
    tick_every, tick_label_rotation, tick_fontsize : tick control
    station_tick_config : :class:`~pycsamt.ai.plot._style.StationTickConfig` or None
    clim_pct : (lo, hi), default ``(2.0, 98.0)``
    title : str
    figsize : (w, h) or None

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    _default_feats = ["snr", "swift_skew", "asym", "phase_xy", "phase_yx"]
    if features is None:
        features = [f for f in _default_feats if f in df.columns]
    else:
        features = [f for f in features if f in df.columns]
    if not features:
        raise ValueError("No valid feature columns found in DataFrame.")

    tick_cfg = _make_tick_config(tick_every, tick_label_rotation,
                                  tick_fontsize, station_tick_config)

    _cmap_defaults = {
        "snr":        "YlGn",
        "swift_skew": "YlOrRd_r",
        "asym":       "RdBu_r",
        "phase_xy":   "RdBu_r",
        "phase_yx":   "RdBu_r",
        "score":      "RdYlGn",
    }
    cmap_map = dict(_cmap_defaults)
    if cmaps:
        cmap_map.update(cmaps)

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(df["station"].unique(),
                          key=lambda x: (str(x).isdigit(), x))

    n_feat = len(features)
    n_st   = len(st_order)
    freqs  = np.sort(df["freq"].unique())[::-1]
    n_f    = freqs.size

    if figsize is None:
        figsize = (max(8.0, n_st * 0.35), 2.8 * n_feat)

    fig, axes = plt.subplots(n_feat, 1, figsize=figsize, sharex=True,
                              gridspec_kw={"hspace": 0.45})
    if n_feat == 1:
        axes = [axes]

    log_T = np.log10(1.0 / freqs)
    if n_f > 1:
        d = np.diff(log_T)
        ye = np.empty(n_f + 1)
        ye[0] = log_T[0] - 0.5 * abs(d[0])
        ye[1:-1] = log_T[:-1] + 0.5 * d
        ye[-1] = log_T[-1] + 0.5 * abs(d[-1])
    else:
        ye = np.array([log_T[0] - 0.5, log_T[0] + 0.5])

    x_cents = np.arange(n_st, dtype=float)
    xe = np.concatenate([
        [x_cents[0] - 0.5],
        0.5 * (x_cents[:-1] + x_cents[1:]),
        [x_cents[-1] + 0.5],
    ]) if n_st > 1 else np.array([-0.5, 0.5])
    X, Y = np.meshgrid(xe, ye)

    st_idx = {s: i for i, s in enumerate(st_order)}
    fr_idx = {f: i for i, f in enumerate(freqs)}

    for k, feat in enumerate(features):
        ax = axes[k]
        mat = np.full((n_f, n_st), np.nan)
        for _, row in df.iterrows():
            si = st_idx.get(row["station"])
            fi = fr_idx.get(row["freq"])
            if si is not None and fi is not None:
                val = row.get(feat, np.nan)
                if pd.notna(val):
                    mat[fi, si] = float(val)

        fin = mat[np.isfinite(mat)]
        if fin.size:
            vmin = float(np.percentile(fin, clim_pct[0]))
            vmax = float(np.percentile(fin, clim_pct[1]))
        else:
            vmin, vmax = 0.0, 1.0

        if feat in ("asym", "phase_xy", "phase_yx"):
            v = max(abs(vmin), abs(vmax))
            vmin, vmax = -v, v

        qm = ax.pcolormesh(X, Y, mat, cmap=cmap_map.get(feat, "RdBu_r"),
                           vmin=vmin, vmax=vmax, shading="flat", rasterized=True)
        if period_up:
            ax.invert_yaxis()

        cb = plt.colorbar(qm, ax=ax, fraction=0.018, pad=0.01)
        cb.set_label(_FEATURE_LABELS.get(feat, feat), fontsize=7)
        cb.ax.tick_params(labelsize=6)

        pos  = np.linspace(log_T.min(), log_T.max(), n_yticks)
        ylab = []
        for v in pos:
            r = round(v)
            ylab.append(f"$10^{{{r}}}$" if abs(r - v) < 0.04
                        else f"$10^{{{v:.1f}}}$")
        ax.set_yticks(pos)
        ax.set_yticklabels(ylab, fontsize=tick_fontsize)
        ax.set_ylabel("Period (s)", fontsize=7)

    # bottom x-axis via StationTickConfig
    tick_cfg.apply(axes[-1], x_cents, st_order,
                   xlabel="Station",
                   xlim=(xe[0], xe[-1]))

    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_score_distribution   (standalone)
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_score_distribution(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    score_threshold: float = 0.5,
    n_bins: int = 20,
    show_kde: bool = True,
    kde_bw: float | str | None = None,
    fill_alpha: float = 0.35,
    line_lw: float = 1.8,
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.3,
    threshold_ls: str = "--",
    xlabel: str = "QC score",
    ylabel: str = "Density",
    title: str = "",
    show_legend: bool = True,
    figsize: tuple[float, float] = (6.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Score histogram and optional KDE per profile.

    Each profile's scores are shown as a semi-transparent histogram with
    an overlaid KDE curve.  The review threshold is marked as a vertical
    line.  Suitable as a standalone figure or embedded in a summary panel.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Same flexible input as :func:`plot_qc_scores`.
    station_labels : list of str or None
    profile_colors : dict, list, or None
    score_threshold : float, default ``0.5``
    n_bins : int, default ``20``
    show_kde : bool, default ``True``
        Overlay a KDE curve on top of each histogram.
    kde_bw : float, str, or None
        KDE bandwidth passed to ``scipy.stats.gaussian_kde``.
        ``None`` → Scott's rule.
    fill_alpha : float, default ``0.35``
    line_lw : float, default ``1.8``
    threshold_color, threshold_lw, threshold_ls : str, float, str
    xlabel, ylabel, title : str
    show_legend : bool, default ``True``
    figsize : (w, h), default ``(6.0, 4.0)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    pscores, _ = _normalise_qc_input(scores, station_labels)
    profiles   = list(pscores.keys())
    pcolors    = _resolve_profile_colors(profiles, profile_colors)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for prof in profiles:
        sc  = pscores[prof]
        col = pcolors[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        label = prof if prof != "_" else "all"

        ax.hist(fin, bins=n_bins, range=(0, 1),
                color=col, alpha=fill_alpha, edgecolor="none",
                density=True, label=label)

        if show_kde and fin.size > 3:
            try:
                from scipy.stats import gaussian_kde
                bw  = kde_bw if kde_bw is not None else "scott"
                kde = gaussian_kde(fin, bw_method=bw)
                xs  = np.linspace(0, 1, 300)
                ax.plot(xs, kde(xs), color=col, lw=line_lw)
            except ImportError:
                pass

    ax.axvline(score_threshold, color=threshold_color,
               lw=threshold_lw, ls=threshold_ls)

    # threshold band annotation
    ax.axvspan(0, score_threshold, color="#fde8e8", alpha=0.25, zorder=0)

    ax.set_xlim(0, 1)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")
    if show_legend and (len(profiles) > 1 or profiles[0] != "_"):
        ax.legend(fontsize=7.5, framealpha=0.85)

    return ax


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_score_spread   (standalone)
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_score_spread(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    score_threshold: float = 0.5,
    kind: str = "violin",
    show_scatter: bool = True,
    scatter_alpha: float = 0.25,
    scatter_size: float = 5.0,
    jitter: float = 0.08,
    violin_alpha: float = 0.70,
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.3,
    threshold_ls: str = "--",
    ylabel: str = "QC score",
    title: str = "",
    show_legend: bool = False,
    figsize: tuple[float, float] = (6.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Per-profile score spread: violin, box, or jitter-strip plot.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Same flexible input as :func:`plot_qc_scores`.
    station_labels : list of str or None
    profile_colors : dict, list, or None
    score_threshold : float, default ``0.5``
    kind : ``"violin"`` | ``"box"`` | ``"strip"``, default ``"violin"``
        Plot style.  ``"strip"`` shows only jittered scatter dots.
    show_scatter : bool, default ``True``
        Overlay jittered scatter dots (ignored when *kind* is ``"strip"``).
    scatter_alpha : float, default ``0.25``
    scatter_size : float, default ``5.0``
    jitter : float, default ``0.08``
        Maximum horizontal jitter width.
    violin_alpha : float, default ``0.70``
    threshold_color, threshold_lw, threshold_ls : str, float, str
    ylabel, title : str
    show_legend : bool, default ``False``
    figsize : (w, h), default ``(6.0, 4.0)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    pscores, _ = _normalise_qc_input(scores, station_labels)
    profiles   = list(pscores.keys())
    pcolors    = _resolve_profile_colors(profiles, profile_colors)

    vio_data:   list[np.ndarray] = []
    vio_labels: list[str]        = []
    vio_colors: list[str]        = []

    for prof in profiles:
        sc  = pscores[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        vio_data.append(fin)
        vio_labels.append(prof if prof != "_" else "all")
        vio_colors.append(pcolors[prof])

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    positions = np.arange(1, len(vio_data) + 1, dtype=float)
    rng       = np.random.default_rng(0)

    if kind == "violin" and vio_data:
        parts = ax.violinplot(vio_data, positions=positions,
                               showmedians=True, showextrema=True)
        for body, col in zip(parts["bodies"], vio_colors):
            body.set_facecolor(col)
            body.set_alpha(violin_alpha)
        for part_key in ("cmedians", "cmins", "cmaxes", "cbars"):
            if part_key in parts:
                parts[part_key].set_color("0.3")
                parts[part_key].set_linewidth(0.9)

    elif kind == "box" and vio_data:
        bp = ax.boxplot(vio_data, positions=positions, patch_artist=True,
                        widths=0.5, notch=False,
                        medianprops=dict(color="0.2", lw=1.5),
                        whiskerprops=dict(lw=0.8),
                        capprops=dict(lw=0.8),
                        flierprops=dict(marker=".", ms=3, alpha=0.4))
        for patch, col in zip(bp["boxes"], vio_colors):
            patch.set_facecolor(col)
            patch.set_alpha(violin_alpha)

    # jitter scatter overlay (always for "strip", optional for others)
    if kind == "strip" or show_scatter:
        for i, (sc, col) in enumerate(zip(vio_data, vio_colors)):
            jit = rng.uniform(-jitter, jitter, sc.size)
            ax.scatter(positions[i] + jit, sc,
                       color=col, alpha=scatter_alpha, s=scatter_size, zorder=3)

    ax.axhline(score_threshold, color=threshold_color,
               lw=threshold_lw, ls=threshold_ls, zorder=4)
    ax.axhspan(0, score_threshold, color="#fde8e8", alpha=0.20, zorder=0)

    ax.set_xticks(positions)
    ax.set_xticklabels(vio_labels, fontsize=8)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel(ylabel, fontsize=8)
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")

    return ax


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_summary
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_summary(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    score_threshold: float = 0.5,
    show_scatter: bool = True,
    n_bins: int = 20,
    show_kde: bool = True,
    kde_bw: float | str | None = None,
    spread_kind: str = "violin",
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    suptitle: str = "",
    figsize: tuple[float, float] = (13.0, 9.5),
) -> plt.Figure:
    """
    Three-panel QC summary figure.

    Layout::

        ┌────────────────────────────────┐
        │  (a) Per-station bar chart     │  full width
        ├──────────────────┬─────────────┤
        │  (b) Score dist  │  (c) Spread │
        │  histogram + KDE │  per profile│
        └──────────────────┴─────────────┘

    Panels (b) and (c) are produced by :func:`plot_qc_score_distribution`
    and :func:`plot_qc_score_spread` respectively, so they can also be
    used as standalone figures.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Same flexible input as :func:`plot_qc_scores`.
    station_labels : list of str or None
    profile_colors : dict, list, or None
    score_threshold : float, default ``0.5``
    show_scatter : bool, default ``True``
    n_bins : int, default ``20``
    show_kde : bool, default ``True``
    kde_bw : float, str, or None
    spread_kind : ``"violin"`` | ``"box"`` | ``"strip"``
    tick_every : int or ``"auto"``
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    station_tick_config : :class:`~pycsamt.ai.plot._style.StationTickConfig` or None
    suptitle : str
    figsize : (w, h), default ``(13.0, 9.5)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles  = list(pscores.keys())
    pcolors   = _resolve_profile_colors(profiles, profile_colors)
    tick_cfg  = _make_tick_config(tick_every, tick_label_rotation,
                                   tick_fontsize, station_tick_config)

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs  = fig.add_gridspec(
        2, 2,
        height_ratios=[1, 0.75],
        hspace=0.45,
        wspace=0.30,
    )
    ax_bar  = fig.add_subplot(gs[0, :])
    ax_hist = fig.add_subplot(gs[1, 0])
    ax_vio  = fig.add_subplot(gs[1, 1])

    # ── (a) bar chart ──────────────────────────────────────────────────────
    _render_bar_chart(
        ax_bar, pscores, plabels, pcolors,
        score_threshold=score_threshold,
        show_scatter=show_scatter,
        scatter_src=(scores if isinstance(scores, pd.DataFrame) else None),
        tick_cfg=tick_cfg,
        show_zone_labels=True,
    )
    ax_bar.set_title("(a) Per-station QC scores",
                     fontsize=9, fontweight="bold")

    # ── (b) score distribution ─────────────────────────────────────────────
    plot_qc_score_distribution.__wrapped__(
        scores,
        station_labels=station_labels,
        profile_colors=profile_colors,
        score_threshold=score_threshold,
        n_bins=n_bins,
        show_kde=show_kde,
        kde_bw=kde_bw,
        title="(b) Score distribution",
        ax=ax_hist,
    )

    # ── (c) score spread ───────────────────────────────────────────────────
    plot_qc_score_spread.__wrapped__(
        scores,
        station_labels=station_labels,
        profile_colors=profile_colors,
        score_threshold=score_threshold,
        kind=spread_kind,
        title=f"(c) Spread by profile  ({spread_kind})",
        ax=ax_vio,
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig
