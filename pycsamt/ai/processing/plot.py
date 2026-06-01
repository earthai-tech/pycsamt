# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Visualisation helpers for :class:`~pycsamt.ai.processing.qc.EMQCScorer` output.

Four public functions cover the main inspection use-cases:

* :func:`plot_qc_scores`        — per-station bar chart, multi-profile aware
* :func:`plot_qc_heatmap`       — station × frequency quality heat-map
* :func:`plot_qc_feature_heatmap` — one panel per QC feature
* :func:`plot_qc_summary`       — 3-panel summary figure
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from ..plot._style import EMStyle, EM_COLORS, EM_CMAPS

__all__ = [
    "plot_qc_scores",
    "plot_qc_heatmap",
    "plot_qc_feature_heatmap",
    "plot_qc_summary",
]

# ── default profile colour cycle ──────────────────────────────────────────── #
_PROFILE_COLORS = [
    "#1f77b4", "#ff7f0e", "#2ca02c",
    "#d62728", "#9467bd", "#8c564b",
    "#e377c2", "#7f7f7f", "#bcbd22",
]

_FEATURE_LABELS: Dict[str, str] = {
    "snr": "SNR",
    "swift_skew": "Swift skew |β|",
    "asym": r"Asymmetry  log$_{10}$(|Z$_{xy}$|/|Z$_{yx}$|)",
    "phase_xy": "Phase XY (°)",
    "phase_yx": "Phase YX (°)",
    "score": "QC score",
}

# ─────────────────────────────────────────────────────────────────────────────
# Internal helpers
# ─────────────────────────────────────────────────────────────────────────────

def _resolve_profile_colors(
    profile_names: List[str],
    user_colors: Optional[Any],
) -> Dict[str, str]:
    if isinstance(user_colors, dict):
        return {p: user_colors.get(p, _PROFILE_COLORS[i % len(_PROFILE_COLORS)])
                for i, p in enumerate(profile_names)}
    if isinstance(user_colors, (list, tuple)):
        return {p: user_colors[i % len(user_colors)]
                for i, p in enumerate(profile_names)}
    # None → default cycle
    return {p: _PROFILE_COLORS[i % len(_PROFILE_COLORS)]
            for i, p in enumerate(profile_names)}


def _normalise_qc_input(
    scores: Any,
    station_labels: Optional[List[str]],
) -> Tuple[Dict[str, np.ndarray], Dict[str, Optional[List[str]]]]:
    """
    Convert flexible input to ``{profile → 1-D score array}`` plus labels.

    Accepted formats
    ----------------
    * ``ndarray`` (n_st,)          single profile
    * ``dict[str, ndarray]``       one entry per profile
    * ``pd.DataFrame``             with ``score`` column; optional
                                   ``profile`` and ``station`` columns
    """
    if isinstance(scores, np.ndarray):
        arr = scores.ravel().astype(float)
        return {"_": arr}, {"_": station_labels}

    if isinstance(scores, dict):
        pscores = {k: np.asarray(v, dtype=float).ravel()
                   for k, v in scores.items()}
        plabels = {k: None for k in pscores}
        return pscores, plabels

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
        f"'scores' must be ndarray, dict, or DataFrame; got {type(scores).__name__}."
    )


def _logT_grid_from_df(
    df: pd.DataFrame,
    station_order: List[str],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pivot a per-(station, freq) DataFrame into a 2-D score matrix.

    Returns
    -------
    score2d : ndarray, shape (n_f, n_st)
    freqs   : ndarray, shape (n_f,) Hz
    x_cents : ndarray, shape (n_st,) station x-positions
    """
    freqs = np.sort(df["freq"].unique())[::-1]   # descending Hz → ascending T
    n_f   = freqs.size
    n_st  = len(station_order)
    mat   = np.full((n_f, n_st), np.nan)
    st_idx = {s: i for i, s in enumerate(station_order)}
    fr_idx = {f: i for i, f in enumerate(freqs)}
    for _, row in df.iterrows():
        si = st_idx.get(row["station"])
        fi = fr_idx.get(row["freq"])
        if si is not None and fi is not None:
            mat[fi, si] = float(row["score"])
    return mat, freqs, np.arange(n_st, dtype=float)


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_scores
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_scores(
    scores: Any,
    *,
    station_labels: Optional[List[str]] = None,
    profile_colors: Optional[Any] = None,
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
    bar_edgecolor: str = "none",
    # scatter overlay (per-frequency scores, when DataFrame input)
    show_scatter: bool = True,
    scatter_alpha: float = 0.22,
    scatter_size: float = 5.0,
    scatter_jitter: float = 0.18,
    # profile annotations
    show_profile_labels: bool = True,
    profile_label_y: float = 1.04,
    profile_label_fontsize: int = 8,
    # separator between profiles
    separator_color: str = "gray",
    separator_lw: float = 0.7,
    separator_ls: str = "--",
    # zone text on right margin
    show_zone_labels: bool = True,
    zone_label_fontsize: int = 8,
    # axes styling
    xlabel: str = "Station",
    ylabel: str = "QC score",
    title: str = "",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    ylim: Optional[Tuple[float, float]] = None,
    show_grid: bool = True,
    figsize: Tuple[float, float] = (10.0, 4.2),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Per-station QC-score bar chart, optionally grouped by profile.

    Supports single-profile 1-D arrays, multi-profile ``dict`` inputs,
    and :class:`~pycsamt.ai.processing.qc.EMQCScorer`\\ ``.score_table()``
    DataFrames.  A configurable threshold line separates *good* from
    *review* stations, and an optional scatter overlay shows the spread
    of individual per-frequency scores behind each bar.

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
        *scores* is a DataFrame.  Defaults to integer indices otherwise.
    profile_colors : dict, list, or None
        Mapping ``{profile_name: color_string}``, ordered list of colours,
        or ``None`` to use the built-in cycle.
    score_threshold : float, default ``0.5``
        Threshold line and rejection zone boundary.
    threshold_color : str, default ``"#c0392b"``
    threshold_lw : float, default ``1.4``
    threshold_ls : str, default ``"--"``
    bad_zone_color, bad_zone_alpha : str, float
        Fill colour and opacity of the rejection shading.
    bar_width : float, default ``0.80``
    bar_alpha : float, default ``0.88``
    bar_edgecolor : str, default ``"none"``
    show_scatter : bool, default ``True``
        When *scores* is a DataFrame with a ``freq`` column, draw
        semi-transparent per-frequency score dots behind each bar.
    scatter_alpha : float, default ``0.22``
    scatter_size : float, default ``5.0``
    scatter_jitter : float, default ``0.18``
        Maximum random horizontal jitter applied to scatter dots.
    show_profile_labels : bool, default ``True``
        Annotate each profile group with its name above the bars.
    profile_label_y : float, default ``1.04``
        Y-axis coordinate (in axes fraction) of the profile labels.
    profile_label_fontsize : int, default ``8``
    separator_color, separator_lw, separator_ls : str, float, str
        Style of the vertical profile-separator line.
    show_zone_labels : bool, default ``True``
        Add "Good ▶" and "◀ Review" text annotations on the right margin.
    zone_label_fontsize : int, default ``8``
    xlabel, ylabel, title : str
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    ylim : (ymin, ymax) or None
        Override y-axis limits.
    show_grid : bool, default ``True``
    figsize : (w, h), default ``(10.0, 4.2)``
        Ignored when *ax* is supplied.
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`

    Examples
    --------
    >>> from pycsamt.ai.processing import EMQCScorer
    >>> from pycsamt.ai.processing.plot import plot_qc_scores
    >>> scorer = EMQCScorer(use_ml=False)
    >>> scorer.fit(np.zeros((1, 5)))
    >>> sc_by_profile = {"L18": scorer.transform(feat18), ...}
    >>> plot_qc_scores(sc_by_profile, score_threshold=0.5)
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles         = list(pscores.keys())
    pcolors          = _resolve_profile_colors(profiles, profile_colors)

    # ── build scatter source (per-freq scores when available) ─────────────
    scatter_src: Dict[str, Optional[Tuple[np.ndarray, np.ndarray]]] = {}
    if show_scatter and isinstance(scores, pd.DataFrame) and "freq" in scores.columns:
        has_prof = "profile" in scores.columns
        has_st   = "station" in scores.columns
        for prof in profiles:
            sub = scores[scores["profile"] == prof] if has_prof else scores
            if not has_st:
                scatter_src[prof] = None
                continue
            st_order  = plabels[prof] or list(sub["station"].unique())
            x_scatter = []
            y_scatter = []
            for si, st in enumerate(st_order):
                rows = sub[sub["station"] == st]["score"].values
                if rows.size:
                    rng = np.random.default_rng(abs(hash(st)) % 2**31)
                    jit = rng.uniform(-scatter_jitter, scatter_jitter, rows.size)
                    x_scatter.append(np.full(rows.size, float(si)) + jit)
                    y_scatter.append(rows)
            if x_scatter:
                scatter_src[prof] = (
                    np.concatenate(x_scatter),
                    np.concatenate(y_scatter),
                )
            else:
                scatter_src[prof] = None
    else:
        scatter_src = {p: None for p in profiles}

    # ── create / reuse axes ────────────────────────────────────────────────
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    # ── rejection zone (behind everything) ────────────────────────────────
    ax.axhspan(0, score_threshold, color=bad_zone_color,
               alpha=bad_zone_alpha, zorder=0)

    # ── draw bars profile by profile ───────────────────────────────────────
    x_offset = 0
    all_tick_positions: List[float] = []
    all_tick_labels:    List[str]   = []
    sep_positions: List[float]      = []

    for i, prof in enumerate(profiles):
        sc    = pscores[prof]
        n     = len(sc)
        col   = pcolors[prof]
        labs  = plabels[prof]
        x_pos = np.arange(x_offset, x_offset + n, dtype=float)

        # scatter first (behind bars)
        if scatter_src[prof] is not None:
            xs, ys = scatter_src[prof]
            # adjust xs relative to offset
            xs_abs = xs + x_offset
            ax.scatter(xs_abs, ys, color=col, alpha=scatter_alpha,
                       s=scatter_size, zorder=1)

        # bars
        sc_plot = np.where(np.isfinite(sc), sc, 0.0)
        ax.bar(
            x_pos, sc_plot,
            color=col, width=bar_width,
            edgecolor=bar_edgecolor, alpha=bar_alpha, zorder=2,
        )

        # per-bar thin indicator line at score value
        for xi, yi in zip(x_pos, sc_plot):
            ax.hlines(yi, xi - bar_width * 0.45, xi + bar_width * 0.45,
                      colors=col, lw=0.8, zorder=3, alpha=0.7)

        # profile label above bars
        if show_profile_labels and prof != "_":
            mid_x = x_offset + (n - 1) / 2.0
            ax.text(
                mid_x, profile_label_y,
                prof,
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=profile_label_fontsize,
                fontweight="bold",
                color=col,
            )

        # separator (after every profile except the last)
        if i < len(profiles) - 1:
            sep_x = x_offset + n - 0.5
            ax.axvline(sep_x, color=separator_color,
                       lw=separator_lw, ls=separator_ls,
                       alpha=0.6, zorder=1)
            sep_positions.append(sep_x)

        # accumulate tick info
        all_tick_positions.extend(x_pos.tolist())
        if labs is not None:
            all_tick_labels.extend(labs)
        else:
            all_tick_labels.extend([str(j) for j in range(x_offset, x_offset + n)])

        x_offset += n

    # ── threshold line ─────────────────────────────────────────────────────
    ax.axhline(score_threshold, color=threshold_color,
               lw=threshold_lw, ls=threshold_ls, zorder=4)

    # ── zone labels on right margin ────────────────────────────────────────
    if show_zone_labels:
        y_good   = (1.0 + score_threshold) / 2.0
        y_review = score_threshold / 2.0
        ax.text(1.002, y_good,   "Good ▶",   transform=ax.transAxes,
                ha="left", va="center", fontsize=zone_label_fontsize,
                color="#2ca02c", fontweight="bold")
        ax.text(1.002, y_review, "◀ Review", transform=ax.transAxes,
                ha="left", va="center", fontsize=zone_label_fontsize,
                color=threshold_color, fontweight="bold")

    # ── x-axis ─────────────────────────────────────────────────────────────
    ax.set_xticks(all_tick_positions)
    ha = "right" if tick_label_rotation > 20 else "center"
    ax.set_xticklabels(
        all_tick_labels,
        rotation=tick_label_rotation, ha=ha,
        fontsize=tick_fontsize,
    )
    ax.set_xlim(-0.8, x_offset - 0.2)

    # ── y-axis ─────────────────────────────────────────────────────────────
    if ylim is None:
        ax.set_ylim(0, 1.15)
    else:
        ax.set_ylim(*ylim)

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")

    if show_grid:
        ax.grid(True, axis="y", ls=":", lw=0.4, color="gray", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)

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
    station_labels: Optional[List[str]] = None,
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
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    figsize: Tuple[float, float] = (10.0, 5.2),
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    """
    Station × frequency QC-score heat-map.

    Renders the full per-(station, frequency) score matrix so that
    individual frequency bands with quality problems are immediately
    visible.  Green = good, red = bad; a white contour can mark the
    configurable threshold boundary.

    Parameters
    ----------
    scores : DataFrame
        Must contain ``station``, ``freq``, and ``score`` columns
        (output of :meth:`~pycsamt.ai.processing.qc.EMQCScorer.score_table`).
        Alternatively, a 2-D ``ndarray`` shape ``(n_st, n_f)`` may be
        passed together with *station_labels* and a ``freqs=`` keyword
        argument (see notes below).
    station_labels : list of str or None
        X-axis tick labels.  Inferred from the ``station`` column when
        *scores* is a DataFrame.
    score_threshold : float, default ``0.5``
        Threshold for the optional contour overlay.
    cmap : str, default ``"RdYlGn"``
        Diverging quality colormap.
    vmin, vmax : float
        Colour scale range.
    period_up : bool, default ``True``
        Long period at the top (MT convention).
    show_threshold_contour : bool, default ``True``
    contour_color, contour_lw : str, float
    n_yticks : int, default ``7``
    colorbar_label, xlabel, ylabel, title : str
    tick_label_rotation : float
    tick_fontsize : int
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
        raise ValueError(
            "DataFrame must have columns 'station', 'freq', 'score'."
        )

    # ── build ordered station list ─────────────────────────────────────────
    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = list(scores["station"].unique())
        try:
            st_order = sorted(st_order)
        except TypeError:
            pass

    mat, freqs, x_cents = _logT_grid_from_df(scores, st_order)
    n_f, n_st = mat.shape

    # ── log-period y grid ──────────────────────────────────────────────────
    log_T  = np.log10(1.0 / freqs)
    d_lT   = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0]    = log_T[0]   - 0.5 * abs(d_lT[0]) if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else []
    y_edge[-1]   = log_T[-1]  + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    x_edge = np.concatenate([
        [x_cents[0] - 0.5],
        0.5 * (x_cents[:-1] + x_cents[1:]),
        [x_cents[-1] + 0.5],
    ]) if n_st > 1 else np.array([-0.5, 0.5])

    X, Y = np.meshgrid(x_edge, y_edge)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    qm = ax.pcolormesh(
        X, Y, mat,
        cmap=cmap, vmin=vmin, vmax=vmax,
        shading="flat", rasterized=True,
    )
    if period_up:
        ax.invert_yaxis()

    # ── threshold contour ──────────────────────────────────────────────────
    if show_threshold_contour:
        try:
            x_c = 0.5 * (x_edge[:-1] + x_edge[1:])
            y_c = 0.5 * (y_edge[:-1] + y_edge[1:])
            Xc, Yc = np.meshgrid(x_c, y_c)
            ax.contour(
                Xc, Yc, mat,
                levels=[score_threshold],
                colors=[contour_color], linewidths=[contour_lw],
            )
        except Exception:
            pass

    # ── colorbar ───────────────────────────────────────────────────────────
    cb = plt.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(colorbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # ── y-axis ticks (log-period) ──────────────────────────────────────────
    pos  = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$")
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    # ── x-axis ticks ──────────────────────────────────────────────────────
    ax.set_xticks(x_cents)
    ha = "right" if tick_label_rotation > 20 else "center"
    ax.set_xticklabels(
        st_order, rotation=tick_label_rotation, ha=ha, fontsize=tick_fontsize,
    )
    ax.set_xlim(x_edge[0], x_edge[-1])
    ax.set_xlabel(xlabel, fontsize=8)

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
    features: Optional[List[str]] = None,
    station_labels: Optional[List[str]] = None,
    cmaps: Optional[Dict[str, str]] = None,
    period_up: bool = True,
    n_yticks: int = 5,
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 6,
    clim_pct: Tuple[float, float] = (2.0, 98.0),
    title: str = "",
    figsize: Optional[Tuple[float, float]] = None,
) -> plt.Figure:
    """
    One panel per QC feature — station × frequency heat-maps.

    Useful for diagnosing which signal-quality metric drives poor scores
    at a given station or frequency band.

    Parameters
    ----------
    df : DataFrame
        Output of :meth:`~pycsamt.ai.processing.qc.EMQCScorer.score_table`.
        Must contain ``station``, ``freq``, and feature columns.
    features : list of str or None
        Feature columns to plot.  Defaults to
        ``["snr", "swift_skew", "asym", "phase_xy", "phase_yx"]``.
    station_labels : list of str or None
    cmaps : dict or None
        ``{feature_name: colormap_string}`` overrides.
        Defaults: SNR → ``"YlGn"``, Swift skew → ``"YlOrRd_r"``,
        others → ``"RdBu_r"``.
    period_up : bool, default ``True``
    n_yticks : int, default ``5``
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``6``
    clim_pct : (lo, hi), default ``(2.0, 98.0)``
        Percentile colour limits for each feature independently.
    title : str
        Figure-level title.
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

    # ── station order ──────────────────────────────────────────────────────
    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = list(df["station"].unique())
        try:
            st_order = sorted(st_order)
        except TypeError:
            pass

    n_feat  = len(features)
    n_st    = len(st_order)
    freqs   = np.sort(df["freq"].unique())[::-1]
    n_f     = freqs.size

    if figsize is None:
        figsize = (max(8.0, n_st * 0.35), 2.8 * n_feat)

    fig, axes = plt.subplots(n_feat, 1, figsize=figsize, sharex=True,
                             gridspec_kw={"hspace": 0.45})
    if n_feat == 1:
        axes = [axes]

    # ── grid edges (shared across panels) ─────────────────────────────────
    log_T  = np.log10(1.0 / freqs)
    if n_f > 1:
        d    = np.diff(log_T)
        ye   = np.empty(n_f + 1)
        ye[0] = log_T[0] - 0.5 * abs(d[0])
        ye[1:-1] = log_T[:-1] + 0.5 * d
        ye[-1]   = log_T[-1] + 0.5 * abs(d[-1])
    else:
        ye = np.array([log_T[0] - 0.5, log_T[0] + 0.5])

    x_cents = np.arange(n_st, dtype=float)
    xe = np.concatenate([
        [x_cents[0] - 0.5],
        0.5 * (x_cents[:-1] + x_cents[1:]),
        [x_cents[-1] + 0.5],
    ]) if n_st > 1 else np.array([-0.5, 0.5])
    X, Y = np.meshgrid(xe, ye)

    # ── per-feature panels ─────────────────────────────────────────────────
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

        # make colour limits symmetric for diverging quantities
        if feat in ("asym", "phase_xy", "phase_yx"):
            v = max(abs(vmin), abs(vmax))
            vmin, vmax = -v, v

        qm = ax.pcolormesh(
            X, Y, mat,
            cmap=cmap_map.get(feat, "RdBu_r"),
            vmin=vmin, vmax=vmax,
            shading="flat", rasterized=True,
        )
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

    # ── bottom x-axis ──────────────────────────────────────────────────────
    axes[-1].set_xticks(x_cents)
    ha = "right" if tick_label_rotation > 20 else "center"
    axes[-1].set_xticklabels(
        st_order, rotation=tick_label_rotation, ha=ha, fontsize=tick_fontsize,
    )
    axes[-1].set_xlim(xe[0], xe[-1])
    axes[-1].set_xlabel("Station", fontsize=7)

    if title:
        fig.suptitle(title, fontsize=10, fontweight="bold", y=1.01)

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# plot_qc_summary
# ─────────────────────────────────────────────────────────────────────────────

@EMStyle()
def plot_qc_summary(
    scores: Any,
    *,
    station_labels: Optional[List[str]] = None,
    profile_colors: Optional[Any] = None,
    score_threshold: float = 0.5,
    show_scatter: bool = True,
    n_hist_bins: int = 20,
    kde_bw: Optional[float] = None,
    show_kde: bool = True,
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    suptitle: str = "",
    figsize: Tuple[float, float] = (13.0, 9.5),
) -> plt.Figure:
    """
    Three-panel QC summary figure.

    Layout::

        ┌────────────────────────────────┐
        │  (a) Per-station bar chart     │
        ├──────────────────┬─────────────┤
        │  (b) Score hist  │  (c) Violin │
        │  + KDE per prof  │  per profile│
        └──────────────────┴─────────────┘

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Same flexible input as :func:`plot_qc_scores`.
    station_labels : list of str or None
    profile_colors : dict, list, or None
    score_threshold : float, default ``0.5``
    show_scatter : bool, default ``True``
        Scatter overlay on the bar chart (requires freq column in DataFrame).
    n_hist_bins : int, default ``20``
    kde_bw : float or None
        Bandwidth for KDE curves.  ``None`` → Scott's rule.
    show_kde : bool, default ``True``
        Overlay a KDE on the histogram.
    tick_label_rotation : float, default ``45.0``
    tick_fontsize : int, default ``7``
    suptitle : str
    figsize : (w, h), default ``(13.0, 9.5)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles  = list(pscores.keys())
    pcolors   = _resolve_profile_colors(profiles, profile_colors)

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs  = fig.add_gridspec(
        2, 2,
        height_ratios=[1, 0.75],
        hspace=0.45,
        wspace=0.30,
    )
    ax_bar  = fig.add_subplot(gs[0, :])   # full-width bar chart
    ax_hist = fig.add_subplot(gs[1, 0])   # histogram + KDE
    ax_vio  = fig.add_subplot(gs[1, 1])   # violin per profile

    # ── (a) bar chart ──────────────────────────────────────────────────────
    # call the bare renderer with pre-normalised input so styles match
    _render_bar_chart(
        ax_bar, pscores, plabels, pcolors,
        score_threshold=score_threshold,
        show_scatter=show_scatter,
        scatter_src=(scores if isinstance(scores, pd.DataFrame) else None),
        tick_label_rotation=tick_label_rotation,
        tick_fontsize=tick_fontsize,
    )
    ax_bar.set_title("(a) Per-station QC scores", fontsize=9, fontweight="bold")

    # ── (b) score histogram / KDE per profile ─────────────────────────────
    all_scores_flat: List[np.ndarray] = []
    for prof in profiles:
        sc  = pscores[prof]
        col = pcolors[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        all_scores_flat.append(fin)
        label = prof if prof != "_" else "all"
        ax_hist.hist(
            fin, bins=n_hist_bins, range=(0, 1),
            color=col, alpha=0.35,
            edgecolor="none", density=True, label=label,
        )
        if show_kde and fin.size > 3:
            try:
                from scipy.stats import gaussian_kde
                bw = kde_bw if kde_bw is not None else "scott"
                kde = gaussian_kde(fin, bw_method=bw)
                xs  = np.linspace(0, 1, 300)
                ax_hist.plot(xs, kde(xs), color=col, lw=1.8)
            except ImportError:
                pass

    ax_hist.axvline(score_threshold, color="#c0392b", lw=1.3, ls="--")
    ax_hist.set_xlabel("QC score", fontsize=8)
    ax_hist.set_ylabel("Density", fontsize=8)
    ax_hist.set_xlim(0, 1)
    ax_hist.set_title("(b) Score distribution", fontsize=9, fontweight="bold")
    if len(profiles) > 1:
        ax_hist.legend(fontsize=7, framealpha=0.8)

    # ── (c) violin / box per profile ──────────────────────────────────────
    vio_data  = []
    vio_labels = []
    vio_colors = []
    for prof in profiles:
        sc  = pscores[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        vio_data.append(fin)
        vio_labels.append(prof if prof != "_" else "all")
        vio_colors.append(pcolors[prof])

    if vio_data:
        parts = ax_vio.violinplot(
            vio_data,
            positions=np.arange(1, len(vio_data) + 1),
            showmedians=True, showextrema=True,
        )
        for body, col in zip(parts["bodies"], vio_colors):
            body.set_facecolor(col)
            body.set_alpha(0.7)
        for part_key in ("cmedians", "cmins", "cmaxes", "cbars"):
            if part_key in parts:
                parts[part_key].set_color("0.3")
                parts[part_key].set_linewidth(0.9)

        # scatter jitter overlay
        rng = np.random.default_rng(0)
        for i, (sc, col) in enumerate(zip(vio_data, vio_colors), start=1):
            jit = rng.uniform(-0.08, 0.08, sc.size)
            ax_vio.scatter(
                np.full(sc.size, i) + jit, sc,
                color=col, alpha=0.25, s=5, zorder=2,
            )

        ax_vio.axhline(score_threshold, color="#c0392b", lw=1.3, ls="--", zorder=3)
        ax_vio.set_xticks(np.arange(1, len(vio_data) + 1))
        ax_vio.set_xticklabels(vio_labels, fontsize=8)
        ax_vio.set_ylim(0, 1.05)
        ax_vio.set_ylabel("QC score", fontsize=8)

    ax_vio.set_title("(c) Score spread by profile", fontsize=9, fontweight="bold")

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold", y=1.01)

    return fig


# ── private bar-chart renderer (shared between plot_qc_scores and summary) ── #

def _render_bar_chart(
    ax: plt.Axes,
    pscores: Dict[str, np.ndarray],
    plabels: Dict[str, Optional[List[str]]],
    pcolors: Dict[str, str],
    *,
    score_threshold: float,
    show_scatter: bool,
    scatter_src: Optional[pd.DataFrame],
    tick_label_rotation: float,
    tick_fontsize: int,
) -> None:
    ax.axhspan(0, score_threshold, color="#fde8e8", alpha=0.30, zorder=0)
    profiles = list(pscores.keys())
    x_offset = 0
    all_tick_pos: List[float] = []
    all_tick_lab: List[str]   = []

    for i, prof in enumerate(profiles):
        sc    = pscores[prof]
        n     = len(sc)
        col   = pcolors[prof]
        labs  = plabels[prof]
        x_pos = np.arange(x_offset, x_offset + n, dtype=float)

        # scatter overlay from df
        if show_scatter and scatter_src is not None and "freq" in scatter_src.columns:
            has_st = "station" in scatter_src.columns
            sub = (scatter_src[scatter_src["profile"] == prof]
                   if "profile" in scatter_src.columns else scatter_src)
            if has_st and labs:
                rng = np.random.default_rng(i)
                for si, st in enumerate(labs):
                    rows = sub[sub["station"] == st]["score"].values
                    if rows.size:
                        jit = rng.uniform(-0.18, 0.18, rows.size)
                        ax.scatter(
                            np.full(rows.size, x_offset + si) + jit,
                            rows, color=col, alpha=0.22, s=5, zorder=1,
                        )

        sc_plot = np.where(np.isfinite(sc), sc, 0.0)
        ax.bar(x_pos, sc_plot, color=col, width=0.80,
               edgecolor="none", alpha=0.88, zorder=2)

        if i < len(profiles) - 1:
            ax.axvline(x_offset + n - 0.5, color="gray",
                       lw=0.7, ls="--", alpha=0.6, zorder=1)

        if prof != "_":
            ax.text(
                x_offset + (n - 1) / 2.0, 1.04, prof,
                transform=ax.get_xaxis_transform(),
                ha="center", va="bottom",
                fontsize=8, fontweight="bold", color=col,
            )

        all_tick_pos.extend(x_pos.tolist())
        if labs:
            all_tick_lab.extend(labs)
        else:
            all_tick_lab.extend([str(j) for j in range(x_offset, x_offset + n)])

        x_offset += n

    ax.axhline(score_threshold, color="#c0392b", lw=1.4, ls="--", zorder=4)
    ax.set_xticks(all_tick_pos)
    ha = "right" if tick_label_rotation > 20 else "center"
    ax.set_xticklabels(all_tick_lab, rotation=tick_label_rotation,
                       ha=ha, fontsize=tick_fontsize)
    ax.set_xlim(-0.8, x_offset - 0.2)
    ax.set_ylim(0, 1.15)
    ax.set_xlabel("Station", fontsize=9)
    ax.set_ylabel("QC score", fontsize=9)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="gray", alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
