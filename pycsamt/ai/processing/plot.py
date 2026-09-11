# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Visualisation helpers for :mod:`pycsamt.ai.processing` output.

QC (:class:`~pycsamt.ai.processing.qc.EMQCScorer`)
    * :func:`plot_qc_scores`              — per-station bar chart
    * :func:`plot_qc_heatmap`             — station × freq. heat-map
    * :func:`plot_qc_feature_heatmap`     — one panel per QC feature
    * :func:`plot_qc_score_distribution`  — score histogram + KDE
    * :func:`plot_qc_score_spread`        — violin / box / strip
    * :func:`plot_qc_summary`             — 3-panel summary figure

Denoising (:class:`~pycsamt.ai.processing.denoise.EMDenoiser`)
    * :func:`plot_denoise_spectra`         — before/after spectra, per site
    * :func:`plot_denoise_noise_reduction` — per-station roughness reduction
    * :func:`plot_denoise_summary`         — combined summary figure

Anomaly detection (:class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`)
    * :func:`plot_anomaly_scores`              — per-station bar chart
    * :func:`plot_anomaly_score_distribution`  — score histogram + KDE
    * :func:`plot_anomaly_summary`             — 2-panel summary figure

Dimensionality (:class:`~pycsamt.ai.processing.classify.\
DimensionalityClassifier`)
    * :func:`plot_dimensionality_map`      — station × freq. class map
    * :func:`plot_predicted_strike_rose`             — 2-D strike rose diagram
    * :func:`plot_dimensionality_summary`  — 3-panel summary figure

Gap filling (:class:`~pycsamt.ai.processing.imputer.EMImputer`)
    * :func:`plot_imputer_gaps`            — station × freq. missing-data map
    * :func:`plot_imputer_validation`      — held-out true-vs-reconstructed
      scatter, faceted by component
    * :func:`plot_imputer_reconstruction`  — per-site spectra, reconstructed
      cells highlighted
    * :func:`plot_imputer_summary`         — 3-panel summary figure

Uncertainty calibration (:class:`~pycsamt.ai.processing.uncertainty.\
UncertaintyCalibrator`)
    * :func:`plot_uncertainty_map`         — station × freq. fractional-error
      heat-map (log colour scale)
    * :func:`plot_uncertainty_validation`  — held-out true-vs-calibrated
      scatter
    * :func:`plot_uncertainty_summary`     — 4-panel summary figure

Distortion triage (:class:`~pycsamt.ai.processing.distortion.\
DistortionTypeClassifier`)
    * :func:`plot_distortion_map`            — per-station regime bar chart
    * :func:`plot_distortion_feature_space`  — feature scatter coloured by
      regime, with rule thresholds
    * :func:`plot_distortion_summary`        — 4-panel summary figure

Time-series denoising (:class:`~pycsamt.ai.processing.tsdenoise.\
TimeSeriesDenoiser`)
    * :func:`plot_ts_denoise_mmf_split`  — raw + low-frequency envelope,
      high-frequency residual
    * :func:`plot_ts_denoise_segments`   — residual with SVM-flagged
      noisy windows shaded
    * :func:`plot_ts_denoise_summary`    — 3-panel before/after summary

Shared
    * :func:`plot_training_history` — train/val loss curve for any of the
      three network-based estimators above.
"""

from __future__ import annotations

import copy
from typing import Any

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import BoundaryNorm, ListedColormap, LogNorm

from ...api.labels import STATION_LABEL
from ...api.station import PYCSAMT_STATION_RENDERING
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
    "plot_denoise_spectra",
    "plot_denoise_noise_reduction",
    "plot_denoise_summary",
    "plot_anomaly_scores",
    "plot_anomaly_score_distribution",
    "plot_anomaly_summary",
    "plot_dimensionality_map",
    "plot_predicted_strike_rose",
    "plot_dimensionality_summary",
    "plot_imputer_gaps",
    "plot_imputer_validation",
    "plot_imputer_reconstruction",
    "plot_imputer_summary",
    "plot_uncertainty_map",
    "plot_uncertainty_validation",
    "plot_uncertainty_summary",
    "plot_distortion_map",
    "plot_distortion_feature_space",
    "plot_distortion_summary",
    "plot_ts_denoise_mmf_split",
    "plot_ts_denoise_segments",
    "plot_ts_denoise_summary",
    "plot_training_history",
]

# ── default profile colour cycle ──────────────────────────────────────────── #
_PROFILE_COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
]

_FEATURE_LABELS: dict[str, str] = {
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
    profile_names: list[str],
    user_colors: Any | None,
) -> dict[str, str]:
    if isinstance(user_colors, dict):
        return {
            p: user_colors.get(p, _PROFILE_COLORS[i % len(_PROFILE_COLORS)])
            for i, p in enumerate(profile_names)
        }
    if isinstance(user_colors, (list, tuple)):
        return {
            p: user_colors[i % len(user_colors)]
            for i, p in enumerate(profile_names)
        }
    return {
        p: _PROFILE_COLORS[i % len(_PROFILE_COLORS)]
        for i, p in enumerate(profile_names)
    }


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
        pscores = {
            k: np.asarray(v, dtype=float).ravel() for k, v in scores.items()
        }
        # station_labels applies per-profile only when its length
        # matches that profile's score count -- e.g. the same station
        # set scored two ways (before/after). A shorter or longer list
        # (different surveys, mismatched profile lengths) falls back
        # to the previous behaviour (plain integer tick labels) rather
        # than silently mis-aligning names to the wrong bars.
        plabels = {
            k: (
                list(station_labels)
                if station_labels is not None
                and len(station_labels) == len(sc)
                else None
            )
            for k, sc in pscores.items()
        }
        return pscores, plabels

    if isinstance(scores, pd.DataFrame):
        if "score" not in scores.columns:
            raise ValueError(
                "DataFrame must contain a 'score' column "
                "(use EMQCScorer.score_table())."
            )
        has_prof = "profile" in scores.columns
        has_st = "station" in scores.columns

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


def _apply_section_station_axis(
    ax: plt.Axes,
    x_cents: np.ndarray,
    st_order: list,
    xlim: tuple[float, float],
) -> None:
    """
    Draw the station axis at the top of a period/station 2-D section
    with pyCSAMT's shared station-rendering convention -- downward
    triangle markers just above the axes, station names above them --
    matching every other section in the package (see
    :data:`~pycsamt.api.station.PYCSAMT_STATION_RENDERING` and, for
    example, :mod:`pycsamt.emtools.resphase_psection`).
    """
    style = copy.deepcopy(
        PYCSAMT_STATION_RENDERING.style_for("pseudosection")
    )
    style.xlabel = STATION_LABEL
    style.apply(ax, x_cents, [str(s) for s in st_order], xlim=xlim)


def _station_freq_grid(
    df: pd.DataFrame,
    station_order: list[str],
    value_col: str = "score",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Pivot a long ``(station, freq, value_col)`` table to a dense
    ``(n_freq, n_station)`` grid, ordered by *station_order* and by
    descending frequency.

    Shared by :func:`plot_qc_heatmap` (``value_col="score"``) and
    :func:`plot_dimensionality_map` (``value_col="dim"``).
    """
    freqs = np.sort(df["freq"].unique())[::-1]
    n_f, n_st = freqs.size, len(station_order)
    mat = np.full((n_f, n_st), np.nan)
    st_idx = {s: i for i, s in enumerate(station_order)}
    fr_idx = {f: i for i, f in enumerate(freqs)}
    for _, row in df.iterrows():
        si = st_idx.get(row["station"])
        fi = fr_idx.get(row["freq"])
        if si is not None and fi is not None:
            val = row.get(value_col, np.nan)
            if pd.notna(val):
                mat[fi, si] = float(val)
    return mat, freqs, np.arange(n_st, dtype=float)


# Backward-compatible alias used by plot_qc_heatmap.
def _logT_grid_from_df(
    df: pd.DataFrame,
    station_order: list[str],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return _station_freq_grid(df, station_order, value_col="score")


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
    # ── generalised for reuse beyond QC scores (anomaly, etc.) ────────────
    bad_zone_side: str = "below",
    ylim: tuple[float, float] = (0.0, 1.15),
    ylabel: str = "QC score",
    zone_label_good: str = "Good ▶",
    zone_label_bad: str = "◀ Review",
    good_color: str = "#2ca02c",
) -> None:
    """
    Render bar chart onto *ax* (no figure creation, no return value).

    *bad_zone_side* selects whether the shaded "reject" band sits below
    the threshold (default -- QC scores, where high = good) or above it
    (anomaly scores, where high = bad).
    """
    if bad_zone_side == "below":
        ax.axhspan(
            ylim[0],
            score_threshold,
            color=bad_zone_color,
            alpha=bad_zone_alpha,
            zorder=0,
        )
    else:
        ax.axhspan(
            score_threshold,
            ylim[1],
            color=bad_zone_color,
            alpha=bad_zone_alpha,
            zorder=0,
        )

    profiles = list(pscores.keys())
    x_offset = 0
    all_positions: list[float] = []
    all_labels: list[str] = []

    for i, prof in enumerate(profiles):
        sc = pscores[prof]
        n = len(sc)
        col = pcolors[prof]
        labs = plabels[prof]
        xpos = np.arange(x_offset, x_offset + n, dtype=float)

        # ── optional per-frequency scatter ────────────────────────────────
        if (
            show_scatter
            and scatter_src is not None
            and "freq" in scatter_src.columns
            and "station" in scatter_src.columns
        ):
            sub = (
                scatter_src[scatter_src["profile"] == prof]
                if "profile" in scatter_src.columns
                else scatter_src
            )
            if labs:
                rng = np.random.default_rng(i)
                for si, st in enumerate(labs):
                    rows = sub[sub["station"] == st]["score"].values
                    if rows.size:
                        jit = rng.uniform(-0.18, 0.18, rows.size)
                        ax.scatter(
                            np.full(rows.size, x_offset + si) + jit,
                            rows,
                            color=col,
                            alpha=0.22,
                            s=5,
                            zorder=1,
                        )

        # ── bars ──────────────────────────────────────────────────────────
        sc_plot = np.where(np.isfinite(sc), sc, 0.0)
        ax.bar(
            xpos,
            sc_plot,
            color=col,
            width=bar_width,
            edgecolor="none",
            alpha=bar_alpha,
            zorder=2,
        )

        # thin top-edge line at score value
        for xi, yi in zip(xpos, sc_plot):
            ax.hlines(
                yi,
                xi - bar_width * 0.44,
                xi + bar_width * 0.44,
                colors=col,
                lw=0.9,
                zorder=3,
                alpha=0.65,
            )

        # ── profile separator ─────────────────────────────────────────────
        if i < len(profiles) - 1:
            ax.axvline(
                x_offset + n - 0.5,
                color=separator_color,
                lw=separator_lw,
                ls=separator_ls,
                alpha=0.6,
                zorder=1,
            )

        # ── profile label above bars ──────────────────────────────────────
        if show_profile_labels and prof != "_":
            ax.text(
                x_offset + (n - 1) / 2.0,
                profile_label_y,
                prof,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=profile_label_fontsize,
                fontweight="bold",
                color=col,
            )

        # ── accumulate tick info ──────────────────────────────────────────
        all_positions.extend(xpos.tolist())
        if labs is not None:
            all_labels.extend(labs)
        else:
            all_labels.extend([str(j) for j in range(x_offset, x_offset + n)])

        x_offset += n

    # ── threshold line ─────────────────────────────────────────────────────
    # A dashed/dotted linestyle with lw<=0 makes matplotlib raise
    # ("at least one value in the dash list must be positive") instead
    # of simply drawing nothing -- skip the line entirely rather than
    # pass a degenerate width through.
    if threshold_lw > 0:
        ax.axhline(
            score_threshold,
            color=threshold_color,
            lw=threshold_lw,
            ls=threshold_ls,
            zorder=4,
        )

    # ── zone labels on right margin ────────────────────────────────────────
    if show_zone_labels:
        if bad_zone_side == "below":
            y_good = (1.0 + score_threshold) / 2.0
            y_bad = score_threshold / 2.0
        else:
            frac_thr = (score_threshold - ylim[0]) / (ylim[1] - ylim[0])
            y_good = frac_thr / 2.0
            y_bad = (1.0 + frac_thr) / 2.0
        kw = dict(
            transform=ax.transAxes,
            va="center",
            fontsize=zone_label_fontsize,
            fontweight="bold",
        )
        ax.text(
            1.002, y_good, zone_label_good, ha="left", color=good_color, **kw
        )
        ax.text(
            1.002,
            y_bad,
            zone_label_bad,
            ha="left",
            color=threshold_color,
            **kw,
        )

    # ── x-axis via StationTickConfig ───────────────────────────────────────
    pos_arr = np.asarray(all_positions)
    tick_cfg.apply(
        ax,
        pos_arr,
        all_labels,
        xlabel=xlabel,
        xlim=(-0.8, x_offset - 0.2),
    )

    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=9)
    if show_grid:
        ax.grid(
            True, axis="y", ls=":", lw=0.4, color="gray", alpha=0.5, zorder=0
        )
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
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
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
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)
    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _render_bar_chart(
        ax,
        pscores,
        plabels,
        pcolors,
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
            for p in profiles
            if p != "_"
        ]
        # Match the threshold line's own visibility guard above --
        # a Line2D legend handle with a dashed style and lw<=0 hits
        # the same matplotlib dash-pattern crash the drawn line does,
        # and there is nothing to label if the line itself is hidden.
        if threshold_lw > 0:
            handles.append(
                plt.Line2D(
                    [],
                    [],
                    color=threshold_color,
                    ls=threshold_ls,
                    lw=threshold_lw,
                    label=f"Review threshold ({score_threshold:.2f})",
                )
            )
        ax.legend(
            handles=handles, loc="lower right", fontsize=7.5, framealpha=0.9
        )

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
    station_markers: bool = True,
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
    station_markers : bool, default ``True``
        Draw the station axis at the top with pyCSAMT's shared
        downward-triangle convention (
        :data:`~pycsamt.api.station.PYCSAMT_STATION_RENDERING`), the
        same rendering used throughout the rest of the package for
        period/station sections. Set to ``False`` to fall back to a
        plain bottom axis controlled by *tick_every* /
        *tick_label_rotation* / *tick_fontsize* / *station_tick_config*
        instead.
    tick_every : int or ``"auto"``
        Only used when ``station_markers=False``.
    tick_label_rotation : float, default ``45.0``
        Only used when ``station_markers=False``.
    tick_fontsize : int, default ``7``
        Only used when ``station_markers=False``.
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
        Only used when ``station_markers=False``.
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

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(
            scores["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )

    mat, freqs, x_cents = _logT_grid_from_df(scores, st_order)
    n_f, n_st = mat.shape

    log_T = np.log10(1.0 / freqs)
    d_lT = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0] = log_T[0] - 0.5 * abs(d_lT[0]) if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else np.array([])
    y_edge[-1] = (
        log_T[-1] + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    )
    x_edge = (
        np.concatenate(
            [
                [x_cents[0] - 0.5],
                0.5 * (x_cents[:-1] + x_cents[1:]),
                [x_cents[-1] + 0.5],
            ]
        )
        if n_st > 1
        else np.array([-0.5, 0.5])
    )

    X, Y = np.meshgrid(x_edge, y_edge)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    qm = ax.pcolormesh(
        X,
        Y,
        mat,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        shading="flat",
        rasterized=True,
    )
    if period_up:
        ax.invert_yaxis()

    if show_threshold_contour:
        try:
            x_c = 0.5 * (x_edge[:-1] + x_edge[1:])
            y_c = 0.5 * (y_edge[:-1] + y_edge[1:])
            Xc, Yc = np.meshgrid(x_c, y_c)
            ax.contour(
                Xc,
                Yc,
                mat,
                levels=[score_threshold],
                colors=[contour_color],
                linewidths=[contour_lw],
            )
        except Exception:
            pass

    cb = plt.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(colorbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    # y-axis (log-period)
    pos = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(
            f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
        )
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    if station_markers:
        _apply_section_station_axis(
            ax, x_cents, st_order, (x_edge[0], x_edge[-1])
        )
    else:
        tick_cfg.apply(
            ax, x_cents, st_order, xlabel=xlabel, xlim=(x_edge[0], x_edge[-1])
        )

    if title:
        ax.set_title(
            title,
            fontsize=9,
            fontweight="bold",
            pad=22.0 if station_markers else None,
        )

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
    station_markers: bool = True,
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
    station_markers : bool, default ``True``
        Draw the station axis once, above the top panel, using
        pyCSAMT's shared downward-triangle convention (see
        :func:`plot_qc_heatmap`). Set to ``False`` for a plain bottom
        axis on the last panel, controlled by *tick_every* /
        *tick_label_rotation* / *tick_fontsize* / *station_tick_config*.
    tick_every, tick_label_rotation, tick_fontsize : tick control
        Only used when ``station_markers=False``.
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
        Only used when ``station_markers=False``.
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

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    _cmap_defaults = {
        "snr": "YlGn",
        "swift_skew": "YlOrRd_r",
        "asym": "RdBu_r",
        "phase_xy": "RdBu_r",
        "phase_yx": "RdBu_r",
        "score": "RdYlGn",
    }
    cmap_map = dict(_cmap_defaults)
    if cmaps:
        cmap_map.update(cmaps)

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(
            df["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )

    n_feat = len(features)
    n_st = len(st_order)
    freqs = np.sort(df["freq"].unique())[::-1]
    n_f = freqs.size

    if figsize is None:
        figsize = (max(8.0, n_st * 0.35), 2.8 * n_feat)

    fig, axes = plt.subplots(
        n_feat,
        1,
        figsize=figsize,
        sharex=True,
        layout="constrained",
        gridspec_kw={
            "hspace": 0.15 if station_markers else 0.45
        },
    )
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
    xe = (
        np.concatenate(
            [
                [x_cents[0] - 0.5],
                0.5 * (x_cents[:-1] + x_cents[1:]),
                [x_cents[-1] + 0.5],
            ]
        )
        if n_st > 1
        else np.array([-0.5, 0.5])
    )
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

        qm = ax.pcolormesh(
            X,
            Y,
            mat,
            cmap=cmap_map.get(feat, "RdBu_r"),
            vmin=vmin,
            vmax=vmax,
            shading="flat",
            rasterized=True,
        )
        if period_up:
            ax.invert_yaxis()

        cb = plt.colorbar(qm, ax=ax, fraction=0.018, pad=0.01)
        cb.set_label(_FEATURE_LABELS.get(feat, feat), fontsize=7)
        cb.ax.tick_params(labelsize=6)

        pos = np.linspace(log_T.min(), log_T.max(), n_yticks)
        ylab = []
        for v in pos:
            r = round(v)
            ylab.append(
                f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
            )
        ax.set_yticks(pos)
        ax.set_yticklabels(ylab, fontsize=tick_fontsize)
        ax.set_ylabel("Period (s)", fontsize=7)

    if station_markers:
        # One shared station axis above the top panel; the bottom
        # panel's own x-axis (already hidden by sharex on every panel
        # but the last) is switched off too, so nothing duplicates it.
        _apply_section_station_axis(
            axes[0], x_cents, st_order, (xe[0], xe[-1])
        )
        axes[-1].tick_params(
            axis="x", which="both", bottom=False, labelbottom=False
        )
    else:
        tick_cfg.apply(
            axes[-1], x_cents, st_order, xlabel="Station", xlim=(xe[0], xe[-1])
        )

    if title:
        if station_markers:
            # constrained_layout positions the suptitle itself, above
            # whatever room the top panel's station axis needs.
            fig.suptitle(title, fontsize=10, fontweight="bold")
        else:
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
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for prof in profiles:
        sc = pscores[prof]
        col = pcolors[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        label = prof if prof != "_" else "all"

        ax.hist(
            fin,
            bins=n_bins,
            range=(0, 1),
            color=col,
            alpha=fill_alpha,
            edgecolor="none",
            density=True,
            label=label,
        )

        if show_kde and fin.size > 3:
            try:
                from scipy.stats import gaussian_kde

                bw = kde_bw if kde_bw is not None else "scott"
                kde = gaussian_kde(fin, bw_method=bw)
                xs = np.linspace(0, 1, 300)
                ax.plot(xs, kde(xs), color=col, lw=line_lw)
            except ImportError:
                pass

    ax.axvline(
        score_threshold,
        color=threshold_color,
        lw=threshold_lw,
        ls=threshold_ls,
    )

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
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)

    vio_data: list[np.ndarray] = []
    vio_labels: list[str] = []
    vio_colors: list[str] = []

    for prof in profiles:
        sc = pscores[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        vio_data.append(fin)
        vio_labels.append(prof if prof != "_" else "all")
        vio_colors.append(pcolors[prof])

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    positions = np.arange(1, len(vio_data) + 1, dtype=float)
    rng = np.random.default_rng(0)

    if kind == "violin" and vio_data:
        parts = ax.violinplot(
            vio_data, positions=positions, showmedians=True, showextrema=True
        )
        for body, col in zip(parts["bodies"], vio_colors):
            body.set_facecolor(col)
            body.set_alpha(violin_alpha)
        for part_key in ("cmedians", "cmins", "cmaxes", "cbars"):
            if part_key in parts:
                parts[part_key].set_color("0.3")
                parts[part_key].set_linewidth(0.9)

    elif kind == "box" and vio_data:
        bp = ax.boxplot(
            vio_data,
            positions=positions,
            patch_artist=True,
            widths=0.5,
            notch=False,
            medianprops=dict(color="0.2", lw=1.5),
            whiskerprops=dict(lw=0.8),
            capprops=dict(lw=0.8),
            flierprops=dict(marker=".", ms=3, alpha=0.4),
        )
        for patch, col in zip(bp["boxes"], vio_colors):
            patch.set_facecolor(col)
            patch.set_alpha(violin_alpha)

    # jitter scatter overlay (always for "strip", optional for others)
    if kind == "strip" or show_scatter:
        for i, (sc, col) in enumerate(zip(vio_data, vio_colors)):
            jit = rng.uniform(-jitter, jitter, sc.size)
            ax.scatter(
                positions[i] + jit,
                sc,
                color=col,
                alpha=scatter_alpha,
                s=scatter_size,
                zorder=3,
            )

    ax.axhline(
        score_threshold,
        color=threshold_color,
        lw=threshold_lw,
        ls=threshold_ls,
        zorder=4,
    )
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
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
    suptitle : str
    figsize : (w, h), default ``(13.0, 9.5)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)
    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.1, 0.9],
        hspace=0.0,
        wspace=0.30,
    )
    ax_bar = fig.add_subplot(gs[0, :])
    ax_hist = fig.add_subplot(gs[1, 0])
    ax_vio = fig.add_subplot(gs[1, 1])

    # ── (a) bar chart ──────────────────────────────────────────────────────
    _render_bar_chart(
        ax_bar,
        pscores,
        plabels,
        pcolors,
        score_threshold=score_threshold,
        show_scatter=show_scatter,
        scatter_src=(scores if isinstance(scores, pd.DataFrame) else None),
        tick_cfg=tick_cfg,
        show_zone_labels=True,
    )
    ax_bar.set_title(
        "(a) Per-station QC scores", fontsize=9, fontweight="bold"
    )

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


# ─────────────────────────────────────────────────────────────────────────────
# plot_training_history
# (shared by EMDenoiser / AnomalyDetector / DimensionalityClassifier)
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_training_history(
    history: Any,
    *,
    train_key: str = "train_loss",
    val_key: str = "val_loss",
    train_label: str = "Train",
    val_label: str = "Validation",
    train_color: str = "#1f77b4",
    val_color: str = "#d62728",
    log_y: bool = False,
    mark_best: bool = True,
    best_color: str = "#2ca02c",
    xlabel: str = "Epoch",
    ylabel: str = "Loss",
    title: str = "",
    show_legend: bool = True,
    figsize: tuple[float, float] = (6.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Train/validation loss curve for a pycsamt.ai.processing estimator.

    Works with the ``history_`` property exposed after :meth:`fit` by
    :class:`~pycsamt.ai.processing.denoise.EMDenoiser`,
    :class:`~pycsamt.ai.processing.anomaly.AnomalyDetector`, and
    :class:`~pycsamt.ai.processing.classify.DimensionalityClassifier` —
    pass either the fitted estimator itself or its ``history_`` dict.

    Parameters
    ----------
    history : estimator or dict
        A fitted estimator exposing ``history_``, or a dict with
        *train_key* / *val_key* entries (lists of per-epoch loss
        values).  Estimators trained through a non-network fallback
        (scipy smoothing, PCA, random forest) have an empty history and
        raise :class:`ValueError`.
    train_key, val_key : str
        Dictionary keys holding the per-epoch loss lists.
    train_label, val_label : str
    train_color, val_color : str
    log_y : bool, default ``False``
        Use a log-scaled y-axis — useful when the loss spans more than
        one order of magnitude early in training.
    mark_best : bool, default ``True``
        Mark the lowest-validation-loss epoch with a vertical dashed
        line and an annotation.
    best_color : str
    xlabel, ylabel, title : str
    show_legend : bool, default ``True``
    figsize : (w, h), default ``(6.0, 4.0)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`

    Examples
    --------
    >>> from pycsamt.ai.processing import EMDenoiser
    >>> from pycsamt.ai.processing.plot import plot_training_history
    >>> den = EMDenoiser().fit(X_clean, epochs=60)  # doctest: +SKIP
    >>> plot_training_history(den, title="Denoiser training")  # doctest: +SKIP
    """
    hist = history.history_ if hasattr(history, "history_") else history
    if not isinstance(hist, dict):
        raise TypeError(
            "'history' must be a fitted estimator exposing 'history_' "
            "or a dict of per-epoch loss lists."
        )

    train_loss = hist.get(train_key, [])
    val_loss = hist.get(val_key, [])
    if not train_loss and not val_loss:
        raise ValueError(
            "Empty training history — the estimator was fitted with a "
            "non-network fallback (no epoch loop) or has not been "
            "fitted yet."
        )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    if train_loss:
        epochs = np.arange(1, len(train_loss) + 1)
        ax.plot(
            epochs,
            train_loss,
            color=train_color,
            lw=1.6,
            label=train_label,
        )
    if val_loss:
        epochs = np.arange(1, len(val_loss) + 1)
        ax.plot(
            epochs, val_loss, color=val_color, lw=1.6, label=val_label
        )

        if mark_best:
            best_ep = int(np.argmin(val_loss)) + 1
            best_val = val_loss[best_ep - 1]
            ax.axvline(
                best_ep, color=best_color, lw=1.0, ls="--", alpha=0.8
            )
            ax.annotate(
                f"best epoch {best_ep}\n"
                f"{val_label.lower()}={best_val:.4g}",
                xy=(best_ep, best_val),
                xytext=(8, 10),
                textcoords="offset points",
                fontsize=7.5,
                color=best_color,
            )

    if log_y:
        ax.set_yscale("log")

    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    if show_legend:
        ax.legend(fontsize=8, framealpha=0.9)

    return ax


# ─────────────────────────────────────────────────────────────────────────────
# EMDenoiser plots
# ─────────────────────────────────────────────────────────────────────────────

_DENOISE_COMP4 = [
    r"log$_{10}$|Z$_{xy}$|",
    r"$\phi_{xy}$ (°)",
    r"log$_{10}$|Z$_{yx}$|",
    r"$\phi_{yx}$ (°)",
]
_DENOISE_COMP8 = _DENOISE_COMP4 + [
    r"log$_{10}$|Z$_{xx}$|",
    r"$\phi_{xx}$ (°)",
    r"log$_{10}$|Z$_{yy}$|",
    r"$\phi_{yy}$ (°)",
]


def _denoise_component_labels(n_components: int) -> list[str]:
    if n_components == 8:
        return list(_DENOISE_COMP8)
    if n_components == 4:
        return list(_DENOISE_COMP4)
    return [f"ch{i}" for i in range(n_components)]


@EMStyle()
def plot_denoise_spectra(
    freq: np.ndarray,
    X_raw: np.ndarray,
    X_denoised: np.ndarray,
    *,
    station_labels: list[str] | None = None,
    component_labels: list[str] | None = None,
    sites: list[int] | None = None,
    n_show: int = 4,
    raw_color: str = "#7f7f7f",
    denoised_color: str = "#d62728",
    raw_alpha: float = 0.75,
    log_freq: bool = True,
    sharey: str = "row",
    suptitle: str = "",
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Before/after spectra for
    :class:`~pycsamt.ai.processing.denoise.EMDenoiser`.

    One column per selected site, one row per feature component,
    following the ``(n_sites, n_components, n_freqs)`` array
    convention used by
    :func:`~pycsamt.ai.processing.denoise.prepare_z_features`.

    Parameters
    ----------
    freq : ndarray, shape (n_freqs,)
        Frequency grid in Hz shared by all sites (the ``freq_ref``
        grid from :func:`~pycsamt.ai.processing.denoise.prepare_z_features`).
    X_raw, X_denoised : ndarray, shape (n_sites, n_components, n_freqs)
        Input and
        :meth:`~pycsamt.ai.processing.denoise.EMDenoiser.transform`
        output.
    station_labels : list of str or None
        One label per site.  Defaults to ``site0, site1, ...``.
    component_labels : list of str or None
        One label per row.  Defaults to the standard
        ``[log|Zxy|, phi_xy, log|Zyx|, phi_yx]`` ordering (or its
        8-component variant) used by :func:`prepare_z_features`.
    sites : list of int or None
        Explicit site indices to display.  Overrides *n_show* when
        given.
    n_show : int, default ``4``
        Number of sites to display when *sites* is ``None`` — the
        first *n_show* sites are used.
    raw_color, denoised_color : str
    raw_alpha : float, default ``0.75``
    log_freq : bool, default ``True``
    sharey : ``"row"`` | ``"none"``, default ``"row"``
        Share the y-axis across sites within each component row so
        amplitudes are directly comparable.
    suptitle : str
    figsize : (w, h) or None

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`

    Examples
    --------
    >>> from pycsamt.ai.processing.plot import plot_denoise_spectra
    >>> plot_denoise_spectra(
    ...     freq, X_noisy, X_denoised, n_show=3,
    ... )  # doctest: +SKIP
    """
    X_raw = np.asarray(X_raw, dtype=float)
    X_denoised = np.asarray(X_denoised, dtype=float)
    n_sites, n_comp, _n_f = X_raw.shape

    idx = (
        list(sites)
        if sites is not None
        else list(range(min(n_show, n_sites)))
    )
    labels = (
        station_labels
        if station_labels is not None
        else [f"site{i}" for i in range(n_sites)]
    )
    comp_labels = component_labels or _denoise_component_labels(n_comp)

    n_cols = len(idx)
    if figsize is None:
        figsize = (max(3.0, 3.1 * n_cols), 2.0 * n_comp)

    # NB: matplotlib's own sharey accepts "row"/"none" directly -- do
    # not coerce to bool, which would collapse to sharey="all" and mix
    # amplitude and phase rows onto one axis scale.
    fig, axes = plt.subplots(
        n_comp,
        n_cols,
        figsize=figsize,
        sharex=True,
        sharey=sharey,
        squeeze=False,
    )

    for c, si in enumerate(idx):
        for r in range(n_comp):
            ax = axes[r, c]
            ax.plot(
                freq,
                X_raw[si, r],
                color=raw_color,
                lw=1.1,
                alpha=raw_alpha,
                label="Input" if (r == 0 and c == 0) else None,
            )
            ax.plot(
                freq,
                X_denoised[si, r],
                color=denoised_color,
                lw=1.3,
                label="Denoised" if (r == 0 and c == 0) else None,
            )
            if log_freq:
                ax.set_xscale("log")
            if r == 0:
                ax.set_title(labels[si], fontsize=9, fontweight="bold")
            if c == 0:
                ax.set_ylabel(comp_labels[r], fontsize=8)
            if r == n_comp - 1:
                ax.set_xlabel("Frequency (Hz)", fontsize=8)
            ax.tick_params(labelsize=7)

    handles, hlabels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(
            handles,
            hlabels,
            loc="upper center",
            ncol=2,
            fontsize=8,
            bbox_to_anchor=(0.5, 1.04),
            frameon=False,
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold", y=1.08)

    fig.tight_layout()
    return fig


def _roughness(x: np.ndarray) -> np.ndarray:
    """
    Mean ``|second difference|`` along the last axis.

    A high-frequency-roughness proxy: large for jagged, noisy spectra,
    small for smooth ones.  Shared by :func:`plot_denoise_noise_reduction`.
    """
    d2 = np.diff(x, n=2, axis=-1)
    return np.nanmean(np.abs(d2), axis=-1)


@EMStyle()
def plot_denoise_noise_reduction(
    X_raw: np.ndarray,
    X_denoised: np.ndarray,
    *,
    station_labels: list[str] | None = None,
    aggregate: str = "mean",
    bar_color: str = "#2166ac",
    negative_color: str = "#c0392b",
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    xlabel: str = "Station",
    ylabel: str = "Roughness reduction (%)",
    title: str = "",
    show_grid: bool = True,
    figsize: tuple[float, float] = (10.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Per-station roughness reduction achieved by
    :class:`~pycsamt.ai.processing.denoise.EMDenoiser`.

    "Roughness" is the mean absolute second difference along the
    frequency axis — large for jagged, noisy spectra and small for
    smooth ones.  The plotted quantity is the percentage reduction in
    roughness from *X_raw* to *X_denoised*, aggregated over feature
    components with *aggregate*.  Positive values mean the denoiser
    removed high-frequency scatter; values near zero or negative flag
    stations where denoising had little effect or over-smoothed
    genuine structure.

    Parameters
    ----------
    X_raw, X_denoised : ndarray, shape (n_sites, n_components, n_freqs)
    station_labels : list of str or None
    aggregate : ``"mean"`` | ``"median"``, default ``"mean"``
        Reduction statistic across components.
    bar_color : str
        Colour for stations with a positive (improved) reduction.
    negative_color : str
        Colour for stations with zero or negative reduction.
    tick_every, tick_label_rotation, tick_fontsize : tick control
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
    xlabel, ylabel, title : str
    show_grid : bool, default ``True``
    figsize : (w, h), default ``(10.0, 4.0)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    X_raw = np.asarray(X_raw, dtype=float)
    X_denoised = np.asarray(X_denoised, dtype=float)
    n_sites = X_raw.shape[0]

    r_raw = _roughness(X_raw)
    r_den = _roughness(X_denoised)
    pct = 100.0 * (r_raw - r_den) / (r_raw + 1e-24)

    stat = np.nanmedian if aggregate == "median" else np.nanmean
    pct_site = stat(pct, axis=1)

    labels = (
        station_labels
        if station_labels is not None
        else [f"site{i}" for i in range(n_sites)]
    )
    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    xpos = np.arange(n_sites, dtype=float)
    colors = [bar_color if v >= 0 else negative_color for v in pct_site]
    ax.bar(xpos, pct_site, color=colors, width=0.8, alpha=0.85)
    ax.axhline(0, color="0.3", lw=0.8)

    tick_cfg.apply(
        ax, xpos, labels, xlabel=xlabel, xlim=(-0.8, n_sites - 0.2)
    )
    ax.set_ylabel(ylabel, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    if show_grid:
        ax.grid(True, axis="y", ls=":", lw=0.4, color="gray", alpha=0.5)
    ax.set_axisbelow(True)

    return ax


@EMStyle()
def plot_denoise_summary(
    freq: np.ndarray,
    X_raw: np.ndarray,
    X_denoised: np.ndarray,
    *,
    history: Any = None,
    station_labels: list[str] | None = None,
    component_labels: list[str] | None = None,
    n_show: int = 3,
    suptitle: str = "",
    figsize: tuple[float, float] = (13.0, 10.0),
) -> plt.Figure:
    """
    Combined denoiser figure: before/after spectra, per-station
    roughness reduction, and (when available) the training curve.

    Layout::

        ┌──────────────────────────────────────────┐
        │  Before / after spectra, n_show sites     │
        ├───────────────────────┬──────────────────┤
        │ Roughness reduction   │ Training history  │
        │ per station           │ (if available)    │
        └───────────────────────┴──────────────────┘

    Parameters
    ----------
    freq : ndarray, shape (n_freqs,)
    X_raw, X_denoised : ndarray, shape (n_sites, n_components, n_freqs)
    history : estimator, dict, or None
        Passed to :func:`plot_training_history` when not empty; the
        bottom-right panel is left blank otherwise.
    station_labels, component_labels : list of str or None
    n_show : int, default ``3``
    suptitle : str
    figsize : (w, h), default ``(13.0, 10.0)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    hist = None
    if history is not None:
        h = (
            history.history_
            if hasattr(history, "history_")
            else history
        )
        if isinstance(h, dict) and (
            h.get("train_loss") or h.get("val_loss")
        ):
            hist = h

    n_comp = X_raw.shape[1]
    n_cols = min(n_show, X_raw.shape[0])
    labels = (
        station_labels
        if station_labels is not None
        else [f"site{i}" for i in range(X_raw.shape[0])]
    )
    comp_labels = component_labels or _denoise_component_labels(n_comp)

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.5, 1], hspace=0.4)
    gs_top = gs[0, :].subgridspec(n_comp, n_cols)

    for c in range(n_cols):
        for r in range(n_comp):
            ax = fig.add_subplot(gs_top[r, c])
            ax.plot(
                freq,
                X_raw[c, r],
                color="#7f7f7f",
                lw=1.0,
                alpha=0.75,
                label="Input" if (r == 0 and c == 0) else None,
            )
            ax.plot(
                freq,
                X_denoised[c, r],
                color="#d62728",
                lw=1.2,
                label="Denoised" if (r == 0 and c == 0) else None,
            )
            ax.set_xscale("log")
            if r == 0:
                ax.set_title(labels[c], fontsize=9, fontweight="bold")
            if c == 0:
                ax.set_ylabel(comp_labels[r], fontsize=7.5)
            ax.tick_params(labelsize=6.5)

    ax_bar = fig.add_subplot(gs[1, 0])
    plot_denoise_noise_reduction.__wrapped__(
        X_raw,
        X_denoised,
        station_labels=station_labels,
        title="Roughness reduction",
        ax=ax_bar,
    )

    if hist is not None:
        ax_hist = fig.add_subplot(gs[1, 1])
        plot_training_history.__wrapped__(
            hist, title="Training history", ax=ax_hist
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# AnomalyDetector plots
# ─────────────────────────────────────────────────────────────────────────────


def _resolve_anomaly_threshold(
    pscores: dict[str, np.ndarray],
    threshold: float | None,
    threshold_percentile: float,
) -> float:
    if threshold is not None:
        return float(threshold)
    allv = np.concatenate(
        [v[np.isfinite(v)] for v in pscores.values()]
    )
    if allv.size == 0:
        return 0.0
    return float(np.percentile(allv, threshold_percentile))


@EMStyle()
def plot_anomaly_scores(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    threshold: float | None = None,
    threshold_percentile: float = 95.0,
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.4,
    threshold_ls: str = "--",
    bad_zone_color: str = "#fde8e8",
    bad_zone_alpha: float = 0.30,
    bar_width: float = 0.80,
    bar_alpha: float = 0.88,
    show_profile_labels: bool = True,
    profile_label_fontsize: int = 8,
    separator_color: str = "gray",
    separator_lw: float = 0.7,
    separator_ls: str = "--",
    show_zone_labels: bool = True,
    zone_label_fontsize: int = 8,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    xlabel: str = "Station",
    ylabel: str = "Anomaly score (reconstruction error)",
    title: str = "",
    ylim: tuple[float, float] | None = None,
    show_grid: bool = True,
    figsize: tuple[float, float] = (10.0, 4.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Per-station anomaly-score bar chart, mirroring :func:`plot_qc_scores`.

    Unlike QC scores, higher values are *worse*: the shaded rejection
    band and the "Anomalous" label sit above the threshold line rather
    than below it.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
        Reconstruction-error scores from
        :meth:`~pycsamt.ai.processing.anomaly.AnomalyDetector.transform`.
        Same flexible input as :func:`plot_qc_scores`: a single 1-D
        array, ``{profile: array}``, or a DataFrame with a ``score``
        column.
    station_labels : list of str or None
    profile_colors : dict, list, or None
    threshold : float or None
        Explicit anomaly threshold.  When ``None``, computed as the
        *threshold_percentile*-th percentile of the pooled scores.
        Pass
        :attr:`AnomalyDetector.threshold_
        <pycsamt.ai.processing.anomaly.AnomalyDetector.threshold_>`
        directly to plot the fitted value.
    threshold_percentile : float, default ``95.0``
    threshold_color, threshold_lw, threshold_ls : str, float, str
    bad_zone_color, bad_zone_alpha : str, float
    bar_width, bar_alpha : float
    show_profile_labels : bool, default ``True``
    profile_label_fontsize : int, default ``8``
    separator_color, separator_lw, separator_ls : str, float, str
    show_zone_labels : bool, default ``True``
    zone_label_fontsize : int, default ``8``
    tick_every, tick_label_rotation, tick_fontsize : tick control
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
    xlabel, ylabel, title : str
    ylim : (ymin, ymax) or None
        Defaults to ``(0, 1.15 * max(scores, threshold))``.
    show_grid : bool, default ``True``
    figsize : (w, h), default ``(10.0, 4.2)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`

    Examples
    --------
    >>> from pycsamt.ai.processing.plot import plot_anomaly_scores
    >>> plot_anomaly_scores(
    ...     scores, threshold=det.threshold_,
    ... )  # doctest: +SKIP
    """
    pscores, plabels = _normalise_qc_input(scores, station_labels)
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)
    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )
    thr = _resolve_anomaly_threshold(
        pscores, threshold, threshold_percentile
    )

    if ylim is None:
        allv = np.concatenate(
            [v[np.isfinite(v)] for v in pscores.values()]
        )
        ymax = max(float(np.nanmax(allv)) if allv.size else thr, thr)
        ylim = (0.0, 1.15 * max(ymax, 1e-12))

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _render_bar_chart(
        ax,
        pscores,
        plabels,
        pcolors,
        score_threshold=thr,
        show_scatter=False,
        scatter_src=None,
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
        bad_zone_side="above",
        ylim=ylim,
        ylabel=ylabel,
        zone_label_good="Normal ▶",
        zone_label_bad="◀ Anomalous",
    )

    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")

    if len(profiles) > 1 or profiles[0] != "_":
        handles = [
            mpatches.Patch(fc=pcolors[p], label=p, alpha=bar_alpha)
            for p in profiles
            if p != "_"
        ]
        handles.append(
            plt.Line2D(
                [],
                [],
                color=threshold_color,
                ls=threshold_ls,
                lw=threshold_lw,
                label=f"Threshold ({thr:.3g})",
            )
        )
        ax.legend(
            handles=handles,
            loc="upper right",
            fontsize=7.5,
            framealpha=0.9,
        )

    return ax


@EMStyle()
def plot_anomaly_score_distribution(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    threshold: float | None = None,
    threshold_percentile: float = 95.0,
    n_bins: int = 20,
    show_kde: bool = True,
    kde_bw: float | str | None = None,
    fill_alpha: float = 0.35,
    line_lw: float = 1.8,
    threshold_color: str = "#c0392b",
    threshold_lw: float = 1.3,
    threshold_ls: str = "--",
    xlabel: str = "Anomaly score",
    ylabel: str = "Density",
    title: str = "",
    show_legend: bool = True,
    figsize: tuple[float, float] = (6.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Anomaly-score histogram and optional KDE, mirroring
    :func:`plot_qc_score_distribution`.

    The shaded rejection band and threshold line sit to the right of
    the threshold value (high score = anomalous) — the opposite
    convention from the QC-score version.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
    station_labels : list of str or None
    profile_colors : dict, list, or None
    threshold : float or None
        See :func:`plot_anomaly_scores`.
    threshold_percentile : float, default ``95.0``
    n_bins : int, default ``20``
    show_kde : bool, default ``True``
    kde_bw : float, str, or None
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
    profiles = list(pscores.keys())
    pcolors = _resolve_profile_colors(profiles, profile_colors)
    thr = _resolve_anomaly_threshold(
        pscores, threshold, threshold_percentile
    )

    allv = np.concatenate([v[np.isfinite(v)] for v in pscores.values()])
    xmax = 1.05 * max(float(np.nanmax(allv)) if allv.size else thr, thr)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    for prof in profiles:
        sc = pscores[prof]
        col = pcolors[prof]
        fin = sc[np.isfinite(sc)]
        if not fin.size:
            continue
        label = prof if prof != "_" else "all"

        ax.hist(
            fin,
            bins=n_bins,
            range=(0, xmax),
            color=col,
            alpha=fill_alpha,
            edgecolor="none",
            density=True,
            label=label,
        )

        if show_kde and fin.size > 3:
            try:
                from scipy.stats import gaussian_kde

                bw = kde_bw if kde_bw is not None else "scott"
                kde = gaussian_kde(fin, bw_method=bw)
                xs = np.linspace(0, xmax, 300)
                ax.plot(xs, kde(xs), color=col, lw=line_lw)
            except ImportError:
                pass

    ax.axvline(
        thr, color=threshold_color, lw=threshold_lw, ls=threshold_ls
    )
    ax.axvspan(thr, xmax, color="#fde8e8", alpha=0.25, zorder=0)

    ax.set_xlim(0, xmax)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel(ylabel, fontsize=8)
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")
    if show_legend and (len(profiles) > 1 or profiles[0] != "_"):
        ax.legend(fontsize=7.5, framealpha=0.85)

    return ax


@EMStyle()
def plot_anomaly_summary(
    scores: Any,
    *,
    station_labels: list[str] | None = None,
    profile_colors: Any | None = None,
    threshold: float | None = None,
    threshold_percentile: float = 95.0,
    n_bins: int = 20,
    show_kde: bool = True,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    suptitle: str = "",
    figsize: tuple[float, float] = (13.0, 4.6),
) -> plt.Figure:
    """
    Two-panel anomaly-detection summary: per-station bar chart and
    score distribution, sharing one resolved threshold.

    Parameters
    ----------
    scores : ndarray, dict, or DataFrame
    station_labels : list of str or None
    profile_colors : dict, list, or None
    threshold : float or None
    threshold_percentile : float, default ``95.0``
    n_bins : int, default ``20``
    show_kde : bool, default ``True``
    tick_every, tick_label_rotation, tick_fontsize : tick control
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
    suptitle : str
    figsize : (w, h), default ``(13.0, 4.6)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    pscores, _ = _normalise_qc_input(scores, station_labels)
    thr = _resolve_anomaly_threshold(
        pscores, threshold, threshold_percentile
    )

    fig, (ax_bar, ax_hist) = plt.subplots(
        1,
        2,
        figsize=figsize,
        gridspec_kw={"width_ratios": [1.6, 1]},
        layout="constrained",
    )

    plot_anomaly_scores.__wrapped__(
        scores,
        station_labels=station_labels,
        profile_colors=profile_colors,
        threshold=thr,
        tick_every=tick_every,
        tick_label_rotation=tick_label_rotation,
        tick_fontsize=tick_fontsize,
        station_tick_config=station_tick_config,
        title="Per-station anomaly scores",
        ax=ax_bar,
    )
    plot_anomaly_score_distribution.__wrapped__(
        scores,
        station_labels=station_labels,
        profile_colors=profile_colors,
        threshold=thr,
        n_bins=n_bins,
        show_kde=show_kde,
        title="Score distribution",
        ax=ax_hist,
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# DimensionalityClassifier plots
# ─────────────────────────────────────────────────────────────────────────────

_DIM_COLORS = ("#2ca02c", "#1f77b4", "#d62728")  # 1-D, 2-D, 3-D
_DIM_LABELS = ("1-D", "2-D", "3-D")


@EMStyle()
def plot_dimensionality_map(
    df: pd.DataFrame,
    *,
    dim_col: str = "dim",
    station_labels: list[str] | None = None,
    class_colors: tuple[str, str, str] = _DIM_COLORS,
    class_labels: tuple[str, str, str] = _DIM_LABELS,
    period_up: bool = True,
    n_yticks: int = 7,
    station_markers: bool = True,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    xlabel: str = "Station",
    ylabel: str = "Period (s)",
    title: str = "",
    figsize: tuple[float, float] = (10.0, 5.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Station × frequency dimensionality-class heat-map.

    Renders the categorical 1-D / 2-D / 3-D output of
    :meth:`DimensionalityClassifier.predict_table
    <pycsamt.ai.processing.classify.DimensionalityClassifier.predict_table>`
    the same way :func:`plot_qc_heatmap` renders continuous QC scores,
    with a discrete three-colour legend in place of a colourbar.

    Parameters
    ----------
    df : DataFrame
        Output of ``DimensionalityClassifier.predict_table()``.  Must
        contain ``station``, ``freq``, and *dim_col*.
    dim_col : str, default ``"dim"``
        Column with integer class labels (0=1-D, 1=2-D, 2=3-D).
    station_labels : list of str or None
    class_colors : (str, str, str)
    class_labels : (str, str, str)
    period_up : bool, default ``True``
    n_yticks : int, default ``7``
    station_markers : bool, default ``True``
        Draw the station axis at the top with pyCSAMT's shared
        downward-triangle convention (see :func:`plot_qc_heatmap`).
        Set to ``False`` for a plain bottom axis controlled by
        *tick_every* / *tick_label_rotation* / *tick_fontsize* /
        *station_tick_config*.
    tick_every, tick_label_rotation, tick_fontsize : tick control
        Only used when ``station_markers=False``.
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
        Only used when ``station_markers=False``.
    xlabel, ylabel, title : str
    figsize : (w, h), default ``(10.0, 5.2)``
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not {"station", "freq", dim_col}.issubset(df.columns):
        raise ValueError(
            f"DataFrame must have columns 'station', 'freq', "
            f"{dim_col!r} (from "
            "DimensionalityClassifier.predict_table())."
        )

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(
            df["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )

    mat, freqs, x_cents = _station_freq_grid(
        df, st_order, value_col=dim_col
    )
    n_f, n_st = mat.shape

    log_T = np.log10(1.0 / freqs)
    d_lT = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0] = log_T[0] - 0.5 * abs(d_lT[0]) if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else np.array([])
    y_edge[-1] = (
        log_T[-1] + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    )
    x_edge = (
        np.concatenate(
            [
                [x_cents[0] - 0.5],
                0.5 * (x_cents[:-1] + x_cents[1:]),
                [x_cents[-1] + 0.5],
            ]
        )
        if n_st > 1
        else np.array([-0.5, 0.5])
    )
    X, Y = np.meshgrid(x_edge, y_edge)

    cmap = ListedColormap(list(class_colors))
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    ax.pcolormesh(
        X, Y, mat, cmap=cmap, norm=norm, shading="flat", rasterized=True
    )
    if period_up:
        ax.invert_yaxis()

    pos = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(
            f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
        )
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    if station_markers:
        _apply_section_station_axis(
            ax, x_cents, st_order, (x_edge[0], x_edge[-1])
        )
    else:
        tick_cfg.apply(
            ax, x_cents, st_order, xlabel=xlabel, xlim=(x_edge[0], x_edge[-1])
        )

    # Below the plot rather than inside it -- an in-axes legend would
    # otherwise sit under the station triangles when station_markers=True.
    handles = [
        mpatches.Patch(fc=c, label=lab)
        for c, lab in zip(class_colors, class_labels)
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.07 if station_markers else -0.12),
        fontsize=7.5,
        framealpha=0.9,
        ncol=3,
    )

    if title:
        ax.set_title(
            title,
            fontsize=9,
            fontweight="bold",
            pad=22.0 if station_markers else None,
        )

    return ax


@EMStyle()
def plot_predicted_strike_rose(
    strike_deg: Any,
    *,
    bins: int = 18,
    fold_180: bool = True,
    color: str = "#1f77b4",
    edgecolor: str = "white",
    alpha: float = 0.85,
    show_mean: bool = True,
    mean_color: str = "#d62728",
    title: str = "",
    figsize: tuple[float, float] = (5.0, 5.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Rose diagram of predicted geoelectric strike directions.

    Parameters
    ----------
    strike_deg : array-like or DataFrame
        Strike angles in degrees, typically the ``strike`` column
        returned by ``DimensionalityClassifier.predict_table()``
        (``NaN`` for non-2-D observations).  A DataFrame is expected
        to carry a ``strike`` column.
    bins : int, default ``18``
        Number of angular bins spanning 0-360 degrees (10-degree bins
        by default).
    fold_180 : bool, default ``True``
        Strike has a 180-degree ambiguity (a line has no direction of
        travel).  When ``True``, every angle is plotted together with
        its 180-degree-rotated counterpart so the rose is symmetric —
        the standard convention for strike roses.
    color, edgecolor, alpha : str, str, float
    show_mean : bool, default ``True``
        Draw a radial line at the circular mean strike.
    mean_color : str
    title : str
    figsize : (w, h), default ``(5.0, 5.0)``
    ax : polar Axes or None
        When given, must have been created with
        ``subplot_kw={"projection": "polar"}``.

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes` (polar)

    Examples
    --------
    >>> from pycsamt.ai.processing.plot import plot_predicted_strike_rose
    >>> plot_predicted_strike_rose(
    ...     table["strike"], title="L18 strike rose",
    ... )  # doctest: +SKIP
    """
    if isinstance(strike_deg, pd.DataFrame):
        strike_deg = strike_deg["strike"]
    vals = np.asarray(strike_deg, dtype=float)
    vals = vals[np.isfinite(vals)]

    if ax is None:
        _, ax = plt.subplots(
            figsize=figsize, subplot_kw={"projection": "polar"}
        )

    if vals.size == 0:
        ax.set_title(title or "No 2-D strike estimates", fontsize=9)
        return ax

    angles = np.deg2rad(vals % 180.0)
    if fold_180:
        angles = np.concatenate([angles, angles + np.pi])

    edges = np.linspace(0, 2 * np.pi, bins + 1)
    counts, _ = np.histogram(angles, bins=edges)
    width = edges[1] - edges[0]
    centers = 0.5 * (edges[:-1] + edges[1:])

    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.bar(
        centers,
        counts,
        width=width,
        color=color,
        edgecolor=edgecolor,
        alpha=alpha,
        linewidth=0.6,
    )

    if show_mean:
        theta = np.deg2rad(vals)
        mean_rad = 0.5 * np.arctan2(
            np.mean(np.sin(2 * theta)), np.mean(np.cos(2 * theta))
        )
        mean_deg = np.rad2deg(mean_rad) % 180.0
        for a in (mean_deg, mean_deg + 180.0):
            ax.plot(
                [np.deg2rad(a)] * 2,
                [0, counts.max()],
                color=mean_color,
                lw=1.6,
                ls="--",
            )

    ax.set_title(title, fontsize=9, fontweight="bold")
    ax.tick_params(labelsize=7)

    return ax


@EMStyle()
def plot_dimensionality_summary(
    df: pd.DataFrame,
    *,
    dim_col: str = "dim",
    strike_col: str = "strike",
    station_labels: list[str] | None = None,
    class_colors: tuple[str, str, str] = _DIM_COLORS,
    class_labels: tuple[str, str, str] = _DIM_LABELS,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    rose_bins: int = 18,
    suptitle: str = "",
    figsize: tuple[float, float] = (13.0, 9.0),
) -> plt.Figure:
    """
    Three-panel dimensionality-classification summary.

    Layout::

        ┌────────────────────────────────────┐
        │  Station × frequency dim map        │  full width
        ├───────────────────┬────────────────┤
        │ Class proportions  │ Strike rose    │
        │                    │ (2-D sites)    │
        └───────────────────┴────────────────┘

    Parameters
    ----------
    df : DataFrame
        Output of ``DimensionalityClassifier.predict_table()``.
    dim_col, strike_col : str
    station_labels : list of str or None
    class_colors, class_labels : (str, str, str)
    tick_every, tick_label_rotation, tick_fontsize : tick control
    station_tick_config : :class:`~pycsamt.ai.plot.StationTickConfig` or None
    rose_bins : int, default ``18``
    suptitle : str
    figsize : (w, h), default ``(13.0, 9.0)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(
        2, 2, height_ratios=[1.05, 0.95], hspace=0.0, wspace=0.22
    )

    ax_map = fig.add_subplot(gs[0, :])
    plot_dimensionality_map.__wrapped__(
        df,
        dim_col=dim_col,
        station_labels=station_labels,
        class_colors=class_colors,
        class_labels=class_labels,
        tick_every=tick_every,
        tick_label_rotation=tick_label_rotation,
        tick_fontsize=tick_fontsize,
        station_tick_config=station_tick_config,
        title="(a) Dimensionality map",
        ax=ax_map,
    )

    ax_bar = fig.add_subplot(gs[1, 0])
    counts = df[dim_col].value_counts().reindex([0, 1, 2], fill_value=0)
    frac = 100.0 * counts / max(int(counts.sum()), 1)
    ax_bar.bar(class_labels, frac.to_numpy(), color=class_colors, alpha=0.85)
    ax_bar.set_ylabel("Share of samples (%)", fontsize=8)
    ax_bar.set_title("(b) Class proportions", fontsize=9, fontweight="bold")
    for x, v in enumerate(frac.to_numpy()):
        ax_bar.text(x, v + 1.0, f"{v:.1f}%", ha="center", fontsize=7.5)
    ax_bar.set_ylim(0, max(float(frac.max()) * 1.2, 10.0))

    ax_rose = fig.add_subplot(gs[1, 1], projection="polar")
    plot_predicted_strike_rose.__wrapped__(
        df[strike_col] if strike_col in df.columns else np.array([]),
        bins=rose_bins,
        title="(c) 2-D strike rose",
        ax=ax_rose,
    )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# EMImputer plots
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_imputer_gaps(
    gap_table: pd.DataFrame,
    *,
    status_col: str = "missing",
    station_labels: list[str] | None = None,
    observed_color: str = "#2166ac",
    missing_color: str = "#d62728",
    period_up: bool = True,
    colorbar_label: str = "Status",
    xlabel: str = "Station",
    ylabel: str = "Period (s)",
    title: str = "",
    station_markers: bool = True,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    n_yticks: int = 7,
    figsize: tuple[float, float] = (10.0, 5.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Station x frequency map of genuinely missing cells for
    :class:`~pycsamt.ai.processing.imputer.EMImputer`.

    Parameters
    ----------
    gap_table : DataFrame
        Must contain ``station``, ``freq``, and *status_col* columns.
        *status_col* is ``1`` where the (station, frequency) cell is
        genuinely missing (any tensor component ``NaN``) and ``0``
        where every component is observed -- build it from a site
        collection with
        ``numpy.isnan(z).any(axis=(1, 2))`` per station, one row per
        frequency (see the ``imputer`` user_guide page for a worked
        example).
    station_labels : list of str or None
    observed_color, missing_color : str
    period_up : bool, default ``True``
    colorbar_label, xlabel, ylabel, title : str
    station_markers : bool, default ``True``
        Draw the station axis at the top with pyCSAMT's shared
        downward-triangle convention (
        :data:`~pycsamt.api.station.PYCSAMT_STATION_RENDERING`). Set
        to ``False`` to fall back to a plain bottom axis controlled by
        *tick_every* / *tick_label_rotation* / *tick_fontsize* /
        *station_tick_config* instead.
    tick_every, tick_label_rotation, tick_fontsize, station_tick_config
        Only used when ``station_markers=False``.
    n_yticks : int, default ``7``
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not isinstance(gap_table, pd.DataFrame):
        raise TypeError(
            "plot_imputer_gaps expects a DataFrame with columns "
            "'station', 'freq', and status_col."
        )
    if not {"station", "freq", status_col}.issubset(gap_table.columns):
        raise ValueError(
            f"DataFrame must have columns 'station', 'freq', {status_col!r}."
        )

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(
            gap_table["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )

    mat, freqs, x_cents = _station_freq_grid(
        gap_table, st_order, value_col=status_col
    )
    n_f, n_st = mat.shape

    log_T = np.log10(1.0 / freqs)
    d_lT = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0] = log_T[0] - 0.5 * abs(d_lT[0]) if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else np.array([])
    y_edge[-1] = (
        log_T[-1] + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    )
    x_edge = (
        np.concatenate(
            [
                [x_cents[0] - 0.5],
                0.5 * (x_cents[:-1] + x_cents[1:]),
                [x_cents[-1] + 0.5],
            ]
        )
        if n_st > 1
        else np.array([-0.5, 0.5])
    )

    X, Y = np.meshgrid(x_edge, y_edge)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    cmap = ListedColormap([observed_color, missing_color])
    norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)
    qm = ax.pcolormesh(
        X, Y, mat, cmap=cmap, norm=norm, shading="flat", rasterized=True
    )
    if period_up:
        ax.invert_yaxis()

    cb = plt.colorbar(qm, ax=ax, fraction=0.02, pad=0.01, ticks=[0, 1])
    cb.ax.set_yticklabels(["Observed", "Missing"], fontsize=7)
    cb.set_label(colorbar_label, fontsize=8)

    pos = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(
            f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
        )
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    if station_markers:
        _apply_section_station_axis(
            ax, x_cents, st_order, (x_edge[0], x_edge[-1])
        )
    else:
        tick_cfg.apply(
            ax, x_cents, st_order, xlabel=xlabel, xlim=(x_edge[0], x_edge[-1])
        )

    if title:
        ax.set_title(
            title,
            fontsize=9,
            fontweight="bold",
            pad=22.0 if station_markers else None,
        )

    return ax


def _validation_scatter_panel(
    ax: plt.Axes,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    color: str = "#2166ac",
    point_size: float = 10.0,
    alpha: float = 0.5,
    show_metrics: bool = True,
    xlabel: str = "True value",
    ylabel: str = "Reconstructed value",
    title: str = "",
) -> None:
    """Shared single-panel true-vs-predicted scatter + 1:1 line +
    RMSE/R^2 annotation, one unit system at a time."""
    ax.scatter(
        y_true, y_pred, s=point_size, alpha=alpha, color=color,
        edgecolors="none",
    )
    if len(y_true):
        lo = float(min(y_true.min(), y_pred.min()))
        hi = float(max(y_true.max(), y_pred.max()))
        pad = 0.03 * (hi - lo + 1e-9)
        ax.plot(
            [lo - pad, hi + pad], [lo - pad, hi + pad],
            color="#333333", lw=1.0, ls="--",
        )
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_ylim(lo - pad, hi + pad)

    if show_metrics and len(y_true) > 1:
        rmse = float(np.sqrt(np.mean((y_true - y_pred) ** 2)))
        ss_res = float(np.sum((y_true - y_pred) ** 2))
        ss_tot = float(np.sum((y_true - y_true.mean()) ** 2)) + 1e-24
        r2 = 1.0 - ss_res / ss_tot
        ax.text(
            0.04, 0.96, f"RMSE = {rmse:.3g}\n$R^2$ = {r2:.3f}",
            transform=ax.transAxes, fontsize=7.5, va="top", ha="left",
            bbox={
                "boxstyle": "round", "fc": "white", "ec": "#888888",
                "alpha": 0.85,
            },
        )

    ax.set_xlabel(xlabel, fontsize=8.5)
    ax.set_ylabel(ylabel, fontsize=8.5)
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")
    ax.tick_params(labelsize=7.5)
    ax.set_aspect("equal", adjustable="box")


@EMStyle()
def plot_imputer_validation(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    component_ids: np.ndarray | None = None,
    component_labels: list[str] | None = None,
    colors: list[str] | None = None,
    n_cols: int = 2,
    point_size: float = 10.0,
    alpha: float = 0.5,
    show_metrics: bool = True,
    xlabel: str = "True value",
    ylabel: str = "Reconstructed value",
    suptitle: str = "",
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Predicted-vs-true scatter for a synthetic held-out-cell validation
    of :class:`~pycsamt.ai.processing.imputer.EMImputer`.

    Genuinely missing cells have no known ground truth to compare
    against, so the standard masked-reconstruction validation
    protocol instead hides a known fraction of *already observed*
    cells, reconstructs them, and compares the reconstruction against
    the real, held-back values -- exactly what this scatter plots.

    Parameters
    ----------
    y_true, y_pred : ndarray, shape (n_points,)
        Flattened true / reconstructed values at synthetically hidden
        cells only.
    component_ids : ndarray of int, shape (n_points,), or None
        Per-point feature-component index. When given, one panel is
        drawn *per component* rather than one shared-axis scatter --
        log-amplitude and phase channels live on very different
        scales, and a single pooled RMSE / :math:`R^2` (or a single
        shared axis range) would be dominated by whichever channel
        has the largest numeric range, the same pooled-metric pitfall
        documented for :class:`~pycsamt.ai.processing.denoise.EMDenoiser`
        in :doc:`/user_guide/ai_processing/denoise`.
    component_labels : list of str or None
        Panel title for each ``component_ids`` value.
    colors : list of str or None
    n_cols : int, default ``2``
        Panel-grid columns when faceting by component.
    point_size : float, default ``10.0``
    alpha : float, default ``0.5``
    show_metrics : bool, default ``True``
        Annotate RMSE and :math:`R^2` on each panel.
    xlabel, ylabel : str
    suptitle : str
    figsize : (w, h) or None

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    finite = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true, y_pred = y_true[finite], y_pred[finite]

    if component_ids is None:
        if figsize is None:
            figsize = (5.0, 5.0)
        fig, ax = plt.subplots(figsize=figsize)
        _validation_scatter_panel(
            ax, y_true, y_pred, point_size=point_size, alpha=alpha,
            show_metrics=show_metrics, xlabel=xlabel, ylabel=ylabel,
        )
        if suptitle:
            fig.suptitle(suptitle, fontsize=10, fontweight="bold")
        fig.tight_layout()
        return fig

    cids = np.asarray(component_ids).ravel()[finite]
    uniq = np.unique(cids)
    n_panels = len(uniq)
    n_cols = max(1, min(n_cols, n_panels))
    n_rows = int(np.ceil(n_panels / n_cols))
    if figsize is None:
        figsize = (3.6 * n_cols, 3.6 * n_rows)

    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=figsize, squeeze=False,
    )

    for i, cid in enumerate(uniq):
        ax = axes[i // n_cols, i % n_cols]
        m = cids == cid
        lbl = (
            component_labels[int(cid)]
            if component_labels is not None
            and int(cid) < len(component_labels)
            else f"component {int(cid)}"
        )
        c = (
            colors[i % len(colors)]
            if colors
            else _PROFILE_COLORS[i % len(_PROFILE_COLORS)]
        )
        _validation_scatter_panel(
            ax, y_true[m], y_pred[m], color=c, point_size=point_size,
            alpha=alpha, show_metrics=show_metrics, xlabel=xlabel,
            ylabel=ylabel, title=lbl,
        )

    for j in range(n_panels, n_rows * n_cols):
        axes[j // n_cols, j % n_cols].set_visible(False)

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")
    fig.tight_layout()
    return fig


@EMStyle()
def plot_imputer_reconstruction(
    freq: np.ndarray,
    X_raw: np.ndarray,
    X_filled: np.ndarray,
    *,
    missing_mask: np.ndarray | None = None,
    station_labels: list[str] | None = None,
    component_labels: list[str] | None = None,
    sites: list[int] | None = None,
    n_show: int = 4,
    line_color: str = "#2166ac",
    fill_marker_color: str = "#d62728",
    fill_marker_size: float = 22.0,
    log_freq: bool = True,
    sharey: str = "row",
    suptitle: str = "",
    figsize: tuple[float, float] | None = None,
) -> plt.Figure:
    """
    Per-station, per-channel spectra with reconstructed cells
    highlighted, for :class:`~pycsamt.ai.processing.imputer.EMImputer`.

    One column per selected site, one row per feature component,
    following the ``(n_sites, n_components, n_freqs)`` array
    convention used by
    :func:`~pycsamt.ai.processing.denoise.prepare_z_features`. The
    continuous curve is *X_filled* -- observed and reconstructed
    values alike -- with a marker overlaid at every genuinely missing
    cell, so the reconstructed segment is visually distinct from the
    real, measured backbone around it.

    Parameters
    ----------
    freq : ndarray, shape (n_freqs,)
        Frequency grid in Hz shared by all sites.
    X_raw : ndarray, shape (n_sites, n_components, n_freqs)
        Data *before* imputation -- only used to infer *missing_mask*
        when not given explicitly (``~numpy.isfinite(X_raw)``).
    X_filled : ndarray, shape (n_sites, n_components, n_freqs)
        :meth:`~pycsamt.ai.processing.imputer.EMImputer.transform`
        output -- the curve actually plotted.
    missing_mask : ndarray of bool, same shape, or None
    station_labels : list of str or None
    component_labels : list of str or None
    sites : list of int or None
        Explicit site indices to display. Overrides *n_show* when
        given.
    n_show : int, default ``4``
        Number of sites to display when *sites* is ``None`` -- sites
        with at least one missing cell are preferred, first *n_show*
        such sites are used (falling back to the first *n_show* sites
        overall when fewer than *n_show* have any gap).
    line_color, fill_marker_color : str
    fill_marker_size : float, default ``22.0``
    log_freq : bool, default ``True``
    sharey : ``"row"`` | ``"none"``, default ``"row"``
    suptitle : str
    figsize : (w, h) or None

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    X_raw = np.asarray(X_raw, dtype=float)
    X_filled = np.asarray(X_filled, dtype=float)
    n_sites, n_comp, _n_f = X_raw.shape
    if missing_mask is None:
        missing_mask = ~np.isfinite(X_raw)
    else:
        missing_mask = np.asarray(missing_mask, dtype=bool)

    if sites is not None:
        idx = list(sites)
    else:
        has_gap = [
            i for i in range(n_sites) if missing_mask[i].any()
        ]
        idx = (has_gap + [
            i for i in range(n_sites) if i not in has_gap
        ])[: min(n_show, n_sites)]

    labels = (
        station_labels
        if station_labels is not None
        else [f"site{i}" for i in range(n_sites)]
    )
    comp_labels = component_labels or _denoise_component_labels(n_comp)

    n_cols = len(idx)
    if figsize is None:
        figsize = (max(3.0, 3.1 * n_cols), 2.0 * n_comp)

    fig, axes = plt.subplots(
        n_comp,
        n_cols,
        figsize=figsize,
        sharex=True,
        sharey=sharey,
        squeeze=False,
    )

    for c, si in enumerate(idx):
        for r in range(n_comp):
            ax = axes[r, c]
            ax.plot(
                freq,
                X_filled[si, r],
                color=line_color,
                lw=1.2,
                label="Filled" if (r == 0 and c == 0) else None,
            )
            miss = missing_mask[si, r]
            if miss.any():
                ax.scatter(
                    freq[miss],
                    X_filled[si, r][miss],
                    color=fill_marker_color,
                    s=fill_marker_size,
                    zorder=5,
                    label="Reconstructed",
                )
            if log_freq:
                ax.set_xscale("log")
            if r == 0:
                ax.set_title(labels[si], fontsize=9, fontweight="bold")
            if c == 0:
                ax.set_ylabel(comp_labels[r], fontsize=8)
            if r == n_comp - 1:
                ax.set_xlabel("Frequency (Hz)", fontsize=8)
            ax.tick_params(labelsize=7)

    handles: list = []
    hlabels: list = []
    for a in axes.ravel():
        h, ll = a.get_legend_handles_labels()
        for hh, lbl in zip(h, ll):
            if lbl not in hlabels:
                handles.append(hh)
                hlabels.append(lbl)
    if handles:
        fig.legend(
            handles,
            hlabels,
            loc="upper center",
            ncol=2,
            fontsize=8,
            bbox_to_anchor=(0.5, 1.04),
            frameon=False,
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold", y=1.08)

    fig.tight_layout()
    return fig


@EMStyle()
def plot_imputer_summary(
    gap_table: pd.DataFrame,
    freq: np.ndarray,
    X_raw: np.ndarray,
    X_filled: np.ndarray,
    *,
    missing_mask: np.ndarray | None = None,
    history: Any = None,
    station_labels: list[str] | None = None,
    component_labels: list[str] | None = None,
    n_show: int = 2,
    suptitle: str = "",
    figsize: tuple[float, float] = (12.0, 11.0),
) -> plt.Figure:
    """
    Combined :class:`~pycsamt.ai.processing.imputer.EMImputer` figure:
    the missing-data map, reconstruction detail for a handful of
    sites, and (when available) the training curve.

    Layout::

        ┌──────────────────────────────────────────────┐
        │  (a) Station x frequency missing-data map     │
        ├──────────────────────────────────────────────┤
        │  (b) Reconstruction detail, n_show sites      │
        ├──────────────────────────────────────────────┤
        │  (c) Training history                         │
        └──────────────────────────────────────────────┘

    The held-out-cell validation scatter is deliberately *not* one of
    these panels -- see :func:`plot_imputer_validation` and pass it
    the same ``y_true``/``y_pred`` separately, faceted one panel per
    feature component, the same reason
    :func:`plot_denoise_spectra` never pools log-amplitude and phase
    into one shared-axis view.

    Parameters
    ----------
    gap_table : DataFrame
        Passed to :func:`plot_imputer_gaps`.
    freq, X_raw, X_filled, missing_mask, station_labels, component_labels
        Passed to :func:`plot_imputer_reconstruction`.
    history : estimator, dict, or None
        Passed to :func:`plot_training_history` when not empty; panel
        (c) is left blank otherwise.
    n_show : int, default ``2``
    suptitle : str
    figsize : (w, h), default ``(12.0, 11.0)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    hist = None
    if history is not None:
        h = history.history_ if hasattr(history, "history_") else history
        if isinstance(h, dict) and (
            h.get("train_loss") or h.get("val_loss")
        ):
            hist = h

    X_raw = np.asarray(X_raw, dtype=float)
    X_filled = np.asarray(X_filled, dtype=float)
    n_comp = X_raw.shape[1]
    n_cols = min(n_show, X_raw.shape[0])
    if missing_mask is None:
        missing_mask = ~np.isfinite(X_raw)
    else:
        missing_mask = np.asarray(missing_mask, dtype=bool)
    labels = (
        station_labels
        if station_labels is not None
        else [f"site{i}" for i in range(X_raw.shape[0])]
    )
    comp_labels = component_labels or _denoise_component_labels(n_comp)

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(3, 1, height_ratios=[1.0, 1.3, 0.85])

    ax_gap = fig.add_subplot(gs[0, 0])
    plot_imputer_gaps.__wrapped__(
        gap_table,
        station_labels=station_labels,
        title="(a) Missing-data map",
        ax=ax_gap,
    )

    has_gap = [i for i in range(X_raw.shape[0]) if missing_mask[i].any()]
    idx = (has_gap + [
        i for i in range(X_raw.shape[0]) if i not in has_gap
    ])[:n_cols]

    gs_mid = gs[1, 0].subgridspec(
        n_comp, n_cols, hspace=0.15, wspace=0.25
    )
    for c, si in enumerate(idx):
        for r in range(n_comp):
            ax = fig.add_subplot(gs_mid[r, c])
            ax.plot(freq, X_filled[si, r], color="#2166ac", lw=1.1)
            miss = missing_mask[si, r]
            if miss.any():
                ax.scatter(
                    freq[miss],
                    X_filled[si, r][miss],
                    color="#d62728",
                    s=18.0,
                    zorder=5,
                )
            ax.set_xscale("log")
            if r == 0:
                title = labels[si]
                if c == 0:
                    title = "(b) Reconstruction detail -- " + title
                ax.set_title(title, fontsize=9, fontweight="bold")
            if c == 0:
                ax.set_ylabel(comp_labels[r], fontsize=7.5)
            ax.tick_params(labelsize=6.5)

    if hist is not None:
        ax_hist = fig.add_subplot(gs[2, 0])
        plot_training_history.__wrapped__(
            hist, title="(c) Training history", ax=ax_hist
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# UncertaintyCalibrator plots
# ─────────────────────────────────────────────────────────────────────────────


@EMStyle()
def plot_uncertainty_map(
    table: pd.DataFrame,
    *,
    value_col: str = "z_err_frac_calibrated",
    station_labels: list[str] | None = None,
    cmap: str = "viridis",
    vmin: float | None = None,
    vmax: float | None = None,
    period_up: bool = True,
    n_yticks: int = 7,
    colorbar_label: str = "Fractional error",
    xlabel: str = "Station",
    ylabel: str = "Period (s)",
    title: str = "",
    station_markers: bool = True,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    figsize: tuple[float, float] = (10.0, 5.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Station x frequency heat-map of a fractional-error column, for
    :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`.

    Colour is log-scaled: fractional error routinely spans more than a
    decade within one survey, and a linear scale would wash out the
    low end.

    Parameters
    ----------
    table : DataFrame
        Output of
        :meth:`UncertaintyCalibrator.predict_table
        <pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.predict_table>`
        (or
        :func:`~pycsamt.ai.processing.uncertainty.build_uncertainty_features_table`).
        Must contain ``station``, ``freq``, and *value_col* columns.
    value_col : str, default ``"z_err_frac_calibrated"``
        Column to plot -- pass ``"z_err_frac"`` instead for the
        original field-processing fractional error, e.g. to build a
        matching before/after pair of figures.
    station_labels : list of str or None
    cmap : str, default ``"viridis"``
    vmin, vmax : float or None
        Log-color-scale bounds; ``None`` infers both from the data.
    period_up : bool, default ``True``
    n_yticks : int, default ``7``
    colorbar_label, xlabel, ylabel, title : str
    station_markers : bool, default ``True``
        See :func:`plot_qc_heatmap`.
    tick_every, tick_label_rotation, tick_fontsize, station_tick_config
        Only used when ``station_markers=False``.
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError(
            "plot_uncertainty_map expects a DataFrame with columns "
            "'station', 'freq', and value_col."
        )
    if not {"station", "freq", value_col}.issubset(table.columns):
        raise ValueError(
            f"DataFrame must have columns 'station', 'freq', {value_col!r}."
        )

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if station_labels is not None:
        st_order = station_labels
    else:
        st_order = sorted(
            table["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )

    mat, freqs, x_cents = _station_freq_grid(
        table, st_order, value_col=value_col
    )
    n_f, n_st = mat.shape

    log_T = np.log10(1.0 / freqs)
    d_lT = np.diff(log_T)
    y_edge = np.empty(n_f + 1)
    y_edge[0] = log_T[0] - 0.5 * abs(d_lT[0]) if n_f > 1 else log_T[0] - 0.5
    y_edge[1:-1] = log_T[:-1] + 0.5 * d_lT if n_f > 1 else np.array([])
    y_edge[-1] = (
        log_T[-1] + 0.5 * abs(d_lT[-1]) if n_f > 1 else log_T[-1] + 0.5
    )
    x_edge = (
        np.concatenate(
            [
                [x_cents[0] - 0.5],
                0.5 * (x_cents[:-1] + x_cents[1:]),
                [x_cents[-1] + 0.5],
            ]
        )
        if n_st > 1
        else np.array([-0.5, 0.5])
    )

    X, Y = np.meshgrid(x_edge, y_edge)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    finite = mat[np.isfinite(mat) & (mat > 0)]
    _vmin = vmin if vmin is not None else (
        float(finite.min()) if finite.size else 1e-3
    )
    _vmax = vmax if vmax is not None else (
        float(finite.max()) if finite.size else 1.0
    )
    norm = LogNorm(vmin=max(_vmin, 1e-6), vmax=max(_vmax, _vmin * 1.01))

    qm = ax.pcolormesh(
        X, Y, mat, cmap=cmap, norm=norm, shading="flat", rasterized=True
    )
    if period_up:
        ax.invert_yaxis()

    cb = plt.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label(colorbar_label, fontsize=8)
    cb.ax.tick_params(labelsize=7)

    pos = np.linspace(log_T.min(), log_T.max(), n_yticks)
    labs = []
    for v in pos:
        r = round(v)
        labs.append(
            f"$10^{{{r}}}$" if abs(r - v) < 0.04 else f"$10^{{{v:.1f}}}$"
        )
    ax.set_yticks(pos)
    ax.set_yticklabels(labs, fontsize=tick_fontsize)
    ax.set_ylabel(ylabel, fontsize=8)

    if station_markers:
        _apply_section_station_axis(
            ax, x_cents, st_order, (x_edge[0], x_edge[-1])
        )
    else:
        tick_cfg.apply(
            ax, x_cents, st_order, xlabel=xlabel, xlim=(x_edge[0], x_edge[-1])
        )

    if title:
        ax.set_title(
            title,
            fontsize=9,
            fontweight="bold",
            pad=22.0 if station_markers else None,
        )

    return ax


@EMStyle()
def plot_uncertainty_validation(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    point_size: float = 10.0,
    alpha: float = 0.5,
    show_metrics: bool = True,
    xlabel: str = "True fractional error",
    ylabel: str = "Calibrated fractional error",
    title: str = "",
    figsize: tuple[float, float] = (5.0, 5.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Predicted-vs-true scatter for
    :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`.

    Same held-out-cell validation idea as
    :func:`plot_imputer_validation`, but for a single unit (fractional
    error) rather than several feature components, so one shared axis
    is appropriate here.

    Parameters
    ----------
    y_true, y_pred : ndarray, shape (n_points,)
        Held-out field ``z_err_frac`` and the calibrator's
        reconstruction of it.
    point_size, alpha : float
    show_metrics : bool, default ``True``
    xlabel, ylabel, title : str
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    y_true = np.asarray(y_true, dtype=float).ravel()
    y_pred = np.asarray(y_pred, dtype=float).ravel()
    finite = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true, y_pred = y_true[finite], y_pred[finite]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _validation_scatter_panel(
        ax,
        y_true,
        y_pred,
        point_size=point_size,
        alpha=alpha,
        show_metrics=show_metrics,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
    )
    return ax


@EMStyle()
def plot_uncertainty_summary(
    table: pd.DataFrame,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    *,
    history: Any = None,
    station_labels: list[str] | None = None,
    suptitle: str = "",
    figsize: tuple[float, float] = (13.0, 10.0),
) -> plt.Figure:
    """
    Combined
    :class:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator`
    figure: the calibrated-error map, a per-station original-vs.-
    calibrated comparison, the held-out validation scatter, and (when
    available) the training curve.

    Layout::

        ┌──────────────────────────────────────────────┐
        │ (a) Station x frequency calibrated-error map  │
        ├───────────────────────┬────────────────────────┤
        │ (b) Original vs.      │ (c) Held-out            │
        │     calibrated,       │     validation          │
        │     per station       │                         │
        ├───────────────────────┴────────────────────────┤
        │ (d) Training history                            │
        └──────────────────────────────────────────────┘

    Parameters
    ----------
    table : DataFrame
        Output of
        :meth:`~pycsamt.ai.processing.uncertainty.UncertaintyCalibrator.predict_table`,
        with both ``z_err_frac`` and ``z_err_frac_calibrated`` columns
        -- passed to :func:`plot_uncertainty_map` for panel (a) and
        aggregated per station (median) for panel (b).
    y_true, y_pred : ndarray
        Passed to :func:`plot_uncertainty_validation` for panel (c).
    history : estimator, dict, or None
        Passed to :func:`plot_training_history` when not empty; panel
        (d) is left blank otherwise.
    station_labels : list of str or None
    suptitle : str
    figsize : (w, h), default ``(13.0, 10.0)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    hist = None
    if history is not None:
        h = history.history_ if hasattr(history, "history_") else history
        if isinstance(h, dict) and (
            h.get("train_loss") or h.get("val_loss")
        ):
            hist = h

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[1.1, 1.0, 0.85])

    ax_map = fig.add_subplot(gs[0, :])
    plot_uncertainty_map.__wrapped__(
        table,
        station_labels=station_labels,
        title="(a) Calibrated fractional error",
        ax=ax_map,
    )

    st_order = (
        station_labels
        if station_labels is not None
        else sorted(
            table["station"].unique(), key=lambda x: (str(x).isdigit(), x)
        )
    )
    agg = (
        table.groupby("station")[["z_err_frac", "z_err_frac_calibrated"]]
        .median()
        .reindex(st_order)
    )
    ymax = float(
        np.nanmax(
            [agg["z_err_frac"].max(), agg["z_err_frac_calibrated"].max()]
        )
    )

    ax_bar = fig.add_subplot(gs[1, 0])
    plot_qc_scores.__wrapped__(
        {
            "Original": agg["z_err_frac"].to_numpy(),
            "Calibrated": agg["z_err_frac_calibrated"].to_numpy(),
        },
        station_labels=st_order,
        show_scatter=False,
        show_zone_labels=False,
        show_profile_labels=False,
        bad_zone_alpha=0.0,
        threshold_lw=0.0,
        ylim=(0.0, ymax * 1.2),
        ylabel="Median fractional error",
        title="(b) Original vs. calibrated",
        ax=ax_bar,
    )

    ax_val = fig.add_subplot(gs[1, 1])
    plot_uncertainty_validation.__wrapped__(
        y_true, y_pred, title="(c) Held-out validation", ax=ax_val
    )

    if hist is not None:
        ax_hist = fig.add_subplot(gs[2, :])
        plot_training_history.__wrapped__(
            hist, title="(d) Training history", ax=ax_hist
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# DistortionTypeClassifier plots
# ─────────────────────────────────────────────────────────────────────────────

_DISTORTION_COLORS = ("#2ca02c", "#ff7f0e", "#d62728")  # clean, SS, distorted
_DISTORTION_LABELS_PRETTY = ("Clean", "Static-shift-only", "Distorted")


@EMStyle()
def plot_distortion_map(
    table: pd.DataFrame,
    *,
    station_labels: list[str] | None = None,
    class_colors: tuple[str, str, str] = _DISTORTION_COLORS,
    class_labels: tuple[str, str, str] = _DISTORTION_LABELS_PRETTY,
    bar_width: float = 0.8,
    bar_alpha: float = 0.88,
    tick_every: int | str = "auto",
    tick_label_rotation: float = 45.0,
    tick_fontsize: int = 7,
    station_tick_config: StationTickConfig | None = None,
    xlabel: str = "Station",
    ylabel: str = "Confidence",
    title: str = "",
    show_grid: bool = True,
    figsize: tuple[float, float] = (10.0, 4.0),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Per-station distortion-regime bar chart, for
    :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`.

    Bar colour encodes the predicted regime (clean / static-shift-only
    / distorted); bar height is the model's confidence for that
    station -- both come straight out of
    :meth:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier.predict_table`.

    Parameters
    ----------
    table : DataFrame
        Output of
        :meth:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier.predict_table`.
        Must contain ``station``, ``regime``, and ``confidence``
        columns.
    station_labels : list of str or None
        Station order; defaults to the table's own row order.
    class_colors, class_labels : (str, str, str)
        Colour and legend label for regimes 0, 1, 2.
    bar_width, bar_alpha : float
    tick_every, tick_label_rotation, tick_fontsize, station_tick_config
        Station-axis tick control (see
        :class:`~pycsamt.ai.plot.StationTickConfig`).
    xlabel, ylabel, title : str
    show_grid : bool, default ``True``
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError(
            "plot_distortion_map expects a DataFrame with columns "
            "'station', 'regime', 'confidence'."
        )
    if not {"station", "regime", "confidence"}.issubset(table.columns):
        raise ValueError(
            "DataFrame must have columns 'station', 'regime', "
            "'confidence'."
        )

    if station_labels is not None:
        order = station_labels
        sub = table.set_index("station").reindex(order)
    else:
        sub = table
        order = sub["station"].tolist()

    tick_cfg = _make_tick_config(
        tick_every, tick_label_rotation, tick_fontsize, station_tick_config
    )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    xpos = np.arange(len(order), dtype=float)
    regimes = sub["regime"].to_numpy()
    conf = sub["confidence"].to_numpy()
    colors = [class_colors[int(r)] for r in regimes]

    ax.bar(xpos, conf, color=colors, width=bar_width, alpha=bar_alpha)

    handles = [
        mpatches.Patch(fc=class_colors[i], label=class_labels[i], alpha=0.85)
        for i in range(3)
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7.5, ncol=3)

    tick_cfg.apply(
        ax, xpos, [str(s) for s in order], xlabel=xlabel,
        xlim=(-0.8, len(order) - 0.2),
    )
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_ylim(0.0, 1.05)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    if show_grid:
        ax.grid(True, axis="y", alpha=0.3)

    return ax


@EMStyle()
def plot_distortion_feature_space(
    table: pd.DataFrame,
    *,
    x_col: str = "delta_log10_rho",
    y_col: str = "twist_deg",
    shift_th: float | None = 0.1,
    twist_th: float | None = 10.0,
    class_colors: tuple[str, str, str] = _DISTORTION_COLORS,
    class_labels: tuple[str, str, str] = _DISTORTION_LABELS_PRETTY,
    point_size: float = 45.0,
    alpha: float = 0.85,
    xlabel: str | None = None,
    ylabel: str | None = None,
    title: str = "",
    show_legend: bool = True,
    figsize: tuple[float, float] = (6.0, 5.5),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Two-feature scatter coloured by predicted distortion regime, for
    :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`.

    Shows *why* a station lands in a given class: with the default
    axes, ``delta_log10_rho`` (resistivity departure from the along-
    line spatial trend) drives clean/static-shift-only, and
    ``twist_deg`` (Groom-Bailey rotation) drives the distorted class,
    with the rule-based self-training thresholds drawn as dashed
    guide lines.

    Parameters
    ----------
    table : DataFrame
        Output of
        :meth:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier.predict_table`.
        Must contain ``regime`` plus *x_col* and *y_col*.
    x_col, y_col : str, default ``"delta_log10_rho"``, ``"twist_deg"``
    shift_th, twist_th : float or None
        Rule-based threshold guide lines drawn on the *x* / *y* axis
        respectively (mirrored for negative values where relevant).
        ``None`` skips that guide line -- pass ``None`` for both when
        *x_col* / *y_col* are not the threshold features.
    class_colors, class_labels : (str, str, str)
    point_size, alpha : float
    xlabel, ylabel : str or None
        Default to *x_col* / *y_col*.
    title : str
    show_legend : bool, default ``True``
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if not isinstance(table, pd.DataFrame):
        raise TypeError(
            "plot_distortion_feature_space expects a DataFrame with "
            "columns 'regime', x_col, y_col."
        )
    if not {"regime", x_col, y_col}.issubset(table.columns):
        raise ValueError(
            f"DataFrame must have columns 'regime', {x_col!r}, {y_col!r}."
        )

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    x = table[x_col].to_numpy(dtype=float)
    y = table[y_col].to_numpy(dtype=float)
    regimes = table["regime"].to_numpy()

    for i in range(3):
        m = regimes == i
        if m.any():
            ax.scatter(
                x[m], y[m], s=point_size, alpha=alpha,
                color=class_colors[i], edgecolors="white", linewidths=0.6,
                label=class_labels[i], zorder=3,
            )

    if shift_th is not None:
        ax.axvline(shift_th, color="#888888", lw=1.0, ls="--", zorder=1)
        ax.axvline(-shift_th, color="#888888", lw=1.0, ls="--", zorder=1)
    if twist_th is not None:
        ax.axhline(twist_th, color="#888888", lw=1.0, ls="--", zorder=1)
        ax.axhline(-twist_th, color="#888888", lw=1.0, ls="--", zorder=1)

    ax.set_xlabel(xlabel or x_col, fontsize=9)
    ax.set_ylabel(ylabel or y_col, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    if show_legend:
        ax.legend(fontsize=7.5, framealpha=0.9)

    return ax


@EMStyle()
def plot_distortion_summary(
    table: pd.DataFrame,
    *,
    history: Any = None,
    station_labels: list[str] | None = None,
    x_col: str = "delta_log10_rho",
    y_col: str = "twist_deg",
    shift_th: float = 0.1,
    twist_th: float = 10.0,
    suptitle: str = "",
    figsize: tuple[float, float] = (12.0, 10.5),
) -> plt.Figure:
    """
    Combined
    :class:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier`
    figure: the per-station regime map, the feature-space scatter that
    explains it, class proportions, and (when available) the training
    curve.

    Layout::

        ┌──────────────────────────────────────────────┐
        │ (a) Per-station regime map                    │
        ├───────────────────────┬────────────────────────┤
        │ (b) Feature space     │ (c) Class proportions  │
        ├───────────────────────┴────────────────────────┤
        │ (d) Training history                            │
        └──────────────────────────────────────────────┘

    Parameters
    ----------
    table : DataFrame
        Output of
        :meth:`~pycsamt.ai.processing.distortion.DistortionTypeClassifier.predict_table`.
    history : estimator, dict, or None
        Passed to :func:`plot_training_history` when not empty; panel
        (d) is left blank otherwise.
    station_labels : list of str or None
    x_col, y_col, shift_th, twist_th
        Passed to :func:`plot_distortion_feature_space`.
    suptitle : str
    figsize : (w, h), default ``(12.0, 10.5)``

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    hist = None
    if history is not None:
        h = history.history_ if hasattr(history, "history_") else history
        if isinstance(h, dict) and (
            h.get("train_loss") or h.get("val_loss")
        ):
            hist = h

    fig = plt.figure(figsize=figsize, layout="constrained")
    gs = fig.add_gridspec(3, 2, height_ratios=[0.9, 1.1, 0.85])

    ax_map = fig.add_subplot(gs[0, :])
    plot_distortion_map.__wrapped__(
        table, station_labels=station_labels,
        title="(a) Per-station regime", ax=ax_map,
    )

    ax_feat = fig.add_subplot(gs[1, 0])
    plot_distortion_feature_space.__wrapped__(
        table, x_col=x_col, y_col=y_col, shift_th=shift_th,
        twist_th=twist_th, title="(b) Feature space", ax=ax_feat,
    )

    ax_bar = fig.add_subplot(gs[1, 1])
    counts = table["regime"].value_counts().reindex([0, 1, 2], fill_value=0)
    frac = 100.0 * counts / max(int(counts.sum()), 1)
    ax_bar.bar(
        list(_DISTORTION_LABELS_PRETTY), frac.to_numpy(),
        color=_DISTORTION_COLORS, alpha=0.85,
    )
    ax_bar.set_ylabel("Share of stations (%)", fontsize=8)
    ax_bar.set_title("(c) Class proportions", fontsize=9, fontweight="bold")
    for x, v in enumerate(frac.to_numpy()):
        ax_bar.text(x, v + 1.0, f"{v:.1f}%", ha="center", fontsize=7.5)
    ax_bar.set_ylim(0, max(float(frac.max()) * 1.2, 10.0))
    ax_bar.tick_params(axis="x", labelrotation=15, labelsize=8)

    if hist is not None:
        ax_hist = fig.add_subplot(gs[2, :])
        plot_training_history.__wrapped__(
            hist, title="(d) Training history", ax=ax_hist
        )

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")

    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Time-series denoising
# (:class:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser`)
# ─────────────────────────────────────────────────────────────────────────────

_TSDENOISE_NOISY_SPAN = "#f4c7c3"
_TSDENOISE_RAW_COLOR = "#7f7f7f"
_TSDENOISE_DENOISED_COLOR = "#d62728"
_TSDENOISE_LOW_COLOR = "#1f77b4"


def _shade_noisy_windows(
    ax: plt.Axes, diagnostics: pd.DataFrame, *, color: str = _TSDENOISE_NOISY_SPAN
) -> None:
    if diagnostics is None or diagnostics.empty:
        return
    noisy = diagnostics[diagnostics["label"] == "noisy"]
    for _, row in noisy.iterrows():
        ax.axvspan(
            row["win_start"], row["win_stop"], color=color, alpha=0.6, lw=0,
            zorder=0,
        )


@EMStyle()
def plot_ts_denoise_mmf_split(
    t: np.ndarray,
    x: np.ndarray,
    low: np.ndarray,
    high: np.ndarray,
    *,
    raw_color: str = _TSDENOISE_RAW_COLOR,
    low_color: str = _TSDENOISE_LOW_COLOR,
    xlabel: str = "Time (s)",
    suptitle: str = "",
    figsize: tuple[float, float] = (9.0, 4.6),
) -> plt.Figure:
    """
    Mathematical-morphological-filter (MMF) low/high split.

    Mirrors Fig. 8(a)-(d) of Gui et al. (2024): the raw channel with
    its extracted low-frequency envelope overlaid, and the
    high-frequency residual carrying any strong interference plus
    genuine high-frequency signal.

    Parameters
    ----------
    t : ndarray, shape (n,)
        Time axis (s).
    x, low, high : ndarray, shape (n,)
        Raw channel and :func:`~pycsamt.ai.processing.tsdenoise.mmf_split`
        output (``low + high == x``).
    raw_color, low_color : str
    xlabel : str
    suptitle : str
    figsize : (w, h)

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`

    See Also
    --------
    pycsamt.ai.processing.tsdenoise.mmf_split
        Producer of *low* / *high*.
    """
    fig, (ax0, ax1) = plt.subplots(
        2, 1, figsize=figsize, sharex=True, layout="constrained",
    )
    ax0.plot(t, x, color=raw_color, lw=0.7, alpha=0.85, label="Raw")
    ax0.plot(t, low, color=low_color, lw=1.4, label="MMF low-frequency")
    ax0.set_ylabel("Amplitude", fontsize=9)
    ax0.set_title("(a) Raw + low-frequency envelope", fontsize=9.5,
                   fontweight="bold")
    ax0.legend(fontsize=7.5, frameon=False, loc="upper right")

    ax1.plot(t, high, color=raw_color, lw=0.7)
    ax1.axhline(0, color="0.3", lw=0.6)
    ax1.set_ylabel("Amplitude", fontsize=9)
    ax1.set_xlabel(xlabel, fontsize=9)
    ax1.set_title("(b) High-frequency residual", fontsize=9.5,
                   fontweight="bold")

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")
    return fig


@EMStyle()
def plot_ts_denoise_segments(
    t: np.ndarray,
    high: np.ndarray,
    diagnostics: pd.DataFrame,
    *,
    color: str = _TSDENOISE_RAW_COLOR,
    noisy_span_color: str = _TSDENOISE_NOISY_SPAN,
    xlabel: str = "Time (s)",
    ylabel: str = "Amplitude",
    title: str = "",
    figsize: tuple[float, float] = (9.5, 3.2),
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    High-frequency residual with SVM window classification overlaid.

    Shades every window
    :class:`~pycsamt.ai.processing.tsdenoise.SignalQualityClassifier`
    labelled noisy -- mirrors Fig. 8(e)-(f)'s SVM-Good / SVM-Noisy
    split, drawn as one trace with shaded spans rather than two
    separate panels.

    Parameters
    ----------
    t : ndarray, shape (n,)
    high : ndarray, shape (n,)
        MMF high-frequency residual (see
        :func:`~pycsamt.ai.processing.tsdenoise.mmf_split`).
    diagnostics : DataFrame
        :attr:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser.\
diagnostics_` (or a subset filtered to one channel) -- columns
        ``win_start``, ``win_stop`` (s), ``label``.
    color, noisy_span_color : str
    xlabel, ylabel, title : str
    figsize : (w, h)
    ax : Axes or None

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)

    _shade_noisy_windows(ax, diagnostics, color=noisy_span_color)
    ax.plot(t, high, color=color, lw=0.7)
    ax.axhline(0, color="0.3", lw=0.6)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")

    handles = [
        mpatches.Patch(color=noisy_span_color, alpha=0.6, label="SVM: noisy"),
    ]
    ax.legend(handles=handles, fontsize=7.5, frameon=False,
              loc="upper right")
    return ax


@EMStyle()
def plot_ts_denoise_summary(
    t: np.ndarray,
    x: np.ndarray,
    denoised: np.ndarray,
    diagnostics: pd.DataFrame,
    *,
    low: np.ndarray | None = None,
    high: np.ndarray | None = None,
    raw_color: str = _TSDENOISE_RAW_COLOR,
    denoised_color: str = _TSDENOISE_DENOISED_COLOR,
    noisy_span_color: str = _TSDENOISE_NOISY_SPAN,
    xlabel: str = "Time (s)",
    suptitle: str = "",
    figsize: tuple[float, float] = (10.0, 7.0),
) -> plt.Figure:
    """
    Combined before/after summary for
    :class:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser`.

    Three panels: raw channel with SVM-flagged noisy windows shaded,
    the same windows on the MMF high-frequency residual, and a final
    raw-vs-denoised overlay -- one figure spanning the paper's
    Fig. 8(b), (e)-(f), and (g).

    Parameters
    ----------
    t : ndarray, shape (n,)
    x, denoised : ndarray, shape (n,)
        Channel before and after
        :meth:`~pycsamt.ai.processing.tsdenoise.TimeSeriesDenoiser.apply`
        / :meth:`~...TimeSeriesDenoiser.transform`.
    diagnostics : DataFrame
        See :func:`plot_ts_denoise_segments`.
    low, high : ndarray or None
        :func:`~pycsamt.ai.processing.tsdenoise.mmf_split` output for
        *x*; when given, the second panel plots *high* (the residual
        actually classified) instead of *x*.
    raw_color, denoised_color, noisy_span_color : str
    xlabel : str
    suptitle : str
    figsize : (w, h)

    Returns
    -------
    fig : :class:`matplotlib.figure.Figure`
    """
    fig, (ax0, ax1, ax2) = plt.subplots(
        3, 1, figsize=figsize, sharex=True, layout="constrained",
    )

    _shade_noisy_windows(ax0, diagnostics, color=noisy_span_color)
    ax0.plot(t, x, color=raw_color, lw=0.7)
    ax0.set_ylabel("Amplitude", fontsize=9)
    ax0.set_title("(a) Raw channel -- SVM-flagged windows shaded",
                   fontsize=9.5, fontweight="bold")

    plot_ts_denoise_segments.__wrapped__(
        t, high if high is not None else x, diagnostics,
        noisy_span_color=noisy_span_color,
        title="(b) High-frequency residual", ax=ax1,
    )

    ax2.plot(t, x, color=raw_color, lw=0.7, alpha=0.7, label="Raw")
    ax2.plot(t, denoised, color=denoised_color, lw=1.0, label="Denoised")
    ax2.set_xlabel(xlabel, fontsize=9)
    ax2.set_ylabel("Amplitude", fontsize=9)
    ax2.set_title("(c) Raw vs. denoised", fontsize=9.5, fontweight="bold")
    ax2.legend(fontsize=7.5, frameon=False, loc="upper right")

    if suptitle:
        fig.suptitle(suptitle, fontsize=11, fontweight="bold")
    return fig
