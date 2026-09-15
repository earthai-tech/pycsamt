"""Quality-control confidence ratios for EM transfer functions.

The composite confidence ratio (CR) used by this module is a bounded,
weighted score:

    CR = sum_k w_k s_k / sum_k w_k,  for finite component scores s_k.

The default components are data coverage, tensor uncertainty,
off-diagonal consistency, diagonal leakage, phase smoothness, and spatial
coherence. Each score is clipped to [0, 1], where 1 is most trustworthy.
The default manuscript classes are CR >= 0.95 (safe), 0.85 <= CR < 0.95
(recoverable/marginal), and CR < 0.85 (reject/review).
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle as _Rect

from ..api.labels import LOG10_PERIOD_LABEL
from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.station import (
    PYCSAMT_STATION_RENDERING,
    StationAxisStyle,
)
from ..api.view import maybe_wrap_frame
from ..compat.numpy import trapz as _trapz
from ._core import (
    _axes_list,
    _get_t_block,
    _get_z_block,
    _iter_items,
    _name,
    _station_positions,
    ensure_sites,
)
from .tensor import build_phase_tensor_table

__all__ = [
    "build_qc_table",
    "confidence_ratio",
    "frequency_confidence_table",
    "export_confidence_map",
    "plot_confidence_band_summary",
    "plot_confidence_before_after",
    "plot_confidence_component_map",
    "plot_confidence_coverage_curve",
    "plot_confidence_distribution",
    "plot_confidence_heatmap",
    "plot_confidence_map",
    "plot_confidence_grid_map",
    "plot_confidence_method_comparison",
    "plot_confidence_profile",
    "plot_confidence_rank",
    "plot_confidence_risk_map",
    "plot_frequency_confidence_psection",
    "plot_station_confidence_dashboard",
    "plot_station_confidence_spectrum",
    "qc_flags",
    "station_confidence_table",
]

DEFAULT_CONFIDENCE_WEIGHTS: dict[str, float] = {
    "coverage": 0.35,
    "uncertainty": 0.20,
    "offdiag": 0.15,
    "diagonal": 0.10,
    "phase": 0.10,
    "spatial": 0.10,
}

DEFAULT_CI_HI = 0.95
DEFAULT_CI_LO = 0.85


# ------------------------------ helpers --------------------------------- #


def _station_map_coordinate(ed: Any) -> tuple[float, float, float, float]:
    """Return ``(longitude, latitude, easting, northing)`` for one site."""
    lon = lat = east = north = np.nan
    candidates = (ed, getattr(ed, "edi", None))

    def _first(names: tuple[str, ...]) -> float:
        for candidate in candidates:
            if candidate is None:
                continue
            for name in names:
                try:
                    value = float(getattr(candidate, name))
                except (AttributeError, TypeError, ValueError):
                    continue
                if np.isfinite(value):
                    return value
        return np.nan

    lon = _first(("longitude", "lon", "long"))
    lat = _first(("latitude", "lat"))
    east = _first(("east", "easting", "x"))
    north = _first(("north", "northing", "y"))
    if not (np.isfinite(lon) and np.isfinite(lat)):
        try:
            from ..site.utils import get_coords

            coords = get_coords(getattr(ed, "edi", ed))
            lat, lon = float(coords.lat), float(coords.lon)
        except Exception:
            coords = getattr(ed, "coords", None)
            try:
                lat, lon = float(coords[0]), float(coords[1])
            except (TypeError, ValueError, IndexError):
                pass
    return lon, lat, east, north


def _order_map_route(
    group: pd.DataFrame, xkey: str, ykey: str
) -> pd.DataFrame:
    """Order one survey line by coordinate-derived principal chainage."""
    if len(group) < 3:
        return group
    xy = group[[xkey, ykey]].to_numpy(dtype=float)
    centred = xy - np.nanmean(xy, axis=0)
    if xkey == "longitude":
        centred[:, 0] *= np.cos(np.deg2rad(np.nanmean(group[ykey])))
    try:
        _, _, vh = np.linalg.svd(centred, full_matrices=False)
    except np.linalg.LinAlgError:
        return group
    chainage = centred @ vh[0]
    return group.iloc[np.argsort(chainage, kind="stable")]


def _resolve_section_style(section: str | SectionStyle) -> SectionStyle:
    """Return a copied section style for EMTools pseudo-sections."""
    if isinstance(section, SectionStyle):
        return section.copy()
    return PYCSAMT_SECTION.style_for(str(section)).copy()


def _row_ok_z(z: np.ndarray) -> np.ndarray:
    y = z.reshape(z.shape[0], -1)
    return np.isfinite(y).all(axis=1)


def _row_ok_t(t: np.ndarray) -> np.ndarray:
    y = t.reshape(t.shape[0], -1)
    return np.isfinite(y).all(axis=1)


def _row_nanmedian(values: np.ndarray) -> np.ndarray:
    """Return row medians without warning for all-NaN rows."""
    out = np.full(values.shape[0], np.nan, dtype=float)
    valid_rows = np.isfinite(values).any(axis=1)
    if valid_rows.any():
        out[valid_rows] = np.nanmedian(values[valid_rows], axis=1)
    return out


def _snr_rows(z: np.ndarray, ze: np.ndarray | None) -> np.ndarray:
    """Return per-row SNR without warning for all-NaN rows.

    Same "no warning" convention as :func:`_row_nanmedian`: a frequency
    row where every tensor component is NaN (a real per-frequency data
    gap, not a bug) yields NaN silently instead of a "Mean of empty
    slice" RuntimeWarning.
    """
    n = z.shape[0]
    if ze is None:
        return np.full(n, np.nan, dtype=float)
    z2 = np.abs(z) ** 2
    ze2 = np.abs(ze) ** 2
    valid_z = np.isfinite(z2).any(axis=(1, 2))
    valid_e = np.isfinite(ze2).any(axis=(1, 2))
    a = np.full(n, np.nan, dtype=float)
    e = np.full(n, np.nan, dtype=float)
    if valid_z.any():
        a[valid_z] = np.sqrt(np.nanmean(z2[valid_z], axis=(1, 2)))
    if valid_e.any():
        e[valid_e] = np.sqrt(np.nanmean(ze2[valid_e], axis=(1, 2)))
    return a / (e + 1e-12)


def _offdiag_logmag(z: np.ndarray) -> np.ndarray:
    m = _row_nanmedian(
        np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
    )
    return np.log10(np.maximum(m, 1e-24))


def _clip01(x: Any) -> float:
    """Return finite scalar clipped to the confidence interval."""
    try:
        value = float(x)
    except (TypeError, ValueError):
        return np.nan
    if not np.isfinite(value):
        return np.nan
    return float(np.clip(value, 0.0, 1.0))


def _weighted_nanmean(
    values: dict[str, float], weights: dict[str, float]
) -> float:
    """Return weighted mean ignoring unavailable metrics."""
    total = 0.0
    weight = 0.0
    for key, value in values.items():
        value = _clip01(value)
        w = float(weights.get(key, 0.0))
        if np.isfinite(value) and w > 0.0:
            total += w * value
            weight += w
    return float(total / weight) if weight > 0.0 else np.nan


def _confidence_error(
    values: dict[str, float], n_freq: int, confidence: float
) -> float:
    """Estimate a compact station-level confidence uncertainty."""
    vals = np.asarray(
        [_clip01(value) for value in values.values()],
        dtype=float,
    )
    vals = vals[np.isfinite(vals)]
    if vals.size > 1:
        return float(np.nanstd(vals, ddof=0))
    confidence = _clip01(confidence)
    if not np.isfinite(confidence):
        return np.nan
    n_freq = max(1, int(n_freq))
    return float(np.sqrt(confidence * (1.0 - confidence) / n_freq))


def confidence_ratio(
    scores: dict[str, float],
    *,
    weights: dict[str, float] | None = None,
    n_freq: int = 1,
    return_error: bool = False,
) -> float | tuple[float, float]:
    r"""Compute the composite confidence ratio from diagnostic scores.

    The confidence ratio is a weighted finite-score mean:

    .. math::

        \mathrm{CR} =
        \frac{\sum_k w_k s_k \mathbf{1}_{s_k\ finite}}
             {\sum_k w_k \mathbf{1}_{s_k\ finite}},
        \qquad 0 \leq s_k \leq 1.

    The default score vector is
    ``coverage, uncertainty, offdiag, diagonal, phase, spatial`` with
    weights ``0.35, 0.20, 0.15, 0.10, 0.10, 0.10``. Missing scores are
    ignored and all finite scores are clipped to ``[0, 1]``.

    The optional error is the population spread of available component
    scores; when only one score is available it falls back to the binomial
    standard error ``sqrt(CR * (1 - CR) / n_freq)``.
    """
    use_weights = {**DEFAULT_CONFIDENCE_WEIGHTS, **(weights or {})}
    cr = _weighted_nanmean(scores, use_weights)
    if return_error:
        return cr, _confidence_error(scores, n_freq, cr)
    return cr


def _relerr_score(
    z: np.ndarray, ze: np.ndarray | None, threshold: float
) -> float:
    """Score tensor uncertainty from median relative error."""
    if ze is None:
        return np.nan
    rel = np.abs(ze) / (np.abs(z) + 1e-24)
    if not np.isfinite(rel).any():
        return np.nan
    med = float(np.nanmedian(rel))
    return _clip01(1.0 - med / max(float(threshold), 1e-12))


def _offdiag_consistency_score(z: np.ndarray, tolerance_log10: float) -> float:
    """Score similarity of ``Zxy`` and ``Zyx`` amplitudes."""
    zxy = np.abs(z[:, 0, 1])
    zyx = np.abs(z[:, 1, 0])
    ratio = np.log10((zxy + 1e-24) / (zyx + 1e-24))
    if not np.isfinite(ratio).any():
        return np.nan
    med = float(np.nanmedian(np.abs(ratio)))
    return _clip01(1.0 - med / max(float(tolerance_log10), 1e-12))


def _diagonal_leakage_score(z: np.ndarray, max_fraction: float) -> float:
    """Score how much diagonal impedance leaks into off-diagonal terms."""
    diag = _row_nanmedian(
        np.stack([np.abs(z[:, 0, 0]), np.abs(z[:, 1, 1])], axis=1),
    )
    off = _row_nanmedian(
        np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
    )
    frac = diag / (off + diag + 1e-24)
    if not np.isfinite(frac).any():
        return np.nan
    med = float(np.nanmedian(frac))
    return _clip01(1.0 - med / max(float(max_fraction), 1e-12))


def _phase_smoothness_score(z: np.ndarray, jump_tolerance_deg: float) -> float:
    """Score abrupt phase jumps in the off-diagonal components."""
    phases = []
    for comp in (z[:, 0, 1], z[:, 1, 0]):
        if not np.isfinite(comp).any():
            continue
        ph = np.unwrap(np.angle(comp))
        if ph.size > 1:
            jumps = np.rad2deg(np.abs(np.diff(ph)))
            if np.isfinite(jumps).any():
                phases.append(jumps)
    if not phases:
        return np.nan
    phase_jumps = np.concatenate(phases)
    if not np.isfinite(phase_jumps).any():
        return np.nan
    med_jump = float(np.nanmedian(phase_jumps))
    return _clip01(1.0 - med_jump / max(float(jump_tolerance_deg), 1e-12))


def _station_spatial_scores(
    med_logrho: np.ndarray,
    tolerance_log10: float,
) -> np.ndarray:
    """Score station coherence against immediate neighboring stations."""
    scores = np.full(med_logrho.size, np.nan, dtype=float)
    for i, value in enumerate(med_logrho):
        neighbors = []
        if i > 0 and np.isfinite(med_logrho[i - 1]):
            neighbors.append(med_logrho[i - 1])
        if i + 1 < med_logrho.size and np.isfinite(med_logrho[i + 1]):
            neighbors.append(med_logrho[i + 1])
        if not neighbors or not np.isfinite(value):
            continue
        ref = float(np.nanmedian(neighbors))
        diff = abs(float(value) - ref)
        scores[i] = _clip01(1.0 - diff / max(float(tolerance_log10), 1e-12))
    return scores


def _frequency_spatial_scores(
    table: pd.DataFrame,
    tolerance_log10: float,
) -> np.ndarray:
    """Score frequency samples against same-frequency neighbor stations."""
    scores = np.full(len(table), np.nan, dtype=float)
    if table.empty:
        return scores
    for _, group in table.groupby("frequency_hz", sort=False):
        order = group.sort_values("distance_m")
        idx = order.index.to_numpy(dtype=int)
        values = order["logrho_proxy"].to_numpy(dtype=float)
        for j, row_index in enumerate(idx):
            neighbors = []
            if j > 0 and np.isfinite(values[j - 1]):
                neighbors.append(values[j - 1])
            if j + 1 < values.size and np.isfinite(values[j + 1]):
                neighbors.append(values[j + 1])
            if not neighbors or not np.isfinite(values[j]):
                continue
            ref = float(np.nanmedian(neighbors))
            diff = abs(float(values[j]) - ref)
            scores[row_index] = _clip01(
                1.0 - diff / max(float(tolerance_log10), 1e-12),
            )
    return scores


def _frequency_phase_jump_score(
    z: np.ndarray,
    jump_tolerance_deg: float,
) -> np.ndarray:
    """Return per-frequency smoothness scores for off-diagonal phase."""
    scores = np.full(z.shape[0], np.nan, dtype=float)
    jumps = []
    for comp in (z[:, 0, 1], z[:, 1, 0]):
        if not np.isfinite(comp).any():
            continue
        phase = np.unwrap(np.angle(comp))
        if phase.size < 2:
            continue
        local = np.full(phase.size, np.nan, dtype=float)
        dphase = np.rad2deg(np.abs(np.diff(phase)))
        local[:-1] = dphase
        prior = local[1:].copy()
        current = dphase
        both = np.isfinite(prior) & np.isfinite(current)
        only_current = ~np.isfinite(prior) & np.isfinite(current)
        prior[both] = np.maximum(prior[both], current[both])
        prior[only_current] = current[only_current]
        local[1:] = prior
        jumps.append(local)
    if not jumps:
        return scores
    jump = _row_nanmedian(np.stack(jumps, axis=1))
    valid = np.isfinite(jump)
    scores[valid] = [
        _clip01(1.0 - value / max(float(jump_tolerance_deg), 1e-12))
        for value in jump[valid]
    ]
    return scores


def _frequency_flags(row: pd.Series, ci_hi: float, ci_lo: float) -> str:
    """Return readable quality flags for one frequency-confidence row."""
    flags = []
    if row["confidence"] < ci_lo:
        flags.append("reject")
    elif row["confidence"] < ci_hi:
        flags.append("recoverable")
    if row["coverage"] < 1.0:
        flags.append("missing")
    if np.isfinite(row["uncertainty"]) and row["uncertainty"] < ci_lo:
        flags.append("high_error")
    if np.isfinite(row["offdiag"]) and row["offdiag"] < ci_lo:
        flags.append("offdiag_mismatch")
    if np.isfinite(row["diagonal"]) and row["diagonal"] < ci_lo:
        flags.append("diagonal_leakage")
    if np.isfinite(row["phase"]) and row["phase"] < ci_lo:
        flags.append("phase_jump")
    if np.isfinite(row["spatial"]) and row["spatial"] < ci_lo:
        flags.append("spatial_outlier")
    return ",".join(flags)


def _y_ticks(yall: np.ndarray, ny: int) -> tuple[np.ndarray, list[str]]:
    yt = np.linspace(0, yall.size - 1, num=min(ny, yall.size))
    yv = np.linspace(yall.min(), yall.max(), num=yt.size)
    lab = [f"{v:.2g}" for v in yv]
    return yt, lab


# ------------------------------ tables ---------------------------------- #


def build_qc_table(
    sites: Any,
    *,
    include_skew: bool = True,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    pt = None
    if include_skew:
        pt = build_phase_tensor_table(
            S,
            recursive=False,
            on_dup=on_dup,
            strict=False,
            verbose=verbose,
        )
    rows: list[dict[str, Any]] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        T, t, ft = _get_t_block(ed)
        if Z is None:
            continue
        n = z.shape[0]
        ko = _row_ok_z(z)
        ze = getattr(Z, "z_err", None)
        snr = _snr_rows(z, ze)
        med_snr = float(np.nanmedian(snr))
        n_ok = int(np.nansum(ko))
        n_t = 0
        n_t_ok = 0
        if T is not None and t is not None:
            n_t = t.shape[0]
            n_t_ok = int(np.nansum(_row_ok_t(t)))
        per = 1.0 / fr
        pmin = float(np.nanmin(per)) if per.size else np.nan
        pmax = float(np.nanmax(per)) if per.size else np.nan
        rec = dict(
            station=st,
            n_freq=int(n),
            n_ok=int(n_ok),
            frac_ok=float(n_ok / max(1, n)),
            n_tip=int(n_t),
            n_tip_ok=int(n_t_ok),
            snr_med=med_snr,
            pmin=pmin,
            pmax=pmax,
        )
        if include_skew and pt is not None:
            sdf = pt[pt["station"] == st]
            if not sdf.empty:
                sb = np.abs(sdf["beta"].to_numpy(dtype=float))
                rec["skew_med"] = float(np.nanmedian(sb))
                rec["skew_iqr"] = float(
                    np.nanpercentile(sb, 75) - np.nanpercentile(sb, 25)
                )
            else:
                rec["skew_med"] = np.nan
                rec["skew_iqr"] = np.nan
        rows.append(rec)
    cols = [
        "station",
        "n_freq",
        "n_ok",
        "frac_ok",
        "n_tip",
        "n_tip_ok",
        "snr_med",
        "pmin",
        "pmax",
    ]
    if include_skew:
        cols += ["skew_med", "skew_iqr"]
    df = pd.DataFrame.from_records(rows, columns=cols)

    return maybe_wrap_frame(
        df,
        api=api,
        name="qc_table",
        kind="emtools.qc",
        source=sites,
        description="Station-level transfer-function quality summary.",
    )


def qc_flags(
    sites: Any,
    *,
    min_frac_ok: float = 0.6,
    min_snr_med: float = 2.0,
    max_skew_med: float = 6.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> pd.DataFrame:
    tb = build_qc_table(
        sites,
        include_skew=True,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if tb.empty:
        return tb
    flags = []
    for _, r in tb.iterrows():
        f = []
        if float(r["frac_ok"]) < float(min_frac_ok):
            f.append("low_coverage")
        if np.isfinite(r["snr_med"]) and r["snr_med"] < min_snr_med:
            f.append("low_snr")
        if (
            np.isfinite(r.get("skew_med", np.nan))
            and r["skew_med"] > max_skew_med
        ):
            f.append("high_skew")
        flags.append(",".join(f))
    out = tb.copy()
    out["flags"] = flags
    return out


def station_confidence_table(
    sites: Any,
    *,
    method: str = "composite",
    weights: dict[str, float] | None = None,
    relerr_threshold: float = 0.20,
    offdiag_tolerance_log10: float = 0.35,
    diagonal_leakage_max: float = 0.35,
    phase_jump_tolerance_deg: float = 90.0,
    spatial_tolerance_log10: float = 0.60,
    spacing_m: float = 200.0,
    force_spacing: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    """Return station-level confidence scores for EM transfer functions.

    ``method="presence"`` reproduces the legacy criterion based only on
    finite tensor rows.  ``method="composite"`` combines several station
    trust indicators: finite data coverage, tensor uncertainty when error
    tensors exist, off-diagonal consistency, diagonal leakage, phase
    smoothness, and spatial coherence with neighboring stations.

    ``distance_m`` in the returned table is the real inter-station
    distance projected along the survey line, derived from EDI
    coordinates (east/north, or lat/lon as a fallback) whenever at least
    two stations carry usable coordinates. ``spacing_m`` is only used as
    a uniform per-station fallback for stations without coordinates, or
    for the whole line when no station has any. Pass
    ``force_spacing=True`` to bypass coordinate lookup entirely and lay
    every station out at uniform ``spacing_m`` steps -- e.g. when the
    available coordinates are known to be unreliable.
    """
    method = str(method).lower()
    if method not in {"presence", "composite"}:
        msg = "method must be 'presence' or 'composite'."
        raise ValueError(msg)
    weights = {**DEFAULT_CONFIDENCE_WEIGHTS, **(weights or {})}
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = list(_iter_items(S))
    positions = _station_positions(
        items, spacing_m, force_spacing=force_spacing
    )
    rows: list[dict[str, Any]] = []
    med_logrho = []
    for i, ed in enumerate(items):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None or z is None or fr is None:
            continue
        ze = getattr(Z, "z_err", None)
        ok = _row_ok_z(z)
        coverage = float(np.nansum(ok) / max(1, z.shape[0]))
        zxy = z[:, 0, 1]
        zyx = z[:, 1, 0]
        rho_proxy = 0.5 * (
            np.log10(np.abs(zxy) ** 2 / np.maximum(fr, 1e-24) + 1e-24)
            + np.log10(np.abs(zyx) ** 2 / np.maximum(fr, 1e-24) + 1e-24)
        )
        med_logrho.append(float(np.nanmedian(rho_proxy)))
        score_parts = {
            "coverage": coverage,
            "uncertainty": _relerr_score(z, ze, relerr_threshold),
            "offdiag": _offdiag_consistency_score(
                z,
                offdiag_tolerance_log10,
            ),
            "diagonal": _diagonal_leakage_score(
                z,
                diagonal_leakage_max,
            ),
            "phase": _phase_smoothness_score(z, phase_jump_tolerance_deg),
        }
        confidence = coverage
        if method == "composite":
            confidence = confidence_ratio(score_parts, weights=weights)
            error_parts = score_parts
        else:
            error_parts = {"coverage": coverage}
        confidence_err = _confidence_error(
            error_parts,
            z.shape[0],
            confidence,
        )
        longitude, latitude, easting, northing = _station_map_coordinate(ed)
        rows.append(
            dict(
                station=st,
                distance_m=float(positions[i])
                if i < positions.size
                else np.nan,
                confidence=float(confidence),
                confidence_err=float(confidence_err),
                method=method,
                n_freq=int(z.shape[0]),
                n_ok=int(np.nansum(ok)),
                coverage=score_parts["coverage"],
                uncertainty=score_parts["uncertainty"],
                offdiag=score_parts["offdiag"],
                diagonal=score_parts["diagonal"],
                phase=score_parts["phase"],
                spatial=np.nan,
                longitude=longitude,
                latitude=latitude,
                easting=easting,
                northing=northing,
            )
        )
    if not rows:
        df = pd.DataFrame(
            columns=[
                "station",
                "distance_m",
                "confidence",
                "method",
                "confidence_err",
                "n_freq",
                "n_ok",
                "coverage",
                "uncertainty",
                "offdiag",
                "diagonal",
                "phase",
                "spatial",
                "longitude",
                "latitude",
                "easting",
                "northing",
            ]
        )

        return maybe_wrap_frame(
            df,
            api=api,
            name="station_confidence_table",
            kind="emtools.qc.station_confidence",
            source=sites,
        )
    spatial_scores = _station_spatial_scores(
        np.asarray(med_logrho, dtype=float),
        spatial_tolerance_log10,
    )
    if method == "composite":
        for i, row in enumerate(rows):
            row["spatial"] = (
                float(spatial_scores[i]) if i < spatial_scores.size else np.nan
            )
            parts = {
                key: row[key]
                for key in (
                    "coverage",
                    "uncertainty",
                    "offdiag",
                    "diagonal",
                    "phase",
                    "spatial",
                )
            }
            row["confidence"] = confidence_ratio(parts, weights=weights)
            row["confidence_err"] = _confidence_error(
                parts,
                row["n_freq"],
                row["confidence"],
            )
    else:
        for i, row in enumerate(rows):
            row["spatial"] = (
                float(spatial_scores[i]) if i < spatial_scores.size else np.nan
            )
    df = pd.DataFrame.from_records(rows)

    return maybe_wrap_frame(
        df,
        api=api,
        name="station_confidence_table",
        kind="emtools.qc.station_confidence",
        source=sites,
        description="Station-level composite confidence scores.",
    )


def frequency_confidence_table(
    sites: Any,
    *,
    method: str = "composite",
    weights: dict[str, float] | None = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    relerr_threshold: float = 0.20,
    offdiag_tolerance_log10: float = 0.35,
    diagonal_leakage_max: float = 0.35,
    phase_jump_tolerance_deg: float = 90.0,
    spatial_tolerance_log10: float = 0.60,
    spacing_m: float = 200.0,
    force_spacing: bool = False,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    api: bool | None = None,
) -> Any:
    """Return frequency-level confidence scores for EM stations.

    The returned table has one row for each station-frequency sample.  It is
    designed as a reusable quality-control source for plots, masking rules,
    and inversion-preparation reports.  ``method="presence"`` scores only
    finite impedance-tensor availability.  ``method="composite"`` combines
    coverage, tensor uncertainty, off-diagonal consistency, diagonal leakage,
    phase smoothness, and same-frequency spatial coherence.

    See :func:`station_confidence_table` for how ``distance_m``,
    ``spacing_m``, and ``force_spacing`` interact.
    """
    method = str(method).lower()
    if method not in {"presence", "composite"}:
        msg = "method must be 'presence' or 'composite'."
        raise ValueError(msg)
    weights = {**DEFAULT_CONFIDENCE_WEIGHTS, **(weights or {})}
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    items = list(_iter_items(S))
    positions = _station_positions(
        items, spacing_m, force_spacing=force_spacing
    )
    rows: list[dict[str, Any]] = []
    for station_index, ed in enumerate(items):
        station = _name(ed, station_index)
        Z, z, fr = _get_z_block(ed)
        if Z is None or z is None or fr is None:
            continue
        ze = getattr(Z, "z_err", None)
        z_abs = np.abs(z)
        coverage = np.isfinite(z.reshape(z.shape[0], -1)).mean(axis=1)
        zxy = z[:, 0, 1]
        zyx = z[:, 1, 0]
        logrho_proxy = 0.5 * (
            np.log10(np.abs(zxy) ** 2 / np.maximum(fr, 1e-24) + 1e-24)
            + np.log10(np.abs(zyx) ** 2 / np.maximum(fr, 1e-24) + 1e-24)
        )
        uncertainty = np.full(z.shape[0], np.nan, dtype=float)
        if ze is not None:
            rel = np.abs(ze) / (z_abs + 1e-24)
            rel_med = _row_nanmedian(rel.reshape(rel.shape[0], -1))
            uncertainty = np.asarray(
                [
                    _clip01(1.0 - value / max(float(relerr_threshold), 1e-12))
                    for value in rel_med
                ],
                dtype=float,
            )
        ratio = np.log10((np.abs(zxy) + 1e-24) / (np.abs(zyx) + 1e-24))
        offdiag = np.asarray(
            [
                _clip01(
                    1.0
                    - abs(value) / max(float(offdiag_tolerance_log10), 1e-12),
                )
                for value in ratio
            ],
            dtype=float,
        )
        diag = _row_nanmedian(
            np.stack([np.abs(z[:, 0, 0]), np.abs(z[:, 1, 1])], axis=1),
        )
        off = _row_nanmedian(
            np.stack([np.abs(z[:, 0, 1]), np.abs(z[:, 1, 0])], axis=1),
        )
        frac = diag / (off + diag + 1e-24)
        diagonal = np.asarray(
            [
                _clip01(
                    1.0 - value / max(float(diagonal_leakage_max), 1e-12),
                )
                for value in frac
            ],
            dtype=float,
        )
        phase = _frequency_phase_jump_score(z, phase_jump_tolerance_deg)
        for freq_index, freq in enumerate(fr):
            parts = {
                "coverage": float(coverage[freq_index]),
                "uncertainty": float(uncertainty[freq_index]),
                "offdiag": float(offdiag[freq_index]),
                "diagonal": float(diagonal[freq_index]),
                "phase": float(phase[freq_index]),
            }
            confidence = parts["coverage"]
            if method == "composite":
                confidence = confidence_ratio(parts, weights=weights)
            error_parts = (
                parts
                if method == "composite"
                else {
                    "coverage": parts["coverage"],
                }
            )
            row = dict(
                station=station,
                station_index=int(station_index),
                distance_m=(
                    float(positions[station_index])
                    if station_index < positions.size
                    else np.nan
                ),
                frequency_hz=float(freq),
                period_s=float(1.0 / freq) if freq else np.nan,
                log10_period=(
                    float(np.log10(1.0 / freq)) if freq > 0 else np.nan
                ),
                confidence=float(confidence),
                confidence_err=_confidence_error(error_parts, 1, confidence),
                method=method,
                n_components=int(np.isfinite(z[freq_index]).sum()),
                coverage=parts["coverage"],
                uncertainty=parts["uncertainty"],
                offdiag=parts["offdiag"],
                diagonal=parts["diagonal"],
                phase=parts["phase"],
                spatial=np.nan,
                logrho_proxy=float(logrho_proxy[freq_index]),
                flags="",
            )
            rows.append(row)
    columns = [
        "station",
        "station_index",
        "distance_m",
        "frequency_hz",
        "period_s",
        "log10_period",
        "confidence",
        "confidence_err",
        "method",
        "n_components",
        "coverage",
        "uncertainty",
        "offdiag",
        "diagonal",
        "phase",
        "spatial",
        "logrho_proxy",
        "flags",
    ]
    if not rows:
        df = pd.DataFrame(columns=columns)

        return maybe_wrap_frame(
            df,
            api=api,
            name="frequency_confidence_table",
            kind="emtools.qc.frequency_confidence",
            source=sites,
        )
    table = pd.DataFrame.from_records(rows, columns=columns)
    spatial = _frequency_spatial_scores(table, spatial_tolerance_log10)
    table["spatial"] = spatial
    if method == "composite":
        for index, row in table.iterrows():
            parts = {
                key: row[key]
                for key in (
                    "coverage",
                    "uncertainty",
                    "offdiag",
                    "diagonal",
                    "phase",
                    "spatial",
                )
            }
            confidence = confidence_ratio(parts, weights=weights)
            table.at[index, "confidence"] = confidence
            table.at[index, "confidence_err"] = _confidence_error(
                parts,
                1,
                confidence,
            )
    table["flags"] = [
        _frequency_flags(row, ci_hi, ci_lo) for _, row in table.iterrows()
    ]

    return maybe_wrap_frame(
        table,
        api=api,
        name="frequency_confidence_table",
        kind="emtools.qc.frequency_confidence",
        source=sites,
        description="Frequency-level transfer-function confidence scores.",
    )


# --------------------------- confidence maps ---------------------------- #


def export_confidence_map(
    sites: Any,
    *,
    csv_path: str | Path | None = None,
    surfer_path: str | Path | None = None,
    method: str = "composite",
    coordinate_system: str = "auto",
    line_labels: Any = None,
    grid_shape: tuple[int, int] = (200, 200),
    max_triangle_edge: float | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> dict[str, Path]:
    """Export station confidence to CSV and/or a Surfer DSAA grid.

    CSV contains station coordinates, confidence, uncertainty, and composite
    components. The Surfer grid contains linearly interpolated confidence
    inside the station convex hull; cells outside it (and optionally across
    triangles longer than ``max_triangle_edge``) use Surfer's blank value.
    ``grid_shape`` is ``(nx, ny)``.
    """
    if csv_path is None and surfer_path is None:
        raise ValueError("provide csv_path, surfer_path, or both.")
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    tb = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    geographic = np.isfinite(tb[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(tb[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    tb.insert(1, "coordinate_system", coordinate_system)
    tb.insert(2, "x", tb[xkey])
    tb.insert(3, "y", tb[ykey])
    if line_labels is not None:
        if isinstance(line_labels, str):
            tb.insert(1, "line", line_labels)
        elif hasattr(line_labels, "get"):
            tb.insert(
                1,
                "line",
                tb["station"].astype(str).map(line_labels),
            )
        else:
            values = list(line_labels)
            if len(values) != len(tb):
                raise ValueError("line_labels must contain one value per row.")
            tb.insert(1, "line", values)

    outputs: dict[str, Path] = {}
    if csv_path is not None:
        out = Path(csv_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        tb.to_csv(out, index=False)
        outputs["csv"] = out

    if surfer_path is not None:
        grid_tb = (
            tb.loc[valid, ["x", "y", "confidence"]]
            .groupby(["x", "y"], as_index=False)["confidence"]
            .mean()
        )
        xy = grid_tb[["x", "y"]].to_numpy(dtype=float)
        if len(xy) < 3 or np.linalg.matrix_rank(xy - xy.mean(axis=0)) < 2:
            raise ValueError(
                "Surfer-grid export requires at least three non-collinear "
                "stations; export CSV for a single survey line."
            )
        from matplotlib.tri import LinearTriInterpolator, Triangulation

        triangulation = Triangulation(xy[:, 0], xy[:, 1])
        if max_triangle_edge is not None:
            limit = float(max_triangle_edge)
            if not np.isfinite(limit) or limit <= 0.0:
                raise ValueError(
                    "max_triangle_edge must be positive and finite."
                )
            vertices = xy[triangulation.triangles]
            edges = np.linalg.norm(
                vertices - np.roll(vertices, -1, axis=1), axis=2
            )
            triangulation.set_mask((edges > limit).any(axis=1))
            if np.all(triangulation.mask):
                raise ValueError(
                    "max_triangle_edge masks every grid triangle."
                )
        nx, ny = (int(grid_shape[0]), int(grid_shape[1]))
        if nx < 2 or ny < 2:
            raise ValueError("grid_shape values must both be at least 2.")
        x_axis = np.linspace(float(xy[:, 0].min()), float(xy[:, 0].max()), nx)
        y_axis = np.linspace(float(xy[:, 1].min()), float(xy[:, 1].max()), ny)
        gx, gy = np.meshgrid(x_axis, y_axis)
        interpolator = LinearTriInterpolator(
            triangulation,
            grid_tb["confidence"].to_numpy(dtype=float),
        )
        grid = np.ma.asarray(interpolator(gx, gy)).filled(np.nan)
        from ..interp.export import to_surfer_dsaa

        outputs["surfer"] = to_surfer_dsaa(
            x_axis,
            y_axis,
            grid,
            surfer_path,
        )
    return outputs


def plot_confidence_map(
    sites: Any,
    *,
    method: str = "composite",
    coordinate_system: str = "auto",
    mode: str = "auto",
    line_labels: Any = None,
    connect: bool = True,
    contour_levels: int | Any = 12,
    colorbar_min: float | str = "auto",
    segmented_colors: bool = True,
    max_triangle_edge: float | None = None,
    show_stations: bool = True,
    show_confidence_values: bool = False,
    confidence_value_step: int | None = None,
    confidence_value_fmt: str = "{:.2f}",
    confidence_value_fontsize: float = 7.0,
    show_contour_lines: bool = False,
    contour_line_levels: Any = None,
    contour_line_colors: Any = "0.25",
    contour_linewidths: Any = 0.65,
    contour_linestyles: Any = "solid",
    contour_labels: bool = True,
    contour_label_fmt: str | Any = "%.2f",
    contour_label_fontsize: float = 7.0,
    contour_label_inline: bool = True,
    show_threshold_contours: bool = True,
    threshold_line_color: Any = "black",
    threshold_linewidth: float = 1.35,
    threshold_linestyle: str = "solid",
    station_labels: bool = False,
    station_label_step: int | None = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    boundary_levels: Any = None,
    cmap: str = "RdYlGn",
    marker_size: float = 72.0,
    figsize: tuple[float, float] = (7.5, 5.5),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Map station-level confidence at geographic or projected coordinates.

    Route and scatter modes do not interpolate between stations. Contour mode
    uses a triangular surface inside the survey's convex hull and rejects
    collinear station layouts, so a single profile cannot accidentally appear
    to provide two-dimensional spatial coverage.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Input stations accepted by :func:`station_confidence_table`.
    method : {"presence", "composite"}
        Confidence scoring method.
    coordinate_system : {"auto", "geographic", "projected"}
        ``"auto"`` prefers longitude/latitude and falls back to
        easting/northing. The selected pair must be available per station.
    mode : {"auto", "scatter", "route", "contour"}
        Map representation. ``"auto"`` selects ``"route"`` when ``connect``
        is true and ``"scatter"`` otherwise. ``"contour"`` requires at
        least three non-collinear station coordinates.
    line_labels : mapping, sequence, str, or None
        Optional survey-line membership. A mapping is keyed by station name;
        a sequence follows table order; one string assigns every station to
        that line. When omitted all stations form one route.
    connect : bool
        Overlay routes within each line. In contour mode this is useful for
        retaining the acquisition geometry above the interpolated surface.
    contour_levels : int or array-like
        Number of discrete filled intervals, or explicit boundaries.
    colorbar_min : float or "auto"
        Lower contour/colorbar boundary. ``"auto"`` rounds down the observed
        minimum: to a 0.05 step when all CR values are at least 0.5, otherwise
        to a 0.1 step. Applies consistently to contour, route, and scatter
        modes. Pass ``0.0`` to retain the original full 0--1 scale.
    segmented_colors : bool
        Use discrete color intervals with explicit breaks at 0.50,
        ``ci_lo``, ``ci_hi``, and 1.00. Set false for a continuous gradient.
    max_triangle_edge : float or None
        Optional maximum triangle-edge length, expressed in the selected map
        units. Triangles crossing a larger unsurveyed gap are masked.
    show_stations : bool
        Draw confidence-coloured station markers over the map.
    show_confidence_values : bool
        Annotate the numeric confidence beside each displayed station.
    confidence_value_step : int or None
        Label every Nth station. ``None`` automatically thins surveys with
        more than 20 stations while retaining the final station value.
    confidence_value_fmt : str
        Python format string for station confidence annotations.
    confidence_value_fontsize : float
        Font size for station confidence annotations.
    show_contour_lines : bool
        Overlay ordinary confidence isolines and, by default, their values.
    contour_line_levels : array-like or None
        Isoline values. ``None`` uses the discrete filled-contour boundaries
        that fall within the observed confidence range.
    contour_line_colors, contour_linewidths, contour_linestyles
        Matplotlib styling passed to :meth:`~matplotlib.axes.Axes.tricontour`.
    contour_labels : bool
        Label ordinary contour lines with their confidence values.
    contour_label_fmt : str, mapping, or callable
        Label formatter passed to :meth:`~matplotlib.axes.Axes.clabel`.
    contour_label_fontsize : float
        Numeric contour-label size.
    contour_label_inline : bool
        Remove the line beneath each numeric contour label.
    show_threshold_contours : bool
        Delineate 0.50, ``ci_lo``, and ``ci_hi`` on the contour surface when
        those values are crossed by the observed data. The colorbar always
        marks every applicable boundary.
    threshold_line_color, threshold_linewidth, threshold_linestyle
        Styling for the emphasized confidence-class borders.
    station_labels : bool
        Annotate station names. ``station_label_step`` controls thinning.
    ci_hi, ci_lo : float
        Safe and recoverable/review confidence thresholds.
    boundary_levels : array-like or None
        Confidence-class borders used by segmented colors, emphasized
        isolines, and colorbar ticks. ``None`` uses
        ``(0.50, ci_lo, ci_hi, 1.00)``. For example, pass
        ``(0.50, 0.85, 0.90, 1.00)`` for boundary-only contours at those
        values.
    cmap, marker_size, figsize
        Matplotlib appearance controls.
    recursive, on_dup, strict, verbose
        Passed through the standard EMTools site-loading API.
    ax : matplotlib.axes.Axes or None
        Existing axes, or ``None`` to create one.

    Returns
    -------
    matplotlib.axes.Axes
        The map axes. Its ``_pycsamt_coordinate_system`` attribute records
        the coordinate system selected by ``"auto"``.
    """
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    mode = str(mode).lower()
    if mode not in {"auto", "scatter", "route", "contour"}:
        raise ValueError(
            "mode must be 'auto', 'scatter', 'route', or 'contour'."
        )
    if mode == "auto":
        mode = "route" if connect else "scatter"
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    tb = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    if tb.empty:
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return ax

    geographic = np.isfinite(tb[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(tb[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    tb = tb.loc[valid].copy()
    ax._pycsamt_coordinate_system = coordinate_system
    if tb.empty:
        ax.text(
            0.5,
            0.5,
            f"no {coordinate_system} coordinates",
            ha="center",
            va="center",
        )
        ax.set_xlabel(
            "Longitude (deg)"
            if coordinate_system == "geographic"
            else "Easting (m)"
        )
        ax.set_ylabel(
            "Latitude (deg)"
            if coordinate_system == "geographic"
            else "Northing (m)"
        )
        return ax

    if line_labels is None:
        labels = pd.Series("survey", index=tb.index, dtype=object)
    elif isinstance(line_labels, str):
        labels = pd.Series(line_labels, index=tb.index, dtype=object)
    elif hasattr(line_labels, "get"):
        labels = tb["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(valid):
            raise ValueError(
                "line_labels must contain one value per confidence-table row."
            )
        labels = pd.Series(values, index=valid.index).loc[tb.index]
    tb["_line"] = labels.astype(str)

    plot_vmin = 0.0
    norm = None
    filled_levels = contour_levels
    threshold_values = np.asarray(
        (0.5, ci_lo, ci_hi, 1.0)
        if boundary_levels is None
        else boundary_levels,
        dtype=float,
    )
    threshold_values = np.unique(threshold_values)
    if (
        threshold_values.ndim != 1
        or threshold_values.size < 2
        or not np.isfinite(threshold_values).all()
        or threshold_values[0] < 0.0
        or threshold_values[-1] > 1.0
    ):
        raise ValueError(
            "boundary_levels must be finite values within [0, 1]."
        )
    spatial_boundaries = threshold_values[~np.isclose(threshold_values, 1.0)]
    observed_min = float(tb["confidence"].min())
    if isinstance(colorbar_min, str):
        if colorbar_min.lower() != "auto":
            raise ValueError("colorbar_min must be a number or 'auto'.")
        step = 0.05 if observed_min >= 0.5 else 0.1
        plot_vmin = max(
            0.0,
            np.floor((observed_min + 1e-12) / step) * step,
        )
        plot_vmin = min(plot_vmin, 0.95)
    else:
        plot_vmin = float(colorbar_min)
        if not np.isfinite(plot_vmin) or not 0.0 <= plot_vmin < 1.0:
            raise ValueError("colorbar_min must satisfy 0 <= value < 1.")
    if segmented_colors or mode == "contour":
        if np.isscalar(contour_levels):
            n_levels = int(contour_levels)
            if n_levels < 2:
                raise ValueError("contour_levels must be at least 2.")
            filled_levels = np.linspace(plot_vmin, 1.0, n_levels + 1)
        else:
            filled_levels = np.asarray(contour_levels, dtype=float)
        filled_levels = np.unique(
            np.r_[
                filled_levels,
                [v for v in threshold_values if plot_vmin <= v <= 1.0],
                plot_vmin,
                1.0,
            ]
        )
        if filled_levels.size < 2 or np.any(np.diff(filled_levels) <= 0.0):
            raise ValueError(
                "contour_levels must define increasing boundaries."
            )
        from matplotlib.colors import BoundaryNorm

        norm = BoundaryNorm(filled_levels, plt.get_cmap(cmap).N, clip=True)
    else:
        from matplotlib.colors import Normalize

        norm = Normalize(plot_vmin, 1.0, clip=True)

    draw_routes = connect and mode != "scatter"
    if mode == "contour" and line_labels is None:
        draw_routes = False
    if draw_routes:
        for _, group in tb.groupby("_line", sort=False):
            group = _order_map_route(group, xkey, ykey)
            ax.plot(
                group[xkey],
                group[ykey],
                color="0.35",
                lw=1.2,
                alpha=0.75,
                zorder=2 if mode == "contour" else 1,
            )
    mappable = None
    if mode == "contour":
        from matplotlib.tri import Triangulation

        xy = tb[[xkey, ykey]].to_numpy(dtype=float)
        if len(xy) < 3 or np.linalg.matrix_rank(xy - xy.mean(axis=0)) < 2:
            raise ValueError(
                "contour mode requires at least three non-collinear stations; "
                "use mode='route' for a single survey line."
            )
        try:
            triangulation = Triangulation(xy[:, 0], xy[:, 1])
        except (RuntimeError, ValueError) as exc:
            raise ValueError(
                "confidence contour triangulation failed; check for duplicate "
                "or nearly collinear station coordinates."
            ) from exc
        if max_triangle_edge is not None:
            max_triangle_edge = float(max_triangle_edge)
            if not np.isfinite(max_triangle_edge) or max_triangle_edge <= 0.0:
                raise ValueError(
                    "max_triangle_edge must be a positive finite value."
                )
            triangles = triangulation.triangles
            vertices = xy[triangles]
            edge_lengths = np.linalg.norm(
                vertices - np.roll(vertices, -1, axis=1), axis=2
            )
            triangulation.set_mask(
                (edge_lengths > max_triangle_edge).any(axis=1)
            )
            if np.all(triangulation.mask):
                raise ValueError(
                    "max_triangle_edge masks every contour triangle; increase "
                    "the limit or use mode='route'."
                )
        mappable = ax.tricontourf(
            triangulation,
            tb["confidence"].to_numpy(dtype=float),
            levels=filled_levels,
            cmap=cmap,
            norm=norm,
            extend="neither",
            zorder=0,
        )
        observed_max = float(tb["confidence"].max())
        if show_contour_lines:
            line_levels = (
                np.asarray(filled_levels, dtype=float)
                if contour_line_levels is None
                else np.asarray(contour_line_levels, dtype=float)
            )
            line_levels = np.unique(
                line_levels[
                    (line_levels >= observed_min)
                    & (line_levels <= observed_max)
                ]
            )
            if show_threshold_contours:
                line_levels = np.asarray(
                    [
                        value
                        for value in line_levels
                        if not any(
                            np.isclose(value, threshold)
                            for threshold in spatial_boundaries
                        )
                    ],
                    dtype=float,
                )
            if line_levels.size:
                isolines = ax.tricontour(
                    triangulation,
                    tb["confidence"].to_numpy(dtype=float),
                    levels=line_levels,
                    colors=contour_line_colors,
                    linewidths=contour_linewidths,
                    linestyles=contour_linestyles,
                    zorder=1,
                )
                if contour_labels:
                    ax.clabel(
                        isolines,
                        fmt=contour_label_fmt,
                        fontsize=contour_label_fontsize,
                        inline=contour_label_inline,
                    )
        crossed = [
            value
            for value in spatial_boundaries
            if observed_min <= value <= observed_max
        ]
        if show_threshold_contours and crossed:
            borders = ax.tricontour(
                triangulation,
                tb["confidence"].to_numpy(dtype=float),
                levels=crossed,
                colors=threshold_line_color,
                linewidths=threshold_linewidth,
                linestyles=threshold_linestyle,
                zorder=1.5,
            )
            if contour_labels:
                ax.clabel(
                    borders,
                    fmt=contour_label_fmt,
                    fontsize=contour_label_fontsize,
                    inline=contour_label_inline,
                )
    if show_stations:
        points = ax.scatter(
            tb[xkey],
            tb[ykey],
            c=tb["confidence"],
            cmap=cmap,
            norm=norm,
            vmin=None,
            vmax=None,
            s=marker_size,
            edgecolors="black",
            linewidths=0.8,
            zorder=3 if mode == "contour" else 2,
        )
        if mappable is None:
            mappable = points
    if mappable is not None:
        colorbar = ax.figure.colorbar(mappable, ax=ax, pad=0.02)
    else:
        from matplotlib.cm import ScalarMappable

        mappable = ScalarMappable(norm=norm, cmap=cmap)
        colorbar = ax.figure.colorbar(mappable, ax=ax, pad=0.02)
    colorbar.set_label("Confidence ratio")
    colorbar_ticks = np.unique(
        np.r_[plot_vmin, threshold_values[threshold_values >= plot_vmin]]
    ).tolist()
    colorbar.set_ticks(colorbar_ticks)
    if segmented_colors:
        for value in colorbar_ticks[1:]:
            colorbar.ax.axhline(value, color="black", lw=0.8, alpha=0.85)

    if station_labels:
        step = max(1, int(station_label_step or np.ceil(len(tb) / 20)))
        keep = list(range(0, len(tb), step))
        if len(tb) - 1 not in keep:
            keep.append(len(tb) - 1)
        for i in keep:
            row = tb.iloc[i]
            ax.annotate(
                str(row["station"]),
                (row[xkey], row[ykey]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7,
            )
    if show_confidence_values:
        value_step = max(
            1,
            int(confidence_value_step or np.ceil(len(tb) / 20)),
        )
        value_indices = list(range(0, len(tb), value_step))
        if len(tb) - 1 not in value_indices:
            value_indices.append(len(tb) - 1)
        for i in value_indices:
            row = tb.iloc[i]
            try:
                label = confidence_value_fmt.format(float(row["confidence"]))
            except (AttributeError, ValueError):
                label = f"{float(row['confidence']):.2f}"
            ax.annotate(
                label,
                (row[xkey], row[ykey]),
                xytext=(4, -6),
                textcoords="offset points",
                ha="left",
                va="top",
                fontsize=confidence_value_fontsize,
                color="0.15",
            )
    if tb["_line"].nunique() > 1:
        for line, group in tb.groupby("_line", sort=False):
            row = group.iloc[len(group) // 2]
            ax.annotate(
                str(line),
                (row[xkey], row[ykey]),
                xytext=(5, -10),
                textcoords="offset points",
                fontsize=8,
                fontweight="bold",
            )
    ax.set_xlabel(
        "Longitude (deg)"
        if coordinate_system == "geographic"
        else "Easting (m)"
    )
    ax.set_ylabel(
        "Latitude (deg)"
        if coordinate_system == "geographic"
        else "Northing (m)"
    )
    ax.set_title(f"Station confidence map ({method}, {mode})", fontsize=10)
    ax.grid(True, ls=":", alpha=0.35)
    return ax


def plot_confidence_grid_map(
    sites: Any,
    *,
    method: str = "composite",
    coordinate_system: str = "auto",
    line_labels: Any = None,
    grid_shape: tuple[int, int] = (180, 180),
    interpolation: str = "linear",
    max_triangle_edge: float | None = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    levels: Any = None,
    cmap: str = "RdYlGn",
    segmented_colors: bool = True,
    show_grid_edges: bool = False,
    grid_edgecolor: str = "white",
    grid_linewidth: float = 0.15,
    show_threshold_contours: bool = True,
    threshold_colors: Any = ("#9a6700", "#1b5e20"),
    threshold_linewidths: Any = (1.4, 1.6),
    threshold_linestyles: Any = ("--", "-"),
    contour_labels: bool = True,
    show_stations: bool = True,
    marker_size: float = 18.0,
    connect: bool = True,
    line_names: bool = True,
    nodata_color: str = "white",
    map_aspect: str = "auto",
    figsize: tuple[float, float] = (8.2, 5.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Render confidence on a regular map grid inside the survey hull.

    Unlike :func:`plot_confidence_map` contour mode, this function displays
    the actual regular cells used by raster and Surfer-style workflows.
    Cells outside the triangulated survey footprint remain blank.
    ``grid_shape`` is ``(nx, ny)`` and ``interpolation`` may be ``'linear'``
    or ``'cubic'``.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    interpolation = str(interpolation).lower()
    if interpolation not in {"linear", "cubic"}:
        raise ValueError("interpolation must be 'linear' or 'cubic'.")
    map_aspect = str(map_aspect).lower()
    if map_aspect not in {"auto", "equal", "geographic"}:
        raise ValueError(
            "map_aspect must be 'auto', 'equal', or 'geographic'."
        )
    nx, ny = int(grid_shape[0]), int(grid_shape[1])
    if nx < 2 or ny < 2:
        raise ValueError("grid_shape values must both be at least 2.")

    table = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    geographic = np.isfinite(table[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(table[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    table = table.loc[valid & np.isfinite(table["confidence"])].copy()
    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax._pycsamt_coordinate_system = coordinate_system
    if table.empty:
        ax.text(0.5, 0.5, "no mapped stations", ha="center", va="center")
        return ax

    if line_labels is None:
        labels = pd.Series("survey", index=table.index, dtype=object)
    elif isinstance(line_labels, str):
        labels = pd.Series(line_labels, index=table.index, dtype=object)
    elif hasattr(line_labels, "get"):
        labels = table["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(valid):
            raise ValueError(
                "line_labels must contain one value per table row."
            )
        labels = pd.Series(values, index=valid.index).loc[table.index]
    table["_line"] = labels.astype(str)

    grouped = (
        table.groupby([xkey, ykey], as_index=False)
        .agg({"confidence": "mean"})
        .copy()
    )
    # ``.to_numpy()`` can hand back a read-only view (pandas copy-on-write);
    # matplotlib's Triangulation / CubicTriInterpolator assume writeable
    # input and raise "assignment destination is read-only" otherwise.
    xy = np.array(grouped[[xkey, ykey]].to_numpy(dtype=float), dtype=float)
    if len(xy) < 3 or np.linalg.matrix_rank(xy - xy.mean(axis=0)) < 2:
        raise ValueError(
            "grid maps require at least three non-collinear station "
            "coordinates; use plot_confidence_map(..., mode='route') for "
            "one survey line."
        )
    from matplotlib.tri import (
        CubicTriInterpolator,
        LinearTriInterpolator,
        Triangulation,
    )

    triangulation = Triangulation(xy[:, 0], xy[:, 1])
    if max_triangle_edge is not None:
        limit = float(max_triangle_edge)
        if not np.isfinite(limit) or limit <= 0.0:
            raise ValueError("max_triangle_edge must be positive and finite.")
        vertices = xy[triangulation.triangles]
        edges = np.linalg.norm(
            vertices - np.roll(vertices, -1, axis=1), axis=2
        )
        triangulation.set_mask((edges > limit).any(axis=1))
        if np.all(triangulation.mask):
            raise ValueError("max_triangle_edge masks every grid triangle.")
    x_axis = np.linspace(float(xy[:, 0].min()), float(xy[:, 0].max()), nx)
    y_axis = np.linspace(float(xy[:, 1].min()), float(xy[:, 1].max()), ny)
    grid_x, grid_y = np.meshgrid(x_axis, y_axis)
    interpolator_class = (
        CubicTriInterpolator
        if interpolation == "cubic"
        else LinearTriInterpolator
    )
    interpolator = interpolator_class(
        triangulation,
        np.array(grouped["confidence"].to_numpy(dtype=float), dtype=float),
    )
    confidence_grid = np.ma.asarray(
        interpolator(np.array(grid_x), np.array(grid_y))
    )
    confidence_grid = np.ma.masked_invalid(
        np.ma.clip(confidence_grid, 0.0, 1.0)
    )
    ax._pycsamt_grid_x = x_axis
    ax._pycsamt_grid_y = y_axis
    ax._pycsamt_confidence_grid = confidence_grid

    from matplotlib.colors import BoundaryNorm, Normalize

    if levels is None:
        levels = sorted({0.0, 0.25, 0.50, 0.75, ci_lo, ci_hi, 1.0})
    levels = np.unique(np.clip(np.asarray(levels, dtype=float), 0.0, 1.0))
    if len(levels) < 2 or levels[0] > 0.0 or levels[-1] < 1.0:
        raise ValueError("levels must span 0 to 1 with at least two values.")
    use_cmap = plt.get_cmap(cmap).copy()
    use_cmap.set_bad(nodata_color)
    norm = (
        BoundaryNorm(levels, use_cmap.N, clip=True)
        if segmented_colors
        else Normalize(0.0, 1.0, clip=True)
    )
    grid_plot = ax.pcolormesh(
        grid_x,
        grid_y,
        confidence_grid,
        shading="auto",
        cmap=use_cmap,
        norm=norm,
        edgecolors=grid_edgecolor if show_grid_edges else "none",
        linewidth=grid_linewidth if show_grid_edges else 0.0,
        rasterized=True,
        zorder=0,
    )
    if show_threshold_contours:
        finite = confidence_grid.compressed()
        boundaries = [
            value
            for value in (ci_lo, ci_hi)
            if finite.size and finite.min() < value < finite.max()
        ]
        if boundaries:
            colors = list(np.atleast_1d(threshold_colors))[: len(boundaries)]
            widths = list(np.atleast_1d(threshold_linewidths))[
                : len(boundaries)
            ]
            styles = list(np.atleast_1d(threshold_linestyles))[
                : len(boundaries)
            ]
            contours = ax.contour(
                grid_x,
                grid_y,
                confidence_grid,
                levels=boundaries,
                colors=colors,
                linewidths=widths,
                linestyles=styles,
                zorder=2,
            )
            if contour_labels:
                ax.clabel(contours, fmt="CR %.2f", fontsize=7, inline=True)
    if connect:
        for _, group in table.groupby("_line", sort=False):
            group = _order_map_route(group, xkey, ykey)
            ax.plot(
                group[xkey],
                group[ykey],
                color="0.25",
                lw=0.7,
                alpha=0.72,
                zorder=3,
            )
    if show_stations:
        ax.scatter(
            table[xkey],
            table[ykey],
            c=table["confidence"],
            cmap=use_cmap,
            norm=norm,
            s=marker_size,
            edgecolors="white",
            linewidths=0.45,
            zorder=4,
        )
    if line_names and table["_line"].nunique() > 1:
        for line, group in table.groupby("_line", sort=False):
            row = group.loc[group[ykey].idxmax()]
            ax.annotate(
                str(line),
                (row[xkey], row[ykey]),
                xytext=(0, 4),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="semibold",
                zorder=5,
            )
    colorbar = ax.figure.colorbar(grid_plot, ax=ax, pad=0.018, fraction=0.045)
    colorbar.set_label("Confidence ratio", fontsize=9)
    colorbar.set_ticks(sorted({0.0, 0.5, ci_lo, ci_hi, 1.0}))
    colorbar.ax.tick_params(labelsize=7.5)
    ax.set_xlabel(
        "Longitude (°)" if coordinate_system == "geographic" else "Easting (m)"
    )
    ax.set_ylabel(
        "Latitude (°)" if coordinate_system == "geographic" else "Northing (m)"
    )
    ax.set_title(
        f"Regular-grid confidence map ({method}; {nx} × {ny})",
        loc="left",
        fontsize=11,
        fontweight="semibold",
    )
    xpad = max(float(x_axis[-1] - x_axis[0]) * 0.025, 1e-12)
    ypad = max(float(y_axis[-1] - y_axis[0]) * 0.035, 1e-12)
    ax.set_xlim(float(x_axis[0]) - xpad, float(x_axis[-1]) + xpad)
    ax.set_ylim(float(y_axis[0]) - ypad, float(y_axis[-1]) + ypad)
    ax.grid(True, ls=":", lw=0.45, color="0.72", alpha=0.55)
    ax.tick_params(labelsize=8, direction="out", length=3)
    if map_aspect == "geographic" and coordinate_system == "geographic":
        mean_lat = float(table[ykey].mean())
        ax.set_aspect(1.0 / max(np.cos(np.deg2rad(mean_lat)), 1e-6))
    elif map_aspect == "equal":
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_aspect("auto")
    return ax


def plot_confidence_component_map(
    sites: Any,
    *,
    method: str = "composite",
    components: Any = None,
    coordinate_system: str = "auto",
    line_labels: Any = None,
    connect: bool = True,
    cmap: str = "RdYlGn",
    marker_size: float = 38.0,
    station_labels: bool = False,
    station_label_step: int | None = None,
    line_names: bool = True,
    ncols: int = 4,
    figsize: tuple[float, float] | None = None,
    map_aspect: str = "auto",
    panel_letters: bool = True,
    colorbar_label: str = "Confidence component score",
    axes: Any = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Plot publication-style maps of confidence and its components.

    Every panel uses the same station geometry, extent, and fixed 0--1 color
    normalization, allowing direct scientific comparison between component
    scores. The default seven panels are overall confidence, coverage,
    uncertainty, off-diagonal consistency, diagonal leakage, phase
    smoothness, and spatial coherence.

    ``map_aspect="auto"`` is the compact publication default. Use
    ``"geographic"`` to preserve longitude/latitude ground proportions or
    ``"equal"`` for equal numeric axis units.
    """
    default_components = (
        "confidence",
        "coverage",
        "uncertainty",
        "offdiag",
        "diagonal",
        "phase",
        "spatial",
    )
    components = tuple(
        default_components if components is None else components
    )
    if not components:
        raise ValueError("components must contain at least one table column.")
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    tb = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    missing = [component for component in components if component not in tb]
    if missing:
        raise ValueError(f"unknown confidence components: {missing}.")
    geographic = np.isfinite(tb[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(tb[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    tb = tb.loc[valid].copy()

    ncols = max(1, min(int(ncols), len(components)))
    nrows = int(np.ceil(len(components) / ncols))
    map_aspect = str(map_aspect).lower()
    if map_aspect not in {"auto", "equal", "geographic"}:
        raise ValueError(
            "map_aspect must be 'auto', 'equal', or 'geographic'."
        )
    colorbar_axes = None
    if axes is None:
        if figsize is None:
            figsize = (2.85 * ncols, 3.05 * nrows)
        fig = plt.figure(figsize=figsize)
        grid = fig.add_gridspec(
            nrows,
            ncols,
            left=0.075,
            right=0.975,
            bottom=0.105,
            top=0.895,
            wspace=0.10,
            hspace=0.16,
        )
        flat_axes = []
        shared = None
        for index in range(len(components)):
            ax = fig.add_subplot(
                grid[index // ncols, index % ncols],
                sharex=shared,
                sharey=shared,
            )
            if shared is None:
                shared = ax
            flat_axes.append(ax)
        flat_axes = np.asarray(flat_axes, dtype=object)
        if len(components) < nrows * ncols:
            host = fig.add_subplot(
                grid[len(components) // ncols, len(components) % ncols]
            )
            host.set_axis_off()
            colorbar_axes = host.inset_axes([0.43, 0.06, 0.14, 0.88])
    else:
        flat_axes = np.asarray(axes, dtype=object).ravel()
        if flat_axes.size < len(components):
            raise ValueError(
                "axes must provide at least one axes per component."
            )
        fig = flat_axes[0].figure
    if tb.empty:
        flat_axes[0].text(
            0.5, 0.5, "no mapped stations", ha="center", va="center"
        )
        return fig

    if line_labels is None:
        labels = pd.Series("survey", index=tb.index, dtype=object)
    elif isinstance(line_labels, str):
        labels = pd.Series(line_labels, index=tb.index, dtype=object)
    elif hasattr(line_labels, "get"):
        labels = tb["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(valid):
            raise ValueError(
                "line_labels must contain one value per table row."
            )
        labels = pd.Series(values, index=valid.index).loc[tb.index]
    tb["_line"] = labels.astype(str)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize

    norm = Normalize(0.0, 1.0, clip=True)
    titles = {
        "confidence": "Overall confidence",
        "coverage": "Data coverage",
        "uncertainty": "Tensor uncertainty",
        "offdiag": "Off-diagonal consistency",
        "diagonal": "Diagonal leakage",
        "phase": "Phase smoothness",
        "spatial": "Spatial coherence",
    }
    letters = "abcdefghijklmnopqrstuvwxyz"
    xpad = max(float(tb[xkey].max() - tb[xkey].min()) * 0.03, 1e-12)
    ypad = max(float(tb[ykey].max() - tb[ykey].min()) * 0.03, 1e-12)
    for panel_index, (ax, component) in enumerate(zip(flat_axes, components)):
        if connect:
            for _, group in tb.groupby("_line", sort=False):
                group = _order_map_route(group, xkey, ykey)
                ax.plot(
                    group[xkey],
                    group[ykey],
                    color="0.48",
                    lw=0.75,
                    alpha=0.75,
                    zorder=1,
                )
        finite = np.isfinite(tb[component].to_numpy(dtype=float))
        if finite.any():
            ax.scatter(
                tb.loc[finite, xkey],
                tb.loc[finite, ykey],
                c=tb.loc[finite, component],
                cmap=cmap,
                norm=norm,
                s=marker_size,
                edgecolors="black",
                linewidths=0.55,
                zorder=2,
            )
        else:
            ax.text(
                0.5,
                0.5,
                "not available",
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="0.35",
            )
        title = titles.get(component, component.replace("_", " ").title())
        if panel_letters:
            letter = (
                letters[panel_index]
                if panel_index < len(letters)
                else str(panel_index + 1)
            )
            title = f"({letter})  {title}"
        ax.set_title(title, loc="left", fontsize=9.2, fontweight="semibold")
        ax.set_xlim(float(tb[xkey].min()) - xpad, float(tb[xkey].max()) + xpad)
        ax.set_ylim(float(tb[ykey].min()) - ypad, float(tb[ykey].max()) + ypad)
        ax.grid(True, ls=":", lw=0.55, color="0.78", alpha=0.75)
        row, col = divmod(panel_index, ncols)
        ax.tick_params(
            labelsize=7.5,
            direction="out",
            length=3,
            labelleft=col == 0,
            labelbottom=row == nrows - 1,
        )
        if map_aspect == "geographic" and coordinate_system == "geographic":
            mean_lat = float(tb[ykey].mean())
            ax.set_aspect(1.0 / max(np.cos(np.deg2rad(mean_lat)), 1e-6))
        elif map_aspect == "equal":
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.set_aspect("auto")
        if station_labels and panel_index == 0:
            step = max(1, int(station_label_step or np.ceil(len(tb) / 18)))
            for i in range(0, len(tb), step):
                row = tb.iloc[i]
                ax.annotate(
                    str(row["station"]),
                    (row[xkey], row[ykey]),
                    xytext=(3, 3),
                    textcoords="offset points",
                    fontsize=6,
                )
        if line_names and panel_index == 0 and tb["_line"].nunique() > 1:
            for line, group in tb.groupby("_line", sort=False):
                row = group.loc[group[ykey].idxmax()]
                ax.annotate(
                    str(line),
                    (row[xkey], row[ykey]),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=6.5,
                    fontweight="bold",
                )
    visible_axes = list(flat_axes[: len(components)])
    fig.supxlabel(
        "Longitude (°)"
        if coordinate_system == "geographic"
        else "Easting (m)",
        fontsize=9,
    )
    fig.supylabel(
        "Latitude (°)"
        if coordinate_system == "geographic"
        else "Northing (m)",
        fontsize=9,
    )
    color_mappable = ScalarMappable(norm=norm, cmap=cmap)
    if colorbar_axes is not None:
        colorbar = fig.colorbar(color_mappable, cax=colorbar_axes)
    else:
        colorbar = fig.colorbar(
            color_mappable,
            ax=visible_axes,
            fraction=0.025,
            pad=0.018,
            aspect=32,
        )
    colorbar.set_label(colorbar_label, fontsize=9)
    colorbar.set_ticks([0.0, 0.25, 0.50, 0.75, 1.0])
    colorbar.ax.tick_params(labelsize=7.5)
    fig.suptitle(
        f"Station confidence components ({method})",
        fontsize=11,
        fontweight="semibold",
    )
    return fig


def plot_confidence_method_comparison(
    sites: Any,
    *,
    coordinate_system: str = "auto",
    line_labels: Any = None,
    connect: bool = True,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    confidence_cmap: str = "RdYlGn",
    difference_cmap: str = "RdBu",
    marker_size: float = 42.0,
    station_labels: bool = False,
    station_label_step: int | None = None,
    line_names: bool = True,
    show_statistics: bool = True,
    difference_limit: float | None = None,
    map_aspect: str = "auto",
    figsize: tuple[float, float] = (11.2, 4.5),
    axes: Any = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Compare presence and composite confidence on matched station maps.

    Panels show presence confidence, composite confidence, and
    ``composite - presence``. The confidence panels share a fixed 0--1 color
    scale; the difference panel uses a symmetric zero-centred scale so score
    gains and penalties remain visually comparable.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    map_aspect = str(map_aspect).lower()
    if map_aspect not in {"auto", "equal", "geographic"}:
        raise ValueError(
            "map_aspect must be 'auto', 'equal', or 'geographic'."
        )
    presence = station_confidence_table(
        sites,
        method="presence",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    composite = station_confidence_table(
        sites,
        method="composite",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    coordinate_columns = ["longitude", "latitude", "easting", "northing"]
    table = composite[
        ["station", "confidence", "confidence_err", *coordinate_columns]
    ].rename(
        columns={
            "confidence": "composite",
            "confidence_err": "composite_err",
        }
    )
    table = table.merge(
        presence[["station", "confidence", "confidence_err"]].rename(
            columns={
                "confidence": "presence",
                "confidence_err": "presence_err",
            }
        ),
        on="station",
        how="inner",
    )
    table["difference"] = table["composite"] - table["presence"]
    geographic = np.isfinite(table[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(table[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    table = table.loc[valid].copy()

    if axes is None:
        fig = plt.figure(figsize=figsize)
        grid = fig.add_gridspec(
            3,
            3,
            height_ratios=(1.0, 0.045, 0.065),
            left=0.075,
            right=0.975,
            bottom=0.075,
            top=0.88,
            wspace=0.08,
            hspace=0.12,
        )
        map_axes = []
        shared = None
        for index in range(3):
            ax = fig.add_subplot(grid[0, index], sharex=shared, sharey=shared)
            if shared is None:
                shared = ax
            map_axes.append(ax)
        longitude_ax = fig.add_subplot(grid[1, :])
        longitude_ax.set_axis_off()
        longitude_ax.text(
            0.5,
            0.5,
            "Longitude (°)"
            if coordinate_system == "geographic"
            else "Easting (m)",
            ha="center",
            va="center",
            fontsize=9,
        )
        confidence_cax = fig.add_subplot(grid[2, :2])
        difference_cax = fig.add_subplot(grid[2, 2])
    else:
        map_axes = list(np.asarray(axes, dtype=object).ravel())
        if len(map_axes) < 3:
            raise ValueError("axes must provide three map axes.")
        map_axes = map_axes[:3]
        fig = map_axes[0].figure
        confidence_cax = difference_cax = None
    if table.empty:
        map_axes[0].text(
            0.5, 0.5, "no mapped stations", ha="center", va="center"
        )
        return fig

    if line_labels is None:
        labels = pd.Series("survey", index=table.index, dtype=object)
    elif isinstance(line_labels, str):
        labels = pd.Series(line_labels, index=table.index, dtype=object)
    elif hasattr(line_labels, "get"):
        labels = table["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(valid):
            raise ValueError(
                "line_labels must contain one value per table row."
            )
        labels = pd.Series(values, index=valid.index).loc[table.index]
    table["_line"] = labels.astype(str)

    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import Normalize, TwoSlopeNorm

    confidence_norm = Normalize(0.0, 1.0, clip=True)
    observed_abs = float(np.nanmax(np.abs(table["difference"])))
    if difference_limit is None:
        difference_limit = max(0.05, np.ceil(observed_abs * 20.0) / 20.0)
    difference_limit = float(difference_limit)
    if not np.isfinite(difference_limit) or difference_limit <= 0.0:
        raise ValueError("difference_limit must be positive and finite.")
    difference_norm = TwoSlopeNorm(
        vmin=-difference_limit,
        vcenter=0.0,
        vmax=difference_limit,
    )
    panel_specs = (
        (
            "presence",
            "(a)  Presence confidence",
            confidence_cmap,
            confidence_norm,
        ),
        (
            "composite",
            "(b)  Composite confidence",
            confidence_cmap,
            confidence_norm,
        ),
        (
            "difference",
            "(c)  Composite − presence",
            difference_cmap,
            difference_norm,
        ),
    )
    xpad = max(float(table[xkey].max() - table[xkey].min()) * 0.03, 1e-12)
    ypad = max(float(table[ykey].max() - table[ykey].min()) * 0.03, 1e-12)
    for panel_index, (ax, (column, title, panel_cmap, norm)) in enumerate(
        zip(map_axes, panel_specs)
    ):
        if connect:
            for _, group in table.groupby("_line", sort=False):
                group = _order_map_route(group, xkey, ykey)
                ax.plot(
                    group[xkey],
                    group[ykey],
                    color="0.48",
                    lw=0.8,
                    alpha=0.75,
                    zorder=1,
                )
        ax.scatter(
            table[xkey],
            table[ykey],
            c=table[column],
            cmap=panel_cmap,
            norm=norm,
            s=marker_size,
            edgecolors="black",
            linewidths=0.6,
            zorder=2,
        )
        ax.set_title(title, loc="left", fontsize=10, fontweight="semibold")
        ax.set_xlim(
            float(table[xkey].min()) - xpad, float(table[xkey].max()) + xpad
        )
        ax.set_ylim(
            float(table[ykey].min()) - ypad, float(table[ykey].max()) + ypad
        )
        ax.grid(True, ls=":", lw=0.55, color="0.78", alpha=0.75)
        ax.tick_params(
            labelsize=7.5,
            direction="out",
            length=3,
            labelleft=panel_index == 0,
        )
        if map_aspect == "geographic" and coordinate_system == "geographic":
            mean_lat = float(table[ykey].mean())
            ax.set_aspect(1.0 / max(np.cos(np.deg2rad(mean_lat)), 1e-6))
        elif map_aspect == "equal":
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.set_aspect("auto")
        if show_statistics:
            values = table[column].to_numpy(dtype=float)
            if column == "difference":
                summary = (
                    f"median = {np.nanmedian(values):+.2f}\n"
                    f"range = {np.nanmin(values):+.2f} to "
                    f"{np.nanmax(values):+.2f}"
                )
            else:
                below = int(np.sum(values < ci_lo))
                summary = (
                    f"median = {np.nanmedian(values):.2f}\n"
                    f"CR < {ci_lo:.2f}: {below}/{len(values)}"
                )
            ax.text(
                0.025,
                0.025,
                summary,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                fontsize=7,
                bbox={
                    "boxstyle": "round,pad=0.25",
                    "facecolor": "white",
                    "edgecolor": "0.55",
                    "alpha": 0.88,
                    "linewidth": 0.6,
                },
                zorder=4,
            )
    if station_labels:
        step = max(1, int(station_label_step or np.ceil(len(table) / 18)))
        for i in range(0, len(table), step):
            row = table.iloc[i]
            map_axes[0].annotate(
                str(row["station"]),
                (row[xkey], row[ykey]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6,
            )
    if line_names and table["_line"].nunique() > 1:
        for line, group in table.groupby("_line", sort=False):
            row = group.loc[group[ykey].idxmax()]
            map_axes[0].annotate(
                str(line),
                (row[xkey], row[ykey]),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                va="bottom",
                rotation=90,
                fontsize=6.5,
                fontweight="bold",
            )
    confidence_mappable = ScalarMappable(
        norm=confidence_norm, cmap=confidence_cmap
    )
    difference_mappable = ScalarMappable(
        norm=difference_norm, cmap=difference_cmap
    )
    if confidence_cax is not None:
        confidence_bar = fig.colorbar(
            confidence_mappable,
            cax=confidence_cax,
            orientation="horizontal",
        )
        difference_bar = fig.colorbar(
            difference_mappable,
            cax=difference_cax,
            orientation="horizontal",
        )
    else:
        confidence_bar = fig.colorbar(
            confidence_mappable,
            ax=map_axes[:2],
            orientation="horizontal",
            fraction=0.05,
            pad=0.14,
        )
        difference_bar = fig.colorbar(
            difference_mappable,
            ax=map_axes[2],
            orientation="horizontal",
            fraction=0.05,
            pad=0.14,
        )
    confidence_bar.set_label("Confidence ratio", fontsize=8.5)
    confidence_bar.set_ticks(sorted({0.0, ci_lo, ci_hi, 1.0}))
    difference_bar.set_label("ΔCR (composite − presence)", fontsize=8.5)
    difference_bar.set_ticks([-difference_limit, 0.0, difference_limit])
    for bar in (confidence_bar, difference_bar):
        bar.ax.tick_params(labelsize=7.5)
    if axes is not None:
        fig.supxlabel(
            "Longitude (°)"
            if coordinate_system == "geographic"
            else "Easting (m)",
            fontsize=9,
        )
    fig.supylabel(
        "Latitude (°)"
        if coordinate_system == "geographic"
        else "Northing (m)",
        fontsize=9,
    )
    fig.suptitle(
        "Confidence-method comparison",
        fontsize=12,
        fontweight="semibold",
    )
    return fig


def plot_confidence_before_after(
    before_sites: Any,
    after_sites: Any | None = None,
    *,
    method: str = "composite",
    before_method: str | None = None,
    after_method: str | None = None,
    line_labels: Any = None,
    before_label: str = "Before",
    after_label: str = "After",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    change_tolerance: float = 0.01,
    show_errorbars: bool = True,
    show_station_labels: bool | str = "auto",
    station_label_step: int | None = None,
    show_line_panel: bool = True,
    improvement_color: str = "#1a9850",
    degradation_color: str = "#d73027",
    stable_color: str = "#7f7f7f",
    delta_limit: float | None = None,
    marker_size: float = 38.0,
    figsize: tuple[float, float] = (12.2, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    axes: Any = None,
) -> plt.Figure:
    """Compare matched station confidence before and after processing.

    ``after_sites`` defaults to ``before_sites`` so scoring-method changes can
    also be audited explicitly with ``before_method`` and ``after_method``.
    Stations are matched by name; unmatched stations are reported on the
    returned figure but are excluded from paired change statistics.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    tolerance = float(change_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("change_tolerance must be non-negative and finite.")
    before_method = str(before_method or method).lower()
    after_method = str(after_method or method).lower()
    for name in (before_method, after_method):
        if name not in {"presence", "composite"}:
            raise ValueError(
                "confidence methods must be presence or composite."
            )
    if after_sites is None:
        after_sites = before_sites
    before = station_confidence_table(
        before_sites,
        method=before_method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    after = station_confidence_table(
        after_sites,
        method=after_method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    )
    before_columns = [
        "station",
        "confidence",
        "confidence_err",
        "n_freq",
        "n_ok",
    ]
    after_columns = [
        "station",
        "confidence",
        "confidence_err",
        "n_freq",
        "n_ok",
    ]
    comparison = (
        before[before_columns]
        .rename(
            columns={
                "confidence": "before",
                "confidence_err": "before_err",
                "n_freq": "before_n_freq",
                "n_ok": "before_n_ok",
            }
        )
        .merge(
            after[after_columns].rename(
                columns={
                    "confidence": "after",
                    "confidence_err": "after_err",
                    "n_freq": "after_n_freq",
                    "n_ok": "after_n_ok",
                }
            ),
            on="station",
            how="inner",
        )
    )
    comparison = comparison.loc[
        np.isfinite(comparison[["before", "after"]]).all(axis=1)
    ].copy()
    comparison["delta"] = comparison["after"] - comparison["before"]
    comparison["status"] = np.select(
        [
            comparison["delta"] > tolerance,
            comparison["delta"] < -tolerance,
        ],
        ["improved", "degraded"],
        default="stable",
    )
    if line_labels is None:
        comparison["_line"] = "Survey"
    elif isinstance(line_labels, str):
        comparison["_line"] = line_labels
    elif hasattr(line_labels, "get"):
        comparison["_line"] = (
            comparison["station"].astype(str).map(line_labels).fillna("Survey")
        )
    else:
        values = list(line_labels)
        if len(values) != len(before):
            raise ValueError(
                "line_labels must contain one value per before-data station."
            )
        mapping = dict(zip(before["station"].astype(str), values))
        comparison["_line"] = comparison["station"].astype(str).map(mapping)

    n_panels = 3 if show_line_panel else 2
    if axes is None:
        widths = (1.0, 1.45, 1.0) if show_line_panel else (1.0, 1.45)
        fig, panel_axes = plt.subplots(
            1,
            n_panels,
            figsize=figsize,
            gridspec_kw={"width_ratios": widths},
            constrained_layout=True,
        )
    else:
        panel_axes = list(np.asarray(axes, dtype=object).ravel())
        if len(panel_axes) < n_panels:
            raise ValueError(f"axes must provide at least {n_panels} axes.")
        panel_axes = panel_axes[:n_panels]
        fig = panel_axes[0].figure
    agreement_ax, delta_ax = panel_axes[:2]
    line_ax = panel_axes[2] if show_line_panel else None
    fig._pycsamt_before_after_table = comparison
    fig._pycsamt_unmatched_before = sorted(
        set(before["station"].astype(str)) - set(after["station"].astype(str))
    )
    fig._pycsamt_unmatched_after = sorted(
        set(after["station"].astype(str)) - set(before["station"].astype(str))
    )
    if comparison.empty:
        agreement_ax.text(
            0.5, 0.5, "no matched stations", ha="center", va="center"
        )
        return fig
    colors = comparison["status"].map(
        {
            "improved": improvement_color,
            "degraded": degradation_color,
            "stable": stable_color,
        }
    )
    if show_errorbars:
        agreement_ax.errorbar(
            comparison["before"],
            comparison["after"],
            xerr=comparison["before_err"],
            yerr=comparison["after_err"],
            fmt="none",
            ecolor="0.65",
            elinewidth=0.55,
            alpha=0.55,
            zorder=1,
        )
    agreement_ax.scatter(
        comparison["before"],
        comparison["after"],
        c=colors,
        s=marker_size,
        edgecolors="0.15",
        linewidths=0.5,
        alpha=0.88,
        zorder=2,
    )
    agreement_ax.plot([0, 1], [0, 1], color="0.25", lw=1.0, ls="--")
    for threshold in (ci_lo, ci_hi):
        agreement_ax.axvline(threshold, color="0.55", lw=0.65, ls=":")
        agreement_ax.axhline(threshold, color="0.55", lw=0.65, ls=":")
    counts = comparison["status"].value_counts()
    summary = (
        f"Improved  {counts.get('improved', 0):3d}\n"
        f"Stable    {counts.get('stable', 0):3d}\n"
        f"Degraded  {counts.get('degraded', 0):3d}\n"
        f"Median Δ  {np.nanmedian(comparison['delta']):+.3f}"
    )
    agreement_ax.text(
        0.025,
        0.975,
        summary,
        transform=agreement_ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.8,
        family="monospace",
        bbox={
            "boxstyle": "round,pad=0.3",
            "facecolor": "white",
            "edgecolor": "0.55",
            "alpha": 0.90,
            "linewidth": 0.6,
        },
    )
    agreement_ax.set_xlim(0.0, 1.02)
    agreement_ax.set_ylim(0.0, 1.02)
    agreement_ax.set_xlabel(f"{before_label} confidence")
    agreement_ax.set_ylabel(f"{after_label} confidence")
    agreement_ax.set_title(
        "(a)  Paired agreement", loc="left", fontsize=10, fontweight="semibold"
    )

    ordered = comparison.sort_values(
        "delta", ascending=True, kind="stable"
    ).copy()
    positions = np.arange(len(ordered)) + 1
    delta_colors = ordered["status"].map(
        {
            "improved": improvement_color,
            "degraded": degradation_color,
            "stable": stable_color,
        }
    )
    delta_ax.bar(
        positions,
        ordered["delta"],
        color=delta_colors,
        edgecolor="none",
        width=0.82,
        zorder=2,
    )
    delta_ax.axhline(0.0, color="0.18", lw=0.9)
    delta_ax.axhspan(
        -tolerance, tolerance, color=stable_color, alpha=0.12, lw=0
    )
    if delta_limit is None:
        observed = float(np.nanmax(np.abs(comparison["delta"])))
        delta_limit = max(tolerance * 1.5, observed * 1.12, 0.05)
    delta_limit = float(delta_limit)
    if not np.isfinite(delta_limit) or delta_limit <= 0.0:
        raise ValueError("delta_limit must be positive and finite.")
    delta_ax.set_ylim(-delta_limit, delta_limit)
    delta_ax.set_xlim(0.25, len(ordered) + 0.75)
    delta_ax.set_xlabel("Stations ordered by confidence change")
    delta_ax.set_ylabel(f"ΔCR ({after_label} − {before_label})")
    delta_ax.set_title(
        "(b)  Station-level change",
        loc="left",
        fontsize=10,
        fontweight="semibold",
    )
    if show_station_labels == "auto":
        annotate = len(ordered) <= 32
    elif isinstance(show_station_labels, (bool, np.bool_)):
        annotate = bool(show_station_labels)
    else:
        raise ValueError("show_station_labels must be bool or 'auto'.")
    if annotate:
        step = max(1, int(station_label_step or 1))
        ticks = list(range(0, len(ordered), step))
        if len(ordered) - 1 not in ticks:
            ticks.append(len(ordered) - 1)
        delta_ax.set_xticks(np.asarray(ticks) + 1)
        delta_ax.set_xticklabels(
            ordered.iloc[ticks]["station"].astype(str), rotation=65, ha="right"
        )
    else:
        delta_ax.set_xticks([])

    if line_ax is not None:
        records = []
        for line, group in comparison.groupby("_line", sort=False):
            q1, median, q3 = np.nanpercentile(group["delta"], [25, 50, 75])
            records.append(
                {
                    "line": str(line),
                    "median_delta": float(median),
                    "q1": float(q1),
                    "q3": float(q3),
                    "improved": int(np.sum(group["status"] == "improved")),
                    "degraded": int(np.sum(group["status"] == "degraded")),
                    "n": len(group),
                }
            )
        line_table = pd.DataFrame(records).sort_values(
            "median_delta", ascending=True, kind="stable"
        )
        fig._pycsamt_before_after_line_table = line_table
        y = np.arange(len(line_table))
        line_colors = np.where(
            line_table["median_delta"] > tolerance,
            improvement_color,
            np.where(
                line_table["median_delta"] < -tolerance,
                degradation_color,
                stable_color,
            ),
        )
        line_ax.barh(
            y,
            line_table["median_delta"],
            color=line_colors,
            height=0.58,
            edgecolor="0.20",
            linewidth=0.5,
            zorder=2,
        )
        lower = line_table["median_delta"] - line_table["q1"]
        upper = line_table["q3"] - line_table["median_delta"]
        line_ax.errorbar(
            line_table["median_delta"],
            y,
            xerr=np.vstack([lower, upper]),
            fmt="none",
            ecolor="0.15",
            elinewidth=1.0,
            capsize=2.5,
            zorder=3,
        )
        for position, row in zip(y, line_table.to_dict("records")):
            line_ax.text(
                delta_limit * 0.98,
                position,
                f"+{row['improved']} / −{row['degraded']}  (n={row['n']})",
                ha="right",
                va="center",
                fontsize=6.5,
            )
        line_ax.axvline(0.0, color="0.18", lw=0.9)
        line_ax.axvspan(
            -tolerance, tolerance, color=stable_color, alpha=0.12, lw=0
        )
        line_ax.set_xlim(-delta_limit, delta_limit)
        line_ax.set_yticks(y)
        line_ax.set_yticklabels(line_table["line"])
        line_ax.set_xlabel("Median ΔCR (IQR)")
        line_ax.set_title(
            "(c)  Survey-line change",
            loc="left",
            fontsize=10,
            fontweight="semibold",
        )
    for ax in panel_axes:
        ax.grid(True, ls=":", lw=0.5, color="0.78", alpha=0.7, zorder=0)
        ax.tick_params(labelsize=7.5, direction="out", length=3)
    unmatched = len(fig._pycsamt_unmatched_before) + len(
        fig._pycsamt_unmatched_after
    )
    suffix = f"; {unmatched} unmatched" if unmatched else ""
    title = (
        f"Confidence before–after assessment "
        f"({len(comparison)} matched{suffix})"
    )
    fig.suptitle(
        title,
        fontsize=12,
        fontweight="semibold",
    )
    return fig


def plot_confidence_risk_map(
    sites: Any,
    *,
    method: str = "composite",
    coordinate_system: str = "auto",
    mode: str = "contour",
    line_labels: Any = None,
    connect: bool = True,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    risk_levels: Any = None,
    cmap: str = "YlOrRd",
    segmented_colors: bool = True,
    show_threshold_contours: bool = True,
    threshold_colors: Any = ("#d98e00", "#a50f15"),
    threshold_linewidths: Any = (1.4, 1.8),
    threshold_linestyles: Any = ("--", "-"),
    contour_labels: bool = True,
    show_stations: bool = True,
    marker_size: float = 38.0,
    station_labels: bool = False,
    station_label_step: int | None = None,
    line_names: bool = True,
    show_summary: bool = True,
    max_triangle_edge: float | None = None,
    map_aspect: str = "auto",
    figsize: tuple[float, float] = (8.2, 5.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot spatial confidence risk, defined as ``1 - confidence``.

    The confidence limits map directly to operational risk classes:
    ``risk <= 1-ci_hi`` is low, ``1-ci_hi < risk <= 1-ci_lo`` is
    moderate, and ``risk > 1-ci_lo`` is high. Contour mode is intended for
    surveys containing at least three non-collinear stations; scatter mode
    remains scientifically honest for a single profile.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    coordinate_system = str(coordinate_system).lower()
    if coordinate_system not in {"auto", "geographic", "projected"}:
        raise ValueError(
            "coordinate_system must be 'auto', 'geographic', or 'projected'."
        )
    mode = str(mode).lower()
    if mode not in {"scatter", "contour"}:
        raise ValueError("mode must be 'scatter' or 'contour'.")
    map_aspect = str(map_aspect).lower()
    if map_aspect not in {"auto", "equal", "geographic"}:
        raise ValueError(
            "map_aspect must be 'auto', 'equal', or 'geographic'."
        )

    table = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    geographic = np.isfinite(table[["longitude", "latitude"]]).all(axis=1)
    projected = np.isfinite(table[["easting", "northing"]]).all(axis=1)
    if coordinate_system == "auto":
        coordinate_system = "geographic" if geographic.any() else "projected"
    valid = geographic if coordinate_system == "geographic" else projected
    xkey, ykey = (
        ("longitude", "latitude")
        if coordinate_system == "geographic"
        else ("easting", "northing")
    )
    table = table.loc[valid & np.isfinite(table["confidence"])].copy()
    table["risk"] = 1.0 - table["confidence"]

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    ax._pycsamt_coordinate_system = coordinate_system
    ax._pycsamt_risk_table = table
    if table.empty:
        ax.text(0.5, 0.5, "no mapped stations", ha="center", va="center")
        return ax

    if line_labels is None:
        labels = pd.Series("survey", index=table.index, dtype=object)
    elif isinstance(line_labels, str):
        labels = pd.Series(line_labels, index=table.index, dtype=object)
    elif hasattr(line_labels, "get"):
        labels = table["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(valid):
            raise ValueError(
                "line_labels must contain one value per table row."
            )
        labels = pd.Series(values, index=valid.index).loc[table.index]
    table["_line"] = labels.astype(str)

    low_limit = 1.0 - ci_hi
    high_limit = 1.0 - ci_lo
    if risk_levels is None:
        risk_levels = sorted(
            {0.0, low_limit, 0.10, high_limit, 0.25, 0.50, 0.75, 1.0}
        )
    levels = np.unique(np.clip(np.asarray(risk_levels, dtype=float), 0.0, 1.0))
    if len(levels) < 2 or levels[0] > 0.0 or levels[-1] < 1.0:
        raise ValueError(
            "risk_levels must span 0 to 1 with at least two values."
        )

    from matplotlib.colors import BoundaryNorm, Normalize

    norm = (
        BoundaryNorm(levels, plt.get_cmap(cmap).N, clip=True)
        if segmented_colors
        else Normalize(0.0, 1.0, clip=True)
    )
    triangulation = None
    surface = None
    if mode == "contour":
        xy = table[[xkey, ykey]].to_numpy(dtype=float)
        if len(xy) < 3 or np.linalg.matrix_rank(xy - xy.mean(axis=0)) < 2:
            raise ValueError(
                "contour risk maps require at least three non-collinear "
                "station coordinates; use mode='scatter' for one line."
            )
        from matplotlib.tri import Triangulation

        triangulation = Triangulation(xy[:, 0], xy[:, 1])
        if max_triangle_edge is not None:
            limit = float(max_triangle_edge)
            if not np.isfinite(limit) or limit <= 0.0:
                raise ValueError(
                    "max_triangle_edge must be positive and finite."
                )
            vertices = xy[triangulation.triangles]
            edges = np.linalg.norm(
                vertices - np.roll(vertices, -1, axis=1), axis=2
            )
            triangulation.set_mask((edges > limit).any(axis=1))
            if np.all(triangulation.mask):
                raise ValueError(
                    "max_triangle_edge masks every grid triangle."
                )
        surface = ax.tricontourf(
            triangulation,
            table["risk"],
            levels=levels,
            cmap=cmap,
            norm=norm,
            extend="neither",
            zorder=0,
        )
        if show_threshold_contours:
            observed = table["risk"].to_numpy(dtype=float)
            boundaries = [
                value
                for value in (low_limit, high_limit)
                if observed.min() < value < observed.max()
            ]
            if boundaries:
                colors = list(np.atleast_1d(threshold_colors))[
                    : len(boundaries)
                ]
                widths = list(np.atleast_1d(threshold_linewidths))[
                    : len(boundaries)
                ]
                styles = list(np.atleast_1d(threshold_linestyles))[
                    : len(boundaries)
                ]
                contours = ax.tricontour(
                    triangulation,
                    table["risk"],
                    levels=boundaries,
                    colors=colors,
                    linewidths=widths,
                    linestyles=styles,
                    zorder=2,
                )
                if contour_labels:
                    fmt = {
                        low_limit: f"low/moderate  R={low_limit:.2f}",
                        high_limit: f"moderate/high  R={high_limit:.2f}",
                    }
                    ax.clabel(contours, fmt=fmt, fontsize=7, inline=True)

    if connect:
        for _, group in table.groupby("_line", sort=False):
            group = _order_map_route(group, xkey, ykey)
            ax.plot(
                group[xkey],
                group[ykey],
                color="0.30",
                lw=0.75,
                alpha=0.72,
                zorder=3,
            )
    points = None
    if show_stations or mode == "scatter":
        points = ax.scatter(
            table[xkey],
            table[ykey],
            c=table["risk"],
            cmap=cmap,
            norm=norm,
            s=marker_size,
            edgecolors="white" if mode == "contour" else "0.15",
            linewidths=0.55,
            zorder=4,
        )
    if station_labels:
        step = max(1, int(station_label_step or np.ceil(len(table) / 20)))
        for index in range(0, len(table), step):
            row = table.iloc[index]
            ax.annotate(
                str(row["station"]),
                (row[xkey], row[ykey]),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=6.2,
                zorder=5,
            )
    if line_names and table["_line"].nunique() > 1:
        for line, group in table.groupby("_line", sort=False):
            row = group.loc[group[ykey].idxmax()]
            ax.annotate(
                str(line),
                (row[xkey], row[ykey]),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                fontsize=7,
                fontweight="semibold",
                zorder=6,
            )

    mappable = surface if surface is not None else points
    colorbar = ax.figure.colorbar(mappable, ax=ax, pad=0.018, fraction=0.045)
    colorbar.set_label("Confidence risk,  R = 1 − CR", fontsize=9)
    colorbar.set_ticks(sorted({0.0, low_limit, high_limit, 0.5, 1.0}))
    colorbar.ax.tick_params(labelsize=7.5)
    for boundary in (low_limit, high_limit):
        colorbar.ax.axhline(boundary, color="0.15", lw=0.8)

    if show_summary:
        risk = table["risk"].to_numpy(dtype=float)
        counts = (
            int(np.sum(risk <= low_limit)),
            int(np.sum((risk > low_limit) & (risk <= high_limit))),
            int(np.sum(risk > high_limit)),
        )
        total = max(len(risk), 1)
        summary = (
            f"Low       {counts[0]:3d}  ({100 * counts[0] / total:4.1f}%)\n"
            f"Moderate  {counts[1]:3d}  ({100 * counts[1] / total:4.1f}%)\n"
            f"High      {counts[2]:3d}  ({100 * counts[2] / total:4.1f}%)"
        )
        ax.text(
            0.018,
            0.022,
            summary,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=7.2,
            family="monospace",
            bbox={
                "boxstyle": "round,pad=0.35",
                "facecolor": "white",
                "edgecolor": "0.40",
                "alpha": 0.90,
                "linewidth": 0.65,
            },
            zorder=7,
        )

    ax.set_xlabel(
        "Longitude (°)" if coordinate_system == "geographic" else "Easting (m)"
    )
    ax.set_ylabel(
        "Latitude (°)" if coordinate_system == "geographic" else "Northing (m)"
    )
    ax.set_title(
        f"Confidence-risk map ({method})",
        loc="left",
        fontsize=11,
        fontweight="semibold",
    )
    xspan = float(table[xkey].max() - table[xkey].min())
    yspan = float(table[ykey].max() - table[ykey].min())
    xpad = max(0.025 * xspan, 1e-12)
    ypad = max(0.035 * yspan, 1e-12)
    ax.set_xlim(
        float(table[xkey].min()) - xpad, float(table[xkey].max()) + xpad
    )
    ax.set_ylim(
        float(table[ykey].min()) - ypad, float(table[ykey].max()) + ypad
    )
    ax.grid(True, ls=":", lw=0.5, color="0.72", alpha=0.65)
    ax.tick_params(labelsize=8, direction="out", length=3)
    if map_aspect == "geographic" and coordinate_system == "geographic":
        mean_lat = float(table[ykey].mean())
        ax.set_aspect(1.0 / max(np.cos(np.deg2rad(mean_lat)), 1e-6))
    elif map_aspect == "equal":
        ax.set_aspect("equal", adjustable="box")
    else:
        ax.set_aspect("auto")
    return ax


def plot_confidence_heatmap(
    sites: Any,
    *,
    method: str = "composite",
    components: Any = None,
    line_labels: Any = None,
    station_order: str = "route",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    cmap: str = "RdYlGn",
    segmented_colors: bool = True,
    annotate: bool | str = "auto",
    annotation_fmt: str = ".2f",
    annotation_fontsize: float = 5.5,
    station_label_step: int | None = None,
    show_line_names: bool = True,
    show_line_separators: bool = True,
    missing_color: str = "#d9d9d9",
    figsize: tuple[float, float] = (13.0, 4.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot station confidence and its diagnostic components as a matrix.

    Stations are columns and metrics are rows. With ``station_order='route'``
    each survey line is ordered using coordinate-derived chainage and lines
    remain contiguous. Missing component scores are shown explicitly rather
    than being assigned a misleading confidence colour.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    labels_by_component = {
        "confidence": "Overall confidence",
        "coverage": "Data coverage",
        "uncertainty": "Tensor uncertainty",
        "offdiag": "Off-diagonal consistency",
        "diagonal": "Diagonal leakage",
        "phase": "Phase smoothness",
        "spatial": "Spatial coherence",
    }
    if components is None:
        components = list(labels_by_component)
    else:
        components = [str(value).lower() for value in components]
    unknown = sorted(set(components) - set(labels_by_component))
    if unknown:
        raise ValueError(f"unknown confidence components: {unknown}")
    if not components:
        raise ValueError("components must contain at least one metric.")
    station_order = str(station_order).lower()
    if station_order not in {"route", "input", "confidence"}:
        raise ValueError(
            "station_order must be 'route', 'input', or 'confidence'."
        )

    table = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    if line_labels is None:
        lines = pd.Series("survey", index=table.index, dtype=object)
    elif isinstance(line_labels, str):
        lines = pd.Series(line_labels, index=table.index, dtype=object)
    elif hasattr(line_labels, "get"):
        lines = table["station"].astype(str).map(line_labels).fillna("survey")
    else:
        values = list(line_labels)
        if len(values) != len(table):
            raise ValueError(
                "line_labels must contain one value per table row."
            )
        lines = pd.Series(values, index=table.index)
    table["_line"] = lines.astype(str)
    table["_input_order"] = np.arange(len(table))

    if station_order == "confidence":
        table = table.sort_values("confidence", ascending=False, kind="stable")
    elif station_order == "route":
        ordered = []
        for _, group in table.groupby("_line", sort=False):
            geographic = np.isfinite(group[["longitude", "latitude"]]).all(
                axis=1
            )
            projected = np.isfinite(group[["easting", "northing"]]).all(axis=1)
            if geographic.all():
                group = _order_map_route(group, "longitude", "latitude")
            elif projected.all():
                group = _order_map_route(group, "easting", "northing")
            elif np.isfinite(group["distance_m"]).all():
                group = group.sort_values("distance_m", kind="stable")
            ordered.append(group)
        if ordered:
            table = pd.concat(ordered, axis=0)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    if table.empty:
        ax.text(0.5, 0.5, "no confidence data", ha="center", va="center")
        return ax
    matrix = table[components].to_numpy(dtype=float).T
    masked = np.ma.masked_invalid(matrix)

    from matplotlib.colors import BoundaryNorm, Normalize

    use_cmap = plt.get_cmap(cmap).copy()
    use_cmap.set_bad(missing_color)
    boundaries = sorted({0.0, 0.25, 0.50, 0.75, ci_lo, ci_hi, 1.0})
    norm = (
        BoundaryNorm(boundaries, use_cmap.N, clip=True)
        if segmented_colors
        else Normalize(0.0, 1.0, clip=True)
    )
    image = ax.imshow(
        masked,
        cmap=use_cmap,
        norm=norm,
        interpolation="nearest",
        aspect="auto",
        origin="upper",
    )
    ax._pycsamt_confidence_table = table
    ax._pycsamt_heatmap_matrix = matrix

    n_stations = len(table)
    if annotate == "auto":
        show_annotations = n_stations * len(components) <= 180
    elif isinstance(annotate, (bool, np.bool_)):
        show_annotations = bool(annotate)
    else:
        raise ValueError("annotate must be bool or 'auto'.")
    if show_annotations:
        for row in range(matrix.shape[0]):
            for column in range(matrix.shape[1]):
                value = matrix[row, column]
                if not np.isfinite(value):
                    continue
                color = "white" if value < 0.30 or value > 0.82 else "black"
                ax.text(
                    column,
                    row,
                    format(value, annotation_fmt),
                    ha="center",
                    va="center",
                    fontsize=annotation_fontsize,
                    color=color,
                )

    step = max(1, int(station_label_step or np.ceil(n_stations / 28)))
    ticks = list(range(0, n_stations, step))
    if n_stations - 1 not in ticks:
        ticks.append(n_stations - 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels(
        table.iloc[ticks]["station"].astype(str), rotation=60, ha="right"
    )
    ax.set_yticks(np.arange(len(components)))
    ax.set_yticklabels([labels_by_component[key] for key in components])
    ax.tick_params(axis="x", labelsize=6.5, length=2)
    ax.tick_params(axis="y", labelsize=8, length=0)
    ax.set_xlabel("Stations ordered along survey lines", fontsize=9)

    if "confidence" in components and len(components) > 1:
        confidence_row = components.index("confidence")
        edge = confidence_row + 0.5
        if confidence_row == len(components) - 1:
            edge = confidence_row - 0.5
        ax.axhline(edge, color="black", lw=1.5, zorder=4)
    groups = []
    start = 0
    line_values = table["_line"].to_numpy(dtype=str)
    for index in range(1, n_stations + 1):
        if index == n_stations or line_values[index] != line_values[start]:
            groups.append((line_values[start], start, index - 1))
            start = index
    for line, first, last in groups:
        if show_line_separators and last < n_stations - 1:
            ax.axvline(last + 0.5, color="white", lw=2.4, zorder=3)
            ax.axvline(last + 0.5, color="0.20", lw=0.55, zorder=4)
        if show_line_names and (len(groups) > 1 or line != "survey"):
            ax.text(
                (first + last) / 2,
                1.015,
                line,
                transform=ax.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=7.5,
                fontweight="semibold",
                clip_on=False,
            )

    colorbar = ax.figure.colorbar(image, ax=ax, pad=0.015, fraction=0.025)
    colorbar.set_label("Confidence score", fontsize=9)
    colorbar.set_ticks(sorted({0.0, 0.5, ci_lo, ci_hi, 1.0}))
    colorbar.ax.tick_params(labelsize=7.5)
    ax.set_title(
        f"Station confidence diagnostic matrix ({method})",
        loc="left",
        fontsize=11,
        fontweight="semibold",
        pad=18 if show_line_names else 8,
    )
    return ax


def plot_confidence_coverage_curve(
    sites: Any,
    *,
    method: str = "composite",
    line_labels: Any = None,
    thresholds: Any = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    show_station_retention: bool = True,
    show_data_retention: bool = True,
    show_route_retention: bool = True,
    show_line_curves: bool = True,
    show_threshold_values: bool = True,
    station_color: str = "#2166ac",
    data_color: str = "#7b3294",
    route_color: str = "#1b7837",
    line_cmap: str = "tab10",
    figsize: tuple[float, float] = (11.5, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    axes: Any = None,
) -> plt.Figure:
    """Plot survey retention as the minimum confidence requirement rises.

    Station retention is the fraction of stations meeting a threshold. Data
    retention weights qualifying stations by their number of valid transfer-
    function rows. Route retention is the fraction of connected survey length
    whose two bounding stations both qualify; it therefore detects spatial
    fragmentation that a station count alone can hide.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    if not any(
        (show_station_retention, show_data_retention, show_route_retention)
    ):
        raise ValueError("enable at least one retention measure.")
    if thresholds is None:
        thresholds = np.linspace(0.0, 1.0, 201)
    thresholds = np.unique(np.asarray(thresholds, dtype=float))
    if (
        len(thresholds) < 2
        or not np.isfinite(thresholds).all()
        or thresholds[0] < 0.0
        or thresholds[-1] > 1.0
    ):
        raise ValueError(
            "thresholds must contain at least two finite values within 0 to 1."
        )
    table = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    table = table.loc[np.isfinite(table["confidence"])].copy()
    if line_labels is None:
        table["_line"] = "Survey"
    elif isinstance(line_labels, str):
        table["_line"] = line_labels
    elif hasattr(line_labels, "get"):
        table["_line"] = (
            table["station"].astype(str).map(line_labels).fillna("Survey")
        )
    else:
        values = list(line_labels)
        if len(values) != len(table):
            raise ValueError("line_labels must contain one value per station.")
        table["_line"] = values

    if axes is None:
        if show_line_curves:
            fig, panel_axes = plt.subplots(
                1,
                2,
                figsize=figsize,
                gridspec_kw={"width_ratios": (1.35, 1.0)},
                constrained_layout=True,
            )
        else:
            fig, survey_ax = plt.subplots(
                figsize=figsize, constrained_layout=True
            )
            panel_axes = [survey_ax]
    else:
        panel_axes = list(np.asarray(axes, dtype=object).ravel())
        required = 2 if show_line_curves else 1
        if len(panel_axes) < required:
            raise ValueError(f"axes must provide at least {required} axes.")
        panel_axes = panel_axes[:required]
        fig = panel_axes[0].figure
    survey_ax = panel_axes[0]
    line_ax = panel_axes[1] if show_line_curves else None
    if table.empty:
        survey_ax.text(
            0.5, 0.5, "no confidence data", ha="center", va="center"
        )
        return fig

    def _ordered_with_segment_lengths(group: pd.DataFrame) -> pd.DataFrame:
        geographic = np.isfinite(group[["longitude", "latitude"]]).all(axis=1)
        projected = np.isfinite(group[["easting", "northing"]]).all(axis=1)
        if geographic.all():
            ordered = _order_map_route(group, "longitude", "latitude").copy()
            xy = ordered[["longitude", "latitude"]].to_numpy(float)
            mean_lat = float(np.nanmean(xy[:, 1]))
            dx = np.diff(xy[:, 0]) * 111_320.0 * np.cos(np.deg2rad(mean_lat))
            dy = np.diff(xy[:, 1]) * 110_540.0
            lengths = np.hypot(dx, dy)
        elif projected.all():
            ordered = _order_map_route(group, "easting", "northing").copy()
            xy = ordered[["easting", "northing"]].to_numpy(float)
            lengths = np.linalg.norm(np.diff(xy, axis=0), axis=1)
        elif np.isfinite(group["distance_m"]).all():
            ordered = group.sort_values("distance_m", kind="stable").copy()
            lengths = np.abs(np.diff(ordered["distance_m"].to_numpy(float)))
        else:
            ordered = group.copy()
            lengths = np.ones(max(0, len(group) - 1), dtype=float)
        ordered.attrs["segment_lengths"] = lengths
        return ordered

    ordered_lines = [
        _ordered_with_segment_lengths(group)
        for _, group in table.groupby("_line", sort=False)
    ]
    confidence = table["confidence"].to_numpy(float)
    n_ok = table["n_ok"].to_numpy(float)
    station_retention = []
    data_retention = []
    route_retention = []
    total_valid = max(float(np.nansum(n_ok)), 1.0)
    total_route = sum(
        float(np.nansum(group.attrs["segment_lengths"]))
        for group in ordered_lines
    )
    for threshold in thresholds:
        keep = confidence >= threshold
        station_retention.append(float(np.mean(keep)))
        data_retention.append(float(np.nansum(n_ok[keep]) / total_valid))
        retained_route = 0.0
        for group in ordered_lines:
            line_keep = group["confidence"].to_numpy(float) >= threshold
            lengths = group.attrs["segment_lengths"]
            if len(lengths):
                retained_route += float(
                    np.nansum(lengths[line_keep[:-1] & line_keep[1:]])
                )
        if total_route > 0.0:
            route_retention.append(retained_route / total_route)
        else:
            route_retention.append(float(np.mean(keep)))
    curve_table = pd.DataFrame(
        {
            "threshold": thresholds,
            "station_retention": station_retention,
            "data_retention": data_retention,
            "route_retention": route_retention,
        }
    )
    fig._pycsamt_coverage_curve_table = curve_table
    specs = (
        (
            show_station_retention,
            "station_retention",
            "Stations",
            station_color,
            "-",
        ),
        (
            show_data_retention,
            "data_retention",
            "Valid TF samples",
            data_color,
            "--",
        ),
        (
            show_route_retention,
            "route_retention",
            "Connected route",
            route_color,
            "-.",
        ),
    )
    for enabled, column, label, color, linestyle in specs:
        if not enabled:
            continue
        values = curve_table[column].to_numpy(float)
        auc = float(_trapz(values, thresholds))
        survey_ax.plot(
            thresholds,
            values,
            color=color,
            lw=2.0,
            ls=linestyle,
            label=f"{label}  (AUC={auc:.2f})",
        )
    for threshold in (ci_lo, ci_hi):
        survey_ax.axvline(threshold, color="0.25", lw=0.85, ls=":")
    survey_ax.axvspan(0.0, ci_lo, color="#f4a6a6", alpha=0.12, lw=0)
    survey_ax.axvspan(ci_lo, ci_hi, color="#f6d78b", alpha=0.16, lw=0)
    survey_ax.axvspan(ci_hi, 1.0, color="#a9d8a0", alpha=0.16, lw=0)
    survey_ax.set_title(
        "(a)  Survey retention",
        loc="left",
        fontsize=10,
        fontweight="semibold",
    )
    survey_ax.set_xlabel("Minimum accepted confidence ratio")
    survey_ax.set_ylabel("Retained survey fraction")
    survey_ax.legend(frameon=False, fontsize=7.5, loc="lower left")

    if line_ax is not None:
        line_records = []
        colors = plt.get_cmap(line_cmap)(
            np.linspace(0.0, 1.0, max(len(ordered_lines), 2))
        )
        for color, group in zip(colors, ordered_lines):
            values = group["confidence"].to_numpy(float)
            retention = np.asarray(
                [np.mean(values >= threshold) for threshold in thresholds]
            )
            line = str(group["_line"].iloc[0])
            auc = float(_trapz(retention, thresholds))
            line_ax.step(
                thresholds,
                retention,
                where="post",
                color=color,
                lw=1.55,
                label=f"{line}  (AUC={auc:.2f})",
            )
            record = {"line": line, "auc": auc, "n_stations": len(values)}
            for threshold in (ci_lo, ci_hi):
                record[f"retained_at_{threshold:.2f}"] = float(
                    np.mean(values >= threshold)
                )
                if show_threshold_values:
                    line_ax.scatter(
                        [threshold],
                        [record[f"retained_at_{threshold:.2f}"]],
                        s=17,
                        color=color,
                        edgecolor="white",
                        linewidth=0.4,
                        zorder=4,
                    )
            line_records.append(record)
        fig._pycsamt_line_coverage_table = pd.DataFrame(line_records)
        for threshold in (ci_lo, ci_hi):
            line_ax.axvline(threshold, color="0.25", lw=0.85, ls=":")
        line_ax.set_title(
            "(b)  Retention by survey line",
            loc="left",
            fontsize=10,
            fontweight="semibold",
        )
        line_ax.set_xlabel("Minimum accepted confidence ratio")
        line_ax.set_ylabel("Retained station fraction")
        line_ax.legend(
            frameon=False,
            fontsize=7,
            loc="lower left",
            ncol=2 if len(ordered_lines) > 4 else 1,
        )
    for ax in panel_axes:
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.02)
        ax.grid(True, ls=":", lw=0.5, color="0.78", alpha=0.7)
        ax.tick_params(labelsize=7.5, direction="out", length=3)
    fig.suptitle(
        f"Confidence coverage curve ({method})",
        fontsize=12,
        fontweight="semibold",
    )
    return fig


def plot_confidence_distribution(
    sites: Any,
    *,
    method: str = "composite",
    metric: str = "confidence",
    line_labels: Any = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    bins: int | Any = 16,
    density: bool = True,
    show_rug: bool = True,
    show_violin: bool = True,
    show_points: bool = True,
    presence_color: str = "#377eb8",
    composite_color: str = "#e6550d",
    metric_color: str = "#5e3c99",
    class_colors: tuple[str, str, str] = (
        "#f4a6a6",
        "#f6d78b",
        "#a9d8a0",
    ),
    figsize: tuple[float, float] = (12.0, 4.3),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    axes: Any = None,
) -> plt.Figure:
    """Summarize confidence distributions, exceedance, and line spread.

    The three panels show a histogram with a lightweight Gaussian density,
    empirical cumulative distributions, and line-wise violin/box summaries.
    ``method='both'`` compares matched presence and composite scores. For a
    diagnostic component, select one method and pass its column as ``metric``.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    method = str(method).lower()
    if method not in {"presence", "composite", "both"}:
        raise ValueError("method must be 'presence', 'composite', or 'both'.")
    valid_metrics = {
        "confidence": "Confidence ratio",
        "coverage": "Data coverage",
        "uncertainty": "Tensor uncertainty",
        "offdiag": "Off-diagonal consistency",
        "diagonal": "Diagonal leakage",
        "phase": "Phase smoothness",
        "spatial": "Spatial coherence",
    }
    metric = str(metric).lower()
    if metric not in valid_metrics:
        raise ValueError(f"unknown confidence metric: {metric!r}")
    if method == "both" and metric != "confidence":
        raise ValueError("method='both' is available only for confidence.")

    methods = ("presence", "composite") if method == "both" else (method,)
    frames = []
    for name in methods:
        frame = station_confidence_table(
            sites,
            method=name,
            recursive=recursive,
            on_dup=on_dup,
            strict=strict,
            verbose=verbose,
            api=False,
        ).copy()
        frame["_method"] = name
        frames.append(frame)
    table = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if line_labels is None:
        table["_line"] = "Survey"
    elif isinstance(line_labels, str):
        table["_line"] = line_labels
    elif hasattr(line_labels, "get"):
        table["_line"] = (
            table["station"].astype(str).map(line_labels).fillna("Survey")
        )
    else:
        values = list(line_labels)
        expected = len(frames[0]) if frames else 0
        if len(values) != expected:
            raise ValueError("line_labels must contain one value per station.")
        mapping = dict(zip(frames[0]["station"].astype(str), values))
        table["_line"] = table["station"].astype(str).map(mapping)
    table = table.loc[np.isfinite(table[metric])].copy()

    if axes is None:
        fig, panel_axes = plt.subplots(
            1,
            3,
            figsize=figsize,
            gridspec_kw={"width_ratios": (1.15, 1.0, 1.35)},
            constrained_layout=True,
        )
    else:
        panel_axes = list(np.asarray(axes, dtype=object).ravel())
        if len(panel_axes) < 3:
            raise ValueError("axes must provide three plotting axes.")
        panel_axes = panel_axes[:3]
        fig = panel_axes[0].figure
    hist_ax, ecdf_ax, violin_ax = panel_axes
    if table.empty:
        hist_ax.text(0.5, 0.5, "no finite scores", ha="center", va="center")
        return fig

    palette = {
        "presence": presence_color,
        "composite": composite_color,
    }
    if method != "both":
        palette[method] = (
            metric_color if metric != "confidence" else composite_color
        )
    class_spans = (
        (0.0, ci_lo, class_colors[0]),
        (ci_lo, ci_hi, class_colors[1]),
        (ci_hi, 1.0, class_colors[2]),
    )
    for ax in (hist_ax, ecdf_ax):
        for lower, upper, color in class_spans:
            ax.axvspan(lower, upper, color=color, alpha=0.20, lw=0, zorder=0)
        for threshold in (ci_lo, ci_hi):
            ax.axvline(threshold, color="0.25", lw=0.8, ls="--", zorder=1)

    if np.isscalar(bins):
        histogram_bins = np.linspace(0.0, 1.0, int(bins) + 1)
    else:
        histogram_bins = np.asarray(bins, dtype=float)
    if len(histogram_bins) < 2:
        raise ValueError("bins must define at least one interval.")
    xgrid = np.linspace(0.0, 1.0, 300)
    for name in methods:
        values = table.loc[table["_method"] == name, metric].to_numpy(float)
        color = palette[name]
        hist_ax.hist(
            values,
            bins=histogram_bins,
            density=density,
            histtype="stepfilled",
            color=color,
            edgecolor=color,
            linewidth=0.8,
            alpha=0.28 if method == "both" else 0.55,
            label=name.capitalize(),
            zorder=2,
        )
        if density and len(values) > 1:
            spread = float(np.nanstd(values, ddof=1))
            bandwidth = max(0.025, 1.06 * spread * len(values) ** (-0.2))
            kernels = np.exp(
                -0.5 * ((xgrid[:, None] - values[None, :]) / bandwidth) ** 2
            )
            kde = kernels.mean(axis=1) / (bandwidth * np.sqrt(2.0 * np.pi))
            hist_ax.plot(xgrid, kde, color=color, lw=1.8, zorder=3)
        if show_rug:
            hist_ax.plot(
                values,
                np.full(len(values), -0.015),
                "|",
                color=color,
                ms=4,
                alpha=0.55,
                transform=hist_ax.get_xaxis_transform(),
                clip_on=False,
            )
        ordered = np.sort(values)
        probability = np.arange(1, len(ordered) + 1) / len(ordered)
        ecdf_ax.step(
            ordered,
            probability,
            where="post",
            color=color,
            lw=1.8,
            label=name.capitalize(),
        )
        median = float(np.nanmedian(values))
        hist_ax.axvline(median, color=color, lw=1.2, ls=":", zorder=3)

    hist_ax.set_title(
        "(a)  Score distribution", loc="left", fontweight="semibold"
    )
    hist_ax.set_ylabel("Probability density" if density else "Station count")
    ecdf_ax.set_title(
        "(b)  Cumulative fraction", loc="left", fontweight="semibold"
    )
    ecdf_ax.set_ylabel("Fraction of stations ≤ score")
    ecdf_ax.set_ylim(0.0, 1.02)
    ecdf_ax.set_yticks(np.linspace(0.0, 1.0, 6))
    if method == "both":
        hist_ax.legend(frameon=False, fontsize=8, loc="upper left")
        ecdf_ax.legend(frameon=False, fontsize=8, loc="lower right")

    line_order = list(dict.fromkeys(table["_line"].astype(str)))
    offsets = (
        np.linspace(-0.18, 0.18, len(methods))
        if len(methods) > 1
        else np.asarray([0.0])
    )
    rng = np.random.default_rng(0)
    for method_index, name in enumerate(methods):
        for line_index, line in enumerate(line_order):
            values = table.loc[
                (table["_method"] == name)
                & (table["_line"].astype(str) == line),
                metric,
            ].to_numpy(float)
            if not len(values):
                continue
            position = line_index + 1 + offsets[method_index]
            if show_violin and len(values) > 1:
                violin = violin_ax.violinplot(
                    [values],
                    positions=[position],
                    widths=0.30 if len(methods) > 1 else 0.58,
                    showextrema=False,
                )
                for body in violin["bodies"]:
                    body.set_facecolor(palette[name])
                    body.set_edgecolor("0.20")
                    body.set_linewidth(0.6)
                    body.set_alpha(0.55)
            q1, median, q3 = np.nanpercentile(values, [25, 50, 75])
            violin_ax.vlines(position, q1, q3, color="0.15", lw=2.0, zorder=4)
            violin_ax.scatter(
                [position],
                [median],
                s=18,
                color="white",
                edgecolor="0.15",
                zorder=5,
            )
            if show_points:
                jitter = rng.uniform(-0.035, 0.035, len(values))
                violin_ax.scatter(
                    np.full(len(values), position) + jitter,
                    values,
                    s=6,
                    color=palette[name],
                    alpha=0.35,
                    linewidths=0,
                    zorder=3,
                )
    violin_ax.axhspan(0.0, ci_lo, color=class_colors[0], alpha=0.12, lw=0)
    violin_ax.axhspan(ci_lo, ci_hi, color=class_colors[1], alpha=0.16, lw=0)
    violin_ax.axhspan(ci_hi, 1.0, color=class_colors[2], alpha=0.16, lw=0)
    for threshold in (ci_lo, ci_hi):
        violin_ax.axhline(threshold, color="0.25", lw=0.8, ls="--")
    violin_ax.set_xticks(np.arange(1, len(line_order) + 1))
    violin_ax.set_xticklabels(line_order, rotation=35, ha="right")
    violin_ax.set_ylim(0.0, 1.02)
    violin_ax.set_ylabel(valid_metrics[metric])
    violin_ax.set_title(
        "(c)  Survey-line spread", loc="left", fontweight="semibold"
    )
    if method == "both":
        from matplotlib.lines import Line2D

        handles = [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=palette[name],
                markeredgecolor="0.20",
                markersize=5,
                label=name.capitalize(),
            )
            for name in methods
        ]
        violin_ax.legend(
            handles=handles,
            frameon=False,
            fontsize=7.5,
            loc="lower right",
        )

    for ax in panel_axes:
        ax.grid(True, ls=":", lw=0.5, color="0.78", alpha=0.7)
        ax.tick_params(labelsize=7.5, direction="out", length=3)
    for ax in (hist_ax, ecdf_ax):
        ax.set_xlim(0.0, 1.0)
        ax.set_xlabel(valid_metrics[metric])
    method_title = "presence vs composite" if method == "both" else method
    fig.suptitle(
        f"Confidence distribution ({method_title})",
        fontsize=12,
        fontweight="semibold",
    )
    return fig


def plot_confidence_rank(
    sites: Any,
    *,
    method: str = "composite",
    metric: str = "confidence",
    line_labels: Any = None,
    order: str = "worst",
    max_stations: int | None = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    cmap: str = "RdYlGn",
    connect: bool = True,
    annotate_stations: bool | str = "auto",
    station_label_step: int | None = None,
    show_values: bool = False,
    show_line_panel: bool = True,
    line_statistic: str = "median",
    show_iqr: bool = True,
    marker_size: float = 36.0,
    figsize: tuple[float, float] = (11.5, 4.8),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    axes: Any = None,
) -> plt.Figure:
    """Rank station scores and survey-line performance.

    Station ranks use the exact selected metric. The line panel ranks either
    the median (default, robust to isolated failures) or mean, and optionally
    shows the interquartile range. Rank 1 follows ``order``: the lowest score
    for ``'worst'`` and the highest score for ``'best'``.
    """
    if not (0.0 <= ci_lo <= ci_hi <= 1.0):
        raise ValueError(
            "confidence thresholds must satisfy 0 <= ci_lo <= ci_hi <= 1."
        )
    metrics = {
        "confidence": "Confidence ratio",
        "coverage": "Data coverage",
        "uncertainty": "Tensor uncertainty",
        "offdiag": "Off-diagonal consistency",
        "diagonal": "Diagonal leakage",
        "phase": "Phase smoothness",
        "spatial": "Spatial coherence",
    }
    metric = str(metric).lower()
    if metric not in metrics:
        raise ValueError(f"unknown confidence metric: {metric!r}")
    order = str(order).lower()
    if order not in {"worst", "best"}:
        raise ValueError("order must be 'worst' or 'best'.")
    line_statistic = str(line_statistic).lower()
    if line_statistic not in {"median", "mean"}:
        raise ValueError("line_statistic must be 'median' or 'mean'.")

    table = station_confidence_table(
        sites,
        method=method,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        api=False,
    ).copy()
    if line_labels is None:
        table["_line"] = "Survey"
    elif isinstance(line_labels, str):
        table["_line"] = line_labels
    elif hasattr(line_labels, "get"):
        table["_line"] = (
            table["station"].astype(str).map(line_labels).fillna("Survey")
        )
    else:
        values = list(line_labels)
        if len(values) != len(table):
            raise ValueError("line_labels must contain one value per station.")
        table["_line"] = values
    table = table.loc[np.isfinite(table[metric])].copy()
    ascending = order == "worst"
    ranked = table.sort_values(
        metric, ascending=ascending, kind="stable"
    ).copy()
    if max_stations is not None:
        max_stations = int(max_stations)
        if max_stations < 1:
            raise ValueError("max_stations must be at least 1.")
        ranked = ranked.iloc[:max_stations].copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)

    if axes is None:
        if show_line_panel:
            fig, panel_axes = plt.subplots(
                1,
                2,
                figsize=figsize,
                gridspec_kw={"width_ratios": (2.25, 1.0)},
                constrained_layout=True,
            )
        else:
            fig, station_ax = plt.subplots(
                figsize=figsize, constrained_layout=True
            )
            panel_axes = [station_ax]
    else:
        panel_axes = list(np.asarray(axes, dtype=object).ravel())
        required = 2 if show_line_panel else 1
        if len(panel_axes) < required:
            raise ValueError(f"axes must provide at least {required} axes.")
        panel_axes = panel_axes[:required]
        fig = panel_axes[0].figure
    station_ax = panel_axes[0]
    line_ax = panel_axes[1] if show_line_panel else None
    fig._pycsamt_rank_table = ranked
    if ranked.empty:
        station_ax.text(0.5, 0.5, "no finite scores", ha="center", va="center")
        return fig

    from matplotlib.colors import BoundaryNorm

    boundaries = sorted({0.0, 0.25, 0.50, 0.75, ci_lo, ci_hi, 1.0})
    norm = BoundaryNorm(boundaries, plt.get_cmap(cmap).N, clip=True)
    if connect:
        station_ax.plot(
            ranked["rank"],
            ranked[metric],
            color="0.48",
            lw=0.8,
            zorder=1,
        )
    points = station_ax.scatter(
        ranked["rank"],
        ranked[metric],
        c=ranked[metric],
        cmap=cmap,
        norm=norm,
        s=marker_size,
        edgecolor="0.15",
        linewidth=0.5,
        zorder=3,
    )
    for lower, upper, color in (
        (0.0, ci_lo, "#f4a6a6"),
        (ci_lo, ci_hi, "#f6d78b"),
        (ci_hi, 1.0, "#a9d8a0"),
    ):
        station_ax.axhspan(
            lower, upper, color=color, alpha=0.16, lw=0, zorder=0
        )
    for threshold, label in (
        (ci_lo, f"CR = {ci_lo:.2f}"),
        (ci_hi, f"CR = {ci_hi:.2f}"),
    ):
        station_ax.axhline(threshold, color="0.25", lw=0.85, ls="--", zorder=2)
        station_ax.text(
            0.995,
            threshold,
            label,
            transform=station_ax.get_yaxis_transform(),
            ha="right",
            va="bottom",
            fontsize=6.5,
            color="0.25",
        )
    if annotate_stations == "auto":
        annotate = len(ranked) <= 35
    elif isinstance(annotate_stations, (bool, np.bool_)):
        annotate = bool(annotate_stations)
    else:
        raise ValueError("annotate_stations must be bool or 'auto'.")
    if annotate:
        step = max(1, int(station_label_step or 1))
        indices = list(range(0, len(ranked), step))
        if len(ranked) - 1 not in indices:
            indices.append(len(ranked) - 1)
        for index in indices:
            row = ranked.iloc[index]
            label = str(row["station"])
            if show_values:
                label += f"  {row[metric]:.2f}"
            station_ax.annotate(
                label,
                (row["rank"], row[metric]),
                xytext=(0, 5),
                textcoords="offset points",
                ha="center",
                va="bottom",
                rotation=70 if len(ranked) > 15 else 45,
                fontsize=6,
                clip_on=True,
            )
    station_ax.set_xlim(0.25, len(ranked) + 0.75)
    station_ax.set_ylim(0.0, 1.025)
    station_ax.set_xlabel(f"Station rank (1 = {order})")
    station_ax.set_ylabel(metrics[metric])
    station_ax.set_title(
        "(a)  Station ranking", loc="left", fontsize=10, fontweight="semibold"
    )
    colorbar = fig.colorbar(points, ax=station_ax, pad=0.012, fraction=0.035)
    colorbar.set_label(metrics[metric], fontsize=8)
    colorbar.set_ticks(sorted({0.0, 0.5, ci_lo, ci_hi, 1.0}))
    colorbar.ax.tick_params(labelsize=7)

    if line_ax is not None:
        records = []
        for line, group in table.groupby("_line", sort=False):
            values = group[metric].to_numpy(float)
            q1, median, q3 = np.nanpercentile(values, [25, 50, 75])
            statistic = (
                median if line_statistic == "median" else np.nanmean(values)
            )
            records.append(
                {
                    "line": str(line),
                    "score": float(statistic),
                    "q1": float(q1),
                    "q3": float(q3),
                    "n": len(values),
                }
            )
        line_table = pd.DataFrame(records).sort_values(
            "score", ascending=ascending, kind="stable"
        )
        fig._pycsamt_line_rank_table = line_table
        positions = np.arange(len(line_table))
        colors = plt.get_cmap(cmap)(norm(line_table["score"].to_numpy(float)))
        line_ax.barh(
            positions,
            line_table["score"],
            color=colors,
            edgecolor="0.20",
            linewidth=0.55,
            height=0.62,
            zorder=2,
        )
        if show_iqr:
            lower = line_table["score"] - line_table["q1"]
            upper = line_table["q3"] - line_table["score"]
            lower = np.maximum(lower, 0.0)
            upper = np.maximum(upper, 0.0)
            line_ax.errorbar(
                line_table["score"],
                positions,
                xerr=np.vstack([lower, upper]),
                fmt="none",
                ecolor="0.12",
                elinewidth=1.1,
                capsize=2.5,
                zorder=3,
            )
        for position, row in zip(positions, line_table.to_dict("records")):
            line_ax.text(
                min(row["score"] + 0.018, 0.98),
                position,
                f"{row['score']:.2f}  (n={row['n']})",
                va="center",
                ha="left" if row["score"] < 0.90 else "right",
                fontsize=6.8,
                color="0.15",
                zorder=4,
            )
        for threshold in (ci_lo, ci_hi):
            line_ax.axvline(threshold, color="0.25", lw=0.8, ls="--", zorder=1)
        line_ax.set_yticks(positions)
        line_ax.set_yticklabels(line_table["line"])
        line_ax.invert_yaxis()
        line_ax.set_xlim(0.0, 1.0)
        line_ax.set_xlabel(f"Line {line_statistic}")
        line_ax.set_title(
            "(b)  Survey-line ranking",
            loc="left",
            fontsize=10,
            fontweight="semibold",
        )
    for ax in panel_axes:
        ax.grid(True, ls=":", lw=0.5, color="0.78", alpha=0.7, zorder=0)
        ax.tick_params(labelsize=7.5, direction="out", length=3)
    fig.suptitle(
        f"Confidence ranking ({method}; {metrics[metric].lower()})",
        fontsize=12,
        fontweight="semibold",
    )
    return fig


# -------------------- confidence profile (Kouadio et al. 2024 Fig. 3) --- #


def plot_confidence_profile(
    sites: Any,
    *,
    method: str = "presence",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    shade_recoverable: bool = True,
    shade_mode: str = "score",
    annotate_low: bool = True,
    annotate_low_step: int | None = None,
    station_labels: bool = True,
    station_label_step: int | None = None,
    show_errorbars: bool = True,
    smart_ylim: bool = True,
    ylim: tuple[float, float] | None = None,
    weights: dict[str, float] | None = None,
    spacing_m: float = 200.0,
    force_spacing: bool = False,
    figsize: tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """
    Profile confidence-ratio (CR) scatter plot along the survey line.

    Reproduces the Fig. 3 style from Kouadio et al. (2024): one dot per
    station coloured green (CR >= ``ci_hi``), pink
    (``ci_lo`` <= CR < ``ci_hi``), or red (CR < ``ci_lo``), with
    dashed threshold lines.

    With ``method="presence"``, CR is the fraction of frequencies with a
    valid finite Z tensor. With ``method="composite"``, CR combines
    coverage, tensor uncertainty, off-diagonal consistency, diagonal
    leakage, phase smoothness, and neighbor coherence.

    Parameters
    ----------
    sites : path, EDI-like, Sites, or iterable
        Input sites.
    ci_hi : float
        Upper CR threshold (default 0.95, "safe", green).
    ci_lo : float
        Lower CR threshold (default 0.85, "recoverable", pink).
    shade_recoverable : bool
        If ``True``, draw an interval cue for stations below ``ci_hi``.
    shade_mode : {"score", "full", "none"}
        ``"score"`` draws compact vertical intervals tied to each station
        point. ``"full"`` preserves the older full-height station shading.
        ``"none"`` disables station interval shading.
    annotate_low : bool
        If ``True``, draw a rotated station-name label above each point
        below ``ci_lo``. Set ``False`` to turn these off entirely -- e.g.
        when ``station_labels`` (the top-axis station ticks) already
        identifies every station and the per-point labels would just
        duplicate it.
    annotate_low_step : int or None
        Gap between labeled low-confidence points, analogous to
        ``station_label_step`` but applied only to the (typically much
        smaller) subset of points below ``ci_lo``. ``None`` auto-thins
        once there are more than 18 low points, the same threshold used
        for the top axis, so a survey where most stations are flagged
        doesn't end up with every single one labeled. ``1`` forces every
        low point to be labeled regardless of count.
    station_label_step : int or None
        Gap between visible station labels on the top axis. ``None`` chooses
        a readable spacing automatically while keeping all station tick marks.
    show_errorbars : bool
        If ``True``, draw the station-level confidence uncertainty returned by
        :func:`station_confidence_table`.
    smart_ylim : bool
        If ``True``, zoom the lower y-limit when every station confidence is
        above ``ci_lo`` so small departures from the safe threshold remain
        visible.
    ylim : tuple of float or None
        Explicit y-axis limits. Overrides ``smart_ylim`` when provided.
    spacing_m : float
        The x-axis is real inter-station distance projected along the
        survey line (from EDI east/north, or lat/lon as a fallback)
        whenever at least two stations carry usable coordinates.
        ``spacing_m`` is only used as a uniform fallback for stations
        without coordinates, or for the whole line when none have any.
    force_spacing : bool
        If ``True``, skip coordinate lookup entirely and lay every
        station out at uniform ``spacing_m`` steps -- e.g. when the
        available coordinates are known to be unreliable and a
        user-supplied spacing should be trusted instead.
    figsize : tuple
        Figure size when a new figure is created.
    recursive, on_dup, strict, verbose
        Passed to :func:`ensure_sites`.
    ax : matplotlib.axes.Axes or None
        Axes to draw on; created if *None*.

    Returns
    -------
    ax : matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    shade_mode = str(shade_mode).lower()
    if shade_mode not in {"score", "full", "none"}:
        msg = "shade_mode must be 'score', 'full', or 'none'."
        raise ValueError(msg)

    tb = station_confidence_table(
        sites,
        method=method,
        weights=weights,
        spacing_m=spacing_m,
        force_spacing=force_spacing,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if tb.empty:
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        ax.set_xlabel("Distance along profile (m)")
        ax.set_ylabel("Confidence ratio")
        return ax

    xs = tb["distance_m"].to_numpy(dtype=float)
    ys = tb["confidence"].to_numpy(dtype=float)
    yerr = tb.get(
        "confidence_err",
        pd.Series(np.nan, index=tb.index),
    ).to_numpy(dtype=float)
    names = tb["station"].astype(str).tolist()
    colors = np.full(len(tb), "#d62728", dtype=object)
    colors[ys >= ci_lo] = "#ff99c8"
    colors[ys >= ci_hi] = "#20b455"

    finite_xs = xs[np.isfinite(xs)]
    if finite_xs.size > 1:
        step_width = float(np.nanmedian(np.diff(np.sort(finite_xs))))
    else:
        step_width = float(spacing_m)
    bar_width = max(step_width * 0.18, 1.0)

    if shade_recoverable and shade_mode == "full":
        order = np.argsort(xs)
        xs_ordered = xs[order]
        for idx, xpos in zip(order, xs_ordered):
            if ci_lo <= ys[idx] < ci_hi:
                if xs_ordered.size > 1:
                    diffs = np.diff(xs_ordered)
                    step = float(np.nanmedian(diffs))
                else:
                    step = spacing_m
                ax.axvspan(
                    xpos - 0.35 * step,
                    xpos + 0.35 * step,
                    ymin=0.0,
                    ymax=1.0,
                    color="#f3a6c9",
                    alpha=0.35,
                    lw=0,
                    zorder=0,
                )
    elif shade_recoverable and shade_mode == "score":
        for x, y in zip(xs, ys):
            if not np.isfinite(x) or not np.isfinite(y):
                continue
            if y >= ci_hi:
                continue
            if y >= ci_lo:
                ax.bar(
                    x,
                    y - ci_lo,
                    bottom=ci_lo,
                    width=bar_width,
                    color="#f3a6c9",
                    alpha=0.45,
                    lw=0,
                    zorder=1,
                )
                ax.bar(
                    x,
                    ci_hi - y,
                    bottom=y,
                    width=bar_width,
                    color="#8fd19e",
                    alpha=0.35,
                    lw=0,
                    zorder=1,
                )
            else:
                ax.bar(
                    x,
                    ci_lo - y,
                    bottom=y,
                    width=bar_width,
                    color="#d62728",
                    alpha=0.30,
                    lw=0,
                    zorder=1,
                )
                ax.bar(
                    x,
                    ci_hi - ci_lo,
                    bottom=ci_lo,
                    width=bar_width,
                    color="#f3a6c9",
                    alpha=0.35,
                    lw=0,
                    zorder=1,
                )

    if len(xs):
        ax.plot(xs, ys, color="black", lw=1.5, zorder=2)
        if show_errorbars and np.isfinite(yerr).any():
            ax.errorbar(
                xs,
                ys,
                yerr=np.clip(yerr, 0.0, 0.5),
                fmt="none",
                ecolor="0.25",
                elinewidth=0.8,
                capsize=2.5,
                alpha=0.65,
                zorder=2,
            )
        ax.scatter(
            xs,
            ys,
            c=colors,
            s=64,
            zorder=3,
            edgecolors="black",
            linewidths=1.0,
        )
    if annotate_low:
        low_idx = np.flatnonzero(ys < ci_lo)
        if low_idx.size:
            if annotate_low_step is None:
                low_step = (
                    max(1, int(np.ceil(low_idx.size / 12)))
                    if low_idx.size > 18
                    else 1
                )
            else:
                low_step = max(1, int(annotate_low_step))
            keep = low_idx[::low_step]
            if low_idx[-1] not in keep:
                keep = np.r_[keep, low_idx[-1]]
            for i in keep:
                ax.text(
                    xs[i],
                    max(ys[i] + 0.04, 0.04),
                    names[i],
                    ha="center",
                    va="bottom",
                    rotation=90,
                    fontsize=7,
                )

    ax.axhline(
        ci_hi,
        ls="--",
        color="black",
        lw=1.1,
        alpha=0.85,
    )
    ax.axhline(
        ci_lo,
        ls="--",
        color="black",
        lw=1.1,
        alpha=0.85,
    )
    handles = [
        plt.Line2D(
            [],
            [],
            marker="o",
            ls="",
            mfc="#20b455",
            mec="black",
            label=f"Conf. >= {ci_hi:.2f}",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            ls="",
            mfc="#ff99c8",
            mec="black",
            label=f"{ci_lo:.2f} <= Conf. < {ci_hi:.2f}",
        ),
        plt.Line2D(
            [],
            [],
            marker="o",
            ls="",
            mfc="#8b0026",
            mec="black",
            label=f"Conf. < {ci_lo:.2f}",
        ),
    ]
    if station_labels:
        top = ax.secondary_xaxis("top")
        top.set_xticks(xs, minor=True)
        top.tick_params(which="minor", length=3)
        if station_label_step is None:
            if len(xs) > 18:
                step = max(1, int(np.ceil(len(xs) / 12)))
            else:
                step = 1
        else:
            step = max(1, int(station_label_step))
        idx = np.arange(0, len(xs), step, dtype=int)
        if len(xs) and len(xs) - 1 not in idx:
            idx = np.r_[idx, len(xs) - 1]
        top.set_xticks(xs[idx])
        top.set_xticklabels(
            [names[i] for i in idx],
            rotation=90,
            fontsize=7,
        )
        top.tick_params(which="major", length=5)
        top.set_xlabel("Station")
    if ylim is not None:
        ax.set_ylim(*ylim)
    elif smart_ylim and np.nanmin(ys) >= ci_lo:
        low = max(0.0, min(ci_lo - 0.05, np.nanmin(ys) - 0.05))
        ax.set_ylim(low, 1.03)
    else:
        low = min(0.0, np.nanmin(ys) - 0.05)
        ax.set_ylim(max(-0.03, low), 1.08)
    ticks = sorted({0.0, ci_lo, ci_hi, 1.0})
    ticks = [
        tick for tick in ticks if ax.get_ylim()[0] <= tick <= ax.get_ylim()[1]
    ]
    if ticks:
        ax.set_yticks(ticks)
    ax.set_xlabel("Distance along profile (m)")
    ax.set_ylabel("Confidence ratio")
    ax.legend(handles=handles, fontsize=8, loc="lower left")
    title = "Station confidence"
    if method != "presence":
        title += f" ({method})"
    ax.set_title(title, fontsize=10)
    ax.grid(True, ls=":", alpha=0.4)
    return ax


def plot_frequency_confidence_psection(
    sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    metric: str = "confidence",
    cmap: str = "RdYlGn",
    section: str | SectionStyle = "dynamic",
    figsize: tuple[float, float] | None = None,
    station_label_step: int | None = None,
    station_preset: str = "pseudosection",
    station_style: StationAxisStyle | None = None,
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot frequency confidence as a station-period pseudo-section."""
    section_style = _resolve_section_style(section)
    # "down" triggers invert_yaxis() so short T (high freq, shallow) is at TOP.
    section_style.axis.y_direction = "down"
    tb = frequency_confidence_table(
        sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        spacing_m=spacing_m,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if tb.empty:
        if ax is None:
            _, ax = plt.subplots(
                figsize=figsize or section_style.figsize_for(),
            )
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return ax
    if metric not in tb.columns:
        msg = f"metric {metric!r} is not available in the confidence table."
        raise ValueError(msg)

    stations = tb.drop_duplicates("station").sort_values("station_index")
    station_names = stations["station"].astype(str).tolist()
    yvals = np.sort(tb["log10_period"].dropna().unique())
    if ax is None:
        _, ax = plt.subplots(
            figsize=figsize
            or section_style.figsize_for(
                n_stations=len(station_names),
                n_y=yvals.size,
                labels=station_names,
                colorbar=True,
            ),
        )
    matrix = np.full((yvals.size, len(station_names)), np.nan, dtype=float)
    for j, station in enumerate(station_names):
        sub = tb[tb["station"] == station]
        lookup = {
            float(row.log10_period): float(row[metric])
            for _, row in sub.iterrows()
            if np.isfinite(row.log10_period)
        }
        for i, yval in enumerate(yvals):
            matrix[i, j] = lookup.get(float(yval), np.nan)

    im = ax.imshow(
        matrix,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
        cmap=cmap,
        vmin=0.0,
        vmax=1.0,
        extent=(-0.5, len(station_names) - 0.5, yvals.min(), yvals.max()),
    )
    ticks = np.arange(len(station_names))
    style = station_style or PYCSAMT_STATION_RENDERING.style_for(
        station_preset or section_style.station_preset,
    )
    if station_label_step is not None:
        style = copy.copy(style)
        style.every = int(station_label_step)
    style.apply(
        ax,
        ticks,
        station_names,
        xlim=(-0.5, len(station_names) - 0.5),
    )
    section_style.apply_axis(
        ax,
        xlabel="Station",
        ylabel=r"$\log_{10}T$ (s)",
        title=f"Frequency confidence ({method})",
    )
    section_style.add_colorbar(
        im,
        ax,
        label=metric.replace("_", " ").title(),
    )
    return ax


def plot_station_confidence_spectrum(
    sites: Any,
    *,
    station: str | None = None,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    figsize: tuple[float, float] = (7.0, 4.0),
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot confidence components versus period for one station."""
    tb = frequency_confidence_table(
        sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        spacing_m=spacing_m,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if tb.empty:
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return ax
    if station is None:
        station = str(tb["station"].iloc[0])
    sub = tb[tb["station"].astype(str) == str(station)].sort_values(
        "log10_period",
    )
    if sub.empty:
        msg = f"station {station!r} is not present in the confidence table."
        raise ValueError(msg)
    x = sub["log10_period"].to_numpy(dtype=float)
    y = sub["confidence"].to_numpy(dtype=float)
    yerr = sub["confidence_err"].to_numpy(dtype=float)
    ax.fill_between(
        x,
        ci_lo,
        ci_hi,
        color="#f3a6c9",
        alpha=0.20,
        zorder=0,
    )
    ax.axhline(ci_hi, color="black", ls="--", lw=1.0)
    ax.axhline(ci_lo, color="black", ls="--", lw=1.0)
    ax.plot(x, y, color="black", lw=1.4, label="confidence")
    if np.isfinite(yerr).any():
        ax.errorbar(
            x,
            y,
            yerr=np.clip(yerr, 0.0, 0.5),
            fmt="none",
            ecolor="0.25",
            elinewidth=0.8,
            capsize=2,
            alpha=0.65,
        )
    colors = np.full(y.size, "#d62728", dtype=object)
    colors[y >= ci_lo] = "#ff99c8"
    colors[y >= ci_hi] = "#20b455"
    ax.scatter(x, y, c=colors, edgecolors="black", s=42, zorder=3)
    for key, color in (
        ("coverage", "#4e79a7"),
        ("offdiag", "#f28e2b"),
        ("diagonal", "#e15759"),
        ("phase", "#76b7b2"),
        ("spatial", "#59a14f"),
    ):
        vals = sub[key].to_numpy(dtype=float)
        if np.isfinite(vals).any():
            ax.plot(x, vals, lw=0.9, alpha=0.70, color=color, label=key)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel(r"$\log_{10}T$ (s)")
    ax.set_ylabel("Confidence")
    ax.set_title(f"{station} frequency confidence", fontsize=10)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(fontsize=7, ncol=2)
    return ax


def _confidence_panel_background(
    ax: plt.Axes, ci_hi: float, ci_lo: float
) -> None:
    """Draw confidence threshold bands for one dashboard axis."""
    ax.axhspan(0.0, ci_lo, color="#d62728", alpha=0.06, lw=0)
    ax.axhspan(ci_lo, ci_hi, color="#f3a6c9", alpha=0.10, lw=0)
    ax.axhspan(ci_hi, 1.0, color="#8fd19e", alpha=0.08, lw=0)
    ax.axhline(ci_hi, color="black", ls="--", lw=0.8, alpha=0.75)
    ax.axhline(ci_lo, color="black", ls="--", lw=0.8, alpha=0.75)


def _confidence_panel_line(
    ax: plt.Axes,
    x: np.ndarray,
    y: np.ndarray,
    *,
    color: str,
    label: str,
    ci_hi: float,
    ci_lo: float,
    yerr: np.ndarray | None = None,
) -> None:
    """Plot one dashboard line with threshold colouring."""
    _confidence_panel_background(ax, ci_hi, ci_lo)
    ax.plot(x, y, color=color, lw=1.35, label=label)
    if yerr is not None and np.isfinite(yerr).any():
        ax.errorbar(
            x,
            y,
            yerr=np.clip(yerr, 0.0, 0.5),
            fmt="none",
            ecolor="0.25",
            elinewidth=0.75,
            capsize=2,
            alpha=0.60,
        )
    marker_colors = np.full(y.size, "#d62728", dtype=object)
    marker_colors[y >= ci_lo] = "#ff99c8"
    marker_colors[y >= ci_hi] = "#20b455"
    ax.scatter(
        x,
        y,
        c=marker_colors,
        edgecolors="black",
        linewidths=0.55,
        s=24,
        zorder=3,
    )
    ax.set_ylim(-0.03, 1.05)
    ax.grid(True, ls=":", alpha=0.35)


def plot_station_confidence_dashboard(
    sites: Any,
    *,
    station: str | None = None,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    axes=None,
    figsize: tuple[float, float] = (10.5, 6.0),
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
) -> plt.Figure:
    """Plot a 2-by-3 confidence dashboard for one station.

    The dashboard separates the final confidence score from the diagnostic
    components used to build it, avoiding the visual crowding of a single
    overlay axis.
    """
    tb = frequency_confidence_table(
        sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        spacing_m=spacing_m,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    axes_given = _axes_list(axes, 6) if axes is not None else None
    if axes_given is None:
        fig, axes_grid = plt.subplots(
            2,
            3,
            figsize=figsize,
            sharex=True,
            sharey=True,
            constrained_layout=True,
        )
        flat_axes = axes_grid.ravel()
    else:
        flat_axes = np.asarray(axes_given, dtype=object)
        fig = flat_axes[0].figure
    if tb.empty:
        flat_axes[0].text(0.5, 0.5, "no stations", ha="center", va="center")
        return fig
    if station is None:
        station = str(tb["station"].iloc[0])
    sub = tb[tb["station"].astype(str) == str(station)].sort_values(
        "log10_period",
    )
    if sub.empty:
        msg = f"station {station!r} is not present in the confidence table."
        raise ValueError(msg)
    x = sub["log10_period"].to_numpy(dtype=float)
    panel_specs = [
        (
            "Overall confidence",
            "confidence",
            "black",
            sub["confidence_err"].to_numpy(dtype=float),
        ),
        ("Data coverage", "coverage", "#4e79a7", None),
        ("Tensor uncertainty", "uncertainty", "#9c755f", None),
        ("Offdiag consistency", "offdiag", "#f28e2b", None),
        ("Diagonal leakage", "diagonal", "#e15759", None),
        ("Phase + spatial coherence", None, "#76b7b2", None),
    ]
    for ax, (title, key, color, yerr) in zip(flat_axes, panel_specs):
        if key is None:
            _confidence_panel_background(ax, ci_hi, ci_lo)
            for sub_key, sub_color in (
                ("phase", "#76b7b2"),
                ("spatial", "#59a14f"),
            ):
                y = sub[sub_key].to_numpy(dtype=float)
                if np.isfinite(y).any():
                    ax.plot(x, y, color=sub_color, lw=1.25, label=sub_key)
                    ax.scatter(
                        x,
                        y,
                        color=sub_color,
                        edgecolors="black",
                        linewidths=0.45,
                        s=20,
                        zorder=3,
                    )
            ax.legend(fontsize=7, loc="lower left")
            ax.set_ylim(-0.03, 1.05)
            ax.grid(True, ls=":", alpha=0.35)
        else:
            y = sub[key].to_numpy(dtype=float)
            if np.isfinite(y).any():
                _confidence_panel_line(
                    ax,
                    x,
                    y,
                    color=color,
                    label=key,
                    ci_hi=ci_hi,
                    ci_lo=ci_lo,
                    yerr=yerr,
                )
            else:
                _confidence_panel_background(ax, ci_hi, ci_lo)
                ax.text(
                    0.5,
                    0.5,
                    "not available",
                    ha="center",
                    va="center",
                    transform=ax.transAxes,
                    color="0.35",
                )
                ax.set_ylim(-0.03, 1.05)
                ax.grid(True, ls=":", alpha=0.35)
        ax.set_title(title, fontsize=9)
    axes_grid = np.asarray(flat_axes, dtype=object).reshape(2, 3)
    for ax in axes_grid[:, 0]:
        ax.set_ylabel("Confidence")
    for ax in axes_grid[-1, :]:
        ax.set_xlabel(r"$\log_{10}T$ (s)")
    fig.suptitle(
        f"{station} frequency-confidence dashboard ({method})",
        fontsize=11,
    )
    return fig


def plot_confidence_band_summary(
    sites: Any,
    *,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    figsize: tuple[float, float] = (8.0, 4.0),
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    """Plot line-wide confidence statistics for each period sample."""
    tb = frequency_confidence_table(
        sites,
        method=method,
        ci_hi=ci_hi,
        ci_lo=ci_lo,
        spacing_m=spacing_m,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if tb.empty:
        ax.text(0.5, 0.5, "no stations", ha="center", va="center")
        return ax
    summary = (
        tb.groupby("log10_period")["confidence"]
        .agg(["median", "mean"])
        .reset_index()
        .sort_values("log10_period")
    )
    bands = (
        tb.assign(
            safe=tb["confidence"] >= ci_hi,
            recoverable=(tb["confidence"] >= ci_lo)
            & (tb["confidence"] < ci_hi),
            reject=tb["confidence"] < ci_lo,
        )
        .groupby("log10_period")[["safe", "recoverable", "reject"]]
        .mean()
        .reset_index()
        .sort_values("log10_period")
    )
    x = summary["log10_period"].to_numpy(dtype=float)
    ax.plot(
        x,
        summary["median"].to_numpy(dtype=float),
        color="black",
        lw=1.5,
        label="median confidence",
    )
    ax.plot(
        x,
        summary["mean"].to_numpy(dtype=float),
        color="0.35",
        lw=1.0,
        ls="--",
        label="mean confidence",
    )
    ax.fill_between(
        x,
        0.0,
        bands["reject"].to_numpy(dtype=float),
        color="#d62728",
        alpha=0.25,
        label="rejected fraction",
    )
    ax.fill_between(
        x,
        bands["reject"].to_numpy(dtype=float),
        (
            bands["reject"].to_numpy(dtype=float)
            + bands["recoverable"].to_numpy(dtype=float)
        ),
        color="#f3a6c9",
        alpha=0.30,
        label="recoverable fraction",
    )
    ax.axhline(ci_hi, color="black", ls="--", lw=1.0)
    ax.axhline(ci_lo, color="black", ls="--", lw=1.0)
    ax.set_ylim(-0.03, 1.05)
    ax.set_xlabel(r"$\log_{10}T$ (s)")
    ax.set_ylabel("Confidence / station fraction")
    ax.set_title(f"Period-band confidence summary ({method})", fontsize=10)
    ax.grid(True, ls=":", alpha=0.4)
    ax.legend(fontsize=7)
    return ax


# ----------------------------- coverage plot ----------------------------- #


def plot_coverage_psection(
    sites: Any,
    *,
    metric: str = "presence",  # presence|snr|offdiag
    alpha_by: str = "none",  # none|snr
    section: str | SectionStyle = "dynamic",
    figsize: tuple[float, float] | None = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    section_style = _resolve_section_style(section)
    # "down" triggers invert_yaxis() so short T (high freq, shallow) is at TOP.
    section_style.axis.y_direction = "down"
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    sts: list[str] = []
    Ys: list[np.ndarray] = []
    Ms: list[np.ndarray] = []
    As: list[np.ndarray] = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        sts.append(st)
        lp = np.log10(np.maximum(1.0 / fr, 1e-9))
        Ys.append(lp)
        if metric == "offdiag":
            M = _offdiag_logmag(z)
        elif metric == "snr":
            ze = Z[3] if isinstance(Z, tuple) else None
            if ze is None:
                _Zobj = getattr(ed, "Z", None) or getattr(
                    getattr(ed, "edi", None), "Z", None
                )
                ze = getattr(_Zobj, "z_err", None)
            M = _snr_rows(z, ze)
        else:
            M = _row_ok_z(z).astype(float)
        Ms.append(M.astype(float))
        if alpha_by == "snr":
            _Zobj = getattr(ed, "Z", None) or getattr(
                getattr(ed, "edi", None), "Z", None
            )
            ze = getattr(_Zobj, "z_err", None)
            A = _snr_rows(z, ze)
        else:
            A = np.ones_like(M, dtype=float)
        As.append(A)
    if not sts:
        if ax is None:
            _, ax = plt.subplots(
                figsize=figsize or section_style.figsize_for(),
            )
        ax.text(0.5, 0.5, "no data", ha="center", va="center")
        return ax
    yall = np.unique(np.concatenate(Ys))
    nx = len(sts)
    Zm = np.zeros((yall.size, nx, 4), dtype=float)
    v = []
    a = []
    for j, (lp, m, al) in enumerate(zip(Ys, Ms, As)):
        i = np.searchsorted(yall, lp)
        i = np.clip(i, 0, yall.size - 1)
        vv = np.nan_to_num(m, nan=np.nan)
        # RGBA alpha must be in [0, 1]; alpha_by="snr" feeds raw SNR
        # ratios (routinely > 1), so clip rather than pass them straight
        # through (imshow silently clips anyway, with a warning).
        aa = np.clip(np.nan_to_num(al, nan=0.0), 0.0, 1.0)
        v.append(vv)
        a.append(aa)
        # map metric to color
        if metric == "presence":
            col = (0.20, 0.60, 0.20)
            Zm[i, j, :3] = col
            Zm[i, j, 3] = aa
        else:
            # use viridis for metric
            pass
    if metric != "presence":
        V = np.concatenate(v)
        V = V[np.isfinite(V)]
        v0 = np.nanpercentile(V, 5) if V.size else 0.0
        v1 = np.nanpercentile(V, 95) if V.size else 1.0
        for j, (lp, m, al) in enumerate(zip(Ys, Ms, As)):
            i = np.searchsorted(yall, lp)
            i = np.clip(i, 0, yall.size - 1)
            sc = (m - v0) / (v1 - v0 + 1e-12)
            sc = np.clip(sc, 0.0, 1.0)
            rgb = plt.cm.viridis(sc)
            Zm[i, j, :3] = rgb[:, :3]
            Zm[i, j, 3] = np.clip(np.nan_to_num(al, nan=0.0), 0.0, 1.0)
    if ax is None:
        _, ax = plt.subplots(
            figsize=figsize
            or section_style.figsize_for(
                n_stations=len(sts),
                n_y=yall.size,
                labels=sts,
                colorbar=False,
            ),
        )
    ax.imshow(
        Zm,
        aspect="auto",
        origin="lower",
        interpolation="nearest",
    )
    section_style.apply_axis(
        ax,
        xlabel="Station",
        ylabel=r"$\log_{10}T$ (s)",
    )
    section_style.apply_stations(
        ax,
        np.arange(nx),
        sts,
        xlim=(-0.5, nx - 0.5),
    )
    yt, yl = _y_ticks(yall, 8)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl)
    return ax


# ----------------------------- SNR histogram ----------------------------- #


def plot_snr_hist(
    sites: Any,
    *,
    bins: int = 40,
    figsize: tuple[float, float] = (7.2, 3.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    vals: list[float] = []
    for _, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        if isinstance(Z, tuple) and len(Z) == 4:
            _, z, fr, ze = Z
        else:
            _Zobj = getattr(ed, "Z", None) or getattr(
                getattr(ed, "edi", None), "Z", None
            )
            ze = getattr(_Zobj, "z_err", None)
        snr = _snr_rows(z, ze)
        vals.extend(list(snr))
    v = np.array(vals, dtype=float)
    v = v[np.isfinite(v)]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if v.size == 0:
        ax.text(
            0.5,
            0.5,
            "SNR histogram requires impedance\n"
            "error data (z_err not available)",
            ha="center",
            va="center",
            fontsize=9,
            color="#888888",
            transform=ax.transAxes,
        )
    else:
        ax.hist(v, bins=int(max(8, bins)))
    ax.set_xlabel("row SNR (|Z|/σ)")
    ax.set_ylabel("count")
    ax.set_title("SNR Histogram")
    return ax


# ----------------------------- quicklook -------------------------------- #


def plot_qc_quicklook(
    sites: Any,
    *,
    axes=None,
    figsize: tuple[float, float] = (10.0, 8.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    axes_given = _axes_list(axes, 3) if axes is not None else None
    if axes_given is None:
        fig = plt.figure(figsize=figsize)
        gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.25)
        ax1 = fig.add_subplot(gs[0, :])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[1, 1])
    else:
        ax1, ax2, ax3 = axes_given
        fig = ax1.figure
    plot_coverage_psection(
        sites,
        metric="presence",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax1,
    )
    plot_coverage_psection(
        sites,
        metric="snr",
        alpha_by="snr",
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax2,
    )
    plot_snr_hist(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
        ax=ax3,
    )
    return fig


# ---------------------- z-block helper (errors) ------------------------- #


def _zblk(ed: Any, need_err: bool = False):
    try:
        return _get_z_block(ed, with_errors=need_err)
    except TypeError:
        try:
            return _get_z_block(ed, with_errors=need_err)
        except TypeError:
            return _get_z_block(ed)


# ----------------------- rho_a + error propagation ---------------------- #


def _rhoa_xy_yx(
    z: np.ndarray, fr: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    c = 0.2 / (fr + 1e-24)
    rxy = c * (np.abs(z[:, 0, 1]) ** 2)
    ryx = c * (np.abs(z[:, 1, 0]) ** 2)
    return rxy, ryx


def _rhoa_ci(
    z: np.ndarray,
    ze: np.ndarray | None,
    fr: np.ndarray,
    *,
    comp: str = "xy",  # xy|yx
    pcts: tuple[float, ...] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    seed: int | None = 0,
) -> np.ndarray:
    a, b = (0, 1) if comp == "xy" else (1, 0)
    zz = z[:, a, b]
    if ze is None:
        c = 0.2 / (fr + 1e-24)
        m = c * (np.abs(zz) ** 2)
        P = [np.zeros_like(m) for _ in pcts]
        return np.vstack([m] + P).T
    ee = ze[:, a, b]
    g = np.isfinite(zz) & np.isfinite(ee)
    if not np.any(g):
        m = np.full(z.shape[0], np.nan, dtype=float)
        P = [np.full_like(m, np.nan) for _ in pcts]
        return np.vstack([m] + P).T
    rng = np.random.default_rng(seed)
    nf = zz.size
    n = int(max(16, n_draws))
    # complex Gaussian, σ equals |ze|
    E = (
        rng.standard_normal((n, nf)) + 1j * rng.standard_normal((n, nf))
    ) / np.sqrt(2.0)
    E = E * ee[None, :]
    Zs = zz[None, :] + E
    c = 0.2 / (fr + 1e-24)
    R = c[None, :] * (np.abs(Zs) ** 2)
    M = np.nanmedian(R, axis=0)
    Q = [np.nanpercentile(R, q, axis=0) for q in pcts]
    return np.vstack([M] + Q).T


def _shade_band(
    ax: plt.Axes,
    x: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    *,
    alpha: float = 0.25,
    color: str = "C0",
):
    ax.fill_between(x, lo, hi, alpha=alpha, color=color)


# ----------------------- 19) Consistency fan chart ---------------------- #


def plot_consistency_fan(
    sites: Any,
    *,
    station: str | None = None,
    other: Any | None = None,  # optional comparison Sites
    comps: tuple[str, str] = ("xy", "yx"),
    pcts: tuple[float, float, float] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    figsize: tuple[float, float] = (8.6, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    ed_map = {}
    for i, ed in enumerate(_iter_items(S)):
        ed_map[_name(ed, i)] = ed
    if not ed_map:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no sites", ha="center", va="center")
        return ax
    if station is None:
        station = sorted(ed_map.keys())[0]
    ed = ed_map.get(station, None)
    if ed is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "station not found", ha="center", va="center")
        return ax
    out = _zblk(ed, need_err=True)
    if len(out) == 4:
        Z, z, fr, ze = out
    else:
        Z, z, fr = out[:3]
        ze = None
    if Z is None:
        if ax is None:
            _, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "no Z", ha="center", va="center")
        return ax
    per = 1.0 / fr
    x = per
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    ax.set_xscale("log")
    cols = {"xy": "C0", "yx": "C2"}
    for c in comps:
        CI = _rhoa_ci(z, ze, fr, comp=c, pcts=pcts, n_draws=n_draws)
        med = CI[:, 0]
        lo = np.minimum(CI[:, 1], CI[:, 2])
        hi = np.maximum(CI[:, 1], CI[:, 2])
        _shade_band(ax, x, lo, hi, color=cols[c], alpha=0.20)
        ax.plot(x, med, "-", lw=2.0, color=cols[c], label=f"ρa_{c}")
    if other is not None:
        So = ensure_sites(other, recursive=False, strict=False)
        # overlay only medians (dashed)
        for i2, edo in enumerate(_iter_items(So)):
            if _name(edo, i2) != station:
                continue
            Z2, z2, fr2 = _zblk(edo)[:3]
            if Z2 is None:
                break
            x2 = 1.0 / fr2
            rxy2, ryx2 = _rhoa_xy_yx(z2, fr2)
            if "xy" in comps:
                ax.plot(
                    x2, rxy2, "--", lw=1.2, color=cols["xy"], label="after xy"
                )
            if "yx" in comps:
                ax.plot(
                    x2, ryx2, "--", lw=1.2, color=cols["yx"], label="after yx"
                )
            break
    ax.set_xlabel("Period (s)")
    ax.set_ylabel("ρa (Ω·m)")
    ax.grid(True, alpha=0.25, which="both")
    ax.set_title(str(station))
    ax.legend(ncol=2, fontsize=8)
    return ax


# ---------------------- 20) XY–YX crossover map ------------------------- #


def plot_xyyx_crossover_map(
    sites: Any,
    *,
    figsize: tuple[float, float] = (9.0, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: plt.Axes | None = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    Xs, Ys, labels = [], [], []
    sts = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i)
        sts.append(st)
        Z, z, fr = _zblk(ed)[:3]
        if Z is None:
            continue
        rxy, ryx = _rhoa_xy_yx(z, fr)
        d = rxy - ryx
        p = 1.0 / fr
        if d.size < 2:
            continue
        s = np.sign(d)
        sc = s[:-1] * s[1:] <= 0.0
        idx = np.where(sc)[0]
        if idx.size == 0:
            continue
        lp = np.log10(np.maximum(p, 1e-9))
        for k in idx:
            w1 = np.abs(d[k])
            w2 = np.abs(d[k + 1])
            t = w1 / (w1 + w2 + 1e-24)
            y = (1.0 - t) * lp[k] + t * lp[k + 1]
            Xs.append(i)
            Ys.append(y)
            labels.append(st)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if not Xs:
        ax.text(0.5, 0.5, "no crossovers", ha="center", va="center")
        return ax
    ax.scatter(Xs, Ys, s=16, c="crimson", alpha=0.8)
    ax.set_ylabel(LOG10_PERIOD_LABEL)
    PYCSAMT_STATION_RENDERING.apply(
        ax,
        np.arange(len(sts), dtype=float),
        sts,
        preset="pseudosection",
        xlim=(-0.5, len(sts) - 0.5),
    )
    # y ticks from Ys
    yall = np.array(Ys, dtype=float)
    yt, yl = _y_ticks(yall, 8)
    ax.set_yticks(yt)
    lo, hi = float(np.nanmin(yall)), float(np.nanmax(yall))
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.05 * (hi - lo))
    if not ax.yaxis_inverted():
        ax.invert_yaxis()
    return ax


# ---------------------- 21) Noise cone overlay -------------------------- #


def overlay_noise_cone(
    ax: plt.Axes,
    period: np.ndarray,
    lo: np.ndarray,
    hi: np.ndarray,
    *,
    color: str = "0.6",
    alpha: float = 0.18,
):
    x = period
    _shade_band(ax, x, lo, hi, color=color, alpha=alpha)


# ---------------------- 22) Spectral hole finder ------------------------ #


def overlay_spectral_holes(
    ax: plt.Axes,
    sites: Any,
    *,
    thresh_dec: float = 0.30,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
):
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    # assumes x = station index, y = log10(period)
    xmap = {}
    for i, ed in enumerate(_iter_items(S)):
        xmap[_name(ed, i)] = i
    for i, ed in enumerate(_iter_items(S)):
        Z, z, fr = _zblk(ed)[:3]
        if Z is None:
            continue
        p = 1.0 / fr
        lp = np.sort(np.log10(np.maximum(p, 1e-9)))
        if lp.size < 2:
            continue
        d = np.diff(lp)
        holes = np.where(d > float(thresh_dec))[0]
        for h in holes:
            y0, y1 = lp[h], lp[h + 1]
            r = _Rect(
                (i - 0.45, y0),
                0.90,
                (y1 - y0),
                facecolor=(0.8, 0.2, 0.2, 0.08),
                edgecolor="none",
                zorder=0.0,
            )
            ax.add_patch(r)


# Public quick-look helpers retain concise summaries for API autosummary.
overlay_noise_cone.__doc__ = (
    "Overlay lower and upper noise envelopes on an existing period axis."
)
overlay_spectral_holes.__doc__ = (
    "Highlight gaps in spectral coverage on an existing QC plot."
)
plot_consistency_fan.__doc__ = (
    "Plot cross-station response consistency as a fan diagram."
)
plot_coverage_psection.__doc__ = (
    "Plot frequency coverage and data availability as a pseudosection."
)
plot_qc_quicklook.__doc__ = (
    "Create a compact multi-panel quality-control summary for a survey."
)
plot_snr_hist.__doc__ = (
    "Plot the distribution of signal-to-noise ratios across survey data."
)
plot_xyyx_crossover_map.__doc__ = (
    "Map XY/YX crossover behaviour across stations and frequencies."
)
