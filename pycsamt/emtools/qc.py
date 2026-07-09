
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
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle as _Rect

from ..api.section import PYCSAMT_SECTION, SectionStyle
from ..api.labels import LOG10_PERIOD_LABEL
from ..api.station import (
    PYCSAMT_STATION_RENDERING,
    StationAxisStyle,
)
from ..api.view import maybe_wrap_frame

from ._core import (
    ensure_sites,
    _axes_list,
    _iter_items,
    _get_z_block,
    _get_t_block,
    _name,
    _station_positions,
)
from .tensor import build_phase_tensor_table

__all__ = [
    "build_qc_table",
    "confidence_ratio",
    "frequency_confidence_table",
    "plot_confidence_band_summary",
    "plot_confidence_profile",
    "plot_frequency_confidence_psection",
    "plot_station_confidence_dashboard",
    "plot_station_confidence_spectrum",
    "qc_flags",
    "station_confidence_table",
]

DEFAULT_CONFIDENCE_WEIGHTS: Dict[str, float] = {
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


def _snr_rows(z: np.ndarray, ze: Optional[np.ndarray]) -> np.ndarray:
    if ze is None:
        return np.full(z.shape[0], np.nan, dtype=float)
    a = np.sqrt(np.nanmean(np.abs(z) ** 2, axis=(1, 2)))
    e = np.sqrt(np.nanmean(np.abs(ze) ** 2, axis=(1, 2)))
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


def _weighted_nanmean(values: Dict[str, float],
                      weights: Dict[str, float]) -> float:
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


def _confidence_error(values: Dict[str, float],
                      n_freq: int,
                      confidence: float) -> float:
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
    scores: Dict[str, float],
    *,
    weights: Optional[Dict[str, float]] = None,
    n_freq: int = 1,
    return_error: bool = False,
) -> float | Tuple[float, float]:
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


def _relerr_score(z: np.ndarray, ze: Optional[np.ndarray],
                  threshold: float) -> float:
    """Score tensor uncertainty from median relative error."""
    if ze is None:
        return np.nan
    rel = np.abs(ze) / (np.abs(z) + 1e-24)
    if not np.isfinite(rel).any():
        return np.nan
    med = float(np.nanmedian(rel))
    return _clip01(1.0 - med / max(float(threshold), 1e-12))


def _offdiag_consistency_score(z: np.ndarray,
                               tolerance_log10: float) -> float:
    """Score similarity of ``Zxy`` and ``Zyx`` amplitudes."""
    zxy = np.abs(z[:, 0, 1])
    zyx = np.abs(z[:, 1, 0])
    ratio = np.log10((zxy + 1e-24) / (zyx + 1e-24))
    if not np.isfinite(ratio).any():
        return np.nan
    med = float(np.nanmedian(np.abs(ratio)))
    return _clip01(1.0 - med / max(float(tolerance_log10), 1e-12))


def _diagonal_leakage_score(z: np.ndarray,
                            max_fraction: float) -> float:
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


def _phase_smoothness_score(z: np.ndarray,
                            jump_tolerance_deg: float) -> float:
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


def _frequency_flags(row: pd.Series,
                     ci_hi: float,
                     ci_lo: float) -> str:
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


def _y_ticks(yall: np.ndarray, ny: int) -> Tuple[np.ndarray, List[str]]:
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
    rows: List[Dict[str, Any]] = []
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
                    np.nanpercentile(sb, 75)
                    - np.nanpercentile(sb, 25)
                )
            else:
                rec["skew_med"] = np.nan
                rec["skew_iqr"] = np.nan
        rows.append(rec)
    cols = [
        "station", "n_freq", "n_ok", "frac_ok", "n_tip",
        "n_tip_ok", "snr_med", "pmin", "pmax",
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
        if np.isfinite(r.get("skew_med", np.nan)) and \
           r["skew_med"] > max_skew_med:
            f.append("high_skew")
        flags.append(",".join(f))
    out = tb.copy()
    out["flags"] = flags
    return out


def station_confidence_table(
    sites: Any,
    *,
    method: str = "composite",
    weights: Optional[Dict[str, float]] = None,
    relerr_threshold: float = 0.20,
    offdiag_tolerance_log10: float = 0.35,
    diagonal_leakage_max: float = 0.35,
    phase_jump_tolerance_deg: float = 90.0,
    spatial_tolerance_log10: float = 0.60,
    spacing_m: float = 200.0,
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
    positions = _station_positions(items, spacing_m)
    rows: List[Dict[str, Any]] = []
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
        rows.append(
            dict(
                station=st,
                distance_m=float(positions[i]) if i < positions.size else np.nan,
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
            )
        )
    if not rows:
        df = pd.DataFrame(
            columns=[
                "station", "distance_m", "confidence", "method",
                "confidence_err", "n_freq", "n_ok", "coverage",
                "uncertainty", "offdiag", "diagonal", "phase",
                "spatial",
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
                float(spatial_scores[i])
                if i < spatial_scores.size
                else np.nan
            )
            parts = {
                key: row[key]
                for key in (
                    "coverage", "uncertainty", "offdiag",
                    "diagonal", "phase", "spatial",
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
                float(spatial_scores[i])
                if i < spatial_scores.size
                else np.nan
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
    weights: Optional[Dict[str, float]] = None,
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    relerr_threshold: float = 0.20,
    offdiag_tolerance_log10: float = 0.35,
    diagonal_leakage_max: float = 0.35,
    phase_jump_tolerance_deg: float = 90.0,
    spatial_tolerance_log10: float = 0.60,
    spacing_m: float = 200.0,
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
    positions = _station_positions(items, spacing_m)
    rows: List[Dict[str, Any]] = []
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
                    1.0 - abs(value)
                    / max(float(offdiag_tolerance_log10), 1e-12),
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
            error_parts = parts if method == "composite" else {
                "coverage": parts["coverage"],
            }
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
                    float(np.log10(1.0 / freq))
                    if freq > 0
                    else np.nan
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
        "station", "station_index", "distance_m", "frequency_hz",
        "period_s", "log10_period", "confidence", "confidence_err",
        "method", "n_components", "coverage", "uncertainty", "offdiag",
        "diagonal", "phase", "spatial", "logrho_proxy", "flags",
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
                    "coverage", "uncertainty", "offdiag",
                    "diagonal", "phase", "spatial",
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
        _frequency_flags(row, ci_hi, ci_lo)
        for _, row in table.iterrows()
    ]

    return maybe_wrap_frame(
        table,
        api=api,
        name="frequency_confidence_table",
        kind="emtools.qc.frequency_confidence",
        source=sites,
        description="Frequency-level transfer-function confidence scores.",
    )


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
    station_labels: bool = True,
    station_label_step: Optional[int] = None,
    show_errorbars: bool = True,
    smart_ylim: bool = True,
    ylim: Optional[Tuple[float, float]] = None,
    weights: Optional[Dict[str, float]] = None,
    spacing_m: float = 200.0,
    figsize: Tuple[float, float] = (9.0, 4.0),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
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
        Fallback station spacing [m] used when no coordinate metadata is
        available on the EDI objects.
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
        for x, y, name in zip(xs, ys, names):
            if y < ci_lo:
                ax.text(
                    x,
                    max(y + 0.04, 0.04),
                    name,
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
            [], [], marker="o", ls="", mfc="#20b455", mec="black",
            label=f"Conf. >= {ci_hi:.2f}",
        ),
        plt.Line2D(
            [], [], marker="o", ls="", mfc="#ff99c8", mec="black",
            label=f"{ci_lo:.2f} <= Conf. < {ci_hi:.2f}",
        ),
        plt.Line2D(
            [], [], marker="o", ls="", mfc="#8b0026", mec="black",
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
    ticks = [tick for tick in ticks if ax.get_ylim()[0] <= tick <= ax.get_ylim()[1]]
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
    figsize: Optional[Tuple[float, float]] = None,
    station_label_step: Optional[int] = None,
    station_preset: str = "pseudosection",
    station_style: Optional[StationAxisStyle] = None,
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
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
            figsize=figsize or section_style.figsize_for(
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
    station: Optional[str] = None,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    figsize: Tuple[float, float] = (7.0, 4.0),
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
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


def _confidence_panel_background(ax: plt.Axes,
                                 ci_hi: float,
                                 ci_lo: float) -> None:
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
    yerr: Optional[np.ndarray] = None,
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
    station: Optional[str] = None,
    method: str = "composite",
    ci_hi: float = DEFAULT_CI_HI,
    ci_lo: float = DEFAULT_CI_LO,
    axes=None,
    figsize: Tuple[float, float] = (10.5, 6.0),
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
    figsize: Tuple[float, float] = (8.0, 4.0),
    spacing_m: float = 200.0,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
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
    alpha_by: str = "none",    # none|snr
    section: str | SectionStyle = "dynamic",
    figsize: Optional[Tuple[float, float]] = None,
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
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
    sts: List[str] = []
    Ys: List[np.ndarray] = []
    Ms: List[np.ndarray] = []
    As: List[np.ndarray] = []
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
                    getattr(ed, "edi", None), "Z", None)
                ze = getattr(_Zobj, "z_err", None)
            M = _snr_rows(z, ze)
        else:
            M = _row_ok_z(z).astype(float)
        Ms.append(M.astype(float))
        if alpha_by == "snr":
            _Zobj = getattr(ed, "Z", None) or getattr(
                getattr(ed, "edi", None), "Z", None)
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
        v.append(vv); a.append(aa)
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
            figsize=figsize or section_style.figsize_for(
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
    figsize: Tuple[float, float] = (7.2, 3.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites,
        recursive=recursive,
        on_dup=on_dup,
        strict=strict,
        verbose=verbose,
    )
    vals: List[float] = []
    for _, ed in enumerate(_iter_items(S)):
        Z, z, fr = _get_z_block(ed)
        if Z is None:
            continue
        if isinstance(Z, tuple) and len(Z) == 4:
            _, z, fr, ze = Z
        else:
            _Zobj = getattr(ed, "Z", None) or getattr(
                getattr(ed, "edi", None), "Z", None)
            ze = getattr(_Zobj, "z_err", None)
        snr = _snr_rows(z, ze)
        vals.extend(list(snr))
    v = np.array(vals, dtype=float)
    v = v[np.isfinite(v)]
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if v.size == 0:
        ax.text(
            0.5, 0.5,
            "SNR histogram requires impedance\nerror data (z_err not available)",
            ha="center", va="center", fontsize=9, color="#888888",
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
    figsize: Tuple[float, float] = (10.0, 8.0),
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

def _rhoa_xy_yx(z: np.ndarray, fr: np.ndarray) -> Tuple[
    np.ndarray, np.ndarray
]:
    c = 0.2 / (fr + 1e-24)
    rxy = c * (np.abs(z[:, 0, 1]) ** 2)
    ryx = c * (np.abs(z[:, 1, 0]) ** 2)
    return rxy, ryx


def _rhoa_ci(
    z: np.ndarray,
    ze: Optional[np.ndarray],
    fr: np.ndarray,
    *,
    comp: str = "xy",        # xy|yx
    pcts: Tuple[float, ...] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    seed: Optional[int] = 0,
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
    E = (rng.standard_normal((n, nf))
         + 1j * rng.standard_normal((n, nf))) / np.sqrt(2.0)
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
    station: Optional[str] = None,
    other: Any | None = None,      # optional comparison Sites
    comps: Tuple[str, str] = ("xy", "yx"),
    pcts: Tuple[float, float, float] = (10.0, 50.0, 90.0),
    n_draws: int = 200,
    figsize: Tuple[float, float] = (8.6, 4.2),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
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
        ax.text(0.5, 0.5, "station not found",
                ha="center", va="center")
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
        CI = _rhoa_ci(
            z, ze, fr, comp=c, pcts=pcts, n_draws=n_draws
        )
        med = CI[:, 0]
        lo = np.minimum(CI[:, 1], CI[:, 2])
        hi = np.maximum(CI[:, 1], CI[:, 2])
        _shade_band(ax, x, lo, hi, color=cols[c], alpha=0.20)
        ax.plot(x, med, "-", lw=2.0, color=cols[c],
                label=f"ρa_{c}")
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
                ax.plot(x2, rxy2, "--", lw=1.2, color=cols["xy"],
                        label="after xy")
            if "yx" in comps:
                ax.plot(x2, ryx2, "--", lw=1.2, color=cols["yx"],
                        label="after yx")
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
    figsize: Tuple[float, float] = (9.0, 4.6),
    recursive: bool = True,
    on_dup: str = "replace",
    strict: bool = False,
    verbose: int = 0,
    ax: Optional[plt.Axes] = None,
) -> plt.Axes:
    S = ensure_sites(
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
    )
    Xs, Ys, labels = [], [], []
    sts = []
    for i, ed in enumerate(_iter_items(S)):
        st = _name(ed, i); sts.append(st)
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
            w1 = np.abs(d[k]); w2 = np.abs(d[k + 1])
            t = w1 / (w1 + w2 + 1e-24)
            y = (1.0 - t) * lp[k] + t * lp[k + 1]
            Xs.append(i)
            Ys.append(y)
            labels.append(st)
    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    if not Xs:
        ax.text(0.5, 0.5, "no crossovers",
                ha="center", va="center")
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
        sites, recursive=recursive, on_dup=on_dup,
        strict=strict, verbose=verbose,
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
