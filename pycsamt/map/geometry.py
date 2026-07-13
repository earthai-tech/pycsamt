# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""Real-geometry line-layout helpers for 3-D survey maps.

Multi-line 3-D views (fence, block, depth-slice) need two things per
survey line: an along-line ("profile distance") axis and a
line-to-line ("offset") position. Placing lines with a synthetic
``index * spacing`` stack ignores where the lines actually sit on the
ground, which staggers panels that should be parallel and evenly
spaced.

This module derives both from real station coordinates instead:
every station is projected onto a shared local (u, v) frame — ``u``
along the survey's dominant strike, ``v`` across it — so parallel
lines stay parallel and the line-to-line spacing reflects their true
geographic separation.

Pure numpy, no Plotly/Dash/pycsamt-app dependencies, so it is safe to
import from both :mod:`pycsamt.map` (the code-first / mapview 3-D
builder) and :mod:`pycsamt.app.web` (the legacy 3-D map callbacks) —
whichever data structure each caller uses (``MapData`` stations, or
plain per-line dicts) is adapted to plain id/lat/lon/line arrays
before calling in here.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

__all__ = [
    "equirect_xy",
    "fit_strike",
    "survey_uv",
    "normalize_offsets",
    "resolve_offset",
]

# Metres per degree, equirectangular approximation — adequate for the
# local extent (tens of km) of a single survey; avoids a hard pyproj
# dependency for what is otherwise a display-layout computation.
_LON_M_PER_DEG = 111_320.0
_LAT_M_PER_DEG = 110_574.0


def equirect_xy(
    lat: Sequence[float],
    lon: Sequence[float],
) -> tuple[np.ndarray, np.ndarray]:
    """Local planar metres for *lat*/*lon*, equirectangular approx.

    Centred implicitly at the mean latitude of *lat*.
    """
    lat_arr = np.asarray(lat, dtype=float)
    lon_arr = np.asarray(lon, dtype=float)
    mlat = np.deg2rad(float(np.nanmean(lat_arr)))
    x = lon_arr * _LON_M_PER_DEG * np.cos(mlat)
    y = lat_arr * _LAT_M_PER_DEG
    return x, y


def fit_strike(
    points: np.ndarray,
    groups: np.ndarray,
) -> np.ndarray | None:
    """Return the shared along-line unit direction for *points*.

    Fits each group's (line's) own best-fit axis — SVD on that
    group's centred points — then averages them, sign-aligned and
    weighted by point count. This reflects "which way each line
    actually runs" rather than the aggregate shape of the whole
    point cloud, which can pick the wrong axis when a survey's
    footprint is roughly as wide (line-to-line) as it is long
    (along-line), e.g. a short/dense grid.

    Returns ``None`` when no direction can be fit (fewer than two
    points overall, or a degenerate/singular point cloud).
    """
    directions: list[np.ndarray] = []
    weights: list[int] = []
    for name in np.unique(groups):
        mask = groups == name
        if int(mask.sum()) < 2:
            continue
        sub = points[mask]
        sub = sub - sub.mean(axis=0)
        try:
            _, singular, vt = np.linalg.svd(sub, full_matrices=False)
        except np.linalg.LinAlgError:
            continue
        if singular[0] <= 0:
            continue
        directions.append(vt[0])
        weights.append(int(mask.sum()))
    if not directions:
        if points.shape[0] < 2:
            return None
        try:
            _, singular, vt = np.linalg.svd(points, full_matrices=False)
        except np.linalg.LinAlgError:
            return None
        return vt[0] if singular[0] > 0 else None
    ref = directions[0]
    total = np.zeros(2)
    for direction, weight in zip(directions, weights):
        if np.dot(direction, ref) < 0:
            direction = -direction
        total += direction * weight
    norm = float(np.linalg.norm(total))
    return total / norm if norm > 0 else ref


def survey_uv(
    ids: Sequence[str],
    lat: Sequence[float],
    lon: Sequence[float],
    line: Sequence[str],
) -> dict[str, tuple[float, float]]:
    """Project every station onto a shared local (u, v) frame, metres.

    ``u`` is the along-strike distance, ``v`` the cross-strike
    offset, both measured from the combined centroid of every station
    passed in. *ids*, *lat*, *lon*, *line* must be the same length and
    aligned (one entry per station).

    Returns ``{}`` when there are fewer than two stations or no
    strike direction could be fit (e.g. all stations coincide).
    """
    ids = list(ids)
    if len(ids) < 2:
        return {}
    x, y = equirect_xy(lat, lon)
    pts = np.column_stack([x, y])
    centered = pts - pts.mean(axis=0)
    strike = fit_strike(centered, np.asarray(line, dtype=object))
    if strike is None:
        return {}
    perp = np.array([-strike[1], strike[0]])
    u = centered @ strike
    v = centered @ perp
    return {sid: (float(ui), float(vi)) for sid, ui, vi in zip(ids, u, v)}


def normalize_offsets(
    offsets: dict[str, float | None],
) -> dict[str, float] | None:
    """Shift *offsets* so the smallest value is zero.

    Returns ``None`` (signalling "no real geometry, fall back to a
    synthetic layout") when *offsets* is empty or any value is
    missing/non-finite — mixing real and synthetic placement in the
    same plot would look inconsistent, so it is all-or-nothing.
    """
    if not offsets:
        return None
    for value in offsets.values():
        if value is None or not np.isfinite(value):
            return None
    base = min(offsets.values())
    return {name: float(value) - base for name, value in offsets.items()}


def resolve_offset(
    name: str,
    idx: int,
    real_offsets: dict[str, float] | None,
    fallback_unit: float,
    spacing: float,
) -> float:
    """Cross-strike placement for one line, in metres.

    Uses the survey's true line-to-line spacing (scaled by
    *spacing*, e.g. a UI "line spacing" multiplier where ``1.0`` means
    true spacing) when *real_offsets* is available, otherwise falls
    back to a synthetic ``idx * spacing * fallback_unit`` stack.
    """
    if real_offsets is not None:
        return real_offsets[name] * spacing
    return idx * spacing * fallback_unit
