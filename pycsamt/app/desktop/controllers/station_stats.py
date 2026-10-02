# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
Per-station quality statistics for the main window's Station Statistics
card (Qt-free).

The card used to show "quality dots" computed as ``n_frequencies // 10``
-- a count, not a quality.  :func:`station_stats` measures it from the
data:

* **completeness** -- share of frequencies whose four impedance
  components are all finite (off-diagonals only for scalar data);
* **SNR** -- median ``|Z| / sigma`` over the off-diagonal components;
* **2-D consistency** -- share of frequencies with phase-tensor skew
  ``|beta| <= 3 deg`` (Caldwell et al., 2004);

and grades the station A-E from their mean score.  Per-frequency relative
errors colour the coverage strip, and apparent resistivity / phase curves
feed the response preview.  :func:`survey_medians` gives the same numbers
over the whole survey, for the "vs survey" comparison.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

__all__ = ["GRADE_COLOURS", "StationStats", "station_stats",
           "survey_medians", "grade_of"]

# Filled badges with >= 4.5:1 contrast on white text (as Pipeline Studio).
GRADE_COLOURS = {"A": "#2a7f3f", "B": "#1864ab", "C": "#9c4f00",
                 "D": "#c2410c", "E": "#c92a2a", "–": "#5f6b7a"}
_MU0 = 4e-7 * math.pi


@dataclass
class StationStats:
    name: str
    freq: np.ndarray = field(default_factory=lambda: np.array([]))
    rel_err: np.ndarray = field(default_factory=lambda: np.array([]))
    period: np.ndarray = field(default_factory=lambda: np.array([]))
    rho_xy: np.ndarray = field(default_factory=lambda: np.array([]))
    rho_yx: np.ndarray = field(default_factory=lambda: np.array([]))
    phi_xy: np.ndarray = field(default_factory=lambda: np.array([]))
    phi_yx: np.ndarray = field(default_factory=lambda: np.array([]))
    completeness: float = float("nan")
    snr: float = float("nan")
    skew_ok: float = float("nan")
    scalar: bool = False
    grade: str = "–"
    score: float = float("nan")

    def breakdown(self) -> list[tuple[str, str, float]]:
        """(label, value text, 0-1 score) rows for the card/tooltips."""
        return [
            ("Completeness", _pct(self.completeness),
             _nz(self.completeness)),
            ("Signal / noise", "–" if not math.isfinite(self.snr)
             else f"{self.snr:.0f}", _snr_score(self.snr)),
            ("2-D consistency", _pct(self.skew_ok), _nz(self.skew_ok)),
        ]


def _pct(v: float) -> str:
    return "–" if not math.isfinite(v) else f"{100 * v:.0f} %"


def _nz(v: float) -> float:
    return 0.0 if not math.isfinite(v) else float(v)


def _snr_score(snr: float) -> float:
    """0 at SNR 1, 1 at SNR >= 100 (log scale)."""
    if not math.isfinite(snr) or snr <= 0:
        return 0.0
    return float(np.clip(math.log10(snr) / 2.0, 0.0, 1.0))


def grade_of(score: float) -> str:
    if not math.isfinite(score):
        return "–"
    for g, lo in (("A", 0.85), ("B", 0.70), ("C", 0.55), ("D", 0.40)):
        if score >= lo:
            return g
    return "E"


def _apparent(z: np.ndarray, f: np.ndarray):
    """rho_a (ohm m) and phase (deg, first quadrant) from Z in
    [mV/km]/[nT] (EDI units: Z_SI = Z * 4e-4 * pi... via 0.2 T |Z|^2)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        rho = 0.2 / f * np.abs(z) ** 2
        phi = np.mod(np.degrees(np.angle(z)), 180.0)
    return rho, phi


def _skew_beta(z: np.ndarray) -> np.ndarray:
    """Phase-tensor skew beta (deg) per frequency, NaN if singular."""
    X, Y = z.real, z.imag
    out = np.full(z.shape[0], np.nan)
    for i in range(z.shape[0]):
        x = X[i]
        if not np.isfinite(z[i]).all():
            continue
        det = np.linalg.det(x)
        if abs(det) < 1e-30:
            continue
        P = np.linalg.solve(x, Y[i])
        tr = P[0, 0] + P[1, 1]
        out[i] = 0.5 * np.degrees(np.arctan2(P[0, 1] - P[1, 0], tr))
    return out


def station_stats(site) -> StationStats:
    """Measured quality of one :class:`~pycsamt.site.base.Site`."""
    name = str(getattr(site, "name", ""))
    try:
        f = np.asarray(site.freq, dtype=float)
        z = np.asarray(site.z, dtype=complex)
    except Exception:
        return StationStats(name)
    if f.ndim != 1 or z.ndim != 3 or z.shape[0] != f.size or f.size == 0:
        return StationStats(name)
    try:
        ze = np.asarray(site.z_err, dtype=float)
        if ze.shape != z.shape:
            ze = None
    except Exception:
        ze = None
    order = np.argsort(f)
    f, z = f[order], z[order]
    ze = ze[order] if ze is not None else None

    xy, yx = z[:, 0, 1], z[:, 1, 0]
    scalar = not np.isfinite(yx).any()
    if scalar:
        complete = np.isfinite(xy)
    else:
        complete = np.isfinite(z).all(axis=(1, 2))
    comps = [(0, 1)] if scalar else [(0, 1), (1, 0)]
    rel = np.full(f.size, np.nan)
    snr_vals = []
    if ze is not None:
        import warnings

        with np.errstate(invalid="ignore", divide="ignore"), \
                warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            r = np.nanmean(np.column_stack([
                ze[:, a, b] / np.abs(z[:, a, b]) for a, b in comps]), axis=1)
        rel = np.where(np.isfinite(r) & (r > 0), r, np.nan)
        for a, b in comps:
            with np.errstate(invalid="ignore", divide="ignore"):
                s = np.abs(z[:, a, b]) / ze[:, a, b]
            snr_vals.append(s[np.isfinite(s) & (s > 0)])
    snr_all = np.concatenate(snr_vals) if snr_vals else np.array([])
    snr = float(np.median(snr_all)) if snr_all.size else float("nan")

    if scalar:
        skew_ok = float("nan")
    else:
        beta = _skew_beta(z)
        good = np.isfinite(beta)
        skew_ok = float(np.mean(np.abs(beta[good]) <= 3.0)) \
            if good.any() else float("nan")

    rho_xy, phi_xy = _apparent(xy, f)
    rho_yx, phi_yx = _apparent(yx, f)
    parts = [float(np.mean(complete)), _snr_score(snr)]
    if math.isfinite(skew_ok):
        parts.append(skew_ok)
    if not math.isfinite(snr):
        parts = [p for i, p in enumerate(parts) if i != 1]
    score = float(np.mean(parts)) if parts else float("nan")
    return StationStats(
        name=name, freq=f, rel_err=rel, period=1.0 / f,
        rho_xy=rho_xy, rho_yx=rho_yx, phi_xy=phi_xy, phi_yx=phi_yx,
        completeness=float(np.mean(complete)), snr=snr, skew_ok=skew_ok,
        scalar=scalar, grade=grade_of(score), score=score,
    )


def survey_medians(sites) -> dict:
    """Median statistics over every station (for "vs survey")."""
    rows = []
    fmin, fmax = [], []
    for s in sites or []:
        st = station_stats(s)
        rows.append((st.freq.size, st.completeness, st.snr, st.skew_ok,
                     st.score))
        if st.freq.size:
            fmin.append(st.freq.min())
            fmax.append(st.freq.max())
    if not rows:
        return {}
    a = np.asarray(rows, dtype=float)

    def med(col):
        v = a[:, col]
        v = v[np.isfinite(v)]
        return float(np.median(v)) if v.size else float("nan")

    return {"n_stations": len(rows), "nfreq": med(0), "completeness": med(1),
            "snr": med(2), "skew_ok": med(3), "score": med(4),
            "fmin": float(min(fmin)) if fmin else float("nan"),
            "fmax": float(max(fmax)) if fmax else float("nan")}
