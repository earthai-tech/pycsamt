# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.methods.emap
====================
EMAP — Electromagnetic Array Profiling helper used by
``pycsamt.utils.em._prepare_sensitivity`` to compute per-station,
per-frequency Bahr skewness (η) and Swift rotational invariant (μ)
matrices on a common frequency grid.
"""
from __future__ import annotations

from typing import Any, List, Optional, Tuple

import numpy as np

from ..emtools._core import _iter_items, _get_z_block, _name, ensure_sites


# ── helpers ──────────────────────────────────────────────────────────────────

def _bahr_skew(z: np.ndarray) -> np.ndarray:
    """Bahr (1988) phase-sensitive skewness η per frequency."""
    Zxx, Zxy = z[:, 0, 0], z[:, 0, 1]
    Zyx, Zyy = z[:, 1, 0], z[:, 1, 1]
    s1, s2 = Zxx + Zyy, Zxy - Zyx
    d1, d2 = Zxx - Zyy, Zxy + Zyx
    num = np.abs(s1) ** 2 + np.abs(s2) ** 2
    den = np.abs(d1) ** 2 + np.abs(d2) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        eta = np.sqrt(np.where(den > 0, num / den, np.nan))
    return eta


def _swift_skew(z: np.ndarray) -> np.ndarray:
    """Swift (1967) skew μ = |Zxx − Zyy| / |Zxy + Zyx| per frequency."""
    Zxx, Zxy = z[:, 0, 0], z[:, 0, 1]
    Zyx, Zyy = z[:, 1, 0], z[:, 1, 1]
    denom = np.abs(Zxy + Zyx)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu = np.where(denom < 1e-30, 0.0, np.abs(Zxx - Zyy) / denom)
    return mu


def _suppress_outliers(arr: np.ndarray, sigma: float = 3.0) -> np.ndarray:
    """Replace values more than *sigma* std from the median with NaN."""
    med = np.nanmedian(arr)
    std = np.nanstd(arr)
    if std == 0 or not np.isfinite(std):
        return arr
    out = arr.copy()
    out[np.abs(arr - med) > sigma * std] = np.nan
    return out


def _interp_to_grid(freqs_src: np.ndarray, vals: np.ndarray,
                    freqs_tgt: np.ndarray) -> np.ndarray:
    """Log-linear interpolation of *vals* from *freqs_src* to *freqs_tgt*."""
    if freqs_src.size == 0 or vals.size == 0:
        return np.full(freqs_tgt.size, np.nan)
    lf_src = np.log10(np.clip(freqs_src, 1e-10, None))
    lf_tgt = np.log10(np.clip(freqs_tgt, 1e-10, None))
    return np.interp(lf_tgt, lf_src, vals, left=np.nan, right=np.nan)


# ── EMAP class ────────────────────────────────────────────────────────────────

class EMAP:
    """
    Electromagnetic Array Profiling processor.

    Fits a collection of MT sites and exposes skewness matrices on a
    common logarithmically-spaced frequency grid.

    Attributes (after :meth:`fit`)
    --------------------------------
    freqs_ : ndarray, shape (n_freqs,)
        Common frequency vector (Hz), sorted descending.
    id : list of str
        Station labels in the order they were fitted.
    ediObjs_ : list
        Raw EDI/site items passed to :meth:`fit`.
    survey_name : str
        Empty string unless set externally.
    """

    def __init__(self) -> None:
        self.freqs_: Optional[np.ndarray] = None
        self.id: List[str] = []
        self.ediObjs_: List[Any] = []
        self.survey_name: str = ""
        self._skew_cols: List[np.ndarray] = []
        self._mu_cols:   List[np.ndarray] = []

    # ── fit ──────────────────────────────────────────────────────────────

    def fit(self, edis_list: Any) -> "EMAP":
        """
        Extract per-station Z tensors and build common-grid skewness arrays.

        Parameters
        ----------
        edis_list : Sites | list of EDI/Z-like objects
            Any collection accepted by :func:`~pycsamt.emtools._core._iter_items`.

        Returns
        -------
        self
        """
        # Coerce to Sites if possible so _iter_items works uniformly
        try:
            sites = ensure_sites(edis_list, recursive=True,
                                 on_dup="replace", strict=False, verbose=0)
        except Exception:
            sites = edis_list   # use as-is; _iter_items handles lists

        # Collect per-station skew/mu and build a union frequency grid
        per_station_freqs: List[np.ndarray] = []
        per_station_eta:   List[np.ndarray] = []
        per_station_mu:    List[np.ndarray] = []
        names: List[str] = []
        raw:   List[Any] = []

        for i, ed in enumerate(_iter_items(sites)):
            _, z, fr = _get_z_block(ed)
            if z is None or fr is None or fr.size == 0:
                continue
            # Sort descending (high freq first, MT convention)
            order = np.argsort(fr)[::-1]
            fr = fr[order]
            z  = z[order]

            eta = _bahr_skew(z)
            mu  = _swift_skew(z)

            per_station_freqs.append(fr)
            per_station_eta.append(eta)
            per_station_mu.append(mu)
            names.append(_name(ed, i))
            raw.append(ed)

        if not per_station_freqs:
            # Nothing usable — return an empty result so callers don't crash
            self.freqs_ = np.array([1.0])
            self.id = []
            self.ediObjs_ = []
            return self

        # Build a common frequency grid: log-spaced from global min/max
        all_freqs = np.concatenate(per_station_freqs)
        f_min = np.nanmin(all_freqs[all_freqs > 0])
        f_max = np.nanmax(all_freqs)
        # Number of grid points ≈ median number of frequencies per station
        n_pts = int(np.median([len(f) for f in per_station_freqs]))
        n_pts = max(n_pts, 4)
        common = np.logspace(np.log10(f_max), np.log10(f_min), n_pts)

        # Interpolate each station to the common grid
        self._skew_cols = [
            _interp_to_grid(fr, eta, common)
            for fr, eta in zip(per_station_freqs, per_station_eta)
        ]
        self._mu_cols = [
            _interp_to_grid(fr, mu, common)
            for fr, mu in zip(per_station_freqs, per_station_mu)
        ]

        self.freqs_   = common
        self.id       = names
        self.ediObjs_ = raw
        return self

    # ── skew ─────────────────────────────────────────────────────────────

    def skew(
        self,
        method: str = "Bahr",
        *,
        suppress_outliers: bool = True,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Return skewness matrices on the common frequency grid.

        Parameters
        ----------
        method : {'Bahr', 'Swift'}
            Primary skewness metric.  The secondary metric (the other
            one) is always returned as the second element.
        suppress_outliers : bool
            Replace values beyond 3σ from the column median with NaN.

        Returns
        -------
        primary_matrix : ndarray, shape (n_freqs, n_stations)
        secondary_matrix : ndarray, shape (n_freqs, n_stations)
        """
        if self.freqs_ is None:
            raise RuntimeError("Call fit() before skew().")

        n_f = self.freqs_.size
        n_s = len(self._skew_cols)

        if n_s == 0:
            empty = np.full((n_f, 0), np.nan)
            return empty, empty

        skew_mat = np.column_stack(self._skew_cols).reshape(n_f, n_s)
        mu_mat   = np.column_stack(self._mu_cols).reshape(n_f, n_s)

        if suppress_outliers:
            for j in range(n_s):
                skew_mat[:, j] = _suppress_outliers(skew_mat[:, j])
                mu_mat[:, j]   = _suppress_outliers(mu_mat[:, j])

        if str(method).lower().startswith("b"):
            return skew_mat, mu_mat
        return mu_mat, skew_mat
