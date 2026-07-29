# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.metadata.frequency
===========================

Named frequency band registry for electromagnetic geophysics.

Every EM method is characterised by its operational frequency range.
:class:`FrequencyBand` captures that range, together with period bounds
and a rule-of-thumb depth-of-investigation (DOI) estimate based on the
skin-depth formula:

    δ (m) ≈ 503 √( ρ / f )     (ρ in Ω·m, f in Hz)

The pre-defined registry :data:`MT_BANDS` covers the eight principal
EM methods used in applied geophysics.

Quick start
-----------
::

    from pycsamt.metadata.frequency import (
        MT_BANDS,
        band_for_frequency,
        REGISTRY,
    )

    b = MT_BANDS["AMT"]
    print(b.f_min, b.f_max)  # 10.0  100000.0
    print(b.period_min, b.period_max)
    print(b.doi_range_m(rho=100.0))  # DOI estimate at 100 Ω·m

    # generate log-spaced frequencies inside the band
    freqs = b.logspace(n=30)

    # which bands contain 1 Hz?
    matches = band_for_frequency(1.0)

    # forward/em1d default
    lo, hi = frequency_range("AMT")
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "FrequencyBand",
    "MT_BANDS",
    "REGISTRY",
    "band_for_frequency",
    "frequency_range",
    "register_band",
]

# Skin-depth reference resistivity used for DOI estimates (Ω·m)
_DOI_REF_RHO = 100.0

# mu_0 / (2 pi) prefactor for skin depth: δ = sqrt(2 ρ / (μ0 ω))
# → δ ≈ 503.3 √(ρ / f) metres
_SKIN_DEPTH_COEFF = 503.3


# ---------------------------------------------------------------------------
# FrequencyBand
# ---------------------------------------------------------------------------


@dataclass
class FrequencyBand:
    """Specification of an EM geophysical frequency band.

    Parameters
    ----------
    name : str
        Short identifier, e.g. ``"AMT"``.
    label : str
        Human-readable name, e.g. ``"Audio-frequency MT"``.
    f_min : float
        Minimum operational frequency in Hz.
    f_max : float
        Maximum operational frequency in Hz.
    method : str
        EM acquisition method: ``"MT"``, ``"AMT"``, ``"CSAMT"``,
        ``"TEM"``, ``"CSEM"``, etc.
    doi_ref_rho : float, default 100.0
        Reference resistivity (Ω·m) used to compute DOI estimates via
        the skin-depth formula.
    notes : str
        Optional annotation.

    Computed attributes
    -------------------
    period_min, period_max : float
        Period bounds derived from ``f_max`` and ``f_min``.
    n_decades : float
        Number of frequency decades spanned by the band.

    Examples
    --------
    ::

        b = MT_BANDS["AMT"]
        print(b.f_min, b.f_max)  # 10.0  100000.0
        print(b.period_range)  # (1e-5, 0.1) s
        print(b.doi_range_m())  # DOI at 100 Ω·m
        print(b.doi_range_m(rho=10.0))  # DOI at 10 Ω·m

        # log-spaced frequencies
        freqs = b.logspace(30)

        # test membership
        assert 1000.0 in b
    """

    name: str
    label: str
    f_min: float
    f_max: float
    method: str
    doi_ref_rho: float = _DOI_REF_RHO
    notes: str = ""

    def __post_init__(self) -> None:
        if self.f_min <= 0 or self.f_max <= 0:
            raise ValueError("f_min and f_max must be > 0")
        if self.f_min > self.f_max:
            self.f_min, self.f_max = self.f_max, self.f_min

    # ------------------------------------------------------------------
    # Computed properties
    # ------------------------------------------------------------------

    @property
    def period_min(self) -> float:
        """Minimum period in seconds (= 1 / f_max)."""
        return 1.0 / self.f_max

    @property
    def period_max(self) -> float:
        """Maximum period in seconds (= 1 / f_min)."""
        return 1.0 / self.f_min

    @property
    def period_range(self) -> tuple[float, float]:
        """``(period_min, period_max)`` in seconds."""
        return (self.period_min, self.period_max)

    @property
    def frequency_range(self) -> tuple[float, float]:
        """``(f_min, f_max)`` in Hz."""
        return (self.f_min, self.f_max)

    @property
    def n_decades(self) -> float:
        """Number of frequency decades in the band."""
        return math.log10(self.f_max / self.f_min)

    # ------------------------------------------------------------------
    # DOI estimation via skin depth
    # ------------------------------------------------------------------

    def skin_depth_m(self, f: float, rho: float | None = None) -> float:
        """Skin depth δ (m) at frequency *f* and resistivity *rho*.

        Formula: δ ≈ 503.3 √(ρ / f)

        Parameters
        ----------
        f : float
            Frequency in Hz.
        rho : float, optional
            Resistivity in Ω·m; defaults to :attr:`doi_ref_rho`.
        """
        r = rho if rho is not None else self.doi_ref_rho
        return _SKIN_DEPTH_COEFF * math.sqrt(max(r, 1e-9) / max(f, 1e-12))

    def doi_range_m(
        self,
        rho: float | None = None,
    ) -> tuple[float, float]:
        """Return the *(shallow, deep)* depth-of-investigation estimate in metres.

        Shallow DOI corresponds to ``f_max``; deep DOI to ``f_min``.

        Parameters
        ----------
        rho : float, optional
            Resistivity in Ω·m; defaults to :attr:`doi_ref_rho`.
        """
        return (
            round(self.skin_depth_m(self.f_max, rho), 1),
            round(self.skin_depth_m(self.f_min, rho), 1),
        )

    # ------------------------------------------------------------------
    # Utilities
    # ------------------------------------------------------------------

    def logspace(self, n: int = 30) -> np.ndarray:
        """Return *n* log-spaced frequencies inside the band (Hz).

        Parameters
        ----------
        n : int, default 30
            Number of frequencies.
        """
        return np.logspace(
            math.log10(self.f_min),
            math.log10(self.f_max),
            num=n,
        )

    def contains(self, f: float) -> bool:
        """Return True when *f* (Hz) lies within the band."""
        return self.f_min <= f <= self.f_max

    def overlaps(self, other: FrequencyBand) -> bool:
        """Return True when *self* and *other* share any frequency."""
        return self.f_min <= other.f_max and other.f_min <= self.f_max

    def intersection(self, other: FrequencyBand) -> tuple[float, float] | None:
        """Return the (f_min, f_max) intersection with *other*, or None."""
        lo = max(self.f_min, other.f_min)
        hi = min(self.f_max, other.f_max)
        return (lo, hi) if lo <= hi else None

    def clip_frequencies(self, freqs: Any) -> np.ndarray:
        """Return only those elements of *freqs* that lie inside the band."""
        f = np.asarray(freqs, dtype=float)
        return f[(f >= self.f_min) & (f <= self.f_max)]

    def __contains__(self, f: Any) -> bool:
        try:
            return self.contains(float(f))
        except (TypeError, ValueError):
            return False

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "label": self.label,
            "f_min": self.f_min,
            "f_max": self.f_max,
            "period_min": self.period_min,
            "period_max": self.period_max,
            "method": self.method,
            "doi_ref_rho": self.doi_ref_rho,
            "doi_shallow_m": self.doi_range_m()[0],
            "doi_deep_m": self.doi_range_m()[1],
            "n_decades": round(self.n_decades, 2),
            "notes": self.notes,
        }

    def __repr__(self) -> str:
        doi = self.doi_range_m()
        return (
            f"FrequencyBand({self.name!r}  "
            f"{self.f_min:.2g}–{self.f_max:.2g} Hz  "
            f"DOI≈{doi[0]:.0f}–{doi[1]:.0f} m  "
            f"[{self.method}])"
        )


# ---------------------------------------------------------------------------
# Pre-defined band registry
# ---------------------------------------------------------------------------

#: Dictionary of all pre-defined frequency bands.
MT_BANDS: dict[str, FrequencyBand] = {
    # ── Passive MT / broadband ──────────────────────────────────────────
    "LMT": FrequencyBand(
        name="LMT",
        label="Long-period MT",
        f_min=1e-5,
        f_max=1e-1,
        method="MT",
        notes="Geomagnetic depth sounding; crustal and mantle targets",
    ),
    "MT": FrequencyBand(
        name="MT",
        label="Magnetotelluric (standard)",
        f_min=1e-4,
        f_max=1e3,
        method="MT",
        notes="Broadband passive-source MT; crustal to sedimentary targets",
    ),
    "BBMT": FrequencyBand(
        name="BBMT",
        label="Broadband MT (wide-range)",
        f_min=1e-4,
        f_max=1e4,
        method="MT",
        notes="Combined MT+AMT coverage; typical modern broadband system",
    ),
    "AMT": FrequencyBand(
        name="AMT",
        label="Audio-frequency MT",
        f_min=1e1,
        f_max=1e5,
        method="AMT",
        notes="High-frequency MT; shallow targets (tens of metres)",
    ),
    "LAMT": FrequencyBand(
        name="LAMT",
        label="Long-period AMT",
        f_min=1.0,
        f_max=1e4,
        method="AMT",
        notes="Intermediate AMT; 1–10 000 Hz",
    ),
    # ── Controlled-source ──────────────────────────────────────────────
    "CSAMT": FrequencyBand(
        name="CSAMT",
        label="Controlled-source AMT",
        f_min=1e0,
        f_max=1e4,
        method="CSAMT",
        notes=(
            "Active source; far-field approximation valid above "
            "a few source–receiver offsets"
        ),
    ),
    "CSEM": FrequencyBand(
        name="CSEM",
        label="Marine CSEM",
        f_min=1e-2,
        f_max=1e0,
        method="CSEM",
        notes="Low-frequency towed-source marine survey; HC reservoir targets",
    ),
    # ── Time-domain ────────────────────────────────────────────────────
    "TEM": FrequencyBand(
        name="TEM",
        label="Time-domain EM (TEM / TDEM)",
        f_min=1e-4,
        f_max=1e2,
        method="TEM",
        notes=(
            "Equivalent frequency range for a central-loop TEM system; "
            "ground/airborne; conductive target detection"
        ),
    ),
}

#: Global mutable registry — start with MT_BANDS; users can add custom bands.
REGISTRY: dict[str, FrequencyBand] = dict(MT_BANDS)


# ---------------------------------------------------------------------------
# Public helpers
# ---------------------------------------------------------------------------


def register_band(band: FrequencyBand) -> None:
    """Add or replace a :class:`FrequencyBand` in the global registry.

    Parameters
    ----------
    band : FrequencyBand
        The band to register under ``band.name``.
    """
    REGISTRY[band.name] = band


def band_for_frequency(
    f: float,
    registry: dict[str, FrequencyBand] | None = None,
) -> list[FrequencyBand]:
    """Return all bands in the registry that contain *f* Hz.

    Parameters
    ----------
    f : float
        Query frequency in Hz.
    registry : dict, optional
        Registry to search; defaults to :data:`REGISTRY`.

    Returns
    -------
    list of FrequencyBand
        Sorted from narrowest to widest (fewest to most decades).
    """
    reg = registry if registry is not None else REGISTRY
    matches = [b for b in reg.values() if b.contains(f)]
    return sorted(matches, key=lambda b: b.n_decades)


def frequency_range(
    method: str,
    registry: dict[str, FrequencyBand] | None = None,
) -> tuple[float, float]:
    """Return ``(f_min, f_max)`` for a named method or band.

    Parameters
    ----------
    method : str
        Band name (e.g. ``"AMT"``) **or** method string (e.g.
        ``"CSAMT"``).  Case-insensitive.  When multiple bands share the
        same method, the union range is returned.
    registry : dict, optional
        Registry to search; defaults to :data:`REGISTRY`.

    Returns
    -------
    tuple[float, float]
        ``(f_min, f_max)`` in Hz.

    Raises
    ------
    KeyError
        When *method* is not found in any band.
    """
    reg = registry if registry is not None else REGISTRY
    key = method.upper()

    # Direct name lookup first
    if key in reg:
        b = reg[key]
        return (b.f_min, b.f_max)

    # Method-level union
    matched = [b for b in reg.values() if b.method.upper() == key]
    if not matched:
        raise KeyError(
            f"No band found for method/name {method!r}.  "
            f"Available: {sorted(reg.keys())}"
        )
    return (
        min(b.f_min for b in matched),
        max(b.f_max for b in matched),
    )


def doi_estimate(
    f: float,
    rho: float = _DOI_REF_RHO,
) -> float:
    """Return the skin-depth DOI estimate in metres at frequency *f*.

    Parameters
    ----------
    f : float
        Frequency in Hz.
    rho : float, default 100.0
        Reference resistivity in Ω·m.
    """
    return _SKIN_DEPTH_COEFF * math.sqrt(max(rho, 1e-9) / max(f, 1e-12))
