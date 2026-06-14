# -*- coding: utf-8 -*-
# Author: LKouadio <etanoyau@gmail.com>
# License: LGPL-3.0
"""
pycsamt.metadata.quality
========================

Data-quality assessment for MT impedance sites.

:class:`DataQuality` is a first-class object that captures the
coverage, validity, and signal-to-noise ratio of every impedance
component at a single station.  :func:`assess_collection` computes
quality for an entire :class:`~pycsamt.site.base.Sites` collection
and returns a DataFrame suitable for reporting, filtering, or
upstream ML pipelines.

Quality thresholds
------------------
+----------+----------------------+---------------------------+
| Flag     | Coverage             | Meaning                   |
+==========+======================+===========================+
| GOOD     | ≥ 0.90               | suitable for inversion    |
+----------+----------------------+---------------------------+
| PARTIAL  | 0.50 – 0.90          | use with care             |
+----------+----------------------+---------------------------+
| POOR     | 0.01 – 0.50          | likely noisy / incomplete |
+----------+----------------------+---------------------------+
| MISSING  | 0.00                 | component absent          |
+----------+----------------------+---------------------------+

Quick start
-----------
::

    from pycsamt.metadata.quality import DataQuality, assess_collection

    # single site
    dq = DataQuality.from_site(site)
    print(dq.overall)          # QualityFlag.GOOD / PARTIAL / POOR / MISSING
    print(dq.summary())

    # collection
    df = assess_collection(sites)    # pandas DataFrame
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..api.view import maybe_wrap_frame

__all__ = [
    "QualityFlag",
    "ComponentQuality",
    "DataQuality",
    "assess_collection",
    "quality_dataframe",
]

# Component ordering
_Z_COMPONENTS = ("Zxx", "Zxy", "Zyx", "Zyy")
_ALL_COMPONENTS = _Z_COMPONENTS + ("Tipper",)

# Coverage thresholds
_THRESH_GOOD    = 0.90
_THRESH_PARTIAL = 0.50


# ---------------------------------------------------------------------------
# QualityFlag
# ---------------------------------------------------------------------------

class QualityFlag(str, Enum):
    """Data-quality classification for one component or a whole station.

    Attributes
    ----------
    GOOD : "good"
        Coverage ≥ 90 %; component is reliable for inversion.
    PARTIAL : "partial"
        50 % ≤ coverage < 90 %; component has significant gaps.
    POOR : "poor"
        0 % < coverage < 50 %; component is largely missing or noisy.
    MISSING : "missing"
        Zero finite values; component is absent.
    """

    GOOD    = "good"
    PARTIAL = "partial"
    POOR    = "poor"
    MISSING = "missing"

    @classmethod
    def from_coverage(cls, coverage: float) -> "QualityFlag":
        """Return the flag that matches *coverage* (0–1)."""
        if coverage <= 0.0:
            return cls.MISSING
        if coverage < _THRESH_PARTIAL:
            return cls.POOR
        if coverage < _THRESH_GOOD:
            return cls.PARTIAL
        return cls.GOOD

    @property
    def rank(self) -> int:
        """Numeric severity rank (0 = worst, 3 = best)."""
        return {"missing": 0, "poor": 1, "partial": 2, "good": 3}[self.value]

    @classmethod
    def worst(cls, flags: "List[QualityFlag]") -> "QualityFlag":
        """Return the worst (lowest-rank) flag in *flags*."""
        if not flags:
            return cls.MISSING
        return min(flags, key=lambda f: f.rank)

    @classmethod
    def best(cls, flags: "List[QualityFlag]") -> "QualityFlag":
        """Return the best (highest-rank) flag in *flags*."""
        if not flags:
            return cls.MISSING
        return max(flags, key=lambda f: f.rank)


# ---------------------------------------------------------------------------
# ComponentQuality
# ---------------------------------------------------------------------------

@dataclass
class ComponentQuality:
    """Quality metrics for a single impedance component.

    Parameters
    ----------
    name : str
        Component label: ``"Zxx"``, ``"Zxy"``, ``"Zyx"``, ``"Zyy"``,
        or ``"Tipper"``.
    coverage : float
        Fraction of finite values in the frequency range (0–1).
    n_valid : int
        Number of frequencies with finite values.
    n_total : int
        Total number of frequencies.
    flag : QualityFlag
        Auto-computed classification from *coverage*.
    snr_mean : float, optional
        Mean signal-to-noise ratio (dB) when available.
    snr_std : float, optional
        Standard deviation of the SNR (dB) when available.

    Examples
    --------
    ::

        cq = ComponentQuality.from_array("Zxy", z_arr)
        print(cq.flag)          # QualityFlag.GOOD
        print(cq.pct_str)       # "100%"
    """

    name: str
    coverage: float
    n_valid: int
    n_total: int
    flag: QualityFlag = field(init=False)
    snr_mean: Optional[float] = None
    snr_std: Optional[float] = None

    def __post_init__(self) -> None:
        self.flag = QualityFlag.from_coverage(self.coverage)

    @property
    def pct_str(self) -> str:
        """Coverage as a formatted percentage string."""
        return f"{self.coverage * 100:.0f}%"

    @classmethod
    def from_array(
        cls,
        name: str,
        arr: Any,
        snr: Optional[Any] = None,
    ) -> "ComponentQuality":
        """Compute quality from a raw array.

        Parameters
        ----------
        name : str
            Component label.
        arr : array-like or None
            Raw complex or real impedance values; shape ``(n,)`` or
            ``(n, 2, 2)`` (only the relevant component is expected).
        snr : array-like, optional
            Per-frequency SNR values in dB (same length as *arr*).
        """
        if arr is None:
            return cls(name=name, coverage=0.0, n_valid=0, n_total=0)

        a = np.asarray(arr)
        if a.size == 0:
            return cls(name=name, coverage=0.0, n_valid=0, n_total=0)

        # For complex arrays check both real and imaginary parts
        if np.iscomplexobj(a):
            finite = np.isfinite(a.real) & np.isfinite(a.imag)
        else:
            finite = np.isfinite(a.astype(float, copy=False))

        n_total = int(finite.size)
        n_valid = int(finite.sum())
        coverage = n_valid / n_total if n_total > 0 else 0.0

        snr_mean = snr_std = None
        if snr is not None:
            s = np.asarray(snr, dtype=float).ravel()
            fs = s[np.isfinite(s)]
            if fs.size:
                snr_mean = float(fs.mean())
                snr_std  = float(fs.std())

        return cls(
            name=name,
            coverage=coverage,
            n_valid=n_valid,
            n_total=n_total,
            snr_mean=snr_mean,
            snr_std=snr_std,
        )

    def to_dict(self) -> Dict[str, Any]:
        return {
            "name":     self.name,
            "coverage": round(self.coverage, 4),
            "n_valid":  self.n_valid,
            "n_total":  self.n_total,
            "flag":     self.flag.value,
            "snr_mean": round(self.snr_mean, 2) if self.snr_mean is not None else None,
            "snr_std":  round(self.snr_std, 2)  if self.snr_std  is not None else None,
        }

    def __repr__(self) -> str:
        snr = (
            f"  SNR={self.snr_mean:.1f} dB"
            if self.snr_mean is not None else ""
        )
        return (
            f"ComponentQuality({self.name!r}  {self.pct_str}"
            f"  [{self.flag.value}]{snr})"
        )


# ---------------------------------------------------------------------------
# DataQuality
# ---------------------------------------------------------------------------

@dataclass
class DataQuality:
    """Data-quality summary for a single MT station.

    Parameters
    ----------
    station : str
        Station identifier.
    n_freq : int
        Total number of frequencies.
    freq_min : float, optional
        Minimum frequency (Hz).
    freq_max : float, optional
        Maximum frequency (Hz).
    components : list of ComponentQuality
        Per-component quality records.
    overall : QualityFlag
        Worst flag across all present components (auto-computed).

    Examples
    --------
    ::

        dq = DataQuality.from_site(site)
        print(dq.overall)
        print(dq.get("Zxy").coverage)
        df_row = dq.to_dict()
    """

    station: str
    n_freq: int
    freq_min: Optional[float] = None
    freq_max: Optional[float] = None
    components: List[ComponentQuality] = field(default_factory=list)
    overall: QualityFlag = field(init=False)

    def __post_init__(self) -> None:
        self.overall = QualityFlag.worst(
            [c.flag for c in self.components if c.n_total > 0]
        ) if self.components else QualityFlag.MISSING

    # ------------------------------------------------------------------
    # Constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_site(cls, site: Any) -> "DataQuality":
        """Build a :class:`DataQuality` from a Site-like object.

        Accepts any object exposing ``.name``, ``.freq``, ``.z``,
        and ``.tipper`` (as per :class:`~pycsamt.site.base.SiteMixin`).
        """
        name  = getattr(site, "name", "?")
        freq  = _safe_array(getattr(site, "freq", None))
        z     = getattr(site, "z",     None)
        tip   = getattr(site, "tipper", None)

        n_freq    = int(freq.size)  if freq  is not None else 0
        freq_min  = float(freq.min()) if freq is not None and freq.size else None
        freq_max  = float(freq.max()) if freq is not None and freq.size else None

        comps: List[ComponentQuality] = []
        if z is not None:
            z_arr = np.asarray(z)
            for idx, cname in enumerate(_Z_COMPONENTS):
                col = _extract_component(z_arr, idx)
                comps.append(ComponentQuality.from_array(cname, col))
        if tip is not None:
            t_arr = np.asarray(tip)
            comps.append(ComponentQuality.from_array("Tipper", t_arr))

        return cls(
            station=name,
            n_freq=n_freq,
            freq_min=freq_min,
            freq_max=freq_max,
            components=comps,
        )

    @classmethod
    def from_edi(cls, edi_path: Any) -> "DataQuality":
        """Build from a spectra/impedance EDI file path."""
        from pycsamt.seg.edi import EDIFile  # noqa: PLC0415
        from pycsamt.site.base import Site   # noqa: PLC0415
        ed = EDIFile(str(edi_path))
        site = Site(ed)
        return cls.from_site(site)

    # ------------------------------------------------------------------
    # Accessors
    # ------------------------------------------------------------------

    def get(self, name: str) -> Optional[ComponentQuality]:
        """Return :class:`ComponentQuality` for *name*, or None."""
        for c in self.components:
            if c.name.lower() == name.lower():
                return c
        return None

    @property
    def z_components(self) -> List[ComponentQuality]:
        """Return only the Z-tensor component records."""
        return [c for c in self.components if c.name in _Z_COMPONENTS]

    @property
    def has_tipper(self) -> bool:
        t = self.get("Tipper")
        return t is not None and t.n_valid > 0

    @property
    def mean_coverage(self) -> float:
        """Mean coverage across all components that have data."""
        vals = [c.coverage for c in self.components if c.n_total > 0]
        return float(np.mean(vals)) if vals else 0.0

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------

    def summary(self) -> str:
        """Return a compact multi-line summary."""
        lines = [
            f"Station : {self.station}",
            f"  Frequencies : {self.n_freq}"
            + (f"  [{self.freq_min:.3g} – {self.freq_max:.3g} Hz]"
               if self.freq_min is not None else ""),
            f"  Overall     : {self.overall.value.upper()}",
        ]
        for c in self.components:
            bar = _bar(c.coverage, width=8)
            lines.append(
                f"  {c.name:<8} {bar} {c.pct_str:>5}  [{c.flag.value}]"
            )
        return "\n".join(lines)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "station":   self.station,
            "n_freq":    self.n_freq,
            "freq_min":  self.freq_min,
            "freq_max":  self.freq_max,
            "overall":   self.overall.value,
            "mean_coverage": round(self.mean_coverage, 4),
            "components": [c.to_dict() for c in self.components],
        }

    def __repr__(self) -> str:
        return (
            f"DataQuality({self.station!r}  {self.n_freq} freq"
            f"  overall={self.overall.value})"
        )


# ---------------------------------------------------------------------------
# Collection utilities
# ---------------------------------------------------------------------------

def assess_collection(sites: Any) -> List[DataQuality]:
    """Compute :class:`DataQuality` for every site in *sites*.

    Parameters
    ----------
    sites :
        Iterable of Site-like objects (e.g.
        :class:`~pycsamt.site.base.Sites`).

    Returns
    -------
    list of DataQuality
    """
    return [DataQuality.from_site(s) for s in sites]


def quality_dataframe(sites: Any, *, api: bool | None = None) -> Any:
    """Return a :class:`pandas.DataFrame` with one quality row per station.

    Columns: ``station``, ``n_freq``, ``freq_min``, ``freq_max``,
    ``overall``, ``mean_coverage``, and one ``cov_<component>`` column
    per impedance component.
    """
    try:
        import pandas as pd  # noqa: PLC0415
    except ImportError as exc:
        raise ImportError("pandas is required for quality_dataframe()") from exc

    rows = []
    for dq in assess_collection(sites):
        row: Dict[str, Any] = {
            "station":       dq.station,
            "n_freq":        dq.n_freq,
            "freq_min":      dq.freq_min,
            "freq_max":      dq.freq_max,
            "overall":       dq.overall.value,
            "mean_coverage": round(dq.mean_coverage, 4),
        }
        for c in dq.components:
            row[f"cov_{c.name}"] = round(c.coverage, 4)
            row[f"flag_{c.name}"] = c.flag.value
        rows.append(row)
    df = pd.DataFrame(rows)

    return maybe_wrap_frame(
        df,
        api=api,
        name="quality_dataframe",
        kind="metadata.quality",
        source=sites,
        description="Per-station data quality assessment table.",
    )


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _safe_array(arr: Any) -> Optional[np.ndarray]:
    if arr is None:
        return None
    a = np.asarray(arr)
    return a if a.size > 0 else None


def _extract_component(z: np.ndarray, idx: int) -> np.ndarray:
    """Extract component *idx* (0–3) from a Z array of any shape."""
    if z.ndim == 3 and z.shape[1:] == (2, 2):
        row, col = divmod(idx, 2)
        return z[:, row, col]
    if z.ndim == 2 and z.shape[1] == 4:
        return z[:, idx]
    return z.ravel()


def _bar(fraction: float, width: int = 8) -> str:
    """Unicode block-char progress bar."""
    fraction = max(0.0, min(1.0, float(fraction)))
    filled = round(fraction * width)
    return "█" * filled + "░" * (width - filled)
