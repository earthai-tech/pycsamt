"""AMT/CSAMT-specific edge quality-control metrics.

The generic :class:`~pycsamt.iot.edge.EdgeProcessor` checks coverage, RMS,
and spikes. Field electromagnetics needs more: powerline-harmonic
contamination, channel SNR, saturation/clipping, contact-resistance
proxies, resolvable frequency coverage, live spectra, impedance
stability, and sensor dropout. Those diagnostics live here.

All spectral routines are implemented with numpy only (a small internal
Welch estimator) so the subpackage imports without SciPy. If SciPy is
installed it is used automatically for a better PSD.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from ..api.property import PyCSAMTObject
from ..api.view import maybe_wrap_frame
from . import _common as _c

__all__ = [
    "HarmonicPeak",
    "PowerlineHarmonics",
    "FrequencyCoverage",
    "ImpedanceStability",
    "detect_powerline_harmonics",
    "estimate_channel_snr",
    "check_channel_saturation",
    "check_contact_resistance",
    "estimate_frequency_coverage",
    "compute_live_spectra",
    "assess_impedance_stability",
    "detect_sensor_dropout",
    "amt_edge_report",
    "amt_edge_table",
]


# ---------------------------------------------------------------------------
# signal preparation helpers
# ---------------------------------------------------------------------------
def _prep_signal(data: Any) -> np.ndarray:
    """Return a finite 1D float signal, interpolating interior NaNs."""
    x = np.asarray(data, dtype=float).ravel()
    if x.size == 0:
        return x
    finite = np.isfinite(x)
    if finite.all():
        return x
    if not finite.any():
        return np.empty(0)
    idx = np.arange(x.size)
    # Interpolate interior gaps; edge NaNs take the nearest finite value.
    x = np.interp(idx, idx[finite], x[finite])
    return x


def _positive_sample_rate(sample_rate: Any) -> float:
    return _c.as_positive(sample_rate, "sample_rate")


def _welch_psd(
    x: np.ndarray,
    fs: float,
    *,
    nperseg: int | None = None,
    detrend: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """One-sided power spectral density via Welch's method (numpy).

    Falls back to :func:`scipy.signal.welch` when SciPy is available.
    Returns ``(frequencies_hz, psd)``. Empty inputs yield empty arrays.
    """
    x = np.asarray(x, dtype=float).ravel()
    n = x.size
    if n < 2:
        return np.empty(0), np.empty(0)
    try:  # Prefer SciPy when present for a well-tested estimator.
        from scipy.signal import welch as _scipy_welch

        seg = int(nperseg or min(n, 256))
        seg = max(8, min(seg, n))
        freqs, psd = _scipy_welch(
            x, fs=fs, nperseg=seg,
            detrend="linear" if detrend else False,
        )
        return np.asarray(freqs), np.asarray(psd)
    except Exception:
        pass

    seg = int(nperseg or min(n, 256))
    seg = max(8, min(seg, n))
    step = max(1, seg // 2)
    window = np.hanning(seg)
    win_power = np.sum(window ** 2)
    starts = range(0, n - seg + 1, step)
    segments = list(starts) or [0]
    psd_acc = np.zeros(seg // 2 + 1)
    count = 0
    for start in segments:
        block = x[start:start + seg]
        if block.size < seg:
            block = np.pad(block, (0, seg - block.size))
        if detrend:
            block = block - np.mean(block)
        spectrum = np.fft.rfft(block * window)
        psd_block = (np.abs(spectrum) ** 2) / (fs * win_power)
        if psd_block.size > 1:
            psd_block[1:-1] *= 2.0  # one-sided scaling
        psd_acc += psd_block
        count += 1
    psd = psd_acc / max(count, 1)
    freqs = np.fft.rfftfreq(seg, d=1.0 / fs)
    return freqs, psd


def _spectrum_df(freqs: np.ndarray) -> float:
    """Return the frequency-bin spacing, or ``nan`` if undefined."""
    if freqs.size < 2:
        return float("nan")
    return float(freqs[1] - freqs[0])


def _band_power(
    freqs: np.ndarray,
    psd: np.ndarray,
    lo: float,
    hi: float,
) -> float:
    """Rectangular band power over ``[lo, hi]``.

    Uses ``sum(psd) * df`` rather than the trapezoidal rule so that a band
    containing only a single FFT bin still returns non-zero power.
    """
    if freqs.size == 0:
        return 0.0
    df = _spectrum_df(freqs)
    if not np.isfinite(df) or df <= 0:
        return 0.0
    mask = (freqs >= lo) & (freqs <= hi)
    if not mask.any():
        return 0.0
    return float(np.sum(psd[mask]) * df)


def _total_power(freqs: np.ndarray, psd: np.ndarray) -> float:
    """Rectangular integral of the whole one-sided PSD."""
    df = _spectrum_df(freqs)
    if not np.isfinite(df) or df <= 0:
        return 0.0
    return float(np.sum(psd) * df)


# ---------------------------------------------------------------------------
# powerline harmonics
# ---------------------------------------------------------------------------
@dataclass
class HarmonicPeak(PyCSAMTObject):
    """Contamination measured at one powerline harmonic."""

    order: int
    frequency_hz: float
    power_ratio: float
    flagged: bool

    def as_dict(self) -> dict[str, Any]:
        return dict(
            order=self.order,
            frequency_hz=self.frequency_hz,
            power_ratio=self.power_ratio,
            flagged=self.flagged,
        )


@dataclass
class PowerlineHarmonics(PyCSAMTObject):
    """Result of :func:`detect_powerline_harmonics`."""

    mains_hz: float
    peaks: list[HarmonicPeak] = field(default_factory=list)
    total_ratio: float = 0.0
    contaminated: bool = False

    @property
    def dominant(self) -> HarmonicPeak | None:
        """Return the strongest harmonic, if any were measured."""
        if not self.peaks:
            return None
        return max(self.peaks, key=lambda p: p.power_ratio)

    def as_dict(self) -> dict[str, Any]:
        return dict(
            mains_hz=self.mains_hz,
            total_ratio=self.total_ratio,
            contaminated=self.contaminated,
            n_flagged=sum(1 for p in self.peaks if p.flagged),
            dominant_hz=(self.dominant.frequency_hz if self.dominant else None),
            peaks=[p.as_dict() for p in self.peaks],
        )


def detect_powerline_harmonics(
    data: Any,
    sample_rate: float,
    *,
    mains_hz: float = 50.0,
    n_harmonics: int = 5,
    bandwidth_hz: float = 1.0,
    threshold_ratio: float = 0.05,
) -> PowerlineHarmonics:
    """Detect mains-frequency harmonics in a time series.

    Parameters
    ----------
    data : array-like
        Single-channel time series.
    sample_rate : float
        Sampling frequency in Hz.
    mains_hz : float
        Powerline fundamental (50 or 60 Hz typically).
    n_harmonics : int
        Number of harmonics (including the fundamental) to test.
    bandwidth_hz : float
        Half-width of the integration band around each harmonic.
    threshold_ratio : float
        Per-harmonic band-power fraction above which a harmonic is
        flagged as contaminating.

    Returns
    -------
    PowerlineHarmonics
    """
    fs = _positive_sample_rate(sample_rate)
    mains_hz = _c.as_positive(mains_hz, "mains_hz")
    n_harmonics = max(1, int(n_harmonics))
    bandwidth_hz = _c.as_positive(bandwidth_hz, "bandwidth_hz")
    threshold_ratio = _c.as_probability(threshold_ratio, "threshold_ratio")

    x = _prep_signal(data)
    freqs, psd = _welch_psd(x, fs)
    result = PowerlineHarmonics(mains_hz=mains_hz)
    if freqs.size == 0:
        return result
    total_power = _total_power(freqs, psd)
    df = _spectrum_df(freqs)
    nyquist = fs / 2.0
    total_ratio = 0.0
    for k in range(1, n_harmonics + 1):
        fc = k * mains_hz
        if fc >= nyquist:
            break
        # Ensure the integration band spans at least one bin either side.
        bw = max(bandwidth_hz, 1.5 * df) if np.isfinite(df) else bandwidth_hz
        band = _band_power(freqs, psd, fc - bw, fc + bw)
        if band == 0.0 and np.isfinite(df):
            idx = int(np.argmin(np.abs(freqs - fc)))
            band = float(psd[idx] * df)
        ratio = float(band / total_power) if total_power > 0 else 0.0
        flagged = ratio >= threshold_ratio
        total_ratio += ratio
        result.peaks.append(
            HarmonicPeak(
                order=k,
                frequency_hz=fc,
                power_ratio=ratio,
                flagged=flagged,
            )
        )
    result.total_ratio = float(total_ratio)
    result.contaminated = any(p.flagged for p in result.peaks)
    return result


# ---------------------------------------------------------------------------
# SNR, saturation, contact resistance
# ---------------------------------------------------------------------------
def estimate_channel_snr(
    data: Any,
    sample_rate: float | None = None,
    *,
    signal_band_hz: tuple[float, float] | None = None,
) -> float:
    """Estimate channel SNR in decibels.

    Two estimators are provided:

    * If ``sample_rate`` and ``signal_band_hz`` are given, SNR is the
      ratio of in-band to out-of-band spectral power.
    * Otherwise a time-domain estimate is used: the signal power is the
      variance of the series and the noise power is derived from the
      variance of first differences (a white-noise proxy).
    """
    x = _prep_signal(data)
    if x.size < 3:
        return float("nan")
    if sample_rate is not None and signal_band_hz is not None:
        fs = _positive_sample_rate(sample_rate)
        lo, hi = signal_band_hz
        lo = _c.as_nonnegative(lo, "signal_band_hz[0]")
        hi = _c.as_positive(hi, "signal_band_hz[1]")
        freqs, psd = _welch_psd(x, fs)
        if freqs.size == 0:
            return float("nan")
        in_band = _band_power(freqs, psd, lo, hi)
        total = _total_power(freqs, psd)
        noise = max(total - in_band, 1e-30)
        if in_band <= 0:
            return float("nan")
        return float(10.0 * np.log10(in_band / noise))
    signal_var = float(np.var(x))
    noise_var = float(np.var(np.diff(x)) / 2.0)
    if noise_var <= 0 or signal_var <= 0:
        return float("nan")
    return float(10.0 * np.log10(signal_var / noise_var))


def check_channel_saturation(
    data: Any,
    *,
    limit: float | None = None,
    max_clip_fraction: float = 0.01,
    tol: float = 1e-9,
) -> dict[str, Any]:
    """Detect ADC clipping / saturation in a channel.

    When ``limit`` is provided, samples with ``abs(x) >= limit`` count as
    saturated. Otherwise, samples equal to the observed min/max (within
    ``tol``) are treated as clipped, which catches rail-to-rail
    saturation without a known full-scale value.
    """
    x = np.asarray(data, dtype=float).ravel()
    finite = x[np.isfinite(x)]
    n = finite.size
    if n == 0:
        return dict(
            n_samples=0, n_clipped=0, clip_fraction=float("nan"),
            saturated=False, limit=limit,
        )
    max_clip_fraction = _c.as_probability(max_clip_fraction, "max_clip_fraction")
    if limit is not None:
        limit = _c.as_positive(limit, "limit")
        clipped = np.abs(finite) >= limit
    else:
        hi = float(np.max(finite))
        lo = float(np.min(finite))
        clipped = (np.abs(finite - hi) <= tol) | (np.abs(finite - lo) <= tol)
        # A non-degenerate signal that never revisits its extremes is fine.
        if hi == lo:
            clipped = np.ones_like(finite, dtype=bool)
    n_clipped = int(np.sum(clipped))
    frac = n_clipped / n
    return dict(
        n_samples=n,
        n_clipped=n_clipped,
        clip_fraction=float(frac),
        saturated=bool(frac > max_clip_fraction),
        limit=limit,
    )


def check_contact_resistance(
    ex: Any,
    ey: Any = None,
    *,
    sample_rate: float | None = None,
    noise_rms_threshold: float | None = None,
) -> dict[str, Any]:
    """Proxy assessment of electrode contact quality.

    True contact resistance requires a current injection measurement. On
    passive AMT electric channels, poor contact manifests as elevated
    low-frequency noise, large DC offset, and drift. This routine reports
    those proxies per channel (``ex`` and optional ``ey``) and flags a
    channel when its high-pass noise RMS exceeds ``noise_rms_threshold``
    (when provided) or when Ex/Ey noise is strongly imbalanced.
    """
    def _stats(sig: Any) -> dict[str, float] | None:
        x = _prep_signal(sig)
        if x.size < 3:
            return None
        dc = float(np.mean(x))
        # Drift: magnitude of a linear trend across the window.
        t = np.arange(x.size, dtype=float)
        slope = float(np.polyfit(t, x, 1)[0]) if np.ptp(t) > 0 else 0.0
        drift = abs(slope) * x.size
        detrended = x - (slope * t + float(np.mean(x)))
        noise_rms = float(np.sqrt(np.mean(detrended ** 2)))
        return dict(dc_offset=dc, drift=drift, noise_rms=noise_rms)

    ex_stats = _stats(ex)
    ey_stats = _stats(ey) if ey is not None else None
    out: dict[str, Any] = dict(ex=ex_stats, ey=ey_stats)

    flags: list[str] = []
    if noise_rms_threshold is not None:
        thr = _c.as_positive(noise_rms_threshold, "noise_rms_threshold")
        for name, st in (("ex", ex_stats), ("ey", ey_stats)):
            if st is not None and st["noise_rms"] > thr:
                flags.append(f"{name}_noise_above_threshold")
    imbalance = float("nan")
    if ex_stats is not None and ey_stats is not None:
        a, b = ex_stats["noise_rms"], ey_stats["noise_rms"]
        denom = min(a, b)
        if denom > 0:
            imbalance = max(a, b) / denom
            if imbalance > 5.0:
                flags.append("ex_ey_noise_imbalance")
    out.update(
        noise_imbalance=imbalance,
        flags=flags,
        ok=(len(flags) == 0 and ex_stats is not None),
    )
    return out


# ---------------------------------------------------------------------------
# frequency coverage
# ---------------------------------------------------------------------------
@dataclass
class FrequencyCoverage(PyCSAMTObject):
    """Result of :func:`estimate_frequency_coverage`."""

    sample_rate_hz: float
    nyquist_hz: float
    f_low_hz: float
    f_high_hz: float
    n_decades: float
    coverage_fraction: float = float("nan")
    missing_bands: list[tuple[float, float]] = field(default_factory=list)

    def as_dict(self) -> dict[str, Any]:
        return dict(
            sample_rate_hz=self.sample_rate_hz,
            nyquist_hz=self.nyquist_hz,
            f_low_hz=self.f_low_hz,
            f_high_hz=self.f_high_hz,
            n_decades=self.n_decades,
            coverage_fraction=self.coverage_fraction,
            missing_bands=[list(b) for b in self.missing_bands],
        )


def estimate_frequency_coverage(
    timeseries: Any,
    sample_rate: float,
    *,
    target_bands: Iterable[tuple[float, float]] | None = None,
    snr_floor_db: float = 6.0,
) -> FrequencyCoverage:
    """Estimate the resolvable frequency band of a recording.

    The PSD noise floor is taken as its median. Frequencies whose power
    exceeds the floor by ``snr_floor_db`` are considered resolved; the
    lowest and highest such frequencies define the covered band. When
    ``target_bands`` are supplied, the fraction that falls inside the
    covered band is reported along with any missing bands.
    """
    fs = _positive_sample_rate(sample_rate)
    snr_floor_db = _c.as_nonnegative(snr_floor_db, "snr_floor_db")
    x = _prep_signal(timeseries)
    nyquist = fs / 2.0
    coverage = FrequencyCoverage(
        sample_rate_hz=fs,
        nyquist_hz=nyquist,
        f_low_hz=float("nan"),
        f_high_hz=float("nan"),
        n_decades=float("nan"),
    )
    freqs, psd = _welch_psd(x, fs)
    if freqs.size == 0:
        return coverage
    positive = freqs > 0
    freqs, psd = freqs[positive], psd[positive]
    if freqs.size == 0 or not np.any(psd > 0):
        return coverage
    floor = float(np.median(psd[psd > 0]))
    ratio_db = 10.0 * np.log10(np.maximum(psd, 1e-30) / max(floor, 1e-30))
    resolved = ratio_db >= snr_floor_db
    if not resolved.any():
        return coverage
    f_lo = float(freqs[resolved][0])
    f_hi = float(freqs[resolved][-1])
    coverage.f_low_hz = f_lo
    coverage.f_high_hz = f_hi
    coverage.n_decades = (
        float(np.log10(f_hi / f_lo)) if f_lo > 0 and f_hi > 0 else float("nan")
    )
    if target_bands is not None:
        bands = [
            (
                _c.as_positive(lo, "target_band_low"),
                _c.as_positive(hi, "target_band_high"),
            )
            for lo, hi in target_bands
        ]
        if bands:
            covered = [b for b in bands if b[0] >= f_lo and b[1] <= f_hi]
            coverage.missing_bands = [b for b in bands if b not in covered]
            coverage.coverage_fraction = len(covered) / len(bands)
    return coverage


def compute_live_spectra(
    data: Any,
    sample_rate: float,
    *,
    nperseg: int | None = None,
) -> dict[str, np.ndarray]:
    """Return ``{"frequency_hz": ..., "psd": ...}`` for live display."""
    fs = _positive_sample_rate(sample_rate)
    x = _prep_signal(data)
    freqs, psd = _welch_psd(x, fs, nperseg=nperseg)
    return {"frequency_hz": freqs, "psd": psd}


# ---------------------------------------------------------------------------
# impedance stability
# ---------------------------------------------------------------------------
@dataclass
class ImpedanceStability(PyCSAMTObject):
    """Result of :func:`assess_impedance_stability`."""

    n_windows: int
    cv_magnitude: float
    phase_std_deg: float
    stable: bool

    def as_dict(self) -> dict[str, Any]:
        return dict(
            n_windows=self.n_windows,
            cv_magnitude=self.cv_magnitude,
            phase_std_deg=self.phase_std_deg,
            stable=self.stable,
        )


def assess_impedance_stability(
    z_windows: Any,
    *,
    max_cv: float = 0.15,
    max_phase_std_deg: float = 10.0,
) -> ImpedanceStability:
    """Assess the stability of per-window impedance estimates.

    Parameters
    ----------
    z_windows : array-like of complex
        Impedance estimates, shape ``(n_windows,)`` or
        ``(n_windows, n_freq)``. Real inputs are treated as magnitudes
        with zero phase.
    max_cv : float
        Maximum coefficient of variation of ``|Z|`` for a stable result.
    max_phase_std_deg : float
        Maximum phase standard deviation (degrees) for a stable result.
    """
    z = np.asarray(z_windows, dtype=complex)
    if z.ndim == 1:
        z = z.reshape(-1, 1)
    if z.size == 0 or z.shape[0] < 2:
        return ImpedanceStability(
            n_windows=int(z.shape[0]) if z.ndim >= 1 else 0,
            cv_magnitude=float("nan"),
            phase_std_deg=float("nan"),
            stable=False,
        )
    mag = np.abs(z)
    phase = np.degrees(np.angle(z))
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_mag = np.mean(mag, axis=0)
        cv = np.where(mean_mag > 0, np.std(mag, axis=0) / mean_mag, np.nan)
    cv_magnitude = float(np.nanmean(cv)) if np.any(np.isfinite(cv)) else float(
        "nan"
    )
    phase_std = float(np.nanmean(np.std(phase, axis=0)))
    stable = bool(
        np.isfinite(cv_magnitude)
        and cv_magnitude <= max_cv
        and phase_std <= max_phase_std_deg
    )
    return ImpedanceStability(
        n_windows=int(z.shape[0]),
        cv_magnitude=cv_magnitude,
        phase_std_deg=phase_std,
        stable=stable,
    )


# ---------------------------------------------------------------------------
# sensor dropout
# ---------------------------------------------------------------------------
def detect_sensor_dropout(
    data: Any,
    *,
    min_flat_run: int = 8,
    flat_tol: float = 1e-12,
) -> dict[str, Any]:
    """Detect NaN gaps and stuck-value (flatline) runs in a channel.

    Returns counts for NaN samples and the longest run of (near-)constant
    consecutive samples, which typically indicates a disconnected or
    stuck sensor.
    """
    x = np.asarray(data, dtype=float).ravel()
    n = x.size
    min_flat_run = max(2, int(min_flat_run))
    if n == 0:
        return dict(
            n_samples=0, n_nan=0, nan_fraction=float("nan"),
            longest_flat_run=0, n_flat_runs=0, dropout=False,
        )
    nan_mask = ~np.isfinite(x)
    n_nan = int(np.sum(nan_mask))

    longest = 1
    current = 1
    n_flat_runs = 0
    flagged_run = False
    for i in range(1, n):
        same = (
            np.isfinite(x[i]) and np.isfinite(x[i - 1])
            and abs(x[i] - x[i - 1]) <= flat_tol
        )
        if same:
            current += 1
        else:
            if current >= min_flat_run:
                n_flat_runs += 1
                flagged_run = True
            current = 1
        longest = max(longest, current)
    if current >= min_flat_run:
        n_flat_runs += 1
        flagged_run = True

    return dict(
        n_samples=n,
        n_nan=n_nan,
        nan_fraction=float(n_nan / n),
        longest_flat_run=int(longest),
        n_flat_runs=int(n_flat_runs),
        dropout=bool(flagged_run or n_nan > 0),
    )


# ---------------------------------------------------------------------------
# aggregation
# ---------------------------------------------------------------------------
def amt_edge_report(
    data: Any,
    sample_rate: float,
    *,
    mains_hz: float = 50.0,
    signal_band_hz: tuple[float, float] | None = None,
) -> dict[str, Any]:
    """Run the core AMT edge diagnostics on one channel and collate them."""
    harmonics = detect_powerline_harmonics(
        data, sample_rate, mains_hz=mains_hz
    )
    coverage = estimate_frequency_coverage(data, sample_rate)
    return dict(
        snr_db=estimate_channel_snr(
            data, sample_rate, signal_band_hz=signal_band_hz
        ),
        powerline=harmonics.as_dict(),
        saturation=check_channel_saturation(data),
        dropout=detect_sensor_dropout(data),
        frequency_coverage=coverage.as_dict(),
    )


def amt_edge_table(
    reports: dict[str, dict[str, Any]] | Iterable[tuple[str, dict[str, Any]]],
    *,
    api: bool | None = None,
) -> Any:
    """Flatten one or more :func:`amt_edge_report` results into a table.

    Accepts a ``{channel: report}`` mapping (or ``(channel, report)``
    pairs) and returns one row per channel with the headline metrics.
    """
    items = (
        list(reports.items()) if isinstance(reports, dict)
        else list(reports)
    )
    rows: list[dict[str, Any]] = []
    for channel, report in items:
        powerline = report.get("powerline", {})
        saturation = report.get("saturation", {})
        dropout = report.get("dropout", {})
        coverage = report.get("frequency_coverage", {})
        rows.append(dict(
            channel=str(channel).lower(),
            snr_db=report.get("snr_db"),
            powerline_contaminated=powerline.get("contaminated"),
            powerline_total_ratio=powerline.get("total_ratio"),
            saturated=saturation.get("saturated"),
            clip_fraction=saturation.get("clip_fraction"),
            dropout=dropout.get("dropout"),
            nan_fraction=dropout.get("nan_fraction"),
            f_low_hz=coverage.get("f_low_hz"),
            f_high_hz=coverage.get("f_high_hz"),
            n_decades=coverage.get("n_decades"),
        ))
    df = pd.DataFrame.from_records(rows)
    return maybe_wrap_frame(
        df,
        api=api,
        name="iot_amt_edge_table",
        kind="iot.edge.amt",
        source=items,
        description="AMT/CSAMT edge quality-control metrics by channel.",
    )
